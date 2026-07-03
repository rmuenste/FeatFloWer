MODULE EL_TRANSFER

  USE ISO_C_BINDING, ONLY: C_DOUBLE, C_INT, C_SHORT
  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE PP3D_MPI, ONLY: myid, master, showid, MPI_COMM_SUBS
  ! Global multigrid level from COMMON /MGPAR/ (def_FEAT). Imported under an alias
  ! so it does not collide with the per-call `ilev` dummy argument. E013Sum3 reads
  ! this global ILEV to pick MGE013(ILEV); it must be set to the EL coupling level
  ! before the feedback assembly (see EL_PARTICLE_MESH_PASS).
  USE def_FEAT, ONLY: el_mg_ilev => ILEV
  ! Global mesh dimensions from COMMON /TRIAD/ (def_FEAT). NDFGL (called from
  ! EL_INTEGRATE_PARTICLE / EL_DEPOSIT_PARTICLE) derives Q2 global DOF indices
  ! from these offsets; like ILEV they are left at whatever level the previous
  ! solver phase set, so they must be pinned to the EL coupling level for the pass
  ! (see EL_PARTICLE_MESH_PASS). Aliased to avoid colliding with the local `nel`.
  USE def_FEAT, ONLY: el_tr_nel => NEL, el_tr_nvt => NVT, &
                      el_tr_net => NET, el_tr_nat => NAT
  USE DEM_QUERY, ONLY: tParticleData
#ifdef HAVE_PE
  USE DEM_QUERY, ONLY: numLocalParticles, getAllParticles, setForcesMapped, &
                       getElSubsteps, getElStepsize
#endif
  USE EL_CONFIG, ONLY: el_eps_f_min, el_eps_f_relax, &
                       el_apply_particle_forces, el_apply_fluid_feedback, &
                       el_write_diagnostics, el_momentum_audit_freq, &
                       el_drag_semi_implicit
  USE EL_FORCES, ONLY: tElForceResult, EL_COMPUTE_PARTICLE_FORCES, &
                       EL_SEMI_IMPLICIT_DRAG_IMPULSE, EL_PARTICLE_VOLUME
  USE EL_FIELDS, ONLY: el_field_data, EL_ENSURE_FIELDS, EL_BEGIN_FIELD_UPDATE, &
                       EL_CAPTURE_FLUID_FEEDBACK_SOURCE
  USE EL_HALO, ONLY: tElParticleRecord, EL_REFRESH_COUPLING_HALO, &
                     EL_REDUCE_TO_OWNERS, EL_BROADCAST_FROM_OWNERS
  USE EL_QUADRATURE, ONLY: EL_SAMPLE_SIZE, EL_SAMPLE_NORMALIZATION, &
                           EL_INTEGRATE_PARTICLE, &
                           EL_DEPOSIT_PARTICLE

  IMPLICIT NONE

#ifdef HAVE_PE
  INTERFACE
    SUBROUTINE step_el_particles_c() BIND(C, NAME='step_el_frozen_trace_')
    END SUBROUTINE step_el_particles_c

    SUBROUTINE sync_el_forces_c() BIND(C, NAME='sync_el_frozen_forces_')
    END SUBROUTINE sync_el_forces_c

    SUBROUTINE set_el_hydro_state_c(bytes, drag_b, u_f, f_other, active) &
      BIND(C, NAME='set_el_hydro_state')
      USE ISO_C_BINDING, ONLY: C_DOUBLE, C_INT, C_SHORT
      INTEGER(C_SHORT), INTENT(IN) :: bytes(8)
      REAL(C_DOUBLE), INTENT(IN) :: drag_b
      REAL(C_DOUBLE), INTENT(IN) :: u_f(3), f_other(3)
      INTEGER(C_INT), INTENT(IN) :: active
    END SUBROUTINE set_el_hydro_state_c
  END INTERFACE
#endif

CONTAINS

  SUBROUTINE EL_ASSERT_PE_TIMESTEP(dt, mfile, context)

    REAL*8, INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: mfile
    CHARACTER(LEN=*), INTENT(IN) :: context

#ifdef HAVE_PE
    REAL*8 :: pe_dt, tol
    INTEGER :: ierr

    pe_dt = getElStepsize()
    tol = 1.0d-10*MAX(ABS(dt),1.0d-300)
    IF (ABS(dt-pe_dt).GT.tol) THEN
      IF (myid.EQ.showid) THEN
        WRITE(*,'(A,A,A,2ES24.16,A,ES14.6)') &
          'Fatal E-L timestep mismatch in ', TRIM(context), &
          ': CFD dt, PE stepsize_ = ', dt, pe_dt, ' tolerance= ', tol
        WRITE(mfile,'(A,A,A,2ES24.16,A,ES14.6)') &
          'Fatal E-L timestep mismatch in ', TRIM(context), &
          ': CFD dt, PE stepsize_ = ', dt, pe_dt, ' tolerance= ', tol
      END IF
      CALL MPI_Abort(MPI_COMM_SUBS, 905, ierr)
    END IF
#endif

  END SUBROUTINE EL_ASSERT_PE_TIMESTEP

  SUBROUTINE EL_PARTICLE_MESH_PASS(mesh, ilev, val_u, val_v, val_w, val_p, &
                                   density, viscosity, gravity, dt, mfile, istep)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, mfile, istep
    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:), val_p(:)
    REAL*8, INTENT(IN) :: density, viscosity, gravity(3), dt

    TYPE(tParticleData), ALLOCATABLE :: owned_particles(:)
    TYPE(tElParticleRecord), ALLOCATABLE :: records(:)
    REAL*8, ALLOCATABLE :: local_samples(:,:), owned_samples(:,:)
    REAL*8, ALLOCATABLE :: owned_result(:,:), record_result(:,:)
    REAL*8 :: expected_volume
#ifdef HAVE_PE
    REAL*8 :: f_other(3)
#endif
    TYPE(tElForceResult) :: force_result
    REAL*8 :: local_deposited, global_deposited, local_expected, global_expected
    REAL*8 :: local_feedback(3), global_feedback(3)
    REAL*8 :: local_drag(3), global_drag(3), local_lift(3), global_lift(3)
    REAL*8 :: drag_force_eff(3), feedback_used(3), m_p
    REAL*8 :: local_drag_impl(3), global_drag_impl(3)
    INTEGER :: n_sub
    REAL*8 :: local_pressure(3), global_pressure(3), local_grav_buoy(3)
    REAL*8 :: global_grav_buoy(3), local_particle_total(3), global_particle_total(3)
    REAL*8 :: local_rhs(3), global_rhs(3), feedback_residual(3)
    REAL*8 :: feedback_residual_norm, volume_rel_error
    INTEGER :: n_owned, n_records, i, ierr, ndof, nel
    INTEGER :: local_record_rank, global_record_ranks
    INTEGER :: global_records, max_records
    INTEGER :: saved_mg_ilev
    INTEGER :: saved_tr_nel, saved_tr_nvt, saved_tr_net, saved_tr_nat
    INTEGER :: global_owned
#ifdef HAVE_PE
    INTEGER :: hydro_active
#endif
    INTEGER :: clipped_local, clipped_global
    LOGICAL :: advance_history
    CHARACTER(LEN=16) :: coupling_mode

    IF (myid.EQ.master) RETURN

#ifdef PE_SERIAL_MODE
    IF (myid.EQ.showid) THEN
      WRITE(*,'(A)') 'Fatal: q2p1_el_pipeflow requires distributed PE ownership.'
      WRITE(*,'(A)') 'PE_SERIAL_MODE duplicates all particles on every CFD rank and is not scalable.'
      WRITE(*,'(A)') 'Reconfigure with USE_PE=ON and USE_PE_SERIAL_MODE=OFF.'
      WRITE(mfile,'(A)') 'Fatal: E-L transfer does not support PE_SERIAL_MODE.'
    END IF
    CALL MPI_Abort(MPI_COMM_SUBS, 904, ierr)
#endif

    nel = mesh%level(ilev)%nel
    ndof = mesh%level(ilev)%nvt + mesh%level(ilev)%net + &
           mesh%level(ilev)%nat + nel
    CALL EL_ENSURE_FIELDS(nel, ndof)

    ! Pin the global /TRIAD/ mesh offsets to the EL coupling level for the whole
    ! pass. NDFGL (inside EL_INTEGRATE_PARTICLE and EL_DEPOSIT_PARTICLE) builds Q2
    ! global DOF indices as e.g. centre = NVT+NET+NAT+iel. If NVT/NET/NAT are left
    ! at another level, those indices no longer match the EL coupling arrays.
    ! Restored at the end of the pass.
    saved_tr_nel = el_tr_nel
    saved_tr_nvt = el_tr_nvt
    saved_tr_net = el_tr_net
    saved_tr_nat = el_tr_nat
    el_tr_nel = nel
    el_tr_nvt = mesh%level(ilev)%nvt
    el_tr_net = mesh%level(ilev)%net
    el_tr_nat = mesh%level(ilev)%nat

    advance_history = istep.LT.0
    CALL EL_BEGIN_FIELD_UPDATE(advance_history)

    n_owned = 0
#ifdef HAVE_PE
    n_owned = numLocalParticles()
    IF (n_owned.GT.0) THEN
      ALLOCATE(owned_particles(n_owned))
      CALL getAllParticles(owned_particles)
    ELSE
      ALLOCATE(owned_particles(0))
    END IF
#else
    ALLOCATE(owned_particles(0))
#endif

    CALL EL_REFRESH_COUPLING_HALO(mesh, ilev, owned_particles, records, mfile)
    n_records = SIZE(records)
    ALLOCATE(local_samples(EL_SAMPLE_SIZE,n_records))
    ALLOCATE(owned_samples(EL_SAMPLE_SIZE,n_owned))
    ALLOCATE(owned_result(5,n_owned), record_result(5,n_records))
    local_samples = 0.0d0
    owned_samples = 0.0d0
    owned_result = 0.0d0
    record_result = 0.0d0
    local_feedback = 0.0d0
    local_drag = 0.0d0
    local_lift = 0.0d0
    local_pressure = 0.0d0
    local_grav_buoy = 0.0d0
    local_particle_total = 0.0d0
    local_drag_impl = 0.0d0
    coupling_mode = 'explicit'
    IF (el_drag_semi_implicit) coupling_mode = 'semi_implicit'
#ifdef HAVE_PE
    IF (advance_history) CALL EL_ASSERT_PE_TIMESTEP(dt, mfile, 'pre-advance pass')
    n_sub = getElSubsteps()
#else
    n_sub = 1
#endif

    DO i=1,n_records
      CALL EL_INTEGRATE_PARTICLE(mesh, ilev, val_u, val_v, val_w, val_p, &
        el_field_data%epsilon_f, records(i), local_samples(:,i))
    END DO
    CALL EL_REDUCE_TO_OWNERS(local_samples, owned_samples)

    DO i=1,n_owned
      IF (owned_samples(1,i).LE.0.0d0) THEN
        WRITE(*,'(A,I0,A,I0,A,3ES14.6)') 'Fatal E-L transfer error on rank ', &
          myid, ': particle ', owned_particles(i)%systemIdx, &
          ' has no positive fluid support at ', owned_particles(i)%position
        CALL MPI_Abort(MPI_COMM_SUBS, 901, ierr)
      END IF
      CALL EL_FORCES_FROM_PARTICLE(owned_particles(i), owned_samples(:,i), &
        density, viscosity, gravity, force_result)
      owned_result(1,i) = owned_samples(EL_SAMPLE_NORMALIZATION,i)
      owned_result(5,i) = force_result%drag_B
      ! Default (post-advance diagnostic pass): explicit drag + lift. The
      ! pre-advance pass below replaces the drag part with the sub-cycled
      ! implicit drag impulse PE actually applies (issue-D fix), then the
      ! finalised feedback_used is deposited after the PE-arming block.
      drag_force_eff = 0.0d0
      feedback_used = force_result%feedback_force
      local_drag = local_drag + force_result%drag
      local_lift = local_lift + force_result%lift
      local_pressure = local_pressure + force_result%pressure
      local_grav_buoy = local_grav_buoy + force_result%grav_buoy
      local_particle_total = local_particle_total + force_result%particle_total
#ifdef HAVE_PE
      ! Only the pre-advance pass (advance_history, istep<0) is allowed to arm
      ! PE: it runs immediately before EL_ADVANCE_PARTICLES. The post-step
      ! diagnostic pass (istep>0) must NOT push a hydro state, otherwise it
      ! would leave PE armed with stale fluid forcing between this step and the
      ! next pre-advance pass. The PE side clears all hydro states at the end of
      ! each PE step, so a body is integrated semi-implicitly only when it was
      ! freshly armed for the step that is about to run.
      IF (el_apply_particle_forces .AND. advance_history) THEN
        owned_particles(i)%force = force_result%particle_total
        owned_particles(i)%torque = 0.0d0
        CALL setForcesMapped(owned_particles(i))
        f_other = force_result%pressure + force_result%lift + &
                  force_result%grav_buoy
        ! Issue-D fix: feed the fluid the drag impulse PE will actually apply over
        ! its n_sub sub-cycled semi-implicit steps (drag_force_eff), plus explicit
        ! lift -- not the explicit old-velocity drag -- so the momentum the fluid
        ! gains equals the drag momentum the particle loses. Lift is in f_other
        ! (shaping the implicit velocity, hence excluded from drag_force_eff) AND
        ! spread explicitly; the Newton pair holds for lift with no double count.
        m_p = EL_PARTICLE_VOLUME(owned_particles(i)%radius) * &
              owned_particles(i)%density
        CALL EL_SEMI_IMPLICIT_DRAG_IMPULSE(owned_particles(i)%velocity, m_p, &
          force_result%drag_B, force_result%u_f, f_other, dt, n_sub, drag_force_eff)
        IF (el_drag_semi_implicit) THEN
          ! Fluid semi-implicit split around PE's actual sub-cycled drag impulse:
          !   target fluid force = -drag_force_eff - lift
          !   target = (B*u_f - drag_force_eff - lift) - B*u_new
          ! EL_DEPOSIT_PARTICLE subtracts feedback_used, so pass
          !   feedback_used = drag_force_eff + lift - B*u_f.
          feedback_used = drag_force_eff + force_result%lift - &
            force_result%drag_B*force_result%u_f
        ELSE
          feedback_used = drag_force_eff + force_result%lift
        END IF
        hydro_active = 1
        CALL set_el_hydro_state_c(owned_particles(i)%bytes, &
          force_result%drag_B, force_result%u_f, f_other, hydro_active)
      END IF
#endif
      owned_result(2:4,i) = feedback_used
      local_feedback = local_feedback + feedback_used
      local_drag_impl = local_drag_impl + drag_force_eff
    END DO

    CALL EL_BROADCAST_FROM_OWNERS(owned_result, record_result)
    DO i=1,n_records
      CALL EL_DEPOSIT_PARTICLE(mesh, ilev, records(i), record_result(1,i), &
        record_result(2:4,i), el_field_data, drag_b=record_result(5,i))
    END DO

    local_rhs(1) = SUM(el_field_data%force_rhs(1,:))
    local_rhs(2) = SUM(el_field_data%force_rhs(2,:))
    local_rhs(3) = SUM(el_field_data%force_rhs(3,:))

    ! Capture the fluid feedback source from the DISTRIBUTED force_rhs, i.e.
    ! BEFORE E013Sum3. The momentum RHS this source is added to (QuadSc%defU,
    ! alongside AddGravForce) is in distributed/additive form and is summed once
    ! by the fluid solver. If we captured after E013Sum3 the source would already
    ! be consistent (summed), and the later solver sum would multiply it
    ! by the sharing multiplicity on partition-boundary DOFs (double-counting).
    ! The diagnostic force_rhs is still summed below for reporting/output.
    IF (advance_history) CALL EL_CAPTURE_FLUID_FEEDBACK_SOURCE()
    ! The parallel assembly (E013Sum3) of the spread force is only needed when the
    ! feedback is actually added to the fluid momentum, which is itself gated on
    ! el_apply_fluid_feedback (QuadSc_main); in one-way mode the summed force_rhs
    ! is never consumed, so skip it.
    !
    ! E013Sum3 picks the communication structure MGE013(ILEV) using the GLOBAL
    ! multigrid level (COMMON /MGPAR/), not our `ilev` argument. On entry the
    ! global level is whatever the previous solver phase left it at (e.g. a coarse
    ! pressure level), so MGE013(ILEV) would not match the force_rhs DOF layout
    ! (sized at mesh%level(ilev)) and the sum would index the wrong/!unallocated
    ! structure. Pin the global level to the EL coupling level for the sum, then
    ! restore it so we do not perturb the surrounding solver state.
    IF (el_apply_fluid_feedback) THEN
      saved_mg_ilev = el_mg_ilev
      el_mg_ilev = ilev
      CALL E013Sum3(el_field_data%force_rhs(1,:),el_field_data%force_rhs(2,:), &
                    el_field_data%force_rhs(3,:))
      el_mg_ilev = saved_mg_ilev
    END IF
    CALL EL_FINALIZE_VOID_FRACTION(dt, clipped_local)
    CALL MPI_Allreduce(clipped_local, clipped_global, 1, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_SUBS, ierr)

    IF (clipped_global.GT.0 .AND. myid.EQ.showid) THEN
      WRITE(*,'(A,I0,A,I0)') 'WARNING: E-L void-fraction clipping occurred in ', &
        clipped_global, ' cells at step ', ABS(istep)
      WRITE(mfile,'(A,I0,A,I0)') 'WARNING: E-L void-fraction clipping occurred in ', &
        clipped_global, ' cells at step ', ABS(istep)
    END IF

    ! Slice dvol to nel: it carries the FEAT total-volume slot at nel+1
    ! (see el_diagnostics volume integrals); alpha_p is nel-sized.
    local_deposited = SUM(el_field_data%alpha_p(1:nel)* &
                          mesh%level(ilev)%dvol(1:nel))
    local_expected = 0.0d0
    DO i=1,n_owned
      local_expected = local_expected + 4.0d0*ACOS(-1.0d0)* &
        owned_particles(i)%radius**3/3.0d0
    END DO
    CALL MPI_Allreduce(local_deposited, global_deposited, 1, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_expected, global_expected, 1, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_feedback, global_feedback, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_drag, global_drag, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_drag_impl, global_drag_impl, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_lift, global_lift, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_pressure, global_pressure, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_grav_buoy, global_grav_buoy, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_particle_total, global_particle_total, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_rhs, global_rhs, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(n_owned, global_owned, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    IF (el_write_diagnostics) THEN
      local_record_rank = 0
      IF (n_records.GT.0) local_record_rank = 1
      CALL MPI_Allreduce(local_record_rank, global_record_ranks, 1, &
        MPI_INTEGER, MPI_SUM, MPI_COMM_SUBS, ierr)
      CALL MPI_Allreduce(n_records, global_records, 1, MPI_INTEGER, &
        MPI_SUM, MPI_COMM_SUBS, ierr)
      CALL MPI_Allreduce(n_records, max_records, 1, MPI_INTEGER, &
        MPI_MAX, MPI_COMM_SUBS, ierr)
    END IF

    feedback_residual = global_rhs + global_feedback
    feedback_residual_norm = SQRT(SUM(feedback_residual**2))

    IF (myid.EQ.showid .AND. el_write_diagnostics .AND. istep.GE.0) THEN
      expected_volume = MAX(global_expected,TINY(1.0d0))
      volume_rel_error = ABS(global_deposited-global_expected)/expected_volume
      WRITE(*,'(A,I0,A,I0,A,3ES14.6)') 'EL step ', istep, &
        ': owned particles=', global_owned, ', volume deposited/expected/error=', &
        global_deposited, global_expected, volume_rel_error
      WRITE(mfile,'(A,I0,A,3ES14.6)') 'EL step ', istep, &
        ': volume deposited/expected/error=', global_deposited, global_expected, &
        volume_rel_error
      WRITE(*,'(A,I0,A,2ES14.6,A,ES14.6)') &
        'EL_VOLUME_CONSERVATION step= ', istep, ' deposited_expected= ', &
        global_deposited, global_expected, ' rel_error= ', volume_rel_error
      WRITE(mfile,'(A,I0,A,2ES14.6,A,ES14.6)') &
        'EL_VOLUME_CONSERVATION step= ', istep, ' deposited_expected= ', &
        global_deposited, global_expected, ' rel_error= ', volume_rel_error
      WRITE(*,'(A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6)') &
        'EL_FORCE_BUDGET step= ', istep, ' drag= ', global_drag, &
        ' pressure= ', global_pressure, ' lift= ', global_lift, &
        ' grav_buoy= ', global_grav_buoy, ' total= ', global_particle_total
      WRITE(mfile,'(A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6)') &
        'EL_FORCE_BUDGET step= ', istep, ' drag= ', global_drag, &
        ' pressure= ', global_pressure, ' lift= ', global_lift, &
        ' grav_buoy= ', global_grav_buoy, ' total= ', global_particle_total
      WRITE(*,'(A,I0,A,I0,A,I0,A,I0)') 'EL_HALO_RECORDS step= ', &
        istep, ' record_ranks= ', global_record_ranks, &
        ' total_records= ', global_records, ' max_records= ', max_records
      WRITE(mfile,'(A,I0,A,I0,A,I0,A,I0)') 'EL_HALO_RECORDS step= ', &
        istep, ' record_ranks= ', global_record_ranks, &
        ' total_records= ', global_records, ' max_records= ', max_records
    END IF

    IF (myid.EQ.showid .AND. istep.GE.0 .AND. &
        MOD(ABS(istep), el_momentum_audit_freq).EQ.0) THEN
      WRITE(*,'(A,I0,A,3ES14.6,A,3ES14.6,A,ES14.6)') &
        'EL_FEEDBACK_CONSERVATION step= ', istep, ' mode= '//TRIM(coupling_mode)// &
        ' rhs= ', global_rhs, ' feedback= ', global_feedback, ' residual= ', &
        feedback_residual_norm
      WRITE(mfile,'(A,I0,A,3ES14.6,A,3ES14.6,A,ES14.6)') &
        'EL_FEEDBACK_CONSERVATION step= ', istep, ' mode= '//TRIM(coupling_mode)// &
        ' rhs= ', global_rhs, ' feedback= ', global_feedback, ' residual= ', &
        feedback_residual_norm
    END IF

    ! Restore the global /TRIAD/ offsets we pinned to the coupling level on entry.
    el_tr_nel = saved_tr_nel
    el_tr_nvt = saved_tr_nvt
    el_tr_net = saved_tr_net
    el_tr_nat = saved_tr_nat

    DEALLOCATE(owned_particles, records, local_samples, owned_samples)
    DEALLOCATE(owned_result, record_result)

  END SUBROUTINE EL_PARTICLE_MESH_PASS

  SUBROUTINE EL_FORCES_FROM_PARTICLE(particle, sample, density, viscosity, &
                                     gravity, result)

    TYPE(tParticleData), INTENT(IN) :: particle
    REAL*8, INTENT(IN) :: sample(EL_SAMPLE_SIZE)
    REAL*8, INTENT(IN) :: density, viscosity, gravity(3)
    TYPE(tElForceResult), INTENT(OUT) :: result
    REAL*8 :: state(8)

    state(1:3) = particle%position
    state(4:6) = particle%velocity
    state(7) = particle%radius
    state(8) = particle%density
    CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, density, viscosity, &
                                    gravity, result)

  END SUBROUTINE EL_FORCES_FROM_PARTICLE

  SUBROUTINE EL_ADVANCE_PARTICLES()

#ifdef HAVE_PE
    IF (myid.NE.master) THEN
      CALL sync_el_forces_c()
      CALL step_el_particles_c()
    END IF
#endif

  END SUBROUTINE EL_ADVANCE_PARTICLES

  SUBROUTINE EL_FINALIZE_VOID_FRACTION(dt, clipped_count)

    REAL*8, INTENT(IN) :: dt
    INTEGER, INTENT(OUT) :: clipped_count
    REAL*8, ALLOCATABLE :: epsilon_new(:)

    clipped_count = 0
    IF (.NOT.ALLOCATED(el_field_data%epsilon_f)) RETURN
    IF (SIZE(el_field_data%epsilon_f).EQ.0) RETURN

    ALLOCATE(epsilon_new(SIZE(el_field_data%epsilon_f)))
    clipped_count = COUNT(1.0d0-el_field_data%alpha_p.LT.el_eps_f_min)
    epsilon_new = MAX(1.0d0-el_field_data%alpha_p,el_eps_f_min)
    el_field_data%epsilon_f = el_eps_f_relax*el_field_data%epsilon_f_old + &
      (1.0d0-el_eps_f_relax)*epsilon_new
    IF (dt.GT.0.0d0) THEN
      el_field_data%deps_f_dt = &
        (el_field_data%epsilon_f-el_field_data%epsilon_f_old)/dt
    ELSE
      el_field_data%deps_f_dt = 0.0d0
    END IF
    DEALLOCATE(epsilon_new)

  END SUBROUTINE EL_FINALIZE_VOID_FRACTION

END MODULE EL_TRANSFER
