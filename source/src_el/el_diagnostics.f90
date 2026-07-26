MODULE EL_DIAGNOSTICS

  USE MPI
  USE PP3D_MPI, ONLY: myid, showid, MPI_COMM_SUBS
  USE DEM_QUERY, ONLY: tParticleData
  USE TYPES, ONLY: tMultiMesh
  USE EL_QUADRATURE, ONLY: EL_INTEGRATE_FLUID_MOMENTUM, EL_INTEGRATE_UZ_LANES
  USE EL_CONFIG, ONLY: el_momentum_audit_freq, el_tavg_window, &
                       el_domain_type, el_cylinder_center, &
                       el_cylinder_radius, el_cylinder_axis, &
                       el_momentum_fix, el_meso_filter_bins, &
                       el_shear_rate
  USE def_FEAT, ONLY: NITNS
#ifdef HAVE_PE
  USE DEM_QUERY, ONLY: numLocalParticles, getAllParticles, &
                       getElContactCount, getElMaxPenetration, &
                       getElStepsize, getElLubricationVirial, &
                       getElLubricationPairs, getElLubricationImpulse
#endif

  IMPLICIT NONE

  LOGICAL :: el_momentum_reference_set = .FALSE.
  REAL*8 :: el_momentum_reference(3) = 0.0d0
  REAL*8 :: el_pair_momentum_before(3) = 0.0d0
  REAL*8 :: el_fluid_pair_before(3) = 0.0d0
  REAL*8 :: el_fluid_pair_mid(3) = 0.0d0
  LOGICAL :: el_fluid_pair_stage_set = .FALSE.
  ! Compensating acceleration (per unit mass, ConstForce convention) for
  ! the measured-leak momentum fix; applied by EL_APPLY_MOMENTUM_FIX via
  ! Grav_QuadSc next step. Recursion C_{n+1} = C_n - mismatch_n keeps the
  ! cumulative fluid-momentum drift bounded by a single step's leak.
  REAL*8 :: el_momfix_accel(3) = 0.0d0
  ! Mesoscale lane-mode filter state: global horizontal domain bounds,
  ! resolved once on first use (Allreduce over vertex coordinates).
  LOGICAL :: el_meso_bounds_set = .FALSE.
  REAL*8 :: el_meso_xmin = 0.0d0, el_meso_xlen = 0.0d0
  REAL*8 :: el_meso_ymin = 0.0d0, el_meso_ylen = 0.0d0
  INTEGER :: el_susp_tavg_count = 0
  REAL*8 :: el_susp_tavg_sig(9) = 0.0d0
  INTEGER :: el_mean_slip_tavg_count = 0
  REAL*8 :: el_mean_slip_tavg_up(3) = 0.0d0
  REAL*8 :: el_mean_slip_tavg_uf_intr(3) = 0.0d0
  REAL*8 :: el_mean_slip_tavg_uf_super(3) = 0.0d0
  REAL*8 :: el_mean_slip_tavg_slip_intr(3) = 0.0d0
  REAL*8 :: el_mean_slip_tavg_eps = 0.0d0

CONTAINS

  SUBROUTINE EL_WRITE_MOMENTUM_DIAGNOSTICS(val_u, val_v, val_w, &
                                           time, mfile, istep, mesh, ilev, &
                                           density, epsilon_f)

    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:)
    REAL*8, INTENT(IN) :: time
    INTEGER, INTENT(IN) :: mfile, istep
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: density, epsilon_f(:)

    REAL*8 :: local_particle(3), global_particle(3)
    REAL*8 :: local_velocity_sum(3), global_velocity_sum(3)
    REAL*8 :: mean_velocity(3), speed, ref_norm
    REAL*8 :: local_vol, global_vol, local_eps_vol, global_eps_vol
    REAL*8 :: uf_intr(3), uf_super(3), slip_intr(3), eps_mean
    REAL*8 :: tavg_up(3), tavg_uf_intr(3), tavg_uf_super(3)
    REAL*8 :: tavg_slip_intr(3), tavg_eps
    REAL*8 :: local_tgran, global_tgran, max_overlap_local, max_overlap
    REAL*8 :: susp_local(9), susp_virial(9), susp_sig(9), susp_tavg(9)
    REAL*8 :: susp_dt, susp_eta
    INTEGER :: lubpairs_local, lubpairs
    REAL*8 :: lubdp_local(3), lubdp(3)
    INTEGER :: tavg_start_step, ncontacts_local, ncontacts
    ! Element-integrated (consistent, no shared-DOF double count) fluid momentum,
    ! plain rho*u and void-fraction-weighted rho*eps_f*u.
    REAL*8 :: local_fei(3), global_fei(3), local_fei_eps(3), global_fei_eps(3)
    REAL*8 :: total_ei(3), drift_ei(3), drift_ei_norm, drift_ei_rel
    REAL*8 :: total_ei_eps(3), drift_ei_eps(3), drift_ei_eps_norm
    INTEGER :: local_count, global_count, ierr, nel_vol

    IF (MOD(ABS(istep), el_momentum_audit_freq).NE.0) RETURN

    CALL EL_LOCAL_PARTICLE_MOMENTUM(local_particle, local_velocity_sum, &
                                    local_count)

    CALL MPI_Allreduce(local_particle, global_particle, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_velocity_sum, global_velocity_sum, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_count, global_count, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_SUBS, ierr)

    ! Element-integrated fluid momentum (each cell once; no shared-DOF double
    ! count). global_fei = integral rho*u; global_fei_eps = integral rho*eps_f*u.
    CALL EL_INTEGRATE_FLUID_MOMENTUM(mesh, ilev, val_u, val_v, val_w, &
                                     epsilon_f, density, local_fei, local_fei_eps)
    CALL MPI_Allreduce(local_fei, global_fei, 3, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_fei_eps, global_fei_eps, 3, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    ! dvol carries the FEAT total-volume slot at index nel+1 (mesh_refine.f90
    ! SETARE convention: DVOL(NEL+1)=SUM) -- an unsliced SUM double-counts the
    ! domain volume and an unsliced epsilon_f*dvol product is nonconforming.
    ! Always slice to the first nel entries.
    nel_vol = mesh%level(ilev)%nel
    local_vol = SUM(mesh%level(ilev)%dvol(1:nel_vol))
    local_eps_vol = SUM(epsilon_f(1:nel_vol)*mesh%level(ilev)%dvol(1:nel_vol))
    CALL MPI_Allreduce(local_vol, global_vol, 1, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_eps_vol, global_eps_vol, 1, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)

    ! Element-integrated totals/drift, measured against the SAME t=0 reference
    ! (at t=0 the fluid is at rest, so the initial particle momentum is enough).
    total_ei = global_fei + global_particle
    IF (.NOT.el_momentum_reference_set) THEN
      el_momentum_reference = total_ei
      el_momentum_reference_set = .TRUE.
    END IF
    ref_norm = SQRT(SUM(el_momentum_reference**2))
    drift_ei = total_ei - el_momentum_reference
    drift_ei_norm = SQRT(SUM(drift_ei**2))
    drift_ei_rel = drift_ei_norm / MAX(ref_norm, 1.0d-30)
    total_ei_eps = global_fei_eps + global_particle
    drift_ei_eps = total_ei_eps - el_momentum_reference
    drift_ei_eps_norm = SQRT(SUM(drift_ei_eps**2))

    mean_velocity = 0.0d0
    IF (global_count.GT.0) THEN
      mean_velocity = global_velocity_sum / DBLE(global_count)
    END IF
    speed = SQRT(SUM(mean_velocity**2))

    uf_intr = 0.0d0
    uf_super = 0.0d0
    IF (global_eps_vol.GT.0.0d0) uf_intr = global_fei_eps / &
      (density*global_eps_vol)
    IF (global_vol.GT.0.0d0) uf_super = global_fei_eps / &
      (density*global_vol)
    slip_intr = mean_velocity - uf_intr
    eps_mean = 0.0d0
    IF (global_vol.GT.0.0d0) eps_mean = global_eps_vol / global_vol

    local_tgran = 0.0d0
#ifdef HAVE_PE
    IF (local_count.GT.0) THEN
      CALL EL_LOCAL_GRANULAR_TEMPERATURE(mean_velocity, local_tgran)
    END IF
#endif
    CALL MPI_Allreduce(local_tgran, global_tgran, 1, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    IF (global_count.GT.0) global_tgran = global_tgran / DBLE(global_count)


    ! Pairwise-lubrication impulse virial (Kroupa closure in PE): the MPI sum
    ! counts each pair exactly once (0.5 weight per locally-owned member on
    ! the PE side); particle-phase stress sigma = -virial/(dt*V). Zero when
    ! lubrication is disabled. The pair count double-counts cross-rank pairs
    ! (informational, like ncontacts).
    susp_local = 0.0d0
    lubpairs_local = 0
    susp_dt = 0.0d0
    lubdp_local = 0.0d0
#ifdef HAVE_PE
    CALL getElLubricationVirial(susp_local)
    CALL getElLubricationImpulse(lubdp_local)
    lubpairs_local = getElLubricationPairs()
    susp_dt = getElStepsize()
#endif
    CALL MPI_Allreduce(lubdp_local, lubdp, 3, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(susp_local, susp_virial, 9, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(lubpairs_local, lubpairs, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    susp_sig = 0.0d0
    IF (susp_dt.GT.0.0d0 .AND. global_vol.GT.0.0d0) &
      susp_sig = -susp_virial / (susp_dt*global_vol)

    tavg_start_step = 0
    IF (NITNS.GT.0) tavg_start_step = INT((1.0d0-el_tavg_window)*DBLE(NITNS))
    IF (ABS(istep).GE.tavg_start_step) THEN
      el_mean_slip_tavg_count = el_mean_slip_tavg_count + 1
      el_mean_slip_tavg_up = el_mean_slip_tavg_up + mean_velocity
      el_mean_slip_tavg_uf_intr = el_mean_slip_tavg_uf_intr + uf_intr
      el_mean_slip_tavg_uf_super = el_mean_slip_tavg_uf_super + uf_super
      el_mean_slip_tavg_slip_intr = el_mean_slip_tavg_slip_intr + slip_intr
      el_mean_slip_tavg_eps = el_mean_slip_tavg_eps + eps_mean
      el_susp_tavg_count = el_susp_tavg_count + 1
      el_susp_tavg_sig = el_susp_tavg_sig + susp_sig
    END IF
    susp_tavg = 0.0d0
    IF (el_susp_tavg_count.GT.0) susp_tavg = el_susp_tavg_sig / &
      DBLE(el_susp_tavg_count)
    IF (el_mean_slip_tavg_count.GT.0) THEN
      tavg_up = el_mean_slip_tavg_up / DBLE(el_mean_slip_tavg_count)
      tavg_uf_intr = el_mean_slip_tavg_uf_intr / DBLE(el_mean_slip_tavg_count)
      tavg_uf_super = el_mean_slip_tavg_uf_super / DBLE(el_mean_slip_tavg_count)
      tavg_slip_intr = el_mean_slip_tavg_slip_intr / DBLE(el_mean_slip_tavg_count)
      tavg_eps = el_mean_slip_tavg_eps / DBLE(el_mean_slip_tavg_count)
    ELSE
      tavg_up = 0.0d0
      tavg_uf_intr = 0.0d0
      tavg_uf_super = 0.0d0
      tavg_slip_intr = 0.0d0
      tavg_eps = 0.0d0
    END IF

#ifdef HAVE_PE
    ncontacts_local = getElContactCount()
    max_overlap_local = getElMaxPenetration()
#else
    ncontacts_local = 0
    max_overlap_local = 0.0d0
#endif
    ! Rank-summed contacts are a stationarity diagnostic; shadow-copy
    ! boundary contacts may be represented on more than one rank.
    CALL MPI_Allreduce(ncontacts_local, ncontacts, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(max_overlap_local, max_overlap, 1, MPI_DOUBLE_PRECISION, &
      MPI_MAX, MPI_COMM_SUBS, ierr)

    IF (myid.EQ.showid) THEN
      WRITE(*,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
        'EL_MOMENTUM_ELEMINT time= ', time, ' step= ', istep, ' fluid= ', &
        global_fei, ' total= ', total_ei, ' drift= ', drift_ei_norm, &
        ' drift_rel= ', drift_ei_rel, ' drift_eps= ', drift_ei_eps_norm
      WRITE(mfile,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
        'EL_MOMENTUM_ELEMINT time= ', time, ' step= ', istep, ' fluid= ', &
        global_fei, ' total= ', total_ei, ' drift= ', drift_ei_norm, &
        ' drift_rel= ', drift_ei_rel, ' drift_eps= ', drift_ei_eps_norm
      WRITE(*,'(A,ES14.6,A,I0,A,I0,A,3ES14.6,A,ES14.6)') &
        'EL_TERMINAL_VEL time= ', time, ' step= ', istep, ' count= ', &
        global_count, ' velocity= ', mean_velocity, ' speed= ', speed
      WRITE(mfile,'(A,ES14.6,A,I0,A,I0,A,3ES14.6,A,ES14.6)') &
        'EL_TERMINAL_VEL time= ', time, ' step= ', istep, ' count= ', &
        global_count, ' velocity= ', mean_velocity, ' speed= ', speed
      WRITE(*,'(A,ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,ES14.6)') &
        'EL_MEAN_SLIP t= ', time, ' up= ', mean_velocity, &
        ' uf_intr= ', uf_intr, ' uf_super= ', uf_super, &
        ' slip_intr= ', slip_intr, ' eps= ', eps_mean
      WRITE(mfile,'(A,ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,ES14.6)') &
        'EL_MEAN_SLIP t= ', time, ' up= ', mean_velocity, &
        ' uf_intr= ', uf_intr, ' uf_super= ', uf_super, &
        ' slip_intr= ', slip_intr, ' eps= ', eps_mean
      WRITE(*,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,ES14.6)') &
        'EL_MEAN_SLIP_TAVG t= ', time, ' samples= ', el_mean_slip_tavg_count, &
        ' up= ', tavg_up, ' uf_intr= ', tavg_uf_intr, &
        ' uf_super= ', tavg_uf_super, ' slip_intr= ', tavg_slip_intr, &
        ' eps= ', tavg_eps
      WRITE(mfile,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,ES14.6)') &
        'EL_MEAN_SLIP_TAVG t= ', time, ' samples= ', el_mean_slip_tavg_count, &
        ' up= ', tavg_up, ' uf_intr= ', tavg_uf_intr, &
        ' uf_super= ', tavg_uf_super, ' slip_intr= ', tavg_slip_intr, &
        ' eps= ', tavg_eps
      WRITE(*,'(A,ES14.6,A,I0,A,ES14.6,A,ES14.6)') &
        'EL_CONTACT_STATS t= ', time, ' ncontacts= ', ncontacts, &
        ' max_overlap= ', max_overlap, ' Tgran= ', global_tgran
      WRITE(mfile,'(A,ES14.6,A,I0,A,ES14.6,A,ES14.6)') &
        'EL_CONTACT_STATS t= ', time, ' ncontacts= ', ncontacts, &
        ' max_overlap= ', max_overlap, ' Tgran= ', global_tgran
      ! sig_xz row-major component (1-based index 3); etaL = sig_xz/G for the
      ! frozen linear-shear box u_x(z) with G = el_shear_rate.
      susp_eta = 0.0d0
      IF (ABS(el_shear_rate).GT.0.0d0) susp_eta = susp_tavg(3)/el_shear_rate
      WRITE(*,'(A,ES14.6,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,9ES14.6)') &
        'EL_SUSP_STRESS t= ', time, ' pairs= ', lubpairs, &
        ' sig_xz= ', susp_sig(3), ' sig_xz_tavg= ', susp_tavg(3), &
        ' etaL_tavg= ', susp_eta, ' sig= ', susp_sig
      WRITE(mfile,'(A,ES14.6,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,9ES14.6)') &
        'EL_SUSP_STRESS t= ', time, ' pairs= ', lubpairs, &
        ' sig_xz= ', susp_sig(3), ' sig_xz_tavg= ', susp_tavg(3), &
        ' etaL_tavg= ', susp_eta, ' sig= ', susp_sig
      WRITE(*,'(A,ES14.6,A,3ES14.6)') &
        'EL_LUB_IMPULSE t= ', time, ' net_dp= ', lubdp
      WRITE(mfile,'(A,ES14.6,A,3ES14.6)') &
        'EL_LUB_IMPULSE t= ', time, ' net_dp= ', lubdp
    END IF

  END SUBROUTINE EL_WRITE_MOMENTUM_DIAGNOSTICS

  SUBROUTINE EL_FLUID_PAIR_BEGIN(val_u, val_v, val_w, mesh, ilev, density, &
                                 epsilon_f)

    ! Snapshot the global element-integrated fluid momentum (plain rho*u)
    ! immediately before the fluid solve. Paired with EL_FLUID_PAIR_END.
    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: density, epsilon_f(:)

    REAL*8 :: local_f(3), local_feps(3)
    INTEGER :: ierr

    CALL EL_INTEGRATE_FLUID_MOMENTUM(mesh, ilev, val_u, val_v, val_w, &
                                     epsilon_f, density, local_f, local_feps)
    CALL MPI_Allreduce(local_f, el_fluid_pair_before, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)

  END SUBROUTINE EL_FLUID_PAIR_BEGIN

  SUBROUTINE EL_FLUID_PAIR_STAGE(val_u, val_v, val_w, mesh, ilev, density, &
                                 epsilon_f)

    ! Snapshot the global element-integrated fluid momentum AFTER the
    ! nonlinear momentum (Burgers) solve but BEFORE the projection
    ! (pressure Poisson + velocity correction). Lets EL_FLUID_PAIR_END
    ! split its mismatch into a momentum-solve part (convection/
    ! stabilization assembly -- all sources enter this stage's rhs) and a
    ! projection part, which must be momentum-neutral in a periodic box
    ! (the discrete pressure gradient annihilates constant test fields).
    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: density, epsilon_f(:)

    REAL*8 :: local_f(3), local_feps(3)
    INTEGER :: ierr

    CALL EL_INTEGRATE_FLUID_MOMENTUM(mesh, ilev, val_u, val_v, val_w, &
                                     epsilon_f, density, local_f, local_feps)
    CALL MPI_Allreduce(local_f, el_fluid_pair_mid, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    el_fluid_pair_stage_set = .TRUE.

  END SUBROUTINE EL_FLUID_PAIR_STAGE

  SUBROUTINE EL_MESO_FILTER_APPLY(val_w, mesh, ilev, coor, ndof)

    ! Subtract from u_z its zero-mean x-slab and y-slab averages (the
    ! gravest lane modes of the periodic-box settling instability; the
    ! wall-suppression surrogate, see el_config.f90). The subtracted
    ! nodal field depends on (x,y) only, so on the axis-aligned periodic
    ! box it is discretely divergence-free, and its global mean is
    ! removed so it carries (up to slab-vs-Q2 representation residue)
    ! zero net momentum. Caller pins ILEV to the coupling level.
    REAL*8, INTENT(INOUT) :: val_w(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, ndof
    REAL*8, INTENT(IN) :: coor(:,:)

    INTEGER :: nbins, i, ib, ierr, ig
    REAL*8 :: gmean, xw, yw
    REAL*8 :: bounds_local(4), bounds(4)
    REAL*8, ALLOCATABLE :: lanes_local(:), lanes(:)
    REAL*8, ALLOCATABLE :: fx(:), fy(:)

    nbins = MAX(2, el_meso_filter_bins)

    IF (.NOT.el_meso_bounds_set) THEN
      bounds_local(1) = MINVAL(mesh%level(ilev)%dcorvg(1,1:mesh%level(ilev)%nvt))
      bounds_local(2) = -MAXVAL(mesh%level(ilev)%dcorvg(1,1:mesh%level(ilev)%nvt))
      bounds_local(3) = MINVAL(mesh%level(ilev)%dcorvg(2,1:mesh%level(ilev)%nvt))
      bounds_local(4) = -MAXVAL(mesh%level(ilev)%dcorvg(2,1:mesh%level(ilev)%nvt))
      CALL MPI_Allreduce(bounds_local, bounds, 4, MPI_DOUBLE_PRECISION, &
        MPI_MIN, MPI_COMM_SUBS, ierr)
      el_meso_xmin = bounds(1)
      el_meso_xlen = -bounds(2) - bounds(1)
      el_meso_ymin = bounds(3)
      el_meso_ylen = -bounds(4) - bounds(3)
      el_meso_bounds_set = .TRUE.
    END IF
    IF (el_meso_xlen.LE.0.0d0 .OR. el_meso_ylen.LE.0.0d0) RETURN

    ALLOCATE(lanes_local(4*nbins), lanes(4*nbins), fx(nbins), fy(nbins))
    CALL EL_INTEGRATE_UZ_LANES(mesh, ilev, val_w, nbins, &
      el_meso_xmin, el_meso_xlen, el_meso_ymin, el_meso_ylen, &
      lanes_local(1:nbins), lanes_local(nbins+1:2*nbins), &
      lanes_local(2*nbins+1:3*nbins), lanes_local(3*nbins+1:4*nbins))
    CALL MPI_Allreduce(lanes_local, lanes, 4*nbins, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)

    IF (SUM(lanes(nbins+1:2*nbins)).LE.0.0d0) THEN
      DEALLOCATE(lanes_local, lanes, fx, fy)
      RETURN
    END IF
    gmean = SUM(lanes(1:nbins))/SUM(lanes(nbins+1:2*nbins))
    DO i=1,nbins
      fx(i) = 0.0d0
      IF (lanes(nbins+i).GT.0.0d0) fx(i) = lanes(i)/lanes(nbins+i) - gmean
      fy(i) = 0.0d0
      IF (lanes(3*nbins+i).GT.0.0d0) fy(i) = lanes(2*nbins+i)/lanes(3*nbins+i) - gmean
    END DO

    ! Periodic wrap in the dof binning so matched dofs on the two copies
    ! of a periodic face receive identical corrections.
    DO ig=1,ndof
      xw = MODULO(coor(1,ig)-el_meso_xmin, el_meso_xlen)
      yw = MODULO(coor(2,ig)-el_meso_ymin, el_meso_ylen)
      ib = MIN(nbins, 1 + INT(xw/el_meso_xlen*DBLE(nbins)))
      val_w(ig) = val_w(ig) - fx(ib)
      ib = MIN(nbins, 1 + INT(yw/el_meso_ylen*DBLE(nbins)))
      val_w(ig) = val_w(ig) - fy(ib)
    END DO

    DEALLOCATE(lanes_local, lanes, fx, fy)

  END SUBROUTINE EL_MESO_FILTER_APPLY

  SUBROUTINE EL_FLUID_PAIR_END(val_u, val_v, val_w, mesh, ilev, density, &
                               epsilon_f, source, ext_force, dt, time, &
                               mfile, istep)

    ! Fluid-side momentum audit: the measured change of the element-
    ! integrated fluid momentum over the fluid solve vs the impulse the
    ! step SHOULD add: dt*(deposited EL feedback + external body force).
    ! Viscous/pressure/convective terms are internal and must be momentum-
    ! neutral in a fully periodic domain, so a persistent mismatch measures
    ! the fluid-internal discrete non-conservation (e.g. the convective
    ! form) plus any source-application defect.
    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: density, epsilon_f(:), source(:,:), ext_force(3)
    REAL*8, INTENT(IN) :: dt, time
    INTEGER, INTENT(IN) :: mfile, istep

    REAL*8 :: local_f(3), local_feps(3), p_after(3), dp(3)
    REAL*8 :: local_src(3), global_src(3), expected(3), mismatch(3)
    REAL*8 :: local_vol, global_vol
    INTEGER :: ierr, nel_vol

    CALL EL_INTEGRATE_FLUID_MOMENTUM(mesh, ilev, val_u, val_v, val_w, &
                                     epsilon_f, density, local_f, local_feps)
    CALL MPI_Allreduce(local_f, p_after, 3, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    ! With the momentum fix the mismatch feeds the compensator and must be
    ! computed EVERY step; without it, non-audit steps stop here.
    IF (.NOT.el_momentum_fix .AND. &
        MOD(ABS(istep), el_momentum_audit_freq).NE.0) RETURN

    ! ext_force is the body force per unit volume; integrate over the domain
    ! (dvol carries the FEAT total-volume slot at nel+1 -- slice!).
    nel_vol = mesh%level(ilev)%nel
    local_vol = SUM(mesh%level(ilev)%dvol(1:nel_vol))
    CALL MPI_Allreduce(local_vol, global_vol, 1, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)

    ! The feedback source is stored in distributed/additive form; the rank
    ! sum of DOF sums is the total deposited force.
    local_src(1) = SUM(source(1,:))
    local_src(2) = SUM(source(2,:))
    local_src(3) = SUM(source(3,:))
    CALL MPI_Allreduce(local_src, global_src, 3, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)

    dp = p_after - el_fluid_pair_before
    ! ext_force excludes the compensating force, so with the fix active
    ! mismatch_n = leak_n + C_n and the recursion below yields
    ! C_{n+1} = -leak_n (deadbeat, one-step lag).
    expected = dt*(global_src + ext_force*global_vol)
    mismatch = dp - expected
    IF (el_momentum_fix .AND. dt.GT.0.0d0 .AND. global_vol.GT.0.0d0) THEN
      el_momfix_accel = el_momfix_accel - &
        mismatch/(dt*density*global_vol)
    END IF
    IF (MOD(ABS(istep), el_momentum_audit_freq).NE.0) THEN
      el_fluid_pair_stage_set = .FALSE.
      RETURN
    END IF
    IF (myid.EQ.showid) THEN
      WRITE(*,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6)') &
        'EL_FLUID_PAIR t= ', time, ' step= ', istep, ' dp= ', dp, &
        ' expected= ', expected, ' mismatch= ', mismatch
      WRITE(mfile,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6)') &
        'EL_FLUID_PAIR t= ', time, ' step= ', istep, ' dp= ', dp, &
        ' expected= ', expected, ' mismatch= ', mismatch
      IF (el_fluid_pair_stage_set) THEN
        ! Stage split: all sources belong to the momentum stage, the
        ! projection stage must be momentum-neutral, so
        ! mismatch == mis_mom + mis_prj.
        WRITE(*,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6)') &
          'EL_FLUID_STAGE t= ', time, ' step= ', istep, &
          ' mis_mom= ', el_fluid_pair_mid - el_fluid_pair_before - &
          expected, ' mis_prj= ', p_after - el_fluid_pair_mid
        WRITE(mfile,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6)') &
          'EL_FLUID_STAGE t= ', time, ' step= ', istep, &
          ' mis_mom= ', el_fluid_pair_mid - el_fluid_pair_before - &
          expected, ' mis_prj= ', p_after - el_fluid_pair_mid
      END IF
    END IF
    el_fluid_pair_stage_set = .FALSE.

  END SUBROUTINE EL_FLUID_PAIR_END

  SUBROUTINE EL_NEWTON_PAIR_BEGIN()

    ! Snapshot the global particle momentum immediately before
    ! EL_ADVANCE_PARTICLES. Paired with EL_NEWTON_PAIR_END.
    REAL*8 :: local_p(3), vsum(3)
    INTEGER :: cnt, ierr

    CALL EL_LOCAL_PARTICLE_MOMENTUM(local_p, vsum, cnt)
    CALL MPI_Allreduce(local_p, el_pair_momentum_before, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)

  END SUBROUTINE EL_NEWTON_PAIR_BEGIN

  SUBROUTINE EL_NEWTON_PAIR_END(expected_impulse, time, mfile, istep)

    ! Third-law audit: the measured particle-momentum change over the PE
    ! advance vs the impulse the CFD side charged to the fluid (mirrored
    ! drag + lift + pressure + grav_buoy). Contact impulses are internal
    ! and cancel in the global sum, so mismatch = Newton-pair violation
    ! between PE integration and the deposited feedback.
    REAL*8, INTENT(IN) :: expected_impulse(3), time
    INTEGER, INTENT(IN) :: mfile, istep

    REAL*8 :: local_p(3), p_after(3), vsum(3), dp(3), mismatch(3)
    INTEGER :: cnt, ierr

    CALL EL_LOCAL_PARTICLE_MOMENTUM(local_p, vsum, cnt)
    CALL MPI_Allreduce(local_p, p_after, 3, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    IF (MOD(ABS(istep), el_momentum_audit_freq).NE.0) RETURN

    dp = p_after - el_pair_momentum_before
    mismatch = dp - expected_impulse
    IF (myid.EQ.showid) THEN
      WRITE(*,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6)') &
        'EL_NEWTON_PAIR t= ', time, ' step= ', istep, ' dp= ', dp, &
        ' expected= ', expected_impulse, ' mismatch= ', mismatch
      WRITE(mfile,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6)') &
        'EL_NEWTON_PAIR t= ', time, ' step= ', istep, ' dp= ', dp, &
        ' expected= ', expected_impulse, ' mismatch= ', mismatch
    END IF

  END SUBROUTINE EL_NEWTON_PAIR_END

  SUBROUTINE EL_WRITE_SS_RADIUS(time, mfile, istep)

    ! Radial-position statistics of all particles about the configured
    ! cylinder axis, normalized by the cylinder radius (Segre-Silberberg
    ! migration monitor). Unlike EL_WRITE_MOMENTUM_DIAGNOSTICS this is also
    ! called in prescribed-field (frozen) runs, where the momentum audit is
    ! meaningless but radial trajectories are the payload.
    REAL*8, INTENT(IN) :: time
    INTEGER, INTENT(IN) :: mfile, istep

    REAL*8 :: local_rsum, global_rsum, local_rmin, global_rmin
    REAL*8 :: local_rmax, global_rmax, rmean, r
    INTEGER :: local_count, global_count, ierr
#ifdef HAVE_PE
    TYPE(tParticleData), ALLOCATABLE :: particles(:)
    REAL*8 :: dx(3)
    INTEGER :: i
#endif

    IF (TRIM(el_domain_type).NE.'cylinder') RETURN
    IF (el_cylinder_radius.LE.0.0d0) RETURN
    IF (MOD(ABS(istep), el_momentum_audit_freq).NE.0) RETURN

    local_rsum = 0.0d0
    local_rmin = HUGE(1.0d0)
    local_rmax = -HUGE(1.0d0)
    local_count = 0

#ifdef HAVE_PE
    local_count = numLocalParticles()
    IF (local_count.GT.0) THEN
      ALLOCATE(particles(local_count))
      CALL getAllParticles(particles)
      DO i = 1, local_count
        dx = particles(i)%position - el_cylinder_center
        SELECT CASE (TRIM(el_cylinder_axis))
        CASE ('x')
          r = SQRT(dx(2)**2 + dx(3)**2)
        CASE ('y')
          r = SQRT(dx(1)**2 + dx(3)**2)
        CASE DEFAULT
          r = SQRT(dx(1)**2 + dx(2)**2)
        END SELECT
        r = r / el_cylinder_radius
        local_rsum = local_rsum + r
        local_rmin = MIN(local_rmin, r)
        local_rmax = MAX(local_rmax, r)
      END DO
      DEALLOCATE(particles)
    END IF
#endif

    CALL MPI_Allreduce(local_rsum, global_rsum, 1, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_rmin, global_rmin, 1, MPI_DOUBLE_PRECISION, &
      MPI_MIN, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_rmax, global_rmax, 1, MPI_DOUBLE_PRECISION, &
      MPI_MAX, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_count, global_count, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_SUBS, ierr)

    IF (global_count.LE.0) RETURN
    rmean = global_rsum / DBLE(global_count)

    IF (myid.EQ.showid) THEN
      WRITE(*,'(A,ES14.6,A,I0,A,I0,A,ES14.6,A,ES14.6,A,ES14.6)') &
        'EL_SS_RADIUS t= ', time, ' step= ', istep, ' count= ', &
        global_count, ' rmean= ', rmean, ' rmin= ', global_rmin, &
        ' rmax= ', global_rmax
      WRITE(mfile,'(A,ES14.6,A,I0,A,I0,A,ES14.6,A,ES14.6,A,ES14.6)') &
        'EL_SS_RADIUS t= ', time, ' step= ', istep, ' count= ', &
        global_count, ' rmean= ', rmean, ' rmin= ', global_rmin, &
        ' rmax= ', global_rmax
    END IF

  END SUBROUTINE EL_WRITE_SS_RADIUS

  SUBROUTINE EL_LOCAL_PARTICLE_MOMENTUM(momentum, velocity_sum, count)

    REAL*8, INTENT(OUT) :: momentum(3), velocity_sum(3)
    INTEGER, INTENT(OUT) :: count

#ifdef HAVE_PE
    TYPE(tParticleData), ALLOCATABLE :: particles(:)
    REAL*8 :: mass
    INTEGER :: i
#endif

    momentum = 0.0d0
    velocity_sum = 0.0d0
    count = 0

#ifdef HAVE_PE
    count = numLocalParticles()
    IF (count.LE.0) RETURN

    ALLOCATE(particles(count))
    CALL getAllParticles(particles)

    DO i = 1, count
      mass = 4.0d0*ACOS(-1.0d0)*particles(i)%radius**3 * &
        particles(i)%density / 3.0d0
      momentum = momentum + mass*particles(i)%velocity
      velocity_sum = velocity_sum + particles(i)%velocity
    END DO

    DEALLOCATE(particles)
#endif

  END SUBROUTINE EL_LOCAL_PARTICLE_MOMENTUM

  SUBROUTINE EL_LOCAL_GRANULAR_TEMPERATURE(mean_velocity, tgran_sum)

    REAL*8, INTENT(IN) :: mean_velocity(3)
    REAL*8, INTENT(OUT) :: tgran_sum

#ifdef HAVE_PE
    TYPE(tParticleData), ALLOCATABLE :: particles(:)
    REAL*8 :: vel_fluct(3)
    INTEGER :: i, count
#endif

    tgran_sum = 0.0d0
#ifdef HAVE_PE
    count = numLocalParticles()
    IF (count.LE.0) RETURN
    ALLOCATE(particles(count))
    CALL getAllParticles(particles)
    DO i=1,count
      vel_fluct = particles(i)%velocity - mean_velocity
      tgran_sum = tgran_sum + SUM(vel_fluct**2)/3.0d0
    END DO
    DEALLOCATE(particles)
#endif

  END SUBROUTINE EL_LOCAL_GRANULAR_TEMPERATURE

  SUBROUTINE EL_CAPTURE_MOMENTUM_REFERENCE(val_u, val_v, val_w, mesh, ilev, &
                                           density, epsilon_f)
    ! Capture the same element-integrated total momentum used by the permanent
    ! audit. Call this before the first coupled update so drift is measured from
    ! the true initial condition rather than from the post-step-1 state.
    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: density, epsilon_f(:)

    REAL*8 :: local_fluid(3), global_fluid(3), local_fluid_eps(3)
    REAL*8 :: local_particle(3), global_particle(3)
    REAL*8 :: local_velocity_sum(3)
    INTEGER :: local_count, ierr

    CALL EL_INTEGRATE_FLUID_MOMENTUM(mesh, ilev, val_u, val_v, val_w, &
                                     epsilon_f, density, local_fluid, &
                                     local_fluid_eps)

    CALL EL_LOCAL_PARTICLE_MOMENTUM(local_particle, local_velocity_sum, &
                                    local_count)

    CALL MPI_Allreduce(local_fluid, global_fluid, 3, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_particle, global_particle, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)

    el_momentum_reference = global_fluid + global_particle
    el_momentum_reference_set = .TRUE.

  END SUBROUTINE EL_CAPTURE_MOMENTUM_REFERENCE

  SUBROUTINE EL_RESET_MOMENTUM_REFERENCE()

    el_momentum_reference_set = .FALSE.
    el_momentum_reference = 0.0d0
    el_mean_slip_tavg_count = 0
    el_mean_slip_tavg_up = 0.0d0
    el_mean_slip_tavg_uf_intr = 0.0d0
    el_mean_slip_tavg_uf_super = 0.0d0
    el_mean_slip_tavg_slip_intr = 0.0d0
    el_mean_slip_tavg_eps = 0.0d0

  END SUBROUTINE EL_RESET_MOMENTUM_REFERENCE

END MODULE EL_DIAGNOSTICS
