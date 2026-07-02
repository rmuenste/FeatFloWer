MODULE EL_DIAGNOSTICS

  USE MPI
  USE PP3D_MPI, ONLY: myid, showid, MPI_COMM_SUBS
  USE DEM_QUERY, ONLY: tParticleData
  USE TYPES, ONLY: tMultiMesh
  USE EL_QUADRATURE, ONLY: EL_INTEGRATE_FLUID_MOMENTUM
  USE EL_CONFIG, ONLY: el_momentum_audit_freq
#ifdef HAVE_PE
  USE DEM_QUERY, ONLY: numLocalParticles, getAllParticles
#endif

  IMPLICIT NONE

  LOGICAL :: el_momentum_reference_set = .FALSE.
  REAL*8 :: el_momentum_reference(3) = 0.0d0

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
    ! Element-integrated (consistent, no shared-DOF double count) fluid momentum,
    ! plain rho*u and void-fraction-weighted rho*eps_f*u.
    REAL*8 :: local_fei(3), global_fei(3), local_fei_eps(3), global_fei_eps(3)
    REAL*8 :: total_ei(3), drift_ei(3), drift_ei_norm, drift_ei_rel
    REAL*8 :: total_ei_eps(3), drift_ei_eps(3), drift_ei_eps_norm
    INTEGER :: local_count, global_count, ierr

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
    END IF

  END SUBROUTINE EL_WRITE_MOMENTUM_DIAGNOSTICS

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

  END SUBROUTINE EL_RESET_MOMENTUM_REFERENCE

END MODULE EL_DIAGNOSTICS
