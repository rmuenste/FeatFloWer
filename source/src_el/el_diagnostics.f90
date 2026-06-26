MODULE EL_DIAGNOSTICS

  USE MPI
  USE PP3D_MPI, ONLY: myid, showid, MPI_COMM_SUBS
  USE DEM_QUERY, ONLY: tParticleData
#ifdef HAVE_PE
  USE DEM_QUERY, ONLY: numLocalParticles, getAllParticles
#endif

  IMPLICIT NONE

  LOGICAL :: el_momentum_reference_set = .FALSE.
  REAL*8 :: el_momentum_reference(3) = 0.0d0

CONTAINS

  SUBROUTINE EL_WRITE_MOMENTUM_DIAGNOSTICS(val_u, val_v, val_w, &
                                           fluid_mass, time, mfile, istep)

    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:), fluid_mass(:)
    REAL*8, INTENT(IN) :: time
    INTEGER, INTENT(IN) :: mfile, istep

    REAL*8 :: local_fluid(3), global_fluid(3)
    REAL*8 :: local_particle(3), global_particle(3)
    REAL*8 :: local_velocity_sum(3), global_velocity_sum(3)
    REAL*8 :: total(3), mean_velocity(3), drift(3), speed, drift_norm
    REAL*8 :: ref_norm, drift_rel
    INTEGER :: local_count, global_count, ndof, ierr

    ndof = MIN(SIZE(fluid_mass), SIZE(val_u), SIZE(val_v), SIZE(val_w))
    IF (ndof.GT.0) THEN
      local_fluid(1) = SUM(fluid_mass(1:ndof)*val_u(1:ndof))
      local_fluid(2) = SUM(fluid_mass(1:ndof)*val_v(1:ndof))
      local_fluid(3) = SUM(fluid_mass(1:ndof)*val_w(1:ndof))
    ELSE
      local_fluid = 0.0d0
    END IF

    CALL EL_LOCAL_PARTICLE_MOMENTUM(local_particle, local_velocity_sum, &
                                    local_count)

    CALL MPI_Allreduce(local_fluid, global_fluid, 3, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_particle, global_particle, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_velocity_sum, global_velocity_sum, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_count, global_count, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_SUBS, ierr)

    total = global_fluid + global_particle
    ! Reference is normally captured before the first coupled update (see
    ! EL_CAPTURE_MOMENTUM_REFERENCE, called from the time loop). Fall back to the
    ! first diagnostic call if it was not pre-set.
    IF (.NOT.el_momentum_reference_set) THEN
      el_momentum_reference = total
      el_momentum_reference_set = .TRUE.
    END IF
    drift = total - el_momentum_reference
    drift_norm = SQRT(SUM(drift**2))
    ref_norm = SQRT(SUM(el_momentum_reference**2))
    drift_rel = drift_norm / MAX(ref_norm, 1.0d-30)

    mean_velocity = 0.0d0
    IF (global_count.GT.0) THEN
      mean_velocity = global_velocity_sum / DBLE(global_count)
    END IF
    speed = SQRT(SUM(mean_velocity**2))

    IF (myid.EQ.showid) THEN
      WRITE(*,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,ES14.6,A,ES14.6)') &
        'EL_MOMENTUM time= ', time, ' step= ', istep, ' fluid= ', &
        global_fluid, ' particle= ', global_particle, ' total= ', total, &
        ' drift= ', drift_norm, ' drift_rel= ', drift_rel
      WRITE(mfile,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,ES14.6,A,ES14.6)') &
        'EL_MOMENTUM time= ', time, ' step= ', istep, ' fluid= ', &
        global_fluid, ' particle= ', global_particle, ' total= ', total, &
        ' drift= ', drift_norm, ' drift_rel= ', drift_rel
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

  SUBROUTINE EL_CAPTURE_MOMENTUM_REFERENCE(val_u, val_v, val_w, fluid_mass)
    ! Capture the total-momentum reference (P_fluid + P_particle) at the CURRENT
    ! state. Call this BEFORE the first coupled update so drift is measured from
    ! the true initial condition rather than from the post-step-1 state (which
    ! would hide the first update's imbalance). The fluid arrays are optional: at
    ! the very first step the lumped mass matrix is not assembled yet, but the
    ! initial fluid is at rest so its momentum is zero and can be omitted.
    REAL*8, INTENT(IN), OPTIONAL :: val_u(:), val_v(:), val_w(:), fluid_mass(:)

    REAL*8 :: local_fluid(3), global_fluid(3)
    REAL*8 :: local_particle(3), global_particle(3)
    REAL*8 :: local_velocity_sum(3)
    INTEGER :: local_count, ndof, ierr

    local_fluid = 0.0d0
    IF (PRESENT(val_u).AND.PRESENT(val_v).AND.PRESENT(val_w).AND. &
        PRESENT(fluid_mass)) THEN
      ndof = MIN(SIZE(fluid_mass), SIZE(val_u), SIZE(val_v), SIZE(val_w))
      IF (ndof.GT.0) THEN
        local_fluid(1) = SUM(fluid_mass(1:ndof)*val_u(1:ndof))
        local_fluid(2) = SUM(fluid_mass(1:ndof)*val_v(1:ndof))
        local_fluid(3) = SUM(fluid_mass(1:ndof)*val_w(1:ndof))
      END IF
    END IF

    CALL EL_LOCAL_PARTICLE_MOMENTUM(local_particle, local_velocity_sum, &
                                    local_count)

    CALL MPI_Allreduce(local_fluid, global_fluid, 3, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(local_particle, global_particle, 3, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)

    el_momentum_reference = global_fluid + global_particle
    el_momentum_reference_set = .TRUE.

  END SUBROUTINE EL_CAPTURE_MOMENTUM_REFERENCE

  SUBROUTINE EL_DEBUG_GZSUM(label, field, n, istep, mfile, mass)
    ! Debug probe for the staged fluid z-momentum budget: prints the global
    ! z-component sum  Sum[ (mass*) field ]  with a label, so we can track where
    ! dt*fed_force stops matching the change in fluid z-momentum across the solve
    ! stages (feedback apply -> BC -> velocity solve -> pressure projection).
    CHARACTER(*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: field(:)
    INTEGER, INTENT(IN) :: n, istep, mfile
    REAL*8, INTENT(IN), OPTIONAL :: mass(:)

    REAL*8 :: loc, glob
    INTEGER :: m, ierr

    m = MIN(n, SIZE(field))
    IF (PRESENT(mass)) THEN
      m = MIN(m, SIZE(mass))
      loc = SUM(mass(1:m)*field(1:m))
    ELSE
      loc = SUM(field(1:m))
    END IF
    CALL MPI_Allreduce(loc, glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    IF (myid.EQ.showid) THEN
      WRITE(*,'(A,A,A,I0,A,ES20.12)') 'EL_BUDGET ', label, ' step= ', istep, &
        ' z= ', glob
      WRITE(mfile,'(A,A,A,I0,A,ES20.12)') 'EL_BUDGET ', label, ' step= ', istep, &
        ' z= ', glob
    END IF
  END SUBROUTINE EL_DEBUG_GZSUM

  SUBROUTINE EL_RESET_MOMENTUM_REFERENCE()

    el_momentum_reference_set = .FALSE.
    el_momentum_reference = 0.0d0

  END SUBROUTINE EL_RESET_MOMENTUM_REFERENCE

END MODULE EL_DIAGNOSTICS
