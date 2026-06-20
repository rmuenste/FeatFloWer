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
    IF (.NOT.el_momentum_reference_set) THEN
      el_momentum_reference = total
      el_momentum_reference_set = .TRUE.
    END IF
    drift = total - el_momentum_reference
    drift_norm = SQRT(SUM(drift**2))

    mean_velocity = 0.0d0
    IF (global_count.GT.0) THEN
      mean_velocity = global_velocity_sum / DBLE(global_count)
    END IF
    speed = SQRT(SUM(mean_velocity**2))

    IF (myid.EQ.showid) THEN
      WRITE(*,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,ES14.6)') &
        'EL_MOMENTUM time= ', time, ' step= ', istep, ' fluid= ', &
        global_fluid, ' particle= ', global_particle, ' total= ', total, &
        ' drift= ', drift_norm
      WRITE(mfile,'(A,ES14.6,A,I0,A,3ES14.6,A,3ES14.6,A,3ES14.6,A,ES14.6)') &
        'EL_MOMENTUM time= ', time, ' step= ', istep, ' fluid= ', &
        global_fluid, ' particle= ', global_particle, ' total= ', total, &
        ' drift= ', drift_norm
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

  SUBROUTINE EL_RESET_MOMENTUM_REFERENCE()

    el_momentum_reference_set = .FALSE.
    el_momentum_reference = 0.0d0

  END SUBROUTINE EL_RESET_MOMENTUM_REFERENCE

END MODULE EL_DIAGNOSTICS
