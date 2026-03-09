MODULE timestep_control

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: SetSimulationTimeStep

CONTAINS

  SUBROUTINE SetSimulationTimeStep(dTimeStep)
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: dTimeStep

    DOUBLE PRECISION :: TSTEP, THETA, THSTEP, TIMENS, EPSNS
    INTEGER :: NITNS, ITNS
    COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS

#ifdef HAVE_PE
    INTERFACE
      SUBROUTINE set_pe_timestep(dTime) BIND(C, name="set_pe_timestep_")
        USE, INTRINSIC :: iso_c_binding, ONLY: c_double
        REAL(c_double), INTENT(IN) :: dTime
      END SUBROUTINE set_pe_timestep
    END INTERFACE
#endif

    TSTEP = dTimeStep
    THSTEP = TSTEP * THETA

#ifdef HAVE_PE
    CALL set_pe_timestep(TSTEP)
#endif
  END SUBROUTINE SetSimulationTimeStep

END MODULE timestep_control
