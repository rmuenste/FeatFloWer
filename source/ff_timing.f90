MODULE ff_timing
! Machine-readable per-call solver timing records for the solver baseline
! harness (Phase 0). Call sites are compiled only with -DFF_SOLVER_TIMING
! (CMake option ENABLE_SOLVER_TIMING).
!
! Records use the same field layout as the HYPRE_TIMING records emitted by
! HypreSolver.f90 so that one summarizer parses both:
!
!   FF_TIMING solver=<label> call=<n> setup_s=<s> solve_s=<s> total_s=<s> &
!             iterations=<n> rel_res=<r>
!
! ff_report_timing is collective over MPI_COMM_WORLD: every rank must call
! it. Times and iteration counts are max-reduced and printed by the master.
USE PP3D_MPI
IMPLICIT NONE

INTEGER :: ff_crs_u_calls = 0
INTEGER :: ff_crs_p_calls = 0
INTEGER :: ff_mg_velocity_calls = 0
INTEGER :: ff_mg_pressure_calls = 0

CONTAINS

SUBROUTINE ff_report_timing(label, call_index, setup_time, solve_time, &
                            total_time, num_iterations, rel_res)
CHARACTER(LEN=*), INTENT(IN) :: label
INTEGER, INTENT(IN) :: call_index, num_iterations
REAL*8, INTENT(IN) :: setup_time, solve_time, total_time, rel_res
#ifdef FF_SOLVER_TIMING
INTEGER :: ierr
REAL*8 :: local_times(3), min_times(3), local_stat(2), max_stat(2)

! Times are MIN-reduced: ranks enter an instrumented routine at different
! moments (the idle master arrives first and waits at the next collective),
! so the smallest per-rank duration is the straggler's pass through the
! solver and excludes that wait.
local_times = (/ setup_time, solve_time, total_time /)
CALL MPI_Reduce(local_times, min_times, 3, MPI_DOUBLE_PRECISION, MPI_MIN, &
                0, MPI_COMM_WORLD, ierr)
local_stat = (/ rel_res, DBLE(num_iterations) /)
CALL MPI_Reduce(local_stat, max_stat, 2, MPI_DOUBLE_PRECISION, MPI_MAX, &
                0, MPI_COMM_WORLD, ierr)

IF (myid.eq.0) THEN
  WRITE(*,'(A,A,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,I0,A,ES14.6)') &
    'FF_TIMING solver=', TRIM(label), ' call=', call_index, &
    ' setup_s=', min_times(1), ' solve_s=', min_times(2), &
    ' total_s=', min_times(3), ' iterations=', NINT(max_stat(2)), &
    ' rel_res=', max_stat(1)
END IF
#endif
END SUBROUTINE ff_report_timing

END MODULE ff_timing
