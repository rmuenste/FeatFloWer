PROGRAM Q2P1_EL_PIPEFLOW

  INCLUDE 'defs_include.h'

  USE SOLUTION_IO, ONLY: postprocessing_app
  USE APP_INITIALIZATION, ONLY: init_q2p1_app
  USE PP3D_MPI, ONLY: myMPI_Barrier
  USE POST_UTILS, ONLY: handle_statistics, print_time, sim_finalize
  USE TRANSPORT_Q2P1, ONLY: Transport_q2p1_UxyzP_el

  INTEGER :: ufile
  REAL :: dout = 0.0
  REAL :: tt0 = 0.0
  REAL :: dtt0 = 0.0

  CALL init_q2p1_app(ufile)

  CALL ZTIME(tt0)
  CALL ZTIME(dtt0)

  dout = REAL(INT(timens/dtgmv)+1)*dtgmv

  DO itns=1,nitns
    itnsr = 0
    timnsh = timens
    dt = tstep
    timens = timens + dt

    CALL Transport_q2p1_UxyzP_el(ufile,inonln_u,itns)

    CALL postprocessing_app(dout, inonln_u, inonln_t, ufile)
    CALL print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)
    CALL handle_statistics(tt0,itns)

    istep_ns = istep_ns + 1
    IF (timemx.LE.(timens+1.0d-10)) EXIT
  END DO

  CALL myMPI_Barrier()
  CALL sim_finalize(tt0,ufile)

END PROGRAM Q2P1_EL_PIPEFLOW
