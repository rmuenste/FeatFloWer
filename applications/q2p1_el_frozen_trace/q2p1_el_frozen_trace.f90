PROGRAM Q2P1_EL_FROZEN_TRACE

  include 'defs_include.h'

  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  use el_frozen_driver, only: el_run_frozen_particle_pass, el_step_frozen_particles

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize

  integer            :: ufile,ilog
  real               :: tt0 = 0.0

  !-------INIT PHASE-------

  call init_q2p1_ext(ufile)

  CALL ZTIME(tt0)

  !-------MAIN LOOP-------
  DO itns=1,nitns

    itnsr=0
    timnsh=timens
    dt=tstep
    timens=timens+dt

    ! Keep the loaded carrier field frozen. This loop advances only the
    ! application time bookkeeping and diagnostics shell.
    call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)
    call el_run_frozen_particle_pass(ufile, itns)
    call el_step_frozen_particles(ufile, itns)
    call handle_statistics(tt0,itns)

    istep_ns = istep_ns + 1
    IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  call sim_finalize(tt0,ufile)

END PROGRAM Q2P1_EL_FROZEN_TRACE
