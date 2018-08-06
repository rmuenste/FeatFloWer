PROGRAM Q2P1_CC

  include 'defs_include.h'

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real*8             :: dout = 0.0
  integer            :: ufile,ilog
  real*8             :: tt0 = 0.0
  real*8             :: dt = 0.0

  !-------INIT PHASE-------

  call init_q2p1_cc(ufile)

  CALL ZTIME(tt0)


  dout = Real(INT(timens/dtgmv)+1)*dtgmv

  IF (ccParams%BDF.ne.0) THEN
  IF (myid.ne.master) THEN
  QuadSc%valU_old1 = QuadSc%valU
  QuadSc%valV_old1 = QuadSc%valV
  QuadSc%valW_old1 = QuadSc%valW
  QuadSc%valU_help = QuadSc%valU
  QuadSc%valV_help = QuadSc%valV
  QuadSc%valW_help = QuadSc%valW
  QuadSc%valU_old2 = QuadSc%valU
  QuadSc%valV_old2 = QuadSc%valV
  QuadSc%valW_old2 = QuadSc%valW
  END IF
  END IF

  !-------MAIN LOOP-------
  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt


  ! Solve Navier-Stokes 
  CALL Transport_q2p1_UxyzP_cc(ufile,inonln_u)

  inonln_t = 2

  call postprocessing_fc_ext(dout, iogmv, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile)

  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  call sim_finalize(tt0,ufile)

END PROGRAM Q2P1_CC
