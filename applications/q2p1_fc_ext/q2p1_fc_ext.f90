PROGRAM Q2P1_FC_EXT

  include 'defs_include.h'

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  logical            :: bstop
  real               :: tout = 0.0
  integer            :: ufile, uterm

  !-------INIT PHASE-------

  call init_q2p1_ext()
  ufile=mfile
  uterm=mterm

  tout = Real(INT(timens/dtgmv)+1)*dtgmv

  CALL ZTIME(ttt0)

  !-------MAIN LOOP-------

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  IF (bstop) STOP

  ! Viscoelastic transport equations
  IF(bViscoElastic)CALL Transport_ViscoScalar(mfile)

  ! Solve Navier-Stokes
  !CALL Transport_QuadScalar_fc_ext(mfile,inonln_u,itns)
  CALL Transport_QuadScalar(ufile,inonln_u,itns)

  IF (bTracer) THEN
    ! Solve transport equation for linear scalar
    CALL Transport_LinScalar(ufile,inonln_t)
  ELSE
    inonln_t = 2
  END IF

  call postprocessing_fc_ext(tout, iogmv, inonln_u, inonln_t,ufile)

  call handle_statistics(ttt0,tttx,itns)

  call print_time(timens,timemx,tstep, itns, nitns, ufile,uterm)

  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

5     FORMAT(104('='))

  call sim_finalize(ttt0,ttt1,mfile)

END PROGRAM Q2P1_FC_EXT
