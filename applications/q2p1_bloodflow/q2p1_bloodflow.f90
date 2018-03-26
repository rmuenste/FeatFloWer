PROGRAM Q2P1_BLOODFLOW

  include 'defs_include.h'

  use solution_io, only: postprocessing_app,write_sol_to_file

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog,iXX
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0

  !-------INIT PHASE-------

  call init_q2p1_app(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dout = Real(INT(timens/dtgmv)+1)*dtgmv
  
!   CALL ReadS3Dfile('_adc/egeplast/rheo.s3d')

  !-------MAIN LOOP-------

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt
  inonln_u = 2
  inonln_t = 2

  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
  CALL Transport_q2p1_UxyzP_fc_ext(ufile,inonln_u,itns)


!   write(*,*) myid, ' reached this stage'
  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  call sim_finalize(tt0,ufile)

END PROGRAM Q2P1_BLOODFLOW
