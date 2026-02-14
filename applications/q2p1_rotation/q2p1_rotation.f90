PROGRAM Q2P1_SPAN

  ! Include definitions
  include 'defs_include.h'

  use solution_io, only: postprocessing_app

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize
  USE var_QuadScalar, ONLY :  myTimer

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  real               :: dtt10 = 0.0

  character(len=100) :: arg
  character(len=100) :: version_string
  character(len=100) :: git_commit_hash_trim
  logical :: show_version

  call init_q2p1_app(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dout = Real(INT(timens/dtgmv)+1)*dtgmv

  !-------MAIN LOOP-------

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
  CALL Transport_q2p1_UxyzP_fc_ext(ufile,inonln_u,itns)

  IF (bTracer) THEN
    ! Solve transport equation for linear scalar
    CALL Transport_LinScalar(ufile,inonln_t)
  ELSE
    inonln_t = 2
  END IF

  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

!   write(*,'(I0,A,6I5,6ES12.4)') myid," : ",myTimer%n,myTimer%t
  CALL Output_MPI_Timings()
  
  call sim_finalize(tt0,ufile)


END PROGRAM Q2P1_SPAN
