PROGRAM Q2P1_CREEP
! hello

  include 'defs_include.h'

  use solution_io, only: postprocessing_app

  use app_initialization, only: init_q2p1_app
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  real               :: dtt10 = 0.0

  !-------INIT PHASE-------

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
!  if (myid.ne.0)then
!    call step_simulation()
!  end if

  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
  CALL Transport_q2p1_UxyzP_fc_ext(ufile,inonln_u,itns)
  !call step_simulation()
!
!  IF (bTracer) THEN
!    ! Solve transport equation for linear scalar
!    CALL Transport_LinScalar(ufile,inonln_t)
!  ELSE
!    inonln_t = 2
!  END IF

  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT
  
  END DO

  call MPI_Barrier(MPI_COMM_WORLD)
  call sim_finalize(tt0,ufile)

END PROGRAM Q2P1_CREEP
