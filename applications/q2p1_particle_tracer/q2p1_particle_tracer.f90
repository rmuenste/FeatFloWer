PROGRAM Q2P1_GEOM_TEST

  include 'defs_include.h'

  use solution_io, only: postprocessing_app

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize

  use fbmaux                         

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  real*8, dimension(8) :: x,y,z
  real*8 :: xx,yy,zz,DPARX,DPARY,DPARZ
  integer :: ielem = 0

  !-------INIT PHASE-------
  x = (/ 1.0000000000000000 , 1.0000000000000000 , -1.0000001192092896 , -0.9999996423721313 , 1.0000004768371582 , 0.9999993443489075 , -1.0000003576278687 , -0.9999999403953552 /)

  y = (/ 0.9999999403953552 , -1.0000000000000000 , -0.9999998211860657 , 1.0000003576278687 , 0.3938853740692139 , -1.0000005960464478 , -1.0654625892639160 , 0.6892766952514648 /)

  z = (/ -1.000000000000000, -1.000000000000000, -1.000000000000000, -1.000000000000000, 1.0000000000000000 , 1.4755718708038330 , 1.0597654581069946 , 0.7847060561180115 /)

  xx =  0.8
  yy = -0.4
  zz =  0.35

  if(fbmaux_PointInHex(xx,yy,zz,x,y,z,DPARX,DPARY,DPARZ,ielem))then
    write(*,*)'in'
  else
    write(*,*)'out'
  end if


!  ! !
!  call init_q2p1_app(ufile)

!  CALL ZTIME(tt0)
!  call ztime(dtt0)

!  dout = Real(INT(timens/dtgmv)+1)*dtgmv

!  !-------MAIN LOOP-------

!  DO itns=1,nitns

!  itnsr=0
!  timnsh=timens
!  dt=tstep
!  timens=timens+dt

!  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
!  CALL Transport_q2p1_UxyzP_fc_ext(ufile,inonln_u,itns)

!  inonln_t = 2

!  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

!  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

!  call handle_statistics(tt0,itns)

!  ! Exit if done
!  IF (timemx.LE.(timens+1D-10)) EXIT

!  END DO

!  call sim_finalize(tt0,ufile)

END PROGRAM Q2P1_GEOM_TEST
