PROGRAM Q2P1_DEVEL

  include 'defs_include.h'

  use solution_io, only: postprocessing_app,write_sol_to_file

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize
                         
  use Transport_Q2P1,  only:        updateFBMGeometry
  USE Sigma_User, ONLY: mySigma,myProcess,mySetup,bKTPRelease
  
 use Transport_Q1,  only: Transport_LinScalar_General,Boundary_LinSc_Val_General,AddSource_General
 use var_QuadScalar, only : dTimeStepEnlargmentFactor,istep_ns

  REAL*8 ViscosityModel
  REAL*8 dCharVisco,dCharSize,dCharVelo,dCharShear,TimeStep
  CHARACTER sTimeStep*(9)

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog,iXX
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0

  mySetup%bAutomaticTimeStepControl = .false.
  
  call init_q2p1_ext(ufile)

#if !defined WIN32
  !IF (bKTPRelease) CALL xSEND_START()
#endif

  CALL ZTIME(tt0)
  call ztime(dtt0)

  timens=0d0
  istep_ns = 1
  
  dout = Real(INT(timens/dtgmv)+1)*dtgmv
  
  mySetup%bAutomaticTimeStepControl = .FALSE.
  dTimeStepEnlargmentFactor = 1d0
  
  !-------MAIN LOOP-------
  
  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt
  inonln_u = 2
  inonln_t = 2
  
  if (itns.eq.1) CALL updateFBMGeometry()

  CALL Transport_LinScalar_General(Boundary_LinSc_Val_General,AddSource_General,ufile,inonln_t)

  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  call sim_finalize(tt0,ufile)

END PROGRAM Q2P1_DEVEL
