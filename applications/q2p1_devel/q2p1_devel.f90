PROGRAM Q2P1_DEVEL

  include 'defs_include.h'

  use solution_io, only: postprocessing_app,write_sol_to_file

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize
                         
  use Transport_Q2P1,  only:        updateFBMGeometry
  USE Sigma_User, ONLY: mySigma,myProcess,mySetup
  
!  use Transport_Q1,  only: Transport_GeneralLinScalar,Boundary_LinSc_Val_Weber,AddSource

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

  CALL ZTIME(tt0)
  call ztime(dtt0)

!=====================================================================================
!=====================================================================================
!=====================================================================================
  IF (myid.eq.1) THEN
   WRITE(MTERM,'(A)') 
   WRITE(MTERM,'(A)') " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
   WRITE(MTERM,'(A)') " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
  END IF
  IF (mySetup%bAutomaticTimeStepControl) THEN
    ! get the characteristic viscosity for characteristic shear rate (10.0[1/s])
   dCharSize      = 0.5d0*myProcess%MaxInflowDiameter
   dCharShear     = 3d0
   dCharVisco     = ViscosityModel(dCharShear,myProcess%T0)
   TimeStep       = 2.5d-3 * (dCharSize/dCharVisco)
   WRITE(sTimeStep,'(ES9.1)') TimeStep
   READ(sTimeStep,*) TimeStep

   IF (myid.eq.1) THEN
    WRITE(MTERM,'(A,5ES12.4,A)') " Characteristic size[cm],shear[1/s]_E/U: ",dCharSize,dCharShear
    WRITE(ufile,'(A,5ES12.4,A)') " Characteristic size[cm],shear[1/s]_E/U: ",dCharSize,dCharShear
    WRITE(MTERM,'(A,2ES12.4,A)') " Characteristic viscosity [Pa.s] and corresponding Timestep [s]: ",0.1d0*dCharVisco,TimeStep
    WRITE(ufile,'(A,2ES12.4,A)') " Characteristic viscosity [Pa.s] and corresponding Timestep [s]: ",0.1d0*dCharVisco,TimeStep
   END IF
   
   CALL AdjustTimeStepping(TimeStep)
  END IF
     
  IF (myid.eq.1) THEN
    WRITE(MTERM,'(A,3ES12.4,I10)') " TSTEP,DTGMV,TIMEMX,NITNS ",TSTEP,DTGMV, TIMEMX, NITNS
    WRITE(ufile,'(A,3ES12.4,I10)') " TSTEP,DTGMV,TIMEMX,NITNS ",TSTEP,DTGMV, TIMEMX, NITNS
    WRITE(MTERM,'(A)') " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
    WRITE(MTERM,'(A)') " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
    WRITE(MTERM,'(A)') 
  END IF
!=====================================================================================
!=====================================================================================
!=====================================================================================

  dout = Real(INT(timens/dtgmv)+1)*dtgmv
  
!   CALL updateFBMGeometry()
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

!   IF (bTracer) THEN
!     ! Solve transport equation for linear scalar
!     CALL Transport_GeneralLinScalar(Boundary_LinSc_Val_Weber,AddSource,ufile,inonln_t)
!   ELSE
!     inonln_t = 2
!   END IF

!   write(*,*) myid, ' reached this stage'
  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  call sim_finalize(tt0,ufile)

END PROGRAM Q2P1_DEVEL
