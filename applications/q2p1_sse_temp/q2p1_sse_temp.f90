PROGRAM Q2P1_DEVEL

  include 'defs_include.h'

  use solution_io, only: write_sol_to_file
  
  use visualization_out, only: viz_output_fields_Simple

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize
                         
  use Transport_Q2P1,  only: TemporalFieldInterpolator,QuadSc,LinSc
  
  USE Sigma_User, ONLY: mySigma,myProcess,mySetup,bKTPRelease,myTransientSolution
  
  use Transport_Q1,  only: Transport_LinScalar_XSE,Boundary_LinSc_Val_XSE,AddSource_XSE,&
                           Assemble_LinScOperators_XSE
                           
  use var_QuadScalar, only : istep_ns,Viscosity, Screw, Shell, Shearrate,&
                      mg_mesh,myExport

  USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI
  
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
  real*8 dTimeStep,dPeriod
  integer iRot,iStep,iSubStep

  call init_q2p1_ext(ufile)

#if !defined WIN32
  !IF (bKTPRelease) CALL xSEND_START()
#endif

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dPeriod = 6d1/myProcess%Umdr
  dTimeStep = dPeriod/DBLE(myProcess%nTimeLevels)
  tstep = dTimeStep/myTransientSolution%nTimeSubStep
  
  timens=0d0
  timemx=dPeriod*nitns
  istep_ns = 1
  itns = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  dout = 0d0
  dtgmv = dPeriod/DBLE(myProcess%nTimeLevels)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  
  !-------MAIN LOOP-------
  
  DO iRot=1,nitns
   DO iStep=0,myProcess%nTimeLevels-1
    DO iSubStep=1,myTransientSolution%nTimeSubStep
 
     itns = itns +1 
   
     itnsr=0
     timnsh=timens
     dt=tstep
     timens=timens+dt
     inonln_u = 2
     inonln_t = 2
   
     CALL TemporalFieldInterpolator(iStep,iSubStep)
     CALL Assemble_LinScOperators_XSE(ufile)
     CALL Transport_LinScalar_XSE(Boundary_LinSc_Val_XSE,AddSource_XSE,ufile,inonln_t)

     call print_time(timens, timemx, tstep, itns, nitns*myProcess%nTimeLevels*myTransientSolution%nTimeSubStep, ufile, uterm)

     call handle_statistics(tt0,itns)

     istep_ns = istep_ns + 1
    
     goto 55
    END DO
    
    myProcess%Angle = MOD(INT(myProcess%dAlpha)*(iStep+1),360/myProcess%Periodicity)
    CALL ZTIME(myStat%t0)
    IF (myProcess%Angle.eq.0d0) THEN
     call viz_output_fields_Simple(myExport, int(myProcess%Angle), QuadSc, LinSc, & 
          Viscosity, Screw, Shell, Shearrate, mg_mesh)
    END IF

    CALL Release_ListFiles_SSE_temp(int(myProcess%Angle))
!    call write_sol_to_file(insavn, timens)
    CALL ZTIME(myStat%t1)
    myStat%tDumpOut = myStat%tDumpOut + (myStat%t1-myStat%t0)
    
    
   END DO

  END DO

55 continue
  
  IF (myid.eq.showid) THEN
    WRITE(*,*) "PP3D_LES has successfully finished. "
    WRITE(filehandle,*) "PP3D_LES has successfully finished. "
  END IF

  CALL Barrier_myMPI()
  CALL MPI_Finalize(ierr)

END PROGRAM Q2P1_DEVEL
