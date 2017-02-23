PROGRAM Q2P1_FC_EXT

  USE def_FEAT
  USE PLinScalar, ONLY : Init_PLinScalar,InitCond_PLinLS, &
    UpdateAuxVariables,Transport_PLinLS,Reinitialize_PLinLS, &
    Reinit_Interphase,dMaxSTF
  USE QuadScalar, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar,Transport_QuadScalar,ProlongateSolution, &
    ResetTimer,bTracer,bViscoElastic,StaticMeshAdaptation
  USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
    Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE LinScalar, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  logical            :: bstop
  real               :: tout = 0.0


  ! Initialization for FEATFLOW
  CALL General_init_ext(79,mfile)

  CALL Init_QuadScalar_Stuctures(mfile)

  IF(bViscoElastic)CALL Init_ViscoScalar_Stuctures(mfile)

  CALL Init_LinScalar

  CALL InitCond_LinScalar()

  IF (ISTART.EQ.0) THEN
    IF (myid.ne.0) CALL CreateDumpStructures(1)
    CALL InitCond_QuadScalar()
    IF(bViscoElastic)CALL IniProf_ViscoScalar()
  ELSE
    IF (ISTART.EQ.1) THEN
      IF (myid.ne.0) CALL CreateDumpStructures(1)
      CALL SolFromFile(CSTART,1)
    ELSE
      IF (myid.ne.0) CALL CreateDumpStructures(0)
      CALL SolFromFile(CSTART,0)
      CALL ProlongateSolution()
      IF (myid.ne.0) CALL CreateDumpStructures(1)
    END IF
  END IF

  tout = Real(INT(timens/dtgmv)+1)*dtgmv
  !      CALL Output_Profiles(0)
  !      PAUSE

  CALL ZTIME(ttt0)

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  IF (bstop) STOP

  ! Viscoelastic transport equations
  IF(bViscoElastic)CALL Transport_ViscoScalar(mfile)

  ! Solve Navier-Stokes
  CALL Transport_QuadScalar(mfile,inonln_u,itns)

  IF (bTracer) THEN
    CALL Transport_LinScalar(mfile,inonln_t)
  ELSE
    inonln_t = 2
  END IF

  ! Output the solution in GMV or GiD format
  IF (itns.eq.1) THEN
    CALL ZTIME(myStat%t0)
    CALL Output_Profiles(0)
    CALL ZTIME(myStat%t1)
    myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
  END IF
  IF(tout.LE.(timens+1e-10)) THEN
    iOGMV = NINT(timens/dtgmv)
    IF (itns.ne.1) THEN
      CALL ZTIME(myStat%t0)
      CALL Output_Profiles(iOGMV)
      CALL ZTIME(myStat%t1)
      myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
    END IF
    tout=tout+dtgmv
    ! Save intermediate solution to a dump file
    IF (insav.NE.0.AND.itns.NE.1) THEN
      IF (MOD(iOGMV,insav).EQ.0) THEN
        CALL ZTIME(myStat%t0)
        CALL SolToFile(-1)
        !          CALL Output_DUMPProfiles()
        CALL FBM_ToFile()
        CALL ZTIME(myStat%t1)
        myStat%tDumpOut = myStat%tDumpOut + (myStat%t1-myStat%t0)
      END IF
    END IF
  END IF

  ! Timestep control
  CALL TimeStepCtrl(tstep,inonln_u,inonln_t,mfile)

  ! Interaction from user
  CALL ProcessControl(mfile,mterm)

  ! Statistics reset
  IF (itns.eq.1) THEN
    CALL ResetTimer()
    CALL ZTIME(ttt0)
  END IF

  IF (MOD(itns,10).EQ.5) THEN
    CALL ZTIME(tttx)
    CALL StatOut(tttx-ttt0,0)
  END IF

  IF (myid.eq.showid) THEN
    write(MTERM,5)
    write(MFILE,5)
    write(MTERM,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
      "time:", timens,"/",timemx," | itns:",itns,"/",nitns,&
      " | dt:",tstep
    write(MFILE,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
      "time:", timens,"/",timemx," | itns:",itns,"/",nitns,&
      " | dt:",tstep
    write(MTERM,5)
    write(MFILE,5)
  END IF
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) GOTO 1

  END DO

  1     CONTINUE

  CALL ZTIME(ttt1)

  CALL StatOut(ttt1-ttt0,MFILE)
  CALL StatOut(ttt1-ttt0,MTERM)

  CALL Finalize(MFILE,MTERM)

  5     FORMAT(104('='))

END PROGRAM Q2P1_FC_EXT
!
!----------------------------------------------
!
SUBROUTINE TimeStepCtrl(dt,inlU,inlT,MFILE)

  USE PP3D_MPI,only :myid,ShowID

  INTEGER IADTIM

  REAL*8  TIMEMX,DTMIN,DTMAX

  COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,IADTIM

  INTEGER :: inlU,inlT,MFILE
  INTEGER :: iMem,nMEm=2
  REAL*8  :: dt, dt_old
  CHARACTER(len=9) :: char_dt
  DATA iMem/0/

  IF (IADTIM.EQ.0) RETURN

  iMem = iMem + 1
  dt_old = dt
  IF (((inlU.GT.3).OR. (inlT.GT.5)).AND.iMem.GE.nMem) THEN
    dt=MAX(dt/1.1d0,DTMIN)
    WRITE(char_dt,'(D9.2)') dt
    READ (char_dt,'(D9.2)') dt
  END IF
  IF (((inlU.LT.3).AND.(inlT.LT.4)).AND.iMem.GE.nMem) THEN
    dt=MIN(1.1d0*dt,DTMAX)
    WRITE(char_dt,'(D9.2)') dt
    READ (char_dt,'(D9.2)') dt
  END IF

  IF (dt.NE.dt_old.AND.myid.eq.ShowID) THEN
    WRITE(MTERM,1) dt_old,dt
    WRITE(MFILE,1) dt_old,dt
    !       WRITE(MTERM,*) 
    !       WRITE(MFILE,*) 
  END IF

  IF (dt.NE.dt_old) iMem = 0

  1     FORMAT('Time step change from ',D9.2,' to ',D9.2)

END SUBROUTINE TimeStepCtrl

