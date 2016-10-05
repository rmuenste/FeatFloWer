      PROGRAM LES

      PARAMETER (NNWORK=0,NNARR=299)
      INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

      INTEGER            :: KWORK(1)
      REAL               :: VWORK(1)
      DOUBLE PRECISION   :: DWORK(NNWORK)

      COMMON       NWORK,IWORK,IWMAX,L,DWORK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))

      ! Default names for output devices etc.
      CALL ZINIT(NNWORK,'_data/feat.msg','_data/les.err',
     *     '_data/les.prt','_data/les.sys','_data/les.trc') 

      CALL Main()

      END


      SUBROUTINE Main()

      USE def_FEAT
      USE PLinScalar, ONLY : Init_PLinScalar,InitCond_PLinLS,
     * UpdateAuxVariables,Transport_PLinLS,Reinitialize_PLinLS,
     * Reinit_Interphase,dMaxSTF
      USE QuadScalar, ONLY : Init_QuadScalar_Stuctures,
     * InitCond_QuadScalar,Transport_QuadScalar,ProlongateSolution,
     * ResetTimer,bTracer,bViscoElastic,StaticMeshAdaptation
      USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures,
     * Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
      USE LinScalar, ONLY : Init_LinScalar,InitCond_LinScalar,
     * Transport_LinScalar
      USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
      USE var_QuadScalar, ONLY : myStat,cFBM_File

      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
      INTEGER iOGMV,iTout
      CHARACTER command*100
      CHARACTER CPP3D*(60)
      LOGICAL bstop
      DATA tout=0.0
C
      ! Initialization for FEATFLOW
      CALL General_init(79,mfile)
C


      CALL StaticMeshAdaptation()
C


      CALL Init_QuadScalar_Stuctures(mfile)



      IF(bViscoElastic)CALL Init_ViscoScalar_Stuctures(mfile)


      CALL Init_LinScalar

      CALL InitCond_LinScalar()

                
C
      IF (ISTART.EQ.0) THEN
       IF (myid.ne.0) CALL CreateDumpStructures(1)
!       CALL FBM_FromFile(cFBM_File)
       CALL InitCond_QuadScalar()
!       Initial profiles for the viscoelastic stresses
        IF(bViscoElastic)CALL IniProf_ViscoScalar()
      ELSE
       IF (ISTART.EQ.1) THEN
        IF (myid.ne.0) CALL CreateDumpStructures(1)
!        CALL FBM_FromFile(cFBM_File)
        CALL SolFromFile(CSTART,1)
!        CALL Load_DUMPProfiles(CSTART)
       ELSE
        IF (myid.ne.0) CALL CreateDumpStructures(0)
!        CALL FBM_FromFile(cFBM_File)
        CALL SolFromFile(CSTART,0)
!        CALL Load_LowDUMPProfiles(CSTART)
        CALL ProlongateSolution()
        IF (myid.ne.0) CALL CreateDumpStructures(1)
       END IF
      END IF

      tout = DBLE(INT(timens/dtgmv)+1)*dtgmv
!      CALL Output_Profiles(0)
!      PAUSE

      CALL ZTIME(ttt0)

      DO itns=1,nitns

!        ! New timestep estimate !
!        CALL TimeStepCtrlST(MFILE,MTERM,dMaxSTF,IPR,tstep)

       itnsr=0
       timnsh=timens
       dt=tstep
       timens=timens+dt
C
       IF (bstop) STOP

       ! Viscoelastic transport equations
       IF(bViscoElastic)CALL Transport_ViscoScalar(mfile)

       ! Solve Navier-Stokes
       CALL Transport_QuadScalar(mfile,inonln_u,itns)
C
       IF (bTracer) THEN
        CALL Transport_LinScalar(mfile,inonln_t)
       ELSE
        inonln_t = 2
       END IF
C
c      ! Output the solution in GMV or GiD format
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
       !call TimeStepCtrlFish(MFILE,MTERM,tstep,timens)

       ! Interaction from user
       CALL ProcessControl(mfile,mterm)

       ! Statistics reset
       IF (itns.eq.1) THEN
!        CALL ZTIME(tttx)
!        CALL StatOut(tttx-ttt0,0)
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
       write(MTERM,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')
     * "time:", timens,"/",timemx," | itns:",itns,"/",nitns,
     * " | dt:",tstep
       write(MFILE,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')
     * "time:", timens,"/",timemx," | itns:",itns,"/",nitns,
     * " | dt:",tstep
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

      END 

      SUBROUTINE TimeStepCtrl(dt,inlU,inlT,MFILE)

      USE PP3D_MPI,only :myid,ShowID

      INTEGER M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8,
     *        IEPSAD,IADIN,IREPIT,IADTIM
      REAL*8  TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *        EPSADU,PRDIF1,PRDIF2
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
      INTEGER inlU,inlT,MFILE
      INTEGER :: iMem,nMEm=2
      REAL*8 dt, dt_old
      CHARACTER char_dt*9
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

      END

      SUBROUTINE TimeStepCtrlST(MFILE,MTERM,dVal,IR,dt)
      USE PP3D_MPI,only :myid,ShowID
      IMPLICIT NONE

      INTEGER MFILE,MTERM,IR
      REAL*8 dVal,dt,DR
      CHARACTER char_dt*9

      !!!!! Time Step function = f(Surface Tension)  !!!!
      dt = MAX(MIN(1d-3,(12d-4 - 1.5d-8*dVal)),3d-4)
      dt = 0.00500d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!! Interphase renewal function = f(Surface Tension)  !!!!
      DR = MAX(MIN(50.0,(70.0 - 0.001*dVal)),10.0)
      DR = 5d1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IR = NINT(DR)

      WRITE(char_dt,'(D9.2)') dt
      READ (char_dt,'(D9.2)') dt

      IF (myid.eq.showid) 
     *   write(MFILE,"(2(A8G12.4),I5)") "dMaxSTF=",dVal,"dt=",dt,IR
      IF (myid.eq.showid) 
     *   write(MTERM,"(2(A8G12.4),I5)") "dMaxSTF=",dVal,"dt=",dt,IR

      END

      SUBROUTINE TimeStepCtrlFish(MFILE,MTERM,dt,time)
      USE PP3D_MPI,only :myid,ShowID
      IMPLICIT NONE

      INTEGER MFILE,MTERM,IR
      REAL*8 time,dt
      CHARACTER char_dt*9

      if((time.gt.2.91d0))then
        dt=0.0001d0
      end if

!      if(time.gt.5.0d0)then
!        dt=0.0005d0
!      end if
   

      WRITE(char_dt,'(D9.2)') dt
      READ (char_dt,'(D9.2)') dt

      IF (myid.eq.showid) 
     *   write(MFILE,"(2(A8G12.4))") "timenvs=",time,"dt=",dt
      IF (myid.eq.showid) 
     *   write(MTERM,"(2(A8G12.4))") "timenvs=",time,"dt=",dt

      END
