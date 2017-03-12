SUBROUTINE mySolToFile(iOutput)
USE def_FEAT
USE QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE LinScalar,ONLY:Tracer
USE PP3D_MPI, ONLY:myid

IMPLICIT NONE
INTEGER iOutput

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------
INTEGER ifilen,iOut,nn
DATA ifilen/0/

IF (iOutput.LT.0) THEN
 ifilen=ifilen+1
 iOut=MOD(ifilen+insavn-1,insavn)+1
ELSE
 iOut = iOutput
END IF

CALL WriteSol_Velo(iOut,0,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)
nn = knel(nlmax)
CALL WriteSol_Pres(iOut,0,LinSc%ValP(NLMAX)%x,LinSc%AuxP(NLMAX)%x(1),&
     LinSc%AuxP(NLMAX)%x(nn+1),LinSc%AuxP(NLMAX)%x(2*nn+1),LinSc%AuxP(NLMAX)%x(3*nn+1),nn)
! CALL WriteSol_Coor(iOut,0,DWORK(L(KLCVG(NLMAX))),QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof)
CALL WriteSol_Time(iOut)

if(bViscoElastic)then
  CALL WriteSol_Visco(iOut,0)
end if

if(myid.eq.1 .and. myFBM%nParticles.gt.0)then
  call writeparticles(iOut)
end if

END SUBROUTINE mySolToFile
!
! ----------------------------------------------
!
SUBROUTINE myOutput_Profiles(iOutput)
USE def_FEAT
USE QuadScalar,ONLY:QuadSc,LinSc,PressureToGMV,&
    Viscosity,Distance,Distamce,mgNormShearStress
USE LinScalar,ONLY:Tracer
USE PP3D_MPI, ONLY:myid,showid,Comm_Summ
USE var_QuadScalar,ONLY:myExport,myFBM,mg_mesh
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
! USE PLinScalar,ONLY:PLinScP1toQ1,OutputInterphase,PLinLS,&
!                dNorm,IntPhaseElem,FracFieldQ1
IMPLICIT NONE
INTEGER iOutput,mfile

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------

IF     (myExport%Format.EQ."GMV") THEN

 IF (myid.NE.0) THEN
  NLMAX = NLMAX + 1
  ILEV = myExport%Level
  CALL SETLEV(2)
  CALL Output_GMV_fields(iOutput,&
    mg_mesh%level(ILEV)%dcorvg,&
    mg_mesh%level(ILEV)%kvert)
  NLMAX = NLMAX - 1
 END IF

ELSEIF (myExport%Format.EQ."VTK") THEN

 IF (myid.NE.0) THEN
  NLMAX = NLMAX + 1
  ILEV = myExport%Level
  CALL SETLEV(2)
  CALL Output_VTK_piece(iOutput,&
    mg_mesh%level(ILEV)%dcorvg,&
    mg_mesh%level(ILEV)%kvert)
  NLMAX = NLMAX - 1
 ELSE
  CALL Output_VTK_main(iOutput)
 END IF

END IF

END SUBROUTINE 
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
!
! ----------------------------------------------
!
subroutine postprocessing_fc_ext(tout, iogmv, inlU,inlT,MFILE)

include 'defs_include.h'

implicit none

integer, intent(inout) :: iogmv
real, intent(inout) :: tout
INTEGER :: inlU,inlT,MFILE

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
CALL TimeStepCtrl(tstep,inlU,inlT,mfile)

! Interaction from user
CALL ProcessControl(mfile,mterm)

end subroutine postprocessing_fc_ext
!
! ----------------------------------------------
!
subroutine handle_statistics(dttt0, dtttx, istepns)

include 'defs_include.h'

implicit none

real, intent(inout) :: dttt0, dtttx
integer, intent(inout) :: istepns

! Statistics reset
IF (istepns.eq.1) THEN
  CALL ResetTimer()
  CALL ZTIME(dttt0)
END IF

IF (MOD(istepns,10).EQ.5) THEN
  CALL ZTIME(dtttx)
  CALL StatOut(dtttx-dttt0,0)
END IF

end subroutine handle_statistics
!
! ----------------------------------------------
!
subroutine print_time(dtimens, dtimemx, dt, istepns, istepmaxns, ufile,uterm)

include 'defs_include.h'
USE PP3D_MPI,only :myid,ShowID

implicit none

real*8, intent(inout) :: dtimens, dtimemx, dt
integer, intent(inout) :: istepns , istepmaxns
integer, intent(inout) :: ufile, uterm

IF (myid.eq.showid) THEN
  write(uTERM,5)
  write(uFILE,5)
  write(uTERM,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
    "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
    " | dt:",dt
  write(uFILE,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
    "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
    " | dt:",dt
  write(uTERM,5)
  write(uFILE,5)
END IF

5 FORMAT(104('='))

end subroutine print_time
!
! ----------------------------------------------
!
subroutine sim_finalize(dttt0, dttt1, mfile)

include 'defs_include.h'

implicit none

real, intent(inout) :: dttt0, dttt1
integer, intent(in) :: mfile

CALL ZTIME(dttt1)

CALL StatOut(dttt1-dttt0,MFILE)
CALL StatOut(dttt1-dttt0,MTERM)

CALL Finalize(MFILE,MTERM)

end subroutine sim_finalize
