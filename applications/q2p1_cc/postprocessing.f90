SUBROUTINE myOutput_Profiles(iOutput)
USE def_FEAT
USE Transport_Q2P1,ONLY:QuadSc,LinSc,PressureToGMV,&
    Viscosity,Distance,Distamce,mgNormShearStress
USE Transport_Q1, ONLY:Tracer
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

!if(myid.eq.1 .and. myFBM%nParticles.gt.0)then
!  call writeparticles(iOutput)
!end if

END SUBROUTINE 
!
!----------------------------------------------
!
SUBROUTINE TimeStepCtrl(dt,inlU,inlT, filehandle)

  USE PP3D_MPI,only :myid,ShowID

  INTEGER IADTIM

  REAL*8  TIMEMX,DTMIN,DTMAX

  COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,IADTIM

  integer, intent(in) :: filehandle

  INTEGER :: inlU,inlT
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
    WRITE(filehandle,1) dt_old,dt
  END IF

  IF (dt.NE.dt_old) iMem = 0

  1  FORMAT('Time step change from ',D9.2,' to ',D9.2)

END SUBROUTINE TimeStepCtrl
!
! ----------------------------------------------
!
subroutine postprocessing_fc_ext(dout, iogmv, inlU,inlT,filehandle)

include 'defs_include.h'

implicit none

integer, intent(in) :: filehandle

integer, intent(inout) :: iogmv
real, intent(inout) :: dout

INTEGER :: inlU,inlT,MFILE

! Output the solution in GMV or GiD format
IF (itns.eq.1) THEN
  CALL ZTIME(myStat%t0)
  CALL myOutput_Profiles(0)
  CALL ZTIME(myStat%t1)
  myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
END IF

IF(dout.LE.(timens+1e-10)) THEN

  iOGMV = NINT(timens/dtgmv)
!   IF (itns.ne.1) THEN
    CALL ZTIME(myStat%t0)
    CALL myOutput_Profiles(iOGMV)
    CALL ZTIME(myStat%t1)
    myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
!   END IF
  dout=dout+dtgmv

  ! Save intermediate solution to a dump file
  IF (insav.NE.0.AND.itns.NE.1) THEN
    IF (MOD(iOGMV,insav).EQ.0) THEN
      CALL ZTIME(myStat%t0)
      CALL mySolToFile(-1)
      !CALL myReleaseSmartDumpFiles(-1)
      CALL ZTIME(myStat%t1)
      myStat%tDumpOut = myStat%tDumpOut + (myStat%t1-myStat%t0)
    END IF
  END IF

END IF

! Timestep control
CALL TimeStepCtrl(tstep,inlU,inlT,filehandle)

! Interaction from user
CALL ProcessControl(filehandle,mterm)

end subroutine postprocessing_fc_ext
!
! ----------------------------------------------
!
subroutine handle_statistics(dttt0, istepns)

include 'defs_include.h'

implicit none

real, intent(inout) :: dttt0
integer, intent(inout) :: istepns
real :: dttx = 0.0

! Statistics reset
IF (istepns.eq.1) THEN
  CALL ResetTimer()
  CALL ZTIME(dttt0)
END IF

IF (MOD(istepns,10).EQ.5) THEN
  CALL ZTIME(dttx)
  CALL StatOut(dttx-dttt0,0)
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
subroutine release_mesh() 
USE PP3D_MPI, ONLY : myid,master,showid
USE var_QuadScalar, only : mg_mesh
implicit none

integer :: i
integer :: maxlevel

maxlevel = mg_mesh%maxlevel

if(associated(mg_mesh%level(maxlevel)%dcorvg))then
  deallocate(mg_mesh%level(maxlevel)%dcorvg)
  mg_mesh%level(maxlevel)%dcorvg => null()
end if

end subroutine release_mesh
!
! ----------------------------------------------
!
SUBROUTINE myStatOut(time_passed,myOutFile)
USE def_FEAT
USE PP3D_MPI, ONLY : myid,master,showid,subnodes
USE var_QuadScalar, ONLY : myStat,bNonNewtonian,myMatrixRenewal
IMPLICIT NONE

Real, intent(in) :: time_passed

REAL*8 daux,daux1,ds
INTEGER myFile,myOutFile,itms,istat
LOGICAL bExist

itms = min(itns-1,nitns-1)
ds = DBLE(subnodes)

IF (myid.eq.showid) THEN

      IF (myOutFile.eq.0) THEN
       myFile = 669
       OPEN (UNIT=myFile, FILE='_data/Statistics.txt',action='write',iostat=istat)
       if(istat .ne. 0)then
         write(*,*)"Could not open file for writing in StatOut(). "
       stop          
       end if
      ELSE
       myFile = myOutFile
      END IF

      daux = myStat%tKMat+myStat%tDMat+myStat%tMMat+myStat%tCMat+myStat%tSMat
      daux1 = myStat%tGMVOut +myStat%tDumpOut
      WRITE(myFile,*) 
      WRITE(myFile,8) " Overall time            ",time_passed
      WRITE(myFile,*)  
      WRITE(myFile,8) " Solving time            ",time_passed-daux-daux1
      WRITE(myFile,*) 
      WRITE(myFile,8) " Operator assembly time  ",daux
      WRITE(myFile,8) "  Convection matrix      ",myStat%tKMat
      WRITE(myFile,8) "  Deformation matrix     ",myStat%tSMat
      WRITE(myFile,8) "  Diffusion matrix       ",myStat%tDMat
      WRITE(myFile,8) "  Mass matrix            ",myStat%tMMat
      WRITE(myFile,8) "  Reactive term          ",myStat%tCMat
      WRITE(myFile,*) 
      WRITE(myFile,8) " Output time             ",daux1
      WRITE(myFile,8) "  VTK/GMV Output time    ",myStat%tGMVOut
      WRITE(myFile,8) "  Dump Output time       ",myStat%tDumpOut

      IF (myOutFile.eq.0) THEN
       CLOSE (myFile)
      END IF

END IF

8  FORMAT(A24,' : ',F12.4,'s')

END SUBROUTINE myStatOut
!
! ----------------------------------------------
!
subroutine sim_finalize(dttt0, filehandle)

USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI
USE var_QuadScalar, ONLY : myStat

real, intent(inout) :: dttt0
integer, intent(in) :: filehandle

integer :: ierr
integer :: terminal = 6
real :: time,time_passed

CALL ZTIME(myStat%t0)
! Save the final solution vector in unformatted form
CALL mySolToFile(-1)
!CALL myReleaseSmartDumpFiles(-1)
CALL ZTIME(myStat%t1)
myStat%tDumpOut = myStat%tDumpOut + (myStat%t1-myStat%t0)

CALL ZTIME(time)

time_passed = time - dttt0
CALL Barrier_myMPI()
CALL myStatOut(time_passed,filehandle)

CALL myStatOut(time_passed,terminal)


IF (myid.eq.showid) THEN
  WRITE(MTERM,*) "CC3D_iso_adaptive has successfully finished. "
  WRITE(filehandle,*) "CC3D_iso_adaptive has successfully finished. "
END IF

call release_mesh()

CALL Barrier_myMPI()
CALL MPI_Finalize(ierr)

end subroutine sim_finalize
!
! ----------------------------------------------
!
SUBROUTINE myCreateDumpStructures(iLevel)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE Transport_Q2P1,ONLY:myDump,mg_mesh

IMPLICIT NONE

INTEGER iLevel
INTEGER JEL,KEL,ivt,jvt,nLengthE,nLengthV
INTEGER iaux,jaux,jj,kv(8),II
LOGICAL ,ALLOCATABLE :: bGot(:)
! -------------- workspace -------------------

! IF (myid.NE.0) THEN

 NLMAX = NLMAX+iLevel

 ILEV = NLMIN
 NEL  = mg_mesh%level(ilev)%nel
 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3
 IF(ALLOCATED(myDump%Elements)) DEALLOCATE(myDump%Elements)
 IF(ALLOCATED(myDump%Vertices)) DEALLOCATE(myDump%Vertices)
 ALLOCATE(myDump%Elements(NEL,nLengthE))
 ALLOCATE(myDump%Vertices(NEL,nLengthV))

 DO IEL = 1,mg_mesh%level(1)%nel

  myDump%Elements(IEL,1) = IEL
  iaux = 1
  DO II=1+1,NLMAX
   jaux = iaux
   DO jel=1,jaux
    kel = myDump%Elements(iel,jel)
    CALL Get8Elem(mg_mesh%level(II)%kadj,&
      kv,kel)
    DO jj = 2,8
     iaux = iaux + 1
     myDump%Elements(iel,iaux) = kv(jj)
    END DO
   END DO
  END DO
 END DO 

 ALLOCATE(bGot(mg_mesh%level(NLMAX)%nvt))
 DO IEL = 1,mg_mesh%level(1)%nel
  bGot = .FALSE.
  iaux = 0
  DO JEL = 1,nLengthE
   KEL = myDump%Elements(IEL,JEL)
   
   CALL getVert(mg_mesh%level(NLMAX)%kvert,&
                kv,KEL)

   DO IVT = 1,8
    JVT = kv(IVT)
    IF (.NOT.bGot(JVT)) THEN
     iaux = iaux + 1
     myDump%Vertices(IEL,iaux) = JVT
     bGot(JVT) = .TRUE.
    END IF
   END DO

  END DO

!  IF (iaux.ne.729) WRITE(*,*) myid,iel,iaux
  
 END DO
 DEALLOCATE(bGot)

 NLMAX = NLMAX - iLevel
!   write(*,*) size(myDump%Vertices),"asdas dasd sad sa",myid
!   pause

! END IF


 CONTAINS

 SUBROUTINE Get8Elem(KADJ,k,el)
 INTEGER KADJ(6,*),k(8),el

 k(1) = el
 k(2) = KADJ(3,k(1))
 k(3) = KADJ(3,k(2))
 k(4) = KADJ(3,k(3))
 k(5) = KADJ(6,k(1))
 k(6) = KADJ(3,k(5))
 k(7) = KADJ(3,k(6))
 k(8) = KADJ(3,k(7))

 END SUBROUTINE Get8Elem

 SUBROUTINE getVert(BigKv,SmallKv,elem)
 INTEGER BigKv(8,*),SmallKv(8),elem

 SmallKv(:) = BigKv(:,elem)

 END SUBROUTINE getVert

END SUBROUTINE myCreateDumpStructures
!
! ----------------------------------------------
!
SUBROUTINE myReleaseSmartDumpFiles(iOutO)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier
USE Transport_Q2P1,ONLY:QuadSc,LinSc,myDump
USE PLinScalar,ONLY:PLinLS
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
CHARACTER cOutFile*20,command*100
INTEGER iOutO
INTEGER iOut,myCoarseElem,ivt,jvt,jel,kel,nLengthE,nLengthV
INTEGER iP,lCmnd
INTEGER ifilen,i
DATA ifilen/0/

IF (myid.NE.0) THEN

 NLMAX = NLMAX+1

 IF (iOutO.LT.0) THEN
  ifilen=ifilen+1
  iOut=MOD(ifilen+insavn-1,insavn)+1
 ELSE
  iOut = iOuto
 END IF

 ILEV = NLMIN

 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3

 IF (myid.eq.1) THEN

  command = 'mkdir _dump'
  WRITE(*,*) 'Executing command:"',TRIM(ADJUSTL(command)),'"'
  CALL system(TRIM(ADJUSTL(command)))
  lCmnd = LEN(TRIM(ADJUSTL(command)))
  IF     (iOut.LT.10    ) THEN
   WRITE(command(lCmnd+1:lCmnd+3),'(A2,I1)') '/0',iOut
  ELSEIF (iOut.LT.100   ) THEN
   WRITE(command(lCmnd+1:lCmnd+3),'(A1,I2)') '/',iOut
  ELSE
    WRITE(*,*) "Decrease the output index!"
    STOP
  END IF
  WRITE(*,*) 'Executing command:"',TRIM(ADJUSTL(command)),'"'
  CALL system(TRIM(ADJUSTL(command)))

 END IF

! pause
 CALL myMPI_Barrier()

 cOutFile = '_dump/00/******.prf'

 IF     (iOut.LT.10    ) THEN
  WRITE(cOutFile(7:8),'(A1,I1)') '0',iOut
 ELSEIF (iOut.LT.100   ) THEN
   WRITE(cOutFile(7:8),'(I2)') iOut
 ELSE
   WRITE(*,*) "Decrease the output index!"
   STOP
 END IF

 IF (myid.EQ.1) THEN
  WRITE(*,'(3A)') "Storing the ",TRIM(ADJUSTL(cOutFile))," series for solution backup ..."
 END IF

 DO IEL = 1,KNEL(1)
  myCoarseElem = coarse%myELEMLINK(IEL)
  IF     (myCoarseElem.LT.10    ) THEN
   WRITE(cOutFile(10:15),'(A5,I1)') '00000',myCoarseElem
  ELSEIF (myCoarseElem.LT.100   ) THEN
   WRITE(cOutFile(10:15),'(A4,I2)') '0000',myCoarseElem
  ELSEIF (myCoarseElem.LT.1000  ) THEN
   WRITE(cOutFile(10:15),'(A3,I3)') '000',myCoarseElem
  ELSEIF (myCoarseElem.LT.10000 ) THEN
   WRITE(cOutFile(10:15),'(A2,I4)') '00',myCoarseElem
  ELSEIF (myCoarseElem.LT.100000) THEN
   WRITE(cOutFile(10:15),'(A1,I5)') '0',myCoarseElem
  ELSE
   WRITE(*,*) "Decrease the problem size!"
   STOP
  END IF

  ! Velocities
  OPEN (FILE=TRIM(ADJUSTL(cOutFile)),UNIT = 547)
  WRITE(547,*)  "Velocities"
  DO ivt=1,nLengthV
   jvt = myDump%Vertices(IEL,ivt)
   WRITE(547,'(3D16.8)') QuadSc%valU(jvt),QuadSc%valV(jvt),QuadSc%valW(jvt)
  END DO

  ! Pressure
  WRITE(547,*)  "Pressure"
  DO jel=1,nLengthE
   kel = myDump%Elements(IEL,jel)
   IF (kel.le.knel(NLMAX-1)) THEN
    iP = 4*(kel-1)
   WRITE(547,'(4D16.8)') LinSc%valP(NLMAX-1)%x(iP+1:iP+4)
   END IF   
  END DO

  WRITE(547,*)  "EOF"
  CLOSE(547)

 END DO 

 IF (myid.eq.1) THEN
  cOutFile = '_dump/00/time.prf'
  IF     (iOut.LT.10    ) THEN
   WRITE(cOutFile(7:8),'(A1,I1)') '0',iOut
  ELSEIF (iOut.LT.100   ) THEN
    WRITE(cOutFile(7:8),'(I2)') iOut
  ELSE
    WRITE(*,*) "Decrease the output index!"
    STOP
  END IF

  OPEN (FILE=TRIM(ADJUSTL(cOutFile)),UNIT = 547)
  WRITE (547,*) "timens"
  WRITE (547,'(D16.8)') timens
  CLOSE(547)
 END IF

 NLMAX = NLMAX-1

END IF

END SUBROUTINE myReleaseSmartDumpFiles
!
!---------------------------------------------------------------------------
!
SUBROUTINE myLoadSmartDumpFiles(cFldrInFile,iLevel)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse
USE Transport_Q2P1,ONLY:QuadSc,LinSc,myDump
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE PLinScalar,ONLY:PLinLS
IMPLICIT NONE
INTEGER iLevel
CHARACTER cInFile*99,cFldrInFile*(*)
INTEGER myCoarseElem,ivt,jvt,jel,kel,nLengthE,nLengthV
INTEGER iP,lStr

IF (myid.NE.0) THEN

 NLMAX = NLMAX+iLevel

 ILEV = NLMIN

 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3

 cInFile = ' '
 WRITE(cInFile(1:),'(A)') TRIM(ADJUSTL(cFldrInFile))
 lStr = LEN(TRIM(ADJUSTL(cInFile)))

 IF (myid.EQ.1) THEN
  WRITE(*,'(3A)') "Reading the ",TRIM(ADJUSTL(cInFile))," series for initialization of profiles ..."
 END IF

 DO IEL = 1,KNEL(1)
  myCoarseElem = coarse%myELEMLINK(IEL)
  IF     (myCoarseElem.LT.10    ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A6,I1,A4)') '/00000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.100   ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A5,I2,A4)') '/0000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.1000  ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A4,I3,A4)') '/000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.10000 ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A3,I4,A4)') '/00',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.100000) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A2,I5,A4)') '/0',myCoarseElem,".prf"
  ELSE
   WRITE(*,*) "Decrease the problem size!"
   STOP
  END IF

  ! Velocities
  OPEN (FILE=TRIM(ADJUSTL(cInFile)),UNIT = 547)
  READ(547,*)  
  DO ivt=1,nLengthV
   jvt = myDump%Vertices(IEL,ivt)
   READ(547,*) QuadSc%valU(jvt),QuadSc%valV(jvt),QuadSc%valW(jvt)
  END DO

  ! Pressure
  READ(547,*)  
  DO jel=1,nLengthE
   kel = myDump%Elements(IEL,jel)
   IF (kel.le.knel(NLMAX-1)) THEN
    iP = 4*(kel-1)
   READ(547,*) LinSc%valP(NLMAX-1)%x(iP+1:iP+4)
   END IF   
  END DO

  CLOSE(547)

 END DO 
END IF

 cInFile = ' '
 WRITE(cInFile(1:),'(A)') TRIM(ADJUSTL(cFldrInFile))
 lStr = LEN(TRIM(ADJUSTL(cInFile)))
 WRITE(cInFile(lStr+1:lStr+9),'(A)') '/time.prf'

 OPEN (FILE=TRIM(ADJUSTL(cInFile)),UNIT = 547)
 READ (547,*) 
 READ (547,*) timens
 CLOSE(547)

IF (myid.NE.0) THEN

 NLMAX = NLMAX-iLevel

END IF

END SUBROUTINE myLoadSmartDumpFiles
!
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!
SUBROUTINE mySolToFile(iOutput)
USE def_FEAT
USE Transport_Q2P1,ONLY:QuadSc,LinSc,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE Transport_Q1,ONLY:Tracer
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

CALL myWriteSol_Velo(iOut,0,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)
nn = knel(nlmax)

CALL myWriteSol_Pres(iOut,0,LinSc%ValP(NLMAX)%x,LinSc%AuxP(NLMAX)%x(1),&
     LinSc%AuxP(NLMAX)%x(nn+1),LinSc%AuxP(NLMAX)%x(2*nn+1),LinSc%AuxP(NLMAX)%x(3*nn+1),nn)

!! myCALL WriteSol_Coor(iOut,0,DWORK(L(KLCVG(NLMAX))),QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof)

CALL myWriteSol_Time(iOut)

if(bViscoElastic)then
  CALL myWriteSol_Visco(iOut,0)
end if

END SUBROUTINE mySolToFile
!
! ----------------------------------------------
!
SUBROUTINE mySolFromFile(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid
USE def_FEAT
USE Transport_Q2P1,ONLY:QuadSc,LinSc,SetUp_myQ2Coor,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE Transport_Q1,ONLY:Tracer
IMPLICIT NONE
INTEGER mfile,iLevel,nn
CHARACTER cInFile*(60)

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
INTEGER nLengthV,nLengthE,LevDif
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

CALL myReadSol_Velo(cInFile,iLevel,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)
nn = knel(nlmax)
CALL myReadSol_Pres(cInFile,iLevel,LinSc%ValP(NLMAX)%x,LinSc%AuxP(NLMAX)%x(1),&
     LinSc%AuxP(NLMAX)%x(nn+1),LinSc%AuxP(NLMAX)%x(2*nn+1),LinSc%AuxP(NLMAX)%x(3*nn+1),nn)
! CALL myReadSol_Coor(cInFile,iLevel,DWORK(L(KLCVG(NLMAX))),QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof)
CALL myReadSol_Time(cInFile)

if(bViscoElastic)then
  CALL myReadSol_Visco(cInFile, iLevel)
end if


END SUBROUTINE mySolFromFile
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myWriteSol_Time(iOut)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,subnodes,RECVDD_myMPI,SENDDD_myMPI

INTEGER iInd
INTEGER pID
character cFile*(40)

 IF (myid.eq.showid) THEN
  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_Time.dmp'
  WRITE(MTERM,*) 'Releasing Time level into: "', TRIM(ADJUSTL(cFile)),'"'
 END IF

 IF (myid.eq.0) THEN
  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_Time.dmp'
  OPEN(321,FILE=TRIM(ADJUSTL(cFile)))
  WRITE(321,*) TIMENS
  CLOSE(321)
 END IF

END SUBROUTINE myWriteSol_Time
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myWriteSol_Velo(iInd,iiLev,Field1,Field2,Field3)
interface
  SUBROUTINE myWriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE Transport_Q2P1,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE myWriteSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*)
INTEGER        iInd,iiLev
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Velocity'

CALL myWriteSol(iInd,iType,iiLev,cFF,Field1,Field2,Field3)

END SUBROUTINE myWriteSol_Velo
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myWriteSol_Visco(iInd,iiLev)
interface
  SUBROUTINE myWriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE Transport_Q2P1,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE myWriteSol
end interface
INTEGER        iInd,iiLev
INTEGER        :: iType= 3
CHARACTER*(20) :: cFF='Stress'

CALL myWriteSol(iInd,iType,iiLev,cFF)

END SUBROUTINE myWriteSol_Visco
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myWriteSol_Pres(iInd,iiLev,Field1,Field2,Field3,Field4,Field5,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE myWriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE Transport_Q2P1,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE myWriteSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*),Field4(*),Field5(*)
INTEGER        i,iInd,nn,iiLev
INTEGER        :: iType= 2
CHARACTER*(20) :: cFF='Pressure'

IF (myid.ne.0) THEN
 DO i=1,nn
  Field2(i) = Field1(4*(i-1)+1)
  Field3(i) = Field1(4*(i-1)+2)
  Field4(i) = Field1(4*(i-1)+3)
  Field5(i) = Field1(4*(i-1)+4)
 END DO
!  WRITE(*,*) Field2(1:nn)
END IF


CALL myWriteSol(iInd,iType,iiLev,cFF,Field2,Field3,Field4,Field5)

END SUBROUTINE myWriteSol_Pres
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myWriteSol_Coor(iInd,iiLev,Field,Field1,Field2,Field3,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE myWriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE Transport_Q2P1,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE myWriteSol
end interface
REAL*8         Field(3,*),Field1(*),Field2(*),Field3(*)
INTEGER        i,iInd,nn,iiLev
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Coordinates'

IF (myid.ne.0) THEN
 DO i=1,nn
  Field1(i) = Field(1,i)
  Field2(i) = Field(2,i)
  Field3(i) = Field(3,i)
 END DO
END IF

CALL myWriteSol(iInd,iType,iiLev,cFF,Field1,Field2,Field3)

END SUBROUTINE myWriteSol_Coor
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myWriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
USE Transport_Q2P1,ONLY:myDump, ViscoSc
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
INTEGER iOut,iType,iiLev
REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
CHARACTER cField*(20)
!----------------------------------------------------
INTEGER i,ivt,jvt,jel,kel,iP,nLengthE,nLengthV
REAL*8,ALLOCATABLE :: Field(:,:),auxField(:,:)
REAL*8 dMaxNel
INTEGER pnel,pID
CHARACTER cFile*(40)

IF (myid.NE.0) THEN


 nLengthE = 8**((NLMAX+iiLev)-1)
 nLengthV = (2**((NLMAX+iiLev))+1)**3

!  WRITE(*,*)  'WRITE :::',nLengthE,nLengthV
END IF

 IF (myid.ne.0) dMaxNel = DBLE(KNEL(1))
 CALL COMM_Maximum(dMaxNel)

 IF (myid.eq.showid) THEN
  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_'//TRIM(ADJUSTL(cField))//'.dmp'
  WRITE(MTERM,*) 'Releasing current '//TRIM(ADJUSTL(cField))//' solution into: "'//ADJUSTL(TRIM(cFile)),'"'
 END IF

 IF (myid.eq.0) THEN
  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_'//TRIM(ADJUSTL(cField))//'.dmp'
  OPEN(321,FILE=TRIM(ADJUSTL(cFile)))
 END IF

 IF (iType.eq.1) THEN
  IF (Present(Field1)) CALL CollectVertField(Field1)
  IF (Present(Field2)) CALL CollectVertField(Field2)
  IF (Present(Field3)) CALL CollectVertField(Field3)
 END IF

 IF (iType.eq.2) THEN
  IF (Present(Field1)) CALL CollectElemField(Field1)
  IF (Present(Field2)) CALL CollectElemField(Field2)
  IF (Present(Field3)) CALL CollectElemField(Field3)
  IF (Present(Field4)) CALL CollectElemField(Field4)
 END IF

 IF (iType.eq.3) THEN
   CALL CollectVertField(ViscoSc%val11)
   CALL CollectVertField(ViscoSc%val22)
   CALL CollectVertField(ViscoSc%val33)
   CALL CollectVertField(ViscoSc%val12)
   CALL CollectVertField(ViscoSc%val13)
   CALL CollectVertField(ViscoSc%val23)
 END IF

 IF (myid.eq.0) THEN
  CLOSE(321)
 END IF

 CONTAINS
! -----------------------------------------------------------------
SUBROUTINE CollectVertField(xField)
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthV,0)
  CALL SENDI_myMPI(KNEL(1),0)
  
  ALLOCATE(Field(nLengthV,KNEL(1))) 

  DO iel = 1,KNEL(1)
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    Field(ivt,iel) = xField(jvt)
   END DO
  END DO
 
  CALL SENDD_myMPI(Field,nLengthV*KNEL(1),0)
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthV,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthV,KNEL(1))) 
    ALLOCATE(auxField(nLengthV,pnel)) 
   END IF
   CALL RECVI_myMPI(pnel,pID)

   CALL RECVD_myMPI(auxField,pnel*nLengthV,pID)

   DO I=1,pnel
   IEL = coarse%pELEMLINK(pID,I)
    DO ivt=1,nLengthV
     Field(ivt,iel) = auxField(ivt,I)
    END DO
   END DO
  END DO
 
  DEALLOCATE(auxField) 
!  WRITE(321,'(A)') '-- - ---'
  DO iel=1,KNEL(1)
   !WRITE(321,'(<nLengthV>ES14.6)') Field(1:nLengthV,iel)
   WRITE(321,*) Field(1:nLengthV,iel)
  END DO
 END IF

DEALLOCATE(Field) 

END SUBROUTINE CollectVertField
! -----------------------------------------------------------------
SUBROUTINE CollectElemField(xField)
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthE,0)
  CALL SENDI_myMPI(KNEL(1),0)
  
  ALLOCATE(Field(nLengthE,KNEL(1))) 

  DO iel = 1,KNEL(1)
   DO ivt=1,nLengthE
    jvt = myDump%Elements(IEL,ivt)
    Field(ivt,iel) = xField(jvt)
   END DO
  END DO
 
  CALL SENDD_myMPI(Field,nLengthE*KNEL(1),0)
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthE,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthE,KNEL(1))) 
    ALLOCATE(auxField(nLengthE,pnel)) 
   END IF
   CALL RECVI_myMPI(pnel,pID)

   CALL RECVD_myMPI(auxField,pnel*nLengthE,pID)

   DO I=1,pnel
   IEL = coarse%pELEMLINK(pID,I)
    DO ivt=1,nLengthE
     Field(ivt,iel) = auxField(ivt,I)
    END DO
   END DO
  END DO
 
  DEALLOCATE(auxField) 
!  WRITE(321,'(A)') '-- - ---'
  DO iel=1,KNEL(1)
   WRITE(321,*) Field(1:nLengthE,iel)
   !WRITE(321,'(<nLengthE>ES14.6)') Field(1:nLengthE,iel)
  END DO
 END IF

DEALLOCATE(Field) 

END SUBROUTINE CollectElemField
! -----------------------------------------------------------------

END SUBROUTINE myWriteSol
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myReadSol_Coor(cInFile,iLevel,Field,Field1,Field2,Field3,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE myReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE Transport_Q2P1,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE myReadSol
end interface
REAL*8         Field(3,*),Field1(*),Field2(*),Field3(*)
INTEGER        iLevel,i,nn
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Coordinates'
CHARACTER*(60) :: cInFile

CALL myReadSol(cInFile,iLevel,iType,cFF,Field1,Field2,Field3)

IF (myid.ne.0) THEN
!  WRITE(*,*) Field2(1:nn)
 DO i=1,nn
  Field(1,i) = Field1(i) 
  Field(2,i) = Field2(i) 
  Field(3,i) = Field3(i) 
 END DO
END IF

END SUBROUTINE myReadSol_Coor
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myReadSol_Velo(cInFile,iLevel,Field1,Field2,Field3)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE myReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE Transport_Q2P1,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE myReadSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*)
INTEGER        iLevel
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Velocity'
CHARACTER*(60) :: cInFile

CALL myReadSol(cInFile,iLevel,iType,cFF,Field1,Field2,Field3)

END SUBROUTINE myReadSol_Velo
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myReadSol_Visco(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE myReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE Transport_Q2P1,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE myReadSol
end interface
INTEGER        iLevel
INTEGER        :: iType= 3
CHARACTER*(20) :: cFF='Stress'
CHARACTER*(60) :: cInFile

CALL myReadSol(cInFile,iLevel,iType,cFF)

END SUBROUTINE myReadSol_Visco
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myReadSol_Pres(cInFile,iLevel,Field1,Field2,Field3,Field4,Field5,nn)
USE PP3D_MPI, ONLY:myid,showid,Comm_Summ
interface
  SUBROUTINE myReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE Transport_Q2P1,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE myReadSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*),Field4(*),Field5(*)
INTEGER        iLevel,i,nn
INTEGER        :: iType= 2
CHARACTER*(20) :: cFF='Pressure'
CHARACTER*(60) :: cInFile

CALL myReadSol(cInFile,iLevel,iType,cFF,Field2,Field3,Field4,Field5)

IF (myid.ne.0) THEN
!  WRITE(*,*) Field2(1:nn)
 DO i=1,nn
  Field1(4*(i-1)+1) = Field2(i) 
  Field1(4*(i-1)+2) = Field3(i) 
  Field1(4*(i-1)+3) = Field4(i) 
  Field1(4*(i-1)+4) = Field5(i) 
 END DO
END IF

END SUBROUTINE myReadSol_Pres
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myReadSol_Time(cInFile)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,subnodes,RECVDD_myMPI,SENDDD_myMPI

CHARACTER cInFile*(60)
INTEGER pID

 IF (myid.eq.showid) THEN
  WRITE(MTERM,*) 'Loading Time level from: "'//ADJUSTL(TRIM(cInFile))//'_Time.dmp','"'
 END IF

 IF (myid.eq.0) THEN
  OPEN(321,FILE=TRIM(ADJUSTL(cInFile))//'_Time.dmp')
  READ(321,*) TIMENS
  CLOSE(321)
 END IF

 IF (myid.NE.0) THEN
  CALL RECVDD_myMPI(timens,0)
 ELSE
  DO pID =1,subnodes
   CALL SENDDD_myMPI(timens,pID)
  END DO
 END IF

END SUBROUTINE myReadSol_Time
!
!-------------------------------------------------------------------------------
!
SUBROUTINE myReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
USE Transport_Q2P1,ONLY:myDump,ViscoSc
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
INTEGER iInd,iType,iLevel
REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
CHARACTER cField*(20),cInFile*(60)
!----------------------------------------------------
INTEGER i,ivt,jvt,jel,kel,iP,nLengthE,nLengthV
REAL*8,ALLOCATABLE :: Field(:,:),auxField(:,:)
REAL*8 dMaxNel
INTEGER pnel,pID
CHARACTER cFile*(40)

IF (myid.NE.0) THEN


 nLengthE = 8**(NLMAX+(iLevel-1)-1)
 nLengthV = (2**(NLMAX+(iLevel-1))+1)**3

!  WRITE(*,*)  'READ  :::',nLengthE,nLengthV
END IF

 IF (myid.ne.0) dMaxNel = DBLE(KNEL(1))
 CALL COMM_Maximum(dMaxNel)

 IF (myid.eq.showid) THEN
  WRITE(MTERM,*) 'Loading dumped '//TRIM(ADJUSTL(cField))//' solution from: "'//ADJUSTL(TRIM(cInFile))//'_'//TRIM(ADJUSTL(cField))//'.dmp','"'
 END IF

 IF (myid.eq.0) THEN
  OPEN(321,FILE=TRIM(ADJUSTL(cInFile))//'_'//TRIM(ADJUSTL(cField))//'.dmp')
 END IF

 IF (iType.eq.1) THEN
  IF (Present(Field1)) CALL DistributeVertField(Field1)
  IF (Present(Field2)) CALL DistributeVertField(Field2)
  IF (Present(Field3)) CALL DistributeVertField(Field3)
 END IF

 IF (iType.eq.2) THEN
  IF (Present(Field1)) CALL DistributeElemField(Field1)
  IF (Present(Field2)) CALL DistributeElemField(Field2)
  IF (Present(Field3)) CALL DistributeElemField(Field3)
  IF (Present(Field4)) CALL DistributeElemField(Field4)
 END IF

 IF (iType.eq.3) THEN
   CALL DistributeVertField(ViscoSc%val11)
   CALL DistributeVertField(ViscoSc%val22)
   CALL DistributeVertField(ViscoSc%val33)
   CALL DistributeVertField(ViscoSc%val12)
   CALL DistributeVertField(ViscoSc%val13)
   CALL DistributeVertField(ViscoSc%val23)
 END IF

 IF (myid.eq.0) THEN
  CLOSE(321)
 END IF

 CONTAINS
! -----------------------------------------------------------------
SUBROUTINE DistributeVertField(xField)
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthV,0)
  CALL SENDI_myMPI(KNEL(1),0)
  
  ALLOCATE(Field(nLengthV,KNEL(1))) 

  CALL RECVD_myMPI(Field,nLengthV*KNEL(1),0)

  DO iel = 1,KNEL(1)
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    xField(jvt) = Field(ivt,iel)
   END DO
  END DO
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthV,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthV,KNEL(1))) 
    ALLOCATE(auxField(nLengthV,pnel)) 

    DO iel=1,KNEL(1)
!    READ(321,'(<nLengthV>ES14.6)') Field(1:nLengthV,iel)
     READ(321,*) Field(1:nLengthV,iel)
    END DO
   END IF
   CALL RECVI_myMPI(pnel,pID)

   DO I=1,pnel
   IEL = coarse%pELEMLINK(pID,I)
    DO ivt=1,nLengthV
     auxField(ivt,I) = Field(ivt,iel)
    END DO
   END DO
 
   CALL SENDD_myMPI(auxField,pnel*nLengthV,pID)

  END DO

  DEALLOCATE(auxField) 

 END IF

DEALLOCATE(Field) 

END SUBROUTINE DistributeVertField
! -----------------------------------------------------------------
SUBROUTINE DistributeElemField(xField)
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthE,0)
  CALL SENDI_myMPI(KNEL(1),0)
  
  ALLOCATE(Field(nLengthE,KNEL(1))) 

  CALL RECVD_myMPI(Field,nLengthE*KNEL(1),0)

  DO iel = 1,KNEL(1)
   DO ivt=1,nLengthE
    jvt = myDump%Elements(IEL,ivt)
    xField(jvt) = Field(ivt,iel)
   END DO
  END DO
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthE,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthE,KNEL(1))) 
    ALLOCATE(auxField(nLengthE,pnel)) 
    DO iel=1,KNEL(1)
!    READ(321,'(<nLengthE>ES14.6)') Field(1:nLengthE,iel)
     READ(321,*) Field(1:nLengthE,iel)
    END DO
   END IF
   CALL RECVI_myMPI(pnel,pID)

   DO I=1,pnel
   IEL = coarse%pELEMLINK(pID,I)
    DO ivt=1,nLengthE
     auxField(ivt,I) = Field(ivt,iel)
    END DO
   END DO

   CALL SENDD_myMPI(auxField,pnel*nLengthE,pID)
  END DO
 
  DEALLOCATE(auxField) 

 END IF

DEALLOCATE(Field) 

END SUBROUTINE DistributeElemField
! -----------------------------------------------------------------

END SUBROUTINE myReadSol

