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
  CALL Output_Profiles(0)
  CALL ZTIME(myStat%t1)
  myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
END IF

IF(dout.LE.(timens+1e-10)) THEN

  iOGMV = NINT(timens/dtgmv)
  IF (itns.ne.1) THEN
    CALL ZTIME(myStat%t0)
    CALL Output_Profiles(iOGMV)
    CALL ZTIME(myStat%t1)
    myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
  END IF
  dout=dout+dtgmv

  ! Save intermediate solution to a dump file
  IF (insav.NE.0.AND.itns.NE.1) THEN
    IF (MOD(iOGMV,insav).EQ.0) THEN
      CALL ZTIME(myStat%t0)
      CALL SolToFile(-1)
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

real, intent(inout) :: dttt0
integer, intent(in) :: filehandle

integer :: ierr
integer :: terminal = 6
real :: time,time_passed


CALL ZTIME(time)

time_passed = time - dttt0
CALL myStatOut(time_passed,filehandle)

CALL myStatOut(time_passed,terminal)

! Save the final solution vector in unformatted form
!CALL mySolToFile(-1)
CALL myReleaseSmartDumpFiles(-1)
CALL myOutput_Profiles(0)

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

