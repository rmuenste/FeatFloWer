SUBROUTINE SolToFile(iOutput)
USE def_FEAT
USE QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic
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


END SUBROUTINE SolToFile
!
! ----------------------------------------------
!
SUBROUTINE SolFromFile(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid
USE def_FEAT
USE QuadScalar,ONLY:QuadSc,LinSc,SetUp_myQ2Coor,bViscoElastic
USE LinScalar,ONLY:Tracer
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

CALL ReadSol_Velo(cInFile,iLevel,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)
nn = knel(nlmax)
CALL ReadSol_Pres(cInFile,iLevel,LinSc%ValP(NLMAX)%x,LinSc%AuxP(NLMAX)%x(1),&
     LinSc%AuxP(NLMAX)%x(nn+1),LinSc%AuxP(NLMAX)%x(2*nn+1),LinSc%AuxP(NLMAX)%x(3*nn+1),nn)
! CALL ReadSol_Coor(cInFile,iLevel,DWORK(L(KLCVG(NLMAX))),QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof)
CALL ReadSol_Time(cInFile)

if(bViscoElastic)then
  CALL ReadSol_Visco(cInFile, iLevel)
end if

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MESH eXchange on coarse level !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IF (myid.EQ.0) THEN
!  CALL CreateDumpStructures(0)
! ELSE
!  LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
!  CALL CreateDumpStructures(LevDif)
! END IF
! 
! ILEV = LinSc%prm%MGprmIn%MedLev
! CALL SETLEV(2)
! nLengthV = (2**(ILEV-1)+1)**3
! nLengthE = KNEL(NLMIN)
! ALLOCATE(SendVect(3,nLengthV,nLengthE))
! CALL SendNodeValuesToCoarse(SendVect,DWORK(L(LCORVG)),KWORK(L(LVERT)),nLengthV,nLengthE,NEL,NVT)
! DEALLOCATE(SendVect)
! 
! IF (myid.ne.0) CALL CreateDumpStructures(1)
! 
! ! ILEV = NLMIN
! ! CALL SETLEV(2)
! ! CALL ExchangeNodeValuesOnCoarseLevel(DWORK(L(LCORVG)),KWORK(L(LVERT)),NVT,NEL)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MESH eXchange on coarse level !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ILEV=NLMAX
! CALL SETLEV(2)
! CALL SetUp_myQ2Coor(DWORK(L(LCORVG)),DWORK(L(LCORAG)),&
!      KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)))
! 
! !!!!!!!!!!!!!!!!!!!!!! ALE !!!!!!!!!!!!!!!!!!!!!!
! ! CALL StoreOrigCoor(DWORK(L(KLCVG(NLMAX))))
! !!!!!!!!!!!!!!!!!!!!!! ALE !!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE SolFromFile
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Time(iOut)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,subnodes,RECVDD_myMPI,SENDDD_myMPI
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
end interface
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

END SUBROUTINE WriteSol_Time
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Velo(iInd,iiLev,Field1,Field2,Field3)
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*)
INTEGER        iInd,iiLev
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Velocity'

CALL WriteSol(iInd,iType,iiLev,cFF,Field1,Field2,Field3)

END SUBROUTINE WriteSol_Velo
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Visco(iInd,iiLev)
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
end interface
INTEGER        iInd,iiLev
INTEGER        :: iType= 3
CHARACTER*(20) :: cFF='Stress'

CALL WriteSol(iInd,iType,iiLev,cFF)

END SUBROUTINE WriteSol_Visco
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Pres(iInd,iiLev,Field1,Field2,Field3,Field4,Field5,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
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


CALL WriteSol(iInd,iType,iiLev,cFF,Field2,Field3,Field4,Field5)

END SUBROUTINE WriteSol_Pres
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Coor(iInd,iiLev,Field,Field1,Field2,Field3,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
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

CALL WriteSol(iInd,iType,iiLev,cFF,Field1,Field2,Field3)

END SUBROUTINE WriteSol_Coor

SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
USE QuadScalar,ONLY:myDump, ViscoSc
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

 ILEV = NLMIN

 nLengthE = 8**((NLMAX+iiLev)-1)
 nLengthV = (2**((NLMAX+iiLev))+1)**3

!  WRITE(*,*)  'WRITE :::',nLengthE,nLengthV
END IF

 IF (myid.ne.0) dMaxNel = DBLE(KNEL(NLMIN))
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
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthV,0)
  CALL SENDI_myMPI(KNEL(NLMIN),0)
  
  ALLOCATE(Field(nLengthV,KNEL(NLMIN))) 

  DO iel = 1,KNEL(NLMIN)
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    Field(ivt,iel) = xField(jvt)
   END DO
  END DO
 
  CALL SENDD_myMPI(Field,nLengthV*KNEL(NLMIN),0)
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthV,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthV,KNEL(NLMIN))) 
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
  DO iel=1,KNEL(NLMIN)
   !WRITE(321,'(<nLengthV>ES14.6)') Field(1:nLengthV,iel)
   WRITE(321,*) Field(1:nLengthV,iel)
  END DO
 END IF

DEALLOCATE(Field) 

END SUBROUTINE CollectVertField
! -----------------------------------------------------------------
SUBROUTINE CollectElemField(xField)
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthE,0)
  CALL SENDI_myMPI(KNEL(NLMIN),0)
  
  ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 

  DO iel = 1,KNEL(NLMIN)
   DO ivt=1,nLengthE
    jvt = myDump%Elements(IEL,ivt)
    Field(ivt,iel) = xField(jvt)
   END DO
  END DO
 
  CALL SENDD_myMPI(Field,nLengthE*KNEL(NLMIN),0)
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthE,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 
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
  DO iel=1,KNEL(NLMIN)
   WRITE(321,*) Field(1:nLengthE,iel)
   !WRITE(321,'(<nLengthE>ES14.6)') Field(1:nLengthE,iel)
  END DO
 END IF

DEALLOCATE(Field) 

END SUBROUTINE CollectElemField
! -----------------------------------------------------------------

END SUBROUTINE WriteSol
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Coor(cInFile,iLevel,Field,Field1,Field2,Field3,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
REAL*8         Field(3,*),Field1(*),Field2(*),Field3(*)
INTEGER        iLevel,i,nn
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Coordinates'
CHARACTER*(60) :: cInFile

CALL ReadSol(cInFile,iLevel,iType,cFF,Field1,Field2,Field3)

IF (myid.ne.0) THEN
!  WRITE(*,*) Field2(1:nn)
 DO i=1,nn
  Field(1,i) = Field1(i) 
  Field(2,i) = Field2(i) 
  Field(3,i) = Field3(i) 
 END DO
END IF

END SUBROUTINE ReadSol_Coor
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Velo(cInFile,iLevel,Field1,Field2,Field3)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*)
INTEGER        iLevel
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Velocity'
CHARACTER*(60) :: cInFile

CALL ReadSol(cInFile,iLevel,iType,cFF,Field1,Field2,Field3)

END SUBROUTINE ReadSol_Velo
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Visco(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
INTEGER        iLevel
INTEGER        :: iType= 3
CHARACTER*(20) :: cFF='Stress'
CHARACTER*(60) :: cInFile

CALL ReadSol(cInFile,iLevel,iType,cFF)

END SUBROUTINE ReadSol_Visco
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Pres(cInFile,iLevel,Field1,Field2,Field3,Field4,Field5,nn)
USE PP3D_MPI, ONLY:myid,showid,Comm_Summ
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*),Field4(*),Field5(*)
INTEGER        iLevel,i,nn
INTEGER        :: iType= 2
CHARACTER*(20) :: cFF='Pressure'
CHARACTER*(60) :: cInFile

CALL ReadSol(cInFile,iLevel,iType,cFF,Field2,Field3,Field4,Field5)

IF (myid.ne.0) THEN
!  WRITE(*,*) Field2(1:nn)
 DO i=1,nn
  Field1(4*(i-1)+1) = Field2(i) 
  Field1(4*(i-1)+2) = Field3(i) 
  Field1(4*(i-1)+3) = Field4(i) 
  Field1(4*(i-1)+4) = Field5(i) 
 END DO
END IF

END SUBROUTINE ReadSol_Pres
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Time(cInFile)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,subnodes,RECVDD_myMPI,SENDDD_myMPI
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
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

END SUBROUTINE ReadSol_Time
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
USE QuadScalar,ONLY:myDump,ViscoSc
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

 ILEV = NLMIN

 nLengthE = 8**(NLMAX+(iLevel-1)-1)
 nLengthV = (2**(NLMAX+(iLevel-1))+1)**3

!  WRITE(*,*)  'READ  :::',nLengthE,nLengthV
END IF

 IF (myid.ne.0) dMaxNel = DBLE(KNEL(NLMIN))
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
  CALL SENDI_myMPI(KNEL(NLMIN),0)
  
  ALLOCATE(Field(nLengthV,KNEL(NLMIN))) 

  CALL RECVD_myMPI(Field,nLengthV*KNEL(NLMIN),0)

  DO iel = 1,KNEL(NLMIN)
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
    ALLOCATE(Field(nLengthV,KNEL(NLMIN))) 
    ALLOCATE(auxField(nLengthV,pnel)) 

    DO iel=1,KNEL(NLMIN)
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
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthE,0)
  CALL SENDI_myMPI(KNEL(NLMIN),0)
  
  ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 

  CALL RECVD_myMPI(Field,nLengthE*KNEL(NLMIN),0)

  DO iel = 1,KNEL(NLMIN)
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
    ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 
    ALLOCATE(auxField(nLengthE,pnel)) 
    DO iel=1,KNEL(NLMIN)
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

END SUBROUTINE ReadSol
!
!-------------------------------------------------------------------------------
!
SUBROUTINE Output_Profiles(iOutput)
USE def_FEAT
USE QuadScalar,ONLY:QuadSc,LinSc,PressureToGMV,&
    Viscosity,Distance,Distamce,mgNormShearStress
USE LinScalar,ONLY:Tracer
USE PP3D_MPI, ONLY:myid,showid,Comm_Summ
USE var_QuadScalar,ONLY:myExport
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
  CALL Output_GMV_fields(iOutput,DWORK(L(KLCVG(ILEV))),KWORK(L(KLVERT(ILEV))))
  NLMAX = NLMAX - 1
 END IF

ELSEIF (myExport%Format.EQ."VTK") THEN

 IF (myid.NE.0) THEN
  NLMAX = NLMAX + 1
  ILEV = myExport%Level
  CALL SETLEV(2)
  CALL Output_VTK_piece(iOutput,DWORK(L(KLCVG(ILEV))),KWORK(L(KLVERT(ILEV))))
  NLMAX = NLMAX - 1
 ELSE
  CALL Output_VTK_main(iOutput)
 END IF

END IF


END
!
! ----------------------------------------------
!
SUBROUTINE Output_Profiles_sub(U,V,W,P,T,kv,dx,&
           iO,nvt,nel,time)
USE PP3D_MPI, ONLY:myid,showid
IMPLICIT NONE
REAL*8 dx(3,*),U(*),V(*),W(*),P(*),T(*),time
REAL*8 DGP,DAUX
INTEGER kv(8,*),nvt,nel,IP
CHARACTER cf*(20),cmm*(20),cm*(20)
INTEGER :: mf=120,ivt,iel
INTEGER iO


  WRITE(cf(13:15),'(A3)') "000"
  WRITE(cf(16:20),'(I1,A4)') iO,".gmv"


 cf="_gmv/               "

 IF(myid.lt.10) THEN
  WRITE(cf(6:12),'(A4,I1,I1,A1)') 'res_',0,myid,'_'
 ELSE
  WRITE(cf(6:12),'(A4,I2,A1)') 'res_',myid,'_'
 END IF

 IF     ((iO.GE.0   ).AND.(iO.LT.10   )) THEN
  WRITE(cf(13:15),'(A3)') "000"
  WRITE(cf(16:20),'(I1,A4)') iO,".gmv"
 ELSEIF ((iO.GE.10  ).AND.(iO.LT.100  )) THEN
  WRITE(cf(13:14),'(A2)') "00"
  WRITE(cf(15:20),'(I2,A4)') iO,".gmv"
 ELSEIF ((iO.GE.100 ).AND.(iO.LT.1000 )) THEN
  WRITE(cf(13:13),'(A1)') "0"
  WRITE(cf(14:20),'(I3,A4)') iO,".gmv"
 ELSEIF ((iO.GE.1000).AND.(iO.LT.10000)) THEN
  WRITE(cf(13:20),'(I4,A4)') iO,".gmv"
 ELSEIF (iO.GE.10000)                       THEN
  STOP
 END IF

  IF(myid.eq.showid) WRITE(*,*) "Outputting gmv file into ",cf

 cmm="msh                 "
 IF(myid.lt.10) WRITE(cmm(4:10),'(A1,I1,I1,A4)') '_',0,myid,".gmv"
 IF(myid.ge.10) WRITE(cmm(4:10),'(A1,I2,A4)') '_',myid,".gmv"

 cm="_gmv/msh             "
 IF(myid.lt.10) WRITE(cm(9:15),'(A1,I1,I1,A4)') '_',0,myid,".gmv"
 IF(myid.ge.10) WRITE(cm(9:15),'(A1,I2,A4)') '_',myid,".gmv"

 IF (iO.EQ.0) CALL XGMVMS(mf,cm,nel,nvt,kv,dx)

 OPEN (UNIT=mf,FILE=cf)
 WRITE(mf,'(A)')'gmvinput ascii'
 WRITE(UNIT=mf,FMT=*) 'nodes fromfile "',TRIM(cmm),'"'
 WRITE(UNIT=mf,FMT=*) 'cells fromfile "',TRIM(cmm),'"'

 WRITE(mf,*)  'velocity 1'
 DO IVT=1,NVT
  WRITE(mf,1000) REAL(U(IVT))
 END DO
 DO IVT=1,NVT
  WRITE(mf,1000) REAL(V(IVT))
 END DO
 DO IVT=1,NVT
  WRITE(mf,1000) REAL(W(IVT))
 END DO

 WRITE(mf,*)  'variable'
 WRITE(mf,*)  'pressure 1'
 DO IVT=1,NVT
  WRITE(mf,1000) REAL(P(IVT))
 END DO

 WRITE(mf,*)  'temperature 1'
 DO IVT=1,NVT
  WRITE(mf,1000) REAL(T(IVT))
 END DO

 WRITE(mf,*)  'endvars'
 WRITE(mf,*)  'probtime ',REAL(time)
 WRITE(mf,*)  'endgmv'
 CLOSE(mf)

1000  FORMAT(E12.5)
1100  FORMAT(8I8)

END
!
! ----------------------------------------------
!
SUBROUTINE Output_DUMPProfiles()
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE QuadScalar,ONLY:QuadSc,LinSc,bTracer
USE PLinScalar,ONLY:PLinLS
USE LinScalar,ONLY:Tracer
IMPLICIT NONE
CHARACTER COFile*15
INTEGER itwx,ifilen,i
DATA ifilen/0/

IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)

 ifilen=ifilen+1
 itwx=MOD(ifilen+insavn-1,insavn)+1
 COFile='_ns/QL        '
 IF (itwx.lt.10) WRITE(COFile(7:9),'(I1,I1,A1)') 0,itwx,'_'
 IF (itwx.ge.10) WRITE(COFile(7:9),'(I2,A1)') itwx,'_'
 IF (myid.eq.showid) &
 WRITE(*,*) "Outputting dump profiles into ",COFile
 IF     (myid.LT.10 ) THEN
  WRITE(COFile(10:12),'(A,I1)') "00",myid
 ELSEIF (myid.LT.100) THEN
  WRITE(COFile(10:12),'(A,I2)') "0",myid
 ELSE
  WRITE(COFile(10:12),'(I3)') myid
 END IF

 OPEN (UNIT=2,FILE=COFile)

 DO I=1,SIZE(QuadSc%ValU)
  WRITE(2,*) QuadSc%ValU(i)
 END DO
 DO I=1,SIZE(QuadSc%ValV)
  WRITE(2,*) QuadSc%ValV(i)
 END DO
 DO I=1,SIZE(QuadSc%ValW)
  WRITE(2,*) QuadSc%ValW(i)
 END DO

 DO I=1,SIZE(LinSc%ValP(NLMAX)%x)
  WRITE(2,*) LinSc%ValP(NLMAX)%x(i)
 END DO

 IF (bTracer) THEN
  DO I=1,SIZE(Tracer%Val(NLMAX+1)%x)
   WRITE(2,*) Tracer%Val(NLMAX+1)%x(i)
  END DO
 END IF
 WRITE(2,*) timens

 CLOSE(2)
END IF

END
!
! ----------------------------------------------
!
SUBROUTINE Load_DUMPProfiles(cdump)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,ShareValueD_myMPI
USE QuadScalar,ONLY:QuadSc,LinSc,bTracer
USE PLinScalar,ONLY:PLinLS
USE LinScalar,ONLY:Tracer
IMPLICIT NONE
CHARACTER CIFile*(30),cdump*(*)
INTEGER IFl,I
REAL*8 ddd
LOGICAL bExist
REAL*8 daux(1)

 CIFile=TRIM(ADJUSTL(cdump))
 IF (myid.ne.0) THEN
  IF (myid.eq.showid) &
  WRITE(*,*) "Loading initial dump profiles from ",CIFile

  INQUIRE(FILE=CIFile,EXIST=bExist)
  IF (.NOT.bExist) THEN
   WRITE(*,*) "File ",CIFile," could not be read ...",myid
   RETURN
  END IF

  OPEN (UNIT=1,FILE=CIFile) 

  DO I=1,SIZE(QuadSc%ValU)
   READ(1,*) QuadSc%ValU(i)
  END DO
  DO I=1,SIZE(QuadSc%ValV)
   READ(1,*) QuadSc%ValV(i)
  END DO
  DO I=1,SIZE(QuadSc%ValW)
   READ(1,*) QuadSc%ValW(i)
  END DO

  DO I=1,SIZE(LinSc%ValP(NLMAX)%x)
   READ(1,*) LinSc%ValP(NLMAX)%x(i)
  END DO

  IF (bTracer) THEN
   DO I=1,SIZE(Tracer%Val(NLMAX+1)%x)
    READ(1,*) Tracer%Val(NLMAX+1)%x(i)
   END DO
  END IF

  READ(1,*) timens
  CLOSE(1)
 END IF

 daux(1) = timens
 CALL ShareValueD_myMPI(daux,1,1)
 timens = daux(1)

END
!
! ----------------------------------------------
!
SUBROUTINE Load_LowDUMPProfiles(cdump)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,ShareValueD_myMPI
USE QuadScalar,ONLY:QuadSc,LinSc,bTracer
USE PLinScalar,ONLY:PLinLS
USE LinScalar,ONLY:Tracer
IMPLICIT NONE
CHARACTER CIFile*(30),cdump*(*)
INTEGER IFl,I,ndof
REAL*8 ddd
LOGICAL bExist
REAL*8 daux(1)

 CIFile=TRIM(ADJUSTL(cdump))
 IF (myid.ne.0) THEN
  ILEV = NLMAX-1
  ndof = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
  IF (myid.eq.showid) &
  WRITE(*,*) "Loading initial dump profiles from ",CIFile

  INQUIRE(FILE=CIFile,EXIST=bExist)
  IF (.NOT.bExist) THEN
   WRITE(*,*) "File ",CIFile," could not be read ...",myid
   RETURN
  END IF

  OPEN (UNIT=1,FILE=CIFile) 

  DO I=1,ndof
   READ(1,*) QuadSc%ValU(i)
  END DO
  DO I=1,ndof
   READ(1,*) QuadSc%ValV(i)
  END DO
  DO I=1,ndof
   READ(1,*) QuadSc%ValW(i)
  END DO

  DO I=1,SIZE(LinSc%ValP(ILEV)%x)
   READ(1,*) LinSc%ValP(ILEV)%x(i)
  END DO

  READ(1,*) timens
  CLOSE(1)
 END IF

 daux(1) = timens
 CALL ShareValueD_myMPI(daux,1,1)
 timens = daux(1)

END
!
! ----------------------------------------------
!
SUBROUTINE FBM_ToFile()
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myFBM
CHARACTER COFile*15
INTEGER itwx,ifilen,i
INTEGER iP
DATA ifilen/0/

IF (myid.eq.showid) THEN

 ifilen=ifilen+1
 itwx=MOD(ifilen+insavn-1,insavn)+1
 COFile='_ns/PT        '
 IF (itwx.lt.10) WRITE(COFile(7:9),'(I1,I1)') 0,itwx
 IF (itwx.ge.10) WRITE(COFile(7:9),'(I2)') itwx
 WRITE(*,*) "Outputting particle data into ",COFile

 OPEN (UNIT=2,FILE=COFile)

 WRITE(2,'(I10)') myFBM%nParticles
 DO iP = 1,myFBM%nParticles
  WRITE(2,'(A)')    myFBM%ParticleNew(iP)%cType
  WRITE(2,'(D12.4)')  myFBM%ParticleNew(iP)%density
  WRITE(2,'(D12.4)')  myFBM%ParticleNew(iP)%sizes(1)
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%Position
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%Velocity
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%Angle
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%AngularVelocity
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%Acceleration
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%FrameVelocity
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%ResistanceForce
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%TorqueForce
 END DO

 CLOSE(2)

END IF

END SUBROUTINE FBM_ToFile
!
! ----------------------------------------------
!
SUBROUTINE FBM_FromFile(CIFile)
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myFBM
 INTEGER iP
 CHARACTER*(*) CIFile

 IF (myid.eq.showid) &
 WRITE(*,*) "Loading particle data from ",CIFile

 OPEN(1,FILE=TRIM(ADJUSTL(CIFile)))

 READ(1,'(I10)') myFBM%nParticles
 ALLOCATE(myFBM%ParticleNew(myFBM%nParticles),myFBM%ParticleOld(myFBM%nParticles))
 ALLOCATE(myFBM%Force(6*myFBM%nParticles))
 DO iP = 1,myFBM%nParticles
  READ(1,'(A)')    myFBM%ParticleOld(iP)%cType
  myFBM%ParticleOld(iP)%cType = TRIM(ADJUSTL(myFBM%ParticleOld(iP)%cType))
  READ(1,'(D12.4)')  myFBM%ParticleOld(iP)%density
  READ(1,'(D12.4)')  myFBM%ParticleOld(iP)%sizes(1)
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%Position
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%Velocity
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%Angle
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%AngularVelocity
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%Acceleration
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%FrameVelocity
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%ResistanceForce
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%TorqueForce
 END DO

 myFBM%ParticleNew = myFBM%ParticleOld
 CLOSE(1)

END SUBROUTINE FBM_FromFile
!
!---------------------------------------------------------------------------
!
SUBROUTINE CreateDumpStructures(iLevel)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE QuadScalar,ONLY:myDump

IMPLICIT NONE

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))

INTEGER iLevel
INTEGER JEL,KEL,ivt,jvt,nLengthE,nLengthV
INTEGER iaux,jaux,jj,kv(8),II
LOGICAL ,ALLOCATABLE :: bGot(:)
! -------------- workspace -------------------

! IF (myid.NE.0) THEN

 NLMAX = NLMAX+iLevel

 ILEV = NLMIN
 NEL  = KNEL(ILEV)
 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3
 IF(ALLOCATED(myDump%Elements)) DEALLOCATE(myDump%Elements)
 IF(ALLOCATED(myDump%Vertices)) DEALLOCATE(myDump%Vertices)
 ALLOCATE(myDump%Elements(NEL,nLengthE))
 ALLOCATE(myDump%Vertices(NEL,nLengthV))

 DO IEL = 1,KNEL(NLMIN)

  myDump%Elements(IEL,1) = IEL
  iaux = 1
  DO II=1+1,NLMAX
   jaux = iaux
   DO jel=1,jaux
    kel = myDump%Elements(iel,jel)
    CALL Get8Elem(KWORK(L(KLADJ(II))),kv,kel)
    DO jj = 2,8
     iaux = iaux + 1
     myDump%Elements(iel,iaux) = kv(jj)
    END DO
   END DO
  END DO
 END DO 

 ALLOCATE(bGot(KNVT(NLMAX)))
 DO IEL = 1,KNEL(NLMIN)
  bGot = .FALSE.
  iaux = 0
  DO JEL = 1,nLengthE
   KEL = myDump%Elements(IEL,JEL)
   
   CALL getVert(KWORK(L(KLVERT(NLMAX))),kv,KEL)

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

END SUBROUTINE CreateDumpStructures
!
!---------------------------------------------------------------------------
!
SUBROUTINE ReleaseSmartDumpFiles(iOutO)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier
USE QuadScalar,ONLY:QuadSc,LinSc,myDump
USE PLinScalar,ONLY:PLinLS
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

 DO IEL = 1,KNEL(NLMIN)
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

END SUBROUTINE ReleaseSmartDumpFiles
!
!---------------------------------------------------------------------------
!
SUBROUTINE LoadSmartDumpFiles(cFldrInFile,iLevel)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse
USE QuadScalar,ONLY:QuadSc,LinSc,myDump
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

 DO IEL = 1,KNEL(NLMIN)
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

END SUBROUTINE LoadSmartDumpFiles
!
!---------------------------------------------------------------------------
!
SUBROUTINE LoadSmartAdaptedMeshFile(dcorvg,cMeshInFile,iLevel)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse
USE QuadScalar,ONLY:myDump
IMPLICIT NONE
REAL*8 dcorvg(3,*)
INTEGER iLevel
CHARACTER cInFile*99,cMeshInFile*(*)
INTEGER myCoarseElem,ivt,jvt,jel,kel,nLengthE,nLengthV
INTEGER iP,lStr

IF (myid.NE.0) THEN

 NLMAX = NLMAX+iLevel

 ILEV = NLMIN

 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3

 cInFile = ' '
 WRITE(cInFile(1:),'(A)') TRIM(ADJUSTL(cMeshInFile))
 lStr = LEN(TRIM(ADJUSTL(cInFile)))

 IF (myid.EQ.1) THEN
  WRITE(*,'(3A)') "Reading the ",TRIM(ADJUSTL(cInFile))," series for preadaptation of the mesh ..."
 END IF

 DO IEL = 1,KNEL(NLMIN)
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
   READ(547,*) dcorvg(:,jvt)
  END DO

  CLOSE(547)

 END DO 

 NLMAX = NLMAX-iLevel

ELSE
  OPEN (FILE=TRIM(ADJUSTL(cMeshInFile))//"/coarse.prf",UNIT = 547)
  DO ivt = 1, NVT
   READ(547,*) dcorvg(:,ivt)
  END DO
  CLOSE (547)
END IF

END SUBROUTINE LoadSmartAdaptedMeshFile


SUBROUTINE Output_VTK_piece(iO,dcoor,kvert)
USE def_FEAT
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE QuadScalar,ONLY: QuadSc,LinSc,Viscosity,Distance,Distamce,mgNormShearStress,myALE
USE QuadScalar,ONLY: MixerKnpr,FictKNPR,ViscoSc
USE LinScalar,ONLY:Tracer
USE var_QuadScalar,ONLY:myExport, Properties, bViscoElastic

IMPLICIT NONE
REAL*8 dcoor(3,*)
INTEGER kvert(8,*),iO,iUnit,ioffset,ive,ivt,iField,i
CHARACTER fileid*(5),filename*(27),procid*(3)
INTEGER NoOfElem,NoOfVert
REAL*8,ALLOCATABLE ::  tau(:,:)
REAL*8 psi(6)

NoOfElem = KNEL(ILEV)
NoOfVert = KNVT(ILEV)

filename=" "
WRITE(filename(1:),'(A,I5.5,A4)') '_vtk/res_node_***.',iO,".vtu"

IF(myid.eq.showid) WRITE(*,'(104("="))') 
IF(myid.eq.showid) WRITE(*,*) "Outputting vtk file into ",filename
WRITE(filename(15:17),'(I3.3)') myid

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",KNVT(ILEV),""" NumberOfCells=""",NoOfElem,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Velocity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Velocity",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(QuadSc%ValU(ivt)),REAL(QuadSc%ValV(ivt)),REAL(QuadSc%ValW(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Stress')
  if(bViscoElastic)then

  ALLOCATE(tau(6,NoOfVert))
  DO i=1,NoOfVert
   psi = [ViscoSc%Val11(i),ViscoSc%Val22(i),ViscoSc%Val33(i),&
          ViscoSc%Val12(i),ViscoSc%Val13(i),ViscoSc%Val23(i)]
   CALL ConvertPsiToTau(psi,tau(:,i))   
     tau(1,i) = (tau(1,i) - 1d0)/Properties%ViscoLambda
     tau(2,i) = (tau(2,i) - 1d0)/Properties%ViscoLambda
     tau(3,i) = (tau(3,i) - 1d0)/Properties%ViscoLambda
     tau(4,i) = (tau(4,i) - 0d0)/Properties%ViscoLambda
     tau(5,i) = (tau(5,i) - 0d0)/Properties%ViscoLambda
     tau(6,i) = (tau(6,i) - 0d0)/Properties%ViscoLambda
  END DO

  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Stress",""" NumberOfComponents=""6"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,6E16.7)')"        ",REAL(tau(1,ivt)),REAL(tau(2,ivt)),REAL(tau(3,ivt)),REAL(tau(4,ivt)),REAL(tau(5,ivt)),REAL(tau(6,ivt))
  end do

  DEALLOCATE(tau)
  write(iunit, *)"        </DataArray>"
  end if
 
 CASE('MeshVelo')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MeshVelocity",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(myALE%MeshVelo(1,ivt)),REAL(myALE%MeshVelo(2,ivt)),REAL(myALE%MeshVelo(3,ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Pressure_V')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_V",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(LinSc%Q2(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Temperature')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Temperature",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Tracer%Val(NLMAX)%x(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Mixer')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Mixer",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(FictKNPR(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Viscosity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Viscosity",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Viscosity(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Monitor')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Monitor",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(myALE%monitor(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 END SELECT 

END DO

write(iunit, '(A)')"    </PointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <CellData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Pressure_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_E",""" format=""ascii"">"
   do ivt=1,NoOfElem
    ive = 4*(ivt-1)+1
    write(iunit, '(A,E16.7)')"        ",REAL(LinSc%ValP(NLMAX-1)%x(ive))
   end do
   write(iunit, *)"        </DataArray>"
  END IF

 CASE('Viscosity_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Viscosity_E",""" format=""ascii"">"
   do ivt=1,NoOfElem
    write(iunit, '(A,E16.7)')"        ",REAL(Viscosity((nvt+net+nat+ivt)))
   end do
   write(iunit, *)"        </DataArray>"
  END IF

 END SELECT

END DO

write(iunit, '(A)')"    </CellData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do ivt=1,NoOfVert
 write(iunit,'(A10,3E16.7)')"          ",REAL(myALE%OldCoor(1,ivt)),REAL(myALE%OldCoor(2,ivt)),REAL(myALE%OldCoor(3,ivt))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",NoOfElem-1,""">"
do ive=1,NoOfElem   
 write(iunit, '(8I10)')kvert(1,ive)-1,kvert(2,ive)-1,kvert(3,ive)-1,kvert(4,ive)-1,&
                       kvert(5,ive)-1,kvert(6,ive)-1,kvert(7,ive)-1,kvert(8,ive)-1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*NoOfElem,""">"
ioffset=NoOfElem/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,NoOfElem
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=NoOfElem/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,NoOfElem
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_VTK_piece



SUBROUTINE Output_VTK_main(iO)
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE var_QuadScalar,ONLY:myExport,bViscoElastic
USE def_FEAT

IMPLICIT NONE
INTEGER iO,iproc,iField
INTEGER :: iMainUnit=555
CHARACTER mainname*(20) 
CHARACTER filename*(26)

! generate the file name
mainname=' '
WRITE(mainname(1:),'(A,I5.5,A5)') '_vtk/main.',iO,'.pvtu'

OPEN (UNIT=imainunit,FILE=mainname)

write(imainunit, *)"<VTKFile type=""PUnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(imainunit, *)"  <PUnstructuredGrid GhostLevel=""0"">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, '(A)')"    <PPointData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Velocity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Velocity",""" NumberOfComponents=""3""/>"
 CASE('Stress')
  if(bViscoElastic)then
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Stress",""" NumberOfComponents=""6""/>"
  end if
 CASE('MeshVelo')
  write(imainunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MeshVelocity",""" NumberOfComponents=""3""/>"
 CASE('Pressure_V')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Pressure_V","""/>"
 CASE('Temperature')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Temperature","""/>"
 CASE('Mixer')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Mixer","""/>"
 CASE('Viscosity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Viscosity","""/>"
 CASE('Monitor')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Monitor","""/>"

 END SELECT
END DO

write(imainunit, '(A)')"    </PPointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, '(A)')"    <PCellData>"
DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Pressure_E')
!  WRITE(*,*) myExport%Level,myExport%LevelMax,myExport%Level.EQ.myExport%LevelMax
  IF (myExport%Level.EQ.myExport%LevelMax) THEN
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Pressure_E","""/>"
  END IF

 CASE('Viscosity_E')
!  WRITE(*,*) myExport%Level,myExport%LevelMax,myExport%Level.EQ.myExport%LevelMax
  IF (myExport%Level.EQ.myExport%LevelMax) THEN
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Viscosity_E","""/>"
  END IF

 END SELECT
END DO
write(imainunit, '(A)')"    </PCellData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, *)"    <PPoints>"
write(imainunit, *)"      <PDataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3""/>"
write(imainunit, *)"    </PPoints>"

do iproc=1,subnodes
 filename=" "
 WRITE(filename(1:),'(A9,I3.3,A1,I5.5,A4)') 'res_node_',iproc,'.',iO,".vtu"
 write(imainunit, '(A,A,A)')"      <Piece Source=""",trim(adjustl(filename)),"""/>"  
end do
write(imainunit, *)"  </PUnstructuredGrid>"
write(imainunit, *)"  </VTKFile>"
close(imainunit)

END SUBROUTINE Output_VTK_main


SUBROUTINE Output_GMV_fields(iO,dcoor,kvert)
USE def_FEAT
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE QuadScalar,ONLY:QuadSc,LinSc,Viscosity,Distance,Distamce,mgNormShearStress,myALE
USE LinScalar,ONLY:Tracer
USE var_QuadScalar,ONLY:myExport

IMPLICIT NONE
REAL*8 dcoor(3,*)
INTEGER kvert(8,*),iO
INTEGER NoOfVert,NoOfElem
CHARACTER cf*(30),cmm*(30),cm*(30)
INTEGER i,j,iField
INTEGER :: iOutUnit=555

 NoOfElem = KNEL(ILEV)
 NoOfVert = KNVT(ILEV)
 cf  = ' '
 cm  = ' '
 cmm = ' '

 WRITE(cf(1:),'(A9,A2,A1,I4.4,A4)') '_gmv/res_','**','_',iO,".gmv"
 WRITE(cmm(1:),'(A4,I2.2,A4)') 'msh_',myid,".gmv"
 WRITE(cm (1:),'(A9,I2.2,A4)') '_gmv/msh_',myid,".gmv"
!  WRITE(cf(1:),'(A9,A3,A1,I5.5,A4)') '_gmv/res_','***','_',iO,".gmv"
!  WRITE(cmm(1:),'(A4,I3.3,A4)') 'msh_',myid,".gmv"
!  WRITE(cm (1:),'(A9,I3.3,A4)') '_gmv/msh_',myid,".gmv"

 IF(myid.eq.showid) WRITE(*,'(104("="))') 
 IF(myid.eq.showid) WRITE(*,*) "Outputting gmv file into ","[",ADJUSTL(TRIM(cf)),"]"
 WRITE(cf(10:11),'(I2.2)') myid
!  WRITE(cf(10:12),'(I3.3)') myid

 IF (iO.EQ.0) THEN
  OPEN (UNIT=iOutUnit,FILE=ADJUSTL(TRIM(cm)))

  WRITE(iOutUnit,'(A)')'gmvinput ascii'
  WRITE(iOutUnit,*)'nodes ',NoOfVert

  DO i=1,NoOfVert
   WRITE(iOutUnit,1200) REAL(dcoor(1,i))
  END DO
  DO i=1,NoOfVert
   WRITE(iOutUnit,1200) REAL(dcoor(2,i))
  END DO
  DO i=1,NoOfVert
   WRITE(iOutUnit,1200) REAL(dcoor(3,i))
  END DO

  WRITE(iOutUnit,*)'cells ',NoOfElem
  DO i=1,NoOfElem
   WRITE(iOutUnit,*)'hex 8'
   WRITE(iOutUnit,1300) (kvert(j,i),j=1,8)
  END DO

  WRITE(iOutUnit,*)  'endgmv'

  CLOSE  (iOutUnit)
 
 END IF
! CALL XGMVMS(iOutUnit,cm,NoOfElem,NoOfVert,kvert,dcoor)

 OPEN (UNIT=iOutUnit,FILE=cf)
 WRITE(iOutUnit,'(A)')'gmvinput ascii'
 WRITE(UNIT=iOutUnit,FMT=*) 'nodes fromfile "',ADJUSTL(TRIM(cmm)),'"'
 WRITE(UNIT=iOutUnit,FMT=*) 'cells fromfile "',ADJUSTL(TRIM(cmm)),'"'

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Velocity')
  WRITE(iOutUnit,*)  'velocity 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(QuadSc%ValU(i))
  END DO
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(QuadSc%ValV(i))
  END DO
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(QuadSc%ValW(i))
  END DO

 END SELECT
END DO
 
WRITE(iOutUnit,*)  'variable'

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Pressure_V')
  WRITE(iOutUnit,*)  'pressure_V 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(LinSc%Q2(i))
  END DO

 CASE('Temperature')
  WRITE(iOutUnit,*)  'temperature 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(Tracer%Val(NLMAX)%x(i))
  END DO

 CASE('Mixer')
  WRITE(iOutUnit,*)  'mixer 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(Distamce(i))
  END DO

 CASE('Monitor')
  WRITE(iOutUnit,*)  'monitor 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(myALE%Monitor(i))
  END DO

 CASE('Viscosity')
  WRITE(iOutUnit,*)  'viscosity 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(Viscosity(i))
  END DO

 CASE('Pressure_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   WRITE(iOutUnit,*)  'pressure_E 0'
   DO i=1,NoOfElem
    j = 4*(i-1) + 1
    WRITE(iOutUnit,1000) REAL(LinSc%ValP(ILEV)%x(j))
   END DO
  END IF

 END SELECT
END DO

 WRITE(iOutUnit,*)  'endvars'
 WRITE(iOutUnit,*)  'probtime ',REAL(timens)
 WRITE(iOutUnit,*)  'endgmv'
 CLOSE(iOutUnit)

1000  FORMAT(E12.5)
1100  FORMAT(8I8)
1200  FORMAT(E12.5)
1300  FORMAT(8I8)

END SUBROUTINE Output_GMV_fields

!----------------------------------------------

SUBROUTINE FBM_SetParticles
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myFBM
 INTEGER iP,iParticles,ipc
 REAL*8 velx,vely,velz,angx,angy,angz,angvelx,angvely,angvelz
 REAL*8 posx,posy,posz

  DO iP = 1,myFBM%nParticles
   ipc=iP-1

   call getpos(posx,posy,posz,ipc)
   call getvel(velx,vely,velz,ipc)
   call getangle(angx,angy,angz,ipc)
   call getangvel(angvelx,angvely,angvelz,ipc)

   myFBM%ParticleNew(iP)%Position(1)=posx
   myFBM%ParticleNew(iP)%Position(2)=posy
   myFBM%ParticleNew(iP)%Position(3)=posz

   myFBM%ParticleNew(iP)%Velocity(1)=velx
   myFBM%ParticleNew(iP)%Velocity(2)=vely
   myFBM%ParticleNew(iP)%Velocity(3)=velz 

   myFBM%ParticleNew(iP)%Angle(1)=angx
   myFBM%ParticleNew(iP)%Angle(2)=angy
   myFBM%ParticleNew(iP)%Angle(3)=angz

   myFBM%ParticleNew(iP)%AngularVelocity(1)=angvelx
   myFBM%ParticleNew(iP)%AngularVelocity(2)=angvely
   myFBM%ParticleNew(iP)%AngularVelocity(3)=angvelz

  END DO

!     IF(myid.eq.1) THEN
!       IP=1
!       write(*,*)'new'
!       WRITE(*,'(A19,3D12.4)') "ResistanceForce: ",&
!           myFBM%particleNew(IP)%ResistanceForce
!       WRITE(*,'(A19,3D12.4)') "TorqueForce: ",&
!           myFBM%particleNew(IP)%TorqueForce
!       WRITE(*,'(A19,3D12.4)') "Position: ",&
!           myFBM%particleNew(IP)%Position
!       WRITE(*,'(A19,3D12.4)') "Velocity: ",&
!           myFBM%particleNew(IP)%Velocity
!       WRITE(*,'(A19,3D12.4)') "Angle: ",&
!           myFBM%particleNew(IP)%Angle
!       WRITE(*,'(A19,3D12.4)') "AngularVelocity: ",&
!           myFBM%particleNew(IP)%AngularVelocity
!     END IF

END SUBROUTINE FBM_SetParticles

!----------------------------------------------

SUBROUTINE FBM_GetParticles
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myFBM
 INTEGER iP,iParticles,ipc
 REAL*8 velx,vely,velz,angx,angy,angz,angvelx,angvely,angvelz
 REAL*8 posx,posy,posz,drad,density
 REAL*8 forcex,forcey,forcez
 REAL*8 torquex,torquey,torquez

 call getnumparticles(iParticles)
 myFBM%nParticles=iParticles-1

 ALLOCATE(myFBM%ParticleNew(myFBM%nParticles),myFBM%ParticleOld(myFBM%nParticles))
 ALLOCATE(myFBM%iel_ug(myFBM%nParticles))

 myFBM%iel_ug=0
 
 ALLOCATE(myFBM%Force(6*myFBM%nParticles))

 if(myid.eq.1)write(*,*)'Number of particles: ',myFBM%nParticles

  DO iP = 1,myFBM%nParticles
   ipc=iP-1
   myFBM%ParticleOld(iP)%cType = 'Sphere'
   
   ! get the density
   call getdensity(density,ipc)
   myFBM%ParticleOld(iP)%density = density
   myFBM%ParticleNew(iP)%density = density
 
   ! get the radius
   call getradius(drad,ipc)
   myFBM%ParticleOld(iP)%sizes(1)=drad
   myFBM%ParticleNew(iP)%sizes(1)=drad

   call getpos(posx,posy,posz,ipc)
   call getvel(velx,vely,velz,ipc)
   call getangle(angx,angy,angz,ipc)
   call getangvel(angvelx,angvely,angvelz,ipc)

   call getforce(forcex,forcey,forcez,ipc)
   call gettorque(torquex,torquey,torquez,ipc)

   myFBM%ParticleOld(iP)%Position(1)=posx
   myFBM%ParticleOld(iP)%Position(2)=posy
   myFBM%ParticleOld(iP)%Position(3)=posz

   myFBM%ParticleNew(iP)%Position(1)=posx
   myFBM%ParticleNew(iP)%Position(2)=posy
   myFBM%ParticleNew(iP)%Position(3)=posz

   myFBM%ParticleOld(iP)%Velocity(1)=velx
   myFBM%ParticleOld(iP)%Velocity(2)=vely
   myFBM%ParticleOld(iP)%Velocity(3)=velz 

   myFBM%ParticleNew(iP)%Velocity(1)=velx
   myFBM%ParticleNew(iP)%Velocity(2)=vely
   myFBM%ParticleNew(iP)%Velocity(3)=velz 

   myFBM%ParticleOld(iP)%Angle(1)=angx
   myFBM%ParticleOld(iP)%Angle(2)=angy
   myFBM%ParticleOld(iP)%Angle(3)=angz

   myFBM%ParticleNew(iP)%Angle(1)=angx
   myFBM%ParticleNew(iP)%Angle(2)=angy
   myFBM%ParticleNew(iP)%Angle(3)=angz

   myFBM%ParticleOld(iP)%AngularVelocity(1)=angvelx
   myFBM%ParticleOld(iP)%AngularVelocity(2)=angvely
   myFBM%ParticleOld(iP)%AngularVelocity(3)=angvelz

   myFBM%ParticleNew(iP)%AngularVelocity(1)=angvelx
   myFBM%ParticleNew(iP)%AngularVelocity(2)=angvely
   myFBM%ParticleNew(iP)%AngularVelocity(3)=angvelz

   myFBM%ParticleOld(iP)%Acceleration(1)=0
   myFBM%ParticleOld(iP)%Acceleration(2)=0
   myFBM%ParticleOld(iP)%Acceleration(3)=0

   myFBM%ParticleNew(iP)%Acceleration(1)=0
   myFBM%ParticleNew(iP)%Acceleration(2)=0
   myFBM%ParticleNew(iP)%Acceleration(3)=0

   myFBM%ParticleOld(iP)%FrameVelocity(1)=0
   myFBM%ParticleOld(iP)%FrameVelocity(2)=0
   myFBM%ParticleOld(iP)%FrameVelocity(3)=0

   myFBM%ParticleOld(iP)%ResistanceForce(1)=forcex
   myFBM%ParticleOld(iP)%ResistanceForce(2)=forcey
   myFBM%ParticleOld(iP)%ResistanceForce(3)=forcez

   myFBM%ParticleNew(iP)%ResistanceForce(1)=forcex
   myFBM%ParticleNew(iP)%ResistanceForce(2)=forcey
   myFBM%ParticleNew(iP)%ResistanceForce(3)=forcez

   myFBM%ParticleOld(iP)%TorqueForce(1)=torquex
   myFBM%ParticleOld(iP)%TorqueForce(2)=torquey
   myFBM%ParticleOld(iP)%TorqueForce(3)=torquez

   myFBM%ParticleNew(iP)%TorqueForce(1)=torquex
   myFBM%ParticleNew(iP)%TorqueForce(2)=torquey
   myFBM%ParticleNew(iP)%TorqueForce(3)=torquez

  END DO

!     IF(myid.eq.1) THEN
!       IP=1
!       write(*,*)'new'
!       WRITE(*,'(A19,3D12.4)') "ResistanceForce: ",&
!           myFBM%particleNew(IP)%ResistanceForce
!       WRITE(*,'(A19,3D12.4)') "TorqueForce: ",&
!           myFBM%particleNew(IP)%TorqueForce
!       WRITE(*,'(A19,3D12.4)') "Position: ",&
!           myFBM%particleNew(IP)%Position
!       WRITE(*,'(A19,3D12.4)') "Velocity: ",&
!           myFBM%particleNew(IP)%Velocity
!       WRITE(*,'(A19,3D12.4)') "Angle: ",&
!           myFBM%particleNew(IP)%Angle
!       WRITE(*,'(A19,3D12.4)') "AngularVelocity: ",&
!           myFBM%particleNew(IP)%AngularVelocity
!       write(*,*)'old'
!       WRITE(*,'(A19,3D12.4)') "ResistanceForce: ",&
!           myFBM%ParticleOld(IP)%ResistanceForce
!       WRITE(*,'(A19,3D12.4)') "TorqueForce: ",&
!           myFBM%ParticleOld(IP)%TorqueForce
!       WRITE(*,'(A19,3D12.4)') "Position: ",&
!           myFBM%ParticleOld(IP)%Position
!       WRITE(*,'(A19,3D12.4)') "Velocity: ",&
!           myFBM%ParticleOld(IP)%Velocity
!       WRITE(*,'(A19,3D12.4)') "Angle: ",&
!           myFBM%ParticleOld(IP)%Angle
!       WRITE(*,'(A19,3D12.4)') "AngularVelocity: ",&
!           myFBM%ParticleOld(IP)%AngularVelocity
!     END IF

END SUBROUTINE FBM_GetParticles
!
!----------------------------------------------
!
SUBROUTINE FBM_ScatterParticles
USE PP3D_MPI, ONLY:myid,showid,SENDD_myMPI,RECVD_myMPI,subnodes
USE var_QuadScalar,ONLY:myFBM
IMPLICIT NONE
INTEGER iP,iParticles,ipc,pID
REAL*8 velx,vely,velz,angx,angy,angz,angvelx,angvely,angvelz
REAL*8 posx,posy,posz
real*8, allocatable :: particleArray(:)

if(myid.eq.0)then

! master copies particle data from gpu
  DO iP = 1,myFBM%nParticles
    ipc=iP-1
    call getpos(posx,posy,posz,ipc)
    call getvel(velx,vely,velz,ipc)
    call getangle(angx,angy,angz,ipc)
    call getangvel(angvelx,angvely,angvelz,ipc)

    myFBM%ParticleNew(iP)%Position(1)=posx
    myFBM%ParticleNew(iP)%Position(2)=posy
    myFBM%ParticleNew(iP)%Position(3)=posz
    myFBM%ParticleNew(iP)%Velocity(1)=velx
    myFBM%ParticleNew(iP)%Velocity(2)=vely
    myFBM%ParticleNew(iP)%Velocity(3)=velz 
  END DO

end if

allocate(particleArray(6*myFBM%nParticles))
particleArray=0d0
! master scatters data to non-master processes
IF(myid.eq.0)THEN

  ! copy particle date into an array for communication
  DO iP = 0,myFBM%nParticles-1
    particleArray(iP*6+1)=myFBM%ParticleNew(iP+1)%Position(1)
    particleArray(iP*6+2)=myFBM%ParticleNew(iP+1)%Position(2)
    particleArray(iP*6+3)=myFBM%ParticleNew(iP+1)%Position(3)
    particleArray(iP*6+4)=myFBM%ParticleNew(iP+1)%Velocity(1)
    particleArray(iP*6+5)=myFBM%ParticleNew(iP+1)%Velocity(2)
    particleArray(iP*6+6)=myFBM%ParticleNew(iP+1)%Velocity(3)
  END DO

  ! send to the other processes
  DO pID=1,subnodes
    CALL SENDD_myMPI(particleArray,myFBM%nParticles*6,pID)  
  END DO
ELSE
 ! non-master processes receive data from master
 CALL RECVD_myMPI(particleArray,myFBM%nParticles*6,0)
  DO iP = 0,myFBM%nParticles-1
    myFBM%ParticleNew(iP+1)%Position(1)=particleArray(iP*6+1)
    myFBM%ParticleNew(iP+1)%Position(2)=particleArray(iP*6+2)
    myFBM%ParticleNew(iP+1)%Position(3)=particleArray(iP*6+3)
    myFBM%ParticleNew(iP+1)%Velocity(1)=particleArray(iP*6+4)
    myFBM%ParticleNew(iP+1)%Velocity(2)=particleArray(iP*6+5)
    myFBM%ParticleNew(iP+1)%Velocity(3)=particleArray(iP*6+6)
  END DO 
END IF
! end scatter  

! executed on the non-master processes
! copy particle data to c++ lib
if(myid.ne.0)then
  DO iP = 1,myFBM%nParticles
    ipc=iP-1
    posx=myFBM%ParticleNew(iP)%Position(1)
    posy=myFBM%ParticleNew(iP)%Position(2)
    posz=myFBM%ParticleNew(iP)%Position(3)

    velx=myFBM%ParticleNew(iP)%Velocity(1)
    vely=myFBM%ParticleNew(iP)%Velocity(2)
    velz=myFBM%ParticleNew(iP)%Velocity(3)

    call setpositionid(posx,posy,posz,ipc)
    call setvelocityid(velx,vely,velz,ipc)
  END DO
end if

END SUBROUTINE FBM_ScatterParticles



