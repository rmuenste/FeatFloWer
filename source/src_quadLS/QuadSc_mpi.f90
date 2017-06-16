#ifdef MUMPS_AVAIL
SUBROUTINE E013_Comm_Master(DCORVG,KVERT,KEDGE,KAREA,NVT,NET,NAT,NEL)

USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY :my_crs_e013_map,my_crs_e010_map

IMPLICIT NONE

REAL*8  DCORVG(3,*)
INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)
INTEGER NAT,NEL,NVT,NET,ndof

INTEGER NeighE(2,12),NeighA(4,6)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
REAL*8 PX,PY,PZ
REAL*8 P1X,P1Y,P1Z,P2X,P2Y,P2Z

INTEGER,ALLOCATABLE :: iCoor(:),iAux(:)
REAL*8, ALLOCATABLE :: dCoor(:,:),rCoor(:,:)

INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,jaux
INTEGER pndof,cndof,pID,myINDEX
LOGICAL bFound

ndof = NVT+NET+NAT+NEL

ALLOCATE (iAux(1:ndof))
ALLOCATE (dCoor(3,1:ndof))
ALLOCATE (iCoor(1:ndof))

DO i=1,nvt
 PX = dcorvg(1,I)
 PY = dcorvg(2,I)
 PZ = dcorvg(3,I)
 dCoor(:,i) = [PX,PY,PZ]
END DO

k=1
DO i=1,nel
 DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
   ivt1 = kvert(NeighE(1,j),i)
   ivt2 = kvert(NeighE(2,j),i)
   PX = 0.5d0*(dcorvg(1,ivt1) + dcorvg(1,ivt2))
   PY = 0.5d0*(dcorvg(2,ivt1) + dcorvg(2,ivt2))
   PZ = 0.5d0*(dcorvg(3,ivt1) + dcorvg(3,ivt2))
   dCoor(:,nvt+k) = [PX,PY,PZ]
   k = k + 1
  END IF
 END DO
END DO

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   PX = 0.25d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2)+dcorvg(1,ivt3)+dcorvg(1,ivt4))
   PY = 0.25d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2)+dcorvg(2,ivt3)+dcorvg(2,ivt4))
   PZ = 0.25d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2)+dcorvg(3,ivt3)+dcorvg(3,ivt4))
   dCoor(:,nvt+net+k) = [PX,PY,PZ]
   k = k + 1
  END IF
 END DO
END DO

DO i=1,nel
 PX = 0d0
 PY = 0d0
 PZ = 0d0
 DO j=1,8
  PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
  PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
  PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
 END DO
 dCoor(:,nvt+net+nat+i) = [PX,PY,PZ]
END DO

IF (myid.eq.0) THEN
 DO pID=1,subnodes
  CALL RECVI_myMPI(pNDOF,pID)
  CALL SENDI_myMPI(ndof,pID)
  CALL SENDD_myMPI(dCoor,3*ndof,pID)
 END DO
ELSE
 CALL SENDI_myMPI(ndof,0)
 CALL RECVI_myMPI(cNDOF,0)
 ALLOCATE (rCoor(3,cndof))
 CALL RECVD_myMPI(rCoor,3*cNdof,0)
END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

jAux = 0
IF (myid.NE.0) THEN
 DO I=1,ndof
  P1X = dCoor(1,I)
  P1Y = dCoor(2,I)
  P1Z = dCoor(3,I)
  bFound = .FALSE.
  DO J=1,cndof
   P2X = rCoor(1,J)
   P2Y = rCoor(2,J)
   P2Z = rCoor(3,J)
   IF ((ABS(P1X-P2X).LT.DEpsPrec).AND.&
       (ABS(P1Y-P2Y).LT.DEpsPrec).AND.&
       (ABS(P1Z-P2Z).LT.DEpsPrec)) THEN
     iCoor(I) = j
     jAux = jAux + 1
     bFound = .TRUE.
    END IF
   END DO
!    IF (.NOT.BFOUND) WRITE(*,*) 'shit!!', myid,i
  END DO
 END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

IF (myid.eq.0) THEN
 ALLOCATE (my_crs_e013_map(subnodes))
 DO pID=1,subnodes
  CALL RECVI_myMPI(pNDOF,pID)
  my_crs_e013_map(pid)%ndof = pNDOF
  ALLOCATE(my_crs_e013_map(pid)%ind(pNDOF))
  CALL RECVK_myMPI(my_crs_e013_map(pid)%ind,pNdof,pID)
 END DO
ELSE
 CALL SENDI_myMPI(ndof,0)
 CALL SENDK_myMPI(iCoor,ndof,0)
END IF

IF (myid.eq.0) THEN
 DO pID=1,subnodes
  jaux = 0
  DO i=1,my_crs_e013_map(pid)%ndof
   jaux = jaux + 3
  END DO
  ALLOCATE(my_crs_e013_map(pid)%indE(jaux))
  ALLOCATE(my_crs_e013_map(pid)%dBuffer(jaux))
  my_crs_e013_map(pid)%cc_ndof = jaux

  jaux = 0
  DO i=1,my_crs_e013_map(pid)%ndof
   myINDEX = my_crs_e013_map(pid)%ind(i)
   my_crs_e013_map(pid)%indE(0*my_crs_e013_map(pid)%ndof + i) = 0*ndof + myINDEX
   my_crs_e013_map(pid)%indE(1*my_crs_e013_map(pid)%ndof + i) = 1*ndof + myINDEX
   my_crs_e013_map(pid)%indE(2*my_crs_e013_map(pid)%ndof + i) = 2*ndof + myINDEX
   IF (myINDEX.gt.NVT+NET+NAT) THEN
    jaux = jaux + 1
    my_crs_e013_map(pid)%indE(3*my_crs_e013_map(pid)%ndof + 4*(jaux-1)+1) = 3*ndof + 4*(myINDEX-(NVT+NET+NAT)-1)+1
    my_crs_e013_map(pid)%indE(3*my_crs_e013_map(pid)%ndof + 4*(jaux-1)+2) = 3*ndof + 4*(myINDEX-(NVT+NET+NAT)-1)+2
    my_crs_e013_map(pid)%indE(3*my_crs_e013_map(pid)%ndof + 4*(jaux-1)+3) = 3*ndof + 4*(myINDEX-(NVT+NET+NAT)-1)+3
    my_crs_e013_map(pid)%indE(3*my_crs_e013_map(pid)%ndof + 4*(jaux-1)+4) = 3*ndof + 4*(myINDEX-(NVT+NET+NAT)-1)+4
   END IF
  END DO

 END DO
END IF


END
!
! ----------------------------------------------
!
! SUBROUTINE Create_GlobalNumbering
! USE PP3D_MPI
! USE def_feat
! USE var_QuadScalar, ONLY : GlobalNumbering,my_crs_e013_map
! IMPLICIT NONE
! INTEGER i,j,ndof,nnQ2,nnP0,nQ2,nP0,iQ2,iP0,pID
! REAL*8, ALLOCATABLE :: D(:)
! REAL*8 daux
! 
! ILEV = NLMIN
! CALL SETLEV(2)
! 
! nnQ2  = NEL+NVT+NAT+NET
! nnP0 =  NEL
! ndof = 3*nnQ2 + 4*nnP0
! 
! D = 0
! 
! nQ2=0
! nP0=0
! j = 0
! 
! IF (myid.eq.0) THEN
! 
!  DO pID=1,subnodes
!   DO i=1,my_crs_e013_map(pid)%cc_ndof
!    j = my_crs_e013_map(pid)%indE(i)
!    my_crs_e013_map(pid)%dBuffer(i) = dble(j)
!   END DO
! !   write(*,*) myid, '<--', my_crs_e013_map(pid)%cc_ndof
!   CALL sendD_myMPI(my_crs_e013_map(pid)%dBuffer,my_crs_e013_map(pid)%cc_ndof,pID)
!  END DO
! 
! ELSE
! 
!   ALLOCATE (GlobalNumbering(ndof),D(ndof))
! !   write(*,*) myid, '-->', ndof  
!   CALL recvD_myMPI(D,ndof,0)
!   GlobalNumbering = INT(D)
!   DEALLOCATE (D)
!   WRITE(*,*) GlobalNumbering
! END IF
! 
! pause
! END SUBROUTINE Create_GlobalNumbering
!
! ----------------------------------------------
!
#endif
SUBROUTINE  Comm_Coor(dc,d1,d2,d3,d4,n)
USE PP3D_MPI
IMPLICIT NONE
REAL*8 dc(3,*),d1(*),d2(*),d3(*),d4(*)
INTEGER i,n

DO i=1,n
 d1(i) = dc(1,i)
 d2(i) = dc(2,i)
 d3(i) = dc(3,i)
 d4(i) = 1d0
END DO

CALL E013Sum(d1)
CALL E013Sum(d2)
CALL E013Sum(d3)
CALL E013Sum(d4)

DO i=1,n
 dc(1,i) = d1(i)/d4(i)
 dc(2,i) = d2(i)/d4(i)
 dc(3,i) = d3(i)/d4(i)
END DO

END SUBROUTINE  Comm_Coor
!
! ----------------------------------------------
!
SUBROUTINE CommunicateSurface()
USE PP3D_MPI
USE var_QuadScalar, ONLY : myTSurf,Properties
!-------------------------------------------------------
!-------------------------------------------------------
IMPLICIT NONE
INTEGER i,j,k,kk,pID,iaux,jaux
REAL*8, ALLOCATABLE :: dpack(:)
INTEGER iIF

do iIF=1,Properties%nInterface

IF (myid.ne.0) THEN
 CALL SENDI_myMPI(myTSurf(iIF)%nT,0)
 ALLOCATE(dpack(myTSurf(iIF)%nT*9*3))
ELSE
 myTSurf(iIF)%nT = 0
 DO pID=1,subnodes
  CALL RECVI_myMPI(iaux,pID)
  myTSurf(iIF)%nT = myTSurf(iIF)%nT + iaux
 END DO
 IF (ALLOCATED(myTSurf(iIF)%T)) DEALLOCATE (myTSurf(iIF)%T)
 ALLOCATE (myTSurf(iIF)%T(myTSurf(iIF)%nT))
 ALLOCATE (dpack(myTSurf(iIF)%nT*9*3))
END IF

!IF (myid.eq.0) WRITE(*,*) "nT = ", myTSurf(iIF)%nT
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

IF (myid.ne.0) THEN
 CALL SENDI_myMPI(myTSurf(iIF)%nT,0)
 kk = 1
 DO i=1,myTSurf(iIF)%nT
  DO j=1,9
   DO k=1,3
    dpack(kk) = myTSurf(iIF)%T(i)%C(k,j)
    kk = kk + 1
   END DO
  END DO
 END DO
 CALL SENDD_myMPI(dpack,myTSurf(iIF)%nT*9*3,0)
ELSE
 jaux = 0
 DO pID=1,subnodes
  CALL RECVI_myMPI(iaux,pID)
  CALL RECVD_myMPI(dpack,iaux*9*3,pID)
  kk = 1
  DO i=jaux+1,jaux + iaux
   DO j=1,9
    DO k=1,3
     myTSurf(iIF)%T(i)%C(k,j) = dpack(kk) 
     kk = kk + 1
    END DO
   END DO
  END DO
  jaux = jaux + iaux
!  myTSurf(iIF)%nT = myTSurf(iIF)%nT + iaux
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

IF (myid.ne.0) THEN

 DEALLOCATE (dpack,myTSurf(iIF)%T)
 CALL RECVI_myMPI(myTSurf(iIF)%nT,0)
 ALLOCATE(myTSurf(iIF)%T(myTSurf(iIF)%nT),dpack(myTSurf(iIF)%nT*9*3))
 CALL RECVD_myMPI(dpack,myTSurf(iIF)%nT*9*3,0)
 
 kk = 1
 DO i=1,myTSurf(iIF)%nT
  DO j=1,9
   DO k=1,3
    myTSurf(iIF)%T(i)%C(k,j) = dpack(kk)
    kk = kk + 1
   END DO
  END DO
 END DO

ELSE

 kk = 1
 DO i=1,myTSurf(iIF)%nT
  DO j=1,9
   DO k=1,3
    dpack(kk) = myTSurf(iIF)%T(i)%C(k,j)
    kk = kk + 1
   END DO
  END DO
 END DO

 DO pID=1,subnodes
  CALL SENDI_myMPI(myTSurf(iIF)%nT,pID)
  CALL SENDD_myMPI(dpack,myTSurf(iIF)%nT*9*3,pID)
 END DO

END IF

DEALLOCATE (dpack)

end do 

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)


  
! IF (myid.eq.1) THEN
!   WRITE(*,'(A,I,A)') "0 ", myTSurf(iIF)%nT*5, "1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"
!  DO i=1,myTSurf(iIF)%nT
!   DO j=1,5
!    WRITE(*,'(3ES12.4)') myTSurf(iIF)%T(i)%C(:,j)
!   END DO
! !    WRITE(*,'(5ES12.4)') myTSurf(iIF)%T(i)%C(1,:)
! !    WRITE(*,'(5ES12.4)') myTSurf(iIF)%T(i)%C(2,:)
! !    WRITE(*,'(5ES12.4)') myTSurf(iIF)%T(i)%C(3,:)
! !    WRITE(*,*) " - - - -  - - - - - - - - - - - - - - - - - - - - - -"
!  END DO
! END IF
! pause

  END SUBROUTINE CommunicateSurface
!
!-------------------------------------------------------  
!
SUBROUTINE ParPresComm_Init(KCOL,KLD,NEQ,NEL,ILEV)
USE PP3D_MPI
!-------------------------------------------------------
!-------------------------------------------------------
IMPLICIT NONE
INTEGER KCOL(*),KLD(*),NEQ,NEL,ILEV
INTEGER pID,pJD,I,J,IDOF,IA,IB,IEL,NNEL,MNEL
INTEGER ELEMS(NEL)

ALLOCATE(MGE013(ILEV)%SP(subnodes))

! Get the number of elements that are needed from the other side 
DO pID=1,subnodes
 MGE013(ILEV)%SP(pID)%nElems = (/0, 0/)
 MGE013(ILEV)%SP(pID)%nEntries = (/0, 0/)
 MGE013(ILEV)%SP(pID)%Num = 0 
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  ELEMS = 0
  DO I=1,MGE013(ILEV)%ST(pID)%Num
   IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
   DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
     IEL = (KCOL(IA)-1)/4+1
     ELEMS(IEL) = 1
   END DO
  END DO
  NNEL = 0
  DO IEL=1,NEL
   IF (ELEMS(IEL).NE.0) NNEL=NNEL+1
  END DO
  MGE013(ILEV)%SP(pID)%nElems(1) = NNEL
  !write(*,*) myid, "nElems=", NNEL
 END IF
END DO

! Send them!
DO pID=1,subnodes
 IF (myid.eq.pID) THEN
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    CALL RECVI_myMPI(NNEL,pJD)
    !write(*,*) "r:",myid,PJD,NNEL    
    MGE013(ILEV)%SP(pJD)%nElems(2) = NNEL
   END IF
  END DO
 ELSE
  NNEL = MGE013(ILEV)%SP(pID)%nElems(1)
  !write(*,*) "s:",myid,PID,NNEL  
  CALL SENDI_myMPI(NNEL ,pID)

 END IF
END DO

! Get the send-list of parallel elements
DO pID=1,subnodes
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  ELEMS = 0
  DO I=1,MGE013(ILEV)%ST(pID)%Num
   IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
   DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
     IEL = (KCOL(IA)-1)/4+1
     ELEMS(IEL) = 1
   END DO
  END DO
  NNEL = 0
  DO IEL=1,NEL
   IF (ELEMS(IEL).NE.0) NNEL=NNEL+1
  END DO
  MGE013(ILEV)%SP(pID)%Num = MAX(MGE013(ILEV)%SP(pID)%nElems(1),MGE013(ILEV)%SP(pID)%nElems(2))
  ALLOCATE(MGE013(ILEV)%SP(pID)%VertLink(2,MGE013(ILEV)%SP(pID)%Num))
  MGE013(ILEV)%SP(pID)%VertLink(:,:)=0
  ALLOCATE(MGE013(ILEV)%SP(pID)%SDVect(4*MGE013(ILEV)%SP(pID)%Num))
  MGE013(ILEV)%SP(pID)%SDVect=0d0
  ALLOCATE(MGE013(ILEV)%SP(pID)%RDVect(4*MGE013(ILEV)%SP(pID)%Num))
  MGE013(ILEV)%SP(pID)%RDVect=0d0
  NNEL = 0
  DO IEL=1,NEL
   IF (ELEMS(IEL).NE.0) THEN
    NNEL = NNEL + 1
    MGE013(ILEV)%SP(pID)%VertLink(1,NNEL) = IEL
   END IF
  END DO
 END IF
END DO

! Get the receive-list of parallel elements
MNEL = 0
DO pID=1,subnodes
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  NNEL = 0
  DO I=1,MGE013(ILEV)%SP(pID)%nElems(2)
   NNEL = NNEL + 1
   MNEL = MNEL + 1
   MGE013(ILEV)%SP(pID)%VertLink(2,NNEL) = MNEL
   !WRITE(*,*) MGE013(ILEV)%SP(pID)%VertLink(1:2,1)
  END DO
 END IF
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)


END
!
!
!
SUBROUTINE Create_ParB_LD_b(KLD_B,KLD,NEQ,ILEV)
USE PP3D_MPI
IMPLICIT NONE
!-------------------------------------------------------
!-------------------------------------------------------
! TYPE LocComm
!  INTEGER, ALLOCATABLE :: a(:)
! END TYPE LocComm
! TYPE(LocComm), ALLOCATABLE :: EntNumR(:),EntNumS(:)

INTEGER KLD(*),KLD_B(*),NEQ,ILEV
INTEGER, ALLOCATABLE :: EntNum(:)
INTEGER pID,pJD,I,J,IDOF,IA,IB,IEL,NNEL
INTEGER iLoc_1,iLoc_2
CHARACTER*11 myFile

WRITE(myFile(1:11),'(A4,I1,A1,I1,A4)') "Parb",myid,"_",ilev,".txt"
OPEN(789,FILE=myFile)

! ALLOCATE(EntNumS(subnodes),EntNumR(subnodes))

! Get the length of parallel KCOL and send it over!
DO pID=1,subnodes
 IF (pID.NE.myid) THEN
  IF (MGE013(ILEV)%ST(pID)%Num.GT.0) THEN
   NNEL  = 0
   ALLOCATE (EntNum(MGE013(ILEV)%ST(pID)%Num))
   EntNum = 0
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
    DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
      EntNum(I) = EntNum(I) + 1
      NNEL = NNEL + 1
    END DO
   END DO

!    MGE013(ILEV)%SP(pID)%nEntries(1) = NNEL
!    CALL SENDI_myMPI(NNEL ,pID)
!    NNEL = MGE013(ILEV)%ST(pID)%Num
     WRITE(789,*) "to - ",pID,NNEL
     WRITE(789,*) EntNum
!    CALL SENDK_myMPI(EntNum,NNEL,pID)

   DEALLOCATE(EntNum)

  END IF
 ELSE
!   DO pJD=1,subnodes
!    IF (myid.NE.pJD) THEN
!     IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
! 
!      CALL RECVI_myMPI(NNEL,pJD)
!      MGE013(ILEV)%SP(pJD)%nEntries(2) = NNEL
!      NNEL = MGE013(ILEV)%ST(pJD)%Num
! !      ALLOCATE (EntNumR(pJD)%a(NNEL))
!      ALLOCATE (EntNum(NNEL))
! !      CALL RECVK_myMPI(EntNumR(pJD)%a,NNEL,pJD)
!      CALL RECVK_myMPI(EntNum,NNEL,pJD)
! !      WRITE(789,*) "from - ",pJD
! !      WRITE(789,*) EntNum
!      DO I=1,MGE013(ILEV)%ST(pJD)%Num
!       IDOF=MGE013(ILEV)%ST(pJD)%VertLink(2,I)
! !       KLD_B(IDOF) = KLD_B(IDOF) + EntNumR(pJD)%a(I)
!       KLD_B(IDOF) = KLD_B(IDOF) + EntNum(I)
!      END DO
! 
!      DEALLOCATE(EntNum)
! 
!     END IF
!    END IF
!   END DO
 END IF
END DO

! iLoc_1 = KLD_B(1)
! KLD_B(1) = 1
! DO IDOF=2,NEQ+1
!  iLoc_2 = KLD_B(IDOF)
!  KLD_B(IDOF) = KLD_B(IDOF-1) + 4*iLoc_1
!  iLoc_1 = iLoc_2
! END DO

! WRITE(789,*) NEQ
! WRITE(789,*) 
! 
! DO IDOF=1,NEQ+1
!  WRITE(789,*) KLD_B(IDOF)
! END DO

CLOSE(789)

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
! DEALLOCATE(EntNumS,EntNumR)

END SUBROUTINE Create_ParB_LD_b
!
!
!
SUBROUTINE Create_ParB_LD(KLD_B,KLD,NEQ,ILEV)
USE PP3D_MPI
IMPLICIT NONE
!-------------------------------------------------------
!-------------------------------------------------------
TYPE LocComm
 INTEGER, ALLOCATABLE :: a(:)
END TYPE LocComm
TYPE(LocComm), ALLOCATABLE :: EntNumR(:),EntNumS(:)

INTEGER KLD(*),KLD_B(*),NEQ,ILEV
INTEGER pID,pJD,I,J,IDOF,IA,IB,IEL,NNEL,nnS
INTEGER iLoc_1,iLoc_2
CHARACTER*11 myFile

! WRITE(myFile(1:11),'(A4,I1,A1,I1,A4)') "ParB",myid,"_",ilev,".txt"
! OPEN(789,FILE=myFile)

ALLOCATE(EntNumS(subnodes),EntNumR(subnodes))

! Get the length of parallel KCOL and send it over!
DO pID=1,subnodes
 IF (pID.NE.myid) THEN
  IF (MGE013(ILEV)%ST(pID)%Num.GT.0) THEN
   NNEL  = 0
   ALLOCATE (EntNumS(pID)%a(MGE013(ILEV)%ST(pID)%Num))
   EntNumS(pID)%a(:) = 0
!    IF (myid.eq.1) WRITE(*,*) "To -",pID
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
    nnS = 0
    DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
      nnS = nnS + 1
    END DO
    NNEL = NNEL + nnS
    EntNumS(pID)%a(I) = EntNumS(pID)%a(I) + nnS
!     IF (myid.eq.1) WRITE(*,*) i,EntNumS(pID)%a(I),NNEL
   END DO

   MGE013(ILEV)%SP(pID)%nEntries(1) = NNEL
   CALL SENDI_myMPI(NNEL ,pID)
   NNEL = MGE013(ILEV)%ST(pID)%Num
!      WRITE(789,*) "to - ",pID,NNEL
!      WRITE(789,*) EntNumS(pID)%a
   CALL SENDK_myMPI(EntNumS(pID)%a,NNEL,pID)

  END IF
 ELSE
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     CALL RECVI_myMPI(NNEL,pJD)
     MGE013(ILEV)%SP(pJD)%nEntries(2) = NNEL
     NNEL = MGE013(ILEV)%ST(pJD)%Num
     ALLOCATE (EntNumR(pJD)%a(NNEL))
     EntNumR(pJD)%a=0
     CALL RECVK_myMPI(EntNumR(pJD)%a,NNEL,pJD)
!      WRITE(789,*) "from - ",pJD
!      WRITE(789,*) EntNumR(pJD)%a
    END IF
   END IF
  END DO
 END IF
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

DO pJD=1,subnodes
 IF (myid.NE.pJD) THEN
  IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
   DO I=1,MGE013(ILEV)%ST(pJD)%Num
    IDOF=MGE013(ILEV)%ST(pJD)%VertLink(2,I)
    KLD_B(IDOF) = KLD_B(IDOF) + EntNumR(pJD)%a(I)
   END DO
  END IF
 END IF
END DO

iLoc_1 = KLD_B(1)
KLD_B(1) = 1
DO IDOF=2,NEQ+1
 iLoc_2 = KLD_B(IDOF)
 KLD_B(IDOF) = KLD_B(IDOF-1) + 4*iLoc_1
 iLoc_1 = iLoc_2
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
DEALLOCATE(EntNumS,EntNumR)

END SUBROUTINE Create_ParB_LD
!
!
!
! The hell knows why but this routine is highy compilaer and amchine dependent ...
! I tried to adjust that more compilers like it but can be a problematic part ...

SUBROUTINE Create_ParB_COLMAT(BXMAT,BYMAT,BZMAT,BXPMAT,BYPMAT,BZPMAT,&
 KLD_B,KCOL_B,KLD_T,KLD,KCOL,NEQ,NEL,pNEL,ILEV)
USE PP3D_MPI
IMPLICIT NONE
!-------------------------------------------------------
!-------------------------------------------------------
REAL*8 BXMAT(*),BYMAT(*),BZMAT(*)
REAL*8 BXPMAT(*),BYPMAT(*),BZPMAT(*)
INTEGER KCOL(*),KLD(*),KCOL_B(*),KLD_B(*),KLD_T(*),NEQ,NEL,pNEL,ILEV
INTEGER pID,pJD,I,J,IDOF,IA,IB,IEL,NNEL
INTEGER ELEMS(NEL)
INTEGER, ALLOCATABLE :: KCOL_E(:),KLD_E(:)
REAL*8, ALLOCATABLE ::  MAT_EX(:),MAT_EY(:),MAT_EZ(:)
CHARACTER*10 myFile

! Fill up the entries of parallel KCOL
DO pID=1,subnodes
 IF (pID.NE.myid) THEN
  IF (MGE013(ILEV)%ST(pID)%Num.GT.0) THEN
   ELEMS = 0
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
    DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
      IEL = (KCOL(IA)-1)/4+1
      ELEMS(IEL) = 1
    END DO
   END DO

   ALLOCATE (KCOL_E(MGE013(ILEV)%SP(pID)%nEntries(1)))
   KCOL_E=0
   ALLOCATE (KLD_E(MGE013(ILEV)%ST(pID)%Num+1))
   KLD_E=0
   ALLOCATE (MAT_EX(4*MGE013(ILEV)%SP(pID)%nEntries(1)))
   MAT_EX=0
   ALLOCATE (MAT_EY(4*MGE013(ILEV)%SP(pID)%nEntries(1)))
   MAT_EY=0
   ALLOCATE (MAT_EZ(4*MGE013(ILEV)%SP(pID)%nEntries(1)))
   MAT_EZ=0

   DO IEL=2,NEL
    ELEMS(IEL) = ELEMS(IEL) + ELEMS(IEL-1)
   END DO

   KLD_E(1)=1
   IB = 0
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
    KLD_E(I+1) = KLD_E(I)
    DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
      IEL = (KCOL(IA)-1)/4+1
      KCOL_E(KLD_E(I+1)) = ELEMS(IEL)
      KLD_E(I+1) = KLD_E(I+1) + 1
    END DO
    DO IA=KLD(IDOF),KLD(IDOF+1)-1
     IB = IB + 1
     MAT_EX(IB) = BXMAT(IA)
     MAT_EY(IB) = BYMAT(IA)
     MAT_EZ(IB) = BZMAT(IA)
    END DO
   END DO

   NNEL = MGE013(ILEV)%ST(pID)%Num+1
   CALL SENDK_myMPI(KLD_E,NNEL,pID)
   NNEL = MGE013(ILEV)%SP(pID)%nEntries(1)
   CALL SENDK_myMPI(KCOL_E,NNEL,pID)
   CALL SENDD_myMPI(MAT_EX,4*NNEL,pID)
   CALL SENDD_myMPI(MAT_EY,4*NNEL,pID)
   CALL SENDD_myMPI(MAT_EZ,4*NNEL,pID)

   DEALLOCATE (KCOL_E,KLD_E)
   DEALLOCATE (MAT_EX,MAT_EY,MAT_EZ)

  END IF
 ELSE
  pNEL = 0
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     NNEL = MGE013(ILEV)%ST(pJD)%Num+1
     ALLOCATE (KLD_E(NNEL))
     KLD_E=0
     CALL RECVK_myMPI(KLD_E,NNEL,pJD)
     NNEL = MGE013(ILEV)%SP(pJD)%nEntries(2)
     ALLOCATE (KCOL_E(NNEL))
     KCOL_E=0
     ALLOCATE (MAT_EX(4*NNEL))
     ALLOCATE (MAT_EY(4*NNEL))
     ALLOCATE (MAT_EZ(4*NNEL))
     MAT_EX=0
     MAT_EY=0
     MAT_EZ=0     
     CALL RECVK_myMPI(KCOL_E,NNEL,pJD)
     CALL RECVD_myMPI(MAT_EX,4*NNEL,pJD)
     CALL RECVD_myMPI(MAT_EY,4*NNEL,pJD)
     CALL RECVD_myMPI(MAT_EZ,4*NNEL,pJD)
     DO I=1,MGE013(ILEV)%ST(pJD)%Num
      IDOF=MGE013(ILEV)%ST(pJD)%VertLink(2,I)
      DO IA=KLD_E(I),KLD_E(I+1)-1
       IB = KLD_T(IDOF)
       KCOL_B(IB+0) = 4*(pNel + KCOL_E(IA)-1)+1
       KCOL_B(IB+1) = 4*(pNel + KCOL_E(IA)-1)+2
       KCOL_B(IB+2) = 4*(pNel + KCOL_E(IA)-1)+3
       KCOL_B(IB+3) = 4*(pNel + KCOL_E(IA)-1)+4

       BXPMat(IB+0) = MAT_EX(4*(IA-1)+1)
       BXPMat(IB+1) = MAT_EX(4*(IA-1)+2)
       BXPMat(IB+2) = MAT_EX(4*(IA-1)+3)
       BXPMat(IB+3) = MAT_EX(4*(IA-1)+4)

       BYPMat(IB+0) = MAT_EY(4*(IA-1)+1)
       BYPMat(IB+1) = MAT_EY(4*(IA-1)+2)
       BYPMat(IB+2) = MAT_EY(4*(IA-1)+3)
       BYPMat(IB+3) = MAT_EY(4*(IA-1)+4)

       BZPMat(IB+0) = MAT_EZ(4*(IA-1)+1)
       BZPMat(IB+1) = MAT_EZ(4*(IA-1)+2)
       BZPMat(IB+2) = MAT_EZ(4*(IA-1)+3)
       BZPMat(IB+3) = MAT_EZ(4*(IA-1)+4)

       KLD_T(IDOF) = KLD_T(IDOF) + 4
      END DO
     END DO

     pNEL = pNel + MGE013(ILEV)%SP(pJD)%nElems(2)
     DEALLOCATE (KCOL_E,KLD_E)
     DEALLOCATE (MAT_EX,MAT_EY,MAT_EZ)

    END IF
   END IF
  END DO
 END IF
END DO

END SUBROUTINE Create_ParB_COLMAT
!
!
!
SUBROUTINE E013_CreateComm(DCORVG,DCORAG,KVERT,KEDGE,KAREA,&
NVT,NET,NAT,NEL,iMedLev)

USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX

IMPLICIT NONE

INTEGER iMedLev
CHARACTER*14 ccf
REAL*8  DCORAG(3,*),DCORVG(3,*)
INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)
INTEGER NAT,NEL,NVT,NET,ndof
REAL*8, ALLOCATABLE :: dCoor(:,:)
INTEGER,ALLOCATABLE :: iCoor(:),iAux(:)
INTEGER :: IAT,JAT,IEL,pID,pJD,I,J,K,nSize,jAux,pNDOF
INTEGER :: IVT(4),IV(4)
INTEGER :: IET(4),IE(4)
INTEGER Neigh(2,12)
DATA Neigh/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
REAL*8 P1X,P1Y,P1Z,P2X,P2Y,P2Z

IF (.NOT.ALLOCATED(MGE013)) ALLOCATE(MGE013(NLMIN:NLMAX))

ndof = NVT+NET+NAT+NEL

IF (myid.EQ.0) GOTO 88

ALLOCATE (iAux(1:ndof))
ALLOCATE (dCoor(3,1:ndof))
ALLOCATE (iCoor(1:ndof))


!IF (ALLOCATED(CoorSP)) DEALLOCATE (CoorSP)
ALLOCATE(CoorSP(subnodes))
 CoorSP(:)%Num = 0

DO pJD=1,subnodes

 nSize = MGE013(ILEV-1)%ST(pJD)%Num
 IF (nSIZE.EQ.0) THEN
  CYCLE
 END IF

 iAux = 0
 jAux = 0

 DO i=1,nSize
   j = MGE013(ILEV-1)%ST(pJD)%VertLink(1,i)
   jAux = jAux + 1
   iAux(j) = 1
 END DO

 k=1
 DO i=1,NEL
  DO j=1,12
   IF (k.eq.KEDGE(j,i)) THEN
    ivt(1) = KVERT(Neigh(1,j),i)
    ivt(2) = KVERT(Neigh(2,j),i)
    IF (iAux(ivt(1)).EQ.1.AND.iAux(ivt(2)).EQ.1) THEN
     iAux(NVT + k) = 1
     jAux = jAux + 1
    END IF
    k = k + 1
   END IF
  END DO
 END DO

 DO pID=1,mg_mpi(ILEV)%NeighNum
  IF (mg_mpi(ILEV)%parST(pID)%Neigh.EQ.pJD) GOTO 1
 END DO

 GOTO 2

1 nSIZE = mg_mpi(ILEV)%parST(pID)%Num
!  pJD   = mg_mpi(ILEV)%parST(pID)%Neigh

 DO I=1,nSIZE
  IEL = mg_mpi(ILEV)%parST(pID)%ElemLink(1,I)
  JAT = mg_mpi(ILEV)%parST(pID)%SideLink(  I)

  IF (JAT.eq.1) THEN
   IV(1) = 1; IV(2) = 2; IV(3) = 3; IV(4) = 4;
   IE(1) = 1; IE(2) = 2; IE(3) = 3; IE(4) = 4;
  END IF

  IF (JAT.eq.2) THEN
   IV(1) = 1; IV(2) = 2; IV(3) = 6; IV(4) = 5;
   IE(1) = 1; IE(2) = 6; IE(3) = 9; IE(4) = 5;
  END IF

  IF (JAT.eq.3) THEN
   IV(1) = 2; IV(2) = 3; IV(3) = 7; IV(4) = 6;
   IE(1) = 2; IE(2) = 7; IE(3) = 10; IE(4) = 6;
  END IF

  IF (JAT.eq.4) THEN
   IV(1) = 3; IV(2) = 4; IV(3) = 8; IV(4) = 7;
   IE(1) = 3; IE(2) = 8; IE(3) = 11; IE(4) = 7;
  END IF

  IF (JAT.eq.5) THEN
   IV(1) = 4; IV(2) = 1; IV(3) = 5; IV(4) = 8;
   IE(1) = 4; IE(2) = 5; IE(3) = 12; IE(4) = 8;
  END IF

  IF (JAT.eq.6) THEN
   IV(1) = 5; IV(2) = 6; IV(3) = 7; IV(4) = 8;
   IE(1) = 9; IE(2) = 10; IE(3) = 11; IE(4) = 12;
  END IF

  IVT(1) = KVERT(IV(1),IEL); IVT(2) = KVERT(IV(2),IEL)
  IVT(3) = KVERT(IV(3),IEL); IVT(4) = KVERT(IV(4),IEL)

  IET(1) = KEDGE(IE(1),IEL); IET(2) = KEDGE(IE(2),IEL)
  IET(3) = KEDGE(IE(3),IEL); IET(4) = KEDGE(IE(4),IEL)

  IAT    = KAREA(JAT,IEL)

  IF (iAUX(IVT(1)).EQ.0) jAux = jAux + 1
  IF (iAUX(IVT(2)).EQ.0) jAux = jAux + 1
  IF (iAUX(IVT(3)).EQ.0) jAux = jAux + 1
  IF (iAUX(IVT(4)).EQ.0) jAux = jAux + 1

  IF (iAUX(NVT+IET(1)).EQ.0) jAux = jAux + 1
  IF (iAUX(NVT+IET(2)).EQ.0) jAux = jAux + 1
  IF (iAUX(NVT+IET(3)).EQ.0) jAux = jAux + 1
  IF (iAUX(NVT+IET(4)).EQ.0) jAux = jAux + 1

  IF (iAUX(NVT+NET+IAT).EQ.0) jAux = jAux + 1

!   write(*,*) IVT(1),IVT(2),IVT(3),IVT(4)
  iAUX(IVT(1)) = 1; iAUX(IVT(2)) = 1
  iAUX(IVT(3)) = 1; iAUX(IVT(4)) = 1

  iAUX(NVT+IET(1)) = 1; iAUX(NVT+IET(2)) = 1
  iAUX(NVT+IET(3)) = 1; iAUX(NVT+IET(4)) = 1

  iAUX(NVT+NET+IAT) = 1;
 END DO

2 CONTINUE

 ALLOCATE(CoorSP(pJD)%dCoor(3,jAux))
 ALLOCATE(CoorSP(pJD)%iCoor(  jAux))
 CoorSP(pJD)%Num = jAux

 jAux = 0
 DO i=1,NVT
  IF (iAux(i).EQ.1) THEN
   jAux = jAux + 1
   CoorSP(pJD)%dCoor(1,jAux) = DCORVG(1,i)
   CoorSP(pJD)%dCoor(2,jAux) = DCORVG(2,i)
   CoorSP(pJD)%dCoor(3,jAux) = DCORVG(3,i)
   CoorSP(pJD)%iCoor(  jAux) = i
  END IF
 END DO

 k=1
 DO i=1,NEL
  DO j=1,12
   IF (k.eq.KEDGE(j,i)) THEN
    IF (iAux(NVT+k).EQ.1) THEN
     jAux = jAux + 1
     ivt(1) = KVERT(Neigh(1,j),i)
     ivt(2) = KVERT(Neigh(2,j),i)
     CoorSP(pJD)%dCoor(1,jAux) = 0.5d0*(DCORVG(1,ivt(1))+DCORVG(1,ivt(2)))
     CoorSP(pJD)%dCoor(2,jAux) = 0.5d0*(DCORVG(2,ivt(1))+DCORVG(2,ivt(2)))
     CoorSP(pJD)%dCoor(3,jAux) = 0.5d0*(DCORVG(3,ivt(1))+DCORVG(3,ivt(2)))
     CoorSP(pJD)%iCoor(  jAux) = NVT+k
    END IF
    k = k + 1
   END IF
  END DO
 END DO

 DO i=1,NAT
  IF (iAux(NVT+NET+i).EQ.1) THEN
   jAux = jAux + 1
   CoorSP(pJD)%dCoor(1,jAux) = DCORAG(1,i)
   CoorSP(pJD)%dCoor(2,jAux) = DCORAG(2,i)
   CoorSP(pJD)%dCoor(3,jAux) = DCORAG(3,i)
   CoorSP(pJD)%iCoor(  jAux) = NVT+NET+i
  END IF
 END DO

END DO

! IF (ALLOCATED(CoorST)) DEALLOCATE (CoorST)
ALLOCATE(CoorST(subnodes))

! write(*,'(I3,2(A,<subnodes>I5))') myid," : ", MGE013(ILEV-1)%ST(:)%Num
! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
! write(*,'(I3,2(A,<subnodes>I5))') myid," : ",CoorSP(1:subnodes)%Num

!write(*,'(I3,2(A,<subnodes>I5))') myid," : ",CoorST(1:subnodes)%Num
! write(*,'(I3,2(A,<subnodes>I5))') myid," : ",CoorSP(1:subnodes)%Num," : ",CoorST(1:subnodes)%Num
! pause

DO pID=1,subnodes
 IF (myid.eq.pID) THEN
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    CALL RECVI_myMPI(pNDOF,pJD)
    CoorST(pJD)%Num = pNDOF
    IF (pNDOF.GT.0) THEN
     ALLOCATE(CoorST(pJD)%dCoor(3,pNDOF))
     ALLOCATE(CoorST(pJD)%iCoor(  pNDOF))
     CALL RECVD_myMPI(CoorST(pJD)%dCoor,3*pNDOF,pJD)
     CALL RECVK_myMPI(CoorST(pJD)%iCoor,  pNDOF,pJD)
    END IF
!    ELSE
!     CoorST(pJD)%Num = jAux
!     ALLOCATE(CoorST(pJD)%dCoor(3,jAux))
!     ALLOCATE(CoorST(pJD)%iCoor(  jAux))
!     CoorST(pJD)%dCoor(:,:) = dCoor(:,1:jAux)
!     CoorST(pJD)%iCoor(  :) = iCoor(  1:jAux)
   END IF
  END DO
 ELSE
  jAux = CoorSP(pID)%Num
  CALL SENDI_myMPI(jAux ,pID)
  IF (jAux.GT.0) THEN
   CALL SENDD_myMPI(CoorSP(pID)%dCoor,3*jAux,pID)
   CALL SENDK_myMPI(CoorSP(pID)%iCoor,  jAux,pID)
  END IF
 END IF

END DO

! WRITE(ccf,'(A,I1.1,A,I3.3,A)') "comm_",ILEV,"_",myid,".txt"
! OPEN (FILE=ccf,UNIT=994)

! write(*,'(I3,2(A,<subnodes>I5))') myid," : ",CoorSP(1:subnodes)%Num
! pause
! if (myid.eq.1) OPEN(987,FILE='VERT1.txt')
! if (myid.eq.2) OPEN(987,FILE='VERT2.txt')
! if (myid.eq.3) OPEN(987,FILE='VERT3.txt')
! 
! DO pID=1,subnodes
! WRITE(987,*) pID
! DO I=1,CoorST(pID)%Num
!  WRITE(987,*) CoorST(pID)%dCoor(:,I)
! END DO
! END DO
! CLOSE(987)

ALLOCATE(MGE013(ILEV)%ST(subnodes))
ALLOCATE(MGE013(ILEV)%UE(NDOF))
ALLOCATE(MGE013(ILEV)%UE11(NDOF))
ALLOCATE(MGE013(ILEV)%UE22(NDOF))
ALLOCATE(MGE013(ILEV)%UE33(NDOF))

DO pID=1,subnodes
 MGE013(ILEV)%ST(pID)%Num = CoorSP(pID)%Num
END DO


! if (myid.eq.1) OPEN(987,ACCESS="APPEND",FILE='VERT1.txt')


DO pID=1,subnodes
 jAux = 0
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  ALLOCATE(MGE013(ILEV)%ST(pID)%VertLink(2,MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SVVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RVVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SDVect  (3*MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RDVect  (3*MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SBVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RBVect  (  MGE013(ILEV)%ST(pID)%Num))
  DO I=1,CoorST(pID)%Num
   P1X = CoorST(pID)%dCoor(1,I)
   P1Y = CoorST(pID)%dCoor(2,I)
   P1Z = CoorST(pID)%dCoor(3,I)
   DO J=1,CoorSP(pID)%Num
    P2X = CoorSP(pID)%dCoor(1,J)
    P2Y = CoorSP(pID)%dCoor(2,J)
    P2Z = CoorSP(pID)%dCoor(3,J)
    IF ((ABS(P1X-P2X).LT.DEpsPrec).AND.&
        (ABS(P1Y-P2Y).LT.DEpsPrec).AND.&
        (ABS(P1Z-P2Z).LT.DEpsPrec)) THEN
     jAux = jAux + 1
     MGE013(ILEV)%ST(pID)%VertLink(1,jAux) = CoorSP(pID)%iCoor(J)
     MGE013(ILEV)%ST(pID)%VertLink(2,jAux) = CoorSP(pID)%iCoor(J)
     EXIT
    END IF
   END DO
  END DO
!  IF (jAux.EQ.MGE013(ILEV)%ST(pID)%Num) write(*,*) "problem ",jAux,MGE013(ILEV)%ST(pID)%Num
  CALL SORT1D(MGE013(ILEV)%ST(pID)%VertLink(1,:),MGE013(ILEV)%ST(pID)%Num)

!   WRITE(994,*) pID,MGE013(ILEV)%ST(pID)%Num,"---"
  DO I=1,MGE013(ILEV)%ST(pID)%Num
!    WRITE(994,*) MGE013(ILEV)%ST(pID)%VertLink(:,I)
  END DO

 END IF
END DO

! CLOSE(994)

! if (myid.eq.1) OPEN(987,ACCESS="APPEND",FILE='VERT1.txt')
! if (myid.eq.2) OPEN(987,ACCESS="APPEND",FILE='VERT2.txt')
! if (myid.eq.3) OPEN(987,ACCESS="APPEND",FILE='VERT3.txt')
! 
! WRITE(987,*) "--ILEV--",ILEV
! DO pID=1,subnodes
! WRITE(987,*) "--pID--",pID
! DO I=1,MGE013(ILEV)%ST(pID)%Num
!  WRITE(987,*) I,MGE013(ILEV)%ST(pID)%VertLink(:,I)
! END DO
! END DO
! CLOSE(987)

DEALLOCATE (iAux,dCoor,iCoor)
DEALLOCATE (CoorSP,CoorST)

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

88 CONTINUE

IF (ILEV.eq.iMedLev) THEN
 IF (myid.eq.MASTER) THEN
  jAUX = 0
  DO pID = 1,subnodes
   if (coarse%pNEL(pID)*(8**(iMedLev-1)).GT.jAUX) jAUX = coarse%pNEL(pID)*(8**(iMedLev-1))
  END DO
  ALLOCATE (MGE013(ILEV)%CRSVect(jAUX*4))
 ELSE
  ALLOCATE (MGE013(ILEV)%CRSVect(NEL*4))
 END IF
END IF

END


! ----------------------------------------------
SUBROUTINE E011_CreateComm(NVT)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE

INTEGER NVT,pID,jAux,i

IF (myid.EQ.0) RETURN

ALLOCATE(E011ST(subnodes))

ILEV = NLMAX

DO pID=1,subnodes
  E011ST(pID)%Num = MGE013(ILEV)%ST(pID)%Num
END DO

ALLOCATE(E011_UE(NVT))

DO pID=1,subnodes
 IF (E011ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  jAux = E011ST(pID)%Num
  ALLOCATE(E011ST(pID)%VertLink(2,E011ST(pID)%Num))
  ALLOCATE(E011ST(pID)%SVVect  (  E011ST(pID)%Num))
  ALLOCATE(E011ST(pID)%RVVect  (  E011ST(pID)%Num))
  ALLOCATE(E011ST(pID)%SDVect  (  E011ST(pID)%Num))
  ALLOCATE(E011ST(pID)%RDVect  (  E011ST(pID)%Num))
  ALLOCATE(E011ST(pID)%SBVect  (  E011ST(pID)%Num))
  ALLOCATE(E011ST(pID)%RBVect  (  E011ST(pID)%Num))
  DO i=1,jAux
   E011ST(pID)%VertLink(1,i) =  MGE013(ILEV)%ST(pID)%VertLink(1,i)
   E011ST(pID)%VertLink(2,i) =  MGE013(ILEV)%ST(pID)%VertLink(2,i)
  END DO
 END IF
END DO

END



SUBROUTINE E013_CreateComm_read(DCORVG,DCORAG,KVERT,KEDGE,KAREA,&
NVT,NET,NAT,NEL,iMedLev)

USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX

IMPLICIT NONE

INTEGER iMedLev
CHARACTER*14 ccf
REAL*8  DCORAG(3,*),DCORVG(3,*)
INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)
INTEGER NAT,NEL,NVT,NET,ndof
REAL*8, ALLOCATABLE :: dCoor(:,:)
INTEGER,ALLOCATABLE :: iCoor(:),iAux(:)
INTEGER :: IAT,JAT,IEL,pID,pJD,I,J,K,nSize,jAux,pNDOF
INTEGER :: IVT(4),IV(4)
INTEGER :: IET(4),IE(4)
INTEGER Neigh(2,12)
DATA Neigh/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
REAL*8 P1X,P1Y,P1Z,P2X,P2Y,P2Z

IF (.NOT.ALLOCATED(MGE013)) ALLOCATE(MGE013(NLMIN:NLMAX))

ndof = NVT+NET+NAT+NEL

IF (myid.EQ.0) GOTO 88

ALLOCATE(MGE013(ILEV)%ST(subnodes))
ALLOCATE(MGE013(ILEV)%UE(NDOF))
ALLOCATE(MGE013(ILEV)%UE11(NDOF))
ALLOCATE(MGE013(ILEV)%UE22(NDOF))
ALLOCATE(MGE013(ILEV)%UE33(NDOF))

! WRITE(ccf,'(A,I1.1,A,I3.3,A)') "comn_",ILEV,"_",myid,".txt"
! OPEN (FILE=ccf,UNIT=994)

MGE013(ILEV)%ST(:)%Num = 0

DO pJD=1,mg_mpi(ILEV)%NeighNum

  jAux = 0

!   READ(994,*) pID,MGE013(ILEV)%ST(pID)%Num

  ALLOCATE(MGE013(ILEV)%ST(pID)%VertLink(2,MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SVVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RVVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SDVect  (3*MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RDVect  (3*MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SBVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RBVect  (  MGE013(ILEV)%ST(pID)%Num))

  DO I=1,MGE013(ILEV)%ST(pID)%Num
!    READ(994,*) MGE013(ILEV)%ST(pID)%VertLink(:,I)
  END DO

END DO

! CLOSE(994)

88 CONTINUE

IF (ILEV.eq.iMedLev) THEN
 IF (myid.eq.MASTER) THEN
  jAUX = 0
  DO pID = 1,subnodes
   if (coarse%pNEL(pID)*(8**(iMedLev-1)).GT.jAUX) jAUX = coarse%pNEL(pID)*(8**(iMedLev-1))
  END DO
  ALLOCATE (MGE013(ILEV)%CRSVect(jAUX*4))
 ELSE
  ALLOCATE (MGE013(ILEV)%CRSVect(NEL*4))
 END IF
END IF

END


SUBROUTINE E013_CreateComm_coarse(DCORVG,DCORAG,KVERT,KEDGE,KAREA,&
NVT,NET,NAT,NEL,iMedLev)

USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX

IMPLICIT NONE

INTEGER iMedLev
CHARACTER*14 ccf
REAL*8  DCORAG(3,*),DCORVG(3,*)
INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)
INTEGER NAT,NEL,NVT,NET,ndof
REAL*8, ALLOCATABLE :: dCoor(:,:)
INTEGER,ALLOCATABLE :: iCoor(:),iAux(:)
INTEGER :: IAT,JAT,IEL,pID,pJD,I,J,K,nSize,jAux,pNDOF
INTEGER :: IVT(4),IV(4)
INTEGER :: IET(4),IE(4)
INTEGER Neigh(2,12)
DATA Neigh/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
REAL*8 P1X,P1Y,P1Z,P2X,P2Y,P2Z

IF (.NOT.ALLOCATED(MGE013)) ALLOCATE(MGE013(NLMIN:NLMAX))

ndof = NVT+NET+NAT+NEL

IF (myid.EQ.0) GOTO 88

ALLOCATE (iAux(1:ndof))
ALLOCATE (dCoor(3,1:ndof))
ALLOCATE (iCoor(1:ndof))

DO pID=1,mg_mpi(ILEV)%NeighNum
 nSIZE = mg_mpi(ILEV)%parST(pID)%Num
 DO I=1,nSIZE
  IEL = mg_mpi(ILEV)%parST(pID)%ElemLink(1,I)
  JAT = mg_mpi(ILEV)%parST(pID)%SideLink(  I)

  IF (JAT.eq.1) THEN
   IV(1) = 1; IV(2) = 2; IV(3) = 3; IV(4) = 4;
   IE(1) = 1; IE(2) = 2; IE(3) = 3; IE(4) = 4;
  END IF

  IF (JAT.eq.2) THEN
   IV(1) = 1; IV(2) = 2; IV(3) = 6; IV(4) = 5;
   IE(1) = 1; IE(2) = 6; IE(3) = 9; IE(4) = 5;
  END IF

  IF (JAT.eq.3) THEN
   IV(1) = 2; IV(2) = 3; IV(3) = 7; IV(4) = 6;
   IE(1) = 2; IE(2) = 7; IE(3) = 10; IE(4) = 6;
  END IF

  IF (JAT.eq.4) THEN
   IV(1) = 3; IV(2) = 4; IV(3) = 8; IV(4) = 7;
   IE(1) = 3; IE(2) = 8; IE(3) = 11; IE(4) = 7;
  END IF

  IF (JAT.eq.5) THEN
   IV(1) = 4; IV(2) = 1; IV(3) = 5; IV(4) = 8;
   IE(1) = 4; IE(2) = 5; IE(3) = 12; IE(4) = 8;
  END IF

  IF (JAT.eq.6) THEN
   IV(1) = 5; IV(2) = 6; IV(3) = 7; IV(4) = 8;
   IE(1) = 9; IE(2) = 10; IE(3) = 11; IE(4) = 12;
  END IF

  IVT(1) = KVERT(IV(1),IEL); IVT(2) = KVERT(IV(2),IEL)
  IVT(3) = KVERT(IV(3),IEL); IVT(4) = KVERT(IV(4),IEL)

  IET(1) = KEDGE(IE(1),IEL); IET(2) = KEDGE(IE(2),IEL)
  IET(3) = KEDGE(IE(3),IEL); IET(4) = KEDGE(IE(4),IEL)

  IAT    = KAREA(JAT,IEL)

!   write(*,*) IVT(1),IVT(2),IVT(3),IVT(4)
  iAUX(IVT(1)) = 1; iAUX(IVT(2)) = 1
  iAUX(IVT(3)) = 1; iAUX(IVT(4)) = 1

  iAUX(NVT+IET(1)) = 1; iAUX(NVT+IET(2)) = 1
  iAUX(NVT+IET(3)) = 1; iAUX(NVT+IET(4)) = 1

  iAUX(NVT+NET+IAT) = 1;
 END DO
END DO

jAux = 0
DO i=1,NVT
 IF (iAux(i).EQ.1) THEN
  jAux = jAux + 1
  dCoor(1,jAux) = DCORVG(1,i)
  dCoor(2,jAux) = DCORVG(2,i)
  dCoor(3,jAux) = DCORVG(3,i)
  iCoor(  jAux) = i
 END IF
END DO

k=1
DO i=1,NEL
 DO j=1,12
  IF (k.eq.KEDGE(j,i)) THEN
   IF (iAux(NVT+k).EQ.1) THEN
    jAux = jAux + 1
    ivt(1) = KVERT(Neigh(1,j),i)
    ivt(2) = KVERT(Neigh(2,j),i)
    dCoor(1,jAux) = 0.5d0*(DCORVG(1,ivt(1))+DCORVG(1,ivt(2)))
    dCoor(2,jAux) = 0.5d0*(DCORVG(2,ivt(1))+DCORVG(2,ivt(2)))
    dCoor(3,jAux) = 0.5d0*(DCORVG(3,ivt(1))+DCORVG(3,ivt(2)))
    iCoor(  jAux) = NVT+k
   END IF
   k = k + 1
  END IF
 END DO
END DO

DO i=1,NAT
 IF (iAux(NVT+NET+i).EQ.1) THEN
  jAux = jAux + 1
  dCoor(1,jAux) = DCORAG(1,i)
  dCoor(2,jAux) = DCORAG(2,i)
  dCoor(3,jAux) = DCORAG(3,i)
  iCoor(  jAux) = NVT+NET+i
 END IF
END DO


IF (ALLOCATED(CoorST)) DEALLOCATE (CoorST)
ALLOCATE(CoorST(subnodes))

DO pID=1,subnodes
 IF (myid.eq.pID) THEN
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    CALL RECVI_myMPI(pNDOF,pJD)
    CoorST(pJD)%Num = pNDOF
    ALLOCATE(CoorST(pJD)%dCoor(3,pNDOF))
    ALLOCATE(CoorST(pJD)%iCoor(  pNDOF))
    CALL RECVD_myMPI(CoorST(pJD)%dCoor,3*pNDOF,pJD)
    CALL RECVK_myMPI(CoorST(pJD)%iCoor,  pNDOF,pJD)
   ELSE
    CoorST(pJD)%Num = jAux
    ALLOCATE(CoorST(pJD)%dCoor(3,jAux))
    ALLOCATE(CoorST(pJD)%iCoor(  jAux))
    CoorST(pJD)%dCoor(:,:) = dCoor(:,1:jAux)
    CoorST(pJD)%iCoor(  :) = iCoor(  1:jAux)
   END IF
  END DO
 ELSE
  CALL SENDI_myMPI(jAux ,pID)
  CALL SENDD_myMPI(dCoor,3*jAux,pID)
  CALL SENDK_myMPI(iCoor,  jAux,pID)
 END IF

END DO

! if (myid.eq.1) OPEN(987,FILE='VERT1.txt')
! if (myid.eq.2) OPEN(987,FILE='VERT2.txt')
! if (myid.eq.3) OPEN(987,FILE='VERT3.txt')
! 
! DO pID=1,subnodes
! WRITE(987,*) pID
! DO I=1,CoorST(pID)%Num
!  WRITE(987,*) CoorST(pID)%dCoor(:,I)
! END DO
! END DO
! CLOSE(987)

ALLOCATE(MGE013(ILEV)%ST(subnodes))
ALLOCATE(MGE013(ILEV)%UE(NDOF))
ALLOCATE(MGE013(ILEV)%UE11(NDOF))
ALLOCATE(MGE013(ILEV)%UE22(NDOF))
ALLOCATE(MGE013(ILEV)%UE33(NDOF))

DO pID=1,subnodes
 jAux = 0
 IF (pID.NE.myid) THEN
  DO I=1,CoorST(pID)%Num
   P1X = CoorST(pID)%dCoor(1,I)
   P1Y = CoorST(pID)%dCoor(2,I)
   P1Z = CoorST(pID)%dCoor(3,I)
   DO J=1,CoorST(myid)%Num
    P2X = CoorST(myid)%dCoor(1,J)
    P2Y = CoorST(myid)%dCoor(2,J)
    P2Z = CoorST(myid)%dCoor(3,J)
    IF ((ABS(P1X-P2X).LT.DEpsPrec).AND.&
        (ABS(P1Y-P2Y).LT.DEpsPrec).AND.&
        (ABS(P1Z-P2Z).LT.DEpsPrec)) THEN
     jAux = jAux + 1
    END IF
   END DO
  END DO
 END IF
!  WRITE(*,*) myid,pID,jAux
 MGE013(ILEV)%ST(pID)%Num = jAux
END DO

! WRITE(ccf,'(A,I1.1,A,I3.3,A)') "comn_",ILEV,"_",myid,".txt"
! OPEN (FILE=ccf,UNIT=994)

DO pID=1,subnodes
 jAux = 0
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  ALLOCATE(MGE013(ILEV)%ST(pID)%VertLink(2,MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SVVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RVVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SDVect  (3*MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RDVect  (3*MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%SBVect  (  MGE013(ILEV)%ST(pID)%Num))
  ALLOCATE(MGE013(ILEV)%ST(pID)%RBVect  (  MGE013(ILEV)%ST(pID)%Num))
  DO I=1,CoorST(pID)%Num
   P1X = CoorST(pID)%dCoor(1,I)
   P1Y = CoorST(pID)%dCoor(2,I)
   P1Z = CoorST(pID)%dCoor(3,I)
   DO J=1,CoorST(myid)%Num
    P2X = CoorST(myid)%dCoor(1,J)
    P2Y = CoorST(myid)%dCoor(2,J)
    P2Z = CoorST(myid)%dCoor(3,J)
    IF ((ABS(P1X-P2X).LT.DEpsPrec).AND.&
        (ABS(P1Y-P2Y).LT.DEpsPrec).AND.&
        (ABS(P1Z-P2Z).LT.DEpsPrec)) THEN
     jAux = jAux + 1
     MGE013(ILEV)%ST(pID)%VertLink(1,jAux) = CoorST(myid)%iCoor(J)
     MGE013(ILEV)%ST(pID)%VertLink(2,jAux) = CoorST(myid)%iCoor(J)
     EXIT
    END IF
   END DO
  END DO
  CALL SORT1D(MGE013(ILEV)%ST(pID)%VertLink(1,:),MGE013(ILEV)%ST(pID)%Num)

!   WRITE(994,*) pID,MGE013(ILEV)%ST(pID)%Num,"---"
  DO I=1,MGE013(ILEV)%ST(pID)%Num
!    WRITE(994,*) MGE013(ILEV)%ST(pID)%VertLink(:,I)
  END DO

 END IF
END DO

! CLOSE(994)

! if (myid.eq.1) OPEN(987,ACCESS="APPEND",FILE='VERT1.txt')
! if (myid.eq.2) OPEN(987,ACCESS="APPEND",FILE='VERT2.txt')
! if (myid.eq.3) OPEN(987,ACCESS="APPEND",FILE='VERT3.txt')
! 
! WRITE(987,*) "--ILEV--",ILEV
! DO pID=1,subnodes
! WRITE(987,*) "--pID--",pID
! DO I=1,MGE013(ILEV)%ST(pID)%Num
!  WRITE(987,*) I,MGE013(ILEV)%ST(pID)%VertLink(:,I)
! END DO
! END DO
! CLOSE(987)

DEALLOCATE (iAux,dCoor,iCoor)
DEALLOCATE (CoorST)


if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

88 CONTINUE

IF (ILEV.eq.iMedLev) THEN
 IF (myid.eq.MASTER) THEN
  jAUX = 0
  DO pID = 1,subnodes
   if (coarse%pNEL(pID)*(8**(iMedLev-1)).GT.jAUX) jAUX = coarse%pNEL(pID)*(8**(iMedLev-1))
  END DO
  ALLOCATE (MGE013(ILEV)%CRSVect(jAUX*4))
 ELSE
  ALLOCATE (MGE013(ILEV)%CRSVect(NEL*4))
 END IF
END IF

END


! In the following I will try to extend the pressure gradient
! matrix with the parallel blocks. Essentially I will try to 
! incorporate those coefficients in the given line of velocity 
! DOF  into the matrix, which are from the parallel block

! First I need to know which are the elements from the parallel 
! blocks, that need to share their coefficients for the parallel 
! velocity DOFS

! I also need to extend the pressure field by the ones from the 
! parallel blocks

SUBROUTINE Extract_ParElems_Orig(BXMAT,BYMAT,BZMAT,KCOL,KLD,NEQ,NEL,pNEL)
USE PP3D_MPI
USE def_QuadScalar, ONLY : mg_qlPMat,mg_BXPMat,mg_BYPMat,mg_BZPMat
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
!-------------------------------------------------------
!-------------------------------------------------------
REAL*8 BXMAT(*),BYMAT(*),BZMAT(*)
INTEGER KCOL(*),KLD(*),NEQ,NEL,pNEL
INTEGER pID,pJD,I,J,IDOF,IA,IB,IEL,NNEL,MNEL
INTEGER ELEMS(NEL)
INTEGER, ALLOCATABLE :: KCOL_E(:),KLD_E(:),EntNum(:),BP_LDC(:)
REAL*8, ALLOCATABLE ::  MAT_EX(:),MAT_EY(:),MAT_EZ(:)
CHARACTER*10 myFile

ALLOCATE(MGE013(ILEV)%SP(subnodes))

! Get the number of elements that are needed from the other side 
DO pID=1,subnodes
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  ELEMS = 0
  DO I=1,MGE013(ILEV)%ST(pID)%Num
   IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
!    write(*,*) (KCOL(IA),IA=KLD(IDOF),KLD(IDOF+1)-1)
   DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
     IEL = (KCOL(IA)-1)/4+1
     ELEMS(IEL) = 1
   END DO
  END DO
  NNEL = 0
  DO IEL=1,NEL
   IF (ELEMS(IEL).NE.0) NNEL=NNEL+1
  END DO
  MGE013(ILEV)%SP(pID)%nElems(1) = NNEL
 END IF
END DO

! Send them!
DO pID=1,subnodes
 IF (myid.eq.pID) THEN
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    CALL RECVI_myMPI(NNEL,pJD)
    MGE013(ILEV)%SP(pJD)%nElems(2) = NNEL
   END IF
  END DO
 ELSE
  NNEL = MGE013(ILEV)%SP(pID)%nElems(1)
  CALL SENDI_myMPI(NNEL ,pID)
 END IF
END DO

! Get the send-list of parallel elements
DO pID=1,subnodes
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  ELEMS = 0
  DO I=1,MGE013(ILEV)%ST(pID)%Num
   IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
   DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
     IEL = (KCOL(IA)-1)/4+1
     ELEMS(IEL) = 1
   END DO
  END DO
  NNEL = 0
  DO IEL=1,NEL
   IF (ELEMS(IEL).NE.0) NNEL=NNEL+1
  END DO
  MGE013(ILEV)%SP(pID)%Num = MAX(MGE013(ILEV)%SP(pID)%nElems(1),MGE013(ILEV)%SP(pID)%nElems(2))
  ALLOCATE(MGE013(ILEV)%SP(pID)%VertLink(2,MGE013(ILEV)%SP(pID)%Num))
  ALLOCATE(MGE013(ILEV)%SP(pID)%SDVect(4*MGE013(ILEV)%SP(pID)%Num))
  ALLOCATE(MGE013(ILEV)%SP(pID)%RDVect(4*MGE013(ILEV)%SP(pID)%Num))
!   ALLOCATE(MGE013(ILEV)%SP(pID)%SVVect(MGE013(ILEV)%SP(pID)%Num))
!   ALLOCATE(MGE013(ILEV)%SP(pID)%RVVect(MGE013(ILEV)%SP(pID)%Num))
  NNEL = 0
  DO IEL=1,NEL
   IF (ELEMS(IEL).NE.0) THEN
    NNEL = NNEL + 1
    MGE013(ILEV)%SP(pID)%VertLink(1,NNEL) = IEL
   END IF
  END DO
 END IF
END DO

! Get the receive-list of parallel elements
MNEL = 0
DO pID=1,subnodes
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
  NNEL = 0
  DO I=1,MGE013(ILEV)%SP(pID)%nElems(2)
   NNEL = NNEL + 1
   MNEL = MNEL + 1
   MGE013(ILEV)%SP(pID)%VertLink(2,NNEL) = MNEL
  END DO
 END IF
END DO

! Prepare and allocate structrures of local KLD_P
mg_qlPMat(ILEV)%nu = NEQ
IF (ALLOCATED(BP_LDC)) DEALLOCATE (BP_LDC)
ALLOCATE(mg_qlPMat(ILEV)%LdA(mg_qlPMat(ILEV)%nu+1))
ALLOCATE(BP_LDC(mg_qlPMat(ILEV)%nu+1))
mg_qlPMat(ILEV)%LdA = 0

! Get the length of parallel KCOL and send it over!
DO pID=1,subnodes
 IF (pID.NE.myid) THEN
  IF (MGE013(ILEV)%ST(pID)%Num.GT.0) THEN
   NNEL  = 0
   ALLOCATE (EntNum(MGE013(ILEV)%ST(pID)%Num))
   EntNum = 0
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
    DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
      EntNum(I) = EntNum(I) + 1
      NNEL = NNEL + 1
    END DO
   END DO

   MGE013(ILEV)%SP(pID)%nEntries(1) = NNEL
   CALL SENDI_myMPI(NNEL ,pID)
   NNEL = MGE013(ILEV)%ST(pID)%Num
   CALL SENDK_myMPI(EntNum,NNEL,pID)
   DEALLOCATE (EntNum)

  END IF
 ELSE
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     CALL RECVI_myMPI(NNEL,pJD)
     MGE013(ILEV)%SP(pJD)%nEntries(2) = NNEL
     NNEL = MGE013(ILEV)%ST(pJD)%Num
     ALLOCATE (EntNum(NNEL))
     CALL RECVK_myMPI(EntNum,NNEL,pJD)
     DO I=1,MGE013(ILEV)%ST(pJD)%Num
      IDOF=MGE013(ILEV)%ST(pJD)%VertLink(2,I)
      mg_qlPMat(ILEV)%LdA(IDOF) = mg_qlPMat(ILEV)%LdA(IDOF) + EntNum(I)
     END DO
     DEALLOCATE (EntNum)

    END IF
   END IF
  END DO
 END IF
END DO

! Prepare and allocate structrures of local KCOL_P
!write(*,*) myid,"--",BP_LD(1)
BP_LDC=mg_qlPMat(ILEV)%LdA
mg_qlPMat(ILEV)%LdA(1)=1
DO IDOF=2,NEQ+1
 mg_qlPMat(ILEV)%LdA(IDOF) = 4*BP_LDC(IDOF-1) + mg_qlPMat(ILEV)%LdA(IDOF-1)
END DO
BP_LDC=mg_qlPMat(ILEV)%LdA
mg_qlPMat(ILEV)%na = 0
DO pID=1,subnodes
 mg_qlPMat(ILEV)%na = mg_qlPMat(ILEV)%na + MGE013(ILEV)%SP(pID)%nEntries(2)
END DO
ALLOCATE(mg_qlPMat(ILEV)%ColA(4*mg_qlPMat(ILEV)%na))
mg_qlPMat(ILEV)%ColA=0
ALLOCATE(mg_BXPMat(ILEV)%a(4*mg_qlPMat(ILEV)%na),mg_BYPMat(ILEV)%a(4*mg_qlPMat(ILEV)%na),mg_BZPMat(ILEV)%a(4*mg_qlPMat(ILEV)%na))


! Fill up the entries of parallel KCOL
DO pID=1,subnodes
 IF (pID.NE.myid) THEN
  IF (MGE013(ILEV)%ST(pID)%Num.GT.0) THEN
   ELEMS = 0
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
    DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
      IEL = (KCOL(IA)-1)/4+1
      ELEMS(IEL) = 1
    END DO
   END DO

   ALLOCATE (KCOL_E(MGE013(ILEV)%SP(pID)%nEntries(1)))
   ALLOCATE (KLD_E(MGE013(ILEV)%ST(pID)%Num+1))
   ALLOCATE (MAT_EX(4*MGE013(ILEV)%SP(pID)%nEntries(1)))
   ALLOCATE (MAT_EY(4*MGE013(ILEV)%SP(pID)%nEntries(1)))
   ALLOCATE (MAT_EZ(4*MGE013(ILEV)%SP(pID)%nEntries(1)))

   DO IEL=2,NEL
    ELEMS(IEL) = ELEMS(IEL) + ELEMS(IEL-1)
   END DO

   KLD_E(1)=1
   IB = 0
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    IDOF=MGE013(ILEV)%ST(pID)%VertLink(1,I)
    KLD_E(I+1) = KLD_E(I)
    DO IA=KLD(IDOF),KLD(IDOF+1)-1,4
      IEL = (KCOL(IA)-1)/4+1
      KCOL_E(KLD_E(I+1)) = ELEMS(IEL)
      KLD_E(I+1) = KLD_E(I+1) + 1
    END DO
    DO IA=KLD(IDOF),KLD(IDOF+1)-1
     IB = IB + 1
     MAT_EX(IB) = BXMAT(IA)
     MAT_EY(IB) = BYMAT(IA)
     MAT_EZ(IB) = BZMAT(IA)
    END DO
   END DO

   NNEL = MGE013(ILEV)%ST(pID)%Num+1
   CALL SENDK_myMPI(KLD_E,NNEL,pID)
   NNEL = MGE013(ILEV)%SP(pID)%nEntries(1)
   CALL SENDK_myMPI(KCOL_E,NNEL,pID)
   CALL SENDD_myMPI(MAT_EX,4*NNEL,pID)
   CALL SENDD_myMPI(MAT_EY,4*NNEL,pID)
   CALL SENDD_myMPI(MAT_EZ,4*NNEL,pID)

   DEALLOCATE (KCOL_E,KLD_E,MAT_EX,MAT_EY,MAT_EZ)

  END IF
 ELSE
  pNEL = 0
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     NNEL = MGE013(ILEV)%ST(pJD)%Num+1
     ALLOCATE (KLD_E(NNEL))
     CALL RECVK_myMPI(KLD_E,NNEL,pJD)
     NNEL = MGE013(ILEV)%SP(pJD)%nEntries(2)
     ALLOCATE (KCOL_E(NNEL))
     ALLOCATE (MAT_EX(4*NNEL))
     ALLOCATE (MAT_EY(4*NNEL))
     ALLOCATE (MAT_EZ(4*NNEL))
     CALL RECVK_myMPI(KCOL_E,NNEL,pJD)
     CALL RECVD_myMPI(MAT_EX,4*NNEL,pJD)
     CALL RECVD_myMPI(MAT_EY,4*NNEL,pJD)
     CALL RECVD_myMPI(MAT_EZ,4*NNEL,pJD)
     DO I=1,MGE013(ILEV)%ST(pJD)%Num
      IDOF=MGE013(ILEV)%ST(pJD)%VertLink(2,I)
      DO IA=KLD_E(I),KLD_E(I+1)-1
       IB = BP_LDC(IDOF)
       mg_qlPMat(ILEV)%ColA(IB+0) = 4*(pNel + KCOL_E(IA)-1)+1
       mg_qlPMat(ILEV)%ColA(IB+1) = 4*(pNel + KCOL_E(IA)-1)+2
       mg_qlPMat(ILEV)%ColA(IB+2) = 4*(pNel + KCOL_E(IA)-1)+3
       mg_qlPMat(ILEV)%ColA(IB+3) = 4*(pNel + KCOL_E(IA)-1)+4

       mg_BXPMat(ILEV)%a(IB+0) = MAT_EX(4*(IA-1)+1)
       mg_BXPMat(ILEV)%a(IB+1) = MAT_EX(4*(IA-1)+2)
       mg_BXPMat(ILEV)%a(IB+2) = MAT_EX(4*(IA-1)+3)
       mg_BXPMat(ILEV)%a(IB+3) = MAT_EX(4*(IA-1)+4)

       mg_BYPMat(ILEV)%a(IB+0) = MAT_EY(4*(IA-1)+1)
       mg_BYPMat(ILEV)%a(IB+1) = MAT_EY(4*(IA-1)+2)
       mg_BYPMat(ILEV)%a(IB+2) = MAT_EY(4*(IA-1)+3)
       mg_BYPMat(ILEV)%a(IB+3) = MAT_EY(4*(IA-1)+4)

       mg_BZPMat(ILEV)%a(IB+0) = MAT_EZ(4*(IA-1)+1)
       mg_BZPMat(ILEV)%a(IB+1) = MAT_EZ(4*(IA-1)+2)
       mg_BZPMat(ILEV)%a(IB+2) = MAT_EZ(4*(IA-1)+3)
       mg_BZPMat(ILEV)%a(IB+3) = MAT_EZ(4*(IA-1)+4)

       BP_LDC(IDOF) = BP_LDC(IDOF) + 4
      END DO
     END DO

     pNEL = pNel + MGE013(ILEV)%SP(pJD)%nElems(2)
     DEALLOCATE (KCOL_E,KLD_E,MAT_EX,MAT_EY,MAT_EZ)

    END IF
   END IF
  END DO
 END IF
END DO

! IF (.TRUE.) THEN
! WRITE(myFile(1:10),'(A4,I2,A4)') "matB",myid,".txt"
!  OPEN(888,FILE=myFile)
!  DO I=1,mg_qlPMat(ILEV)%nu
!   write(888,'(I6,A1,50I12)') I,"|", (mg_qlPMat(ILEV)%ColA(IA),IA=mg_qlPMat(ILEV)%LdA(I),mg_qlPMat(ILEV)%LdA(I+1)-1)
!   write(888,'(I6,A1,50D12.4)') I,"|", (mg_BXPMat(ILEV)%a(IA),IA=mg_qlPMat(ILEV)%LdA(I),mg_qlPMat(ILEV)%LdA(I+1)-1)
!   write(888,'(I6,A1,50D12.4)') I,"|", (mg_BYPMat(ILEV)%a(IA),IA=mg_qlPMat(ILEV)%LdA(I),mg_qlPMat(ILEV)%LdA(I+1)-1)
!   write(888,'(I6,A1,50D12.4)') I,"|", (mg_BZPMat(ILEV)%a(IA),IA=mg_qlPMat(ILEV)%LdA(I),mg_qlPMat(ILEV)%LdA(I+1)-1)
!   write(888,"(200('-'))")
!  END DO
!  CLOSE(888)
! END IF

! DO pID=1,subnodes
!   IF (MGE013(ILEV)%ST(pID)%Num.GT.0) THEN
!    WRITE(*,'(A1,7I8,A1)') "-",myid,pID,MGE013(ILEV)%SP(pID)%nEntries(:),&
!                           MGE013(ILEV)%SP(pID)%nElems(:),pNEl,"-"
!   END IF
! END DO

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
!  WRITE(*,'(A,I2,A)') "Hello, ",myid, " is done!"
!  pause

END



SUBROUTINE GetParPressure(P,PP)
USE PP3D_MPI
USE def_feat, ONLY: ILEV

IMPLICIT NONE
INTEGER I,pID,pJD,nSize,II,JJ
REAL*8 P(*),PP(*)
CHARACTER*10 myFile

IF (myid.ne.MASTER) THEN

DO pID=1,subnodes
 IF (myid.eq.pID) THEN
  DO pJD=1,subnodes
   IF (myid.NE.pJD.AND.MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSize = MGE013(ILEV)%SP(pJD)%nElems(2)
     !write(*,*) "r", myid, pJD, nsize     
     CALL RECVD_myMPI(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,pJD)
     DO I=1,nSize
      II = 4*(I-1)+1
      JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(2,I)-1)+1
      PP(JJ+0) = MGE013(ILEV)%SP(pJD)%RDVect(II+0)
      PP(JJ+1) = MGE013(ILEV)%SP(pJD)%RDVect(II+1)
      PP(JJ+2) = MGE013(ILEV)%SP(pJD)%RDVect(II+2)
      PP(JJ+3) = MGE013(ILEV)%SP(pJD)%RDVect(II+3)
     END DO

   END IF
  END DO
 ELSE
  IF (MGE013(ILEV)%SP(pID)%Num.GT.0) THEN
   nSize = MGE013(ILEV)%SP(pID)%nElems(1)
   DO I=1,nSize
    II = 4*(I-1)+1
    JJ = 4*(MGE013(ILEV)%SP(pID)%VertLink(1,I)-1)+1
    MGE013(ILEV)%SP(pID)%SDVect(II+0) = P(JJ+0)
    MGE013(ILEV)%SP(pID)%SDVect(II+1) = P(JJ+1)
    MGE013(ILEV)%SP(pID)%SDVect(II+2) = P(JJ+2)
    MGE013(ILEV)%SP(pID)%SDVect(II+3) = P(JJ+3)
   END DO
   !write(*,*) "s", myid, pID, nsize
   CALL SENDD_myMPI(MGE013(ILEV)%SP(pID)%SDVect,4*nSIZE,pID)
  END IF
 END IF
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)


! JJ = 0
! DO pID=1,subnodes
!  JJ = JJ + MGE013(ILEV)%SP(pID)%nElems(2)
! END DO
! ! 
! WRITE(myFile(1:10),'(A4,I2,A4)') "file",myid,".txt"
! OPEN(987,FILE=myFile)
!  DO I=1,JJ
!  WRITE(987,'(I8,4D12.4)') 4*(i-1)+1,PP((4*(I-1)+1):(4*(I-1)+4))
!  END DO
! CLOSE(987)
! PAUSE

END IF

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013Sum(FX) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV

REAL*8  FX(*)
INTEGER I,pID,pJD,nSIZE,nEIGH

IF (myid.ne.MASTER) THEN

 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num

     DO I=1,nSIZE
      MGE013(ILEV)%ST(pJD)%SDVect(I) = FX(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
     END DO

     CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,pID)

    END IF
   END DO
  ELSE
   DO pJD=1,subnodes
     IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
      nSIZE = MGE013(ILEV)%ST(pJD)%Num

      CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,pJD)

     END IF
   END DO
  END IF
 END DO

 DO pJD=1,subnodes
   IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num

     DO I=1,nSIZE
       FX(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(I)
     END DO

   END IF
 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E013Sum
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWSum(FX) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel

REAL*8  FX(*)
INTEGER I,pID,pJD,nSIZE,nEIGH
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3

LEQ = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
LEQ1 =0
LEQ2 =LEQ
LEQ3 =2*LEQ

IF (myid.ne.MASTER) THEN

 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

!      write(*,'(10I10)')myid,leq1,leq2,leq3,meq1,meq2,meq3
!      pause
     DO I=1,nSIZE
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = FX(LEQ1+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = FX(LEQ2+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = FX(LEQ3+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
     END DO

     CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,pID)

    END IF
   END DO
  ELSE
   DO pJD=1,subnodes
     IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
      nSIZE = MGE013(ILEV)%ST(pJD)%Num

      CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,pJD)

     END IF
   END DO
  END IF
 END DO

 DO pJD=1,subnodes
   IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     DO I=1,nSIZE
       FX(LEQ1+MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX(LEQ1+MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ1+I)
       FX(LEQ2+MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX(LEQ2+MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ2+I)
       FX(LEQ3+MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX(LEQ3+MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ3+I)
     END DO

   END IF
 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E013UVWSum
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWMAT(A11,A22,A33,KLDA,NU) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV

REAL*8 A11(*),A22(*),A33(*)
INTEGER KLDA(*),NU
INTEGER I,J
INTEGER pID,pJD,nSIZE
INTEGER MEQ,MEQ1,MEQ2,MEQ3

IF (myid.ne.MASTER) THEN

 DO I=1,NU
  MGE013(ILEV)%UE11(I)=A11(KLDA(I))
  MGE013(ILEV)%UE22(I)=A22(KLDA(I))
  MGE013(ILEV)%UE33(I)=A33(KLDA(I))
 ENDDO

 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     DO I=1,nSIZE
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = MGE013(ILEV)%UE11(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = MGE013(ILEV)%UE22(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = MGE013(ILEV)%UE33(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
     END DO

     CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,pID)

    END IF
   END DO
  ELSE
   DO pJD=1,subnodes
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,pJD)

    END IF

   END DO
  END IF
 END DO

 DO pJD=1,subnodes

   IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     DO I=1,nSIZE
       MGE013(ILEV)%UE11(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       MGE013(ILEV)%UE11(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) + MGE013(ILEV)%ST(pJD)%RDVect(MEQ1+I)
       MGE013(ILEV)%UE22(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       MGE013(ILEV)%UE22(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) + MGE013(ILEV)%ST(pJD)%RDVect(MEQ2+I)
       MGE013(ILEV)%UE33(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       MGE013(ILEV)%UE33(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) + MGE013(ILEV)%ST(pJD)%RDVect(MEQ3+I)
     END DO

   END IF

 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E013UVWMAT
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013MAT(A,KLDA,NU) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV

REAL*8 A(*)
INTEGER KLDA(*),NU
INTEGER I,J
INTEGER pID,pJD,nSIZE

IF (myid.ne.MASTER) THEN

 DO I=1,NU
  MGE013(ILEV)%UE(I)=A(KLDA(I))
 ENDDO

 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num

     DO I=1,nSIZE
      MGE013(ILEV)%ST(pJD)%SDVect(I) = MGE013(ILEV)%UE(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
     END DO

     CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,pID)

    END IF
   END DO
  ELSE
   DO pJD=1,subnodes
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,pJD)

    END IF

   END DO
  END IF
 END DO

 DO pJD=1,subnodes

   IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     nSIZE = MGE013(ILEV)%ST(pJD)%Num

     DO I=1,nSIZE
       MGE013(ILEV)%UE(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       MGE013(ILEV)%UE(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) + MGE013(ILEV)%ST(pJD)%RDVect(I)
     END DO

   END IF

 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E013MAT

SUBROUTINE E012DISTR_L1(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J

!  WRITE(*,*) myid, "is here", ILEV
 IF (myid.NE.0) THEN
  pN = N
  CALL SENDI_myMPI(pN,0)
  CALL SENDD_myMPI(D,4*N,0)
 ELSE
  DO pID=1,subnodes
   CALL RECVI_myMPI(pN ,pID)
   CALL RECVD_myMPI(MGE013(ILEV)%CRSVect,4*pN,pID)

   DO I=1,pN
    J = coarse%pELEMLINK(pID,I)
    D(4*(J-1)+1) = MGE013(ILEV)%CRSVect(4*(I-1)+1)
    D(4*(J-1)+2) = MGE013(ILEV)%CRSVect(4*(I-1)+2)
    D(4*(J-1)+3) = MGE013(ILEV)%CRSVect(4*(I-1)+3)
    D(4*(J-1)+4) = MGE013(ILEV)%CRSVect(4*(I-1)+4)
   END DO

  END DO

!   OPEN (885,FILE ='gfgfgf.txt')
!   DO I=1,N
!   WRITE(885,*) I,D(I)
!   END DO
!   CLOSE(885)
 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E012DISTR_L1

SUBROUTINE E012GATHR_L1(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J 

!  WRITE(*,*) myid, "is here", ILEV
 IF (myid.NE.0) THEN
  pN = N
  CALL SENDI_myMPI(pN,0)
  CALL RECVD_myMPI(D,4*N,0)
 ELSE
  DO pID=1,subnodes
   CALL RECVI_myMPI(pN ,pID)

   DO I=1,pN
    J = coarse%pELEMLINK(pID,I)
    MGE013(ILEV)%CRSVect(4*(I-1)+1) = D(4*(J-1)+1)
    MGE013(ILEV)%CRSVect(4*(I-1)+2) = D(4*(J-1)+2)
    MGE013(ILEV)%CRSVect(4*(I-1)+3) = D(4*(J-1)+3)
    MGE013(ILEV)%CRSVect(4*(I-1)+4) = D(4*(J-1)+4)
   END DO
   CALL SENDD_myMPI(MGE013(ILEV)%CRSVect,4*pN,pID)
  END DO

!   OPEN (885,FILE ='gfgfgf.txt')
!   DO I=1,N
!   WRITE(885,*) I,D(I)
!   END DO
!   CLOSE(885)
 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E012GATHR_L1

SUBROUTINE E012DISTR_L2(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J,JJ,II,iu,ppN

!  WRITE(*,*) myid, "is here", ILEV
 IF (myid.NE.0) THEN
  pN = N
  CALL SENDI_myMPI(pN,0)
  CALL SENDD_myMPI(D,4*N,0)
 ELSE

  ppn = N

  DO pID=1,subnodes
   CALL RECVI_myMPI(pN ,pID)
   CALL RECVD_myMPI(MGE013(ILEV)%CRSVect,4*pN,pID)

   DO II=1,pN/8
    JJ = coarse%pELEMLINK(pID,II)

    I = II
    J = JJ
    D(4*(J-1)+1) = MGE013(ILEV)%CRSVect(4*(I-1)+1)
    D(4*(J-1)+2) = MGE013(ILEV)%CRSVect(4*(I-1)+2)
    D(4*(J-1)+3) = MGE013(ILEV)%CRSVect(4*(I-1)+3)
    D(4*(J-1)+4) = MGE013(ILEV)%CRSVect(4*(I-1)+4)

    DO iu = 1,7

     J = ppN/8 + (JJ-1)*7 + iu
     I =  pN/8 + (II-1)*7 + iu

     D(4*(J-1)+1) = MGE013(ILEV)%CRSVect(4*(I-1)+1)
     D(4*(J-1)+2) = MGE013(ILEV)%CRSVect(4*(I-1)+2)
     D(4*(J-1)+3) = MGE013(ILEV)%CRSVect(4*(I-1)+3)
     D(4*(J-1)+4) = MGE013(ILEV)%CRSVect(4*(I-1)+4)
    END DO

   END DO

  END DO

!   OPEN (885,FILE ='gfgfgf.txt')
!   DO I=1,N
!   WRITE(885,*) I,D(I)
!   END DO
!   CLOSE(885)
 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E012DISTR_L2


SUBROUTINE E012GATHR_L2(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J,II,JJ,iu,ppN

!  WRITE(*,*) myid, "is here", ILEV
 IF (myid.NE.0) THEN
  pN = N
  CALL SENDI_myMPI(pN,0)
  CALL RECVD_myMPI(D,4*N,0)
 ELSE
  DO pID=1,subnodes

   ppN = N

   CALL RECVI_myMPI(pN ,pID)

   DO II=1,pN/8
    JJ = coarse%pELEMLINK(pID,II)

    I = II
    J = JJ

    MGE013(ILEV)%CRSVect(4*(I-1)+1) = D(4*(J-1)+1)
    MGE013(ILEV)%CRSVect(4*(I-1)+2) = D(4*(J-1)+2)
    MGE013(ILEV)%CRSVect(4*(I-1)+3) = D(4*(J-1)+3)
    MGE013(ILEV)%CRSVect(4*(I-1)+4) = D(4*(J-1)+4)

    DO iu = 1,7

     J = ppN/8 + (JJ-1)*7 + iu
     I = pN/8 + (II-1)*7 + iu

    MGE013(ILEV)%CRSVect(4*(I-1)+1) = D(4*(J-1)+1)
    MGE013(ILEV)%CRSVect(4*(I-1)+2) = D(4*(J-1)+2)
    MGE013(ILEV)%CRSVect(4*(I-1)+3) = D(4*(J-1)+3)
    MGE013(ILEV)%CRSVect(4*(I-1)+4) = D(4*(J-1)+4)
    END DO

   END DO
   CALL SENDD_myMPI(MGE013(ILEV)%CRSVect,4*pN,pID)
  END DO

!   OPEN (885,FILE ='gfgfgf.txt')
!   DO I=1,N
!   WRITE(885,*) I,D(I)
!   END DO
!   CLOSE(885)
!   pause
 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E012GATHR_L2

SUBROUTINE E012DISTR_L3(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J,JJ,II,III,JJJ,iu,ppN,iw

!  WRITE(*,*) myid, "is here", ILEV
 IF (myid.NE.0) THEN
  pN = N
  CALL SENDI_myMPI(pN,0)
  CALL SENDD_myMPI(D,4*N,0)
 ELSE

  ppn = N

  DO pID=1,subnodes
   CALL RECVI_myMPI(pN ,pID)
   CALL RECVD_myMPI(MGE013(ILEV)%CRSVect,4*pN,pID)

   DO II=1,pN/64
    JJ = coarse%pELEMLINK(pID,II)

    I = II
    J = JJ
    D(4*(J-1)+1) = MGE013(ILEV)%CRSVect(4*(I-1)+1)
    D(4*(J-1)+2) = MGE013(ILEV)%CRSVect(4*(I-1)+2)
    D(4*(J-1)+3) = MGE013(ILEV)%CRSVect(4*(I-1)+3)
    D(4*(J-1)+4) = MGE013(ILEV)%CRSVect(4*(I-1)+4)

    DO iw = 1,7
     J = ppN/8 + (JJ-1)*7 + iw
     I =  pN/8 + (II-1)*7 + iw

     D(4*(J-1)+1) = MGE013(ILEV)%CRSVect(4*(I-1)+1)
     D(4*(J-1)+2) = MGE013(ILEV)%CRSVect(4*(I-1)+2)
     D(4*(J-1)+3) = MGE013(ILEV)%CRSVect(4*(I-1)+3)
     D(4*(J-1)+4) = MGE013(ILEV)%CRSVect(4*(I-1)+4)
    END DO

    DO iu = 1,7

     JJJ = ppN/64 + (JJ-1)*7 + iu
     III =  pN/64 + (II-1)*7 + iu
     I = III
     J = JJJ

     D(4*(J-1)+1) = MGE013(ILEV)%CRSVect(4*(I-1)+1)
     D(4*(J-1)+2) = MGE013(ILEV)%CRSVect(4*(I-1)+2)
     D(4*(J-1)+3) = MGE013(ILEV)%CRSVect(4*(I-1)+3)
     D(4*(J-1)+4) = MGE013(ILEV)%CRSVect(4*(I-1)+4)

     DO iw = 1,7

      J = ppN/8 + (JJJ-1)*7 + iw
      I =  pN/8 + (III-1)*7 + iw
 
      D(4*(J-1)+1) = MGE013(ILEV)%CRSVect(4*(I-1)+1)
      D(4*(J-1)+2) = MGE013(ILEV)%CRSVect(4*(I-1)+2)
      D(4*(J-1)+3) = MGE013(ILEV)%CRSVect(4*(I-1)+3)
      D(4*(J-1)+4) = MGE013(ILEV)%CRSVect(4*(I-1)+4)
     END DO

    END DO

   END DO

  END DO

!   OPEN (885,FILE ='gfgfgf.txt')
!   DO I=1,N
!   WRITE(885,*) I,D(I)
!   END DO
!   CLOSE(885)
 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E012DISTR_L3


SUBROUTINE E012GATHR_L3(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J,II,JJ,III,JJJ,iu,ppN,iw

!  WRITE(*,*) myid, "is here", ILEV
 IF (myid.NE.0) THEN
  pN = N
  CALL SENDI_myMPI(pN,0)
  CALL RECVD_myMPI(D,4*N,0)
 ELSE
  DO pID=1,subnodes

   ppN = N

   CALL RECVI_myMPI(pN ,pID)

   DO II=1,pN/64
    JJ = coarse%pELEMLINK(pID,II)

    I = II
    J = JJ

    MGE013(ILEV)%CRSVect(4*(I-1)+1) = D(4*(J-1)+1)
    MGE013(ILEV)%CRSVect(4*(I-1)+2) = D(4*(J-1)+2)
    MGE013(ILEV)%CRSVect(4*(I-1)+3) = D(4*(J-1)+3)
    MGE013(ILEV)%CRSVect(4*(I-1)+4) = D(4*(J-1)+4)

    DO iw = 1,7

     J = ppN/8 + (JJ-1)*7 + iw
     I = pN/8 + (II-1)*7 + iw

     MGE013(ILEV)%CRSVect(4*(I-1)+1) =D(4*(J-1)+1)
     MGE013(ILEV)%CRSVect(4*(I-1)+2) =D(4*(J-1)+2)
     MGE013(ILEV)%CRSVect(4*(I-1)+3) =D(4*(J-1)+3)
     MGE013(ILEV)%CRSVect(4*(I-1)+4) =D(4*(J-1)+4)
    END DO

    DO iu = 1,7

     JJJ = ppN/64 + (JJ-1)*7 + iu
     III = pN/64 + (II-1)*7 + iu
     I = III
     J = JJJ

    MGE013(ILEV)%CRSVect(4*(I-1)+1) = D(4*(J-1)+1)
    MGE013(ILEV)%CRSVect(4*(I-1)+2) = D(4*(J-1)+2)
    MGE013(ILEV)%CRSVect(4*(I-1)+3) = D(4*(J-1)+3)
    MGE013(ILEV)%CRSVect(4*(I-1)+4) = D(4*(J-1)+4)

    DO iw = 1,7

     J = ppN/8 + (JJJ-1)*7 + iw
     I = pN/8 + (III-1)*7 + iw

     MGE013(ILEV)%CRSVect(4*(I-1)+1) =D(4*(J-1)+1)
     MGE013(ILEV)%CRSVect(4*(I-1)+2) =D(4*(J-1)+2)
     MGE013(ILEV)%CRSVect(4*(I-1)+3) =D(4*(J-1)+3)
     MGE013(ILEV)%CRSVect(4*(I-1)+4) =D(4*(J-1)+4)
    END DO

    END DO

   END DO
   CALL SENDD_myMPI(MGE013(ILEV)%CRSVect,4*pN,pID)
  END DO

!   OPEN (885,FILE ='gfgfgf.txt')
!   DO I=1,N
!   WRITE(885,*) I,D(I)
!   END DO
!   CLOSE(885)
!   pause
 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E012GATHR_L3

SUBROUTINE E013SendK(iP,jP,Num)
USE PP3D_MPI
IMPLICIT NONE
Integer iP,jP,Num

IF (myid.eq.iP) CALL SENDI_myMPI(Num,jP)
IF (myid.eq.jP) CALL RECVI_myMPI(Num ,iP)

END SUBROUTINE E013SendK

SUBROUTINE CommBarrier()
USE PP3D_MPI
IMPLICIT NONE

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE CommBarrier
