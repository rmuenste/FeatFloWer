!
! ----------------------------------------------
!
SUBROUTINE CommunicateInterface(lIfc,gIfc)
USE PP3D_MPI
USE var_QuadScalar, ONLY : tTetInterface

IMPLICIT NONE
TYPE(tTetInterface) lIfc,gIfc
INTEGER nG,nT,nnG,nnT,pID

IF (myid.ne.0) THEN
 nG = lIfc%nG
 nT = lIfc%nT
 CALL SENDI_myMPI(nG,0)
 CALL SENDI_myMPI(nT,0)
 CALL RECVI_myMPI(nnG,0)
 CALL RECVI_myMPI(nnT,0)
ELSE
 nnG = 1
 nnT = 0
 DO pID=1,subnodes
  CALL RECVI_myMPI(nG,pID)
  CALL RECVI_myMPI(nT,pID)
  nnG = nnG + nG
  nnT = nnT + nT
 END DO

 DO pID=1,subnodes
  CALL SENDI_myMPI(nnG,pID)
  CALL SENDI_myMPI(nnT,pID)
 END DO

END IF

IF (ALLOCATED(gIfc%X)) DEALLOCATE(gIfc%X)
IF (ALLOCATED(gIfc%L)) DEALLOCATE(gIfc%L)
gIFC%nG = nnG
gIFC%nT = nnT
ALLOCATE(gIFC%X(3,3*nnT))
ALLOCATE(gIFC%L(nnG))

IF (myid.ne.0) THEN
 nG = lIfc%nG
 nT = lIfc%nT
 CALL SENDI_myMPI(nG,0)
 CALL SENDK_myMPI(lIfc%L,nG,0)

 CALL SENDI_myMPI(nT,0)
 CALL SENDD_myMPI(lIfc%X,9*nT,0)

 CALL RECVK_myMPI(gIfc%L,nnG,0)
 CALL RECVD_myMPI(gIfc%X,9*gIFC%nT,0)

ELSE
 nnG = 1
 nnT = 0
 DO pID=1,subnodes
  CALL RECVI_myMPI(nG,pID)
  CALL RECVK_myMPI(gIfc%L(nnG+1:nnG+nG),nG,pID)
  nnG = nnG + nG

  CALL RECVI_myMPI(nT,pID)
  CALL RECVD_myMPI(gIfc%X(:,nnT+1:nnT+3*nT),9*nT,pID)
  nnT = nnT + 3*nT
 END DO

 DO pID=1,subnodes
  CALL SENDK_myMPI(gIfc%L,nnG,pID)
  CALL SENDD_myMPI(gIfc%X,9*gIFC%nT,pID)
 END DO

END IF

END
!
! ----------------------------------------------
!
subroutine E010_CollectCoarseVector(X,nP)
USE PP3D_MPI
implicit none
REAL*8 X(*)
integer nP
integer nnP,i,j,pID
real*8, allocatable :: buffer(:)

if (myid.ne.0) then

  nnP = nP  
  call SENDI_myMPI(nnP,0)
  CALL SENDD_myMPI(X,nnP,0)
  
else

 allocate(buffer(np))
 do pID=1,subnodes 
  call RECVI_myMPI(nnP,pID)
  call RECVD_myMPI(buffer,nnP,pID)
  DO i=1,nnP
   j = coarse%pELEMLINK(pID,I)
   X(j) = buffer(i)
  end do
 end do
 deallocate(buffer)
 
end if


END subroutine E010_CollectCoarseVector
!
! ----------------------------------------------
!
subroutine SendPressBCElemsToCoarse(iVect,nP)
USE PP3D_MPI
implicit none
integer iVect(*),nP
integer, allocatable :: iSend(:)
integer nnp,i,iP

if (myid.ne.0) then
  call RECVI_myMPI(nnP,0)
  allocate(iSend(nnp))
  iSend = 0
  do i=1,np
   if (iVect(i).ne.0) iSend(coarse%myELEMLINK(i)) = 1   
  end do
else
 do iP=1,subnodes 
  call SENDI_myMPI(nP,IP)
 end do
 allocate(iSend(np))
end if

if (myid.ne.0) then
  call SENDK_myMPI(iSend,nnP,0)
else
 nnP = nP
 do iP=1,subnodes 
  call RECVK_myMPI(iSend,nnP,iP)
  do i=1,nnP
   if (iSend(i).eq.1) iVect(i) = 1
  end do
 end do
 do i=1,np
!    if (iVect(i).ne.0) write(*,*) i  
 end do
end if

deallocate(iSend)

end subroutine SendPressBCElemsToCoarse
!
! ----------------------------------------------
!
SUBROUTINE  Comm_Solution(d1,d2,d3,d4,n)
USE PP3D_MPI
IMPLICIT NONE
REAL*8 d1(*),d2(*),d3(*),d4(*)
INTEGER i,n

DO i=1,n
 d4(i) = 1d0
END DO

CALL E013Sum(d1)
CALL E013Sum(d2)
CALL E013Sum(d3)
CALL E013Sum(d4)

DO i=1,n
 d1(i) = d1(i)/d4(i)
 d2(i) = d2(i)/d4(i)
 d3(i) = d3(i)/d4(i)
END DO

END SUBROUTINE  Comm_Solution
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE Create_GlobalP1CommNumbering(lScalar,lPScalar,DCORVG,KVERT,KEDGE,KAREA,NVT,NET,NAT,NEL)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY : mg_kVector,HlobalParallelList1,HlobalParallelList2,&
    HlobalParallelBufferIn,HlobalParallelBufferOut,HlobalNList,HlobalNBuffer,&
    HlobalParallelList3,TLinScalar,TParLinScalar
implicit none
TYPE(TLinScalar), INTENT(INOUT), TARGET :: lScalar
TYPE(TParLinScalar), INTENT(INOUT), TARGET ::  lPScalar
REAL*8  DCORVG(3,*)
INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)
INTEGER NAT,NEL,NVT,NET,ndof

INTEGER NeighE(2,12),NeighA(4,6)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
REAL*8 PX,PY,PZ
REAL*8 P1X,P1Y,P1Z,P2X,P2Y,P2Z
REAL*8 Point1(3),Point2(3),dist

INTEGER,ALLOCATABLE :: iAux(:),offset(:)
INTEGER iVertices,jVertices,pID,nLenght,nPLength
INTEGER i,j,k,iSend,iRecv
 
 IF (.not.Allocated(HlobalParallelList1)) ALLOCATE(HlobalParallelList1(NLMIN:NLMAX))
 IF (.not.Allocated(HlobalParallelList2)) ALLOCATE(HlobalParallelList2(NLMIN:NLMAX))
 IF (.not.Allocated(HlobalParallelList3)) ALLOCATE(HlobalParallelList3(NLMIN:NLMAX))
 IF (.not.Allocated(HlobalNBuffer)) ALLOCATE(HlobalNBuffer(NLMIN:NLMAX))
 IF (.not.Allocated(HlobalNList)) ALLOCATE(HlobalNList(NLMIN:NLMAX))

 ndof = nel
 allocate(iaux(ndof))
 iaux = 0

 DO pID=1,subnodes
   DO I=1,MGE013(ILEV)%SP(pID)%nElems(1)
    j=MGE013(ILEV)%SP(pID)%VertLink(1,I)
    iaux(j) = 1
   END DO
 END DO

 iVertices = 0
 DO i=1,ndof
  if (iaux(i).eq.1) iVertices  = iVertices + 1
 END DO
 
 HlobalNList(ILEV)   = iVertices
 IF (.not.Allocated(HlobalParallelList1(ilev)%x)) ALLOCATE(HlobalParallelList1(ilev)%x(iVertices))
 IF (.not.Allocated(HlobalParallelList2(ilev)%x)) ALLOCATE(HlobalParallelList2(ilev)%x(iVertices))

 ALLOCATE(offset(subnodes+1))

 offset(1) = 0
 DO pID=1,subnodes
  iSend = 0
  if (myid.eq.pID) iSend = iVertices
  CALL MPI_Allreduce(iSend,iRecv,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SUBS,iERR)
  iSend = iRecv
  offset(pID+1) = offset(pID) + iSend 
 END DO

 nLenght=offset(subnodes+1)
 HlobalNBuffer(ILEV) = nLenght
 
 k = 0
 lScalar%valP(ilev)%x= 0
 
 DO i=1,ndof
  if (iaux(i).eq.1) THEN
   k = k + 1
   HlobalParallelList1(ilev)%x(k) = offset(myid)+k
   HlobalParallelList2(ilev)%x(k) = i
   lScalar%valP(ilev)%x(i) = DBLE(offset(myid)+k) !DBLE((myid))!
  end if
 END DO

 nPLength = 0
 DO pID=1,subnodes
  nPLength = nPLength + MGE013(ILEV)%SP(pID)%nElems(2)
 END DO
 
 CALL GetParPresEntries(lScalar%valP(ilev)%x,lPScalar%val)
 IF (.not.Allocated(HlobalParallelList3(ilev)%x)) ALLOCATE(HlobalParallelList3(ilev)%x(nPLength))
  
 k = 0
 do i=1,nPLength
   k = k + 1
   HlobalParallelList3(ilev)%x(k) = INT(lPScalar%val(i))
 end do

 if (ilev.eq.nlmax) ALLOCATE(HlobalParallelBufferIn(4*nLenght))
 if (ilev.eq.nlmax) ALLOCATE(HlobalParallelBufferOut(4*nLenght))
! !  WRITE(*,*) offset

! if (ilev.eq.nlmax) then
!  write(*,*) myid, nLenght,nPLength,' : ', HlobalParallelList3(ilev)%x 
! END IF
! if (ilev.eq.nlmax) pause
 
 deallocate(iaux,offset)

END SUBROUTINE Create_GlobalP1CommNumbering
!
! ----------------------------------------------
!
SUBROUTINE Create_GlobalQ2CommNumbering(DCORVG,KVERT,KEDGE,KAREA,NVT,NET,NAT,NEL)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY :mg_kVector,GlobalParallelList1,GlobalParallelList2,&
    GlobalParallelBufferIn,GlobalParallelBufferOut,GlobalNList,GlobalNBuffer

REAL*8  DCORVG(3,*)
INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)
INTEGER NAT,NEL,NVT,NET,ndof

INTEGER NeighE(2,12),NeighA(4,6)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
REAL*8 PX,PY,PZ
REAL*8 P1X,P1Y,P1Z,P2X,P2Y,P2Z
REAL*8 Point1(3),Point2(3),dist

INTEGER,ALLOCATABLE :: iAux(:),offset(:)
REAL*8 ,ALLOCATABLE :: dAux(:,:),dAuxO(:,:)
INTEGER iVertices,jVertices,pID,nLenght

 IF (.not.Allocated(GlobalParallelList1)) ALLOCATE(GlobalParallelList1(NLMIN:NLMAX))
 IF (.not.Allocated(GlobalParallelList2)) ALLOCATE(GlobalParallelList2(NLMIN:NLMAX))
 IF (.not.Allocated(GlobalNBuffer)) ALLOCATE(GlobalNBuffer(NLMIN:NLMAX))
 IF (.not.Allocated(GlobalNList)) ALLOCATE(GlobalNList(NLMIN:NLMAX))


 ndof = nvt+net+nat+nel
 allocate(iaux(ndof))

 iaux = 2**10

 DO pID=1,subnodes
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    j = MGE013(ILEV)%ST(pID)%VertLink(1,i)
    iaux(j) = min(iaux(j),pID)
   END DO
 END DO

 iVertices = 0
 DO i=1,ndof
  if (iaux(i).le.subnodes) iVertices  = iVertices + 1
 END DO
 IF (.not.Allocated(GlobalParallelList1(ilev)%x)) ALLOCATE(GlobalParallelList1(ilev)%x(iVertices))
 IF (.not.Allocated(GlobalParallelList2(ilev)%x)) ALLOCATE(GlobalParallelList2(ilev)%x(iVertices))
 
 jVertices = 0
 DO i=1,ndof
  if (iaux(i).le.subnodes.and.iaux(i).gt.myid) jVertices  = jVertices + 1
 END DO


ALLOCATE(offset(subnodes))

offset(1) = 0
DO pID=1,subnodes-1
 iSend = 0
 if (myid.eq.pID) iSend = jVertices
 CALL MPI_Allreduce(iSend,iRecv,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SUBS,iERR)
 iSend = iRecv
 offset(pID+1) = offset(pID) + iSend 
END DO

 nLenght=offset(subnodes)
 ALLOCATE(dAux(3,nLenght),dAuxO(3,nLenght))

 daux = 0d0
 iaux = 2**10
 k = offset(myid)
 DO pID=1,subnodes
   DO I=1,MGE013(ILEV)%ST(pID)%Num
    j = MGE013(ILEV)%ST(pID)%VertLink(1,i)
    if (iaux(j).gt.subnodes.and.pID.gt.myid) then
     k = k + 1
     daux(1,k) = dcorvg(1,j)
     daux(2,k) = dcorvg(2,j)
     daux(3,k) = dcorvg(3,j)
    end if
    iaux(j) = min(iaux(j),pID)
   END DO
 END DO

 CALL MPI_Allreduce(dAux,dAuxO,3*nLenght,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SUBS,iERR)

!  if (myid.eq.5) write(*,*)iaux

k = 0
 DO i=1,ndof
  inode = iaux(i)
  IF (inode.le.subnodes) then
   Point1 = dcorvg(:,i)
   IF (inode.lt.myid) THEN
    DO j=offset(inode)+1,offset(inode+1)
     Point2 = dauxO(:,j)
     dist = SQRT((Point1(1)-Point2(1))**2d0+(Point1(2)-Point2(2))**2d0+(Point1(3)-Point2(3))**2d0)
     IF (dist.lt.dEpsPrec) THEN
      k = k + 1
      GlobalParallelList1(ilev)%x(k) = j
      GlobalParallelList2(ilev)%x(k) = i
!        IF (jnode.ne.myid.and.jnode.ne.inode) WRITE(*,*) 'Found ! ',i,myid,inode,jnode
      GOTO 1 
     END IF
    END DO
   ELSE
    DO j=offset(myid)+1,offset(myid+1)
     Point2 = dauxO(:,j)
     dist = SQRT((Point1(1)-Point2(1))**2d0+(Point1(2)-Point2(2))**2d0+(Point1(3)-Point2(3))**2d0)
     IF (dist.lt.dEpsPrec) THEN
      k = k + 1
      GlobalParallelList1(ilev)%x(k) = j
      GlobalParallelList2(ilev)%x(k) = i
!        IF (jnode.ne.myid.and.jnode.ne.inode) WRITE(*,*) 'Found ! ',i,myid,inode,jnode
      GOTO 1 
     END IF
    END DO
   END IF   
!    DO jnode = 1,subnodes-1
!    DO j=offset(jnode)+1,offset(jnode+1)
!     Point2 = dauxO(:,j)
!     dist = SQRT((Point1(1)-Point2(1))**2d0+(Point1(2)-Point2(2))**2d0+(Point1(3)-Point2(3))**2d0)
!     IF (dist.lt.dEpsPrec) THEN
!       IF (jnode.ne.myid.and.jnode.ne.inode) WRITE(*,*) 'Found ! ',i,myid,inode,jnode
!      GOTO 1 
!     END IF
!    END DO
!    END DO
   WRITE(*,*) 'not found ... ',i,myid,inode
1  CONTINUE
  END IF
 END DO


 GlobalNBuffer(ILEV) = nLenght
 GlobalNList(ILEV)   = iVertices
 if (ilev.eq.nlmax) ALLOCATE(GlobalParallelBufferIn(3*nLenght))
 if (ilev.eq.nlmax) ALLOCATE(GlobalParallelBufferOut(3*nLenght))

 if (myid.eq.1) write(*,*) ilev,'offset: ',offset
 DEALLOCATE(iaux)
 DEALLOCATE(offset, daux, dauxO)

END SUBROUTINE Create_GlobalQ2CommNumbering
!
! ----------------------------------------------
!
SUBROUTINE ExchangeVelocitySolutionCoarseSub(U,V,W,dU,dV,dW,N)
USE PP3D_MPI

IMPLICIT NONE
REAL*8 U(*),V(*),W(*),dU(*),dV(*),dW(*)
INTEGER N
!!!!!!!!!!!!!!!!!
INTEGER i,j,pN,pID

IF (myid.eq.0) THEN

 DO pID=1,subnodes
  CALL RECVI_myMPI(pN ,pID)
  CALL RECVD_myMPI(dU,pN,pID)
  CALL RECVD_myMPI(dV,pN,pID)
  CALL RECVD_myMPI(dW,pN,pID)

  DO i=1,pN
   J = mQ2(pID)%x(I)
   U(J) = dU(i)
   V(J) = dV(i)
   W(J) = dW(i)
  END DO
 END DO

ELSE

 CALL SENDI_myMPI(N ,0)
 CALL SENDD_myMPI(U,N,0)
 CALL SENDD_myMPI(V,N,0)
 CALL SENDD_myMPI(W,N,0)

END IF

END SUBROUTINE ExchangeVelocitySolutionCoarseSub
!
! ----------------------------------------------
!
SUBROUTINE Create_GlobalNumbering_CC
USE PP3D_MPI
USE def_feat
USE var_QuadScalar, ONLY : GlobalNumbering,my_crs_e013_map
IMPLICIT NONE
INTEGER i,j,ndof,nnQ2,nnP0,nQ2,nP0,iQ2,iP0,pID
REAL*8, ALLOCATABLE :: D(:)
REAL*8 daux

ILEV = NLMIN
CALL SETLEV(2)

nnQ2  = NEL+NVT+NAT+NET
nnP0 =  NEL
ndof = 3*nnQ2 + 4*nnP0

nQ2=0
nP0=0
j = 0

IF (myid.eq.0) THEN

 DO pID=1,subnodes
  DO i=1,my_crs_e013_map(pid)%cc_ndof

   j = my_crs_e013_map(pid)%indE(i)
   my_crs_e013_map(pid)%dBuffer(i) = dble(j)
  END DO

  CALL sendD_myMPI(my_crs_e013_map(pid)%dBuffer,my_crs_e013_map(pid)%cc_ndof,pID)
 END DO

ELSE

  ALLOCATE (GlobalNumbering(ndof))
  ALLOCATE (D(ndof))

  CALL recvD_myMPI(D,ndof,0)
  GlobalNumbering = INT(D)
  DEALLOCATE (D)
END IF

END SUBROUTINE Create_GlobalNumbering_CC
!
! ----------------------------------------------
!
SUBROUTINE E013_Comm_Master(DCORVG,KVERT,KEDGE,KAREA,NVT,NET,NAT,NEL)

USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY :my_crs_e013_map

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
   IF ((ABS(P1X-P2X).LT.1D-5).AND.&
       (ABS(P1Y-P2Y).LT.1D-5).AND.&
       (ABS(P1Z-P2Z).LT.1D-5)) THEN
     iCoor(I) = j
     jAux = jAux + 1
     bFound = .TRUE.
    END IF
   END DO
!    IF (.NOT.BFOUND) WRITE(*,*) 'pair was not found!', myid,i
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
   IF (my_crs_e013_map(pid)%ind(i).gt.NVT+NET+NAT) THEN
    jaux = jaux + 4
   END IF
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

END SUBROUTINE E013_Comm_Master
!
! ----------------------------------------------
!
SUBROUTINE Create_GlobalNumbering(DCORVG,KVERT,KEDGE,KAREA,NVT,NET,NAT,NEL)

USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY :GlobalNumberingQ2,GlobalNumberingP1,myGlobalNumberingMap,myGlobal_ndof

IMPLICIT NONE

REAL*8  DCORVG(3,*)
INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)
INTEGER NAT,NEL,NVT,NET,ndof

INTEGER NeighE(2,12),NeighA(4,6)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
REAL*8 PX,PY,PZ
REAL*8 P1X,P1Y,P1Z,P2X,P2Y,P2Z,dist

INTEGER,ALLOCATABLE :: iCoor(:),iAux(:)
REAL*8, ALLOCATABLE :: dCoor(:,:),rCoor(:,:),D(:)

INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,jaux
INTEGER nnQ2,nnP0,pndof,cndof,pID,myINDEX
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
 CALL InitOctTree(rCoor,cndof)
 DO I=1,ndof
  CALL FindInOctTree(rCoor,cndof,dCoor(:,i),J,dist)
  IF (J.lt.0) then
   WRITE(*,*) I,"PROBLEM of Q2 dof assignement ..."
  end if
  IF (DIST.LT.DEpsPrec) THEN 
   iCoor(I) = j
   jAux = jAux + 1
  END IF
  CALL FindInPeriodicOctTree(rCoor,cndof,dCoor(:,i),J,dist,dPeriodicity)
  IF (DIST.LT.DEpsPrec) THEN 
   iCoor(I) = j
   jAux = jAux + 1
  END IF
!   P1X = dCoor(1,I)
!   P1Y = dCoor(2,I)
!   P1Z = dCoor(3,I)
!   bFound = .FALSE.
!   DO J=1,cndof
!    P2X = rCoor(1,J)
!    P2Y = rCoor(2,J)
!    P2Z = rCoor(3,J)
!     IF (((ABS(P1X-P2X).LT.DEpsPrec).OR.(ABS(ABS(P1X-P2X)-dPeriodicity(1)).LT.DEpsPrec)).AND.&
!         ((ABS(P1Y-P2Y).LT.DEpsPrec).OR.(ABS(ABS(P1Y-P2Y)-dPeriodicity(2)).LT.DEpsPrec)).AND.& 
!         ((ABS(P1Z-P2Z).LT.DEpsPrec).OR.(ABS(ABS(P1Z-P2Z)-dPeriodicity(3)).LT.DEpsPrec))) THEN
!      iCoor(I) = j
!      jAux = jAux + 1
!      bFound = .TRUE.
!     END IF
!    END DO
!    IF (.NOT.BFOUND) WRITE(*,*) 'pair was not found!', myid,i
  END DO
 CALL FreeOctTree()
 END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

IF (myid.eq.0) THEN
 ALLOCATE (myGlobalNumberingMap(subnodes))
 DO pID=1,subnodes
  CALL RECVI_myMPI(pNDOF,pID)
  myGlobalNumberingMap(pid)%ndof = pNDOF
  ALLOCATE(myGlobalNumberingMap(pid)%ind(pNDOF))
  CALL RECVK_myMPI(myGlobalNumberingMap(pid)%ind,pNdof,pID)
 END DO
ELSE
 CALL SENDI_myMPI(ndof,0)
 CALL SENDK_myMPI(iCoor,ndof,0)
END IF

IF (myid.eq.0) THEN
 DO pID=1,subnodes
  jaux = 0
  DO i=1,myGlobalNumberingMap(pid)%ndof
   jaux = jaux + 3
  END DO
  ALLOCATE(myGlobalNumberingMap(pid)%indQ2(jaux))
  ALLOCATE(myGlobalNumberingMap(pid)%dBufferQ2(jaux))
  myGlobalNumberingMap(pid)%ndof_Q2 = jaux

  jaux = 0
  DO i=1,myGlobalNumberingMap(pid)%ndof
   myINDEX = myGlobalNumberingMap(pid)%ind(i)
   myGlobalNumberingMap(pid)%indQ2(0*myGlobalNumberingMap(pid)%ndof + i) = 0*ndof + myINDEX
   myGlobalNumberingMap(pid)%indQ2(1*myGlobalNumberingMap(pid)%ndof + i) = 1*ndof + myINDEX
   myGlobalNumberingMap(pid)%indQ2(2*myGlobalNumberingMap(pid)%ndof + i) = 2*ndof + myINDEX
  END DO

 END DO
END IF

IF (myid.eq.0) THEN
 DO pID=1,subnodes
  jaux = 0
  DO i=1,myGlobalNumberingMap(pid)%ndof
   IF (myGlobalNumberingMap(pid)%ind(i).gt.NVT+NET+NAT) THEN
    jaux = jaux + 4
   END IF
  END DO
  ALLOCATE(myGlobalNumberingMap(pid)%indP1(jaux))
  ALLOCATE(myGlobalNumberingMap(pid)%dBufferP1(jaux))
  myGlobalNumberingMap(pid)%ndof_P1 = jaux

  jaux = 0
  DO i=1,myGlobalNumberingMap(pid)%ndof
   myINDEX = myGlobalNumberingMap(pid)%ind(i)
   IF (myINDEX.gt.NVT+NET+NAT) THEN
    jaux = jaux + 1
    myGlobalNumberingMap(pid)%indP1(4*(jaux-1)+1) = 4*(myINDEX-(NVT+NET+NAT)-1)+1
    myGlobalNumberingMap(pid)%indP1(4*(jaux-1)+2) = 4*(myINDEX-(NVT+NET+NAT)-1)+2
    myGlobalNumberingMap(pid)%indP1(4*(jaux-1)+3) = 4*(myINDEX-(NVT+NET+NAT)-1)+3
    myGlobalNumberingMap(pid)%indP1(4*(jaux-1)+4) = 4*(myINDEX-(NVT+NET+NAT)-1)+4
   END IF
  END DO

 END DO
END IF

DEALLOCATE (iAux)
DEALLOCATE (dCoor)
DEALLOCATE (iCoor)

nnQ2  = NEL+NVT+NAT+NET
nnP0 =  NEL

j = 0

ndof = 3*nnQ2
IF (myid.eq.0) THEN
 DO pID=1,subnodes
  DO i=1,myGlobalNumberingMap(pid)%ndof_Q2
   j = myGlobalNumberingMap(pid)%indQ2(i)
   myGlobalNumberingMap(pid)%dBufferQ2(i) = dble(j)
  END DO
  CALL sendD_myMPI(myGlobalNumberingMap(pid)%dBufferQ2,myGlobalNumberingMap(pid)%ndof_Q2,pID)
 END DO
ELSE
  ALLOCATE (GlobalNumberingQ2(ndof),D(ndof))
  CALL recvD_myMPI(D,ndof,0)
  GlobalNumberingQ2 = INT(D)
  DEALLOCATE (D)
END IF

ndof = 4*nnP0
IF (myid.eq.0) THEN
 DO pID=1,subnodes
  DO i=1,myGlobalNumberingMap(pid)%ndof_P1
   j = myGlobalNumberingMap(pid)%indP1(i)
   myGlobalNumberingMap(pid)%dBufferP1(i) = dble(j)
  END DO
!   write(*,*) myid, '<--', myGlobalNumberingMap(pid)%cc_ndof
  CALL sendD_myMPI(myGlobalNumberingMap(pid)%dBufferP1,myGlobalNumberingMap(pid)%ndof_P1,pID)
 END DO

ELSE

  ALLOCATE (GlobalNumberingP1(ndof),D(ndof))
!   write(*,*) myid, '-->', ndof  
  CALL recvD_myMPI(D,ndof,0)
  GlobalNumberingP1 = INT(D)
  DEALLOCATE (D)
!   WRITE(*,*) GlobalNumberingP1
END IF


IF (myid.eq.0) THEN
 ndof = nnQ2
 myGlobal_ndof = ndof
 DO pID=1,subnodes
  CALL sendI_myMPI(ndof,pID)
 END DO

ELSE
  CALL recvI_myMPI(ndof,0)
  myGlobal_ndof  = ndof
END IF

! pause

END SUBROUTINE Create_GlobalNumbering
!
! ----------------------------------------------
!
SUBROUTINE GetQ2CoarseMapper(DCORVG,KVERT,KEDGE,KAREA,NVT,NET,NAT,NEL)
USE PP3D_MPI

IMPLICIT NONE
REAL*8  DCORVG(3,*)
INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)
INTEGER NVT,NEL,NAT,NET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
INTEGER pID
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
REAL*8 PX,PY,PZ,DIST
INTEGER ppC(5),mppC(5)
INTEGER NeighE(2,12),NeighA(4,6)
LOGICAL bFound

integer :: ndof

REAL*8, ALLOCATABLE :: dCoor(:,:)
REAL*8, ALLOCATABLE :: pdCoor(:,:)

DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

ndof = NVT+NET+NAT+NEL

ALLOCATE (dCoor(3,1:ndof))
ALLOCATE (pdCoor(3,1:ndof))

DO I=1,NVT
 dcoor(1,I)=DCORVG(1,I)
 dcoor(2,I)=DCORVG(2,I)
 dcoor(3,I)=DCORVG(3,I)
END DO

k=1
DO i=1,nel
 DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
   ivt1 = kvert(NeighE(1,j),i)
   ivt2 = kvert(NeighE(2,j),i)
   PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
   PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
   PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
   dcoor(:,nvt+k) = [PX,PY,PZ]
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
   dcoor(:,nvt+net+k) = [PX,PY,PZ]
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
 dcoor(:,nvt+net+nat+i) = [PX,PY,PZ]
END DO

ALLOCATE (mQ2(subnodes))

IF (myid.eq.0) THEN

mppC(1) = 0
mppC(2) = NVT
mppC(3) = NVT+NET
mppC(4) = NVT+NET+NAT
mppC(5) = NVT+NET+NAT+NEL
!  coarse%pVERTLINK = 0
! 
 DO pID=1,subnodes
  ppC(1) = 0
  CALL RECVI_myMPI(ppC(2) ,pID)
  CALL RECVI_myMPI(ppC(3) ,pID)
  CALL RECVI_myMPI(ppC(4) ,pID)
  CALL RECVI_myMPI(ppC(5) ,pID)
  CALL RECVD_myMPI(pdCoor,3*ppC(5),pID)

  ALLOCATE(mQ2(pID)%x(ppC(5)))

!   coarse%pNVT(pID)=pNVT
! 
  DO k=1,4
   DO I=ppC(k)+1,ppC(k+1)
    bFound = .FALSE.
    DO J=mppC(k)+1,mppC(k+1)

     DIST = SQRT((pdcoor(1,I)-dcoor(1,J))**2d0 + &
                 (pdcoor(2,I)-dcoor(2,J))**2d0 + &
                 (pdcoor(3,I)-dcoor(3,J))**2d0)
     IF (DIST.LT.DEpsPrec) THEN

         mQ2(pID)%x(I)=J 
         bFound = .TRUE.
         GOTO 1

     END IF

    END DO
    IF (.NOT.bFound) WRITE(*,*) "problem with pid=",pID,ppC(k+1)-ppC(k),mppC(k+1)-mppC(k),i,j
1   CONTINUE
   END DO
  END DO

 END DO

ELSE
 CALL SENDI_myMPI(NVT ,0)
 CALL SENDI_myMPI(NVT+NET ,0)
 CALL SENDI_myMPI(NVT+NET+NAT ,0)
 CALL SENDI_myMPI(NVT+NET+NAT+NEL ,0)
 CALL SENDD_myMPI(dcoor,3*(NVT+NET+NAT+NEL),0)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE GetQ2CoarseMapper
!
! ----------------------------------------------
!
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
   IF (myid.NE.pJD.and.VerticeCommunicationScheme(pJD).ne.0) THEN
    CALL RECVI_myMPI(NNEL,pJD)
    !write(*,*) "r:",myid,PJD,NNEL    
    MGE013(ILEV)%SP(pJD)%nElems(2) = NNEL
   END IF
  END DO
 ELSE
  IF (VerticeCommunicationScheme(pID).ne.0) THEN
   NNEL = MGE013(ILEV)%SP(pID)%nElems(1)
   !write(*,*) "s:",myid,PID,NNEL  
   CALL SENDI_myMPI(NNEL ,pID)
  END IF
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

CALL MemoryPrint(1,'s','COMM0')

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

CALL MemoryPrint(1,'s','COMM1')

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
   IF (myid.NE.pJD.and.VerticeCommunicationScheme(pJD).gt.0) THEN
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
  IF (VerticeCommunicationScheme(pID).gt.0) THEN
   jAux = CoorSP(pID)%Num
   CALL SENDI_myMPI(jAux ,pID)
   IF (jAux.GT.0) THEN
    CALL SENDD_myMPI(CoorSP(pID)%dCoor,3*jAux,pID)
    CALL SENDK_myMPI(CoorSP(pID)%iCoor,  jAux,pID)
   END IF
  END IF
 END IF

END DO

CALL MemoryPrint(1,'s','COMM2')

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

CALL MemoryPrint(1,'s','COMM3')

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
    IF (((ABS(P1X-P2X).LT.DEpsPrec).OR.(ABS(ABS(P1X-P2X)-dPeriodicity(1)).LT.DEpsPrec)).AND.&
        ((ABS(P1Y-P2Y).LT.DEpsPrec).OR.(ABS(ABS(P1Y-P2Y)-dPeriodicity(2)).LT.DEpsPrec)).AND.&
        ((ABS(P1Z-P2Z).LT.DEpsPrec).OR.(ABS(ABS(P1Z-P2Z)-dPeriodicity(3)).LT.DEpsPrec))) THEN
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

CALL MemoryPrint(1,'s','COMM4')
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

! CALL MemoryPrint(1,'s','COMM5')

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

CALL MemoryPrint(1,'s','COMM-done')

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
USE OctTreeSearch

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
REAL*8 P1X,P1Y,P1Z,P2X,P2Y,P2Z,dist

IF (.NOT.ALLOCATED(MGE013)) ALLOCATE(MGE013(NLMIN:NLMAX))

ndof = NVT+NET+NAT+NEL

IF (myid.EQ.0) GOTO 88

ALLOCATE (iAux(1:ndof))
ALLOCATE (dCoor(3,1:ndof))
ALLOCATE (iCoor(1:ndof))

iAux = 0
dCoor = 0d0
iCoor = 0

CALL MemoryPrint(1,'s','comm0')

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

CALL MemoryPrint(1,'s','comm1')

DO pID=1,subnodes
 IF (myid.eq.pID) THEN
  DO pJD=1,subnodes
   IF (myid.NE.pJD) THEN
    IF (VerticeCommunicationScheme(pJD).gt.0) THEN
     CALL RECVI_myMPI(pNDOF,pJD)
     CoorST(pJD)%Num = pNDOF
     ALLOCATE(CoorST(pJD)%dCoor(3,pNDOF))
     ALLOCATE(CoorST(pJD)%iCoor(  pNDOF))
     CALL RECVD_myMPI(CoorST(pJD)%dCoor,3*pNDOF,pJD)
     CALL RECVK_myMPI(CoorST(pJD)%iCoor,  pNDOF,pJD)
    END IF
   ELSE
    CoorST(pJD)%Num = jAux
    ALLOCATE(CoorST(pJD)%dCoor(3,jAux))
    ALLOCATE(CoorST(pJD)%iCoor(  jAux))
    CoorST(pJD)%dCoor(:,:) = dCoor(:,1:jAux)
    CoorST(pJD)%iCoor(  :) = iCoor(  1:jAux)
   END IF
  END DO
 ELSE
  IF (VerticeCommunicationScheme(pID).gt.0) THEN
   CALL SENDI_myMPI(jAux ,pID)
   CALL SENDD_myMPI(dCoor,3*jAux,pID)
   CALL SENDK_myMPI(iCoor,  jAux,pID)
  END IF
 END IF

END DO

CALL MemoryPrint(1,'s','comm2')

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

CALL MemoryPrint(1,'s','comm3')

CALL InitOctTree(CoorST(myid)%dCoor,CoorST(myid)%Num)

DO pID=1,subnodes
 jAux = 0
 IF (pID.NE.myid.and.VerticeCommunicationScheme(pID).gt.0) THEN
  DO I=1,CoorST(pID)%Num
!    CALL FindInOctTree(CoorST(myid)%dCoor,CoorST(myid)%Num,CoorST(pID)%dCoor(:,I),J,dist)
!    IF (J.gt.0.and.J.le.CoorST(myid)%Num) then
!     IF (DIST.LT.DEpsPrec) THEN 
!      jAux = jAux + 1
!     END IF
!    END IF
!   CALL FindInPeriodicOctTree(CoorST(myid)%dCoor,CoorST(myid)%Num,CoorST(pID)%dCoor(:,I),J,dist,dPeriodicity)
!    IF (J.gt.0.and.J.le.CoorST(myid)%Num) then
!     IF (DIST.LT.DEpsPrec) THEN 
!      jAux = jAux + 1
!     END IF
!    END IF

   P1X = CoorST(pID)%dCoor(1,I)
   P1Y = CoorST(pID)%dCoor(2,I)
   P1Z = CoorST(pID)%dCoor(3,I)
   DO J=1,CoorST(myid)%Num
    P2X = CoorST(myid)%dCoor(1,J)
    P2Y = CoorST(myid)%dCoor(2,J)
    P2Z = CoorST(myid)%dCoor(3,J)
    IF (((ABS(P1X-P2X).LT.DEpsPrec).OR.(ABS(ABS(P1X-P2X)-dPeriodicity(1)).LT.DEpsPrec)).AND.&
        ((ABS(P1Y-P2Y).LT.DEpsPrec).OR.(ABS(ABS(P1Y-P2Y)-dPeriodicity(2)).LT.DEpsPrec)).AND.& 
        ((ABS(P1Z-P2Z).LT.DEpsPrec).OR.(ABS(ABS(P1Z-P2Z)-dPeriodicity(3)).LT.DEpsPrec))) THEN
     jAux = jAux + 1
    END IF
   END DO
  END DO
 END IF
!  WRITE(*,*) myid,pID,jAux
 MGE013(ILEV)%ST(pID)%Num = jAux
END DO

! WRITE(*,'(A,I0,A,1000(I0," "))')  "pID=",myid," :: ", MGE013(ILEV)%ST(:)%Num
! WRITE(*,'(A,I0,A,1000(I0," "))')  "pID=",myid," :: ", VerticeCommunicationScheme
! pause

CALL MemoryPrint(1,'s','comm4')

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
!    CALL FindInOctTree(CoorST(myid)%dCoor,CoorST(myid)%Num,CoorST(pID)%dCoor(:,I),J,dist)
!    IF (J.gt.0.and.J.le.CoorST(myid)%Num) then
!     IF (DIST.LT.DEpsPrec) THEN 
!      jAux = jAux + 1
!      MGE013(ILEV)%ST(pID)%VertLink(1,jAux) = CoorST(myid)%iCoor(J)
!      MGE013(ILEV)%ST(pID)%VertLink(2,jAux) = CoorST(myid)%iCoor(J)
!     END IF
!    END IF
!    CALL FindInPeriodicOctTree(CoorST(myid)%dCoor,CoorST(myid)%Num,CoorST(pID)%dCoor(:,I),J,dist,dPeriodicity)
!    IF (J.gt.0.and.J.le.CoorST(myid)%Num) then
!     IF (DIST.LT.DEpsPrec) THEN 
!      jAux = jAux + 1
!      MGE013(ILEV)%ST(pID)%VertLink(1,jAux) = CoorST(myid)%iCoor(J)
!      MGE013(ILEV)%ST(pID)%VertLink(2,jAux) = CoorST(myid)%iCoor(J)
!     END IF
!    END IF
   P1X = CoorST(pID)%dCoor(1,I)
   P1Y = CoorST(pID)%dCoor(2,I)
   P1Z = CoorST(pID)%dCoor(3,I)
   DO J=1,CoorST(myid)%Num
    P2X = CoorST(myid)%dCoor(1,J)
    P2Y = CoorST(myid)%dCoor(2,J)
    P2Z = CoorST(myid)%dCoor(3,J)
    IF (((ABS(P1X-P2X).LT.DEpsPrec).OR.(ABS(ABS(P1X-P2X)-dPeriodicity(1)).LT.DEpsPrec)).AND.&
        ((ABS(P1Y-P2Y).LT.DEpsPrec).OR.(ABS(ABS(P1Y-P2Y)-dPeriodicity(2)).LT.DEpsPrec)).AND.& 
        ((ABS(P1Z-P2Z).LT.DEpsPrec).OR.(ABS(ABS(P1Z-P2Z)-dPeriodicity(3)).LT.DEpsPrec))) THEN
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

CALL MemoryPrint(1,'s','comm5')

CALL FreeOctTree()

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

CALL MemoryPrint(1,'w','comm6')

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

CALL MemoryPrint(1,'w','comm-done')


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
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE GetParPresEntries(P,PP)
USE PP3D_MPI
USE def_feat, ONLY: ILEV

IMPLICIT NONE
INTEGER I,J,pID,pJD,nSize
REAL*8 P(*),PP(*)

IF (myid.ne.MASTER) THEN

DO pID=1,subnodes
 IF (myid.eq.pID) THEN
  DO pJD=1,subnodes
   IF (myid.NE.pJD.AND.MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSize = MGE013(ILEV)%SP(pJD)%nElems(2)
     CALL RECVD_myMPI(MGE013(ILEV)%SP(pJD)%RDVect,nSIZE,pJD)
     DO I=1,nSize
      J = MGE013(ILEV)%SP(pJD)%VertLink(2,I)
      PP(J) = MGE013(ILEV)%SP(pJD)%RDVect(I)
     END DO

   END IF
  END DO
 ELSE
  IF (MGE013(ILEV)%SP(pID)%Num.GT.0) THEN
   nSize = MGE013(ILEV)%SP(pID)%nElems(1)
   DO I=1,nSize
    J = MGE013(ILEV)%SP(pID)%VertLink(1,I)
    MGE013(ILEV)%SP(pID)%SDVect(I) = P(J)
   END DO
   CALL SENDD_myMPI(MGE013(ILEV)%SP(pID)%SDVect,nSIZE,pID)
  END IF
 END IF
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END IF

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE GetHYPREParPressureIndices(Ind_PP)
USE var_QuadScalar, ONLY :  iCommSwitch
IMPLICIT NONE
INTEGER Ind_PP(*)

if (iCommSwitch.eq.1) THEN
 WRITE(*,*) 'Subroutine not yet implemented ... '
 STOP
END IF
if (iCommSwitch.eq.2) THEN
 WRITE(*,*) 'Subroutine not yet implemented ... '
 STOP
END IF

if (iCommSwitch.eq.3) CALL GetHYPREParPressureIndicesSUPER(Ind_PP)

if (iCommSwitch.eq.4) CALL GetHYPREParPressureIndicesSUPER(Ind_PP)

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE GetParPressure(P,PP)
USE var_QuadScalar, ONLY :  iCommSwitch
IMPLICIT NONE
REAL*8 P(*),PP(*)

if (iCommSwitch.eq.1) CALL GetParPressureNEW(P,PP)
if (iCommSwitch.eq.2) CALL GetParPressureOLD(P,PP)
if (iCommSwitch.eq.3) CALL GetParPressureSUPER(P,PP)
if (iCommSwitch.eq.4) CALL GetParPressureRec(P,PP) 

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013SUM(FX)
USE var_QuadScalar, ONLY :  iCommSwitch
IMPLICIT NONE
REAL*8 FX(*)

if (iCommSwitch.eq.1) CALL E013SumNEW(FX)
if (iCommSwitch.eq.2) CALL E013SumOLD(FX)
if (iCommSwitch.eq.3) CALL E013SumSUPER(FX)
if (iCommSwitch.eq.4) CALL E013SumRec(FX)

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013Sum3(FX1,FX2,FX3)
USE var_QuadScalar, ONLY :  iCommSwitch
IMPLICIT NONE
REAL*8 FX1(*),FX2(*),FX3(*)

if (iCommSwitch.eq.2) CALL E013Sum3OLD(FX1,FX2,FX3)
if (iCommSwitch.eq.3) CALL E013Sum3SUPER(FX1,FX2,FX3)
if (iCommSwitch.eq.4) CALL E013Sum3Rec(FX1,FX2,FX3)

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWSum(FX)
USE var_QuadScalar, ONLY :  iCommSwitch
IMPLICIT NONE
REAL*8 FX(*)

if (iCommSwitch.eq.1) CALL E013UVWSumNEW(FX)
if (iCommSwitch.eq.2) CALL E013UVWSumOLD(FX)
if (iCommSwitch.eq.3) CALL E013UVWSumSUPER(FX)
if (iCommSwitch.eq.4) CALL E013UVWSumRec(FX)

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWMAT(A11,A22,A33,KLDA,NU) !ok
USE var_QuadScalar, ONLY :  iCommSwitch
IMPLICIT NONE
REAL*8 A11(*),A22(*),A33(*)
INTEGER KLDA(*),NU

if (iCommSwitch.eq.2) CALL E013UVWMAT_OLD(A11,A22,A33,KLDA,NU)
if (iCommSwitch.eq.3) CALL E013UVWMAT_SUPER(A11,A22,A33,KLDA,NU)
if (iCommSwitch.eq.4) CALL E013UVWMAT_Rec(A11,A22,A33,KLDA,NU)

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013MAT(A,KLDA,NU) !ok
USE var_QuadScalar, ONLY :  iCommSwitch
IMPLICIT NONE
REAL*8 A(*)
INTEGER KLDA(*),NU

if (iCommSwitch.eq.2) CALL E013MAT_OLD(A,KLDA,NU)
if (iCommSwitch.eq.3) CALL E013MAT_SUPER(A,KLDA,NU)
if (iCommSwitch.eq.4) CALL E013MAT_Rec(A,KLDA,NU)

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013SumNEW(FX) !ok
USE var_QuadScalar, ONLY : GlobalParallelList1,GlobalParallelList2,GlobalParallelBufferIn,&
    GlobalParallelBufferOut,GlobalNList,GlobalNBuffer,myStat
USE PP3D_MPI
USE def_feat, ONLY: ILEV

REAL*8  FX(*)
INTEGER i,nLenght
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)
 
 nLenght = GlobalNBuffer(ILEV)

 GlobalParallelBufferIn(1:nLenght) = 0d0

 DO i = 1,GlobalNList(ILEV)
  GlobalParallelBufferIn(GlobalParallelList1(ILEV)%x(i)) = FX(GlobalParallelList2(ILEV)%x(i))
 END DO

 CALL MPI_Allreduce(GlobalParallelBufferIn,GlobalParallelBufferOut,nLenght,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SUBS,iERR)

 DO i = 1,GlobalNList(ILEV)
  FX(GlobalParallelList2(ILEV)%x(i)) = GlobalParallelBufferOut(GlobalParallelList1(ILEV)%x(i))
 END DO

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)

END IF

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013PerSum(FX1,FX2,DC) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV

REAL*8  FX1(*),FX2(*),DC(*)
INTEGER I,pID,pJD,nSIZE,nEIGH
INTEGER MEQ,MEQ1,MEQ2,MEQ3


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

     DO I=1,nSIZE
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = FX1(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = FX2(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) =  DC(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
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

     jCount = 0
     DO I=1,nSIZE
       DC_my   = DC(MGE013(ILEV)%ST(pJD)%VertLink(2,I))
       DC_recv = MGE013(ILEV)%ST(pJD)%RDVect(MEQ3+I)
       IF (ABS(DC_my-DC_recv).lt.dEpsPrec) THEN

        FX1(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
        FX1(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ1+I)
        FX2(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
        FX2(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ2+I)
      ELSE
!        jCount =      jCount + 1
      END IF
     END DO
!      IF (jCount.GT.0)  WRITE(*,'(A,3I)') "na ja ...", myid,pjd,jCount

   END IF
 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E013PerSum
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWSumNEW(FX) !ok
USE var_QuadScalar, ONLY : GlobalParallelList1,GlobalParallelList2,GlobalParallelBufferIn,&
    GlobalParallelBufferOut,GlobalNList,GlobalNBuffer,knvt,knet,knat,knel,myStat
USE PP3D_MPI
USE def_feat, ONLY: ILEV

IMPLICIT NONE
REAL*8  FX(*)
INTEGER i,nLenght,ndof,ii,jj
INTEGER MEQ1,MEQ2,MEQ3
INTEGER LEQ1,LEQ2,LEQ3
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

 nLenght = GlobalNBuffer(ILEV)
 ndof = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)

 LEQ1 =0
 LEQ2 =ndof
 LEQ3 =2*ndof
 MEQ1 = 0
 MEQ2 = nLenght
 MEQ3 = 2*nLenght

 GlobalParallelBufferIn(1:3*nLenght) = 0d0

 DO i = 1,GlobalNList(ILEV)
  ii = GlobalParallelList1(ILEV)%x(i)
  jj = GlobalParallelList2(ILEV)%x(i)
  GlobalParallelBufferIn(MEQ1+ii) = FX(LEQ1+jj)
  GlobalParallelBufferIn(MEQ2+ii) = FX(LEQ2+jj)
  GlobalParallelBufferIn(MEQ3+ii) = FX(LEQ3+jj)
 END DO
 
 CALL MPI_Allreduce(GlobalParallelBufferIn,GlobalParallelBufferOut,3*nLenght,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SUBS,iERR)

 DO i = 1,GlobalNList(ILEV)
  ii = GlobalParallelList1(ILEV)%x(i)
  jj = GlobalParallelList2(ILEV)%x(i)
  FX(LEQ1+jj) = GlobalParallelBufferOut(MEQ1+ii) 
  FX(LEQ2+jj) = GlobalParallelBufferOut(MEQ2+ii) 
  FX(LEQ3+jj) = GlobalParallelBufferOut(MEQ3+ii) 
 END DO

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)

END IF

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE GetParPressureNEW(P,PP)
USE PP3D_MPI , ONLY: myid,master,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SUBS
USE def_feat, ONLY: ILEV
USE var_QuadScalar, ONLY :HlobalParallelList1,HlobalParallelList2,&
    HlobalParallelBufferIn,HlobalParallelBufferOut,HlobalNList,HlobalNBuffer,&
    HlobalParallelList3,myStat

IMPLICIT NONE
INTEGER ierr
INTEGER I,pID,pJD,nLenght,iP,jP
REAL*8 P(*),PP(*)
CHARACTER*80 cFile
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)
 
 nLenght = HlobalNBuffer(ILEV)

 HlobalParallelBufferIn  = 0d0

 DO i=1,HlobalNList(ILEV)
  iP = 4*(HlobalParallelList1(ilev)%x(i)-1) + 1
  jP = 4*(HlobalParallelList2(ilev)%x(i)-1) + 1
  HlobalParallelBufferIn(iP:iP+3) = P(jP:jP+3)
 END DO

 CALL MPI_Allreduce(HlobalParallelBufferIn,HlobalParallelBufferOut,4*nLenght,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SUBS,iERR)

 DO i=1,SIZE(HlobalParallelList3(ilev)%x)
  iP = 4*(i-1) + 1
  jP = 4*(HlobalParallelList3(ilev)%x(i)-1) + 1
  PP(iP:iP+3) = HlobalParallelBufferOut(jP:jP+3) 
 END DO

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommP = myStat%tCommP + (tt1-tt0)
 
END IF

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE GetParPressureSUPER(P,PP)
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar, ONLY :CommOrder,myStat,BaSynch

IMPLICIT NONE
! INTEGER ierr
INTEGER I,iRound,pJD,nLenght,iP,jP,nSize,II,JJ
REAL*8 P(*),PP(*)
CHARACTER*80 cFile
REAL*4 tt0,tt1
INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)

IF (myid.ne.MASTER) THEN

 send_req = MPI_REQUEST_NULL
 recv_req = MPI_REQUEST_NULL

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0) THEN
   IF (MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%SP(pJD)%Num

     IF (myid.lt.pJD) THEN
      DO I=1,nSize
       II = 4*(I-1)+1
       JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(1,I)-1)+1
       if(jj.gt.0)then
         MGE013(ILEV)%SP(pJD)%SDVect(II+0) = P(JJ+0)
         MGE013(ILEV)%SP(pJD)%SDVect(II+1) = P(JJ+1)
         MGE013(ILEV)%SP(pJD)%SDVect(II+2) = P(JJ+2)
         MGE013(ILEV)%SP(pJD)%SDVect(II+3) = P(JJ+3)
       end if
      END DO
     END IF
     
     IF (myid.gt.pJD) THEN
      DO I=1,nSize
       II = 4*(I-1)+1
       ! Here is something going on
       ! JJ is 0 sometimes
       JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(1,I)-1)+1
       if(jj .gt. 0)then
         MGE013(ILEV)%SP(pJD)%SDVect(II+0) = P(JJ+0)
         MGE013(ILEV)%SP(pJD)%SDVect(II+1) = P(JJ+1)
         MGE013(ILEV)%SP(pJD)%SDVect(II+2) = P(JJ+2)
         MGE013(ILEV)%SP(pJD)%SDVect(II+3) = P(JJ+3)
       end if
      END DO
     END IF
   END IF
  END IF
 END DO

! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0) THEN
   IF (MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%SP(pJD)%Num

     !!!!     sends pID ----> pJD
     IF (myid.lt.pJD) THEN
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%SP(pJD)%SDVect,4*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%SP(pJD)%SDVect,4*nSIZE,pJD)
      END IF
     ELSE
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,pJD)
      END IF
     END IF
     
     !!!!     sends pJD ----> pID
     IF (myid.lt.pJD) THEN
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,pJD)
      END IF
     ELSE
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%SP(pJD)%SDVect,4*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%SP(pJD)%SDVect,4*nSIZE,pJD)
      END IF
     END IF

   END IF
  END IF
 END DO

 IF (BaSynch) then
  DO pJD=1,subnodes
   IF ((MGE013(ILEV)%SP(pJD)%Num.GT.0)) THEN
    CALL MPI_Wait(recv_req(pJD),STATUS, IERR )
    CALL MPI_Wait(send_req(pJD),STATUS, IERR )
   END IF
  END DO
 ELSE
  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 END IF
 CALL ztime(tt1)
 myStat%tCommP = myStat%tCommP + (tt1-tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0) THEN
   IF (MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%SP(pJD)%Num

     IF (myid.gt.pJD) THEN
      DO I=1,nSize
       II = 4*(I-1)+1
       JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(2,I)-1)+1
       if(jj.gt.0)then
        PP(JJ+0) = MGE013(ILEV)%SP(pJD)%RDVect(II+0)
        PP(JJ+1) = MGE013(ILEV)%SP(pJD)%RDVect(II+1)
        PP(JJ+2) = MGE013(ILEV)%SP(pJD)%RDVect(II+2)
        PP(JJ+3) = MGE013(ILEV)%SP(pJD)%RDVect(II+3)
       end if
      END DO
     END IF
     
     IF (myid.lt.pJD) THEN
      DO I=1,nSize
       II = 4*(I-1)+1
       JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(2,I)-1)+1
       if(jj.gt.0)then
        PP(JJ+0) = MGE013(ILEV)%SP(pJD)%RDVect(II+0)
        PP(JJ+1) = MGE013(ILEV)%SP(pJD)%RDVect(II+1)
        PP(JJ+2) = MGE013(ILEV)%SP(pJD)%RDVect(II+2)
        PP(JJ+3) = MGE013(ILEV)%SP(pJD)%RDVect(II+3)
       end if
      END DO
     END IF
    END IF
   END IF
 END DO

END IF

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE GetHYPREParPressureIndicesSUPER(Ind_PP)
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar, ONLY :CommOrder,myStat,BaSynch,myHYPRE

IMPLICIT NONE
! INTEGER ierr
INTEGER I,iRound,pJD,nLenght,iP,jP,nSize,II,JJ
INTEGER Ind_PP(*)
CHARACTER*80 cFile
REAL*4 tt0,tt1
INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)

IF (myid.ne.MASTER) THEN

 send_req = MPI_REQUEST_NULL
 recv_req = MPI_REQUEST_NULL

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0) THEN
   IF (MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%SP(pJD)%Num

     IF (myid.lt.pJD) THEN
      DO I=1,nSize
       II = 4*(I-1)+1
       JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(1,I)-1)+1
       if(jj.gt.0)then
         MGE013(ILEV)%SP(pJD)%SDVect(II+0) = dble(myHYPRE%Numbering(JJ+0))
         MGE013(ILEV)%SP(pJD)%SDVect(II+1) = dble(myHYPRE%Numbering(JJ+1))
         MGE013(ILEV)%SP(pJD)%SDVect(II+2) = dble(myHYPRE%Numbering(JJ+2))
         MGE013(ILEV)%SP(pJD)%SDVect(II+3) = dble(myHYPRE%Numbering(JJ+3))
       end if
      END DO
     END IF
     
     IF (myid.gt.pJD) THEN
      DO I=1,nSize
       II = 4*(I-1)+1
       ! Here is something going on
       ! JJ is 0 sometimes
       JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(1,I)-1)+1
       if(jj .gt. 0)then
         MGE013(ILEV)%SP(pJD)%SDVect(II+0) = dble(myHYPRE%Numbering(JJ+0))
         MGE013(ILEV)%SP(pJD)%SDVect(II+1) = dble(myHYPRE%Numbering(JJ+1))
         MGE013(ILEV)%SP(pJD)%SDVect(II+2) = dble(myHYPRE%Numbering(JJ+2))
         MGE013(ILEV)%SP(pJD)%SDVect(II+3) = dble(myHYPRE%Numbering(JJ+3))
       end if
      END DO
     END IF
   END IF
  END IF
 END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0) THEN
   IF (MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%SP(pJD)%Num

     !!!!     sends pID ----> pJD
     IF (myid.lt.pJD) THEN
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%SP(pJD)%SDVect,4*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%SP(pJD)%SDVect,4*nSIZE,pJD)
      END IF
     ELSE
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,pJD)
      END IF
     END IF
     
     !!!!     sends pJD ----> pID
     IF (myid.lt.pJD) THEN
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,pJD)
      END IF
     ELSE
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%SP(pJD)%SDVect,4*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%SP(pJD)%SDVect,4*nSIZE,pJD)
      END IF
     END IF

   END IF
  END IF
 END DO

 IF (BaSynch) then
  DO pJD=1,subnodes
   IF ((MGE013(ILEV)%SP(pJD)%Num.GT.0)) THEN
    CALL MPI_Wait(recv_req(pJD),STATUS, IERR )
    CALL MPI_Wait(send_req(pJD),STATUS, IERR )
   END IF
  END DO
 ELSE
  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 END IF
 CALL ztime(tt1)
 myStat%tCommP = myStat%tCommP + (tt1-tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0) THEN
   IF (MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%SP(pJD)%Num

     IF (myid.gt.pJD) THEN
      DO I=1,nSize
       II = 4*(I-1)+1
       JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(2,I)-1)+1
       if(jj.gt.0)then
        Ind_PP(JJ+0) = INT(MGE013(ILEV)%SP(pJD)%RDVect(II+0))
        Ind_PP(JJ+1) = INT(MGE013(ILEV)%SP(pJD)%RDVect(II+1))
        Ind_PP(JJ+2) = INT(MGE013(ILEV)%SP(pJD)%RDVect(II+2))
        Ind_PP(JJ+3) = INT(MGE013(ILEV)%SP(pJD)%RDVect(II+3))
       end if
      END DO
     END IF
     
     IF (myid.lt.pJD) THEN
      DO I=1,nSize
       II = 4*(I-1)+1
       JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(2,I)-1)+1
       if(jj.gt.0)then
        Ind_PP(JJ+0) = INT(MGE013(ILEV)%SP(pJD)%RDVect(II+0))
        Ind_PP(JJ+1) = INT(MGE013(ILEV)%SP(pJD)%RDVect(II+1))
        Ind_PP(JJ+2) = INT(MGE013(ILEV)%SP(pJD)%RDVect(II+2))
        Ind_PP(JJ+3) = INT(MGE013(ILEV)%SP(pJD)%RDVect(II+3))
       end if
      END DO
     END IF
    END IF
   END IF
 END DO

END IF

END SUBROUTINE GetHYPREParPressureIndicesSUPER
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE ExtraxtParallelPattern() !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMAX
USE var_QuadScalar, ONLY : CommOrder

IMPLICIT NONE
INTEGER I,J,pID,pJD,nSIZE,nEIGH
REAL*4 tt0,tt1
INTEGER, ALLOCATABLE :: iTableS(:,:),iTableR(:,:)

ILEV = NLMAX

ALLOCATE(iTableS(subnodes,subnodes))
ALLOCATE(iTableR(subnodes,subnodes))
iTableS = 0

IF (myid.ne.MASTER) THEN
 
 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     iTableS(myid,pJD) = 1

    END IF
   END DO
  END IF
 END DO
 
CALL MPI_Allreduce(iTableS,iTableR,subnodes*subnodes,MPI_INTEGER,MPI_SUM,MPI_COMM_SUBS,iERR)

END IF

CALL OrganizeComm(iTableR,subnodes)
! WRITE(*,*) CommOrder(myid)%x
! pause

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013SumSUPER(FX) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar, ONLY : CommOrder,myStat,BaSynch

REAL*8  FX(*)
INTEGER I,pJD,nSIZE,nEIGH,iRound
REAL*4 tt0,tt1
INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)


!if (myid.eq.1) write(*,*) 'here it goes...'

IF (myid.ne.MASTER) THEN

 send_req = MPI_REQUEST_NULL
 recv_req = MPI_REQUEST_NULL
 
!  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)

  !IF (pJD.ne.0.and.MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
  IF (pJD.ne.0) THEN
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

       nSIZE = MGE013(ILEV)%ST(pJD)%Num

       !!!!     sends pID ----> pJD
       IF (myid.lt.pJD) THEN
        DO I=1,nSIZE
         MGE013(ILEV)%ST(pJD)%SDVect(I) = FX(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
        END DO
        
        IF (BaSynch) then
         CALL MPI_ISEND(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
        ELSE
         CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,pJD)
        END IF
       ELSE
        IF (BaSynch) then
         CALL MPI_IRECV(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
        ELSE
         CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,pJD)
        END IF
       END IF
       
       !!!!     sends pJD ----> pID
       IF (myid.lt.pJD) THEN
        IF (BaSynch) then
         CALL MPI_IRECV(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
        ELSE
         CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,pJD)
        END IF
       ELSE
        DO I=1,nSIZE
         MGE013(ILEV)%ST(pJD)%SDVect(I) = FX(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
        END DO      
        IF (BaSynch) then
         CALL MPI_ISEND(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
        ELSE
         CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,pJD)
        END IF
       END IF

    END IF
  END IF

 END DO

 IF (BaSynch) then
  DO pJD=1,subnodes
   IF ((MGE013(ILEV)%ST(pJD)%Num.GT.0)) THEN
    CALL MPI_Wait(send_req(pJD),STATUS, IERR )
    CALL MPI_Wait(recv_req(pJD),STATUS, IERR )
   END IF
  END DO
 ELSE
  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 END IF
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)
 
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

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013SumOLD(FX) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar, ONLY : myStat

REAL*8  FX(*)
INTEGER I,pID,pJD,nSIZE,nEIGH
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)
 
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

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)

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

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013Max_SUPER(FX) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar, ONLY : CommOrder,myStat

INTEGER  FX(*)
INTEGER I,pJD,nSIZE,nEIGH,iRound
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_CdOMM_SUBS,IERR)
 CALL ztime(tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)

  !IF (pJD.ne.0.and.MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
  IF (pJD.ne.0) THEN
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

       nSIZE = MGE013(ILEV)%ST(pJD)%Num

       !!!!     sends pID ----> pJD
       IF (myid.lt.pJD) THEN
        DO I=1,nSIZE
         MGE013(ILEV)%ST(pJD)%SDVect(I) = DBLE(FX(MGE013(ILEV)%ST(pJD)%VertLink(1,I)))
        END DO
        
        CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,pJD)
       ELSE
        CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,pJD)
       END IF
       
       !!!!     sends pJD ----> pID
       IF (myid.lt.pJD) THEN
        CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,pJD)
       ELSE
        DO I=1,nSIZE
         MGE013(ILEV)%ST(pJD)%SDVect(I) = FX(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
        END DO      
        CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,pJD)
       END IF

    END IF
  END IF

 END DO

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)
 
 DO pJD=1,subnodes
   IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     DO I=1,nSIZE
       FX(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       NINT(MAX(DBLE(FX(MGE013(ILEV)%ST(pJD)%VertLink(2,I))),MGE013(ILEV)%ST(pJD)%RDVect(I)))
     END DO
   END IF
  END DO

END IF

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWSumSUPER(FX) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:CommOrder,knvt,knet,knat,knel,myStat,BaSynch

REAL*8  FX(*)
INTEGER I,pJD,nSIZE,nEIGH
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3
REAL*4 tt0,tt1
INTEGER recv_req(numnodes),send_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)

LEQ = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
LEQ1 =0
LEQ2 =LEQ
LEQ3 =2*LEQ

req = MPI_REQUEST_NULL

IF (myid.ne.MASTER) THEN

!  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0) THEN
    IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     !!!!     sends pID ----> pJD
     IF (myid.lt.pJD) THEN
      DO I=1,nSIZE
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = FX(LEQ1+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = FX(LEQ2+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = FX(LEQ3+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      END DO

      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,pJD)
      END IF
     ELSE
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,pJD)
      END IF
     END IF
     
     !!!!     sends pJD ----> pID
     IF (myid.lt.pJD) THEN
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,pJD)
      END IF
     ELSE
      DO I=1,nSIZE
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = FX(LEQ1+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = FX(LEQ2+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = FX(LEQ3+MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      END DO
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,pJD)
      END IF
     END IF

    END IF
  ELSE
!   req(pJD) = 5
  END IF
 END DO

 IF (BaSynch) then
  DO pJD=1,subnodes
   IF ((MGE013(ILEV)%ST(pJD)%Num.GT.0)) THEN
    CALL MPI_Wait(recv_req(pJD),STATUS, IERR )
    CALL MPI_Wait(send_req(pJD),STATUS, IERR )
   END IF
  END DO
 ELSE
  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 END IF
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)

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

! WRITE(*,*) "Total time for COMM :",   myStat%tCommV
! pause

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013Sum3SUPER(FX1,FX2,FX3) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:CommOrder,knvt,knet,knat,knel,myStat,BaSynch

REAL*8  FX1(*),FX2(*),FX3(*)
INTEGER I,pJD,nSIZE,nEIGH
INTEGER MEQ,MEQ1,MEQ2,MEQ3
REAL*4 tt0,tt1
INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)

IF (myid.ne.MASTER) THEN

 send_req = MPI_REQUEST_NULL
 recv_req = MPI_REQUEST_NULL
 
!  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0) THEN
   IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     !!!!     sends pID ----> pJD
     IF (myid.lt.pJD) THEN
      DO I=1,nSIZE
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = FX1(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = FX2(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = FX3(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      END DO
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,pJD)
      END IF
     ELSE
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,pJD)
      END IF
     END IF
     
     !!!!     sends pJD ----> pID
     IF (myid.lt.pJD) THEN
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,pJD)
      END IF
     ELSE
      DO I=1,nSIZE
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = FX1(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = FX2(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = FX3(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      END DO
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,pJD)
      END IF
     END IF

    END IF
   END IF
 END DO

 IF (BaSynch) then
  DO pJD=1,subnodes
   IF ((MGE013(ILEV)%ST(pJD)%Num.GT.0)) THEN
    CALL MPI_Wait(send_req(pJD),STATUS, IERR )
    CALL MPI_Wait(recv_req(pJD),STATUS, IERR )
   END IF
  END DO
 ELSE
  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 END IF
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)

 DO pJD=1,subnodes
   IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     DO I=1,nSIZE
       FX1(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX1(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ1+I)
       FX2(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX2(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ2+I)
       FX3(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX3(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ3+I)
     END DO

   END IF
 END DO
 
END IF

END 
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013Sum3OLD(FX1,FX2,FX3) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat

REAL*8  FX1(*),FX2(*),FX3(*)
INTEGER I,pID,pJD,nSIZE,nEIGH
INTEGER MEQ,MEQ1,MEQ2,MEQ3
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

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
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = FX1(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = FX2(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = FX3(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
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

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)
 
 DO pJD=1,subnodes
   IF (MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     DO I=1,nSIZE
       FX1(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX1(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ1+I)
       FX2(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX2(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ2+I)
       FX3(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) = &
       FX3(MGE013(ILEV)%ST(pJD)%VertLink(2,I)) +  MGE013(ILEV)%ST(pJD)%RDVect(MEQ3+I)
     END DO

   END IF
 END DO

END IF

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWSumOLD(FX) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat

REAL*8  FX(*)
INTEGER I,pID,pJD,nSIZE,nEIGH
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3
REAL*4 tt0,tt1


LEQ = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
LEQ1 =0
LEQ2 =LEQ
LEQ3 =2*LEQ

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

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

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)
 
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

END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWMAT_OLD(A11,A22,A33,KLDA,NU) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:myStat

REAL*8 A11(*),A22(*),A33(*)
INTEGER KLDA(*),NU
INTEGER I,J
INTEGER pID,pJD,nSIZE
INTEGER MEQ,MEQ1,MEQ2,MEQ3
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 DO I=1,NU
  MGE013(ILEV)%UE11(I)=A11(KLDA(I))
  MGE013(ILEV)%UE22(I)=A22(KLDA(I))
  MGE013(ILEV)%UE33(I)=A33(KLDA(I))
 ENDDO

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

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

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)

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

END SUBROUTINE E013UVWMAT_OLD
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWMAT_SUPER(A11,A22,A33,KLDA,NU) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:CommOrder,myStat,BaSynch

REAL*8 A11(*),A22(*),A33(*)
INTEGER KLDA(*),NU
INTEGER I,J
INTEGER pID,pJD,nSIZE
INTEGER MEQ,MEQ1,MEQ2,MEQ3
REAL*4 tt0,tt1
INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)

IF (myid.ne.MASTER) THEN

 send_req = MPI_REQUEST_NULL
 recv_req = MPI_REQUEST_NULL

 DO I=1,NU
  MGE013(ILEV)%UE11(I)=A11(KLDA(I))
  MGE013(ILEV)%UE22(I)=A22(KLDA(I))
  MGE013(ILEV)%UE33(I)=A33(KLDA(I))
 ENDDO

!  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0)then
   if(MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN

     nSIZE = MGE013(ILEV)%ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     !!!!     sends pID ----> pJD
     IF (myid.lt.pJD) THEN
      DO I=1,nSIZE
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = MGE013(ILEV)%UE11(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = MGE013(ILEV)%UE22(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = MGE013(ILEV)%UE33(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      END DO
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,pJD)
      END IF
     ELSE
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,pJD)
      END IF
     END IF
     
     !!!!     sends pJD ----> pID
     IF (myid.lt.pJD) THEN
      IF (BaSynch) then
       CALL MPI_IRECV(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,recv_req(pJD),IERR)
      ELSE
       CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,3*nSIZE,pJD)
      END IF
     ELSE
      DO I=1,nSIZE
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ1+I) = MGE013(ILEV)%UE11(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ2+I) = MGE013(ILEV)%UE22(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
       MGE013(ILEV)%ST(pJD)%SDVect(MEQ3+I) = MGE013(ILEV)%UE33(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      END DO
      IF (BaSynch) then
       CALL MPI_ISEND(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,MPI_DOUBLE_PRECISION,pJD,101,MPI_COMM_WORLD,send_req(pJD),IERR)
      ELSE
       CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,3*nSIZE,pJD)
      END IF
     END IF

   END IF
  END IF

 END DO

 IF (BaSynch) then
  DO pJD=1,subnodes
   IF ((MGE013(ILEV)%ST(pJD)%Num.GT.0)) THEN
    CALL MPI_Wait(send_req(pJD),STATUS, IERR )
    CALL MPI_Wait(recv_req(pJD),STATUS, IERR )
   END IF
  END DO
 ELSE
  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 END IF
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)

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

END SUBROUTINE E013UVWMAT_SUPER
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013MAT_OLD(A,KLDA,NU) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:myStat

REAL*8 A(*)
INTEGER KLDA(*),NU
INTEGER I,J
INTEGER pID,pJD,nSIZE
REAL*4 tt0,tt1

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

END SUBROUTINE E013MAT_OLD
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013MAT_SUPER(A,KLDA,NU) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:CommOrder,myStat

REAL*8 A(*)
INTEGER KLDA(*),NU
INTEGER I,J
INTEGER pID,pJD,nSIZE
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

 DO I=1,NU
  MGE013(ILEV)%UE(I)=A(KLDA(I))
 ENDDO

 DO iRound=1,SIZE(CommOrder(myid)%x)
  pJD = CommOrder(myid)%x(iRound)
  IF (pJD.ne.0.and.MGE013(ILEV)%ST(pJD)%Num.GT.0) THEN
     nSIZE = MGE013(ILEV)%ST(pJD)%Num

     !!!!     sends pID ----> pJD
     IF (myid.lt.pJD) THEN
      DO I=1,nSIZE
       MGE013(ILEV)%ST(pJD)%SDVect(I) = MGE013(ILEV)%UE(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      END DO
      
      CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,pJD)
     ELSE
      CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,pJD)
     END IF
     
     !!!!     sends pJD ----> pID
     IF (myid.lt.pJD) THEN
      CALL RECVD_myMPI(MGE013(ILEV)%ST(pJD)%RDVect,nSIZE,pJD)
     ELSE
      DO I=1,nSIZE
       MGE013(ILEV)%ST(pJD)%SDVect(I) = MGE013(ILEV)%UE(MGE013(ILEV)%ST(pJD)%VertLink(1,I))
      END DO      
      CALL SENDD_myMPI(MGE013(ILEV)%ST(pJD)%SDVect,nSIZE,pJD)
     END IF

    END IF
 END DO

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommV = myStat%tCommV + (tt1-tt0)

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

END SUBROUTINE E013MAT_SUPER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E013GATHR_L1(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY :GlobalNumberingQ2,GlobalNumberingP1,myGlobalNumberingMap
IMPLICIT NONE
REAL*8 D(*)
INTEGER N,ndof
INTEGER pID,pN,I,J

IF (myid.eq.0) THEN
 DO pID=1,subnodes
  DO i=1,myGlobalNumberingMap(pid)%ndof_Q2
   j = myGlobalNumberingMap(pid)%indQ2(i)
   myGlobalNumberingMap(pid)%dBufferQ2(i) = D(j)
  END DO
 END DO
END IF

IF (myid.eq.0) THEN
 DO pID=1,subnodes
  ndof = myGlobalNumberingMap(pid)%ndof_Q2
  CALL sendD_myMPI(myGlobalNumberingMap(pid)%dBufferQ2,ndof,pID)
 END DO
ELSE
  ndof = 3*N
  CALL recvD_myMPI(D,ndof,0)
END IF

END SUBROUTINE E013GATHR_L1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E013DISTR_L1(D,N)
USE var_QuadScalar, ONLY :GlobalNumberingQ2,GlobalNumberingP1,myGlobalNumberingMap
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N,ndof
INTEGER pID,pN,I,J

IF (myid.eq.0) THEN
 DO pID=1,subnodes
  ndof = myGlobalNumberingMap(pid)%ndof_Q2
  CALL recvD_myMPI(myGlobalNumberingMap(pid)%dBufferQ2,ndof,pID)
 END DO
ELSE
  ndof = 3*N
  CALL sendD_myMPI(D,ndof,0)
END IF


IF (myid.eq.0) THEN
 DO pID=1,subnodes
  DO i=1,myGlobalNumberingMap(pid)%ndof_Q2
   j = myGlobalNumberingMap(pid)%indQ2(i)
   D(j) = D(j) + myGlobalNumberingMap(pid)%dBufferQ2(i)
  END DO
 END DO
END IF

END SUBROUTINE E013DISTR_L1
!
!
!
SUBROUTINE E013_PUTMAT(A11,A22,A33,KLDA,NU) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV

REAL*8 A11(*),A22(*),A33(*)
REAL*8 D11,D22,D33
INTEGER KLDA(*),NU
INTEGER I,J

IF (myid.ne.MASTER) THEN

 DO I=1,NU
  D11 = A11(KLDA(I))
  D22 = A22(KLDA(I))
  D33 = A33(KLDA(I))
  A11(KLDA(I)) = MGE013(ILEV)%UE11(I)
  A22(KLDA(I)) = MGE013(ILEV)%UE22(I)
  A33(KLDA(I)) = MGE013(ILEV)%UE33(I)
  MGE013(ILEV)%UE11(I)=D11
  MGE013(ILEV)%UE22(I)=D22
  MGE013(ILEV)%UE33(I)=D33
 ENDDO

END IF

END SUBROUTINE E013_PUTMAT
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013_GETMAT(A11,A22,A33,KLDA,NU) !ok
USE PP3D_MPI
USE def_feat, ONLY: ILEV

REAL*8 A11(*),A22(*),A33(*)
REAL*8 D11,D22,D33
INTEGER KLDA(*),NU
INTEGER I,J

IF (myid.ne.MASTER) THEN

 DO I=1,NU
  D11 = MGE013(ILEV)%UE11(I)
  D22 = MGE013(ILEV)%UE22(I)
  D33 = MGE013(ILEV)%UE33(I)
  MGE013(ILEV)%UE11(I) = A11(KLDA(I))
  MGE013(ILEV)%UE22(I) = A22(KLDA(I))
  MGE013(ILEV)%UE33(I) = A33(KLDA(I))
  A11(KLDA(I))=D11
  A22(KLDA(I))=D22
  A33(KLDA(I))=D33
 ENDDO

END IF

END SUBROUTINE E013_GETMAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!
!----------------------------------------------
!
SUBROUTINE E013Min(FX) !ok
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
       MIN(FX(MGE013(ILEV)%ST(pJD)%VertLink(2,I)),MGE013(ILEV)%ST(pJD)%RDVect(I))
     END DO

   END IF
 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E013Min

SUBROUTINE COMM_cc_def(d,n)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY :my_crs_e013_map

implicit none
REAL*8 d(*)
INTEGER pid,i,j,n

IF (myid.eq.0) THEN

 d(1:n) = 0d0

 DO pID=1,subnodes
  CALL recvD_myMPI(my_crs_e013_map(pid)%dBuffer,my_crs_e013_map(pid)%cc_ndof,pID)
  DO i=1,my_crs_e013_map(pid)%cc_ndof
   j = my_crs_e013_map(pid)%indE(i)
   d(j) = d(j) + my_crs_e013_map(pid)%dBuffer(i)
  END DO
 END DO

ELSE

  CALL sendD_myMPI(d,n,0)

END IF

END SUBROUTINE COMM_cc_def



SUBROUTINE COMM_cc_sol(d,n)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY :my_crs_e013_map

implicit none
REAL*8 d(*)
INTEGER pid,i,j,n

IF (myid.eq.0) THEN

 DO pID=1,subnodes
  DO i=1,my_crs_e013_map(pid)%cc_ndof
   j = my_crs_e013_map(pid)%indE(i)
   my_crs_e013_map(pid)%dBuffer(i) = d(j)
  END DO
  CALL sendD_myMPI(my_crs_e013_map(pid)%dBuffer,my_crs_e013_map(pid)%cc_ndof,pID)
 END DO

ELSE

  CALL recvD_myMPI(d,n,0)

END IF

END SUBROUTINE COMM_cc_sol
!
!
!
SUBROUTINE OrganizeComm(T,nT)
USE PP3D_MPI, ONLY : myid,subnodes
USE var_QuadScalar, ONLY : CommOrder,nSubCoarseMesh
IMPLICIT NONE
INTEGER nT
INTEGER T(nT,nT)
TYPE tPairs
 INTEGER d(2)
 LOGICAL p,o
END TYPE tPairs
TYPE(tPairs), ALLOCATABLE :: Pairs(:)
TYPE(tPairs), ALLOCATABLE :: Comm(:,:)
INTEGER nPairs,MaxPairs,UsedMaxPairs,MaxComm,LenComm,i,j,k

ALLOCATE(CommOrder(nT))

nPairs = 0
DO i=1,nT
 DO j=1,nT
  IF (T(i,j).ne.0.and.i.lt.j) then
   nPairs = nPairs + 1
  END iF
 END DO
END DO

ALLOCATE(Pairs(nPairs))

nPairs = 0
MaxPairs = 0
DO i=1,nT
 k=0
 DO j=1,nT
  IF (T(i,j).ne.0) k = k + 1 
  IF (T(i,j).ne.0.and.i.lt.j) then
   nPairs = nPairs + 1
   Pairs(nPairs)%d = [i,j]
  END iF
 END DO
 IF (MaxPairs.lt.k) THEN
  MaxPairs = k
  MaxComm  = i
 END IF
END DO

j=mod(nT,2)
LenComm=INT(0.5*(nT-j)) 

UsedMaxPairs = MaxPairs*2
ALLOCATE(Comm(UsedMaxPairs,LenComm))

! write(*,*) "MaxPairs","nPairs","LComm","MComm","sComm"
! write(*,*) MaxPairs,nPairs,LenComm,MaxComm,SIZE(Comm)

CALL PackIt()

CALL OutputComm()

! pause
Call ExtractSingleComm()

 CONTAINS

SUBROUTINE ExtractSingleComm
INTEGER i1,i2

 DO i=1,nT
  ALLOCATE(CommOrder(i)%x(UsedMaxPairs))
  CommOrder(i)%x(:) = 0
 END DO
 
 DO i=1,UsedMaxPairs
  DO j=1,LenComm
   if (comm(i,j)%p) then
     i1 = Comm(i,j)%d(1)
     i2 = Comm(i,j)%d(2)
     CommOrder(i1)%x(i) = i2
     CommOrder(i2)%x(i) = i1
   END IF
  END DO
 END DO

if (myid.eq.1) then
!  WRITE(*,*) 
!  DO i=1,nT
!   WRITE(*,'(I,A$)') i," : "
!   DO j=1,UsedMaxPairs
!    WRITE(*,'(I3$)') CommOrder(i)%x(j)
!   END DO
!   WRITE(*,*) " |"
!  END DO
END if

END 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OutputComm
INTEGER i1,i2
CHARACTER cX*(1)

if (myid.eq.1) then
 WRITE(*,*) 'LenComm: ',LenComm,'UsedMaxPairs: ',UsedMaxPairs
 DO i=1,LenComm
  DO j=1,UsedMaxPairs
   IF (Comm(j,i)%d(1).ne.0) THEN
    if (Comm(j,i)%o) then
     WRITE(*,'(A1,I4,A1,I4,A1$)') '{',Comm(j,i)%d(1),',',Comm(j,i)%d(2),'}'
    else
     WRITE(*,'(A1,I4,A1,I4,A1$)') '[',Comm(j,i)%d(1),',',Comm(j,i)%d(2),']'
    end if
   ELSE
    WRITE(*,'(A$)') '         '
   END IF
  END DO
  write(*,*) ' |'
 END DO
end if

return
if (myid.eq.1) then
 DO i=1,UsedMaxPairs
  DO j=1,LenComm
   if (Comm(i,j)%o) then
    cX=' '
   else
    cX='*'
   end if
   WRITE(*,'(A2,I3,A1,I3,A2,A$)') ' {',Comm(i,j)%d(1),',',Comm(i,j)%d(2),' }',cX
  END DO
  write(*,*) ' '
 END DO
end if

! pause
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Packit
LOGICAL, ALLOCATABLE :: f(:)
INTEGER i1,i2,kSubPart

 allocate(f(nT))

 DO i=1,UsedMaxPairs
  DO j=1,LenComm
   Comm(i,j)%d = [0,0]
   Comm(i,j)%p = .FALSE.
   Comm(i,j)%o = .TRUE.
  END DO
 END DO
 
 Pairs(:)%p = .FALSE.
 DO i=1,SIZE(Pairs)
    kSubPart = FLOOR(DBLE(nT)/DBLE(nSubCoarseMesh)-1d-10)+1
    i1 = FLOOR(DBLE(subnodes-Pairs(i)%d(1)+1)/DBLE(kSubPart)-1d-10)+1
    i2 = FLOOR(DBLE(subnodes-Pairs(i)%d(2)+1)/DBLE(kSubPart)-1d-10)+1
!     i1 = nSubCoarseMesh - i1 + 1
!     i2 = nSubCoarseMesh - i2 + 1
!   i1 = INT(nSubCoarseMesh*Pairs(i)%d(1)/nT)
!   i2 = INT(nSubCoarseMesh*Pairs(i)%d(2)/nT)
  if (i1.eq.i2) then
   Pairs(i)%o = .TRUE.
!   IF (myid.eq.1) WRITE(*,*) i1,' - ', Pairs(i)%d,nSubCoarseMesh
  else
   Pairs(i)%o = .FALSE.
  end if
 END DO
  
 DO i=1,UsedMaxPairs
  f = .false.
  DO j=1,LenComm
   DO k=1,nPairs
    IF ((.not.Pairs(k)%p).and.(.not.Pairs(k)%o)) THEN
     i1 = Pairs(k)%d(1)
     i2 = Pairs(k)%d(2)
     if ((.not.f(i1)).and.(.not.f(i2))) then
      Comm(i,j)%d = Pairs(k)%d
      Comm(i,j)%o = Pairs(k)%o
      Comm(i,j)%p = .true.
      Pairs(k)%p  = .true.
      f(i1)       = .true.
      f(i2)       = .true.
      exit
     end if
    END IF
   END DO
  END DO
 END DO
 
 DO i=1,UsedMaxPairs
  f = .false.
  DO j=1,LenComm
     i1 = Comm(i,j)%d(1)
     i2 = Comm(i,j)%d(2)
     if (i1.ne.0.and.i2.ne.0) then
     f(i1) = .true.
     f(i2) = .true.
     end if
  END DO
  
  DO j=1,LenComm
   IF (.not.Comm(i,j)%p) then
   DO k=1,nPairs
    IF (.not.Pairs(k)%p) THEN
     i1 = Pairs(k)%d(1)
     i2 = Pairs(k)%d(2)
     if ((.not.f(i1)).and.(.not.f(i2))) then
      Comm(i,j)%d = Pairs(k)%d
      Comm(i,j)%o = Pairs(k)%o
      Comm(i,j)%p = .true.
      Pairs(k)%p  = .true.
      f(i1)       = .true.
      f(i2)       = .true.
      exit
     end if
    END IF
   END DO
   end if
  END DO
 END DO

 DO i=1,UsedMaxPairs
  k = 0
  DO j=1,LenComm
   IF (Comm(i,j)%p) k = k + 1
  END DO
  if (k.eq.0) GOTO 1
 END DO
 1 continue
 
 UsedMaxPairs = i-1
 
!  write(*,*) Pairs(:)%p
!  write(*,*) "UsedMaxPairs:",UsedMaxPairs
END 

END 

subroutine GetMyHYPRENumberingLimits(iLow,iUp,N)
USE PP3D_MPI
implicit none
integer iLow,iUp,N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer piVAL,iVAL
INTEGER pID

IF (myid.eq.MASTER) THEN
  iVAL = 0
  DO pID=1,subnodes
   CALL RECVI_myMPI(piVAL,pID)
   iVAL = iVAL + piVAL
   CALL SENDI_myMPI(iVAL,pID)
  END DO
ELSE
  piVAL=N
  CALL SENDI_myMPI(piVAL,0)
  CALL RECVI_myMPI(iVAL,0)
  iUp  = 4*iVAL
  iLow = iUp - 4*N + 1
END IF

end subroutine GetMyHYPRENumberingLimits


SUBROUTINE E010GATHR_L1(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J

 IF (myid.NE.0) THEN
  pN = N
  CALL SENDI_myMPI(pN,0)
  CALL SENDD_myMPI(D,N,0)
 ELSE
!  WRITE(*,*) 'a',allocated(MGE013)
!  IF (allocated(MGE013)) WRITE(*,*) 'a',allocated(MGE013(ILEV)%CRSVect)
!  pause

 DO pID=1,subnodes
   CALL RECVI_myMPI(pN ,pID)
   CALL RECVD_myMPI(MGE013(ILEV)%CRSVect,pN,pID)

   DO I=1,pN
    J = coarse%pELEMLINK(pID,I)
    D(J) = MGE013(ILEV)%CRSVect(I)
   END DO

  END DO

 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E010GATHR_L1

SUBROUTINE E010GATHR_L2(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J,JJ,II,iu,ppN


 IF (myid.NE.0) THEN
  pN = N
  CALL SENDI_myMPI(pN,0)
  CALL SENDD_myMPI(D,N,0)
 ELSE

  ppn = N

  DO pID=1,subnodes
   CALL RECVI_myMPI(pN ,pID)
   CALL RECVD_myMPI(MGE013(ILEV)%CRSVect,pN,pID)

   DO II=1,pN/8
    JJ = coarse%pELEMLINK(pID,II)

    I = II
    J = JJ
    D(J) = MGE013(ILEV)%CRSVect(I)

    DO iu = 1,7

     J = ppN/8 + (JJ-1)*7 + iu
     I =  pN/8 + (II-1)*7 + iu

     D(J) = MGE013(ILEV)%CRSVect(I)
    END DO

   END DO

  END DO

 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E010GATHR_L2

SUBROUTINE E010GATHR_L3(D,N)
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
  CALL SENDD_myMPI(D,N,0)
 ELSE

  ppn = N

  DO pID=1,subnodes
   CALL RECVI_myMPI(pN ,pID)
   CALL RECVD_myMPI(MGE013(ILEV)%CRSVect,pN,pID)

   DO II=1,pN/64
    JJ = coarse%pELEMLINK(pID,II)

    I = II
    J = JJ
    D(J) = MGE013(ILEV)%CRSVect(I)

    DO iw = 1,7
     J = ppN/8 + (JJ-1)*7 + iw
     I =  pN/8 + (II-1)*7 + iw

     D(J) = MGE013(ILEV)%CRSVect(I)
    END DO

    DO iu = 1,7

     JJJ = ppN/64 + (JJ-1)*7 + iu
     III =  pN/64 + (II-1)*7 + iu
     I = III
     J = JJJ

     D(J) = MGE013(ILEV)%CRSVect(I)

     DO iw = 1,7

      J = ppN/8 + (JJJ-1)*7 + iw
      I =  pN/8 + (III-1)*7 + iw
 
      D(J) = MGE013(ILEV)%CRSVect(I)
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

END SUBROUTINE E010GATHR_L3
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE GetParPressureOLD(P,PP)
USE PP3D_MPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar, ONLY : myStat

IMPLICIT NONE
! INTEGER ierr
INTEGER I,pID,pJD,nLenght,iP,jP,nSize,II,JJ
REAL*8 P(*),PP(*)
CHARACTER*80 cFile
REAL*4 tt0,tt1

IF (myid.ne.MASTER) THEN

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt0)

DO pID=1,subnodes
 IF (myid.eq.pID) THEN
  DO pJD=1,subnodes
   IF (myid.NE.pJD.AND.MGE013(ILEV)%SP(pJD)%Num.GT.0) THEN
     nSize = MGE013(ILEV)%SP(pJD)%nElems(2)
     CALL RECVD_myMPI(MGE013(ILEV)%SP(pJD)%RDVect,4*nSIZE,pJD)
     DO I=1,nSize
      II = 4*(I-1)+1
      JJ = 4*(MGE013(ILEV)%SP(pJD)%VertLink(2,I)-1)+1
      ! Here is something going on
      ! JJ is 0 sometimes
      if(jj.gt.0)then
        PP(JJ+0) = MGE013(ILEV)%SP(pJD)%RDVect(II+0)
        PP(JJ+1) = MGE013(ILEV)%SP(pJD)%RDVect(II+1)
        PP(JJ+2) = MGE013(ILEV)%SP(pJD)%RDVect(II+2)
        PP(JJ+3) = MGE013(ILEV)%SP(pJD)%RDVect(II+3)
      end if
     END DO
   END IF
  END DO
 ELSE
  IF (MGE013(ILEV)%SP(pID)%Num.GT.0) THEN
   nSize = MGE013(ILEV)%SP(pID)%nElems(1)
   DO I=1,nSize
    II = 4*(I-1)+1
    JJ = 4*(MGE013(ILEV)%SP(pID)%VertLink(1,I)-1)+1
    if(jj.gt.0)then
      MGE013(ILEV)%SP(pID)%SDVect(II+0) = P(JJ+0)
      MGE013(ILEV)%SP(pID)%SDVect(II+1) = P(JJ+1)
      MGE013(ILEV)%SP(pID)%SDVect(II+2) = P(JJ+2)
      MGE013(ILEV)%SP(pID)%SDVect(II+3) = P(JJ+3)
    end if
   END DO
   CALL SENDD_myMPI(MGE013(ILEV)%SP(pID)%SDVect,4*nSIZE,pID)
  END IF
 END IF
END DO

 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
 CALL ztime(tt1)
 myStat%tCommP = myStat%tCommP + (tt1-tt0)

END IF

END 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! So far only the DISTR/GATHER L1 routines are asynchronous communication !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE E012DISTR_L1(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY : UNF_P_CrsGrid

IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J
INTEGER send_req(numnodes),recv_req(numnodes),my_recv_req,my_send_req
INTEGER STATUS(MPI_STATUS_SIZE)

send_req = MPI_REQUEST_NULL
recv_req = MPI_REQUEST_NULL
my_send_req = MPI_REQUEST_NULL
my_recv_req = MPI_REQUEST_NULL

! not yet initialized
 IF (.not.allocated(UNF_P_CrsGrid)) THEN
  IF (myid.NE.0) THEN
   allocate(UNF_P_CrsGrid(1))
   UNF_P_CrsGrid(1) = N
   CALL SENDI_myMPI(pN,0)
  ELSE
   allocate(UNF_P_CrsGrid(subnodes))
   DO pID=1,subnodes
    CALL RECVI_myMPI(UNF_P_CrsGrid(pID) ,pID)
   END DO
  END IF
 END IF

!  WRITE(*,*) myid, "is here", ILEV
 IF (myid.NE.0) THEN
  CALL MPI_ISEND(D,4*N,MPI_DOUBLE_PRECISION,0,1101,MPI_COMM_WORLD,my_send_req,IERR)
!  CALL SENDD_myMPI(D,4*N,0)
 ELSE
  DO pID=1,subnodes
   pN = UNF_P_CrsGrid(pID)
   CALL MPI_IRECV(MGE013(ILEV)%CRSVect,4*pN,MPI_DOUBLE_PRECISION,pID,1101,MPI_COMM_WORLD,recv_req(pID),IERR)
!   CALL RECVD_myMPI(MGE013(ILEV)%CRSVect,4*pN,pID)

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

 IF (myid.NE.0) THEN
      CALL MPI_Wait(my_send_req,STATUS, IERR )
 ELSE
  DO pID=1,subnodes
      CALL MPI_Wait(recv_req(pID),STATUS, IERR )
  END DO
 END IF
 
! CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE E012DISTR_L1

SUBROUTINE E012GATHR_L1(D,N)
USE PP3D_MPI
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, ONLY : UNF_P_CrsGrid
IMPLICIT NONE
REAL*8 D(*)
INTEGER N
INTEGER pID,pN,I,J 

INTEGER send_req(numnodes),recv_req(numnodes),my_recv_req,my_send_req
INTEGER STATUS(MPI_STATUS_SIZE)

send_req = MPI_REQUEST_NULL
recv_req = MPI_REQUEST_NULL
my_send_req = MPI_REQUEST_NULL
my_recv_req = MPI_REQUEST_NULL

! not yet initialized
 IF (.not.allocated(UNF_P_CrsGrid)) THEN
  IF (myid.NE.0) THEN
   allocate(UNF_P_CrsGrid(1))
   UNF_P_CrsGrid(1) = N
   CALL SENDI_myMPI(pN,0)
  ELSE
   allocate(UNF_P_CrsGrid(subnodes))
   DO pID=1,subnodes
    CALL RECVI_myMPI(UNF_P_CrsGrid(pID) ,pID)
   END DO
  END IF
 END IF
 
!  WRITE(*,*) myid, "is here", ILEV
 IF (myid.NE.0) THEN
!  CALL RECVD_myMPI(D,4*N,0)
  CALL MPI_IRECV(D,4*N,MPI_DOUBLE_PRECISION,0,1102,MPI_COMM_WORLD,my_recv_req,IERR)
 ELSE
  DO pID=1,subnodes

   pN = UNF_P_CrsGrid(pID)
   DO I=1,pN
    J = coarse%pELEMLINK(pID,I)
    MGE013(ILEV)%CRSVect(4*(I-1)+1) = D(4*(J-1)+1)
    MGE013(ILEV)%CRSVect(4*(I-1)+2) = D(4*(J-1)+2)
    MGE013(ILEV)%CRSVect(4*(I-1)+3) = D(4*(J-1)+3)
    MGE013(ILEV)%CRSVect(4*(I-1)+4) = D(4*(J-1)+4)
   END DO
   
   CALL MPI_ISEND(MGE013(ILEV)%CRSVect,4*pN,MPI_DOUBLE_PRECISION,pID,1102,MPI_COMM_WORLD,send_req(pID),IERR)
!   CALL SENDD_myMPI(MGE013(ILEV)%CRSVect,4*pN,pID)
  END DO

!   OPEN (885,FILE ='gfgfgf.txt')
!   DO I=1,N
!   WRITE(885,*) I,D(I)
!   END DO
!   CLOSE(885)
 END IF

 IF (myid.NE.0) THEN
      CALL MPI_Wait(my_recv_req,STATUS, IERR )
 ELSE
  DO pID=1,subnodes
      CALL MPI_Wait(send_req(pID),STATUS, IERR )
  END DO
 END IF
! CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

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
