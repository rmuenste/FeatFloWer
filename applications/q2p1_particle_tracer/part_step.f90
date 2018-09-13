module particle_step

contains

!========================================================================================
!                            Sub: GetPointFromElement
!========================================================================================
SUBROUTINE GetPointFromElement(P8,P,iStart,BF,PR,iP)
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
PARAMETER (Q8=0.125D0)
LOGICAL BF
logical :: ok_flag
REAL*8  P8(3,8), P(3),PR(3)
REAL*8 DJACI(3,3),DJAC(3,3)

DJ11=( P8(1,1)+P8(1,2)+P8(1,3)+P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ12=( P8(2,1)+P8(2,2)+P8(2,3)+P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ13=( P8(3,1)+P8(3,2)+P8(3,3)+P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ21=(-P8(1,1)+P8(1,2)+P8(1,3)-P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ22=(-P8(2,1)+P8(2,2)+P8(2,3)-P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ23=(-P8(3,1)+P8(3,2)+P8(3,3)-P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ31=(-P8(1,1)-P8(1,2)+P8(1,3)+P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ32=(-P8(2,1)-P8(2,2)+P8(2,3)+P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ33=(-P8(3,1)-P8(3,2)+P8(3,3)+P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ41=(-P8(1,1)-P8(1,2)-P8(1,3)-P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ42=(-P8(2,1)-P8(2,2)-P8(2,3)-P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ43=(-P8(3,1)-P8(3,2)-P8(3,3)-P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ51=( P8(1,1)-P8(1,2)+P8(1,3)-P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ52=( P8(2,1)-P8(2,2)+P8(2,3)-P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ53=( P8(3,1)-P8(3,2)+P8(3,3)-P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ61=( P8(1,1)-P8(1,2)-P8(1,3)+P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ62=( P8(2,1)-P8(2,2)-P8(2,3)+P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ63=( P8(3,1)-P8(3,2)-P8(3,3)+P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ71=( P8(1,1)+P8(1,2)-P8(1,3)-P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ72=( P8(2,1)+P8(2,2)-P8(2,3)-P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ73=( P8(3,1)+P8(3,2)-P8(3,3)-P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ81=(-P8(1,1)+P8(1,2)-P8(1,3)+P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ82=(-P8(2,1)+P8(2,2)-P8(2,3)+P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ83=(-P8(3,1)+P8(3,2)-P8(3,3)+P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8

IF (istart.eq.1) THEN
 XI1 = -1d0; XI2 = -1d0; XI3 = -1d0
END IF
IF (istart.eq.2) THEN
 XI1 = +1d0; XI2 = -1d0; XI3 = -1d0
END IF
IF (istart.eq.3) THEN
 XI1 = +1d0; XI2 = +1d0; XI3 = -1d0
END IF
IF (istart.eq.4) THEN
 XI1 = -1d0; XI2 = +1d0; XI3 = -1d0
END IF
IF (istart.eq.5) THEN
 XI1 = -1d0; XI2 = -1d0; XI3 = 1d0
END IF
IF (istart.eq.6) THEN
 XI1 = +1d0; XI2 = -1d0; XI3 = 1d0
END IF
IF (istart.eq.7) THEN
 XI1 = +1d0; XI2 = +1d0; XI3 = 1d0
END IF
IF (istart.eq.8) THEN
 XI1 = -1d0; XI2 = +1d0; XI3 = 1d0
END IF

!XI1 = 0d0; XI2 = 0d0; XI3 = 0d0

DO iter = 1,80

  DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
  DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
  DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
  DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
  DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
  DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
  DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
  DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
  DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
!   DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
!        -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
!        +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))

  XX = DJ11 + DJAC(1,1)*XI1 + DJ31*XI2 + DJ41*XI3 + DJ71*XI2*XI3
  YY = DJ12 + DJ22*XI1 + DJAC(2,2)*XI2 + DJ42*XI3 + DJ62*XI1*XI3
  ZZ = DJ13 + DJ23*XI1 + DJ33*XI2 + DJAC(3,3)*XI3 + DJ53*XI1*XI2

  CALL M33INV (DJAC, DJACI, OK_FLAG)

!  dFact = 0.33d0/ABS(DETJ)
  dFact = -0.33d0
  diff1 = -dFact*(DJACI(1,1)*(P(1)-XX) + DJACI(1,2)*(P(2)-YY) + DJACI(1,3)*(P(3)-ZZ))
  diff2 = -dFact*(DJACI(2,1)*(P(1)-XX) + DJACI(2,2)*(P(2)-YY) + DJACI(2,3)*(P(3)-ZZ))
  diff3 = -dFact*(DJACI(3,1)*(P(1)-XX) + DJACI(3,2)*(P(2)-YY) + DJACI(3,3)*(P(3)-ZZ))
  XI1 = XI1 + diff1
  XI2 = XI2 + diff2
  XI3 = XI3 + diff3
  daux = diff1*diff1 + diff2*diff2 + diff3*diff3

!   WRITE(*,'(I3,9ES12.2)') iter,XI1,XI2,XI3, XX,YY,ZZ, diff1,diff2,diff3
  IF (iter.gt.4.and.(ABS(XI1).GT.1.5d0.OR.ABS(XI2).GT.1.5d0.OR.ABS(XI3).GT.1.5d0)) GOTO 2
  IF (iter.gt.4.and.daux.lt.0.00001d0) GOTO 1

END DO

1 CONTINUE
IF (ABS(XI1).GT.1.01d0.OR.ABS(XI2).GT.1.01d0.OR.ABS(XI3).GT.1.01d0) GOTO 2

BF = .TRUE.
PR = [XI1,XI2,XI3]
RETURN

2 CONTINUE
! IF (iP.EQ.142) WRITE(*,*) '-- ',iter,XI1,XI2,XI3,daux
!pause

END SUBROUTINE GetPointFromElement

!========================================================================================
!                                Sub: GetMeshBounds
!========================================================================================
SUBROUTINE GetMeshBounds(dcorvg,nvtx,xmin,xmax,ymin,ymax,zmin,zmax)
USE PP3D_MPI, ONLY : myid, showid, master, COMM_Minimum, COMM_Maximum, myMPI_Barrier

  IMPLICIT NONE
  REAL*8  dcorvg(3,*)
  real*8  xmin,xmax,ymin,ymax,zmin,zmax
  integer  i,nvtx

  xmin =1D99
  ymin = xmin
  zmin = xmin

  xmax = -1D99
  ymax = xmax
  zmax = xmax

  if (myid .ne. master ) then
    ! Get the local minima/maxima
    xmin = dcorvg(1,1)
    xmax = dcorvg(1,1)

    ymin = dcorvg(2,1)
    ymax = dcorvg(2,1)

    zmin = dcorvg(3,1)
    zmax = dcorvg(3,1)

    do i=1,nvtx
      xmin = min(xmin,dcorvg(1,i))
      xmax = max(xmax,dcorvg(1,i))

      ymin = min(ymin,dcorvg(2,i))
      ymax = max(ymax,dcorvg(2,i))

      zmin = min(zmin,dcorvg(3,i))
      zmax = max(zmax,dcorvg(3,i))
    end do
  end if

  call COMM_Maximum(xmax)
  call COMM_Maximum(ymax)
  call COMM_Maximum(zmax)

  call COMM_Minimum(xmin)
  call COMM_Minimum(ymin)
  call COMM_Minimum(zmin)

end subroutine GetMeshBounds

!========================================================================================
!                                    Sub: M33Inv
!========================================================================================
subroutine M33INV (A, AINV, OK_FLAG)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
LOGICAL, INTENT(OUT) :: OK_FLAG

DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
DOUBLE PRECISION :: DET
DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


DET =   A(1,1)*A(2,2)*A(3,3)  &
      - A(1,1)*A(2,3)*A(3,2)  &
      - A(1,2)*A(2,1)*A(3,3)  &
      + A(1,2)*A(2,3)*A(3,1)  &
      + A(1,3)*A(2,1)*A(3,2)  &
      - A(1,3)*A(2,2)*A(3,1)

 IF (ABS(DET) .LE. EPS) THEN
    AINV = 0.0D0
    OK_FLAG = .FALSE.
    RETURN
 END IF

 COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
 COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
 COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
 COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
 COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
 COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
 COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
 COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
 COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

 AINV = TRANSPOSE(COFACTOR) / DET

 OK_FLAG = .TRUE.

 RETURN

end subroutine M33INV

!========================================================================================
!                                Sub: Extract_Particle
!========================================================================================
SUBROUTINE Extract_Particle(dcorvg,kvert,kedge,karea,kel,nkel,&
                           VelU,VelV,VelW,nvt,net,nat,nel,tDelta)
USE PP3D_MPI, ONLY : myid,master,showid,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,RECVK_myMPI,SENDK_myMPI,MPI_COMM_WORLD

USE types, ONLY : myCompleteSet,myActiveSet,nActiveSet,myExchangeSet,myLostSet,&
                         nExchangeSet,nStartActiveSet,nLostSet

IMPLICIT NONE
REAL*8 dcorvg(3,*),VelU(*),VelV(*),VelW(*),point(3),tLevel,tDelta
INTEGER nkel,nvt,net,nat,nel
INTEGER kel(nkel,*),kvert(8,*),kedge(12,*),karea(6,*)
INTEGER KDFG(27),KDFL(27)
INTEGER i,iPoint,jPoint,ivt,jel,iel,iFoundElem,pID,IERR
REAL*8 dist,P8(3,8),pointR(3),LocVel(3,27)
LOGICAL :: bFound,bOut=.FALSE.
INTEGER :: nIter = 100, iIter, iLoc,iParticel,iActiveParticel
INTEGER iMonitor(40),iAux
REAL*8  dMonitor(40),dAux
INTEGER nXX,kk,iMon
PARAMETER (nXX = 40)
REAL*8 , ALLOCATABLE :: dSendVector(:),dRecvVector(:)
INTEGER, ALLOCATABLE :: iRecvVector(:)

ALLOCATE(dSendVector(nExchangeSet),iRecvVector(nExchangeSet),dRecvVector(nExchangeSet))
dSendVector = 1d30
iRecvVector = 0

IF (myid.ne.0) THEN

DO iParticel = 1,nExchangeSet

iIter = 0

point  = myExchangeSet(iParticel)%coor
tLevel = myExchangeSet(iParticel)%time

55 CONTINUE

bFound = .FALSE.

dMonitor = 1d30
iMonitor = 0

DO i=1,nvt

 dist = sqrt((dcorvg(1,i)-point(1))**2d0 + (dcorvg(2,i)-point(2))**2d0 + (dcorvg(3,i)-point(3))**2d0)
 IF (dist.lt.dMonitor(nXX)) THEN
  dMonitor(nXX) = dist
  iMonitor(nXX) = i
  DO kk=nXX-1,1,-1
   IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
    iAux = iMonitor(kk)
    dAux = dMonitor(kk)
    iMonitor(kk) = iMonitor(kk+1)
    dMonitor(kk) = dMonitor(kk+1)
    iMonitor(kk+1) = iAux
    dMonitor(kk+1) = dAux
   END IF
  END DO
 END IF

END DO


DO iMon = 1,nXX

iPoint = iMonitor(iMon)

DO iel = 1,nkel

 jel = kel(iel,iPoint)

 if (jel.ne.0) then
  DO i=1,8
   ivt = kvert(i,jel)
   P8(:,i) = dcorvg(:,ivt)
  END DO
 ELSE
  EXIT
 END IF

 dist = 1d30
 jPoint = 0

 DO i=1,8
  daux = sqrt((P8(1,i)-point(1))**2d0 + (P8(2,i)-point(2))**2d0 + (P8(3,i)-point(3))**2d0)
  IF (daux.lt.dist) THEN
   dist = daux
   jPoint = i
  END IF
 END DO

 CALL GetPointFromElement(P8,point,jPoint,bFound,pointR,iParticel)

 IF (bFound) THEN
  GOTO 1
 END IF

END DO

END DO

IF (.not.bFound) GOTO 2

1 CONTINUE

dSendVector(iParticel) = MAX(pointR(1),pointR(2),pointR(3))

2 CONTINUE

END DO

END IF

IF (myid.eq.0) THEN
 DO pID=1,subnodes
   CALL RECVD_myMPI(dRecvVector,nExchangeSet,pID)
   DO i=1,nExchangeSet
    IF (dRecvVector(i).lt.dSendVector(i)) THEN
     dSendVector(i) = dRecvVector(i)
     iRecvVector(i) = pID
    END IF
   END DO
 END DO

ELSE
  CALL SENDD_myMPI(dSendVector,nExchangeSet,0)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF (myid.eq.0) THEN
 DO pID=1,subnodes
   CALL SENDK_myMPI(iRecvVector,nExchangeSet,pID)
 END DO
ELSE
  CALL RECVK_myMPI(iRecvVector,nExchangeSet,0)
END IF

nStartActiveSet = nActiveSet
iActiveParticel = 0

DO iParticel = 1,nExchangeSet
 IF (iRecvVector(iParticel).eq.myid) THEN
  iActiveParticel = iActiveParticel + 1
  myActiveSet(nStartActiveSet+iActiveParticel) = myExchangeSet(iParticel)
 END IF
 IF (iRecvVector(iParticel).eq.0) THEN
  nLostSet = nLostSet + 1
  myLostSet(nLostSet) = myExchangeSet(iParticel)
 END IF
END DO

DEALLOCATE(dSendVector,iRecvVector,dRecvVector)

nActiveSet = nStartActiveSet+iActiveParticel
IF (bOut) WRITE(*,*) "My Active Set of Particles after exchange: ", nActiveSet
! pause

END SUBROUTINE Extract_Particle

!========================================================================================
!                            Sub: Exchange_Particle
!========================================================================================
SUBROUTINE Exchange_Particle(nSum)
USE PP3D_MPI, ONLY : myid,master,showid,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,RECVK_myMPI,SENDK_myMPI,MPI_COMM_WORLD
USE types, ONLY : myActiveSet,myExchangeSet,nActiveSet,nExchangeSet,nStartActiveSet

IMPLICIT NONE
INTEGER nSum
INTEGER iN,pID,i,ierr
REAL*8,  ALLOCATABLE :: dAux(:)
INTEGER, ALLOCATABLE :: iAux(:)

ALLOCATE(dAux(4*nSum))
ALLOCATE(iAux(1*nSum))

!  WRITE(*,*) "myid,nsum ",myid,nsum

IF (myid.eq.MASTER) THEN
 nExchangeSet = 0
 DO pID=1,subnodes
  CALL RECVI_myMPI(iN,pID)

  IF (iN.ne.0) THEN
   CALL RECVK_myMPI(iAux,  iN,pID)
   CALL RECVD_myMPI(dAux,4*iN,pID)
   DO i=1,iN
    myExchangeSet(nExchangeSet+i)%coor(1) = dAux(4*(i-1)+1)
    myExchangeSet(nExchangeSet+i)%coor(2) = dAux(4*(i-1)+2)
    myExchangeSet(nExchangeSet+i)%coor(3) = dAux(4*(i-1)+3)
    myExchangeSet(nExchangeSet+i)%time    = dAux(4*(i-1)+4)
    myExchangeSet(nExchangeSet+i)%indice  = iAux(i)
   END DO
   nExchangeSet = nExchangeSet + iN

  END IF
 END DO
ELSE
  CALL SENDI_myMPI(nExchangeSet,0)
  IF (nExchangeSet.ne.0) THEN
   DO i=1,nExchangeSet
    dAux(4*(i-1)+1:4*(i-1)+4) = [myExchangeSet(i)%coor(1),myExchangeSet(i)%coor(2),myExchangeSet(i)%coor(3),myExchangeSet(i)%time]
    iAux(i)                   =  myExchangeSet(i)%indice
   END DO
  CALL SENDK_myMPI(iAux,  nExchangeSet,0)
  CALL SENDD_myMPI(dAux,4*nExchangeSet,0)

  END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF (myid.eq.MASTER) THEN

 DO i=1,nSum
  dAux(4*(i-1)+1:4*(i-1)+4) = [myExchangeSet(i)%coor(1),myExchangeSet(i)%coor(2),myExchangeSet(i)%coor(3),myExchangeSet(i)%time]
  iAux(i)                   =  myExchangeSet(i)%indice
 END DO

 DO pID=1,subnodes
  CALL SENDD_myMPI(dAux,4*nSum,pID)
  CALL SENDK_myMPI(iAux,  nSum,pID)
 END DO

ELSE

  CALL RECVD_myMPI(dAux,4*nSum,0)
  CALL RECVK_myMPI(iAux,  nSum,0)

  DO i=1,nSum
   myExchangeSet(i)%coor(1) = dAux(4*(i-1)+1)
   myExchangeSet(i)%coor(2) = dAux(4*(i-1)+2)
   myExchangeSet(i)%coor(3) = dAux(4*(i-1)+3)
   myExchangeSet(i)%time    = dAux(4*(i-1)+4)
   myExchangeSet(i)%indice  = iAux(i)
  END DO
  nExchangeSet = nSum

END IF

DEALLOCATE(dAux,iAux)

END SUBROUTINE Exchange_Particle

!========================================================================================
!                            Sub: Gather_Particle
!========================================================================================
SUBROUTINE Gather_Particle(nSum)
USE PP3D_MPI, ONLY : myid,master,showid,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,RECVK_myMPI,SENDK_myMPI,MPI_COMM_WORLD
  USE types, ONLY : myActiveSet,myExchangeSet,nActiveSet,nExchangeSet,myCompleteSet,nCompleteSet

IMPLICIT NONE
INTEGER nSum
INTEGER iN,pID,i,ierr
REAL*8,  ALLOCATABLE :: dAux(:)
INTEGER, ALLOCATABLE :: iAux(:)

ALLOCATE(dAux(4*nSum))
ALLOCATE(iAux(1*nSum))

!  WRITE(*,*) "myid,nsum ",myid,nsum

IF (myid.eq.MASTER) THEN
 nCompleteSet = 0
 DO pID=1,subnodes
  CALL RECVI_myMPI(iN,pID)

  IF (iN.ne.0) THEN
   CALL RECVK_myMPI(iAux,  iN,pID)
   CALL RECVD_myMPI(dAux,4*iN,pID)
   DO i=1,iN
    myCompleteSet(nCompleteSet+i)%coor(1) = dAux(4*(i-1)+1)
    myCompleteSet(nCompleteSet+i)%coor(2) = dAux(4*(i-1)+2)
    myCompleteSet(nCompleteSet+i)%coor(3) = dAux(4*(i-1)+3)
    myCompleteSet(nCompleteSet+i)%time    = dAux(4*(i-1)+4)
    myCompleteSet(nCompleteSet+i)%indice  = iAux(i)
   END DO
   nCompleteSet = nCompleteSet + iN

  END IF
 END DO
ELSE
  CALL SENDI_myMPI(nActiveSet,0)
  IF (nActiveSet.ne.0) THEN
   DO i=1,nActiveSet
    dAux(4*(i-1)+1:4*(i-1)+4) = [myActiveSet(i)%coor(1),myActiveSet(i)%coor(2),myActiveSet(i)%coor(3),myActiveSet(i)%time]
    iAux(i)                   =  myActiveSet(i)%indice
   END DO
  CALL SENDK_myMPI(iAux,  nActiveSet,0)
  CALL SENDD_myMPI(dAux,4*nActiveSet,0)
  END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
! if (myid.eq.1) WRITE(*,*) 'so far so good ....'


IF (myid.eq.MASTER) THEN

 DO i=1,nSum
  dAux(4*(i-1)+1:4*(i-1)+4) = [myCompleteSet(i)%coor(1),myCompleteSet(i)%coor(2),myCompleteSet(i)%coor(3),myCompleteSet(i)%time]
  iAux(i)                   =  myCompleteSet(i)%indice
 END DO

 DO pID=1,subnodes
  CALL SENDD_myMPI(dAux,4*nSum,pID)
  CALL SENDK_myMPI(iAux,  nSum,pID)
 END DO

ELSE

  CALL RECVD_myMPI(dAux,4*nSum,0)
  CALL RECVK_myMPI(iAux,  nSum,0)

  DO i=1,nSum
   myCompleteSet(i)%coor(1) = dAux(4*(i-1)+1)
   myCompleteSet(i)%coor(2) = dAux(4*(i-1)+2)
   myCompleteSet(i)%coor(3) = dAux(4*(i-1)+3)
   myCompleteSet(i)%time    = dAux(4*(i-1)+4)
   myCompleteSet(i)%indice  = iAux(i)
  END DO
  nCompleteSet = nSum

END IF

DEALLOCATE(dAux,iAux)

END SUBROUTINE Gather_Particle

!========================================================================================
!                           Sub: Search_Particle
!========================================================================================
SUBROUTINE Search_Particle(dcorvg,kvert,kedge,karea,kel,nkel,&
                           VelU,VelV,VelW,nvt,net,nat,nel,tDelta)
USE types, ONLY : myCompleteSet,myActiveSet,nActiveSet

IMPLICIT NONE
REAL*8 dcorvg(3,*),VelU(*),VelV(*),VelW(*),point(3),tLevel,tDelta
INTEGER nkel,nvt,net,nat,nel
INTEGER kel(nkel,*),kvert(8,*),kedge(12,*),karea(6,*)
INTEGER KDFG(27),KDFL(27)
INTEGER i,iPoint,jPoint,ivt,jel,iel,iFoundElem
REAL*8 dist,P8(3,8),pointR(3),LocVel(3,27)
LOGICAL :: bFound,bOut=.FALSE.
INTEGER :: nIter = 100, iIter, iLoc,iParticel,iActiveParticel
INTEGER iMonitor(10),iAux
REAL*8  dMonitor(10),dAux
INTEGER nXX,kk,iMon
PARAMETER (nXX = 10)

iActiveParticel = 0

DO iParticel = 1,SIZE(myCompleteSet)

iIter = 0

point  = myCompleteSet(iParticel)%coor
tLevel = myCompleteSet(iParticel)%time

55 CONTINUE

dist = 1d30
iPoint = 0
bFound = .FALSE.

! DO jel=1,nel
!
!   DO i=1,8
!    ivt = kvert(i,jel)
!    P8(:,i) = dcorvg(:,ivt)
!    daux = sqrt((P8(1,i)-point(1))**2d0 + (P8(2,i)-point(2))**2d0 + (P8(3,i)-point(3))**2d0)
!    IF (daux.lt.dist) THEN
!     dist = daux
!     jPoint = i
!    END IF
!   END DO

dMonitor = 1d30
iMonitor = 0

DO i=1,nvt

 dist = sqrt((dcorvg(1,i)-point(1))**2d0 + (dcorvg(2,i)-point(2))**2d0 + (dcorvg(3,i)-point(3))**2d0)
 IF (dist.lt.dMonitor(nXX)) THEN
  dMonitor(nXX) = dist
  iMonitor(nXX) = i
  DO kk=nXX-1,1,-1
   IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
    iAux = iMonitor(kk)
    dAux = dMonitor(kk)
    iMonitor(kk) = iMonitor(kk+1)
    dMonitor(kk) = dMonitor(kk+1)
    iMonitor(kk+1) = iAux
    dMonitor(kk+1) = dAux
   END IF
  END DO
 END IF

END DO


DO iMon = 1,nXX

iPoint = iMonitor(iMon)

DO iel = 1,nkel

 jel = kel(iel,iPoint)

 if (jel.ne.0) then
  DO i=1,8
   ivt = kvert(i,jel)
   P8(:,i) = dcorvg(:,ivt)
  END DO
 ELSE
  EXIT
 END IF

 dist = 1d30
 jPoint = 0

 DO i=1,8
  daux = sqrt((P8(1,i)-point(1))**2d0 + (P8(2,i)-point(2))**2d0 + (P8(3,i)-point(3))**2d0)
  IF (daux.lt.dist) THEN
   dist = daux
   jPoint = i
  END IF
 END DO

 CALL GetPointFromElement(P8,point,jPoint,bFound,pointR,iParticel)

 IF (bFound) THEN
!   WRITE(*,*) 'yes',jel,kvert(jPoint,jel)
  GOTO 1
 END IF

END DO !iel
END DO !iMon

IF (.not.bFound) GOTO 2

1 CONTINUE

iActiveParticel = iActiveParticel + 1
myActiveSet(iActiveParticel) = myCompleteSet(iParticel)

2 CONTINUE

END DO

nActiveSet = iActiveParticel
! WRITE(*,*) "My Active Set of Particles: ", iActiveParticel, bFound,iPoint
! WRITE(*,*) "rel elems ", kel(1:nkel,iPoint)
! pause

END SUBROUTINE Search_Particle

!========================================================================================
!                           Sub: Search_Particle
!========================================================================================
subroutine Move_Particle(dcorvg,kvert,kedge,karea,kel,nkel,&
                         Vel0U,Vel0V,Vel0W,Vel1U,Vel1V,Vel1W,nvt,net,nat,nel,tDelta,tStart)
USE PP3D_MPI, ONLY : myid
USE types, ONLY : myActiveSet,myExchangeSet,nActiveSet,nExchangeSet,nStartActiveSet

IMPLICIT NONE
REAL*8 dcorvg(3,*),point(3),tLevel,tDelta,tStart
REAL*8 Vel0U(*),Vel0V(*),Vel0W(*),Vel1U(*),Vel1V(*),Vel1W(*)
INTEGER nkel,nvt,net,nat,nel
INTEGER kel(nkel,*),kvert(8,*),kedge(12,*),karea(6,*)
INTEGER KDFG(27),KDFL(27)
INTEGER i,iPoint,jPoint,ivt,jel,iel,iFoundElem
REAL*8 dist,P8(3,8),pointR(3),LocVel0(3,27),LocVel1(3,27)
LOGICAL :: bFound,bOut=.FALSE.
INTEGER :: nIter = 100, iIter, iLoc,iParticel,iLostParticel,iActiveParticel
INTEGER iMonitor(40),iAux
REAL*8  dMonitor(40),dAux
INTEGER nXX,kk,iMon
PARAMETER (nXX = 40)

iLostParticel = 0
iActiveParticel = 0

DO iParticel = nStartActiveSet+1,nActiveSet

iIter = 0

point  = myActiveSet(iParticel)%coor
tLevel = myActiveSet(iParticel)%time

55 CONTINUE

dist = 1d30
iPoint = 0
bFound = .FALSE.

dMonitor = 1d30
iMonitor = 0

DO i=1,nvt

 dist = sqrt((dcorvg(1,i)-point(1))**2d0 + (dcorvg(2,i)-point(2))**2d0 + (dcorvg(3,i)-point(3))**2d0)
 IF (dist.lt.dMonitor(nXX)) THEN
  dMonitor(nXX) = dist
  iMonitor(nXX) = i
  DO kk=nXX-1,1,-1
   IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
    iAux = iMonitor(kk)
    dAux = dMonitor(kk)
    iMonitor(kk) = iMonitor(kk+1)
    dMonitor(kk) = dMonitor(kk+1)
    iMonitor(kk+1) = iAux
    dMonitor(kk+1) = dAux
   END IF
  END DO
 END IF

END DO


DO iMon = 1,nXX

iPoint = iMonitor(iMon)

DO iel = 1,nkel

 jel = kel(iel,iPoint)

 IF (jel.ne.0) then
  DO i=1,8
   ivt = kvert(i,jel)
   P8(:,i) = dcorvg(:,ivt)
  END DO
 ELSE
  EXIT
 END IF

 dist = 1d30
 jPoint = 0

 DO i=1,8
  daux = sqrt((P8(1,i)-point(1))**2d0 + (P8(2,i)-point(2))**2d0 + (P8(3,i)-point(3))**2d0)
  IF (daux.lt.dist) THEN
   dist = daux
   jPoint = i
  END IF
 END DO

 CALL GetPointFromElement(P8,point,jPoint,bFound,pointR,iParticel)

 IF (bFound) THEN

  CALL NDFGL(JEL,1,13,KVERT,KEDGE,KAREA,KDFG,KDFL)
  DO i=1,27
   ivt  = KDFG(i)
   iLoc = KDFL(i)
   LocVel0(:,iLoc) = [Vel0U(ivt),Vel0V(ivt),Vel0W(ivt)]
   LocVel1(:,iLoc) = [Vel1U(ivt),Vel1V(ivt),Vel1W(ivt)]
  END DO

  IF (bOut) WRITE(*,'(A,I0,3E14.4,A)') '------------   point ', myActiveSet(iParticel)%indice, point, '  ------------'
  CALL MovePointInElement(P8,point,pointR,LocVel0,LocVel1,tLevel,tDelta,tStart,iParticel,jel)
  iFoundElem = jel
  GOTO 1
 END IF

END DO
END DO

IF (.not.bFound) THEN
! myActiveSet(iParticel)%coor = point
! myActiveSet(iParticel)%time = tLevel
 iLostParticel = iLostParticel + 1
 myExchangeSet(iLostParticel)%coor   = point
 myExchangeSet(iLostParticel)%time   = tLevel
 myExchangeSet(iLostParticel)%indice = myActiveSet(iParticel)%indice
!  IF (myExchangeSet(iLostParticel)%coor(2).GT.0d0) THEN
!   WRITE(*,'(A,I,4E16.8)') 'problem',myActiveSet(iParticel)%indice,sqrt((myExchangeSet(iLostParticel)%coor(1))**2d0 + (myExchangeSet(iLostParticel)%coor(2)-16.7d0)**2d0),myExchangeSet(iLostParticel)%coor
!  ELSE
!   WRITE(*,'(A,I,4E16.8)') 'problem',myActiveSet(iParticel)%indice,sqrt((myExchangeSet(iLostParticel)%coor(1))**2d0 + (myExchangeSet(iLostParticel)%coor(2)+16.7d0)**2d0),myExchangeSet(iLostParticel)%coor
!  END IF
 GOTO 2
END IF

1 CONTINUE

iIter = iIter + 1
IF (tLevel.lt.tDelta.AND.iIter.LT.nIter) GOTO 55

IF (bFound.AND.bOut) WRITE(*,'(A,I0,3E14.4,A)') '------------   point ', myActiveSet(iParticel)%indice, point, '  ------------'
IF (bFound.AND.bOut) WRITE(*,'(A,I0,I0,3E14.4,A)') '-------  point ', myActiveSet(iParticel)%indice,iIter,point, '  has been successfully treated -----------'

iActiveParticel = iActiveParticel + 1
myActiveSet(nStartActiveSet+iActiveParticel)%coor   = point
myActiveSet(nStartActiveSet+iActiveParticel)%time   = tLevel
myActiveSet(nStartActiveSet+iActiveParticel)%indice = myActiveSet(iParticel)%indice

2 CONTINUE

END DO

nActiveSet      = nStartActiveSet+iActiveParticel
nExchangeSet    = iLostParticel
nStartActiveSet = nActiveSet

IF (bOut) WRITE(*,*) 'nActiveSet,nExchangeSet: ',nActiveSet,nExchangeSet

! pause
! DO iParticel=nActiveSet+1,SIZE(myActiveSet)
!  myActiveSet(iActiveParticel)%coor   = 0d0
!  myActiveSet(iActiveParticel)%time   = 0d0
!  myActiveSet(iActiveParticel)%indice = 0
! END DO

! pause

END SUBROUTINE Move_Particle

!========================================================================================
!                           Sub: Search_Particle
!========================================================================================
SUBROUTINE MovePointInElement(P8,P,PR,DV0,DV1,tLevel,tEnd,tStart,iP,iE)
!USE Particle, ONLY : myParticleParam
USE PP3D_MPI, ONLY : myid
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
PARAMETER (Q8=0.125D0)
LOGICAL :: BOUT=.FALSE.
REAL*8  P8(3,8), P(3),PR(3),DV0(3,27),DV1(3,27),tLevel,tEnd,tStart
REAL*8 DJACI(3,3),DJAC(3,3)
REAL*8 DBAS(27),XI1,XI2,XI3,DV_Loc(3),P_new(3)
REAL*8 :: DV_Rot(3),dVV(3)
REAL*8 :: DeltaT=1d0,new_timestep
REAL*8 :: PI=3.141592654d0,dAlpha,XB,YB,ZB
REAL*8 :: RK_Velo(3,4)
logical :: ok_flag
dPI = 4d0*dATAN(1d0)

!dFreq = (myParticleParam%f/60d0)
!
!iMove = 0
!dVelo = 0d0
!
!XI1 = PR(1)
!XI2 = PR(2)
!XI3 = PR(3)
!
!IF (bOut) WRITE(*,'(A,7E14.4)') 'my original position: ', XI1,XI2,XI3,P,tLevel
!
!CALL RETURN_Velo()
!
!DJ11=( P8(1,1)+P8(1,2)+P8(1,3)+P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
!DJ12=( P8(2,1)+P8(2,2)+P8(2,3)+P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
!DJ13=( P8(3,1)+P8(3,2)+P8(3,3)+P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
!DJ21=(-P8(1,1)+P8(1,2)+P8(1,3)-P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
!DJ22=(-P8(2,1)+P8(2,2)+P8(2,3)-P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
!DJ23=(-P8(3,1)+P8(3,2)+P8(3,3)-P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
!DJ31=(-P8(1,1)-P8(1,2)+P8(1,3)+P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
!DJ32=(-P8(2,1)-P8(2,2)+P8(2,3)+P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
!DJ33=(-P8(3,1)-P8(3,2)+P8(3,3)+P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
!DJ41=(-P8(1,1)-P8(1,2)-P8(1,3)-P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
!DJ42=(-P8(2,1)-P8(2,2)-P8(2,3)-P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
!DJ43=(-P8(3,1)-P8(3,2)-P8(3,3)-P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
!DJ51=( P8(1,1)-P8(1,2)+P8(1,3)-P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
!DJ52=( P8(2,1)-P8(2,2)+P8(2,3)-P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
!DJ53=( P8(3,1)-P8(3,2)+P8(3,3)-P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8
!DJ61=( P8(1,1)-P8(1,2)-P8(1,3)+P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
!DJ62=( P8(2,1)-P8(2,2)-P8(2,3)+P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
!DJ63=( P8(3,1)-P8(3,2)-P8(3,3)+P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
!DJ71=( P8(1,1)+P8(1,2)-P8(1,3)-P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
!DJ72=( P8(2,1)+P8(2,2)-P8(2,3)-P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
!DJ73=( P8(3,1)+P8(3,2)-P8(3,3)-P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
!DJ81=(-P8(1,1)+P8(1,2)-P8(1,3)+P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
!DJ82=(-P8(2,1)+P8(2,2)-P8(2,3)+P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
!DJ83=(-P8(3,1)+P8(3,2)-P8(3,3)+P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8
!
! DJAC(1,1)=DJ21
! DJAC(1,2)=DJ31
! DJAC(1,3)=DJ41
! DJAC(2,1)=DJ22
! DJAC(2,2)=DJ32
! DJAC(2,3)=DJ42
! DJAC(3,1)=DJ23
! DJAC(3,2)=DJ33
! DJAC(3,3)=DJ43
! DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
!      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
!      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
!
! dElemSize = (ABS(DETJ)**0.3333d0)
! dVelo = SQRT(DV_Loc(1)*DV_Loc(1) + DV_Loc(2)*DV_Loc(2) + DV_Loc(3)*DV_Loc(3))
! DeltaT = myParticleParam%hSize*dElemSize/dVelo
!! DeltaT = 0.001d0
! IF ((tLevel + DeltaT).GT.tEnd) THEN
!  DeltaT = tEnd - tLevel
! END IF
! dFracTime = (tLevel + 0.5d0*DeltaT - tStart)/(tEnd - tStart)
!
!55 CONTINUE
!
!RK_Velo(:,1) = DV_loc
!
!CALL FindPoint(2d0/2d0,DV_loc)  ! --> XI1,XI2,XI3
!CALL RETURN_Velo()              ! --> DV_loc
!! RK_Velo(:,2) = DV_loc
!!
!! CALL FindPoint(1d0/2d0,DV_loc)  ! --> XI1,XI2,XI3
!! CALL RETURN_Velo()              ! --> DV_loc
!! RK_Velo(:,3) = DV_loc
!!
!! CALL FindPoint(2d0/2d0,DV_loc)  ! --> XI1,XI2,XI3
!! CALL RETURN_Velo()              ! --> DV_loc
!! RK_Velo(:,4) = DV_loc
!!
!! DV_Loc = (1d0/6d0)*(RK_Velo(:,1) + 2d0*RK_Velo(:,2) + 2d0*RK_Velo(:,3) + RK_Velo(:,4))
!! ! DV_Loc = (0.25d0*RK_Velo(:,1) + 0.75d0*RK_Velo(:,2))
!!
!! CALL FindPoint(1d0/1d0,DV_loc)  ! --> XI1,XI2,XI3
!P(1) = P_New(1)
!P(2) = P_New(2)
!P(3) = P_New(3)
!! CALL RETURN_Velo()              ! --> DV_loc
!
!tLevel = tLevel + DeltaT
!iMove  = iMove + 1
!
!IF ((ABS(XI1).GT.1.01d0.OR.ABS(XI2).GT.1.01d0.OR.ABS(XI3).GT.1.01d0).OR.(tLevel.ge.tEnd)) THEN
! ! Reached the end of the element -> leave the routine
! GOTO 5
!ELSE
!
!! IF ((tLevel + DeltaT).GT.tEnd) THEN
!!  DeltaT = tEnd - tLevel
!! END IF
!! We are still in the element - go back and continue marching in the element
!GOTO 55
!
!END IF
!
!5 CONTINUE
!
!! pause
!
!XX = P(1)
!YY = P(2)
!ZZ = P(3)
!
!dist = SQRT(YY*YY + XX*XX)
!
!IF (dist.ge.myParticleParam%dEps1*myParticleParam%D_Out*0.5d0) THEN
! dFrac = dist/(myParticleParam%dEps2*myParticleParam%D_Out*0.5d0)
! XX = XX/dFrac
! YY = YY/dFrac
!
! P(1) = XX
! P(2) = YY
!! WRITE(*,*) '... problem is being fixed .. '
!END IF
!
!! IF (ZZ.ge.640.0d0) THEN
!!  P(3) = 1d3
!! END IF
!
CONTAINS

SUBROUTINE RETURN_Velo()

! Loop not neccesary because nothing depends on i?
!DO i=1,27
 DBAS( 1)=-Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
 DBAS( 2)= Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
 DBAS( 3)=-Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
 DBAS( 4)=Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
 DBAS( 5)=Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
 DBAS( 6)=-Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
 DBAS( 7)=Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
 DBAS( 8)=-Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
 DBAS( 9)=Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
 DBAS(10)=-Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
 DBAS(11)=-Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
 DBAS(12)=Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
 DBAS(13)=Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
 DBAS(14)=-Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
 DBAS(15)=Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
 DBAS(16)=-Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
 DBAS(17)=-Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
 DBAS(18)=Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
 DBAS(19)=Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
 DBAS(20)=-Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
 DBAS(21)= -Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
 DBAS(22)= -Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
 DBAS(23)= Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
 DBAS(24)= Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
 DBAS(25)= -Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
 DBAS(26)= Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
 DBAS(27)= Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
!END DO

DV_loc = 0d0
DO I=1,27
 DV_loc(1) = DV_loc(1) + DBAS(I)*DV0(1,I)
 DV_loc(2) = DV_loc(2) + DBAS(I)*DV0(2,I)
 DV_loc(3) = DV_loc(3) + DBAS(I)*DV0(3,I)
END DO

END SUBROUTINE RETURN_Velo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE FindPoint(fFact,dLocalVelo)
REAL*8 fFact,dLocalVelo(3)

DV_Rot = [-2d0*dPI*dFreq*P(2), 2d0*dPI*dFreq*P(1),0d0]
dtheta = datan(P(2)/P(1))
IF (P(1).lt.0d0) THEN
 dtheta = dtheta + dPI
END IF

dRR    = DSQRT(P(1)**2d0 + P(2)**2d0)

dVV  = dLocalVelo + DV_Rot

dU_r     =  dVV(1)*cos(dtheta) + dVV(2)*sin(dtheta)
dU_theta =(-dVV(1)*sin(dtheta) + dVV(2)*cos(dtheta))/dRR
dU_z     =  dVV(3)

dAlpha = fFact*DeltaT*dU_theta
dRho   = fFact*DeltaT*dU_r
dZ     = fFact*DeltaT*dU_z

!!!! Angular transformation
XB = P(1)
YB = P(2)
ZB = P(3)

P_New(1) = XB*cos(dAlpha) - YB*sin(dAlpha)
P_New(2) = XB*sin(dAlpha) + YB*cos(dAlpha)

!!!! Radial transformation
XB = P_New(1)
YB = P_New(2)

dRho0 = DSQRT(XB**2d0 + YB**2d0)
dRho1 = dRho0 + dRho

P_New(1) = XB*dRho1/dRho0
P_New(2) = YB*dRho1/dRho0

!!!! Axial transformation
P_New(3) = ZB + dZ

DO iter = 1,80

  DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
  DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
  DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
  DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
  DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
  DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
  DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
  DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
  DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2

  XX = DJ11 + DJAC(1,1)*XI1 + DJ31*XI2 + DJ41*XI3 + DJ71*XI2*XI3
  YY = DJ12 + DJ22*XI1 + DJAC(2,2)*XI2 + DJ42*XI3 + DJ62*XI1*XI3
  ZZ = DJ13 + DJ23*XI1 + DJ33*XI2 + DJAC(3,3)*XI3 + DJ53*XI1*XI2

  CALL M33INV (DJAC, DJACI, OK_FLAG)

  dFact = -0.33d0
  diff1 = -dFact*(DJACI(1,1)*(P_New(1)-XX) + DJACI(1,2)*(P_New(2)-YY) + DJACI(1,3)*(P_New(3)-ZZ))
  diff2 = -dFact*(DJACI(2,1)*(P_New(1)-XX) + DJACI(2,2)*(P_New(2)-YY) + DJACI(2,3)*(P_New(3)-ZZ))
  diff3 = -dFact*(DJACI(3,1)*(P_New(1)-XX) + DJACI(3,2)*(P_New(2)-YY) + DJACI(3,3)*(P_New(3)-ZZ))
  XI1 = XI1 + diff1
  XI2 = XI2 + diff2
  XI3 = XI3 + diff3
  daux = diff1*diff1 + diff2*diff2 + diff3*diff3

  IF (iter.gt.3.and.daux.lt.0.00001d0) GOTO 1

END DO

1 CONTINUE

end subroutine FindPoint

end subroutine MovePointInElement

end module particle_step
