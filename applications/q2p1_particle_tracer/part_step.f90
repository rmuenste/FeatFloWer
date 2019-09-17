module particle_step

USE var_QuadScalar, only : myParticleParam
USE geometry_processing, only : GetDistToSTL

REAL*8 dist_CGAL,cdx,cdy,cdz

real*8 DJ(8,3),q8
PARAMETER (Q8=0.125D0)
REAL*8 DJACI(3,3),DJAC(3,3),DBAS(27),DV0(3,27),DV1(3,27)
REAL*8 XI1,XI2,XI3,XX,YY,ZZ,diff1,diff2,diff3
REAL*8 :: dPI = 4d0*dATAN(1d0)
integer iter
logical :: ok_flag
REAL*8  DV_Loc(3),dFreq
REAL*8  P8(3,8), P(3), PR(3),P_New(3)
REAL*8 :: DeltaT


contains

!========================================================================================
!                            Sub: GetPointFromElement
!========================================================================================
SUBROUTINE GetPointFromElement(P8,P,iStart,BF,PR,iP)
LOGICAL BF

REAL*8  P8(3,8), P(3),PR(3)

DJ(1,1)=( P8(1,1)+P8(1,2)+P8(1,3)+P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(1,2)=( P8(2,1)+P8(2,2)+P8(2,3)+P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(1,3)=( P8(3,1)+P8(3,2)+P8(3,3)+P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(2,1)=(-P8(1,1)+P8(1,2)+P8(1,3)-P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(2,2)=(-P8(2,1)+P8(2,2)+P8(2,3)-P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(2,3)=(-P8(3,1)+P8(3,2)+P8(3,3)-P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(3,1)=(-P8(1,1)-P8(1,2)+P8(1,3)+P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(3,2)=(-P8(2,1)-P8(2,2)+P8(2,3)+P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(3,3)=(-P8(3,1)-P8(3,2)+P8(3,3)+P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(4,1)=(-P8(1,1)-P8(1,2)-P8(1,3)-P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(4,2)=(-P8(2,1)-P8(2,2)-P8(2,3)-P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(4,3)=(-P8(3,1)-P8(3,2)-P8(3,3)-P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(5,1)=( P8(1,1)-P8(1,2)+P8(1,3)-P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(5,2)=( P8(2,1)-P8(2,2)+P8(2,3)-P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(5,3)=( P8(3,1)-P8(3,2)+P8(3,3)-P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(6,1)=( P8(1,1)-P8(1,2)-P8(1,3)+P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(6,2)=( P8(2,1)-P8(2,2)-P8(2,3)+P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(6,3)=( P8(3,1)-P8(3,2)-P8(3,3)+P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(7,1)=( P8(1,1)+P8(1,2)-P8(1,3)-P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(7,2)=( P8(2,1)+P8(2,2)-P8(2,3)-P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(7,3)=( P8(3,1)+P8(3,2)-P8(3,3)-P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(8,1)=(-P8(1,1)+P8(1,2)-P8(1,3)+P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(8,2)=(-P8(2,1)+P8(2,2)-P8(2,3)+P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(8,3)=(-P8(3,1)+P8(3,2)-P8(3,3)+P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8

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

  DJAC(1,1)=DJ(2,1)+DJ(5,1)*XI2+DJ(6,1)*XI3+DJ(8,1)*XI2*XI3
  DJAC(1,2)=DJ(3,1)+DJ(5,1)*XI1+DJ(7,1)*XI3+DJ(8,1)*XI1*XI3
  DJAC(1,3)=DJ(4,1)+DJ(6,1)*XI1+DJ(7,1)*XI2+DJ(8,1)*XI1*XI2
  DJAC(2,1)=DJ(2,2)+DJ(5,2)*XI2+DJ(6,2)*XI3+DJ(8,2)*XI2*XI3
  DJAC(2,2)=DJ(3,2)+DJ(5,2)*XI1+DJ(7,2)*XI3+DJ(8,2)*XI1*XI3
  DJAC(2,3)=DJ(4,2)+DJ(6,2)*XI1+DJ(7,2)*XI2+DJ(8,2)*XI1*XI2
  DJAC(3,1)=DJ(2,3)+DJ(5,3)*XI2+DJ(6,3)*XI3+DJ(8,3)*XI2*XI3
  DJAC(3,2)=DJ(3,3)+DJ(5,3)*XI1+DJ(7,3)*XI3+DJ(8,3)*XI1*XI3
  DJAC(3,3)=DJ(4,3)+DJ(6,3)*XI1+DJ(7,3)*XI2+DJ(8,3)*XI1*XI2
!   DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
!        -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
!        +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))

  XX = DJ(1,1) + DJAC(1,1)*XI1 + DJ(3,1)*XI2 + DJ(4,1)*XI3 + DJ(7,1)*XI2*XI3
  YY = DJ(1,2) + DJ(2,2)*XI1 + DJAC(2,2)*XI2 + DJ(4,2)*XI3 + DJ(6,2)*XI1*XI3
  ZZ = DJ(1,3) + DJ(2,3)*XI1 + DJ(3,3)*XI2 + DJAC(3,3)*XI3 + DJ(5,3)*XI1*XI2

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
SUBROUTINE Extract_Particle(dcorvg,kvert,kedge,karea,kel_LdA,kel_ColA,&
                           VelU,VelV,VelW,nvt,net,nat,nel,tDelta,d_CorrDist)
USE PP3D_MPI, ONLY : myid,master,showid,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,RECVK_myMPI,SENDK_myMPI,MPI_COMM_WORLD

USE types, ONLY : myCompleteSet,myActiveSet,nActiveSet,myExchangeSet,myLostSet,&
                         nExchangeSet,nStartActiveSet,nLostSet

IMPLICIT NONE
REAL*8 dcorvg(3,*),VelU(*),VelV(*),VelW(*),point(3),tLevel,tDelta,d_CorrDist
INTEGER nkel,nvt,net,nat,nel
INTEGER kel_LdA(*),kel_ColA(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
INTEGER KDFG(27),KDFL(27)
INTEGER i,iPoint,jPoint,ivt,jel,iel,iFoundElem,pID,IERR
REAL*8 dist,pointR(3),LocVel(3,27)
LOGICAL :: bFound,bOut=.FALSE.
INTEGER :: nIter = 100, iIter, iLoc,iParticel,iActiveParticel
INTEGER iMonitor(40),iAux
REAL*8  dMonitor(40),dAux
INTEGER nXX,kk,iMon
PARAMETER (nXX = 40)
REAL*8 , ALLOCATABLE :: dSendVector(:),dRecvVector(:)
INTEGER, ALLOCATABLE :: iRecvVector(:)
REAL*8 cpx,cpy,cpz,cnormal(3)

ALLOCATE(dSendVector(nExchangeSet),iRecvVector(nExchangeSet),dRecvVector(nExchangeSet))
dSendVector = 1d30
iRecvVector = 0

! if (myid.eq.1) then
! DO iParticel = 1,nExchangeSet
!  write(*,*) myExchangeSet(iParticel)%coor(1:3)
! end do
! end if
! pause
! 

IF (myid.ne.0) THEN

DO iParticel = 1,nExchangeSet

iIter = 0

point  = myExchangeSet(iParticel)%coor
tLevel = myExchangeSet(iParticel)%time

  cdx = 1d1*point(1)
  cdy = 1d1*point(2)
  cdz = 1d1*point(3)
  
  CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
  if ((     myParticleParam%bRotationalMovement.and.dist_CGAL.lt. d_CorrDist*0.5d0).or. &
      (.not.myParticleParam%bRotationalMovement.and.dist_CGAL.gt.-d_CorrDist*0.5d0)) then
 ! if (dist_CGAL.lt.+d_CorrDist*0.5d0) then
!  if (dist_CGAL.gt.-d_CorrDist*0.5d0) then
!     write(*,*) 'Particel: ', iParticel 
!   write(*,*) 'before: ', P 

    call getclosestpointid(cdx,cdy,cdz,cpx,cpy,cpz,daux,0);
!   call projectonboundaryid(cdx,cdy,cdz,cpx,cpy,cpz,daux,0)
   
   cnormal = [cpx-cdx,cpy-cdy,cpz-cdz]
   daux = SQRT(cnormal(1)**2d0 + cnormal(2)**2d0 + cnormal(3)**2d0)
   cnormal = cnormal/daux

   if (myParticleParam%bRotationalMovement) then
    cdx = cpx - cnormal(1)*d_CorrDist*0.5d0!/dist_CGAL
    cdy = cpy - cnormal(2)*d_CorrDist*0.5d0!/dist_CGAL
    cdz = cpz - cnormal(3)*d_CorrDist*0.5d0!/dist_CGAL
   else
    cdx = cpx + cnormal(1)*d_CorrDist*0.5d0!/dist_CGAL
    cdy = cpy + cnormal(2)*d_CorrDist*0.5d0!/dist_CGAL
    cdz = cpz + cnormal(3)*d_CorrDist*0.5d0!/dist_CGAL
   end if
   
   P=0.1d0*[cdx,cdy,cdz]
!   write(*,*) 'after: ', P   
!   pause
   point = [P(1),P(2),P(3)]
  end if
  
!  point = [4.6d0, 0d0, 0d0]
  
  dist = (point(1)**2d0 + point(2)**2d0)**0.5d0
  dist_CGAL = myParticleParam%D_Out*0.5d0 - dist
  if (myParticleParam%bRotationalMovement.and.dist_CGAL.lt.+0.1d0*d_CorrDist*1d0) then
   daux = myParticleParam%D_Out*0.5d0 - 0.1d0*d_CorrDist*1.0d0
   cdx = point(1)*daux/dist
   cdy = point(2)*daux/dist
   cdz = point(3)
   if (myid.eq.0) write(*,'(A,3Es12.4,I0)') 'corrected:', myExchangeSet(iParticel)%coor, myExchangeSet(iParticel)%indice
   point = [cdx,cdy,cdz]
  end if
  
!  pause
  
  myExchangeSet(iParticel)%coor = point

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


! write(*,*) myid,point
! write(*,*) myid,dMonitor
! pause

DO iMon = 1,nXX

iPoint = iMonitor(iMon)

if (iPoint.ge.1.and.iPoint.le.nvt) then

DO iel = kel_LdA(iPoint),kel_LdA(iPoint+1)-1
 
 jel = kel_ColA(iel)
!DO iel = 1,nkel

! jel = kel(iel,iPoint)

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

!  write(*,*) '-->',myid,iMon,dist,nkel
!  pause
 
 CALL GetPointFromElement(P8,point,jPoint,bFound,pointR,iParticel)

 IF (bFound) THEN
  GOTO 1
 END IF

END DO

end if
! write(*,*) '-->',myid,iMon,iMon,dist,bFound
END DO

IF (.not.bFound) GOTO 2

1 CONTINUE

dSendVector(iParticel) = MAX(pointR(1),pointR(2),pointR(3))

2 CONTINUE

! write(*,*) '-->',myid,iMon,dist,bFound
! pause

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
  if (myid.eq.0) write(*,'(A,3Es12.4,I0)') 'lost:',myExchangeSet(iParticel)%coor, myExchangeSet(iParticel)%indice
 END IF
END DO

DEALLOCATE(dSendVector,iRecvVector,dRecvVector)

nActiveSet = nStartActiveSet+iActiveParticel
IF (bOut) WRITE(*,*) "My Active Set of Particles after exchange: ", nActiveSet

! write(*,*) myid,"My Active Set of Particles after exchange: ", nActiveSet,nkel
! pause

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
!                           Sub: Move_Particle
!========================================================================================
subroutine Move_Particle(dcorvg,kvert,kedge,karea,kel_LdA,kel_ColA,&
                         Vel0U,Vel0V,Vel0W,Vel1U,Vel1V,Vel1W,nvt,net,nat,nel,tDelta,tStart,d_CorrDist)
USE PP3D_MPI, ONLY : myid
USE types, ONLY : myActiveSet,myExchangeSet,nActiveSet,nExchangeSet,nStartActiveSet

IMPLICIT NONE
REAL*8 dcorvg(3,*),point(3),tLevel,tDelta,tStart,d_CorrDist
REAL*8 Vel0U(*),Vel0V(*),Vel0W(*),Vel1U(*),Vel1V(*),Vel1W(*)
INTEGER nvt,net,nat,nel
INTEGER kel_LdA(*),kel_ColA(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
INTEGER KDFG(27),KDFL(27)
INTEGER i,iPoint,jPoint,ivt,jel,iel,iFoundElem
REAL*8 dist,pointR(3),LocVel0(3,27),LocVel1(3,27)
LOGICAL :: bFound,bOut=.FALSE.
INTEGER :: nIter = 100, iIter, iLoc,iParticel,iLostParticel,iActiveParticel
INTEGER iMonitor(40),iAux
REAL*8  dMonitor(40),dAux
INTEGER nXX,kk,iMon
PARAMETER (nXX = 40)
REAL*8 cpx,cpy,cpz,cnormal(3)

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

if (iPoint.ge.1.and.iPoint.le.nvt) then

DO iel = kel_LdA(iPoint),kel_LdA(iPoint+1)-1
 
 jel = kel_ColA(iel)
! DO iel = 1,nkel
! 
!  jel = kel(iel,iPoint)

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
!    LocVel0(:,iLoc) = [1d-3,0d0,0d0]
!    LocVel1(:,iLoc) = [1d-3,0d0,0d0]
   LocVel0(:,iLoc) = [Vel0U(ivt),Vel0V(ivt),Vel0W(ivt)]
   LocVel1(:,iLoc) = [Vel1U(ivt),Vel1V(ivt),Vel1W(ivt)]
  END DO

  IF (bOut) WRITE(*,'(A,I0,3E14.4,A)') '------------   point ', myActiveSet(iParticel)%indice, point, '  ------------'

  DV0 = LocVel0
  DV1 = LocVel1
  P   = point 
  PR  = pointR
  
  CALL MovePointInElement(tLevel,tDelta,tStart,iParticel,jel)
  
  cdx = 1d1*P(1)
  cdy = 1d1*P(2)
  cdz = 1d1*P(3)
  
  CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
  if ((     myParticleParam%bRotationalMovement.and.dist_CGAL.lt. d_CorrDist*0.5d0).or. &
      (.not.myParticleParam%bRotationalMovement.and.dist_CGAL.gt.-d_CorrDist*0.5d0)) then
!   if (dist_CGAL.lt.d_CorrDist*0.5d0) then
!  if (dist_CGAL.gt.-d_CorrDist*0.5d0) then
!    write(*,*) 'Particel: ', iParticel 
!   write(*,*) 'before: ', P   
   call getclosestpointid(cdx,cdy,cdz,cpx,cpy,cpz,daux,0)
!   call projectonboundaryid(cdx,cdy,cdz,cpx,cpy,cpz,daux,0)
   
   cnormal = [cpx-cdx,cpy-cdy,cpz-cdz]
   daux = SQRT(cnormal(1)**2d0 + cnormal(2)**2d0 + cnormal(3)**2d0)
   cnormal = cnormal/daux
   
   if (myParticleParam%bRotationalMovement) then
    cdx = cpx - cnormal(1)*d_CorrDist*0.5d0!/dist_CGAL
    cdy = cpy - cnormal(2)*d_CorrDist*0.5d0!/dist_CGAL
    cdz = cpz - cnormal(3)*d_CorrDist*0.5d0!/dist_CGAL
   else
    cdx = cpx + cnormal(1)*d_CorrDist*0.5d0!/dist_CGAL
    cdy = cpy + cnormal(2)*d_CorrDist*0.5d0!/dist_CGAL
    cdz = cpz + cnormal(3)*d_CorrDist*0.5d0!/dist_CGAL
   end if
   
!    cdx = cpx - cnormal(1)*d_CorrDist*0.5d0!/dist_CGAL
!    cdy = cpy - cnormal(2)*d_CorrDist*0.5d0!/dist_CGAL
!    cdz = cpz - cnormal(3)*d_CorrDist*0.5d0!/dist_CGAL
   
   P=0.1d0*[cdx,cdy,cdz]
!   write(*,*) 'after: ', P   
!   pause
  end if
  
  dist = (P(1)**2d0 + P(2)**2d0)**0.5d0
  dist_CGAL = myParticleParam%D_Out*0.5d0 - dist
  if (myParticleParam%bRotationalMovement.and.dist_CGAL.lt.+d_CorrDist*1.0d0) then
   daux = myParticleParam%D_Out*0.5d0 - d_CorrDist*1.0d0
   cdx = P(1)*daux/dist
   cdy = P(2)*daux/dist
   cdz = P(3)
   P = [cdx,cdy,cdz]
  end if

  point = [P(1),P(2),P(3)]
!  point = [abs(P(1)),P(2),P(3)]
  
  iFoundElem = jel
  GOTO 1
 END IF

END DO

end if

END DO

IF (.not.bFound) THEN
!   cdx = 1d1*point(1)
!   cdy = 1d1*point(2)
!   cdz = 1d1*point(3)
! 
!   CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
  
!   write(*,*) point,dist_CGAL
!   if (dist_CGAL.gt.-d_CorrDist*0.5d0) then
! !    write(*,*) 'Particel: ', iParticel 
! !   write(*,*) 'before: ', P   
!    call projectonboundaryid(cdx,cdy,cdz,cpx,cpy,cpz,daux,0)
!    
!    cnormal = [cpx-cdx,cpy-cdy,cpz-cdz]
!    daux = SQRT(cnormal(1)**2d0 + cnormal(2)**2d0 + cnormal(3)**2d0)
!    cnormal = cnormal/daux
!    
!    cdx = cpx - cnormal(1)*d_CorrDist*0.5d0!/dist_CGAL
!    cdy = cpy - cnormal(2)*d_CorrDist*0.5d0!/dist_CGAL
!    cdz = cpz - cnormal(3)*d_CorrDist*0.5d0!/dist_CGAL
!    
!    point=0.1d0*[cdx,cdy,cdz]
!    
!    CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
!    call projectonboundaryid(cdx,cdy,cdz,cpx,cpy,cpz,daux,0)
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
SUBROUTINE MovePointInElement(tLevel,tEnd,tStart,iP,iE)
USE var_QuadScalar, ONLY : myParticleParam
USE PP3D_MPI, ONLY : myid
implicit none
LOGICAL :: BOUT=.FALSE.
REAL*8  tLevel,tEnd,tStart
REAL*8 :: dAlpha,XB,YB,ZB
REAL*8 :: RK_Velo(3,4),DETJ,dVelo,dElemSize,dFracTime,dist,dfrac
integer imove,ip,ie
logical :: ok_flag

! write(*,*) myid,'is trying ... ',tLevel,tEnd,tStart
! 
dFreq = (myParticleParam%f/60d0)

iMove = 0
dVelo = 0d0

XI1 = PR(1)
XI2 = PR(2)
XI3 = PR(3)

! WRITE(*,'(A,7E14.4)') 'my original position: ', P,tLevel
IF (bOut) WRITE(*,'(A,7E14.4)') 'my original position: ', XI1,XI2,XI3,P,tLevel

CALL RETURN_Velo()

DJ(1,1)=( P8(1,1)+P8(1,2)+P8(1,3)+P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(1,2)=( P8(2,1)+P8(2,2)+P8(2,3)+P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(1,3)=( P8(3,1)+P8(3,2)+P8(3,3)+P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(2,1)=(-P8(1,1)+P8(1,2)+P8(1,3)-P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(2,2)=(-P8(2,1)+P8(2,2)+P8(2,3)-P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(2,3)=(-P8(3,1)+P8(3,2)+P8(3,3)-P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(3,1)=(-P8(1,1)-P8(1,2)+P8(1,3)+P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(3,2)=(-P8(2,1)-P8(2,2)+P8(2,3)+P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(3,3)=(-P8(3,1)-P8(3,2)+P8(3,3)+P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(4,1)=(-P8(1,1)-P8(1,2)-P8(1,3)-P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(4,2)=(-P8(2,1)-P8(2,2)-P8(2,3)-P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(4,3)=(-P8(3,1)-P8(3,2)-P8(3,3)-P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(5,1)=( P8(1,1)-P8(1,2)+P8(1,3)-P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(5,2)=( P8(2,1)-P8(2,2)+P8(2,3)-P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(5,3)=( P8(3,1)-P8(3,2)+P8(3,3)-P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(6,1)=( P8(1,1)-P8(1,2)-P8(1,3)+P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(6,2)=( P8(2,1)-P8(2,2)-P8(2,3)+P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(6,3)=( P8(3,1)-P8(3,2)-P8(3,3)+P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(7,1)=( P8(1,1)+P8(1,2)-P8(1,3)-P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(7,2)=( P8(2,1)+P8(2,2)-P8(2,3)-P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(7,3)=( P8(3,1)+P8(3,2)-P8(3,3)-P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(8,1)=(-P8(1,1)+P8(1,2)-P8(1,3)+P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(8,2)=(-P8(2,1)+P8(2,2)-P8(2,3)+P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(8,3)=(-P8(3,1)+P8(3,2)-P8(3,3)+P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8

DJAC(1,1)=DJ(2,1)
DJAC(1,2)=DJ(3,1)
DJAC(1,3)=DJ(4,1)
DJAC(2,1)=DJ(2,2)
DJAC(2,2)=DJ(3,2)
DJAC(2,3)=DJ(4,2)
DJAC(3,1)=DJ(2,3)
DJAC(3,2)=DJ(3,3)
DJAC(3,3)=DJ(4,3)
DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))

dElemSize = (ABS(DETJ)**0.3333d0)

! write(*,*) myid,'Elemsize... ',dElemSize,DetJ
! !write(*,*) myid,'P8... ',P8
! pause

dVelo = SQRT(DV_Loc(1)*DV_Loc(1) + DV_Loc(2)*DV_Loc(2) + DV_Loc(3)*DV_Loc(3))
DeltaT = myParticleParam%hSize*dElemSize/dVelo
! DeltaT = 0.001d0
IF ((tLevel + DeltaT).GT.tEnd) THEN
 DeltaT = tEnd - tLevel
END IF
dFracTime = (tLevel + 0.5d0*DeltaT - tStart)/(tEnd - tStart)

55 CONTINUE

RK_Velo(:,1) = DV_Loc

! WRITE(*,'(A,7E14.4,I)') 'position0: ', P,P_New

CALL FindPoint(2d0/2d0,DV_loc)  ! --> XI1,XI2,XI3
! WRITE(*,'(A,7E14.4,I)') 'position1: ', P,P_New
CALL RETURN_Velo()              ! --> DV_loc
! WRITE(*,'(A,7E14.4,I)') 'position2: ', P,P_New
! RK_Velo(:,2) = DV_loc
!
! CALL FindPoint(1d0/2d0,DV_loc)  ! --> XI1,XI2,XI3
! CALL RETURN_Velo()              ! --> DV_loc
! RK_Velo(:,3) = DV_loc
!
! CALL FindPoint(2d0/2d0,DV_loc)  ! --> XI1,XI2,XI3
! CALL RETURN_Velo()              ! --> DV_loc
! RK_Velo(:,4) = DV_loc
!
! DV_Loc = (1d0/6d0)*(RK_Velo(:,1) + 2d0*RK_Velo(:,2) + 2d0*RK_Velo(:,3) + RK_Velo(:,4))
! ! DV_Loc = (0.25d0*RK_Velo(:,1) + 0.75d0*RK_Velo(:,2))
!
! CALL FindPoint(1d0/1d0,DV_loc)  ! --> XI1,XI2,XI3
P(1) = P_New(1)
P(2) = P_New(2)
P(3) = P_New(3)
! CALL RETURN_Velo()              ! --> DV_loc

tLevel = tLevel + DeltaT
iMove  = iMove + 1

! !WRITE(*,'(A,7E14.4,I)') 'position: ', P

IF ((ABS(XI1).GT.1.01d0.OR.ABS(XI2).GT.1.01d0.OR.ABS(XI3).GT.1.01d0).OR.(tLevel.ge.tEnd)) THEN
! Reached the end of the element -> leave the routine
GOTO 5
ELSE

! IF ((tLevel + DeltaT).GT.tEnd) THEN
!  DeltaT = tEnd - tLevel
! END IF
! We are still in the element - go back and continue marching in the element
! write(*,*) myid,'Elemsize... ',dElemSize,DetJ
! write(*,*) myid,'is trying to here ... ',iMove,tLevel,DeltaT,tEnd,myParticleParam%hSize,dElemSize,dVelo
! pause
GOTO 55

END IF

5 CONTINUE

! write(*,*) myid,'is trying to here ... ',tLevel,tEnd,tStart
! pause

XX = P(1)
YY = P(2)
ZZ = P(3)

dist = SQRT(YY*YY + XX*XX)


! WRITE(*,'(A,7E14.4,I)') 'my velocity: ', DV_loc,dVelo,DeltaT
! WRITE(*,'(A,4E14.4,I)') 'my final position: ', P,tLevel,iMove
! pause

!!!!!!!!!!!!!!!!!!!!! CORRECTION !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! CORRECTION !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! CORRECTION !!!!!!!!!!!!!!!!!!!!!!!!
! IF (dist.ge.myParticleParam%dEps1*myParticleParam%D_Out*0.5d0) THEN
! dFrac = dist/(myParticleParam%dEps2*myParticleParam%D_Out*0.5d0)
! XX = XX/dFrac
! YY = YY/dFrac
! 
! P(1) = XX
! P(2) = YY
! ! WRITE(*,*) '... problem is being fixed .. '
! END IF
!!!!!!!!!!!!!!!!!!!!! CORRECTION !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! CORRECTION !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! CORRECTION !!!!!!!!!!!!!!!!!!!!!!!!



end subroutine MovePointInElement

SUBROUTINE RETURN_Velo()
implicit none
!!!!!!!!!!!!!!!!!!!!!!
integer i

! Loop not neccesary because nothing depends on i?
DO i=1,27
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
END DO

DV_loc = 0d0
DO i=1,27
 DV_loc(1) = DV_loc(1) + DBAS(I)*DV0(1,I)
 DV_loc(2) = DV_loc(2) + DBAS(I)*DV0(2,I)
 DV_loc(3) = DV_loc(3) + DBAS(I)*DV0(3,I)
END DO

END SUBROUTINE RETURN_Velo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE FindPoint(fFact,dLocalVelo)
implicit none
REAL*8 fFact,dLocalVelo(3)
!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 :: DV_Rot(3),dVV(3)
REAL*8 :: dAlpha,XB,YB,ZB,dtheta
REAL*8 :: RK_Velo(3,4)
REAL*8 :: dRR,dU_r,dU_theta,dU_z,dRho,dRho0,dRho1,dZ,daux,dFact
!LOGICAL :: bRotationalMovement=.false.

if (myParticleParam%bRotationalMovement) then
 DV_Rot = [-2d0*dPI*dFreq*P(2), 2d0*dPI*dFreq*P(1),0d0]
 dtheta = datan(P(2)/P(1))
 IF (P(1).lt.0d0) THEN
  dtheta = dtheta + dPI
 END IF

 dRR    = DSQRT(P(1)**2d0 + P(2)**2d0)

 dVV  = dLocalVelo + DV_Rot

 dU_r     =  dVV(1)*dcos(dtheta) + dVV(2)*dsin(dtheta)
 dU_theta =(-dVV(1)*dsin(dtheta) + dVV(2)*dcos(dtheta))/dRR
 dU_z     =  dVV(3)

 dAlpha = fFact*DeltaT*dU_theta
 dRho   = fFact*DeltaT*dU_r
 dZ     = fFact*DeltaT*dU_z

 !!!! Angular transformation
 XB = P(1)
 YB = P(2)
 ZB = P(3)

 P_New(1) = XB*dcos(dAlpha) - YB*dsin(dAlpha)
 P_New(2) = XB*dsin(dAlpha) + YB*dcos(dAlpha)

 !!!! Radial transformation
 XB = P_New(1)
 YB = P_New(2)

 dRho0 = DSQRT(XB**2d0 + YB**2d0)
 dRho1 = dRho0 + dRho

 P_New(1) = XB*dRho1/dRho0
 P_New(2) = YB*dRho1/dRho0
 !!!! Axial transformation
 P_New(3) = ZB + dZ
ELSE

 P_New = P + fFact*DeltaT * dLocalVelo

end if

DO iter = 1,80

  DJAC(1,1)=DJ(2,1)+DJ(5,1)*XI2+DJ(6,1)*XI3+DJ(8,1)*XI2*XI3
  DJAC(1,2)=DJ(3,1)+DJ(5,1)*XI1+DJ(7,1)*XI3+DJ(8,1)*XI1*XI3
  DJAC(1,3)=DJ(4,1)+DJ(6,1)*XI1+DJ(7,1)*XI2+DJ(8,1)*XI1*XI2
  DJAC(2,1)=DJ(2,2)+DJ(5,2)*XI2+DJ(6,2)*XI3+DJ(8,2)*XI2*XI3
  DJAC(2,2)=DJ(3,2)+DJ(5,2)*XI1+DJ(7,2)*XI3+DJ(8,2)*XI1*XI3
  DJAC(2,3)=DJ(4,2)+DJ(6,2)*XI1+DJ(7,2)*XI2+DJ(8,2)*XI1*XI2
  DJAC(3,1)=DJ(2,3)+DJ(5,3)*XI2+DJ(6,3)*XI3+DJ(8,3)*XI2*XI3
  DJAC(3,2)=DJ(3,3)+DJ(5,3)*XI1+DJ(7,3)*XI3+DJ(8,3)*XI1*XI3
  DJAC(3,3)=DJ(4,3)+DJ(6,3)*XI1+DJ(7,3)*XI2+DJ(8,3)*XI1*XI2

  XX = DJ(1,1) + DJAC(1,1)*XI1 + DJ(3,1)*XI2 + DJ(4,1)*XI3 + DJ(7,1)*XI2*XI3
  YY = DJ(1,2) + DJ(2,2)*XI1 + DJAC(2,2)*XI2 + DJ(4,2)*XI3 + DJ(6,2)*XI1*XI3
  ZZ = DJ(1,3) + DJ(2,3)*XI1 + DJ(3,3)*XI2 + DJAC(3,3)*XI3 + DJ(5,3)*XI1*XI2

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

end module particle_step
