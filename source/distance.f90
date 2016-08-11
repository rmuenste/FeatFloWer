MODULE DTP
IMPLICIT NONE

TYPE RetVal
 Real*8 normal(3),distance
END TYPE RetVal
TYPE Cube
 Real*8 dCoor(3,8),dFunc(4),dMid(3)
 Real*8 dVals(8),dEdge(3,12),dPoly(3,12)
 INTEGER nEdge
 LOGICAL bEdge(12)
END TYPE Cube
TYPE(Cube) myCube
TYPE(RetVal) myRetVal
REAL*8 Point(3),ProjP(3),ProjPdist,SignedPdist
LOGICAL :: bInOut,bWrite=.FALSE.

CONTAINS

SUBROUTINE GetPolygon(dCoor,dFunc,dMid,dPoly,iN)
Real*8 dCoor(3,8),dFunc(4),dPoly(3,12),dMid(3)
INTEGER iN

bWrite = .FALSE.
myCube%dMid  = dMid
myCube%dCoor = dCoor
myCube%dFunc = dFunc

! Interpolate the distance function into the nodes !
CALL DTP_ComputeVertices
! Search for edges intersected by the interphase !
CALL DTP_GetInterphase
! Reconstruct the correctly numbered polygon  !
CALL DTP_GetPolygon

 dPoly = myCube%dPoly
 iN = myCube%nEdge
! END IF

END SUBROUTINE GetPolygon

SUBROUTINE DistanceToPolygon(dPoint,dGrad,dDist)
Real*8 dPoint(3),dDist,dGrad(3),ds

Point        = dPoint

CALL DTP_ProjectPoint
! Decide if the point is inside or outside of the polygon  !
CALL DTP_InOut
! Compute the shortest distance between the point to the polygon !
CALL DTP_ReturnValues

ds = SQRT(myRetVal%normal(1)**2d0+myRetVal%normal(2)**2d0+&
          myRetVal%normal(3)**2d0)

dDist = myRetVal%distance
dGrad = myRetVal%normal/ds

END SUBROUTINE DistanceToPolygon

SUBROUTINE DistanceToPlane(dCoor,dFunc,dMid,dPoint,dGrad,dDist,b)
LOGICAL b
Real*8 dCoor(3,8),dFunc(4),dPoint(3),dDist,dGrad(3),dMid(3),ds

bWrite = b
myCube%dMid  = dMid
myCube%dCoor = dCoor
myCube%dFunc = dFunc
Point        = dPoint

! Interpolate the distance function into the nodes !
CALL DTP_ComputeVertices
! Search for edges intersected by the interphase !
CALL DTP_GetInterphase
! Reconstruct the correctly numbered polygon  !
CALL DTP_GetPolygon
! Project the point to the plane of the polygon  !
CALL DTP_ProjectPoint
! Decide if the point is inside or outside of the polygon  !
CALL DTP_InOut
! Compute the shortest distance between the point to the polygon !
CALL DTP_ReturnValues

ds = SQRT(myRetVal%normal(1)**2d0+myRetVal%normal(2)**2d0+&
          myRetVal%normal(3)**2d0)

! IF (myRetVal%distance.LT.0d0) THEN
!  dDist = -myRetVal%distance
!  dGrad = -myRetVal%normal/ds
! ELSE
 dDist = myRetVal%distance
 dGrad = myRetVal%normal/ds
! END IF

END SUBROUTINE DistanceToPlane
!
! ------------------------------------------------------------------
!
SUBROUTINE DTP_ComputeVertices
INTEGER i
REAL*8 X,Y,Z,A,B,C,D,M1,M2,M3

IF (bWrite) Write(*,*)
IF (bWrite) Write(*,*) "Computation of nodal distance values ..."

M1 = myCube%dMid(1)
M2 = myCube%dMid(2)
M3 = myCube%dMid(3)
DO i=1,8
 X = myCube%dCoor(1,i)
 Y = myCube%dCoor(2,i)
 Z = myCube%dCoor(3,i)
 A = myCube%dFunc(1)
 B = myCube%dFunc(2)
 C = myCube%dFunc(3)
 D = myCube%dFunc(4)
 myCube%dVals(i) = A + B*(X-M1) + C*(Y-M2) + D*(Z-M3)
 IF (bWrite) Write(*,'(8G12.4)') i,myCube%dVals(i)
END DO

END SUBROUTINE DTP_ComputeVertices
!
! ------------------------------------------------------------------
!
SUBROUTINE DTP_GetInterphase
INTEGER :: iEdge(2,12)
INTEGER i,k,ivt1,ivt2
REAL*8 X1(3),X2(3),dPar
DATA iEdge/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/

IF (bWrite) Write(*,*)
IF (bWrite) Write(*,*) "Get the interphase points ..."

k=0
DO i=1,12
 ivt1 = iEdge(1,i)
 ivt2 = iEdge(2,i)
 myCube%bEdge(i) = .FALSE.
 IF (myCube%dVals(ivt1)*myCube%dVals(ivt2).LE.0d0) THEN
  IF (myCube%dVals(ivt1)*myCube%dVals(ivt2).EQ.0d0) THEN
   Write(*,*) "has to be generalized ... "
   STOP
  END IF
  k=k+1
  X1(1) = myCube%dCoor(1,ivt1)
  X1(2) = myCube%dCoor(2,ivt1)
  X1(3) = myCube%dCoor(3,ivt1)
  X2(1) = myCube%dCoor(1,ivt2)
  X2(2) = myCube%dCoor(2,ivt2)
  X2(3) = myCube%dCoor(3,ivt2)
  dPar  =-myCube%dVals(ivt1)/(myCube%dVals(ivt2)-myCube%dVals(ivt1))
  myCube%dEdge(:,i) = X1(:) + dPar*(X2(:)-X1(:))
  myCube%bEdge(i  ) = .TRUE.
  IF (bWrite) Write(*,'(8G12.4)') k,i,myCube%dEdge(:,i)
 END IF
END DO

myCube%nEdge = k

END SUBROUTINE DTP_GetInterphase
!
! ------------------------------------------------------------------
!
SUBROUTINE DTP_GetPolygon
INTEGER :: iEdge(2,12),iFace(4,6)
DATA iEdge/1,2,1,3,1,4,1,5,5,2,2,3,3,4,4,5,6,2,6,3,6,4,6,5/
DATA iFace/1,2,3,4,1,6,9,5,2,7,10,6,3,8,11,7,4,5,12,8,9,10,11,12/
INTEGER i,j,iE1,iE2,iA1,iA2,iE

IF (bWrite) Write(*,*)
IF (bWrite) Write(*,*) "Reconstruction of the polygon ..."

DO i=1,12
 IF (myCube%bEdge(i)) THEN
  iE2 = i
  iA2 = iEdge(1,iE2)
  myCube%dPoly(:,1) = myCube%dEdge(:,iE2)
  IF (bWrite) Write(*,'(8G12.4)') iE2,iA2,myCube%dEdge(:,iE2)
  EXIT
 END IF
END DO

DO i=2,myCube%nEdge
 iE1=iE2
 iA1=iA2
 DO j=1,4
  iE = iFace(j,iA1)
  IF(myCube%bEdge(iE).AND.iE.NE.iE1) THEN
   iE2=iE
   myCube%dPoly(:,i) = myCube%dEdge(:,iE2)
   IF (iA1.NE.iEdge(1,iE2)) iA2=iEdge(1,iE2)
   IF (iA1.NE.iEdge(2,iE2)) iA2=iEdge(2,iE2)
   IF (bWrite) Write(*,'(8G12.4)') iE2,iA2,myCube%dEdge(:,iE2)
   EXIT
  END IF
 END DO
END DO


! DO i=1,myCube%nEdge
!   IF (bWrite) Write(*,'(8G12.4)') i,myCube%dPoly(:,i)
! END DO

END SUBROUTINE DTP_GetPolygon
!
! ------------------------------------------------------------------
!
SUBROUTINE DTP_ProjectPoint
REAL*8 X,Y,Z,A,B,C,D,M1,M2,M3,daux

 IF (bWrite) Write(*,*)
 IF (bWrite) Write(*,*) "Projection ..."
 M1 = myCube%dMid(1)
 M2 = myCube%dMid(2)
 M3 = myCube%dMid(3)
 X = Point(1)
 Y = Point(2)
 Z = Point(3)
 A = myCube%dFunc(1)
 B = myCube%dFunc(2)
 C = myCube%dFunc(3)
 D = myCube%dFunc(4)

 daux = SQRT(B*B+C*C+D*D)
 ProjPdist = (A+B*(X-M1)+C*(Y-M2)+D*(Z-M3))/daux
 SignedPdist = SIGN(ProjPdist,1d0)
 ProjP(1) = Point(1) - ProjPdist*B
 ProjP(2) = Point(2) - ProjPdist*C
 ProjP(3) = Point(3) - ProjPdist*D
 IF (bWrite) Write(*,'(8G12.4)') ProjP,":",ProjPdist

END SUBROUTINE DTP_ProjectPoint
!
! ------------------------------------------------------------------
!
SUBROUTINE DTP_InOut
REAL*8 :: d2PI = 6.283185307179586476925287
REAL*8 dAngle,dCosTetha,m1,m2,m1m2,P1(3),P2(3)
INTEGER i,j

bInOut=.FALSE.
dAngle=0d0
DO i=1,myCube%nEdge
 P1(:) = myCube%dPoly(:,i  )-ProjP(:)
 j=i+1
 IF (i.eq.myCube%nEdge) j=1
 P2(:) = myCube%dPoly(:,j)-ProjP(:)

 m1 = SQRT(P1(1)*P1(1)+P1(2)*P1(2)+P1(3)*P1(3))
 m2 = SQRT(P2(1)*P2(1)+P2(2)*P2(2)+P2(3)*P2(3))

 m1m2 = m1*m2
 IF (m1m2.LE.1d-8) THEN
  bInOut=.TRUE.
  ProjPdist=0d0
  GOTO 1
 ELSE
  dCosTetha=(P1(1)*P2(1)+P1(2)*P2(2)+P1(3)*P2(3))/m1m2
  dAngle = dAngle + acos(dCosTetha)
 END IF

END DO

IF (ABS(ABS(dAngle)/d2PI-1d0).lt.1d-4) bInOut=.TRUE.

1 CONTINUE

IF (bInOut) THEN
 IF (bWrite) Write(*,*) "The point is inside of the polygon",dAngle!,ABS(dAngle/d2PI-1d0).lt.1d-4
ELSE
 IF (bWrite) Write(*,*) "The point is NOT inside of the polygon",dAngle!,ABS(dAngle/d2PI-1d0).lt.1d-4
END IF

END SUBROUTINE DTP_InOut
!
! ------------------------------------------------------------------
!
SUBROUTINE DTP_ReturnValues
REAL*8 dist, DistMin,P(3),P1(3),P3(3),m1,m2,P2(3),dd,ds
INTEGER i,iNB1,iNB2,iNB3

IF (bWrite) Write(*,*)
IF (bWrite) Write(*,*) "Closest point results..."
IF (bInOut) THEN
 IF (bWrite) Write(*,'(8G12.4)') "proj. point"
 myRetVal%normal(1) = myCube%dFunc(2)
 myRetVal%normal(2) = myCube%dFunc(3)
 myRetVal%normal(3) = myCube%dFunc(4)
 myRetVal%distance  = ProjPdist
ELSE
 DistMin=1d30
 DO i=1,myCube%nEdge
  P = myCube%dPoly(:,i)-ProjP(:)
  dist = SQRT(P(1)*P(1)+P(2)*P(2)+P(3)*P(3))
  IF (DistMin.GT.dist) THEN
   iNB2 = i
   DistMin = dist
!    IF (bWrite) Write(*,*) i, dist
  END IF
 END DO

 IF (iNB2.EQ.1) THEN
  iNB1=myCube%nEdge
 ELSE
  iNB1=iNB2-1
 END IF

 IF (iNB2.EQ.myCube%nEdge) THEN
  iNB3=1
 ELSE
  iNB3=iNB2+1
 END IF

 IF (bWrite) Write(*,'(8G12.4)') "segments: ",iNB1,iNB2," and ",iNB2,iNB3

 P2(:) = ProjP(:) - myCube%dPoly(:,iNB2)

 P1(:) = myCube%dPoly(:,iNB1) - myCube%dPoly(:,iNB2)
 m1 = P1(1)*P2(1)+P1(2)*P2(2)+P1(3)*P2(3)
 IF (m1.LE.0d0) GOTO 1
 m2 = P1(1)*P1(1)+P1(2)*P1(2)+P1(3)*P1(3)
 IF (m2.LE.m1) GOTO 1
 dd= m1/m2
 P2(:)=myCube%dPoly(:,iNB2)+dd*P1(:)
 IF (bWrite) Write(*,'(8G12.4)') "segment 1"
 GOTO 3

1 CONTINUE

 P1(:) = myCube%dPoly(:,iNB3) - myCube%dPoly(:,iNB2)
 m1 = P1(1)*P2(1)+P1(2)*P2(2)+P1(3)*P2(3)
 IF (m1.LE.0d0) GOTO 2
 m2 = P1(1)*P1(1)+P1(2)*P1(2)+P1(3)*P1(3)
 IF (m2.LE.m1) GOTO 2
 dd= m1/m2
 P2(:)=myCube%dPoly(:,iNB2)+dd*P1(:)
 IF (bWrite) Write(*,'(8G12.4)') "segment 2"
 GOTO 3

2 CONTINUE

 P2(:)=myCube%dPoly(:,iNB2)
 IF (bWrite) Write(*,'(8G12.4)') "point 2"

3 CONTINUE

 ProjP(:) = P2(:)
 P2(:) = Point(:)-ProjP(:)
 ds = P2(1)*myCube%dFunc(2)+P2(2)*myCube%dFunc(3)+&
      P2(3)*myCube%dFunc(4)
 ds = ds/ABS(ds)
 ProjPdist = ds*(SQRT(P2(1)*P2(1)+P2(2)*P2(2)+P2(3)*P2(3)))
 myRetVal%distance  = ProjPdist
 myRetVal%normal(1) = P2(1)/ProjPdist
 myRetVal%normal(2) = P2(2)/ProjPdist
 myRetVal%normal(3) = P2(3)/ProjPdist

!  P1(:) = myCube%dPoly(:,iNB1) - myCube%dPoly(:,iNB2)
!  P3(:) = myCube%dPoly(:,iNB3) - myCube%dPoly(:,iNB2)
!  P1(:) = myCube%dPoly(:,iNB2) + 1d-3*P1(:)
!  P3(:) = myCube%dPoly(:,iNB2) + 1d-3*P3(:)
!  P1(:) = P1(:)-ProjP(:)
!  P3(:) = P3(:)-ProjP(:)
! !  IF (bWrite) Write(*,*) "-----------------"
! !  IF (bWrite) Write(*,*) P1-P3
! !  IF (bWrite) Write(*,*) "-----------------"
!  m1 = SQRT(P1(1)*P1(1)+P1(2)*P1(2)+P1(3)*P1(3))
!  m2 = SQRT(P3(1)*P3(1)+P3(2)*P3(2)+P3(3)*P3(3))
!  IF (m1.LT.DistMin.OR.m2.lt.DistMin) THEN
!   IF (m1.LT.DistMin) THEN
!    IF (bWrite) Write(*,'(8G12.4)') "line 1 :", iNB1,iNB2,m1,m2,DistMin
!    
!   ELSE
!    IF (bWrite) Write(*,'(8G12.4)') "line 3 :", iNB2,iNB3,m1,m2,DistMin
!   END IF
!  ELSE
!   IF (bWrite) Write(*,'(8G12.4)') "point 2 :", iNB2
!  END IF

END IF

IF (bWrite) Write(*,'(8G12.4)') myRetVal%normal,":",myRetVal%distance
END SUBROUTINE DTP_ReturnValues
!
! ------------------------------------------------------------------
!
END MODULE DTP

