SUBROUTINE GetInitVal(X,Y,Z,Val)
IMPLICIT NONE
REAL*8 X,Y,Z,Val(4)
REAL*8 :: RY1 = -0.05d0,RY2 = 0.05d0
REAL*8 :: dx,dy,dz,dr,dF,dA=1d0,dB=1d0,dC=2d0

REAL*8 :: RX=0.0D0,RY=0.0D0,RZ= 0.020D0
! REAL*8 :: RX=0.0D0,RY=0.0D0,RZ=-0.29D0
REAL*8 :: RAD

! IF (Z.GT.RZ) THEN
!  dx = X-RX+0.01
!  dy = Y-RY+0.01
!  dz = Z-(0.08d0)
!  RAD = 0.040d0
! ELSE
 dx = X-0d0
 dy = Y-0d0
 dz = Z-(0.04d0)
 RAD = 0.060d0
! END IF

! REAL*8 :: RAD=0.0450d0
! 
! IF (X-Y.GT.0d0) THEN
!  dx = X-0.050d0
!  dy = Y+0.050d0
!  dz = Z+0.050d0
! ELSE
!  dx = X+0.050d0
!  dy = Y-0.050d0
!  dz = Z+0.050d0
! END IF
! IF (Z.LT.-0.3) THEN
!  dx = X-RX
!  dy = Y-RY
!  dz = Z-(-0.42)
! END IF

! IF (Z.GT.RZ) THEN
!  dx = X-RX
!  dy = Y-RY
!  dz = Z-(-0.258)
! ELSE
!  dx = X-RX
!  dy = Y-RY
!  dz = Z-(-0.32)
! END IF

dr = SQRT(dx*dx/dA + dy*dy/dB+ dz*dz/dC)
Val(1) = dr - RAD
Val(2) = dx/(dA*dr)
Val(3) = dy/(dB*dr)
Val(4) = dz/(dC*dr)

RETURN

! REAL*8 :: RX=0.0D0,RY=0.0D0,RZ1=0.24D0,RZ2=0.301
! REAL*8 :: RAD=0.015d0
! 
! Z = Z + 0.045
! IF (Z.GT.RZ1+0.002) THEN
!  IF (Z.GT.RZ2) THEN
!   dx = (X-RX)
!   dy = (Y-RY)
!   dr = SQRT(dx**2d0 + dy**2d0)
!   Val(1) = dr - RAD
!   Val(2) = dx/dr
!   Val(3) = dy/dr
!   Val(4) = 0d0
! ELSE
!   dF = 1d0-0.333333d0*(0.3-Z)/(0.3-RZ1)
!   dx = dF*(X-RX)
!   dy = dF*(Y-RY)
!   dz = 0d0
!   dr = SQRT(dx**2d0 + dy**2d0)
!   Val(1) = dr - RAD
!   Val(2) = dx/dr
!   Val(3) = dy/dr
!   Val(4) = (0.333333d0/(0.3-RZ1))*SQRT(dx**2d0 + dy**2d0)
!  END IF
! ELSE
!  dF = 1d0-0.333333d0
!  dx = (X-RX)
!  dy = (Y-RY)
!  dz = (Z-RZ1)
!  dr = SQRT(dx**2d0 + dy**2d0 + dz**2d0)
!  Val(1) = dr - 0.0218d0
!  Val(2) = dx/dr
!  Val(3) = dy/dr
!  Val(4) = dz/dr
! 
! END IF


RETURN

! REAL*8 :: RX=0.0D0,RY=0.0D0,RZ=0.3D0
! REAL*8 :: RAD=0.045d0
! dx = X-RX
! dy = Y-RY
! dz = Z-RZ
! 
! dr = SQRT(dx*dx + dy*dy+ dz*dz)
! Val(1) = dr - RAD
! Val(2) = dx/dr
! Val(3) = dy/dr
! Val(4) = dz/dr
! 
! RETURN

! REAL*8 :: RX=0.0D0,RY=0.0D0,RZ=0.05D0
! REAL*8 :: RAD=0.05d0
! IF (Y.GE.RY1.AND.Y.LE.RY2) THEN
!  dx = X-RX
!  dy = 0d0
!  dz = Z-RZ
! ELSE
!  IF (Y.LT.RY1) THEN
!   dx = X-RX
!   dy = Y-RY1
!   dz = Z-RZ
!  ELSE
!   dx = X-RX
!   dy = Y-RY2
!   dz = Z-RZ
!  END IF
! END IF
! 
! dr = SQRT(dx*dx + dy*dy+ dz*dz)
! Val(1) = dr - RAD
! Val(2) = dx/dr
! Val(3) = dy/dr
! Val(4) = dz/dr
! 
! RETURN

END
!
! ----------------------------------------------
!
FUNCTION GetBCVal(X,Y,Z,t) RESULT(Val)
IMPLICIT NONE
REAL*8 X,Y,Z,Val,t,dr
REAL*8 PXt,PYt,PZt
REAL*8 :: PX=0.25D0,PY=0.25D0,PZ=0.25D0,RAD=0.015d0

dr = SQRT(X*X + Y*Y)
Val = dr - RAD

! PXt = PX + r*COS(1d0*t)
! PYt = PY + r*SIN(1d0*t)
! PZt = PZ
! 
! dr = SQRT((X-PXt)*(X-PXt)+(Y-PYt)*(Y-PYt)+(Z-PZt)*(Z-PZt))
! Val = (dr - RAD)

RETURN

END
!
! ----------------------------------------------
!
SUBROUTINE RI_BoundCond(Val,iVal,nel)
REAL*8 Val(*)
INTEGER iVAl(*),nel,ip

DO iel=1,nel
 IF (iVal(iel).EQ.0) THEN
  ip=4*(iel-1)
  Val(ip+1) = 0d0
  Val(ip+2) = 0d0
  Val(ip+3) = 0d0
  Val(ip+4) = 0d0
 END IF
END DO

END SUBROUTINE RI_BoundCond

