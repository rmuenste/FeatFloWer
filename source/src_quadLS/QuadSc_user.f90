SUBROUTINE GetPresInitVal(X,Y,Z,Val)
REAL*8 X,Y,Z,Val(4)
REAL*8 :: dScale = (8.55/5.35)*0.008d0*0.45D0/(0.41D0**2d0)

! Val(1) = dScale*(2.5d0-X)
! Val(2) = -dScale
! Val(3) = 0d0
! Val(4) = 0d0!-1d1

Val(1) = 0d0!-1d1*X
Val(2) = 0d0!-1d1
Val(3) = 0d0
Val(4) = 0d0

RETURN

END SUBROUTINE GetPresInitVal   !InitCond
!------------------------------------------------
SUBROUTINE GetCylKnpr(X,Y,Z,iBndr,bCyl)
REAL*8 X,Y,Z
INTEGER iBndr
LOGICAL bCyl
REAL*8 :: PX = 0.5d0,PY = 0.2d0, RAD = 0.05d0

bCyl = .FALSE.
IF (iBndr.EQ.0) RETURN

DIST = SQRT((X-PX)**2d0+(Y-PY)**2d0) - RAD
IF (DIST.LT.1d-4) THEN
 bCyl = .TRUE.
END IF

END SUBROUTINE GetCylKnpr
!--------------------------------------------------
SUBROUTINE GetVeloInitVal(X,Y,Z,ValU,ValV,ValW)
USE var_QuadScalar, ONLY : bViscoElastic
REAL*8 X,Y,Z,ValU,ValV,ValW
REAL*8 PX,PY,dScale
INTEGER iC
REAL*8 :: RX = 0.0d0,RY = 0.0d0,RZ = 0.0d0, RAD = 0.245d0
REAL*8 :: R_inflow=4d0
REAL*8 :: PI=3.141592654d0

!  dScale = 16D0*0.45D0/(0.41D0**4) 
!  Val = dScale*Y*(0.41D0-Y)*Z*(0.41D0-Z)
ValU = 0d0
ValV = 0d0
if(bViscoElastic)then
  ValW = 1d0
else
  ValW = 0d0
end if

! if (Y.lt.0.5d0) then
!  ValU = tanh(3d1*(Y-0.25d0))
! else
!  ValU = tanh(3d1*(0.75d0-Y))
! end if
! 
! ValV = 0.05d0*sin(2d0*PI*X)

! dist = SQRT(x*x + y*y)
! dScale = 5d0/R_inflow*R_inflow
! IF (dist.lt.R_inflow) THEN
!  ValW= dScale*(dist+R_inflow)*(R_inflow-dist)
! END IF

RETURN

END SUBROUTINE GetVeloInitVal
!---------------------------------------------------
SUBROUTINE GetVeloBCVal(X,Y,Z,ValU,ValV,ValW,iT,t)
REAL*8 X,Y,Z,ValU,ValV,ValW,t
REAL*8 :: tt=4d0
INTEGER iT
REAL*8 :: RX = 0.0d0,RY = 0.0d0,RZ = 0.0d0, RAD1 = 0.25d0
REAL*8 :: RY1 = -0.123310811d0,RY2 = 0.123310811d0,Dist1,Dist2,R_In = 0.1d0
REAL*8  dScale
REAL*8 DIST
REAL*8 :: PI=3.141592654d0
REAL*8 :: R_inflow=4d0

ValU = 0d0
ValV = 0d0
ValW = 0d0

IF (iT.EQ.1) THEN
   dist = SQRT(z*z + y*y)
   dScale = 5.0d0/R_inflow*R_inflow
   IF (dist.lt.R_inflow) THEN
    ValU= dScale*(dist+R_inflow)*(R_inflow-dist)
   END IF
END IF

IF (iT.EQ.2) THEN
 dScale=0.2d0*(3d0/2d0)/(0.205d0*0.205d0)
 ValU=dScale*Y*(0.41d0-Y)
 ValV= 0d0
 ValW= 0d0
END IF

IF (iT.EQ.6) THEN
 dScale=(9d0/4d0)/(0.205d0*0.205d0*0.205d0*0.205d0)*sin(t*PI/8d0)
 ValU=dScale*Y*(0.41d0-Y)*Z*(0.41d0-Z)
 ValV= 0d0
 ValW= 0d0
END IF

IF (iT.EQ.8.OR.IT.EQ.9) THEN
 ValW = 1d0
END IF

RETURN

END SUBROUTINE GetVeloBCVal
!------------------------------------------------------
SUBROUTINE GetMixerKnpr(X,Y,Z,iBndr,inpr,D,t)
IMPLICIT NONE
REAL*8 X,Y,Z,t,d
REAL*8 :: PI = 6.283185307179586232
INTEGER :: inpr,iBndr
REAL*8 :: RX0=0d0,RY0=0d0
REAL*8 :: DA=0.268d0,DB=1.d0
REAL*8 :: RAD = 0.145d0,dAlpha,XT,YT,ZT,XB,YB,ZB
REAL*8 :: dBeta,XTT,YTT,ZTT,dist
REAL*8 :: mBottom=0.15d0,mTop=0.55d0,mThickness=0.2d0
REAL*8 dScale

inpr = 0
RETURN
!IF (iBndr.NE.0) RETURN
d=1d0

! First screw
XB = X
YB = Y-0.125
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = -t*PI
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

! Next the point needs to be transformed back to Z = 0 level
dBeta = -(ZT-mBottom)/(mThickness)*PI
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

IF (ZTT.LE.mTop.AND.ZTT.GE.mBottom) THEN
 DIST = SQRT((XTT-RX0)*(XTT-RX0)/DA+(YTT-RY0)*(YTT-RY0)/DB)
 d=dist-rad
 IF (DIST.LT.RAD) THEN
  inpr = 101
 END IF
END IF

! Second screw
XB = X
YB = Y+0.125
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = -t*PI+PI/4d0
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

! Next the point needs to be transformed back to Z = 0 level
dBeta = -(ZT-mBottom)/(mThickness)*PI
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

IF (ZTT.LE.mTop.AND.ZTT.GE.mBottom) THEN
 DIST = SQRT((XTT-RX0)*(XTT-RX0)/DA+(YTT-RY0)*(YTT-RY0)/DB)
 d=min(d,dist-rad)
 IF (DIST.LT.RAD) THEN
  inpr = 102
 END IF
END IF

RETURN

END SUBROUTINE GetMixerKnpr
!------------------------------------------------------------
SUBROUTINE GetVeloMixerVal(X,Y,Z,ValU,ValV,ValW,iP,t)
IMPLICIT NONE
REAL*8 :: PI=6.283185307179586232
INTEGER iP
REAL*8 X,Y,Z,ValU,ValV,ValW,t

IF (iP.EQ.101) THEN
 ValU =  -PI*(Y-0.125d0)
ELSE
 ValU =  -PI*(Y+0.125d0)
END IF
ValV = PI*X
ValW = 0d0


END SUBROUTINE GetVeloMixerVal


