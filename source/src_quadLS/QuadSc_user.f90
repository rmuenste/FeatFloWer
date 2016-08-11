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
ValW = 0d0

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
 ValU = -1d0*(120d0/60d0)*2d0*PI*Y
 ValV = +1d0*(120d0/60d0)*2d0*PI*X
END IF 

! IF (iT.EQ.1) THEN
!    dist = SQRT(x*x + y*y)
!    dScale = 5d0/R_inflow*R_inflow
!    IF (dist.lt.R_inflow) THEN
!     ValW= dScale*(dist+R_inflow)*(R_inflow-dist)
!    END IF
! END IF
! 
! IF (iT.EQ.6) THEN
!  dScale=(9d0/4d0)/(0.205d0*0.205d0*0.205d0*0.205d0)*sin(t*PI/8d0)
!  ValU=dScale*Y*(0.41d0-Y)*Z*(0.41d0-Z)
!  ValV= 0d0
!  ValW= 0d0
! END IF

RETURN

END SUBROUTINE GetVeloBCVal
!------------------------------------------------------
SUBROUTINE GetMixerKnpr(X,Y,Z,myCenterY,iBndr,inpr,D,t)
use QuadScalar, only: myFish
IMPLICIT NONE
REAL*8 X,Y,Z,t,d
REAL*8 :: PI = 4d0*ATAN(1d0)
REAL*8 :: XT,YT,ZT,XB,YB,ZB,dAlpha
REAL*8 :: myCenterY, myCenterRotY
REAL*8 :: dist1,dist2,t_i
INTEGER :: iBndr,inpr,iRot

inpr = 0

myCenterRotY = myCenterY+myFish%centerRotDisp(2)

dist1 = sqrt(x*x + y*y + (z-myCenterY)*(z-myCenterY)) - myFish%radius

IF (dist1.lt.0d0) THEN
 inpr = 1
 D=dist1
END IF

XB = X - myFish%centerRotDisp(1)
YB = Y - myCenterRotY
ZB = Z

! !
! iRot = INT(t/myFish%period)
! t_i  = t - iRot*myFish%period
! 
! IF (t_i.le.myFish%t_open) THEN
!  dAlpha =  0.601d0*(t_i/myFish%t_open)*PI + 0.025d0*PI
! ! ELSEIF(t_i.le.myFish%t_pause)then
! !  dAlpha =  0.626d0*PI
! !elseif(t_i.le.myFish%t_close)then
! else if(t_i.le.myFish%t_close)then
!  dAlpha = -0.601d0*((t_i-myFish%t_pause)/(myFish%t_close-myFish%t_pause))*PI+ 0.626d0*PI
! ELSE
!  dAlpha = 0.025d0*PI 
! End IF
! 
! XT = XB*cos(dAlpha) - YB*sin(dAlpha) + myFish%centerRotDisp(1)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha) + myCenterRotY
! ZT = ZB

! XT = XB
! YT = XB
! ZT = ZB
! 
! 
! dist2 = -min(xT-myFish%barBot,-xT+myFish%barTop,yT-myCenterY,-yT+(myFish%barLength + myFish%centerRotDisp(1) + myCenterY)) 
! IF (dist2.lt.0d0) THEN
!  IF (inpr.ne.101) inpr = 102
! END IF
! 
! D = MIN(dist1,dist2)

RETURN

END SUBROUTINE GetMixerKnpr
!------------------------------------------------------------
SUBROUTINE GetVeloMixerVal(X,Y,Z,ValU,ValV,ValW,iP,t)
use QuadScalar, only: myFish
IMPLICIT NONE
REAL*8 :: PI = 4d0*ATAN(1d0),t_i
INTEGER iP,iRot
REAL*8 X,Y,Z,ValU,ValV,ValW,t
REAL*8 :: myCenterRotY


return
ValU = 0d0
ValV = 0d0
ValW = 0d0

! write(*,*) iP
myCenterRotY = myFish%myCenterY + myFish%centerRotDisp(2)

IF(ip.ne.0)THEN
 ValW=myFish%myVelY
end if

! IF (iP.EQ.101) THEN
! 
! ValV=myFish%myVelY
! 
! ELSE
! 
!  iRot = INT(t/myFish%period)
!  t_i  = t - iRot*myFish%period
! 
!  IF (t_i.le.myFish%t_open) THEN
!   ValU = +0.21d0*PI*(Y-(myFish%myCenterY+myFish%centerRotDisp(2)))
!   ValV = -0.21d0*PI*(X-myFish%centerRotDisp(1)) + myFish%myVelY
! !  ELSEIF(t_i.le.myFish%t_pause)then
! !   ValU = 0.0d0
! !   ValV = myFish%myVelY
! ! elseif(t_i.le.myFish%t_close)then
!  else if(t_i.le.myFish%t_close)then 
!   !ValU = -0.9717d0*PI*(Y-(myFish%myCenterY+myFish%centerRotDisp(2)))
!   !ValV = +0.9717d0*PI*(X-myFish%centerRotDisp(1)) + myFish%myVelY
!   ValU = -1.0d0*PI*(Y-(myFish%myCenterY+myFish%centerRotDisp(2)))
!   ValV = +1.0d0*PI*(X-myFish%centerRotDisp(1)) + myFish%myVelY  
!   else
!    ValU = 0.0d0
!    ValV = myFish%myVelY  
!  END IF

!END IF

END SUBROUTINE GetVeloMixerVal


