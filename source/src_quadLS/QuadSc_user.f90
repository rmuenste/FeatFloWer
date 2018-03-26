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
USE var_QuadScalar, ONLY : bViscoElastic,bViscoElasticFAC
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

if(bViscoElasticFAC)then
  ValV = 0d0
  ValW = 0d0
  dScale=3d0/2d0
  ValU=dScale*(1d0-0.25d0*Y*Y)
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
use var_QuadScalar, only : myFBM
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess
USE PP3D_MPI, ONLY:myid
implicit none

REAL*8 X,Y,Z,ValU,ValV,ValW,t
REAL*8 :: tt=4d0
INTEGER iT
REAL*8 :: RX = 0.0d0,RY = 0.0d0,RZ = 0.0d0, RAD1 = 0.25d0
REAL*8 :: RY1 = -0.123310811d0,RY2 = 0.123310811d0,Dist1,Dist2,R_In = 0.1d0
REAL*8  dScale
REAL*8 dInnerRadius,dOuterRadius,dVolFlow,daux
REAL*8 DIST
REAL*8 :: PI=dATAN(1d0)*4d0
REAL*8 :: R_inflow=4d0,dNx,dNy,dNz,dNorm


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

IF (iT.EQ.4.OR.iT.EQ.5) THEN
 dScale=3d0/2d0
 ValU=dScale*(1d0-0.25d0*Y*Y)
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

if(it.eq.10)then

  dInnerRadius = myProcess%MinInflowDiameter*0.5d0
  dOuterRadius = myProcess%MaxInflowDiameter*0.5d0

  dVolFlow = (1e3/3.6d3)*myProcess%Massestrom/(myThermodyn%density) ! mm3/s
  daux = (PI/6d0)*(dInnerRadius+dOuterRadius)*((dOuterRadius-dInnerRadius)**3d0)
  dScale = (dVolFlow/1d0)/daux

  DIST = SQRT(X**2d0+Y**2d0)

  IF (DIST.GT.dInnerRadius.AND.DIST.LT.dOuterRadius) THEN
   ValW= dScale*(DIST-dInnerRadius)*(dOuterRadius-DIST)
  END IF
end if


if(it.eq.11)then
  ValU =  -DBLE(myProcess%ind)*PI*Y*(myProcess%Umdr/3d1)
  ValV =   DBLE(myProcess%ind)*PI*X*(myProcess%Umdr/3d1)
  ValW =   0d0
end if

!! InitBoundaryStructure: changed from I1 to I2 !!
! used for pipe/pipesphere.prj
IF (iT.EQ.25) THEN
 ValU= 1d0
 ValV= 0d0
 ValW= 0d0
END IF
! used for FAC3Ds
IF (it.EQ.26) THEN
 dScale=0.2d0*(9d0/4d0)/(0.205d0*0.205d0*0.205d0*0.205d0)
 ValU=dScale*(0.41d0-Y)*Y*(0.41d0-Z)*Z
 ValV=0d0
 ValW=0d0
END IF
IF (it.EQ.27) THEN
 dScale=(9d0/4d0)/(0.205d0*0.205d0*0.205d0*0.205d0)
 ValU=dScale*(0.41d0-Y)*Y*(0.41d0-Z)*Z
 ValV=0d0
 ValW=0d0
END IF

IF (iT.EQ.21) THEN
 R_inflow = 0.55d0
 dist = SQRT((x-0.47)**2d0 + (y+0.60)**2d0 + (z-3.15)**2d0)
 dScale = 0.1d0/R_inflow*R_inflow
 IF (dist.lt.R_inflow) THEN
  ValW= -dScale*(dist+R_inflow)*(R_inflow-dist)
 END IF
END IF

IF (iT.EQ.22) THEN
 R_inflow = 0.3d0
 dist = SQRT((x-5.5)**2d0 + (y-0.01)**2d0 + (z+0.84)**2d0)
 IF (dist.lt.R_inflow) THEN
  dScale = 5.0d0/R_inflow*R_inflow
  ValU= -0.5d0*dScale*(dist+R_inflow)*(R_inflow-dist)
  ValV=  0d0
  ValU= -0.5d0*dScale*(dist+R_inflow)*(R_inflow-dist)
 END IF
END IF

!Aneurysm Case 1
IF (iT.EQ.33) THEN
 R_inflow = 1.0d0
 dist = SQRT((x+6.0)**2d0 + (y-12.9)**2d0 + (z-5.2)**2d0)
 IF (dist.lt.R_inflow) THEN
  dScale = 1.0d0/R_inflow*R_inflow
  dScale = dScale*(dist+R_inflow)*(R_inflow-dist)
  ValU=  1.0084836*dScale
  ValV= -2.0737041*dScale
  ValW= -1.0046438*dScale
 END IF
END IF

!Aneurysm Case 2
IF (iT.EQ.34) THEN
 dNx = 3.0218865
 dNY = 1.3312501
 dNz =-2.3385094
 dNorm = SQRT(dNX*dNx +dNy*dNy + dNz*dNz)
 dNx = dNx/dNorm
 dNY = dNy/dNorm
 dNz = dNz/dNorm
 dist = SQRT((X-(-4.72736))*(X-(-4.72736)) + (Y-(-1.301935))*(Y-(-1.301935)) + (Z-(2.571679))*(Z-(2.571679)))
!  dVectMag = 2d0*my_InFlow%VectMag*(1.12d0-dist*dist)
 dScale = 5d0*(2d0*max(0d0,(1.12d0-dist*dist)))
!  WRITE(*,*) dist,dVectMag
 ValU = dNx*dScale
 ValV = dNy*dScale
 ValW = dNz*dScale
END IF
 
IF (iT.EQ.41) THEN
 ValW=RotParabolicVelo(0d0,0d0,100d0,-1d0,0.495d0)
END IF

IF (iT.EQ.45) THEN
 ValW=RotParabolicVelo(0d0,0d0,1449d0,1d0,2.495d0)
END IF

  
IF (iT.EQ.51) THEN
  ValW=RotParabolicVelo(0d0,0d0,1270d0,1d0,2.495d0)
!   ValW=RotParabolicVelo(0d0,0d0,2127d0,1d0,2.495d0) !
END IF
IF (iT.EQ.52) THEN
  ValW=RotParabolicVelo(0d0,6d0,56d0,1d0,1.245d0)
!  ValW=RotParabolicVelo(0d0,6d0,54d0,1d0,1.245d0)
END IF
IF (iT.EQ.53) THEN
  ValW=RotParabolicVelo(0d0,-6d0,56d0,1d0,1.245d0)
!  ValW=RotParabolicVelo(0d0,-6d0,54d0,1d0,1.245d0)
END IF
IF (iT.EQ.54) THEN
  ValW=RotParabolicVelo(0d0,6d0,67d0,1d0,1.245d0)
!   ValW=RotParabolicVelo(0d0,6d0,59d0,1d0,1.245d0)
END IF

! centroplast
IF (iT.EQ.61) THEN
  ValW=RotParabolicVelo(0d0,0d0,25d0,1d0,1.245d0)
!   ValW=RotParabolicVelo(0d0,6d0,59d0,1d0,1.245d0)
END IF

! egeplast
IF (iT.EQ.62) THEN
  ValW=RotParabolicVelo(2.4d0,-4.155d0,20d0,1d0,0.39d0)
!   ValW=RotParabolicVelo(0d0,6d0,59d0,1d0,1.245d0)
END IF

IF (iT.EQ.99) THEN
 ValW = -myFBM%ParticleNew(1)%Velocity(3)
END IF

RETURN

 CONTAINS
 REAL*8 function RotParabolicVelo(XR,YR,DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux,YR,XR
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO) ! mm3/s
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux
  DIST = SQRT((X-XR)**2d0+(Y-YR)**2d0)
 IF (DIST.LT.dR) THEN
  RotParabolicVelo = dScale*(DR - DIST)*(DR + DIST)
 ELSE
  RotParabolicVelo = 0d0
 END IF
 END 

END SUBROUTINE GetVeloBCVal
!------------------------------------------------------
SUBROUTINE GetMixerKnpr(X,Y,Z,iBndr,inpr,D,t)
IMPLICIT NONE
REAL*8 X,Y,Z,t,d
REAL*8 :: PI = 6.283185307179586232
INTEGER :: inpr,iBndr
REAL*8 :: XT,YT,ZT,XB,YB,ZB
REAL*8 dScale

integer ipc,isin
real*8 d_temp

inpr = 0
ipc = 0
isin = 0

xt = 1d1*x
yt = 1d1*y
zt = 1d1*z

call isinelementid(xt,yt,zt,ipc,isin)
call getclosestpointid(xt,yt,zt,xb,yb,zb,d_temp,ipc)

if (isin.gt.0) then
 D = +abs(1d-1*d_temp)
else
 inpr = 103
 D = -abs(1d-1*d_temp)
end if

return

END SUBROUTINE GetMixerKnpr
!------------------------------------------------------------
SUBROUTINE GetVeloMixerVal(X,Y,Z,ValU,ValV,ValW,iP,t)
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess
IMPLICIT NONE
REAL*8 :: myPI
INTEGER iP
REAL*8 X,Y,Z,ValU,ValV,ValW,t

myPI = dATAN(1d0)*4d0

!  write(*,*)'vahhal'

SELECT CASE(iP)
 CASE (100) ! Y positiv
  ValU = 0d0
  ValV = 0d0
  ValW = 0d0
 CASE (101) ! Y positiv
  ValU =  -DBLE(myProcess%ind)*myPI*(Y-mySigma%a/2d0)*(myProcess%Umdr/3d1)
  ValV =   DBLE(myProcess%ind)*myPI*X*(myProcess%Umdr/3d1)
  ValW = 0d0
 CASE (102) ! Y negativ
  ValU =  -DBLE(myProcess%ind)*myPI*(Y+mySigma%a/2d0)*(myProcess%Umdr/3d1)
  ValV =   DBLE(myProcess%ind)*myPI*X*(myProcess%Umdr/3d1)
  ValW = 0d0
 CASE (103) ! Y negativ
  ValU =  -DBLE(myProcess%ind)*myPI*Y*(myProcess%Umdr/3d1)
  ValV =   DBLE(myProcess%ind)*myPI*X*(myProcess%Umdr/3d1)
!   write(*,*)'case 103 valu, valv:',valu,valv
  ValW =   0d0
 CASE (104) ! Y negativ
  ValU =  0d0
  ValV =  0d0
  ValW =  00d0/60d0
END SELECT


END SUBROUTINE GetVeloMixerVal
!
!------------------------------------------------------------
!
FUNCTION ViscosityModel(NormShearSquare)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myRheology
IMPLICIT NONE

real*8 :: ViscosityModel
real*8, intent (in) :: NormShearSquare
REAL*8 :: dStrs, aT = 1d0,dLimStrs
REAL*8 :: VNN
REAL*8 :: dN
! REAL*8 :: VCN,VNN,frac

dStrs = NormShearSquare**0.5d0
dLimStrs = MIN(1d5,MAX(1d-2,ABS(dStrs)))

! Paderborn Carreau
IF (myRheology%Equation.EQ.1) THEN
 VNN = (1d1*myRheology%A)*aT*(1d0+myRheology%B*aT*dLimStrs)**(-myRheology%C)
END IF

!PowerLaw
IF (myRheology%Equation.EQ.2) THEN
 VNN =(1d1*myRheology%K)*aT*(dLimStrs)**(-(1d0-myRheology%n))
END IF

!POLYFLOW carreau
IF (myRheology%Equation.EQ.3) THEN
 VNN = (1d1*myRheology%A)*(1d0+((myRheology%B*dLimStrs)**2d0))**(0.5d0*(myRheology%C-1d0)) 
END IF

!Ellis
IF (myRheology%Equation.EQ.4) THEN
 VNN = (1d1*myRheology%A)/(1d0+(dLimStrs/myRheology%B)**myRheology%C) 
END IF

! HogenPowerLaw
IF (myRheology%Equation.EQ.5) THEN
 dN = Properties%PowerLawExp-1d0
 VNN = Properties%Viscosity(1)*(1d-4 + NormShearSquare)**dN
END IF


dLimStrs = MIN(myRheology%ViscoMax,MAX(myRheology%ViscoMin,ABS(dStrs)))
! WRITE(*,*) dLimStrs,myRheology%eta_max,myRheology%eta_min
ViscosityModel = VNN


RETURN

end function ViscosityModel
