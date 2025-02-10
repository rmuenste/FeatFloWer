SUBROUTINE GetPresInitVal(X,Y,Z,Val)
USE Sigma_User, ONLY : myProcess,mySigma
REAL*8 X,Y,Z,Val(4)
REAL*8 :: dScale = (8.55/5.35)*0.008d0*0.45D0/(0.41D0**2d0)

Val(1) = 0d0
Val(2) = 0d0
Val(3) = 0d0
Val(4) = 0d0

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".OR.&
    ADJUSTL(TRIM(mySigma%cType)).EQ."DIE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."NETZSCH") THEN
 IF (ADJUSTL(TRIM(myProcess%pTYPE)).EQ."PRESSUREDROP") THEN
  Val(1) = myProcess%dPress/mySigma%L*Z
  Val(2) = 0d0!-1d1
  Val(3) = 0d0
  Val(4) = myProcess%dPress/mySigma%L
 END IF
END IF

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
USE var_QuadScalar, ONLY : b2DViscoBench,b3DViscoBench, referenceVelocity
use fbm, only: fbm_updateFBM, fbm_velBCTest, fbm_velValue
REAL*8 X,Y,Z,ValU,ValV,ValW
REAL*8 PX,PY,dScale
INTEGER iC
REAL*8 :: RX = 0.0d0,RY = 0.0d0,RZ = 0.0d0, RAD = 0.245d0
REAL*8 :: R_inflow=4d0
REAL*8 :: PI=3.141592654d0

ValU = 0d0
ValV = 0d0
ValW = 0d0
return

!ValU = 2.0d0 * Y
! dScale = (1d0) * (3d0/2d0)/(0.5d0**2)
! ValU = dScale * (y**2 - 0.25d0)
 dScale = (1d0) * (3d0/2d0)/(0.5d0**2)
 ValU = dScale * (0.5d0-Y) * (0.5d0+Y) - referenceVelocity
 ValV = 0d0
 ValW = 0d0


if(b2DViscoBench)then
  ValV = 0d0
  ValW = 0d0
  dScale=3d0/2d0
  ValU=dScale*(1d0-0.25d0*Y*Y)
end if

if(b3DViscoBench)then
  ValU = 0d0
  ValV = 0d0
  ValW = 1d0
end if

! ValU = 0d0
! ValV = 0d0
! ValW = 0d0

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
use var_QuadScalar, only : myFBM, referenceVelocity
use fbm, only: fbm_updateFBM, fbm_velBCTest, fbm_velValue
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,myMultiMat
USE PP3D_MPI, ONLY:myid
implicit none

REAL*8 X,Y,Z,ValU,ValV,ValW,t
REAL*8 :: tt=4d0
INTEGER iT,j,iInflow,iMat
REAL*8 :: RX = 0.0d0,RY = 0.0d0,RZ = 0.0d0, RAD1 = 0.25d0
REAL*8 :: RY1 = -0.123310811d0,RY2 = 0.123310811d0,Dist1,Dist2,R_In = 0.1d0
REAL*8  dScale,XX,YY,ZZ
REAL*8 dInnerRadius,dOuterRadius,dMassFlow,dVolFlow,daux,dInnerInflowRadius,dDensity,dDensitySlope,&
       dTemperature
REAL*8 DIST
REAL*8 :: PI=dATAN(1d0)*4d0,myTwoPI=2d0*dATAN(1d0)*4d0
REAL*8 :: R_inflow=4d0,dNx,dNy,dNz,dNorm,dCenter(3),dNormal(3),dProfil(3)
REAL*8 :: dMidpointA(3), dMidpointB(3)
REAL*8 :: U_bar, h, normalizedTime, val,dFact
real*8, dimension(11) :: x_arr, y_arr, CC, DD, MM
integer iSubInflow
real*8 :: dRPM 

ValU = 0d0
ValV = 0d0
ValW = 0d0

IF (iT.EQ.23) THEN
 ValU= 1.0d0
 ValV= 0d0
 ValW= 0d0
END IF

IF (iT.EQ.202) THEN
 dScale = (1d0) * (3d0/2d0)/(0.5d0**2)
 !ValU = dScale * (y**2 - 0.25d0)
 ValU = dScale * (0.5d0-Y) * (0.5d0+Y) - referenceVelocity 
 ValV = 0d0
 ValW = 0d0
END IF

IF (iT.EQ.203) THEN
 ValU = 0.0d0 - referenceVelocity
 ValV = 0d0
 ValW = 0d0
END IF

! This is actually 771
! But we hack it to be 773
!IF (iT.EQ.771) THEN
IF (iT.EQ.773) THEN
  dRPM = 4d0
  ValU =  -myTwoPI*Y*(dRPM/6d1)
  ValV =   myTwoPI*X*(dRPM/6d1)
  ! one rotation takes 1min=60s ==> in one roatation the translation is 0.193*4=0.772cm ==> translation velocity is 0.772cm/min = 0.772cm/60s
  ValW =   -0.77d0*(dRPM/60d0)
END IF
 
IF (iT.EQ.772) THEN
  dRPM = 8d0
  ValU =  -myTwoPI*Y*(dRPM/6d1)
  ValV =   myTwoPI*X*(dRPM/6d1)
  ! one rotation takes 1min=60s ==> in one roatation the translation is 0.193*4=0.772cm ==> translation velocity is 0.772cm/min = 0.772cm/60s
  ValW =   -0.77d0*(dRPM/60d0)
END IF
 
! This is actually 773
! But we hack it to be faster
!IF (iT.EQ.773) THEN
IF (iT.EQ.771) THEN
  !dRPM = 12d0
  dRPM = 40d0
  ValU =  -myTwoPI*Y*(dRPM/6d1)
  ValV =   myTwoPI*X*(dRPM/6d1)
  ! one rotation takes 1min=60s ==> in one roatation the translation is 0.193*4=0.772cm ==> translation velocity is 0.772cm/min = 0.772cm/60s
  ValW =   -0.77d0*(dRPM/60d0)
END IF


RETURN

 CONTAINS
include '../include/ProfileFunctions.f90'

END SUBROUTINE GetVeloBCVal
!------------------------------------------------------------
SUBROUTINE GetVeloMixerVal(X,Y,Z,ValU,ValV,ValW,iP,t)
use, intrinsic :: ieee_arithmetic
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess
IMPLICIT NONE
REAL*8 :: myTwoPI=2d0*dATAN(1d0)*4d0
INTEGER iP
REAL*8 X,Y,Z,ValU,ValV,ValW,t,yb,xb,zb,dYShift,dYShiftdt,timeLevel
integer iScrewType,iSeg

IF (iP.lt.200) then

 SELECT CASE(iP)
  CASE (100) ! Y positiv
   ValU = 0d0
   ValV = 0d0
   ValW = 0d0
  CASE (101) ! Y positiv
   IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
    ValU =  -DBLE(myProcess%iInd*myProcess%ind)*myTwoPI*(Y-mySigma%a/2d0)*(myProcess%Umdr/6d1)
    ValV =   DBLE(myProcess%iInd*myProcess%ind)*myTwoPI*(X)*(myProcess%Umdr/6d1)
    ValW = 0d0
   ELSE
    CALL TransformPointToNonparallelRotAxis(x,y,z,XB,YB,ZB,-1d0)
    ValU =  -DBLE(myProcess%iInd*myProcess%ind)*myTwoPI*YB*(myProcess%Umdr/6d1)
    ValV =   DBLE(myProcess%iInd*myProcess%ind)*myTwoPI*XB*(myProcess%Umdr/6d1)
   END IF
  CASE (102) ! Y negativ
   IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
    ValU =  -DBLE(myProcess%ind)*myTwoPI*(Y+mySigma%a/2d0)*(myProcess%Umdr/6d1)
    ValV =   DBLE(myProcess%ind)*myTwoPI*(X)*(myProcess%Umdr/6d1)
    ValW = 0d0
   ELSE
    CALL TransformPointToNonparallelRotAxis(x,y,z,XB,YB,ZB,+1d0)
    ValU =  -DBLE(myProcess%ind)*myTwoPI*YB*(myProcess%Umdr/6d1)
    ValV =   DBLE(myProcess%ind)*myTwoPI*XB*(myProcess%Umdr/6d1)
   END IF
  CASE (103) ! Y negativ
   ValU =  -DBLE(myProcess%ind)*myTwoPI*Y*(myProcess%Umdr/6d1)
   ValV =   DBLE(myProcess%ind)*myTwoPI*X*(myProcess%Umdr/6d1)
   ValW =   0d0
!    ValU =  0d0
!    ValV =  0d0 
!    ValW =  // Forward 8d0 // Backward -1d0 !cm/s      ///!-DBLE(myProcess%ind)*myProcess%Umdr/1d1 !mm/s
  CASE (104) ! Y negativ
   dYShift   = 2d0*0.88d0*dcos(myProcess%Angle*myTwoPI/360d0)
   
   timeLevel = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)
 !  (myProcess%Angle/360d0) = timeLevel*(myProcess%Umdr/60d0)
   dYShiftdt  = 2d0*0.88d0*dcos(timeLevel*(myProcess%Umdr/60d0)*myTwoPI)
   
   

   dYShift    = 2d0*0.88d0*dcos(myProcess%Angle*myTwoPI/360d0)
   dYShiftdt  = (myProcess%Umdr/60d0)*myTwoPI*2d0*0.88d0*dsin(timeLevel*(myProcess%Umdr/60d0)*myTwoPI)
   
   YB = Y + dYShift
   XB = X 
   ValU =  -DBLE(myProcess%ind)*myTwoPI*YB*(myProcess%Umdr/6d1)
   ValV =   DBLE(myProcess%ind)*myTwoPI*XB*(myProcess%Umdr/6d1)
   ValW =   0d0
   ValV =   ValV + dYShiftdt
   
 !  write(*,*) dYShiftdt
   
 END SELECT

ELSE
 IF (iP.gt.300) THEN
    iSeg = iP - 300
 
    valU = mySigma%mySegment(iSeg)%FBMVeloBC(1)
    valV = mySigma%mySegment(iSeg)%FBMVeloBC(2)
    valW = mySigma%mySegment(iSeg)%FBMVeloBC(3)
    
 ELSE
    iSeg = iP - 200
    
    IF (ieee_is_finite(mySigma%mySegment(iSeg)%SegRotFreq)) then
     ValU =  -DBLE(myProcess%ind)*myTwoPI*Y*(mySigma%mySegment(iSeg)%SegRotFreq/6d1)
     ValV =   DBLE(myProcess%ind)*myTwoPI*X*(mySigma%mySegment(iSeg)%SegRotFreq/6d1)
     ValW =   0d0
    else
     ValU =  -DBLE(myProcess%ind)*myTwoPI*Y*(myProcess%Umdr/6d1)
     ValV =   DBLE(myProcess%ind)*myTwoPI*X*(myProcess%Umdr/6d1)
     ValW =   0d0
    end if
  END IF
  
END IF

END SUBROUTINE GetVeloMixerVal
!
!------------------------------------------------------------
!
FUNCTION ViscosityModel(NormShearSquare,Temperature)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
IMPLICIT NONE

real*8 :: ViscosityModel
real*8, intent (in) :: NormShearSquare
! integer, intent (in) :: iMat
real*8, intent (in), optional :: Temperature
REAL*8 :: dStrs, aT, dLimStrs
REAL*8 :: VNN,daux
REAL*8 :: dN
INTEGER iMat
TYPE(tRheology), POINTER :: myRheology

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
aT = 1d0

iMat = myMultiMat%initMaterial
myRheology => myMultiMat%Mat(iMat)%Rheology

! C1C2
if (present(Temperature)) then
 IF (myRheology%AtFunc.EQ.2) THEN
  daux = - myRheology%C1*(Temperature-myRheology%Tb)/(myRheology%C2 + Temperature- myRheology%Tb)
  aT = EXP(daux)
 END IF

 ! TBTS
 IF (myRheology%AtFunc.EQ.3) THEN
  daux = myRheology%C1*(myRheology%TB-myRheology%TS)/(myRheology%C2 + myRheology%TB - myRheology%TS)
  daux = daux - myRheology%C1*(Temperature-myRheology%TS)/(myRheology%C2 + Temperature- myRheology%TS)
!   aT = EXP(daux)
  aT = 1d1**daux
 END IF
 
!ETB
 IF (myRheology%AtFunc.EQ.4) THEN
  daux = (myRheology%E/8.314d0)*( 1d0/(Temperature+273.15d0) - 1d0/(myRheology%TB+273.15d0))
  aT = EXP(daux)
 END IF
 
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

dStrs = (2d0*NormShearSquare)**0.5d0

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
 VNN = (1d1*myRheology%A)/(1d0+(dLimStrs/myRheology%B)**(myRheology%C-1d0)) 
END IF

! HogenPowerLaw
IF (myRheology%Equation.EQ.5) THEN
 dN = Properties%PowerLawExp-1d0
 VNN = Properties%Viscosity(1)*(1d-4 + NormShearSquare)**dN
END IF

! Bingham
IF (myRheology%Equation.EQ.6) THEN
 VNN = (1d1*myRheology%A)*aT*(1d0+myRheology%B*aT*dLimStrs)**(-myRheology%C)
END IF

! WRITE(*,*) dLimStrs,myRheology%eta_max,myRheology%eta_min
ViscosityModel = VNN

! IF (myRheology%Equation.EQ.2) THEN
!  ViscosityModel = MIN(1d1*myRheology%ViscoMax,MAX(1d1*myRheology%ViscoMin,VNN))
! else
!  ViscosityModel = MIN(1d1*1d9,MAX(1d1*myRheology%ViscoMin,VNN))
! END IF

myRheology => NULL()

RETURN

end function ViscosityModel
!
!------------------------------------------------------------
!
FUNCTION AlphaViscosityMatModel(NormShearSquare,iMat,Temperature)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
USE PP3D_MPI, ONLY:myid
IMPLICIT NONE

real*8 :: AlphaViscosityMatModel
real*8, intent (in) :: NormShearSquare
! real*8, intent (in) :: dAlpha
integer, intent (in)  :: iMAt
real*8, intent (in), optional :: Temperature
REAL*8 :: dStrs, aT,log_aT,dLimStrs,MF
REAL*8 :: VNN,daux
REAL*8 :: dN
TYPE(tRheology), POINTER :: myRheology


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
aT = 1d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! if (dAlpha.gt.0.5d0) then
!  iMat=1
! else
!  iMat=2
! end if
 
myRheology => myMultiMat%Mat(iMat)%Rheology

! C1C2
if (present(Temperature)) then
 IF (myRheology%AtFunc.EQ.2) THEN
  daux = - myRheology%C1*(Temperature-myRheology%Tb)/(myRheology%C2 + Temperature- myRheology%Tb)
  aT = EXP(daux)
 END IF

 ! TBTS
 IF (myRheology%AtFunc.EQ.3) THEN
  daux = myRheology%C1*(myRheology%TB-myRheology%TS)/(myRheology%C2 + myRheology%TB - myRheology%TS) - &
         myRheology%C1*(Temperature-myRheology%TS)/(myRheology%C2 + Temperature- myRheology%TS)
  aT = 1d1**daux
 END IF
 
 ! MeltTBTS
 IF (myRheology%AtFunc.EQ.4) THEN

  CALL MeltFunction_MF(MF,Temperature)
  
  log_aT = myRheology%C1*(myRheology%TB-myRheology%TS)/(myRheology%C2 + myRheology%TB - myRheology%TS) - &
           myRheology%C1*(Temperature-myRheology%TS)/(myRheology%C2 + Temperature- myRheology%TS)

  aT = 1d1**((1d0-MF)*myRheology%log_aT_Tilde_Max + MF*log_aT)
  
 END IF

 !ETB myRheology%E is in J/mol
 IF (myRheology%AtFunc.EQ.5) THEN
  daux = (myRheology%E/8.314d0)*( 1d0/(Temperature+273.15d0) - 1d0/(myRheology%TB+273.15d0))
  aT = EXP(daux)
 END IF
 
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

dStrs = (2d0*NormShearSquare)**0.5d0

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
 VNN = (1d1*myRheology%A)/(1d0+(dLimStrs/myRheology%B)**(myRheology%C-1d0)) 
END IF

!MAS
IF (myRheology%Equation.EQ.7) THEN
 VNN = (1d1*myRheology%A)*(1d0+((myRheology%B*dLimStrs)**myRheology%D))**((myRheology%C-1d0)/myRheology%D) 
END IF

!Yasuda
IF (myRheology%Equation.EQ.8) THEN
 VNN = (1d1*myRheology%A)*aT*(1d0+((myRheology%B*aT*dLimStrs)**myRheology%D))**((myRheology%C-1d0)/myRheology%D) 
END IF

! HogenPowerLaw
IF (myRheology%Equation.EQ.5) THEN
 dN = Properties%PowerLawExp-1d0
 VNN = Properties%Viscosity(1)*(1d-4 + NormShearSquare)**dN
END IF

! Bingham
IF (myRheology%Equation.EQ.6) THEN
 VNN = 1d1*(myRheology%A + myRheology%C/(myRheology%B + dLimStrs) )

! VNN = (1d1*myRheology%A)*aT*(1d0+myRheology%B*aT*dLimStrs)**(-myRheology%C)
END IF

AlphaViscosityMatModel = VNN

myRheology => NULL()

RETURN

end function AlphaViscosityMatModel
!
!
!
SUBROUTINE TransformPointToNonparallelRotAxis(x1,y1,z1,x2,y2,z2,dS)
USE Sigma_User, ONLY: mySigma
REAL*8 x1,y1,z1,x2,y2,z2,dS
! REAL*8 :: dAlpha=0.9d0*datan(1d0)/45d0, RotCenter = 501.82d0
REAL*8 xb,yb,zb,xt,yt,zt

XB = x1
YB = y1
ZB = z1 - mySigma%RotAxisCenter

XT = XB
YT = YB*cos(dS*mySigma%RotAxisAngle) - ZB*sin(dS*mySigma%RotAxisAngle)
ZT = YB*sin(dS*mySigma%RotAxisAngle) + ZB*cos(dS*mySigma%RotAxisAngle)

x2 = XT
y2 = YT
z2 = ZT + mySigma%RotAxisCenter

END SUBROUTINE TransformPointToNonparallelRotAxis
!
!
!
SUBROUTINE SetInitialTemperature(T,Coor,ndof)
USE Sigma_User, ONLY: myProcess
implicit none
integer ndof
real*8 T(*),coor(3,*)
integer i
real*8 X,Y,Z,distance

DO i=1,ndof

 X = coor(1,i)
 Y = coor(2,i)
 Z = coor(3,i)
 
  !!!!!!!!!!!!!!!!!!!! VEKA !!!!!!!!!!!!!!!!!!!!!!!
!  distance = ( 2d0*(X+0d0)**2d0 + 6d0*(Y+0d0)**2d0 + (Z-39.0d0)**2d0 )**0.5d0
!  distance = max(13.5d0,distance)
!  T(i) = 200d0 - 15d0 * (distance-13.5d0)/42.7d0
 
 T(i) = myProcess%T0 + Z*myProcess%T0_slope
 
 
end do

! Temperature = myProcess%T0

END SUBROUTINE SetInitialTemperature
!
! ----------------------------------------------
!
SUBROUTINE MeltFunction_MF(MF,T)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
real*8, intent(out) :: MF
real*8, intent (in) :: T
!---------------------------------------------------
REAL*8 :: TM=403.15d0-273.15d0,TS=0.44d0,MS=0.15d0

 MF = (0.5d0*((tanh((T-TM)*TS)) + 1d0))**MS

END SUBROUTINE MeltFunction_MF
!
! ----------------------------------------------
!
SUBROUTINE MeltFunction_Rho(Rho,T,MF)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
real*8, intent(out) :: Rho
real*8, intent (in) :: T,MF
!---------------------------------------------------
REAL*8 :: mV_s=7.97d-7, nV_s=7.83d-4
REAL*8 :: mV_l=1.06d-6, nV_l=8.03d-4
REAL*8 :: ES=0.5d0
REAL*8 :: V_s,V_l,V

 V_s = mV_s*(T - 273.15d0) + nV_s
 V_l = mV_l*(T - 273.15d0) + nV_l
 
 V = V_s*(1d0 - MF**(ES)) + V_s*(MF)**(ES)
 
 Rho = 1d0/V ! kg/m3 

END SUBROUTINE MeltFunction_Rho
!
! ----------------------------------------------
!
SUBROUTINE MeltFunction_Cp(Cp,T,MF)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
real*8, intent(out) :: Cp
real*8, intent (in) :: T,MF
!---------------------------------------------------
REAL*8 :: Cp_s=1.92d3, nH_s=-515.3d3
REAL*8 :: Cp_l=2.89d3, nH_l=-690.0d3

 ! (Cp_s*T + n_s)*(1-MF) + (Cp_l*T + n_l)*MF =
 ! (1-MF)*Cp_s*T + (1-MF)*n_s + MF*Cp_l*T + MF*n_l = 
 ! [(1-MF)*Cp_s + MF*Cp_l]*T + [(1-MF)*n_s + MF*n_l] =
 ! Cp*T + H_ref
 
 Cp = (1d0-MF)*Cp_s + MF*Cp_l ! J/g/K

END SUBROUTINE MeltFunction_Cp
!
! ----------------------------------------------
!
SUBROUTINE MeltFunction_Lambda(Lambda,T,MF)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
real*8, intent(out) :: Lambda
real*8, intent (in) :: T,MF
!---------------------------------------------------
REAL*8 :: HS=0.40 !W/(m*K)
REAL*8 :: m=0.0 !W/(m*K*K)
REAL*8 :: n=0.260 !W/(m*K)
REAL*8 :: T_0=273.15 !K
REAL*8 :: WS=1.39D-5 !1/(K*K)
REAL*8 :: e=2.7182818284d0
REAL*8 :: Lambda_l,Lambda_s

 Lambda_s = HS*(e**(-WS*((T-273.15)-T_0)**2d0))
 Lambda_l = m*(T-273.15) + n
 
 Lambda = (1d0-MF)*Lambda_s + MF*Lambda_l ! W/(m*K)

END SUBROUTINE MeltFunction_Lambda
!
! ----------------------------------------------
!
FUNCTION WallSlip(dShell,dScrew,iMat,Tau)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
USE PP3D_MPI, ONLY:myid
IMPLICIT NONE

real*8 :: WallSlip
real*8, intent (in) :: dShell,dScrew,Tau
integer, intent(in) :: iMat

real*8 d,d_min,d_max,tau_min,tau_max,tau_factor,d_factor,slip

d_min = 0d0
d_max = 1d1*myMultiMat%Mat(iMat)%Rheology%WS_d          ! scaling to mm 
tau_min = 1d1*myMultiMat%Mat(iMat)%Rheology%WS_TauMin   ! scaling from Pa to 10Pa
tau_max = 1d1*myMultiMat%Mat(iMat)%Rheology%WS_TauMax   ! scaling from Pa to 10Pa
slip = myMultiMat%Mat(iMat)%Rheology%WS_SlipFactor

! d = max(dShell,dScrew)                                  ! mm
d = dShell

d_factor   = 1d0-((max(min(d,d_max),d_min)-d_min)/(d_max-d_min))
tau_factor = ((max(min(tau,tau_max),tau_min)-tau_min)/(tau_max-tau_min))

WallSlip = 1d0-slip*d_factor*tau_factor
! write(*,*) WallSlip 

END FUNCTION WallSlip
