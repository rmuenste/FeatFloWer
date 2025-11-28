! -----------------------------------------------------
!
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
USE var_QuadScalar, ONLY : b2DViscoBench,b3DViscoBench, referenceVelocity, GammaDot
use fbm, only: fbm_updateFBM, fbm_velBCTest, fbm_velValue
REAL*8 X,Y,Z,ValU,ValV,ValW
REAL*8 PX,PY,dScale
INTEGER iC
REAL*8 :: RX = 0.0d0,RY = 0.0d0,RZ = 0.0d0, RAD = 0.245d0
REAL*8 :: R_inflow=4d0
REAL*8 :: PI=3.141592654d0

!ValU = Z * 10.0 
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
use var_QuadScalar, only : myFBM, referenceVelocity, GammaDot
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
real*8 dRPM
REAL*8 ValUT,ValVT,ValWT

ValU = 0d0
ValV = 0d0
ValW = 0d0

IF (iT.EQ.23) THEN
 ValU= 0.0d0 
 ValV= 0d0
 ValW= 0d0
END IF

IF (iT.lt.0) THEN
  
 iInflow = ABS(iT)
 
 IF (ALLOCATED(myProcess%myInflow)) then
  IF (SIZE(myProcess%myInflow).ge.iInflow) then
 
   IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
   
    IF (myProcess%myInflow(iInflow)%iBCType.eq.1) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     dDensity      = myMultiMat%Mat(iMat)%Thermodyn%density
!      dTemperature  = myProcess%myInflow(iInflow)%Temperature
!      dDensitySlope = myMultiMat%Mat(iMat)%Thermodyn%densitySteig
!      dDensity      = dDensity - dTemperature*dDensitySlope
     douterradius  = myProcess%myInflow(iInflow)%outerradius
     dinnerradius  = myProcess%myInflow(iInflow)%innerradius
     dProfil = RotParabolicVelo3D(dMassFlow,dDensity,dOuterRadius)
     ValU = dProfil(1)
     ValV = dProfil(2)
     ValW = dProfil(3)
    END IF

    IF (myProcess%myInflow(iInflow)%iBCType.eq.2) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     dDensity      = myMultiMat%Mat(iMat)%Thermodyn%density
!      dTemperature  = myProcess%myInflow(iInflow)%Temperature
!      dDensitySlope = myMultiMat%Mat(iMat)%Thermodyn%densitySteig
!      dDensity      = dDensity - dTemperature*dDensitySlope
     douterradius  = myProcess%myInflow(iInflow)%outerradius
     dinnerradius  = myProcess%myInflow(iInflow)%innerradius
     dProfil = RotDoubleParabolicVelo3D(dMassFlow,dDensity,dInnerRadius,dOuterRadius)
     ValU = dProfil(1)
     ValV = dProfil(2)
     ValW = dProfil(3)
    END IF
    
    IF (myProcess%myInflow(iInflow)%iBCType.eq.3) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     dDensity      = myMultiMat%Mat(iMat)%Thermodyn%density
!      dTemperature  = myProcess%myInflow(iInflow)%Temperature
!      dDensitySlope = myMultiMat%Mat(iMat)%Thermodyn%densitySteig
!      dDensity      = dDensity - dTemperature*dDensitySlope
     douterradius  = myProcess%myInflow(iInflow)%outerradius
     dProfil = FlatVelo3D(dMassFlow,dDensity,dOuterRadius)
     ValU = dProfil(1)
     ValV = dProfil(2)
     ValW = dProfil(3)
    END IF
    
    IF (myProcess%myInflow(iInflow)%iBCType.eq.4) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     dDensity      = myMultiMat%Mat(iMat)%Thermodyn%density
!      dTemperature  = myProcess%myInflow(iInflow)%Temperature
!      dDensitySlope = myMultiMat%Mat(iMat)%Thermodyn%densitySteig
!      dDensity      = dDensity - dTemperature*dDensitySlope
     douterradius  = myProcess%myInflow(iInflow)%outerradius
     dinnerradius  = myProcess%myInflow(iInflow)%innerradius
     dProfil = CurvedFlatVelo3D(dMassFlow,dDensity,dInnerRadius,dOuterRadius)
     ValU = dProfil(1)
     ValV = dProfil(2)
     ValW = dProfil(3)
    END IF

    IF (myProcess%myInflow(iInflow)%iBCType.eq.5) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     dDensity      = myMultiMat%Mat(iMat)%Thermodyn%density
!      dTemperature  = myProcess%myInflow(iInflow)%Temperature
!      dDensitySlope = myMultiMat%Mat(iMat)%Thermodyn%densitySteig
!      dDensity      = dDensity - dTemperature*dDensitySlope
     dMidpointA    = myProcess%myInflow(iInflow)%midpointA
     dMidpointB    = myProcess%myInflow(iInflow)%midpointB

     dProfil = RectangleVelo3D(dMassFlow,dDensity,dMidpointA,dMidpointB)
     ValU = dProfil(1)
     ValV = dProfil(2)
     ValW = dProfil(3)
    END IF

    IF (myProcess%myInflow(iInflow)%iBCType.eq.6) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     dDensity      = myMultiMat%Mat(iMat)%Thermodyn%density
!      dTemperature  = myProcess%myInflow(iInflow)%Temperature
!      dDensitySlope = myMultiMat%Mat(iMat)%Thermodyn%densitySteig
!      dDensity      = dDensity - dTemperature*dDensitySlope
     dMidpointA    = myProcess%myInflow(iInflow)%midpointA
     dMidpointB    = myProcess%myInflow(iInflow)%midpointB

     dProfil = CurvedRectangleVelo3D(dMassFlow,dDensity,dMidpointA,dMidpointB)
     ValU = dProfil(1)
     ValV = dProfil(2)
     ValW = dProfil(3)
    END IF
    
   ELSE
   
    DO iSubInflow=1,myProcess%myInflow(iInflow)%nSubInflows
    
     IF (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCType.eq.1) then
      dCenter       = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%center
      dNormal       = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%normal
      dMassFlow     = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate
      iMat          = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Material
      dDensity      = myMultiMat%Mat(iMat)%Thermodyn%density
!       dTemperature  = myProcess%myInflow(iInflow)%Temperature
!       dDensitySlope = myMultiMat%Mat(iMat)%Thermodyn%densitySteig
!       dDensity      = dDensity - dTemperature*dDensitySlope
      douterradius  = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius
      dinnerradius  = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius
      dProfil = RotParabolicVelo3D(dMassFlow,dDensity,dOuterRadius)
      daux = dProfil(1)**2d0 + dProfil(2)**2d0 + dProfil(3)**2d0
      if (daux.ne.0d0) then
       ValU = dProfil(1)
       ValV = dProfil(2)
       ValW = dProfil(3)
      END IF
     END IF
     
    END DO
    
   END IF
  else
   write(*,*) 'Inflow array is not allocated!!', iInflow,SIZE(myProcess%myInflow)
   stop
  END IF
 else
  write(*,*) 'Size of allocated inflow array is smaller then requred index:',iInflow
  stop
 END IF
END IF

IF (iT.EQ.111) THEN
   dist = SQRT((x-27)**2d0 + y*y)
   dScale = 5.0d0/0.245**2d0
   IF (dist.lt.R_inflow) THEN
    ValW= -dScale*(dist+0.245)*(0.245-dist)
   END IF
END IF

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

IF (iT.EQ.4.OR.iT.EQ.5) THEN
 dScale=3d0/2d0
 ValU=dScale*(1d0-0.25d0*Y*Y)
 ValV= 0d0
 ValW= 0d0
END IF

! 2D FAC  with Slip
IF (iT.EQ.6) THEN
    dCenter       = [0.05d0, 0.055d0, 0.0d0]
    dNormal       = [1.0d0, 1.0d0, 0.0d0]
    dMassFlow     = 0.5d0
    ddensity      = 1.0d0
    douterradius  = 0.200d0
    dProfil = FlatVelo3D(dMassFlow,dDensity,dOuterRadius)
    ValU = dProfil(1)
    ValV = dProfil(2)
    ValW = dProfil(3)
END IF

IF (iT.EQ.770) THEN
  ValU =  -myTwoPI*Y*(1d0/6d1)
  ValV =   myTwoPI*X*(1d0/6d1)
  ValW =   0.0d0
END IF

IF (iT.EQ.666) THEN
  dCenter=[-0.85d0,0d0,0d0]
  DIST = SQRT((X-dCenter(1))**2d0+(Y-dCenter(2))**2d0)
  DIST = MIN(DIST,0.18d0)
  DAUX= 2d0*(DIST+0.18d0)*(0.18d0-DIST)/(0.18d0*0.18d0)
  ValU =   0.3d0*DAUX
  ValV =   0.0d0*DAUX
  ValW =   0.3d0*DAUX
END IF

IF (iT.EQ.778) THEN
  dScale=1.0d3*(9d0/4d0)/(0.15d0*0.0425d0)**2d0
  XX = X-9.6d0
  ZZ = Z
  ValU =   0.0d0
  ValV =   dScale*ZZ*(0.085d0-ZZ)*XX*(0.300d0-XX)
  ValW =   0.0d0
END IF

IF (iT.EQ.771) THEN
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

IF (iT.EQ.8.OR.IT.EQ.9) THEN
 ValW = 1d0
END IF

if(it.eq.10)then

  IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
   dInnerRadius = 0.5d0*mySigma%Dz_In

  IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
   dVolFlow = (1e3/3.6d3)*myProcess%Massestrom/(myThermodyn%density) ! cm3/s
   daux = (PI/6d0)*(dInnerRadius+mySigma%a/2d0)*((mySigma%a/2d0-dInnerRadius)**3d0)
   dScale = (dVolFlow/2d0)/daux
   IF (Y.LT.0) THEN
     DIST = SQRT(X**2d0+(Y+mySigma%a/2d0)**2d0)
    ELSE
     DIST = SQRT(X**2d0+(Y-mySigma%a/2d0)**2d0)
    END IF
    IF (DIST.GT.dInnerRadius.AND.DIST.LT.mySigma%a/2d0) THEN
     ValW= dScale*(DIST-dInnerRadius)*(mySigma%a/2d0-DIST)
    END IF
  ELSE
   CALL TransformPointToNonparallelRotAxis(0d0,0d0,z,XX,YY,ZZ,+1d0)
   dInnerInflowRadius = abs(YY)
   dVolFlow = (1e3/3.6d3)*myProcess%Massestrom/(myThermodyn%density) ! cm3/s
   daux = (PI/6d0)*(dInnerRadius+dInnerInflowRadius)*((dInnerInflowRadius-dInnerRadius)**3d0)
   dScale = (dVolFlow/2d0)/daux
   
   IF (Y.LT.0) CALL TransformPointToNonparallelRotAxis(x,y,z,XX,YY,ZZ,+1d0)
   IF (Y.GT.0) CALL TransformPointToNonparallelRotAxis(x,y,z,XX,YY,ZZ,-1d0)
    
   DIST = SQRT(XX**2d0 + YY**2d0)
   
   IF (DIST.GT.dInnerRadius.AND.DIST.LT.dInnerInflowRadius) THEN
    ValW= dScale*(DIST-dInnerRadius)*(dInnerInflowRadius-DIST)
   END IF
  END IF

   
  END IF

  IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE") THEN
   dInnerRadius = myProcess%MinInflowDiameter*0.5d0
   dOuterRadius = myProcess%MaxInflowDiameter*0.5d0

   dVolFlow = (1e3/3.6d3)*myProcess%Massestrom/(myThermodyn%density) ! cm3/s
   daux = (PI/6d0)*(dInnerRadius+dOuterRadius)*((dOuterRadius-dInnerRadius)**3d0)
   dScale = (dVolFlow/1d0)/daux

   DIST = SQRT(X**2d0+Y**2d0)
   
   IF (DIST.GT.dInnerRadius.AND.DIST.LT.dOuterRadius) THEN
    ValW= dScale*(DIST-dInnerRadius)*(dOuterRadius-DIST)
   END IF
  END IF
end if

! QUESTIONABLE CASE WHEN THE INFLOW FOR TSE/SSE experiences some strange things

! if(it.eq.100)then
!   IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
!    dInnerRadius = 0.5d0*mySigma%Dz_In
!
!   IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
!    dVolFlow = (1e3/3.6d3)*myProcess%Massestrom/(myThermodyn%density) ! cm3/s
!    daux = (PI/6d0)*(dInnerRadius+mySigma%a/2d0)*((mySigma%a/2d0-dInnerRadius)**3d0)
!    dScale = (dVolFlow/2d0)/daux
!    IF (Y.LT.0) THEN
!      DIST = SQRT(X**2d0+(Y+mySigma%a/2d0)**2d0)
!     ELSE
!      DIST = SQRT(X**2d0+(Y-mySigma%a/2d0)**2d0)
!     END IF
!     IF (DIST.GT.dInnerRadius.AND.DIST.LT.mySigma%a/2d0) THEN
!      ValW= dScale*(DIST-dInnerRadius)*(mySigma%a/2d0-DIST)
!     END IF
!     IF (DIST.LT.dInnerRadius.and.Y.gt.0d0) THEN
!      ValU= -myTwoPI*(Y-mySigma%a/2d0)*(myProcess%Umdr/6d1)*REAL(myProcess%iInd*myProcess%ind)
!      ValV = myTwoPI*X*(myProcess%Umdr/6d1)*REAL(myProcess%iInd*myProcess%ind)
!     END IF
!     IF (DIST.LT.dInnerRadius.and.Y.lt.0d0) THEN
!      ValU= -myTwoPI*(Y+mySigma%a/2d0)*(myProcess%Umdr/6d1)*REAL(myProcess%iInd*myProcess%ind)
!      ValV = myTwoPI*X*(myProcess%Umdr/6d1)*REAL(myProcess%iInd*myProcess%ind)
!     END IF
!   ELSE
!    CALL TransformPointToNonparallelRotAxis(0d0,0d0,z,XX,YY,ZZ,+1d0)
!    dInnerInflowRadius = abs(YY)
!    dVolFlow = (1e3/3.6d3)*myProcess%Massestrom/(myThermodyn%density) ! cm3/s
!    daux = (PI/6d0)*(dInnerRadius+dInnerInflowRadius)*((dInnerInflowRadius-dInnerRadius)**3d0)
!    dScale = (dVolFlow/2d0)/daux
!
!    IF (Y.LT.0) CALL TransformPointToNonparallelRotAxis(x,y,z,XX,YY,ZZ,+1d0)
!    IF (Y.GT.0) CALL TransformPointToNonparallelRotAxis(x,y,z,XX,YY,ZZ,-1d0)
!
!    DIST = SQRT(XX**2d0 + YY**2d0)
!
!    IF (DIST.GT.dInnerRadius.AND.DIST.LT.dInnerInflowRadius) THEN
!     ValW= dScale*(DIST-dInnerRadius)*(dInnerInflowRadius-DIST)
!    END IF
!    IF (DIST.LT.dInnerRadius) THEN
! ! ! ! ! !     IF (DIST.LT.dInnerRadius.and.Y.gt.0d0) THEN
! ! ! ! ! !      ValU=
! ! ! ! ! !      ValV =
! ! ! ! ! !     END IF
! ! ! ! ! !     IF (DIST.LT.dInnerRadius.and.Y.lt.0d0) THEN
! ! ! ! ! !      ValU=
! ! ! ! ! !      ValV =
! ! ! ! ! !     END IF
!    END IF
!   END IF
!
!
!   END IF
!
!   IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE") THEN
!    dInnerRadius = myProcess%MinInflowDiameter*0.5d0
!    dOuterRadius = myProcess%MaxInflowDiameter*0.5d0
!
!    dVolFlow = (1e3/3.6d3)*myProcess%Massestrom/(myThermodyn%density) ! cm3/s
!    daux = (PI/6d0)*(dInnerRadius+dOuterRadius)*((dOuterRadius-dInnerRadius)**3d0)
!    dScale = (dVolFlow/1d0)/daux
!
!    DIST = SQRT(X**2d0+Y**2d0)
!
!    IF (DIST.GT.dInnerRadius.AND.DIST.LT.dOuterRadius) THEN
!     ValW= dScale*(DIST-dInnerRadius)*(dOuterRadius-DIST)
!    END IF
!    IF (DIST.LT.dInnerRadius) THEN
!     ValU =  -DBLE(myProcess%ind)*myTwoPI*Y*(myProcess%Umdr/6d1)
!     ValV =   DBLE(myProcess%ind)*myTwoPI*X*(myProcess%Umdr/6d1)
!    END IF
!   END IF
! end if

if(it.eq.11)then
  IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
   ValU= -myTwoPI*(Y-mySigma%a/2d0)*(myProcess%Umdr/6d1)*REAL(myProcess%iInd*myProcess%ind)
   ValV = myTwoPI*X*(myProcess%Umdr/6d1)*REAL(myProcess%iInd*myProcess%ind)
   ValW = 0d0
  ELSE
   CALL TransformPointToNonparallelRotAxis(x,y,z,XX,YY,ZZ,-1d0)
   ValU =  -DBLE(myProcess%iInd*myProcess%ind)*myTwoPI*YY*(myProcess%Umdr/6d1)
   ValV =   DBLE(myProcess%iInd*myProcess%ind)*myTwoPI*XX*(myProcess%Umdr/6d1)
   ValW = 0d0
   CALL TransformVelocityToNonparallelRotAxis(ValU,ValV,ValW,ValUT,ValVT,ValWT,-1d0)
   ValU = ValUT
   ValV = ValVT
   ValW = ValWT
  END IF
end if

if(it.eq.12)then
  IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
   ValU= -myTwoPI*(Y+mySigma%a/2d0)*(myProcess%Umdr/6d1)*REAL(myProcess%ind)
   ValV = myTwoPI*X*(myProcess%Umdr/6d1)*REAL(myProcess%ind)
   ValW = 0d0
  ELSE
   CALL TransformPointToNonparallelRotAxis(x,y,z,XX,YY,ZZ,+1d0)
   ValU =  -DBLE(myProcess%ind)*myTwoPI*YY*(myProcess%Umdr/6d1)
   ValV =   DBLE(myProcess%ind)*myTwoPI*XX*(myProcess%Umdr/6d1)
   ValW = 0d0
   CALL TransformVelocityToNonparallelRotAxis(ValU,ValV,ValW,ValUT,ValVT,ValWT,+1d0)
   ValU = ValUT
   ValV = ValVT
   ValW = ValWT
  END IF
end if

if(it.eq.13)then
  ValU =  -DBLE(myProcess%ind)*myTwoPI*Y*(myProcess%Umdr/6d1)
  ValV =   DBLE(myProcess%ind)*myTwoPI*X*(myProcess%Umdr/6d1)
  ValW =   0d0
end if

IF (iT.EQ.20) THEN

   dOuterRadius = myProcess%MaxInflowDiameter*0.5d0 !cm
   dVolFlow = myProcess%Massestrom ! kg/h
   dDensity = myThermodyn%density ! g/cm3
   ValU =  -DBLE(myProcess%ind)*myTwoPI*Y*(myProcess%Umdr/6d1)
   ValV =   DBLE(myProcess%ind)*myTwoPI*X*(myProcess%Umdr/6d1)
   if (x.gt.0d0) then
    ValW=RotParabolicVelo2Dz(+0.395d0,-0.065d0,dVolFlow,dDensity,dOuterRadius)
   else
    ValW=RotParabolicVelo2Dz(-0.395d0,+0.065d0,dVolFlow,dDensity,dOuterRadius)
   end if
   
END IF

!! InitBoundaryStructure: changed from I1 to I2 !!
! used for Channel / BDF(2) validation
IF (iT.EQ.21) THEN
 dScale=0.2d0*(3d0/2d0)/(0.205d0*0.205d0)*t
 ValU=dScale*Y*(0.41d0-Y)
 ValV= 0d0
 ValW= 0d0
END IF
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
IF (iT.EQ.28) THEN
 dScale=(9d0/4d0)/(0.205d0*0.205d0*0.205d0*0.205d0)*sin(t*PI/8d0)
 ValU=dScale*Y*(0.41d0-Y)*Z*(0.41d0-Z)
 ValV= 0d0
 ValW= 0d0
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
  dScale = (0.1d0+0.9d0*abs(sin(t*2d0*PI/8d0)))/(R_inflow*R_inflow)
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

!STENT inflow 1
IF (iT.EQ.38) THEN
! R_inflow = 0.250d0
! dist = SQRT((x-(-2.27d0))**2d0 + (y+0.00d0)**2d0 + (z-(-0.07d0))**2d0)
! IF (dist.lt.R_inflow) THEN
!  dScale = 1d3
!   dScale = (0.1d0+0.9d0*abs(sin(t*2d0*PI/8d0)))/(R_inflow*R_inflow)
!  dScale = dScale*(dist+R_inflow)*(R_inflow-dist)
!  ValU= +1.0000d0*dScale
!  ValV=  0.0000d0*dScale
!  ValW=  0.0000d0*dScale

  dist = pulsativeProfile(t)

  dCenter=[-2.27, 0.0, -0.07]
  dNormal=[1.0, 0.0, 0.0]
  dProfil = RotParabolicVelo3D(dist * 60d0,1d0,0.250d0)

  ValU = dProfil(1)
  ValV = dProfil(2)
  ValW = dProfil(3)

! END IF
END IF

!STENT inflow 2
IF (iT.EQ.39) THEN
! R_inflow = 0.175d0
! dist = SQRT((x-(-1.72d0))**2d0 + (y-0.02d0)**2d0 + (z-1.31)**2d0)
!!  write(*,*) dist
! IF (dist.lt.R_inflow) THEN
!  dScale = 1d3
!!   dScale = (0.1d0+0.9d0*abs(sin(t*2d0*PI/8d0)))/(R_inflow*R_inflow)
!  dScale = dScale*(dist+R_inflow)*(R_inflow-dist)
!  ValU= +0.9650d0*dScale
!  ValV= -0.0000d0*dScale
!  ValW= -0.2588d0*dScale
! END IF

  dist = pulsativeProfile(t)

  dCenter=[-1.72, 0.02, 1.31]
  dNormal=[0.965, 0.0, -0.2588]
  dProfil = RotParabolicVelo3D(dist * 20.0d0,1d0,0.175d0)

  ValU = dProfil(1)
  ValV = dProfil(2)
  ValW = dProfil(3)
END IF
 
IF (iT.EQ.444) THEN
 dDensity = myThermodyn%density
 IF (abs(y).gt.2.5d0) THEN
  daux = (3.29d0/3600d0)/(dDensity*1d3)
  dScale = daux/(0.2322d0*0.0012)*1d2 ! cm/s ==> u_mean
  if (y.lt.0d0) THEN
   ValW = max(0d0,3d0/2d0*dScale*(1d0 - (2d0*(y + 5.053d0)/0.12d0)**2d0))
  ELSE
   ValW = max(0d0,3d0/2d0*dScale*(1d0 - (2d0*(y - 5.053d0)/0.12d0)**2d0))
  END IF
 else
  daux = (306.45d0/3600d0)/(dDensity*1d3)
  dScale = 3d0/2d0*daux/(0.2322d0*0.0024)*1d2
  ValW = max(0d0,dScale*(-(2d0*(y + 0.12d0)/0.24d0)*(2d0*(y - 0.12d0)/0.24d0)))
 end if
END IF

IF (iT.EQ.45) THEN
 ValW=RotParabolicVelo2Dz(0d0,0d0,1449d0,1d0,2.495d0)
END IF

IF (iT.EQ.51) THEN
  ValW=RotParabolicVelo2Dz(0d0,0d0,1270d0,1d0,2.495d0)
END IF
IF (iT.EQ.52) THEN
  ValW=RotParabolicVelo2Dz(0d0,6d0,56d0,1d0,1.245d0)
END IF
IF (iT.EQ.53) THEN
  ValW=RotParabolicVelo2Dz(0d0,-6d0,56d0,1d0,1.245d0)
END IF
IF (iT.EQ.54) THEN
  ValW=RotParabolicVelo2Dz(0d0,6d0,67d0,1d0,1.245d0)
END IF
IF (iT.EQ.31) THEN
 dScale=0.2d0*(3d0/2d0)/(0.205d0*0.205d0)*sin(t*PI/8d0)
 ValU=dScale*Y*(0.41d0-Y)
 ValV= 0d0
 ValW= 0d0
END IF

! centroplast
IF (iT.EQ.61) THEN
  ValW=RotParabolicVelo2Dz(0d0,0d0,25d0,1d0,1.245d0)
END IF

! egeplast
IF (iT.EQ.62) THEN
  ValW=RotParabolicVelo2Dz(2.4d0,-4.155d0,7.5d0,1d0,0.39d0)
END IF

! Weber/G1  --> Middle layer
IF (iT.EQ.63) THEN
  ValW=RotParabolicVelo2Dz(-8.46d0,-20.43d0,-400d0,1d0,2.95d0)
END IF

! Weber/G2  --> Outer layer
IF (iT.EQ.64) THEN
  ValW=RotParabolicVelo2Dz(+8.52d0,+20.57d0,-200d0,1d0,1.45d0)
END IF

! Weber/G3  --> Inner layer
IF (iT.EQ.65) THEN
  ValW=RotParabolicVelo2Dz(+0.0d0,+0d0,-200d0,1d0,1.49d0)
!  dCenter=[7.89d0,18.99d0,72.76d0]
!  dNormal=[-0.346835,-0.837317,-0.422617]
!  dProfil = RotParabolicVelo3D(200d0,1d0,1.65d0)
!  ValU = dProfil(1)
!  ValV = dProfil(2)
!  ValW = dProfil(3)
END IF

IF (iT.EQ.70) THEN
 U_bar = 2.0
 dScale=((3d0/2d0)/(0.20d0*0.20d0)) * U_bar
 ValU=dScale*Y*(0.4d0-Y)
 ValV= 0d0
 ValW= 0d0
END IF

! Centroplast - Vollstab
IF (iT.EQ.71) THEN
  ValW=RotParabolicVelo2Dz(+0.0d0,+0.0d0,-24d0,1d0,1.25d0)
END IF
IF (iT.EQ.72) THEN
  ValW=RotParabolicVelo2Dz(+0.0d0,+0.0d0,+24d0,1d0,1.25d0)
END IF

! M+S --> for the meshes prepared by Jens and Raphael
IF (iT.EQ.81) THEN
   dOuterRadius = myProcess%MaxInflowDiameter*0.5d0 !cm
   dVolFlow = myProcess%Massestrom ! kg/h
   dDensity = myThermodyn%density ! g/cm3
   ValW=RotParabolicVelo2Dz(+0d0,+0d0,dVolFlow,dDensity,dOuterRadius)
END IF

! IDE
IF (iT.EQ.82) THEN
  ValV=RotParabolicVelo2Dy(+0.0d0,+0.0d0,-100d0,1d0,3.5d0)
END IF

! Reuter
IF (iT.EQ.86) THEN
  ValW=RotParabolicVelo2Dz(-0.4d0,+0.6d0,-43d0,1d0,2.1d0)
END IF
IF (iT.EQ.87) THEN
 dCenter=[23.6d0,0.6d0,26.85d0]
 dNormal=[-0.965923, 3.10658e-06, -0.25883]
 dProfil = RotParabolicVelo3D(53d0,1d0,2.1d0)
 ValU = dProfil(1)
 ValV = dProfil(2)
 ValW = dProfil(3)

!  dNx = -0.965923
!  dNY = +3.10658e-06     
!  dNz = -0.25883
!  dNorm = SQRT(dNX*dNx +dNy*dNy + dNz*dNz)
!  dNx = dNx/dNorm
!  dNY = dNy/dNorm
!  dNz = dNz/dNorm
!  dist = SQRT((X-23.6d0)*(X-23.6d0) + (Y-0.6d0)*(Y-0.6d0) + (Z-26.85d0)*(Z-26.85d0))
! !  dVectMag = 2d0*my_InFlow%VectMag*(1.12d0-dist*dist)
!  dScale = 5d0*(2d0*max(0d0,(1.12d0-dist*dist)))
!  ValU = dNx*dScale
!  ValV = dNy*dScale
!  ValW = dNz*dScale
END IF

IF (iT.EQ.171) THEN
   dScale = 2.25d0*(Y)*(Y-0.1d0)*(Z)*(Z-0.1d0)/(0.05d0**4d0)
   ValU= +dScale*900d0
END IF

IF (iT.EQ.172) THEN
   dScale = 2.25d0*(Y)*(Y-0.1d0)*(Z)*(Z-0.1d0)/(0.05d0**4d0)
   ValU= -dScale*1300d0
END IF


!!!!!!!!!!!!!!!!!!!!!!!!!! INNOSPIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (iT.EQ.222) THEN
 ValW=-1d0
END IF

IF (iT.EQ.223) THEN
 ValW=RotParabolicVelo2Dz(0.0d0,0.0d0,-1d0,1d0,0.55d0)
END IF

IF (iT.EQ.224) THEN
 ValW=RotParabolicVelo2Dz(1.0d0,0.0d0,-1d0,1d0,0.25d0)
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!! INNOSPIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!

! PP-Weber
IF (iT.EQ.91) THEN
  ValV=RotParabolicVelo2Dy(+0.0d0,+35.2d0,-100d0,1d0,0.5d0)
END IF
IF (iT.EQ.92) THEN
  ValW=RotParabolicVelo2Dz(+0.0d0,+0.0d0,-900d0,1d0,2.5d0)
END IF

! RAIN CARBON
IF (iT.EQ.93) THEN
  ValU=RotParabolicVelo2Dx(+0.0d0,+0.0d0,9217d0,1d0,5.0d0)
END IF

IF (iT.EQ.99) THEN
 ValW = -myFBM%ParticleNew(1)%Velocity(3)
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
REAL*8 dOmega
integer iScrewType,iSeg
REAL*8 ValUT,ValVT,ValWT

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
    CALL TransformVelocityToNonparallelRotAxis(ValU,ValV,ValW,ValUT,ValVT,ValWT,-1d0)
    ValU = ValUT
    ValV = ValVT
    ValW = ValWT
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
    CALL TransformVelocityToNonparallelRotAxis(ValU,ValV,ValW,ValUT,ValVT,ValWT,+1d0)
    ValU = ValUT
    ValV = ValVT
    ValW = ValWT
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

    dOmega = sqrt(mySigma%mySegment(iSeg)%FBMOmegaBC(1)**2d0 + &
                  mySigma%mySegment(iSeg)%FBMOmegaBC(2)**2d0 + &
                  mySigma%mySegment(iSeg)%FBMOmegaBC(3)**2d0)
    if (dOmega.gt.0d0) then
     if (mySigma%mySegment(iSeg)%FBMOmegaBC(1).ne.0d0.and.mySigma%mySegment(iSeg)%FBMOmegaBC(2).eq.0d0.and.mySigma%mySegment(iSeg)%FBMOmegaBC(3).eq.0d0) THEN ! xAxis-rotation
     end if

     if (mySigma%mySegment(iSeg)%FBMOmegaBC(2).ne.0d0.and.mySigma%mySegment(iSeg)%FBMOmegaBC(1).eq.0d0.and.mySigma%mySegment(iSeg)%FBMOmegaBC(3).eq.0d0) THEN ! yAxis-rotation
     end if

     if (mySigma%mySegment(iSeg)%FBMOmegaBC(3).ne.0d0.and.mySigma%mySegment(iSeg)%FBMOmegaBC(2).eq.0d0.and.mySigma%mySegment(iSeg)%FBMOmegaBC(1).eq.0d0) THEN ! zAxis-rotation
      ValU =  -myTwoPI*(Y-mySigma%mySegment(iSeg)%FBMOffsetBC(2))*(mySigma%mySegment(iSeg)%FBMOmegaBC(3)/6d1)
      ValV =   myTwoPI*(X-mySigma%mySegment(iSeg)%FBMOffsetBC(1))*(mySigma%mySegment(iSeg)%FBMOmegaBC(3)/6d1)
      ValW =   0d0
     end if
    end if

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
USE viscosity_model, ONLY : mu_eff,set_mu_f
IMPLICIT NONE

real*8 :: AlphaViscosityMatModel
real*8, intent (in) :: NormShearSquare
! real*8, intent (in) :: dAlpha
integer, intent (in)  :: iMAt
real*8, intent (in), optional :: Temperature
REAL*8 :: dStrs, aT,log_aT,dLimStrs,MF,dLimTemperature
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

 dLimTemperature = MIN(myRheology%TemperatureMax,MAX(myRheology%TemperatureMin,Temperature))

 IF (myRheology%AtFunc.EQ.2) THEN
  daux = - myRheology%C1*(dLimTemperature-myRheology%Tb)/(myRheology%C2 + dLimTemperature- myRheology%Tb)
  aT = EXP(daux)
 END IF

 ! TBTS
 IF (myRheology%AtFunc.EQ.3) THEN
  daux = myRheology%C1*(myRheology%TB-myRheology%TS)/(myRheology%C2 + myRheology%TB - myRheology%TS) - &
         myRheology%C1*(dLimTemperature-myRheology%TS)/(myRheology%C2 + dLimTemperature- myRheology%TS)
  aT = 1d1**daux
 END IF
 
 ! MeltTBTS
 IF (myRheology%AtFunc.EQ.4) THEN

  CALL MeltFunction_MF(MF,dLimTemperature)
  
  log_aT = myRheology%C1*(myRheology%TB-myRheology%TS)/(myRheology%C2 + myRheology%TB - myRheology%TS) - &
           myRheology%C1*(dLimTemperature-myRheology%TS)/(myRheology%C2 + dLimTemperature- myRheology%TS)

  aT = 1d1**((1d0-MF)*myRheology%log_aT_Tilde_Max + MF*log_aT)
  
 END IF

 !ETB myRheology%E is in J/mol
 IF (myRheology%AtFunc.EQ.5) THEN
  daux = (myRheology%E/8.314d0)*( 1d0/(dLimTemperature+273.15d0) - 1d0/(myRheology%TB+273.15d0))
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
!  CALL set_mu_f(Properties%Viscosity(1))
!  VNN = mu_eff(Temperature, dStrs)
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
SUBROUTINE TransformVelocityToNonparallelRotAxis(x1,y1,z1,x2,y2,z2,dS)
USE Sigma_User, ONLY: mySigma
REAL*8 x1,y1,z1,x2,y2,z2,dS

X2 = X1
Y2 = Y1*cos(dS*mySigma%RotAxisAngle) - Z1*sin(dS*mySigma%RotAxisAngle)
Z2 = Y1*sin(dS*mySigma%RotAxisAngle) + Z1*cos(dS*mySigma%RotAxisAngle)

END SUBROUTINE TransformVelocityToNonparallelRotAxis
!
!
!
SUBROUTINE SetInitialTemperature(T,Coor,ndof)
USE Sigma_User, ONLY: myProcess,mySigma
USE PP3D_MPI, ONLY:myid
implicit none
integer ndof
real*8 T(*),coor(3,*)
integer i,iN
real*8 X,Y,Z,distance,dN,daux

DO i=1,ndof

 X = coor(1,i)
 Y = coor(2,i)
 Z = coor(3,i)
 
  !!!!!!!!!!!!!!!!!!!! VEKA !!!!!!!!!!!!!!!!!!!!!!!
!  distance = ( 2d0*(X+0d0)**2d0 + 6d0*(Y+0d0)**2d0 + (Z-39.0d0)**2d0 )**0.5d0
!  distance = max(13.5d0,distance)
!  T(i) = 200d0 - 15d0 * (distance-13.5d0)/42.7d0
 
!   !!!!!!!!!!!!!!!!!!!! KMB Zahradpumpe !!!!!!!!!!!!!!!!!!!!!!!
!    if (abs(z).lt.5.3d0) then
!     T(i) = myProcess%T0
!    else
!     T(i) = myProcess%T0 + 30d0
!    end if
!   !!!!!!!!!!!!!!!!!!!! KMB Zahradpumpe !!!!!!!!!!!!!!!!!!!!!!!

 IF (myProcess%T0_N.gt.0) THEN

  daux = DBLE(myProcess%T0_N-1)*Z/mySigma%L
  iN = floor(1d0+daux)
  if (iN.eq.myProcess%T0_N) then
   iN = iN-1
   dN = 1d0
  ELSE
   dN = daux - floor(daux)
  end if
!  iN = MAX(1,MIN(myProcess%T0_N,INT(daux+0.5d0)))

!  dN = dble(iN) - daux

   T(i) = myProcess%T0_T(iN) + dN*(myProcess%T0_T(iN+1)-myProcess%T0_T(iN))

!    WRITE(*,*) Z,daux,iN,nint(daux),dN,T(i)


 ELSE
  T(i) = myProcess%T0 + Z*myProcess%T0_slope
 END IF

 
end do

!   pause
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

!d = 1d1*dScrew                                  ! mm
d = 1d1*min(dShell,dScrew)                                  ! mm
!d = 1d1*dShell

d_factor   = 1d0-((max(min(d,d_max),d_min)-d_min)/(d_max-d_min))
tau_factor = ((max(min(tau,tau_max),tau_min)-tau_min)/(tau_max-tau_min))

WallSlip = 1d0-slip*d_factor*tau_factor
! write(*,*) WallSlip 

END FUNCTION WallSlip
