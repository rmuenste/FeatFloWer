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
INTEGER iT,iInflow
REAL*8 :: RX = 0.0d0,RY = 0.0d0,RZ = 0.0d0, RAD1 = 0.25d0
REAL*8 :: RY1 = -0.123310811d0,RY2 = 0.123310811d0,Dist1,Dist2,R_In = 0.1d0
REAL*8  dScale,XX,YY,ZZ
REAL*8 dInnerRadius,dOuterRadius,dMassFlow,dVolFlow,daux,dInnerInflowRadius,dDensity
REAL*8 DIST
REAL*8 :: PI=dATAN(1d0)*4d0
REAL*8 :: R_inflow=4d0,dNx,dNy,dNz,dNorm,dCenter(3),dNormal(3),dProfil(3)
REAL*8 :: U_bar


!  INTEGER iBCtype
!  REAL*8  massflowrate, density,outerradius,innerradius
!  REAL*8  center(3),normal(3)

ValU = 0d0
ValV = 0d0
ValW = 0d0

IF (iT.lt.0) THEN
  
 iInflow = ABS(iT)
 
 IF (ALLOCATED(myProcess%myInflow)) then
  IF (SIZE(myProcess%myInflow).le.iInflow) then
 
   IF (myProcess%myInflow(iInflow)%iBCType.eq.1) then
    dCenter       = myProcess%myInflow(iInflow)%center
    dNormal       = myProcess%myInflow(iInflow)%normal
    dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
    ddensity      = myProcess%myInflow(iInflow)%density
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
    ddensity      = myProcess%myInflow(iInflow)%density
    douterradius  = myProcess%myInflow(iInflow)%outerradius
    dinnerradius  = myProcess%myInflow(iInflow)%innerradius
    dProfil = RotDoubleParabolicVelo3D(dMassFlow,dDensity,dInnerRadius,dOuterRadius)
    ValU = dProfil(1)
    ValV = dProfil(2)
    ValW = dProfil(3)
   END IF
  else
   write(*,*) 'Inflow array is not allocated!!'
   stop
  END IF
 else
  write(*,*) 'Size of allocated inflow array is smaller then requred index:',iInflow
  stop
 END IF
 
END IF
! ! PP-Weber
! IF (iT.EQ.-1) THEN
!   ValV=RotParabolicVelo2Dy(+0.0d0,+35.2d0,-50d0,1d0,0.5d0)
! END IF

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
   CALL TransformPointToNonparallelRotAxis(0d0,0d0,0d0,XX,YY,ZZ,+1d0)
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

if(it.eq.11)then
   ValU= -PI*(Y-mySigma%a/2d0)*(myProcess%Umdr/3d1)*REAL(myProcess%iInd*myProcess%ind)
   ValV = PI*X*(myProcess%Umdr/3d1)*REAL(myProcess%iInd*myProcess%ind)
   ValW = 0d0 !myProcess%FrameVelocity
end if
if(it.eq.12)then
   ValU= -PI*(Y+mySigma%a/2d0)*(myProcess%Umdr/3d1)*REAL(myProcess%ind)
   ValV = PI*X*(myProcess%Umdr/3d1)*REAL(myProcess%ind)
   ValW = 0d0
end if
if(it.eq.13)then
  ValU =  -DBLE(myProcess%ind)*PI*Y*(myProcess%Umdr/3d1)
  ValV =   DBLE(myProcess%ind)*PI*X*(myProcess%Umdr/3d1)
  ValW =   0d0
end if


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

  dCenter=[-2.27, 0.0, -0.07]
  dNormal=[1.0, 0.0, 0.0]
  dProfil = RotParabolicVelo3D(60d0,1d0,0.250d0)
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

  dCenter=[-1.72, 0.02, 1.31]
  dNormal=[0.965, 0.0, -0.2588]
  dProfil = RotParabolicVelo3D(30d0,1d0,0.175d0)
  ValU = dProfil(1)
  ValV = dProfil(2)
  ValW = dProfil(3)

END IF
 
IF (iT.EQ.41) THEN
 ValW=RotParabolicVelo2Dz(0d0,0d0,100d0,-1d0,0.495d0)
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

! SKZ_Dietl_1
IF (iT.EQ.73) THEN
 ValU=RotParabolicVelo2Dx(0d0,12d0,-30d0,1d0,1.6d0)
END IF
! SKZ_Dietl_2
IF (iT.EQ.74) THEN
 ValW=RotParabolicVelo2Dz(0d0,0d0,-6d0,1d0,0.9d0)
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

! IDE 201810
IF (iT.EQ.88) THEN
  ValW=RotParabolicVelo2Dz(-0.0d0,+0.0d0,-102d0,1d0,2.45d0)
END IF

!Weber 201810
IF (iT.EQ.89) THEN
  ValW=RotParabolicVelo2Dz(-0.0d0,+0.0d0,-150d0,1d0,1.95d0)
END IF

! PP-Weber
IF (iT.EQ.91) THEN
  ValV=RotParabolicVelo2Dy(+0.0d0,+35.2d0,-50d0,1d0,0.5d0)
END IF
IF (iT.EQ.92) THEN
  ValW=RotParabolicVelo2Dz(+0.0d0,+0.0d0,-450d0,1d0,2.5d0)
END IF

! RAIN CARBON
IF (iT.EQ.93) THEN
  ValU=RotParabolicVelo2Dx(+0.0d0,+0.0d0,4600d0,1d0,5.0d0)
END IF

IF (iT.EQ.99) THEN
 ValW = -myFBM%ParticleNew(1)%Velocity(3)
END IF


RETURN

 CONTAINS
 REAL*8 function RotParabolicVelo2Dx(YR,ZR,DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux,YR,ZR
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO)
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux
  DIST = SQRT((Y-YR)**2d0+(Z-ZR)**2d0)
 IF (DIST.LT.dR) THEN
  RotParabolicVelo2Dx = dScale*(DR - DIST)*(DR + DIST)
 ELSE
  RotParabolicVelo2Dx = 0d0
 END IF
 END 
!------------------------------------------------------------------------------
 REAL*8 function RotParabolicVelo2Dy(XR,ZR,DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux,ZR,XR
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO)
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux
  DIST = SQRT((X-XR)**2d0+(Z-ZR)**2d0)
 IF (DIST.LT.dR) THEN
  RotParabolicVelo2Dy = dScale*(DR - DIST)*(DR + DIST)
 ELSE
  RotParabolicVelo2Dy = 0d0
 END IF
 END 
!------------------------------------------------------------------------------
 REAL*8 function RotParabolicVelo2Dz(XR,YR,DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux,YR,XR
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO)
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux
  DIST = SQRT((X-XR)**2d0+(Y-YR)**2d0)
 IF (DIST.LT.dR) THEN
  RotParabolicVelo2Dz = dScale*(DR - DIST)*(DR + DIST)
 ELSE
  RotParabolicVelo2Dz = 0d0
 END IF
 END 
!------------------------------------------------------------------------------
 function RotParabolicVelo3D(DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux
 REAL*8  dNRM
 REAL*8, dimension(3) :: RotParabolicVelo3D
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO) 
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux

  dNRM = SQRT(dNormal(1)*dNormal(1) +dNormal(2)*dNormal(2) + dNormal(3)*dNormal(3))
  dNormal(1) = dNormal(1)/dNRM
  dNormal(2) = dNormal(2)/dNRM
  dNormal(3) = dNormal(3)/dNRM
  dist = SQRT((X-(dCenter(1)))**2d0 + (Y-(dCenter(2)))**2d0 + (Z-(dCenter(3)))**2d0)

!     write(*,*) 'normal: ',dnormal
!     write(*,*) 'dcenter: ',dCenter
!     write(*,*) 'applied c/n: ',dScale,dist,dr
!     

  IF (DIST.LT.dR) THEN
   RotParabolicVelo3D(:) = dNormal(:)*dScale*(DR - DIST)*(DR + DIST)
  ELSE
   RotParabolicVelo3D(:) = 0d0
  END IF
  
  END 
!------------------------------------------------------------------------------
 function RotDoubleParabolicVelo3D(DM,DRHO,dR1,dR2)
 REAL*8  dVolFlow,DM,DRHO,dR1,dR2,daux
 REAL*8  dNRM
 REAL*8, dimension(3) :: RotDoubleParabolicVelo3D
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO) 
  daux = (PI/6d0)*(dR1+dR2)*((dR2-dR1)**3d0)
  dScale = (dVolFlow/1d0)/daux

  dNRM = SQRT(dNormal(1)*dNormal(1) +dNormal(2)*dNormal(2) + dNormal(3)*dNormal(3))
  dNormal(1) = dNormal(1)/dNRM
  dNormal(2) = dNormal(2)/dNRM
  dNormal(3) = dNormal(3)/dNRM
  dist = SQRT((X-(dCenter(1)))**2d0 + (Y-(dCenter(2)))**2d0 + (Z-(dCenter(3)))**2d0)

  IF (DIST.LT.dR2.and.DIST.GT.dR1) THEN
   RotDoubleParabolicVelo3D(:) = dNormal(:)*dScale*(DR2 - DIST)*(DIST - DR1)
  ELSE
   RotDoubleParabolicVelo3D(:) = 0d0
  END IF
  
  END 

END SUBROUTINE GetVeloBCVal
!------------------------------------------------------------
SUBROUTINE GetVeloMixerVal(X,Y,Z,ValU,ValV,ValW,iP,t)
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess
IMPLICIT NONE
REAL*8 :: myPI
INTEGER iP
REAL*8 X,Y,Z,ValU,ValV,ValW,t,yb,xb,zb

myPI = dATAN(1d0)*4d0

!  write(*,*)'vahhal'

SELECT CASE(iP)
 CASE (100) ! Y positiv
  ValU = 0d0
  ValV = 0d0
  ValW = 0d0
 CASE (101) ! Y positiv
  ValU =  -DBLE(myProcess%iInd*myProcess%ind)*myPI*(Y-mySigma%a/2d0)*(myProcess%Umdr/3d1)
  ValV =   DBLE(myProcess%iInd*myProcess%ind)*myPI*X*(myProcess%Umdr/3d1)
  ValW = 0d0
 CASE (102) ! Y negativ
  ValU =  -DBLE(myProcess%ind)*myPI*(Y+mySigma%a/2d0)*(myProcess%Umdr/3d1)
  ValV =   DBLE(myProcess%ind)*myPI*X*(myProcess%Umdr/3d1)
  ValW = 0d0
!  CASE (101) ! Y positiv
!   IF (ADJUSTL(TRIM(mySigma%RotationAxis)).NE."PARALLEL") THEN
!    CALL TransformPointToNonparallelRotAxis(x,y,z,XB,YB,ZB,-1d0)
!    ValU =  -DBLE(myProcess%iInd*myProcess%ind)*myPI*YB*(myProcess%Umdr/3d1)
!    ValV =   DBLE(myProcess%iInd*myProcess%ind)*myPI*XB*(myProcess%Umdr/3d1)
!   END IF
!  CASE (102) ! Y negativ
!   IF (ADJUSTL(TRIM(mySigma%RotationAxis)).NE."PARALLEL") THEN
!    CALL TransformPointToNonparallelRotAxis(x,y,z,XB,YB,ZB,+1d0)
!    ValU =  -DBLE(myProcess%ind)*myPI*YB*(myProcess%Umdr/3d1)
!    ValV =   DBLE(myProcess%ind)*myPI*XB*(myProcess%Umdr/3d1)
!   END IF
 CASE (103) ! Y negativ
  ValU =  -DBLE(myProcess%ind)*myPI*Y*(myProcess%Umdr/3d1)
  ValV =   DBLE(myProcess%ind)*myPI*X*(myProcess%Umdr/3d1)
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

! Bingham
IF (myRheology%Equation.EQ.6) THEN
 VNN = (1d1*myRheology%A)*aT*(1d0+myRheology%B*aT*dLimStrs)**(-myRheology%C)
END IF

dLimStrs = MIN(myRheology%ViscoMax,MAX(myRheology%ViscoMin,ABS(dStrs)))
! WRITE(*,*) dLimStrs,myRheology%eta_max,myRheology%eta_min
ViscosityModel = VNN


RETURN

end function ViscosityModel
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
