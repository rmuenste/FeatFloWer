MODULE Sigma_User
USE PP3D_MPI, ONLY:myid,showid,subnodes,dZPeriodicLength,dPeriodicity
USE var_QuadScalar ,ONLY : bNoOutflow,activeFBM_Z_Position,dTimeStepEnlargmentFactor

IMPLICIT NONE

TYPE tSegment
  INTEGER :: nOFFfilesR=0,nOFFfilesL=0,nOFFfiles=0
  CHARACTER*200, ALLOCATABLE :: OFFfilesR(:),OFFfilesL(:),OFFfiles(:)
  INTEGER, ALLOCATABLE :: idxCgal(:)
  CHARACTER*20 ObjectType,Unit
  CHARACTER*10 name
  CHARACTER*99 ::  cF
  CHARACTER*8 ::  ART
  INTEGER ::    KNETz,N
  REAL*8  :: Ds,s,delta,Dss,excentre
  REAL*8, ALLOCATABLE :: Zknet(:)
  REAL*8 :: t,D,Alpha,StartAlpha ! t=Gangsteigung
  REAL*8 :: Min, Max,L
  REAL*8 :: ZME_DiscThick,ZME_gap_SG, ZME_gap_SS 
  INTEGER :: ZME_N
  REAL*8  :: SecProf_W, SecProf_D,SecProf_L
  INTEGER :: SecProf_N, SecProf_I
  !!!!!!!!!!!!!!!!!!!!! EWIKON !!!!!!!!!!!!!!!!!!!!!
  INTEGER :: MatInd
  REAL*8 :: HeatSourceMax,HeatSourceMin,UseHeatSource
  REAL*8 :: InitTemp,Volume
  CHARACTER*200 :: TemperatureBC
   
END TYPE tSegment
!------------------------------------------------------------
TYPE tSigma
!   REAL*8 :: Dz_out,Dz_in, a, L, Ds, s, delta,SegmentLength, DZz,W
  CHARACTER cType*(50),cZwickel*(50),RotationAxis*(50)
  REAL*8 :: RotAxisCenter,RotAxisAngle
  REAL*8 :: Dz_out,Dz_in, a, L, SegmentLength, DZz,W
  REAL*8 :: SecStr_W,SecStr_D
  INTEGER :: NumberOfMat,NumberOfSeg, GANGZAHL,STLSeg=0
  TYPE (tSegment), ALLOCATABLE :: mySegment(:)
  
END TYPE tSigma

TYPE(tSigma) :: mySigma
!------------------------------------------------------------
TYPE tInflow
 INTEGER iBCtype
 REAL*8  massflowrate, density,outerradius,innerradius
 REAL*8  center(3),normal(3)
END TYPE tInflow

TYPE tProcess
   REAL*8 :: Umdr, Ta, Ti, T0=0d0, Massestrom, Dreh, Angle, dPress
   REAL*8 :: MinInflowDiameter,MaxInflowDiameter
   INTEGER :: Periodicity,Phase
   REAL*8 :: dAlpha
   CHARACTER*6 :: Rotation !RECHT, LINKS
   CHARACTER*50 :: pTYPE !RECHT, LINKS
   INTEGER :: ind,iInd
  !!!!!!!!!!!!!!!!!!!!! EWIKON !!!!!!!!!!!!!!!!!!!!!
   REAL*8 :: AmbientTemperature,HeatTransferCoeff,ConductiveGradient,ConductiveLambda
   integer   nOfInflows
   TYPE (tInflow), allocatable :: myInflow(:)
   REAL*8 :: TemperatureSensorRadius=0d0, TemperatureSensorCoor(3)=[0d0,0d0,0d0]
END TYPE tProcess

TYPE(tProcess) :: myProcess
!------------------------------------------------------------
TYPE tRheology
   INTEGER :: Equation = 5 !-->> Hogen-Powerlaw
   INTEGER :: AtFunc = 1 !-->> isotherm
   REAL*8 :: A, B, C ! Carreau Parameter
   REAL*8 :: n, K ! Power Law
   REAL*8 :: eta_max, eta_min 
   REAL*8 :: Ts, Tb, C1, C2, E! WLF Parameter
   REAL*8 :: ViscoMin = 1e-4
   REAL*8 :: ViscoMax = 1e10
END TYPE tRheology

TYPE(tRheology) :: myRheology

!------------------------------------------------------------
TYPE tThermodyn
   CHARACTER*60 :: DensityModel='NO'
   REAL*8 :: density=0d0, densitySteig=0d0
   REAL*8 :: lambda=0d0, Cp=0d0, lambdaSteig=0d0,CpSteig=0d0
   REAL*8 :: Alpha=0d0, Beta=0d0, Gamma=0d0
END TYPE tThermodyn

TYPE(tThermodyn) :: myThermodyn
TYPE(tThermodyn), Allocatable  :: myMaterials(:)
!------------------------------------------------------------

TYPE tSetup
 LOGICAL :: bPressureFBM = .FALSE.
 LOGICAL :: bAutomaticTimeStepControl = .TRUE.
 REAL*8 :: CharacteristicShearRate=1d1
 CHARACTER*200 cMeshPath
 CHARACTER*20 cMesher
 INTEGER MeshResolution,nSolutions
 INTEGER m_nT,m_nT1,m_nT2,m_nR,m_nZ,m_nP
 LOGICAL :: bGeoTest=.FALSE.,bSendEmail=.TRUE.
END TYPE tSetup
TYPE(tSetup) :: mySetup

TYPE tOutput
 INTEGER :: nOf1DLayers=16
 INTEGER :: nOfHistogramBins=16
 REAL*8 ::  HistogramShearMax=1e6,HistogramShearMin=1e-2,HistogramViscoMax=1e6,HistogramViscoMin=1e0
 REAL*8  :: CutDtata_1D=0.04d0
END TYPE tOutput

TYPE(tOutput) :: myOutput

LOGICAL bCtrl
REAL*8 DistTolerance
REAL*8 :: dMinOutputPressure = 0d0
LOGICAL :: bKTPRelease = .TRUE.

CONTAINS

! !********************GEOMETRY****************************
!
!********************GEOMETRY****************************
!
! SUBROUTINE DistanceToScrewShell(X,Y,Z,d1,d2,d3,t)
! REAL*8 X,Y,Z,t,d1,d2,d3
! REAL*8 dSeg1,dSeg2,tt
! INTEGER k,iaux
! 
! d1 = 1d-1*DistTolerance
! d2 = 1d-1*DistTolerance
! tt = t
! 
! !----------------------------------------------------------
! DO k=1, mySigma%NumberOfSeg
!  IF (mySigma%mySegment(k)%ART.EQ.'KNET' ) CALL KNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'SKNET') CALL SKNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'FOERD') CALL FOERD_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'SME'  ) CALL SME_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'ZME'  ) CALL ZME_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  d1 = MIN(dSeg1,d1)
!  d2 = MIN(dSeg2,d2)
! END DO
! !-------------------------------------------------------
! 
! CALL Shell_dist(X,Y,Z,d3)
! 
! return
! 
! END SUBROUTINE DistanceToScrewShell
! !
!********************GEOMETRY****************************
!
SUBROUTINE Shell_dist(X,Y,Z,d)
REAL*8 X,Y,Z,d
REAL*8 PX1,PX2,dZ_max,dd1,dd2,dd3

dZ_max = 0.5d0*mySigma%Dz_out



IF (Y.GT.0d0) THEN

PX1 = +SQRT(dZ_max*dZ_max - (mySigma%a/2d0)*(mySigma%a/2d0))
PX2 = -SQRT(dZ_max*dZ_max - (mySigma%a/2d0)*(mySigma%a/2d0))

dd1 = dZ_max - SQRT(X*X + (Y-mySigma%a/2d0)*(Y-mySigma%a/2d0))
dd2 =          SQRT((X-PX1)*(X-PX1) + Y*Y)
dd3 =          SQRT((X-PX2)*(X-PX2) + Y*Y)

ELSE

PX1 = +SQRT(dZ_max*dZ_max - (mySigma%a/2d0)*(mySigma%a/2d0))
PX2 = -SQRT(dZ_max*dZ_max - (mySigma%a/2d0)*(mySigma%a/2d0))

dd1 = dZ_max - SQRT(X*X + (Y+mySigma%a/2d0)*(Y+mySigma%a/2d0))
dd2 =          SQRT((X-PX1)*(X-PX1) + Y*Y)
dd3 =          SQRT((X-PX2)*(X-PX2) + Y*Y)

END IF

d = max(-DistTolerance,min(DistTolerance,dd1,dd2,dd3))

END SUBROUTINE Shell_dist
!
!********************GEOMETRY****************************
!
SUBROUTINE ZME_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
IMPLICIT NONE
REAL*8 X,Y,Z,t
REAL*8 :: dKnetMin,dKnetMax
INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
REAL*8 :: dBeta,XTT,YTT,ZTT,dist
REAL*8 :: ZMIN,ZMAX
REAL*8 dScale, daux,dist1, dist2,dEps,dCut,InnerDist,dInnerRadius
REAL*8 myPI
REAL*8 :: ZME_SegThick

dEps = mySigma%a/1d5

dInnerRadius = 0.5d0*mySigma%Dz_in
! IF (mySigma%GANGZAHL.NE.0) THEN
!  dInnerRadius = 1d0*(mySigma%A - 0.5d0*mySigma%Ds - 1d0*mySigma%s)
! ELSE
!  dInnerRadius = 1d0*0.5d0*mySigma%Ds
! END IF

ZME_SegThick = 2d0*(mySigma%mySegment(iSeg)%ZME_DiscThick+mySigma%mySegment(iSeg)%ZME_gap_SS)

myPI = dATAN(1d0)*4d0

IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
 DO l=1, mySigma%mySegment(iSeg)%ZME_N
  dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick-dEps
  dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*ZME_SegThick+dEps
  IF (Z.GE.dKnetMin.AND.Z.LT.dKnetMax) THEN
   lKnet = l
   EXIT
  END IF
 END DO
ELSE
 IF (Z.lt.mySigma%mySegment(iSeg)%Min) lknet = 1
 IF (Z.gt.mySigma%mySegment(iSeg)%Max) lknet = mySigma%mySegment(iSeg)%ZME_N
END IF

! First screw
dist1 = DistTolerance
XB = X
YB = Y-mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = 0d0 - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

! daux = -5d0
DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%ZME_N,lKnet+1)
 
ZMIN = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick + 0.5d0*mySigma%mySegment(iSeg)%ZME_gap_SS + 0d0*mySigma%mySegment(iSeg)%ZME_DiscThick
ZMAX = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick + 0.5d0*mySigma%mySegment(iSeg)%ZME_gap_SS + 1d0*mySigma%mySegment(iSeg)%ZME_DiscThick

CALL ZME_disc(XT,YT,ZT,ZMIN,ZMAX,iSeg,daux) 
InnerDist = sqrt(xt*xt+yt*yt)-dInnerRadius
dist1 =MIN(daux,dist1,InnerDist)

END DO

CALL SecondaryProfile(XT,YT,ZT,dCut,iSeg)
dist1 = max(dist1,-dCut,-DistTolerance)

IF (dist1.LT.0d0) THEN
 inpr = 101
END IF

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Second screw
dist2 = DistTolerance
XB = X
YB = Y+mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = 0d0 - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

! daux = -5d0
DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%ZME_N,lKnet+1)
! l=1

ZMIN = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick + 1.5d0*mySigma%mySegment(iSeg)%ZME_gap_SS + 1d0*mySigma%mySegment(iSeg)%ZME_DiscThick
ZMAX = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick + 1.5d0*mySigma%mySegment(iSeg)%ZME_gap_SS + 2d0*mySigma%mySegment(iSeg)%ZME_DiscThick

CALL ZME_disc(XT,YT,ZT,ZMIN,ZMAX,iSeg,daux) 
InnerDist = sqrt(xt*xt+yt*yt)-dInnerRadius

! dist2 = max(-DistTolerance,min(DistTolerance,daux,dist2,InnerDist))
 dist2 =MIN(daux,dist2,InnerDist)

END DO

CALL SecondaryProfile(XT,YT,ZT,dCut,iSeg)
dist2 = max(dist2,-dCut,-DistTolerance)
! dist2 = max(dist2,-dCut)

IF (dist2.LT.0d0) THEN
 inpr = 102
END IF

END SUBROUTINE ZME_elem
!
!********************GEOMETRY****************************
!

SUBROUTINE SKNET_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
IMPLICIT NONE
REAL*8 X,Y,Z,t
REAL*8 :: dKnetMin,dKnetMax,dKnetMed
INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
REAL*8 :: dBeta1,dBeta2,XTT,YTT,ZTT,dist
REAL*8 :: dZ
REAL*8 daux1,daux2
REAL*8 dScale, daux,dist1, dist2,dEps,dCut
REAL*8 myPI

dEps = mySigma%a/1d5

myPI = dATAN(1d0)*4d0

IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
 DO l=1, mySigma%mySegment(iSeg)%N
  dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
  dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
  IF (Z.GE.dKnetMin.AND.Z.LT.dKnetMax) THEN
   lKnet = l
   EXIT
  END IF
 END DO
ELSE
 IF (Z.lt.mySigma%mySegment(iSeg)%Min) lknet = 1
 IF (Z.gt.mySigma%mySegment(iSeg)%Max) lknet = mySigma%mySegment(iSeg)%N
END IF

! First screw
dist1 = DistTolerance
XB = X
YB = Y-mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)

dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
dKnetMed = 0.5d0*(dKnetMin+dKnetMax)
dBeta1 = myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2
dBeta2 = myPI*DBLE(l-2)*mySigma%mySegment(iSeg)%Alpha/1.8d2

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta1) - YT*sin(dBeta1)
YTT = XT*sin(dBeta1) + YT*cos(dBeta1)

CALL Schnecke_nGang(XTT,YTT,iSeg,daux1)

XTT = XT*cos(dBeta2) - YT*sin(dBeta2)
YTT = XT*sin(dBeta2) + YT*cos(dBeta2)

CALL Schnecke_nGang(XTT,YTT,iSeg,daux2)

dZ = MIN(ABS(ZT-dKnetMed),ABS(ZT-dKnetMax))
IF (ZT.GE.dKnetMed.and.ZT.LE.dKnetMax) THEN
 IF (daux1.ge.0d0) THEN
  daux=daux1
 ELSE
  daux=Max(daux1,-dZ)
 END IF
END IF

IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMed) THEN
 daux=MAX(daux1,daux2)
 IF (daux.le.0d0) THEN
   dZ = MIN(ABS(ZT-dKnetMin),ABS(ZT-dKnetMed))
   daux=max(daux,-dZ)
 ELSE
  dZ = ABS(ZT-dKnetMed)
  IF (daux1.lt.0d0) THEN
   daux=min(dZ,daux)
  ELSE
   daux=min(daux,SQRT(daux1**2d0+dZ**2d0))
  END IF
 END IF
END IF

IF (ZT.GT.dKnetMax) THEN
 dZ = ABS(ZT-dKnetMax)
 daux=daux1
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

IF (ZT.LT.dKnetMin) THEN
 dZ = ABS(ZT-dKnetMin)
 daux=MAX(daux1,daux2)
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

dist1 = max(-DistTolerance,min(daux,dist1))

END DO

IF (dist1.LT.0d0) THEN
 inpr = 101
END IF


! Second screw
dist2 = DistTolerance
XB = X
YB = Y+mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind

! First the point needs to be transformed back to time = 0
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)

dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
dKnetMed = 0.5d0*(dKnetMin+dKnetMax)
dBeta1 = myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2
dBeta2 = myPI*DBLE(l-2)*mySigma%mySegment(iSeg)%Alpha/1.8d2

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta1) - YT*sin(dBeta1)
YTT = XT*sin(dBeta1) + YT*cos(dBeta1)

CALL Schnecke_nGang(XTT,YTT,iSeg,daux1)

XTT = XT*cos(dBeta2) - YT*sin(dBeta2)
YTT = XT*sin(dBeta2) + YT*cos(dBeta2)

CALL Schnecke_nGang(XTT,YTT,iSeg,daux2)

dZ = MIN(ABS(ZT-dKnetMed),ABS(ZT-dKnetMax))
IF (ZT.GE.dKnetMed.and.ZT.LE.dKnetMax) THEN
 IF (daux1.ge.0d0) THEN
  daux=daux1
 ELSE
  daux=Max(daux1,-dZ)
 END IF
END IF

IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMed) THEN
 daux=MAX(daux1,daux2)
 IF (daux.le.0d0) THEN
   dZ = MIN(ABS(ZT-dKnetMin),ABS(ZT-dKnetMed))
   daux=max(daux,-dZ)
 ELSE
  dZ = ABS(ZT-dKnetMed)
  IF (daux1.lt.0d0) THEN
   daux=min(dZ,daux)
  ELSE
   daux=min(daux,SQRT(daux1**2d0+dZ**2d0))
  END IF
 END IF
END IF

IF (ZT.GT.dKnetMax) THEN
 dZ = ABS(ZT-dKnetMax)
 daux=daux1
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

IF (ZT.LT.dKnetMin) THEN
 dZ = ABS(ZT-dKnetMin)
 daux=MAX(daux1,daux2)
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

dist2 = max(-DistTolerance,min(daux,dist2))
! dist2 =MIN(daux,dist2)

END DO

IF (dist2.LT.0d0) THEN
 inpr = 102
END IF

END SUBROUTINE SKNET_elem
!
!********************GEOMETRY****************************
!
SUBROUTINE SME_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
IMPLICIT NONE
REAL*8 X,Y,Z,t
REAL*8 :: dKnetMin,dKnetMax
INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
REAL*8 :: dBeta,XTT,YTT,ZTT,dist
REAL*8 :: dZ
REAL*8 dScale, daux,dist1, dist2,dEps,dCut
REAL*8 myPI

dEps = mySigma%a/1d5

myPI = dATAN(1d0)*4d0

IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
 dBeta=2d0*(Z-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
ELSE
 IF (Z.lt.mySigma%mySegment(iSeg)%Min) dBeta=2d0*(mySigma%mySegment(iSeg)%min-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
 IF (Z.gt.mySigma%mySegment(iSeg)%Max) dBeta=2d0*(mySigma%mySegment(iSeg)%max-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
END IF

! First screw
! dist1 = 5d0
XB = X
YB = Y-mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

CALL Schnecke_nGang(XTT,YTT,iSeg,daux)
CALL SecondaryProfile(XT,YT,ZT,dCut,iSeg)
daux = max(daux,-dCut)

dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
 IF (daux.ge.0d0) THEN
  daux=daux
 ELSE
  daux=Max(daux,-dZ)
 END IF
ELSE
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

dist1 = max(-DistTolerance,min(DistTolerance,daux))
! dist1 =MIN(daux,dist1)

IF (dist1.LT.0d0) THEN
 inpr = 101
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First screw
! dist2 = 5d0
XB = X
YB = Y+mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind

! First the point needs to be transformed back to time = 0
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

CALL Schnecke_nGang(XTT,YTT,iSeg,daux)
CALL SecondaryProfile(XT,YT,ZT,dCut,iSeg)
daux = max(daux,-dCut)

dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
 IF (daux.ge.0d0) THEN
  daux=daux
 ELSE
  daux=Max(daux,-dZ)
 END IF
ELSE
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

dist2 = max(-DistTolerance,min(DistTolerance,daux))
! dist2 =MIN(daux,dist2)

IF (dist2.LT.0d0) THEN
 inpr = 102
END IF

END SUBROUTINE SME_elem
!
!********************GEOMETRY****************************
!
SUBROUTINE FOERD_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
IMPLICIT NONE
REAL*8 X,Y,Z,t
REAL*8 :: dKnetMin,dKnetMax
INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
REAL*8 :: dBeta,XTT,YTT,ZTT,dist
REAL*8 :: dZ
REAL*8 dScale, daux,dist1, dist2,dEps,dCut
REAL*8 myPI

dEps = mySigma%a/1d5

myPI = dATAN(1d0)*4d0

IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
 dBeta=2d0*(Z-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
ELSE
 IF (Z.lt.mySigma%mySegment(iSeg)%Min) dBeta=2d0*(mySigma%mySegment(iSeg)%min-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
 IF (Z.gt.mySigma%mySegment(iSeg)%Max) dBeta=2d0*(mySigma%mySegment(iSeg)%max-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
END IF

! First screw
! dist1 = DistTolerance
XB = X
YB = Y-mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

CALL Schnecke_nGang(XTT,YTT,iSeg,daux)

dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
 IF (daux.ge.0d0) THEN
  daux=daux
 ELSE
  daux=Max(daux,-dZ)
 END IF
ELSE
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

! dist1 =MIN(daux,dist1)
dist1 = max(-DistTolerance,min(DistTolerance,daux))

IF (dist1.LT.0d0) THEN
 inpr = 101
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First screw
! dist2 = 5d0
XB = X
YB = Y+mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind

! First the point needs to be transformed back to time = 0
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

CALL Schnecke_nGang(XTT,YTT,iSeg,daux)

dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
 IF (daux.ge.0d0) THEN
  daux=daux
 ELSE
  daux=Max(daux,-dZ)
 END IF
ELSE
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

! dist2 =MIN(daux,dist2)
dist2 = max(-DistTolerance,min(DistTolerance,daux))

IF (dist2.LT.0d0) THEN
 inpr = 102
END IF

END SUBROUTINE FOERD_elem
!
!********************GEOMETRY****************************
!
SUBROUTINE KNET_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
IMPLICIT NONE
REAL*8 X,Y,Z,t
REAL*8 :: dKnetMin,dKnetMax
INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
REAL*8 :: dBeta,XTT,YTT,ZTT,dist
REAL*8 :: dZ
REAL*8 dScale, daux,dist1, dist2,dEps,dCut
REAL*8 myPI

dEps = mySigma%a/1d5
myPI = dATAN(1d0)*4d0

IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
 DO l=1, mySigma%mySegment(iSeg)%N
  dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
  dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
  IF (Z.GE.dKnetMin.AND.Z.LT.dKnetMax) THEN
   lKnet = l
   EXIT
  END IF
 END DO
ELSE
 IF (Z.lt.mySigma%mySegment(iSeg)%Min) lknet = 1
 IF (Z.gt.mySigma%mySegment(iSeg)%Max) lknet = mySigma%mySegment(iSeg)%N
END IF

! First screw
dist1 = DistTolerance
XB = X
YB = Y-mySigma%a/2d0
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)

dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
dBeta =  myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

CALL Schnecke_nGang(XTT,YTT,iSeg,daux)

dZ = MIN(ABS(ZT-dKnetMin),ABS(ZT-dKnetMax))
IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMax) THEN
 IF (daux.ge.0d0) THEN
  daux=daux
 ELSE
  daux=Max(daux,-dZ)
 END IF
ELSE
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

dist1 = max(-DistTolerance,min(daux,dist1))

END DO

IF (dist1.LT.0d0) THEN
 inpr = 101
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Second screw
dist2 = DistTolerance
XB = X
YB = Y+mySigma%a/2d0 
ZB = Z

! First the point needs to be transformed back to time = 0
IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind

! First the point needs to be transformed back to time = 0
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)

dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
dBeta =  myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

CALL Schnecke_nGang(XTT,YTT,iSeg,daux)

 dZ = MIN(ABS(ZT-dKnetMin),ABS(ZT-dKnetMax))
IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMax) THEN
 IF (daux.ge.0d0) THEN
  daux=daux
 ELSE
  daux=max(daux,-dZ)
 END IF
ELSE
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

dist2 = max(-DistTolerance,min(daux,dist2))

END DO

IF (dist2.LT.0d0) THEN
 inpr = 102
END IF

END SUBROUTINE KNET_elem
!
!********************GEOMETRY****************************
!
SUBROUTINE EKNET_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
IMPLICIT NONE
REAL*8 X,Y,Z,t
REAL*8 :: dKnetMin,dKnetMax
INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
REAL*8 :: dBeta,XTT,YTT,ZTT,dist
REAL*8 :: dZ, exc
REAL*8 dScale, daux,dist1, dist2,dEps,dCut
REAL*8 myPI

dEps = mySigma%a/1d5
exc = mySigma%mySegment(iSeg)%excentre
myPI = dATAN(1d0)*4d0

IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
 DO l=1, mySigma%mySegment(iSeg)%N
  dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
  dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
  IF (Z.GE.dKnetMin.AND.Z.LT.dKnetMax) THEN
   lKnet = l
   EXIT
  END IF
 END DO
ELSE
 IF (Z.lt.mySigma%mySegment(iSeg)%Min) lknet = 1
 IF (Z.gt.mySigma%mySegment(iSeg)%Max) lknet = mySigma%mySegment(iSeg)%N
END IF

! First screw
dist1 = DistTolerance
XB = X
YB = Y-mySigma%a/2d0 !+exc
ZB = Z

! First the point needs to be transformed back to time = 0
dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

YT = YT + exc

DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)

dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
dBeta =  myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

CALL Schnecke_nGang(XTT,YTT,iSeg,daux)

dZ = MIN(ABS(ZT-dKnetMin),ABS(ZT-dKnetMax))
IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMax) THEN
 IF (daux.ge.0d0) THEN
  daux=daux
 ELSE
  daux=Max(daux,-dZ)
 END IF
ELSE
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

dist1 = max(-DistTolerance,min(daux,dist1))

END DO

IF (dist1.LT.0d0) THEN
 inpr = 101
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Second screw
dist2 = DistTolerance
XB = X
YB = Y+mySigma%a/2d0 !+exc
ZB = Z

! First the point needs to be transformed back to time = 0
IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind

! First the point needs to be transformed back to time = 0
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

YT = YT + exc

DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)

dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
dBeta =  myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2

! Next the point needs to be transformed back to Z = 0 level
XTT = XT*cos(dBeta) - YT*sin(dBeta)
YTT = XT*sin(dBeta) + YT*cos(dBeta)
ZTT = ZT

CALL Schnecke_nGang(XTT,YTT,iSeg,daux)

 dZ = MIN(ABS(ZT-dKnetMin),ABS(ZT-dKnetMax))
IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMax) THEN
 IF (daux.ge.0d0) THEN
  daux=daux
 ELSE
  daux=max(daux,-dZ)
 END IF
ELSE
 IF (daux.le.0d0) THEN
  daux=dZ
 ELSE
  daux=SQRT(daux**2d0+dZ**2d0)
 END IF
END IF

dist2 = max(-DistTolerance,min(daux,dist2))

END DO

IF (dist2.LT.0d0) THEN
 inpr = 102
END IF

END SUBROUTINE EKNET_elem
!
!-----------------------------------------------------------
!
SUBROUTINE ZME_disc(XP,YP,ZP,ZPMIN,ZPMAX,iSeg,Dif)
REAL*8 XP,YP,ZP, Dif,ZPMIN,ZPMAX
REAL*8  Dz, a, delta, R_out
REAL*8 d1,d2,d3 !,dShift
INTEGER iSeg
! INTEGER iSide

Dz=mySigma%Dz_out
a=mySigma%a
delta=mySigma%mySegment(iSeg)%ZME_gap_SG

R_out = Dz/2D0 - delta

! IF (iSide.EQ.1) dShift = -10.0d0
! IF (iSide.EQ.2) dShift = +10.0d0

d1 = dSQRT(XP*XP + YP*YP)-R_out
d2 = ZP - (ZPMAX)
d3 = (ZPMIN) - ZP

Dif = MAX(d1,d2,d3)

END SUBROUTINE ZME_disc
!
!-----------------------------------------------------------
!
SUBROUTINE Schnecke_nGang(XP,YP,iSeg,Dif)
!Dz= Zylinderdurchmesser, a=Achsabstand
!delta= Spiel Schneckenkamm-Gehaeusewand
!s= Spiel Schnecke-Schnecke (from Sigma!)
INTEGER iSeg
REAL*8  Dz, a, delta, s, Ds, exc    
REAL*8 dTheta, XP, YP, Dif, XB, YB
REAL*8 R1, R2, R3, dAlpha,dGamma,dCutAngle,dRadius,dN,dCrit1,dCrit2
REAL*8  :: myPI = dATAN(1d0)*4d0
INTEGER iAux

dN=DBLE(mySigma%GANGZAHL)
Dz=mySigma%Dz_out
a=mySigma%a
delta=mySigma%mySegment(iSeg)%delta
s=mySigma%mySegment(iSeg)%s

IF (mySigma%mySegment(iSeg)%ART .EQ. "EKNET") THEN
   exc=mySigma%mySegment(iSeg)%excentre
ELSE
   exc=0d0
END IF

Dz= Dz -2d0*exc
Ds= Dz - 2d0*delta

R2 = Dz/2D0 - delta
R1 = a-s-R2
R3 = a-s

dRadius = SQRT(XP*XP + YP*YP)
IF (XP.LT.0d0) THEN 
 dGamma  = 0d0
ELSE
 dGamma = myPi
END IF

dGamma = dGamma + dATAN(YP/XP)
dCutAngle = 2d0*myPi/(dN)
iAux = FLOOR(dGamma/dCutAngle)
dGamma = dGamma - DBLE(iAux)*dCutAngle - 0.5d0*dCutAngle
dGamma = ABS(dGamma)
XB = (dRadius*COS(dGamma))
YB = (dRadius*SIN(dGamma))

dTheta = 0.5d0*(myPI/dN - 2d0*ACOS(R3/(2d0*R2)))
!dTheta = 0.5d0*(myPI/dN - 2d0*ACOS((a-s)/(mySigma%Dz_out-2d0*exc-2d0*delta)))

IF (delta.lt.0) THEN
 if (myid.eq.1) then
  WRITE(*,*) "the screw could not be created : wrong delta"
  WRITE(*,*) "theta [deg]", 180d0*dtheta/myPi
  WRITE(*,*) "Ds: ", Ds
  WRITE(*,*) "delta: ", delta
  WRITE(*,*) "s: ", s
  WRITE(*,*) "exc: ", exc
  WRITE(*,*) "a: ", a
  WRITE(*,*) "Dz,Dz*: ", mySigma%Dz_out,Dz
  WRITE(*,*) "R2,R3: ", R2,R3
!   WRITE(*,*) "crit1", (mySigma%Dz_out-2d0*exc-2d0*delta)/(a-s)
!   WRITE(*,*) "crit2", 1d0/cos(mypi/(2d0*dn))
 END IF
 stop
END IF

IF (dTheta.lt.0d0) THEN
 if (myid.eq.1) then
  WRITE(*,*) "the screw could not be created : wrong angle"
  WRITE(*,*) "theta [deg]", 180d0*dtheta/myPi
  WRITE(*,*) "Ds: ", Ds
  WRITE(*,*) "delta: ", delta
  WRITE(*,*) "s: ", s
  WRITE(*,*) "exc: ", exc
  WRITE(*,*) "a: ", a
  WRITE(*,*) "Dz,Dz*: ", mySigma%Dz_out,Dz
  WRITE(*,*) "R2,R3: ", R2,R3
!   WRITE(*,*) "crit1", (mySigma%Dz_out-2d0*exc-2d0*delta)/(a-s)
!   WRITE(*,*) "crit2", 1d0/cos(mypi/(2d0*dn))
 END IF
 stop
ENDIF

Dif = dSQRT(XB*XB+YB*YB)-R2

IF ( dGamma.LE.dTheta .AND. dGamma.GE.0d0 ) THEN
  Dif = dSQRT(XB*XB+YB*YB)-R1
END IF

IF ( dGamma.GT.dTheta .AND. dGamma.LT.(myPI/dN-dTheta)) THEN
  Dif = dSQRT((XB + R2*dCOS(dTheta))**2d0 +(YB + R2*dSIN(dTheta))**2d0)-R3
END IF




END SUBROUTINE Schnecke_nGang
!------------------------------------------------------------
SUBROUTINE SecondaryProfile(XB,YB,ZB,d,k)
REAL*8 XB,YB,ZB,d, XT,YT,XP,YP
INTEGER k
REAL*8 dGamma,dRadius,dCutAngle,dCutRadius,dBeta,dStartAngle
! REAL*8  :: dCutDepth=15d0,dCutWidth=28d0,dCutPitch !,dCutAngleF=0.33333333333333d0 !2d0*0.66666666667d0
! INTEGER :: nCut = 4,SSE_Type = 3
REAL*8  :: dCutDepth,dCutWidth,dCutPitch !,dCutAngleF=0.33333333333333d0 !2d0*0.66666666667d0
INTEGER :: nCut,SSE_Type
INTEGER iAux
REAL*8  :: myPI = dATAN(1d0)*4d0

nCut      = mySigma%mySegment(k)%SecProf_N
SSE_Type  = mySigma%mySegment(k)%SecProf_I
dCutDepth = mySigma%mySegment(k)%SecProf_D
dCutWidth = mySigma%mySegment(k)%SecProf_W
dCutPitch = mySigma%mySegment(k)%SecProf_L

dCutRadius = 0.5d0*mySigma%mySegment(k)%Ds !mySigma%Dz_out/2D0 !- mySigma%delta

! if (k.eq.1) dStartAngle =  myPi* 0d0/180d0
! if (k.eq.2) dStartAngle = 2d0*dCutAngleF*myPI/mySigma%mySegment(1)%t !+myPi*45d0/180d0

! XT = XB
! YT = YB
! dBeta= -2d0*dCutAngleF*(ZB-mySigma%mySegment(k)%Min)*myPI/40d0 !mySigma%mySegment(k)%t

dBeta=2d0*(ZB-mySigma%mySegment(k)%min)*myPI/dCutPitch
XT = XB*cos(dBeta) - YB*sin(dBeta)
YT = XB*sin(dBeta) + YB*cos(dBeta)

dRadius = SQRT(XT*XT + YT*YT)
IF (XT.LT.0d0) THEN 
 dGamma  = 0d0
ELSE
 dGamma = myPi
END IF

dGamma = dGamma + dATAN(YT/XT)
dCutAngle = 2d0*myPi/DBLE(nCut)
iAux = FLOOR(dGamma/dCutAngle)
dGamma = dGamma - DBLE(iAux)*dCutAngle - 0.5d0*dCutAngle
XB = (dRadius*COS(dGamma))
YB = (dRadius*SIN(dGamma))

IF (SSE_Type.eq.1) THEN
 d = SQRT((XB-dCutRadius)*(XB-dCutRadius) + YB*YB) - dCutDepth
END IF

IF (SSE_Type.eq.2) THEN
 d = MAX(-dCutDepth-(XB-dCutRadius),-0.5d0*dCutWidth-YB,-0.5d0*dCutWidth+YB) !XB-dCutRadius+dCutDepth !MAX(XB-dCutRadius+dCutDepth, YB-0.5d0*dCutWidth,0.5d0*dCutWidth-YB)
ENDIF

IF (SSE_Type.eq.3) THEN
 d = FindDist(YB,XB-dCutRadius)
ENDIF

! IF (SSE_Type.eq.2) THEN
!  d = SQRT(dCutWidth*dCutWidth*(XB-dCutRadius)*(XB-dCutRadius) + dCutDepth*dCutDepth*YB*YB)-dCutWidth*dCutDepth
!  d=0.1d0*d
! END IF

 CONTAINS

 FUNCTION FindDist(ppx, ppy)
 REAL*8  ppx, ppy
 REAL*8  FindDist
 REAL*4  px, py, pz, ax, ay, az, bx, by, bz, tol
 INTEGER np,noptfd
 REAL*4  dpmin, fdmin, cx, cy, cz
 INTEGER nlim, itrun, nerr
 INTEGER i
 REAL*8 ::  A,B,C,dSGN,dF,dDist
 REAL*8 :: memory(2,33)

 C=-dCutDepth
 A= -C/(dCutWidth*dCutWidth)*4d0
 B= 0d0

 tol = 1e-6
 np =1
 noptfd = 0

 px = REAL(ppx)
 py = REAL(ppy)
 pz = 0e0

 dF = A*px*px + B*px + C- py
! dpmin = dF
 dSGN = (dF)/ABS(dF)

 ax =  0e0
 ay =  REAL(C)
 az =  0e0

 bx = REAL(ppx)
 by = REAL(A*px*px + B*px + C)
 bz = 0e0

 DO i=1,33
  cx = 0.5d0*(ax + bx)
  cy = REAL(A*cx*cx + B*cx + C)

  memory(:,i) = [cx,cy]

  IF (ABS(ay).lt.ABS(by)) THEN
   bx = cx
   by = cy
  ELSE
   ax = cx
   ay = cy
  END IF

  IF (((ax-bx)**2d0+(ay-by)**2d0).lt. 1d-4*mySigma%Dz_out) EXIT

 END DO

 dDist = SQRT((px-cx)**2d0+(py-cy)**2d0)

 FindDist = dSGN*ABS(dDist) !max(-dF,ppx-mySigma%Dz_out/2D0)

 END FUNCTION FindDist

END SUBROUTINE SecondaryProfile
!
!-----------------------------------------------------------
!
!  End of geometry
!-----------------------------------------------------------
SUBROUTINE Sigma_AdjustTimeParameters(dt,tmax,dtOut)

REAL*8 dt,tmax,dtOut
REAL*8 :: tFact = 10d0

dt    = 1d0/(myProcess%Umdr/6d1)/(tFact*1d3)
dtOut = tFact*24.9999d0*dt
tmax  = myProcess%Dreh/(myProcess%Umdr/6d1)

! IF (bGeoTest) dt = 125*dt

END SUBROUTINE Sigma_AdjustTimeParameters

END MODULE Sigma_User
