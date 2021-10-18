module types
!-------------------------------------------------------------------------------------------------
! A module that contains several variants of the Laplacian 
! smoother 'Umbrella'. Besides the standard version of this 
! smoother type, versions with user-defined weighting functions
! are available
!-------------------------------------------------------------------------------------------------
! No implicit variables in this module
implicit none

!================================================================================================
!  A structure for a bounding box 
!================================================================================================
type tBoundingBox
  double precision, dimension(3,2)  :: vertices
end type tBoundingBox

!================================================================================================
!  A structure that maintains multiple bounding boxes 
!================================================================================================
type tmgBoundingBox
  type(tBoundingBox), dimension(:), allocatable :: bb
end type tmgBoundingBox

type tPostprocessingParams
  real*8 :: U_mean = 0.2d0
  ! This is the setting for 2D FAC, for
  ! the full 3D FAC the value is H = 0.205d0
  real*8 :: H = 0.05d0
  real*8 :: D = 0.1d0
  real*8 :: Sc_U = 1d0
  real*8 :: Sc_Mu = 1d0
  real*8 :: Sc_a = 1d0
end type tPostprocessingParams

!================================================================================================
!  A structure for a single mesh instance 
!================================================================================================
TYPE tMesh
  ! Mesh integer parameters
  integer :: NEL = 0
  integer :: NVT = 0
  integer :: NET = 0
  integer :: NAT = 0
  integer :: NBCT = 0
  integer :: NVBD = 0
  integer :: NEBD = 0
  integer :: NABD = 0

  integer :: NVE = 8
  integer :: NEE = 12
  integer :: NAE = 6

  integer :: NVEL = 0

  integer :: NNelAtVertex = 0

  ! Conntectivity arrays

  ! A pointer to the mesh coordinates, in a mesh hierarchy this array
  ! will be filled with the actual coordinates only on the finest
  ! level of the hierarchical mesh while the lower levels
  ! contain only pointers to the array on the finest level
  real*8,  pointer, dimension(:,:) :: dcorvg => null()

  real*8,  allocatable, dimension(:,:) :: dcorag

  real*8,  allocatable, dimension(:) :: dvol

  ! Coordinate arrays, _old is for mesh deformation 
  real*8,  allocatable, dimension(:,:) :: dcorvg_old

  ! The vertices at an element
  integer, allocatable, dimension(:,:) :: kvert

  ! The elements at a vertex 
  integer, allocatable, dimension(:,:) :: kvel

  ! The edges at a vertex
  integer, allocatable, dimension(:,:) :: kved

  ! 
  integer, allocatable, dimension(:,:) :: kvar

  ! The edge array with its vertex indices
  integer, allocatable, dimension(:,:) :: kedge

  ! The neighbors at an element
  integer, allocatable, dimension(:,:) :: kadj

  integer, allocatable, dimension(:,:) :: kadj2

  ! The vertices at a face
  integer, allocatable, dimension(:,:) :: karea

  integer, allocatable, dimension(:) :: knpr

  ! The elements at a vertex 
  ! elementsAtVertexIdx=array [1..NVT+1] of integer.
  ! Index array for elementsAtVertex of length NVT+1 for describing the
  ! elements adjacent to a corner vertex. for vertex IVT, the array
  ! elementsAtVertex contains the numbers of the elements around this
  ! vertex at indices
  !     elementsAtVertexIdx(IVT)..elementsAtVertexIdx(IVT+1)-1.
  ! By subtracting
  !     elementsAtVertexIdx(IVT+1)-elementsAtVertexIdx(IVT)
  ! One can get the number of elements adjacent to a vertex IVT.
  integer, allocatable, dimension(:) :: elementsAtVertexIdx

  ! Array containing the Elements Adjacent to a Vertex.
  ! Handle to
  !       elementsAtVertex = array(1..*) of integer
  ! elementsAtVertex ( elementsAtVertexIdx(IVT)..elementsAtVertexIdx(IVT+1)-1 )
  ! contains the number of the adjacent element in a vertex.
  ! This replaces the old KVEL array.
  integer, allocatable, dimension(:) :: elementsAtVertex


END TYPE

TYPE tBndryNone
 LOGICAL :: bOuterPoint=.FALSE.
 LOGICAL :: ParamTypes(4)
 INTEGER :: nPoint=0,nLine=0,nSurf=0,nVolume=0
 INTEGER, ALLOCATABLE :: P(:),L(:),S(:),V(:)
 real*8    :: x(3)
 logical   :: bx=.false.
END TYPE tBndryNone

!================================================================================================
!  A structure for a mesh hierarchy 
!================================================================================================
TYPE tMultiMesh

  integer :: nlmin

  ! maximum level used for computation
  integer :: nlmax

  ! this is usually nlmax + 1
  integer :: maxlevel
  
  type(tBndryNone), allocatable, dimension(:) :: BndryNodes

  type(tMesh), allocatable, dimension(:) :: level

END TYPE

TYPE tParticle
 REAL*8 time
 REAL*8 coor(3),velo(3)
 INTEGER indice
 INTEGER :: ID=1
END TYPE tParticle

! Careful this is a type from PLinSc
!TYPE tParam
! REAL*8  thetaRI,AvgGradPhiTol
! REAL*8  defCrit,epsCrit,MinDef,theta
! INTEGER NLminRI,NLmaxRI,nRI
! INTEGER NLmin,NLmax
! INTEGER SolvIter,SolvType,iMass
! INTEGER StabScheme,SrfCubF,VolCubF
! LOGICAL SlopeLimit
!END TYPE tParam

TYPE tMGParamIn
 INTEGER MinLev,MedLev,MaxLev,MaxDifLev,MinIterCycle,MaxIterCycle,nSmootherSteps,nIterCoarse,CrsSolverType,SmootherType
 integer :: vanka
 REAL*8  DefImprCoarse,Criterion1,Criterion2,RLX
 REAL*8 :: CrsRelaxPrm=2d0/3d0,CrsRelaxParPrm=1d0/3d0
 CHARACTER*1 CycleType
 CHARACTER*10 MGProlongation
END TYPE tMGParamIn

TYPE tMGParamOut
 REAL*8 DefInitial,DefFinal
 REAL*8 RhoMG1,RhoMG2
 REAL*8 u_rel(6),p_rel(2)
 INTEGER UsedIterCycle,nIterCoarse
END TYPE tMGParamOut

TYPE tParam
 Character*20 :: cEquation
 Character*20, allocatable :: cField(:)
 integer :: nOfFields

 REAL*8  defCrit,epsCrit,MinDef

 INTEGER NLmin,NLmax

 INTEGER SolvIter,SolvType

 LOGICAL AFC
 TYPE(tMGParamIn) :: MGprmIn
 TYPE(tMGParamOut):: MGprmOut(3)
END TYPE tParam


TYPE mg_vector
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: x
END TYPE mg_vector

TYPE tAFC
 INTEGER :: iedge,nedge,nu
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: inod,jnod,iaux,isep
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: aedge
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: pp,pm,qm,qp
END TYPE
TYPE(tAFC) AFC

TYPE TlMatrix
 INTEGER :: nu,na
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: ColA,LdA
END TYPE


! Careful this type is from PLinScalar and
! has the same name as the one from LinScalar
!TYPE lScalar
! CHARACTER cName*7
! INTEGER :: ndof,na,nel,iNumFace
! REAL*8  , DIMENSION(:)  , ALLOCATABLE :: aux,rhs,def,val_old,val_ad
! REAL*8  , DIMENSION(:)  , ALLOCATABLE :: Q1
! INTEGER , DIMENSION(:,:), ALLOCATABLE :: iParFace
! REAL*8  , DIMENSION(:,:), ALLOCATABLE :: dParFace,dParMidC
! TYPE(mg_vector), DIMENSION(:),ALLOCATABLE :: val
! TYPE(tParam) :: prm
!END TYPE

TYPE lScalar
 CHARACTER cName*7
 INTEGER :: ndof,na
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: knpr
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: aux,rhs,def,val_old,oldSol
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: src,snk
 TYPE(mg_vector), DIMENSION(:),ALLOCATABLE :: val
 TYPE(tParam) :: prm
END TYPE

TYPE tStatistics
 INTEGER :: iNonLin=0,iLinUVW=0,iLinP=0
 REAL  :: tMGUVW=0d0,tMGP=0d0,tDefUVW=0d0,tDefP=0d0,tCorrUVWP=0d0
 REAL  :: tGMVOut=0d0,tDumpOut=0d0
 REAL  :: tSmat=0d0,tKmat=0d0,tDmat=0d0,tMmat=0d0,tCmat=0d0
 REAL  :: tRestUVW=0d0,tProlUVW=0d0,tSmthUVW=0d0,tSolvUVW=0d0
 REAL  :: tRestP=0d0,tProlP=0d0,tSmthP=0d0,tSolvP=0d0
 REAL  :: tCommV = 0d0
 REAL  :: tCommP = 0d0
 REAL  :: tCommS = 0d0
 REAL  :: t0,t1
END TYPE tStatistics
TYPE (tStatistics),save :: myStat

TYPE tProperties
 CHARACTER cName*7
 CHARACTER Material(2)*10
 REAL*8 :: Gravity(3)
 REAL*8, dimension(6) :: ForceScale = (/1d0, 1d0, 1d0, 1d0, 1d0, 1d0/)
 REAL*8 Density(2),Viscosity(2),DiffCoeff(2),Sigma,DiracEps,PowerLawExp
 REAL*8 :: ViscoLambda
 REAL*8 :: ViscoAlphaExp   =-0.1d0, ViscoAlphaImp   =+0.1d0
 REAL*8 :: NS_StabAlpha_Exp=-0.1d0, NS_StabAlpha_Imp=+0.1d0
 INTEGER :: ViscoModel,nTPSubSteps,nTPIterations,nTPFSubSteps
 INTEGER :: nInterface
END TYPE tProperties

TYPE tParamV
 REAL*8  defCrit,MinDef
 Real*8 :: Alpha
 INTEGER SolvType
 INTEGER iMass,NLmin,NLmax
 TYPE(tMGParamIn) :: MGprmIn
 TYPE(tMGParamOut):: MGprmOut(3)
END TYPE tParamV

TYPE tParamP
 TYPE(tMGParamIn) :: MGprmIn
 TYPE(tMGParamOut):: MGprmOut
END TYPE tParamP


TYPE mg_dVector
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: x
END TYPE mg_dVector

TYPE mg_kVector
 INTEGER  , DIMENSION(:) , ALLOCATABLE :: x
END TYPE mg_kVector

TYPE tGradient
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: x,y,z
END TYPE tGradient

TYPE TViscoScalar
 CHARACTER cName*7
 INTEGER :: ndof,na
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: val11,val22,val33,val12,val13,val23,rhs0,ValOld,diag
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: knpr
 TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: def,aux,rhs,sol
 TYPE(tGradient) :: grad11,grad22,grad33,grad12,grad13,grad23
 TYPE(tParamV) :: prm
 LOGICAL :: bProlRest=.FALSE.
END TYPE

TYPE TLinScalar
 CHARACTER cName*7
 INTEGER :: ndof,na
 TYPE(mg_kVector), DIMENSION(:),ALLOCATABLE :: knprP
! TYPE(mg_kVector), DIMENSION(:),ALLOCATABLE :: knpr
! INTEGER , DIMENSION(:)  , ALLOCATABLE :: knpr
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valP_old,P_old,P_new
 REAL*8  , DIMENSION(:)  , ALLOCATABLE ::  ST_P,Q2
 TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: valP,defP,auxP,rhsP,dvalP
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valP_GMV
 TYPE(tParamP) :: prm
 LOGICAL :: bProlRest=.FALSE.
END TYPE

TYPE TQuadScalar
 CHARACTER cName*7
 INTEGER :: ndof,na
 TYPE(mg_kVector), DIMENSION(:),ALLOCATABLE :: knprU
 TYPE(mg_kVector), DIMENSION(:),ALLOCATABLE :: knprV
 TYPE(mg_kVector), DIMENSION(:),ALLOCATABLE :: knprW
 TYPE(mg_kVector), DIMENSION(:),ALLOCATABLE :: knprT
! INTEGER , DIMENSION(:)  , ALLOCATABLE :: knpr
! In this part also valX_old1/2 and valX_help are added for the BDF time-stepping (locally allocated in q2p1_cc)
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: auxU,rhsU,defU,valU_old,valU,valU_old1,valU_help,valU_old2
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: auxV,rhsV,defV,valV_old,valV,valV_old1,valV_help,valV_old2
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: auxW,rhsW,defW,valW_old,valW,valW_old1,valW_help,valW_old2
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valUx,valUy,valUz
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valVx,valVy,valVz
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valWx,valWy,valWz
 TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: def,aux,rhs,sol,dsol
 TYPE(tParamV) :: prm
 LOGICAL :: bProlRest=.FALSE.
END TYPE

type fieldPtr
  real*8, dimension(:), pointer ::p
end type fieldPtr


TYPE TParLinScalar
 INTEGER :: ndof
 REAL*8  , DIMENSION(:)  , ALLOCATABLE::  Val
END TYPE

TYPE TMatrix
 INTEGER :: nu,na
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: ColA,LdA
END TYPE

TYPE TcrsStructure
! INTEGER, ALLOCATABLE :: A_Ld(:),A_Col(:)
 REAL*8, ALLOCATABLE  :: A_Mat(:),A_Rhs(:),A_Sol(:)
 TYPE(TMatrix) A
! INTEGER nu,na
END TYPE

TYPE lScalar3
 CHARACTER cName*7
 INTEGER :: ndof,na
 LOGICAL :: bProlRest=.FALSE. 

 INTEGER , DIMENSION(:)  , ALLOCATABLE :: knprX,knprY,knprZ
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valX_old,valY_old,valZ_old
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valX,valY,valZ
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: defX,defY,defZ
 TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: def,sol,aux,rhs
 TYPE(tParam) :: prm
 
END TYPE

TYPE lScalarField
 CHARACTER cName*7
 integer ID
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: knpr
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: val_old
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: val
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: rhs
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: def
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: aux
END TYPE

TYPE lScalarGen
 CHARACTER cName*7
 INTEGER :: ndof,na,nOfFields
 LOGICAL :: bProlRest=.FALSE. 
 TYPE(lScalarField), DIMENSION(:)  , ALLOCATABLE :: Fld
 TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: def,sol,aux,rhs
 INTEGER, DIMENSION(:),ALLOCATABLE :: knpr
 TYPE(tParam) :: prm
END TYPE

TYPE tVelo
 REAL*8, ALLOCATABLE :: x(:)
 REAL*8, ALLOCATABLE :: y(:)
 REAL*8, ALLOCATABLE :: z(:)
END TYPE tVelo

TYPE(tParticle), ALLOCATABLE :: myCompleteSet(:)
TYPE(tParticle), ALLOCATABLE :: myActiveSet(:)
TYPE(tParticle), ALLOCATABLE :: myExchangeSet(:)
TYPE(tParticle), ALLOCATABLE :: myLostSet(:)
TYPE(tVelo), ALLOCATABLE :: myVelo(:)

INTEGER nCompleteSet,nActiveSet,nExchangeSet,nStartActiveSet,nLostSet


TYPE tParticleInflow
 REAL*8 :: Center(3), Radius
END TYPE tParticleInflow

TYPE tPhysParticles
 REAL*8 :: rho_l, rho_p, d_p, mu_l
 REAL*8 :: gravity(3) 
END TYPE tPhysParticles

TYPE tParticleParam
 REAL*8 dEps1,dEps2, D_Out,D_in, f, Z_seed,Epsilon,hSize,d_CorrDist,OutflowZPos
 REAL*8 :: minFrac
 INTEGER :: nTimeLevels,nParticles,nRotation,Raster,dump_in_file,nPeriodicity
 ! What kind of starting procedure?
 INTEGER :: inittype
 ! If we start from a sourcefile this char stores the name of it
 CHARACTER(len=256) :: sourcefile = ''
 ! Scaling factors for input and output of units
 ! The dump-files (aka the flow-field and coordinates) the code reads
 ! in are in cm, so the particle-code internally will work in cm, too.
 ! Therefore, these factors are applied when reading in the parameters
 ! or when writing out the particle-files
 ! We also need to consider that the sourcefile we read in as particle-seed
 ! can now be in another unit as the cm of the particle-code!
 REAL*8 :: dFacUnitIn = 1.0d0
 REAL*8 :: dFacUnitOut = 1.0d0
 REAL*8 :: dFacUnitSourcefile = 1.0d0

 ! We can make z-cutplanes at as many positions as we want.
 ! However, we will limit it to 99
 integer :: nZposCutplanes
 
 integer :: DumpFormat=1
 real*8, allocatable :: cutplanePositions(:)

 ! Now many segments do we take for coloring?
 real*8 :: numberSegments
 
 LOGICAL :: bRotationalMovement=.true.

 LOGICAL :: bBacktrace=.false.,bPhysParticles=.false.
 
 INTEGER :: NumberOfInflowRegions=0
 TYPE(tParticleInflow), ALLOCATABLE :: InflowRegion(:)
 
 TYPE (tPhysParticles) :: PhysParticles
 
 !!!! Seeding 
 INTEGER  :: Plane = 0, PlaneParticles = 0, VolumeParticles = 0
 REAL*8   :: PlaneOffset=0d0
 
END TYPE tParticleParam

! Define parameters to find out where the particle-seed comes from
! From the parameterfile (aka RTD_param.dat) by construction
integer, parameter, public :: ParticleSeed_Parameterfile = 0
! Read in from a source.csv that contains 3 columns: X,Y,Z-Coordinates. First line is a header
integer, parameter, public :: ParticleSeed_CSVFILE = 1
! From old_output.csv. It contains 4 columns: X,Y,Z-Coordinates + the index. First line is a header
! This matches the description of the outputfiles this code produces - therefore this is good
! for benchmarks.
integer, parameter, public :: ParticleSeed_OUTPUTFILE = 2

! This setting makes possible to seed the particles in the element center of the hex mesh on the finest level
integer, parameter, public :: ParticleSeed_ELEMCENTER = 3

! This defines a planar seeding
integer, parameter, public :: ParticleSeed_PLANE = 8

! This defines a planar seeding
integer, parameter, public :: ParticleSeed_VOLUME = 9

TYPE tMeshInfoParticle
  real*8 xmin, xmax, ymin, ymax, zmin, zmax
END TYPE tMeshInfoParticle
type(tMeshInfoParticle) :: myMeshInfo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SIGMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
TYPE tPID
 REAL*8 T_set,sumI,e_old
 REAL*8 P,I,D,omega_P,omega_I,omega_D,PID
END TYPE tPID

TYPE tSensor
  REAL*8 :: Radius=0d0, Coor(3)=[0d0,0d0,0d0], Volume, CurrentTemperature
  LOGICAL :: HeatingStatus = .true.
  REAL*8 :: MinRegValue= 0d0, MaxRegValue= 0d0
END TYPE tSensor

TYPE tConvergenceDetector
  REAL*8  :: Condition=0.01d0
  INTEGER :: Counter  = 0
  INTEGER :: Limit  = 250
  LOGICAL :: Converged=.FALSE.
END TYPE tConvergenceDetector

TYPE tSegment
  INTEGER :: nOFFfilesR=0,nOFFfilesL=0,nOFFfiles=0
  CHARACTER*200, ALLOCATABLE :: OFFfilesR(:),OFFfilesL(:),OFFfiles(:)
  INTEGER, ALLOCATABLE :: idxCgal(:),idxCgalL(:),idxCgalR(:)
  REAL*8 offsetAngle
  CHARACTER*20 ObjectType,Unit
  CHARACTER*10 name
  CHARACTER*99 ::  cF
  CHARACTER*8 ::  ART
  INTEGER ::    KNETz,N
  REAL*8  :: Ds,s,delta,Dss,excentre,DiscFrac=0.05d0
  REAL*8, ALLOCATABLE :: Zknet(:)
  REAL*8 :: t,D,Alpha,StartAlpha ! t=Gangsteigung
  REAL*8 :: Min, Max,L
  REAL*8 :: ZME_DiscThick,ZME_gap_SG, ZME_gap_SS 
  REAL*8 :: SegRotFreq
  INTEGER :: ZME_N
  REAL*8  :: SecProf_W, SecProf_D,SecProf_L
  INTEGER :: SecProf_N, SecProf_I
  !!!!!!!!!!!!!!!!!!!!! EWIKON !!!!!!!!!!!!!!!!!!!!!
  INTEGER :: MatInd
  REAL*8 :: HeatSourceMax,HeatSourceMin,UseHeatSource
  character*64 :: regulation="SIMPLE"
  TYPE(tSensor) TemperatureSensor
  TYPE(tPID) PID_ctrl
  TYPE(tConvergenceDetector) ConvergenceDetector
  REAL*8 :: InitTemp,Volume
  CHARACTER*200 :: TemperatureBC
  !!!!!!!!!!!!!!!!!!!
  INTEGER GANGZAHL
   
END TYPE tSegment

TYPE tSigma
!   REAL*8 :: Dz_out,Dz_in, a, L, Ds, s, delta,SegmentLength, DZz,W
  CHARACTER cType*(50),cZwickel*(50),RotationAxis*(50)
  LOGICAL :: ScrewCylinderRendering=.TRUE.
  REAL*8 :: RotAxisCenter,RotAxisAngle
  REAL*8 :: Dz_out,Dz_in, a, L, L0, SegmentLength, DZz,W
  REAL*8 :: SecStr_W,SecStr_D
  INTEGER :: NumberOfMat,NumberOfSeg, GANGZAHL,STLSeg=0
  TYPE (tSegment), ALLOCATABLE :: mySegment(:)
  INTEGER :: InnerDiamNParam=0
  REAL*8,ALLOCATABLE ::  InnerDiamDParam(:),InnerDiamZParam(:)
  LOGICAL :: bOnlyBarrelAdaptation=.false., bAnalyticalShearRateRestriction=.false.
  
END TYPE tSigma

TYPE tRheology
   INTEGER :: Equation = 5 !-->> Hogen-Powerlaw
   INTEGER :: AtFunc = 1 !-->> isotherm
   REAL*8 :: A, B, C, D ! Carreau Parameter
   REAL*8 :: n, K ! Power Law
   REAL*8 :: eta_max, eta_min 
   REAL*8 :: Ts, Tb, C1, C2, E! WLF Parameter
   REAL*8 :: ViscoMin = 1e-4
   REAL*8 :: ViscoMax = 1e10
END TYPE tRheology

TYPE tSubInflow
 INTEGER iBCtype,Material
 REAL*8  massflowrate, outerradius,innerradius,temperature
 REAL*8  center(3),normal(3)
END TYPE tSubInflow

TYPE tInflow
 INTEGER :: nSubInflows=0
 TYPE (tSubInflow), dimension(:), allocatable :: mySubInflow
 
 INTEGER iBCtype,Material
 REAL*8  massflowrate, outerradius,innerradius,temperature
 REAL*8  center(3),normal(3)
 real*8, allocatable :: PressureEvolution(:)
END TYPE tInflow

TYPE tSegThermoPhysProp
 real*8 rho,cp,lambda,T_const
 character*256 :: cConstTemp
 logical :: bConstTemp=.false.
ENDTYPE tSegThermoPhysProp

TYPE tProcess
   REAL*8 :: Umdr, Ta, Ti, T0=0d0, T0_Slope=0d0, Massestrom, Dreh, Angle, dPress
   REAL*8 :: HeatFluxThroughBarrelWall_kWm2=0d0
   REAL*8 :: MinInflowDiameter,MaxInflowDiameter
   INTEGER :: Periodicity,Phase, nTimeLevels=36, nPeriodicity=1
   REAL*8 :: dAlpha
   REAL*8 :: ExtrusionGapSize,ExtrusionSpeed
   CHARACTER*6 :: Rotation !RECHT, LINKS
   CHARACTER*50 :: pTYPE !RECHT, LINKS
   INTEGER :: ind,iInd
   REAL*8 :: FBMVeloBC(3)=[0d0,0d0,0d0]
   integer   nOfInflows
   TYPE (tInflow), dimension(:), allocatable :: myInflow
   LOGICAL :: SegmentThermoPhysProps=.false.
   TYPE(tSegThermoPhysProp), allocatable :: SegThermoPhysProp(:)
  !!!!!!!!!!!!!!!!!!!!! EWIKON !!!!!!!!!!!!!!!!!!!!!
   REAL*8 :: AmbientTemperature=280d0,MeltInflowTemperature = 290d0
   REAL*8 :: WorkBenchThickness = 5d0, CoolingWaterTemperature = 55d0, ConductiveLambda = 21d0

!    REAL*8 :: TemperatureSensorRadius=0d0, TemperatureSensorCoor(3)=[0d0,0d0,0d0]
END TYPE tProcess

TYPE tThermodyn
   CHARACTER*60 :: DensityModel='NO'
   REAL*8 :: densityT0=0d0, density=0d0, densitySteig=0d0
   REAL*8 :: lambda=0d0, Cp=0d0, lambdaSteig=0d0,CpSteig=0d0
   REAL*8 :: Alpha=0d0, Beta=0d0, Gamma=0d0
END TYPE tThermodyn

TYPE tSingleMat
   CHARACTER*256 :: cMatNAme='UnknownMaterial'
   TYPE(tRheology)  :: Rheology
   TYPE(tThermodyn) :: Thermodyn
END TYPE tSingleMat

TYPE tMultiMat
   Integer :: nOfMaterials=1,initMaterial=1
   TYPE(tSingleMat) , Allocatable  :: Mat(:)
END TYPE tMultiMat

TYPE tTransientSolution
 INTEGER :: nTimeSubStep = 6, DumpFormat=2 ! LST
 TYPE(mg_dVector), ALLOCATABLE :: Velo(:,:)
 TYPE(mg_dVector), ALLOCATABLE :: Coor(:,:)
 TYPE(mg_dVector), ALLOCATABLE :: Dist(:)
 TYPE(mg_dVector), ALLOCATABLE :: Temp(:)
 TYPE(mg_dVector), ALLOCATABLE :: iSeg(:)
END TYPE tTransientSolution

TYPE tSetup
 LOGICAL :: bPressureFBM = .FALSE.
 REAL*8  :: PressureConvergenceTolerance 
 Logical :: bPressureConvergence= .FALSE.
 LOGICAL :: bAutomaticTimeStepControl = .TRUE.,bRotationalFramOfReference=.FALSE.
 LOGICAL :: bConvergenceEstimator=.FALSE.
 REAL*8 :: CharacteristicShearRate=1d1
 CHARACTER*200 cMeshPath
 CHARACTER*20 cMesher
 INTEGER MeshResolution
 INTEGER m_nT,m_nT1,m_nT2,m_nR,m_nZ,m_nP,m_nX,m_nY,nBoxElem
 REAL*8 m_box(3,2)
 LOGICAL :: bGeoTest=.FALSE.,bSendEmail=.TRUE.
END TYPE tSetup

TYPE tOutput
 INTEGER :: nOf1DLayers=16
 INTEGER :: nOfHistogramBins=16
 REAL*8 ::  HistogramShearMax=1e6,HistogramShearMin=1e-2,HistogramViscoMax=1e6,HistogramViscoMin=1e0
 REAL*8  :: CutDtata_1D=0.04d0
END TYPE tOutput

TYPE tErrorCodes
 INTEGER :: WRONG_SCREW=88
 INTEGER :: DIVERGENCE_U=55
 INTEGER :: DIVERGENCE_T=56
 INTEGER :: SIGMA_READER=44
END TYPE

TYPE tMGSteps
 INTEGER :: m=2
 INTEGER :: i,j,iaux
 real*8  :: daux
 INTEGER, allocatable :: n(:)
 REAL*8, allocatable :: r(:)
END TYPE tMGSteps

end module types
