MODULE var_QuadScalar

use def_FEAT
use types

IMPLICIT NONE

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------

CHARACTER*200 :: ApplicationString=&
"  |                                                                                                  |"

CHARACTER*200 :: VersionString=&
"  |                                                          Version:22.01  Date:2022.02.02          |"

CHARACTER*200 :: myDataFile="_data/q2p1_param.dat"

! There seems to be a problem with this setting 
! when the intel compiler is used, mesh deformation is
! used in every time step (i.e. setting M2D2C2) AND
! the game-of-thrones servers are used. It seems to 
! run fine on BTTF and the LIDO-servers
INTEGER :: iCommSwitch=3
LOGICAL :: BaSynch=.false.
LOGICAL :: bParallel=.true.
LOGICAL :: bMasterTurnedOn=.true. ! because of hypre
LOGICAL :: bMultiMat=.false.
LOGICAL :: DivergedSolution=.false., ConvergedSolution = .false., bAlphaConverged=.false.
REAL*8  :: AlphaControl=0d0

TYPE tSSE_covergence
 REAL*8, allocatable :: Monitor(:)
 REAL*8              :: average,std_dev,dCharVisco
 integer             :: iC,nC,start
END TYPE
TYPE (tSSE_covergence) :: mySSE_covergence

LOGICAL :: SSE_HAS_ANGLE=.false.
real*8  :: extruder_angle = 0.0

integer :: MaxLevelKnownToMaster

TYPE tMatrixRenewal
INTEGER K,D,M,S,C
END TYPE tMatrixRenewal
TYPE(tMatrixRenewal) myMatrixRenewal
LOGICAL :: bNonNewtonian=.TRUE.
LOGICAL :: bNoOutflow,bTracer,bViscoElastic,bRefFrame
LOGICAL :: b2DViscoBench=.FALSE.,b3DViscoBench=.FALSE.

LOGICAL :: bSteadyState =.FALSE.
LOGICAL :: bBoundaryCheck=.FALSE.
LOGICAL :: bNS_Stabilization=.FALSE.

! Integer parameter for terminal output
integer, parameter :: uterm = 6

REAL*8  :: dCGALtoRealFactor = 1d0
INTEGER, PARAMETER :: Giesekus = 0
INTEGER, PARAMETER :: OldroydB = 1
INTEGER :: ProlongationDirection = 0
REAL*8 :: activeFBM_Z_Position=-1d9
REAL*8 :: dTimeStepEnlargmentFactor=1d0
INTEGER :: iTimeStepEnlargmentFactor=1

TYPE tTransform
 INTEGER :: ILINT=2
END TYPE
TYPE (tTransform) Transform


TYPE(TQuadScalar), target :: QuadSc
! TODO: Move these to var
TYPE(TLinScalar)    LinSc
TYPE(TViscoScalar)  ViscoSc
TYPE(TParLinScalar) PLinSc

TYPE(lScalar) Tracer
TYPE(lScalar3) Tracer3
TYPE(lScalarGen) GenLinScalar

TYPE(TcrsStructure) crsSTR

INTEGER , DIMENSION(:)  , ALLOCATABLE :: lqK

TYPE(TMatrix), POINTER :: qlMat,qlPMat,lqMat,lMat,lPMat,qMat

TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_qlMat,mg_qlPMat,mg_lqMat
TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_lMat,mg_lPMat,mg_qMat

REAL*8, DIMENSION(:), POINTER :: BXMat,BYMat,BZMat
REAL*8, DIMENSION(:), POINTER :: BTXMat,BTYMat,BTZMat
REAL*8, DIMENSION(:), POINTER :: BXPMat,BYPMat,BZPMat
REAL*8, DIMENSION(:), POINTER :: Mmat,MlMat,MlPmat,MlRhomat,MlRhoPmat
REAL*8, DIMENSION(:), POINTER :: DMat,Kmat,A11mat,A22mat,A33mat,ConstDMat,hDMat
REAL*8, DIMENSION(:), POINTER :: A12mat,A13mat,A23mat,A21mat,A31mat,A32mat
REAL*8, DIMENSION(:), POINTER :: S11mat,S22mat,S33mat
REAL*8, DIMENSION(:), POINTER :: S12mat,S13mat,S23mat,S21mat,S31mat,S32mat
REAL*8, DIMENSION(:), POINTER :: Cmat,CPMat
REAL*8, DIMENSION(:), POINTER :: VisMat_11,VisMat_22,VisMat_33
REAL*8, DIMENSION(:), POINTER :: VisMat_12,VisMat_13,VisMat_23

TYPE(TMatrix)          :: UMF_lMat
REAL*8 , ALLOCATABLE   :: UMF_CMat(:)

TYPE tGlobalNumberingMap
 INTEGER  :: ndof,ndof_Q2,ndof_P1
 INTEGER , allocatable :: ind(:)
 INTEGER , allocatable :: indQ2(:),indP1(:)
 REAL*8, allocatable   :: dBufferQ2(:),dBufferP1(:)
END TYPE tGlobalNumberingMap
TYPE(tGlobalNumberingMap), ALLOCATABLE :: myGlobalNumberingMap(:)
INTEGER, ALLOCATABLE :: GlobalNumberingQ2(:),GlobalNumberingP1(:)
INTEGER myGlobal_ndof

TYPE(mg_kVector),ALLOCATABLE :: GlobalParallelList1(:),GlobalParallelList2(:)
REAL*8, ALLOCATABLE          :: GlobalParallelBufferOut(:),GlobalParallelBufferIn(:)
INTEGER, ALLOCATABLE         :: GlobalNList(:),GlobalNBuffer(:)

TYPE(mg_kVector),ALLOCATABLE :: HlobalParallelList1(:),HlobalParallelList2(:),HlobalParallelList3(:)
REAL*8, ALLOCATABLE          :: HlobalParallelBufferOut(:),HlobalParallelBufferIn(:)
INTEGER, ALLOCATABLE         :: HlobalNList(:),HlobalNBuffer(:)

TYPE(mg_kVector),ALLOCATABLE :: CommOrder(:)

TYPE mg_Matrix
 REAL*8  , DIMENSION(:)  , ALLOCATABLE  :: a
END TYPE mg_Matrix

TYPE tMGFldMatrix
 TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE  :: fld
END TYPE

REAL*8 :: NewtonForBurgers=0d0
REAL*8  , DIMENSION(:)  , POINTER :: barM11mat,barM22mat,barM33mat,barM12mat,barM13mat,barM23mat,barM21mat,barM31mat,barM32mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM11mat,mg_barM22mat,mg_barM33mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM12mat,mg_barM13mat,mg_barM23mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM21mat,mg_barM31mat,mg_barM32mat

TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_BXMat,mg_BYMat,mg_BZMat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_BTXMat,mg_BTYMat,mg_BTZMat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_BXPMat,mg_BYPMat,mg_BZPMat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_DMat,mg_KMat,mg_ConstDMat,mg_hDMat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_A11mat,mg_A22mat,mg_A33mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_A12mat,mg_A13mat,mg_A23mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_A21mat,mg_A31mat,mg_A32mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_S11mat,mg_S22mat,mg_S33mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_S12mat,mg_S13mat,mg_S23mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_S21mat,mg_S31mat,mg_S32mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_MMat,mg_MlMat,mg_MlPMat,mg_MlRhomat,mg_MlRhoPmat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_CMat,mg_CPMat,mg_P1MMat,mg_P1iMMat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_VisMat_11,mg_VisMat_22,mg_VisMat_33
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_VisMat_12,mg_VisMat_13,mg_VisMat_23

TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_E012Prol,mg_E013Prol,mg_E013Rest
TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_E013ProlM,mg_E013RestM

TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_E011Prol,mg_E011Rest
TYPE(TMatrix)   , DIMENSION(:)  , ALLOCATABLE,  TARGET :: mg_E011ProlM,mg_E011RestM

TYPE tMultiGrid
 CHARACTER*20 cVariable
 CHARACTER*1 CycleType

 CHARACTER*10 MGProlongation

 LOGICAL, POINTER :: bProlRest
 INTEGER, DIMENSION(:), POINTER::  KNPRU,KNPRV,KNPRW
 INTEGER, DIMENSION(:), POINTER::  KNPR
 TYPE(mg_kVector), DIMENSION(:), POINTER::  KNPRP

 TYPE(mg_dVector), DIMENSION(:), POINTER::  X_u,dX_u,D_u,A_u,B_u
 TYPE(mg_dVector), DIMENSION(:), POINTER::  X_p,dX_p,D_p,A_p,B_p
 TYPE (mg_Matrix), DIMENSION(:), POINTER :: BX,BY,BZ,BTX,BTY,BTZ
 TYPE(TMatrix), DIMENSION(:),  POINTER   :: Lq,Lql,Llq


 TYPE(mg_dVector), DIMENSION(:), POINTER::  X,D,AUX,B
 TYPE (mg_Matrix), DIMENSION(:), POINTER :: A,AP
 TYPE (mg_Matrix), DIMENSION(:), POINTER :: A11,A22,A33,A12,A13,A23,A21,A31,A32
 TYPE (tMGFldMatrix), DIMENSION(:), POINTER :: AXX
 TYPE(TMatrix), DIMENSION(:),  POINTER :: L,LP
 REAL*8  , DIMENSION(:)  , POINTER :: XP
 INTEGER MinLev,MaxLev,MedLev,MaxDifLev,MinIterCycle,MaxIterCycle,nIterCoarse,nSmootherSteps,CrsSolverType,SmootherType
 INTEGER :: nOfFields,nOfSubsystemEqs
 integer :: vanka
 REAL*8  DefImprCoarse
 REAL*8  Criterion1,Criterion2,RLX,CrsRelaxPrm,CrsRelaxParPrm
! MGOutputs
 REAL*8  RhoMG1,RhoMG2,DefInitial,DefFinal
 INTEGER UsedIterCycle

END TYPE tMultiGrid

TYPE (tMultiGrid) :: myMG

TYPE tMultiGrid_cc
 CHARACTER*10 cVariable
 CHARACTER*1  CycleType
 CHARACTER*10 MGProlongation
 LOGICAL, POINTER :: bProlRest
 INTEGER, DIMENSION(:), POINTER::  KNPRU,KNPRV,KNPRW

 TYPE(mg_dVector), DIMENSION(:), POINTER::  X_u,dX_u,D_u,A_u,B_u
 TYPE(mg_dVector), DIMENSION(:), POINTER::  X_p,dX_p,D_p,A_p,B_p
 TYPE (mg_Matrix), DIMENSION(:), POINTER :: BX,BY,BZ,BTX,BTY,BTZ
 TYPE(TMatrix), DIMENSION(:),  POINTER   :: Lq,Lql,Llq

 TYPE(mg_dVector), DIMENSION(:), POINTER::  X,D,AUX,B
 TYPE (mg_Matrix), DIMENSION(:), POINTER :: A,AP
 TYPE (mg_Matrix), DIMENSION(:), POINTER :: A11,A22,A33,A12,A13,A23,A21,A31,A32
 TYPE(TMatrix), DIMENSION(:),  POINTER :: L,LP
 REAL*8  , DIMENSION(:)  , POINTER :: XP
 INTEGER MinLev,MaxLev,MedLev,MinIterCycle,MaxIterCycle,nIterCoarse,nSmootherSteps,VANKA
 REAL*8  DefImprCoarse
 REAL*8  Criterion1,Criterion2,RLX
! MGOutputs
 REAL*8  RhoMG1,RhoMG2,DefInitial,DefFinal
 INTEGER UsedIterCycle

END TYPE tMultiGrid_cc
TYPE (tMultiGrid_cc) :: myMG_cc

TYPE(mg_kVector), ALLOCATABLE :: MaterialDistribution(:)

TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: mgDensity(:),mgDiffCoeff(:)
TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: mgNormShearStress(:)

type(tMultiMesh),save :: mg_mesh

INTEGER, ALLOCATABLE :: ParKNPR(:)
INTEGER, ALLOCATABLE :: FictKNPR(:),MixerKnpr(:)
REAL*8, ALLOCATABLE :: Distance(:),Distamce(:),Screw(:),Shell(:),ScrewDist(:,:)
REAL*8, ALLOCATABLE :: Viscosity(:), Shearrate(:),Temperature(:),ElemSizeDist(:),MaxShearRate(:)
REAL*8, ALLOCATABLE :: Temperature_AVG(:)
INTEGER :: iTemperature_AVG = 0

TYPE tCGALObjects
 REAL*8, ALLOCATABLE :: Block(:)
 REAL*8, ALLOCATABLE :: Wire(:)
 REAL*8, ALLOCATABLE :: Channel(:)
 INTEGER, ALLOCATABLE :: Sensor(:)
 INTEGER, ALLOCATABLE :: Segment(:)
END TYPE tCGALObjects
TYPE (tCGALObjects) myHeatObjects

TYPE tParticleFBM
 CHARACTER cTYPE*10
 REAL*8 sizes(20),density
 REAL*8 ResistanceForce(3),TorqueForce(3)
 REAL*8 Position(3),Velocity(3),Angle(3),AngularVelocity(3)
 REAL*8 Acceleration(3),FrameVelocity(3)
END TYPE tParticleFBM

TYPE tFBM
 INTEGER nParticles
 REAL*8,ALLOCATABLE :: Force(:)
 TYPE (tParticleFBM), ALLOCATABLE :: ParticleOld(:),ParticleNew(:)
 integer, allocatable, dimension(:) :: iel_ug   
END TYPE tFBM

CHARACTER cFBM_File*30
TYPE (tFBM) myFBM

LOGICAL, ALLOCATABLE :: BndrForce(:)
REAL*8, ALLOCATABLE :: myQ2Coor(:,:), BoundaryNormal(:,:)

TYPE tBoundary
 LOGICAL, ALLOCATABLE :: nWall(:)   ,nInflown(:),  nOutflow(:)   ,nSymmetry(:)
 LOGICAL, ALLOCATABLE :: bWall(:),bOutflow(:),bSymmetry(:,:),bSlip(:)
 INTEGER, ALLOCATABLE :: iInflow(:), iPhase(:),iTemperature(:)
 INTEGER, ALLOCATABLE :: LS_zero(:)
 LOGICAL, ALLOCATABLE :: bDisp_DBC(:) 
END TYPE tBoundary
TYPE (tBoundary) myBoundary

TYPE tDump
 INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Elements,Vertices
END TYPE tDump
TYPE(tDump) :: myDump

REAL*8, ALLOCATABLE :: dPeriodicVector(:)
REAL*8, allocatable :: mySegmentIndicator(:,:)

TYPE tExport
 INTEGER :: Level,LevelMax
 CHARACTER*(3) :: Format
 CHARACTER*(40),ALLOCATABLE, DIMENSION(:) :: Fields
 CHARACTER*(4) :: cFileName(2)=['res ','main']
END TYPE tExport
TYPE(tExport) :: myExport

INTEGER :: iOutput=0
CHARACTER cProjectGridFile*220,cGridFileName*200,cProjectFile*200,cProjectFolder*200,cProjectNumber*4
INTEGER nSubCoarseMesh

CHARACTER*13 :: outfile="OutFile  .txt"

!                            Changes for CC
!---------------------------------------------------------------------------
TYPE tElementMatrix
 REAL*8 A(85,85)
 INTEGER H(2)
END TYPE

TYPE tMGElementMatrix
 TYPE(tElementMatrix), ALLOCATABLE :: E(:)
END TYPE

TYPE(tmgElementMatrix), ALLOCATABLE :: CC_EMat(:)

TYPE(TMatrix)          :: CC_crs_lMat
REAL*8 , ALLOCATABLE   :: CC_crs_AMat(:)
INTEGER CC_H(2)

TYPE crs_e013_map
 INTEGER  :: ndof,cc_ndof
 INTEGER , allocatable :: ind(:)
 INTEGER , allocatable :: indE(:)
 REAL*8, allocatable   :: dBuffer(:)
END TYPE crs_e013_map

TYPE(crs_e013_map), ALLOCATABLE :: my_crs_e013_map(:)

TYPE CC_Elem
 INTEGER, ALLOCATABLE :: pairE(:,:),pairV(:,:)
 INTEGER, ALLOCATABLE :: E_qq(:),E_lq(:),E_ql(:)
 REAL*8 , ALLOCATABLE :: a(:)
 INTEGEr sym,num
END TYPE

TYPE mg_CCPiece
 INTEGER, ALLOCATABLE :: a(:,:)
 INTEGER, ALLOCATABLE :: LdA_qq(:),ColA_qq(:)
 INTEGER, ALLOCATABLE :: LdA_lq(:),ColA_lq(:)
 INTEGER, ALLOCATABLE :: LdA_ql(:),ColA_ql(:)
 TYPE(TMatrix):: MPatch,MPatchCopy
 TYPE(CC_Elem), ALLOCATABLE :: E(:)
 INTEGER :: nu_qq,na_qq,nu_ql,na_ql,nu_lq,na_lq
END TYPE mg_CCPiece
TYPE(mg_CCPiece),ALLOCATABLE :: my_mg_CCPiece(:)
INTEGER, ALLOCATABLE :: GlobalNumbering(:)

TYPE tCoarseMat
 INTEGER na,nu
 INTEGER, dimension(:), allocatable   :: Row
 INTEGER, dimension(:), allocatable   :: Col

 REAL*8 , ALLOCATABLE   :: A(:)
 REAL*8 , ALLOCATABLE   :: D(:)
END TYPE tCoarseMat

TYPE(tCoarseMat), target :: myCrsMat
!---------------------------------------------------------------------------

integer :: istep_ns = 1

TYPE tTriangle
 REAL*8 :: C(3,9)
END TYPE tTriangle

TYPE tTriangulation
 INTEGER nT
 TYPE(tTriangle), ALLOCATABLE :: T(:)
END TYPE tTriangulation

TYPE(tTriangulation), ALLOCATABLE :: myTSurf(:)

LOGICAL :: bMeshAdaptation = .FALSE.
CHARACTER*100 cAdaptedMeshFile
INTEGER nUmbrellaSteps,nInitUmbrellaSteps


integer :: nUmbrellaStepsLvl(9) = (/0, 0, 0, 0, 0, 0, 0, 0, 0/)
integer :: nMainUmbrellaSteps = 0
REAL*8 :: dIntegralHeat = 0d0, dNozzlePosition

TYPE tALE
 REAL*8, ALLOCATABLE :: Monitor(:)
 REAL*8, ALLOCATABLE :: OldCoor(:,:), NewCoor(:,:),MeshVelo(:,:), OrigCoor(:,:)
 REAL*8, ALLOCATABLE :: Q2coor_old(:,:)
 LOGICAL :: bUseFrameVelocity = .false.
 REAL*8 dFrameVelocity(3),dFrameVelocityChange(3)
END TYPE tALE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! interface structure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE tTetInterface
 INTEGER nG,nT
 REAL*8  , ALLOCATABLE :: X(:,:),Y(:,:)
 INTEGER , ALLOCATABLE :: L(:)
END TYPE tTetInterface

TYPE(tTetInterface) lInterface,gInterface
REAL*8, allocatable :: myFracField(:),myDistance(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! interface structure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TYPE(tALE),save :: myALE

TYPE(tBoundingBox), dimension(:), allocatable :: mgBoundingBox

type(tPostprocessingParams) :: postParams

TYPE t1DOutput
 REAL*8, ALLOCATABLE :: dMean(:),dMin(:),dMax(:),dLoc(:)
 CHARACTER cName*20
END TYPE t1DOutput
TYPE(t1DOutput), TARGET :: my1DOut(11)
REAL*8, ALLOCATABLE :: my1DIntervals(:,:),my1DWeight(:),my1DTorque(:,:)
INTEGER my1DOut_nol

TYPE tViscFunc
 REAL*8 :: shear_rate(66),visc(66)
END TYPE tViscFunc

TYPE(tViscFunc) :: myViscFunc

TYPE(tProperties),save :: Properties

TYPE (tParticleParam) :: myParticleParam

TYPE (tErrorCodes) ::  myErrorCode

TYPE (tMGSteps) MGSteps

contains 
integer function KNEL(ilevel)
  implicit none
  integer, intent(in) :: ilevel

   KNEL = mg_mesh%level(ilevel)%nel
  return 
end function KNEL

integer function KNET(ilevel)
  implicit none
  integer, intent(in) :: ilevel

   KNET = mg_mesh%level(ilevel)%net
  return 
end function KNET

integer function KNAT(ilevel)
  implicit none
  integer, intent(in) :: ilevel

   KNAT = mg_mesh%level(ilevel)%nat
  return 
end function KNAT

integer function KNVT(ilevel)
  implicit none
  integer, intent(in) :: ilevel

   KNVT = mg_mesh%level(ilevel)%nvt
  return 
end function KNVT

END MODULE var_QuadScalar

