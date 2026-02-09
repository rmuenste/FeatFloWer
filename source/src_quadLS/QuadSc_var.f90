MODULE var_QuadScalar

  use def_FEAT
  use types

  IMPLICIT NONE

  ! Workspace & metadata
  ! Add new workspace buffers and identifying strings within this block.
  INTEGER  NNWORK
  PARAMETER (NNWORK=1)
  INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)
  INTEGER            :: KWORK(1)
  REAL               :: VWORK(1)
  DOUBLE PRECISION   :: DWORK(NNWORK)
  COMMON       NWORK,IWORK,IWMAX,L,DWORK
  EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
  CHARACTER*200 :: ApplicationString=&
  "  |                                                                                                  |"
  CHARACTER*200 :: VersionString=&
  "  |                                                          Version:22.01  Date:2022.02.02          |"
  CHARACTER*200 :: myDataFile="_data/q2p1_param.dat"

  ! Runtime control
  ! Extend runtime flags, timers, and control scalars here.
  ! There seems to be a problem with this setting
  ! when the intel compiler is used, mesh deformation is
  ! used in every time step (i.e. setting M2D2C2) AND
  ! the game-of-thrones servers are used. It seems to
  ! run fine on BTTF and the LIDO-servers
  INTEGER :: iCommSwitch=3
  LOGICAL :: BaSynch=.false.
  LOGICAL :: bParallel=.true.
  LOGICAL :: bMemoryPrint=.true.
  LOGICAL :: bMasterTurnedOn=.true. ! because of hypre
  LOGICAL :: bMultiMat=.false.
  LOGICAL :: DivergedSolution=.false., ConvergedSolution = .false., bAlphaConverged=.false.
  REAL*8  :: AlphaControl=0d0
  REAL*8  :: total_lubrication = 0.0
  TYPE tTimer
    integer :: n = 0
    real*8  :: t = 0d0
  END TYPE tTimer
  TYPE (tTimer) :: myTimer(6),MYCOMMTIMER(8)
  TYPE tSSE_covergence
    REAL*8, allocatable :: Monitor(:)
    REAL*8              :: average,std_dev,dCharVisco
    integer             :: iC,nC,start
  END TYPE
  TYPE (tSSE_covergence) :: mySSE_covergence
  LOGICAL :: SSE_HAS_ANGLE=.false.
  real*8  :: extruder_angle = 0.0
  integer :: MaxLevelKnownToMaster
  real*8 :: referenceVelocity = 0.0
  TYPE tMatrixRenewal
    INTEGER K,D,M,S,C
  END TYPE tMatrixRenewal
  TYPE(tMatrixRenewal) :: myMatrixRenewal
  LOGICAL :: bNonNewtonian=.TRUE.
  LOGICAL :: bNoOutflow,bTracer,bViscoElastic,bRefFrame
  LOGICAL :: b2DViscoBench=.FALSE.,b3DViscoBench=.FALSE.
  LOGICAL :: bSteadyState =.FALSE.
  LOGICAL :: bBoundaryCheck=.FALSE.
  LOGICAL :: bNS_Stabilization=.FALSE.
  LOGICAL :: skipFBMForce=.FALSE.,skipFBMDynamics=.FALSE.
  LOGICAL :: bBinaryVtkOutput=.FALSE.
  REAL*8 :: Gamma = 0d0
  integer, parameter :: uterm = 6
  REAL*8  :: dCGALtoRealFactor = 1d0
  INTEGER, PARAMETER :: Giesekus = 0
  INTEGER, PARAMETER :: OldroydB = 1
  INTEGER :: ProlongationDirection = 0
  REAL*8 :: activeFBM_Z_Position=-1d9
  REAL*8 :: dTimeStepEnlargmentFactor=1d0
  INTEGER :: iTimeStepEnlargmentFactor=1
  real*8 :: GammaDot = 0.0d0
  real*8 :: AlphaRelax = 0.0d0
  real*8 :: RadParticle = 1.0d0
  real*8 :: RPM = 12.0d0
  TYPE tTransform
    INTEGER :: ILINT=2
  END TYPE
  TYPE (tTransform) :: Transform
  integer :: istep_ns = 1

  ! Solver matrices
  ! Place solver handles, matrix pointers, and assembled blocks below.
  ! Coupled Q2 system descriptor hosting assembled velocity matrices
  TYPE(TQuadScalar), target :: QuadSc
  ! P1 scalar container storing pressure-related matrices
  TYPE(TLinScalar)    :: LinSc
  ! Viscoelastic scalar block collecting stress/extra field matrices
  TYPE(TViscoScalar)  :: ViscoSc
  ! Distributed-memory linear scalar helper for parallel solves
  TYPE(TParLinScalar) :: PLinSc
  ! Passive scalar matrix container for tracer transport
  TYPE(lScalar) :: Tracer
  ! Vector-valued tracer matrices storing three coupled scalar fields
  TYPE(lScalar3) :: Tracer3
  ! Generic scalar descriptor used by auxiliary linear subsystems
  TYPE(lScalarGen) :: GenLinScalar
  ! Sparse CRS structure reused when building linearised operators
  TYPE(TcrsStructure) :: crsSTR
  ! Level-to-level coupling pattern between linear and quadratic grids
  INTEGER , DIMENSION(:)  , ALLOCATABLE :: lqK
  ! Pointers to active system matrices for l<->q couplings and velocity blocks
  TYPE(TMatrix), POINTER :: qlMat,qlPMat,lqMat,lMat,lPMat,qMat
  ! Multigrid hierarchy copies of ql / qlP / lq transfer matrices
  TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_qlMat,mg_qlPMat,mg_lqMat
  ! Multigrid hierarchy copies of linear and quadratic system matrices
  TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_lMat,mg_lPMat,mg_qMat
  ! Coupling operators mapping velocity divergence into pressure space (B)
  REAL*8, DIMENSION(:), POINTER :: BXMat,BYMat,BZMat
  ! Transposed coupling operators providing B^T contributions
  REAL*8, DIMENSION(:), POINTER :: BTXMat,BTYMat,BTZMat
  ! Coupling operators for pressure mass matrices when prolongating
  REAL*8, DIMENSION(:), POINTER :: BXPMat,BYPMat,BZPMat
  ! Mass matrices for velocity, lumped mass, and density-weighted variants
  REAL*8, DIMENSION(:), POINTER :: Mmat,MlMat,MlPmat,MlRhomat,MlRhoPmat
  ! Diffusion, stiffness, and stabilisation matrices (diagonal tensor blocks)
  REAL*8, DIMENSION(:), POINTER :: DMat,Kmat,A11mat,A22mat,A33mat,ConstDMat,hDMat
  ! Off-diagonal tensor blocks capturing cross-derivative couplings
  REAL*8, DIMENSION(:), POINTER :: A12mat,A13mat,A23mat,A21mat,A31mat,A32mat
  ! Symmetric stress matrix diagonal entries for viscoelastic models
  REAL*8, DIMENSION(:), POINTER :: S11mat,S22mat,S33mat
  ! Off-diagonal stress matrix entries for viscoelastic models
  REAL*8, DIMENSION(:), POINTER :: S12mat,S13mat,S23mat,S21mat,S31mat,S32mat
  ! Antisymmetric vorticity matrix diagonal components
  REAL*8, DIMENSION(:), POINTER :: W11mat,W22mat,W33mat
  ! Antisymmetric vorticity matrix off-diagonal components
  REAL*8, DIMENSION(:), POINTER :: W12mat,W13mat,W23mat,W21mat,W31mat,W32mat
  ! Constraint, pressure-mass, and P1 interpolation matrices
  REAL*8, DIMENSION(:), POINTER :: Cmat,CPMat,P1MMat,P1iMMat
  ! Viscosity matrix diagonal components (normal stresses)
  REAL*8, DIMENSION(:), POINTER :: VisMat_11,VisMat_22,VisMat_33
  ! Viscosity matrix off-diagonal components (shear stresses)
  REAL*8, DIMENSION(:), POINTER :: VisMat_12,VisMat_13,VisMat_23
  ! UMFPACK factorised linear matrix used for direct solves
  TYPE(TMatrix)          :: UMF_lMat
  ! Dense storage of UMFPACK factor values corresponding to UMF_lMat
  REAL*8 , ALLOCATABLE   :: UMF_CMat(:)
  ! UMFPACK CRS grid permutation describing pivoting pattern
  INTEGER, allocatable   :: UNF_P_CrsGrid(:)
  TYPE tGlobalNumberingMap
    ! Total number of dofs for mixed Q2/P1 system and its sub-blocks
    INTEGER  :: ndof,ndof_Q2,ndof_P1
    ! Global index mappings for the combined and split finite element spaces
    INTEGER , allocatable :: ind(:)
    INTEGER , allocatable :: indQ2(:),indP1(:)
    ! Workspace buffers used when scattering between numbering schemes
    REAL*8, allocatable   :: dBufferQ2(:),dBufferP1(:)
  END TYPE tGlobalNumberingMap
  ! Global numbering per subdomain required when assembling level operators
  TYPE(tGlobalNumberingMap), ALLOCATABLE :: myGlobalNumberingMap(:)
  ! Flat global numbering for velocity (Q2) and pressure (P1) unknowns
  INTEGER, ALLOCATABLE :: GlobalNumberingQ2(:),GlobalNumberingP1(:)
  ! Combined total number of degrees of freedom in the active system
  INTEGER :: myGlobal_ndof
  TYPE mg_Matrix
    ! Values of a matrix stored for one multigrid level
    REAL*8  , DIMENSION(:)  , ALLOCATABLE  :: a
  END TYPE mg_Matrix
  TYPE tMGFldMatrix
    ! Field-wise matrix collections (one mg_Matrix per extra field)
    TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE  :: fld
  END TYPE
  ! Newton update scaling factor used in Burgers-type linearisations
  REAL*8 :: NewtonForBurgers=0d0
  ! Barred mass matrices storing time-averaged values for each component
  REAL*8  , DIMENSION(:)  , POINTER :: barM11mat,barM22mat,barM33mat,barM12mat,barM13mat,barM23mat,barM21mat,barM31mat,barM32mat
  ! Multigrid level storage for barred mass matrix diagonal components
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM11mat,mg_barM22mat,mg_barM33mat
  ! Multigrid level storage for barred mass matrix off-diagonal components
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM12mat,mg_barM13mat,mg_barM23mat
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM21mat,mg_barM31mat,mg_barM32mat
  ! Multigrid copies of B operators used in coarse transfer steps
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_BXMat,mg_BYMat,mg_BZMat
  ! Multigrid copies of transposed B operators
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_BTXMat,mg_BTYMat,mg_BTZMat
  ! Multigrid copies of pressure mass matrices for prolongation phases
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_BXPMat,mg_BYPMat,mg_BZPMat
  ! Multigrid diffusion and stiffness matrices (diagonal components)
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_DMat,mg_KMat,mg_ConstDMat,mg_hDMat
  ! Multigrid stiffness tensor diagonal storage (A11/A22/A33)
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_A11mat,mg_A22mat,mg_A33mat
  ! Multigrid stiffness tensor off-diagonal storage (A12...A32)
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_A12mat,mg_A13mat,mg_A23mat
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_A21mat,mg_A31mat,mg_A32mat
  ! Multigrid viscoelastic stress tensor diagonal blocks
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_S11mat,mg_S22mat,mg_S33mat
  ! Multigrid viscoelastic stress tensor off-diagonal blocks
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_S12mat,mg_S13mat,mg_S23mat
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_S21mat,mg_S31mat,mg_S32mat
  ! Multigrid vorticity matrix diagonal entries for rotational smoothing
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_W11mat,mg_W22mat,mg_W33mat
  ! Multigrid vorticity matrix off-diagonal entries
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_W12mat,mg_W13mat,mg_W23mat
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_W21mat,mg_W31mat,mg_W32mat
  ! Multigrid velocity mass and density-weighted mass matrices
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_MMat,mg_MlMat,mg_MlPMat,mg_MlRhomat,mg_MlRhoPmat
  ! Multigrid constraint, pressure mass, and interpolation matrices
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_CMat,mg_CPMat,mg_P1MMat,mg_P1iMMat
  ! Multigrid viscosity tensor diagonal entries for viscoelastic solves
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_VisMat_11,mg_VisMat_22,mg_VisMat_33
  ! Multigrid viscosity tensor off-diagonal entries for viscoelastic solves
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_VisMat_12,mg_VisMat_13,mg_VisMat_23
  ! Prolongation/restriction matrices for E012/E013 element couplings
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_E012Prol,mg_E013Prol,mg_E013Rest
  ! Explicit matrix objects for E013 prolongation/restriction operators
  TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_E013ProlM,mg_E013RestM
  ! Prolongation/restriction matrices for E011 element couplings
  TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_E011Prol,mg_E011Rest
  ! Explicit matrix objects backing the E011 prolongation/restriction operators
  TYPE(TMatrix)   , DIMENSION(:)  , ALLOCATABLE,  TARGET :: mg_E011ProlM,mg_E011RestM
  ! Multigrid smoothing step configuration and statistics
  TYPE (tMGSteps) :: MGSteps

  ! Multigrid hierarchy
  ! Multigrid structures, coarse matrices, and related helpers go here.
  TYPE tMultiGrid
    ! Human-readable identifier for the multigrid solve (e.g. "Velocity")
    CHARACTER*20 :: cVariable
    ! Multigrid cycle flavour requested ('V','W','F',...)
    CHARACTER*1 :: CycleType
    ! Selected prolongation strategy name passed down to setup routines
    CHARACTER*10 :: MGProlongation
    ! Points to flag array that enables prolongation/restriction per level
    LOGICAL, POINTER :: bProlRest
    ! Dirichlet markers for velocity components on each grid level
    INTEGER, DIMENSION(:), POINTER::  KNPRU,KNPRV,KNPRW
    ! Dirichlet markers for pressure / scalar unknowns per level
    INTEGER, DIMENSION(:), POINTER::  KNPR
    ! Element-wise masks used to skip constrained cells in transfer ops
    TYPE(mg_kVector), DIMENSION(:), POINTER::  KNPRP
    ! Level-wise solution, increments, residuals, rhs, and smoother buffers for velocity
    TYPE(mg_dVector), DIMENSION(:), POINTER::  X_u,dX_u,D_u,A_u,B_u
    ! Level-wise solution, increments, residuals, rhs, and smoother buffers for pressure
    TYPE(mg_dVector), DIMENSION(:), POINTER::  X_p,dX_p,D_p,A_p,B_p
    ! Coupled Q2/P1 block matrices (B and BT operators) on each level
    TYPE (mg_Matrix), DIMENSION(:), POINTER :: BX,BY,BZ,BTX,BTY,BTZ
    ! Mixed velocity-pressure coupling matrices assembled per level
    TYPE(TMatrix), DIMENSION(:),  POINTER   :: Lq,Lql,Llq
    ! Generic mg vectors for field solves (solution, defect, auxiliary, rhs)
    TYPE(mg_dVector), DIMENSION(:), POINTER::  X,D,AUX,B
    ! Primary system matrices and their preconditioned variants across levels
    TYPE (mg_Matrix), DIMENSION(:), POINTER :: A,AP
    ! Split tensor blocks storing the directional stiffness contributions
    TYPE (mg_Matrix), DIMENSION(:), POINTER :: A11,A22,A33,A12,A13,A23,A21,A31,A32
    ! Field-specific block matrices (e.g. visco-elastic extra fields)
    TYPE (tMGFldMatrix), DIMENSION(:), POINTER :: AXX
    ! Constraint matrices and pressure blocks used by coarse solvers
    TYPE(TMatrix), DIMENSION(:),  POINTER :: L,LP
    ! Diagonal or smoothing workspace shared by stationary smoothers
    REAL*8  , DIMENSION(:)  , POINTER :: XP
    ! Book-keeping for hierarchy layout and smoothing parameters
    INTEGER :: MinLev,MaxLev,MedLev,MaxDifLev,MinIterCycle,MaxIterCycle,nIterCoarse,&
               & nSmootherSteps,CrsSolverType,SmootherType
    ! Number of coupled fields and sub-systems handled by this multigrid run
    INTEGER :: nOfFields,nOfSubsystemEqs
    ! Vanka smoother configuration flag (0 disables, >0 selects variant)
    integer :: vanka
    ! Measured defect reduction on the coarse grid during the current cycle
    REAL*8  :: DefImprCoarse
    ! Convergence criteria and relaxation factors steering the cycles
    REAL*8  :: Criterion1,Criterion2,RLX,CrsRelaxPrm,CrsRelaxParPrm
    ! Spectral radius estimates and defect history for convergence reports
    REAL*8  :: RhoMG1,RhoMG2,DefInitial,DefFinal
    ! Keeps track of iteration counts used in the most recent cycle
    INTEGER :: UsedIterCycle
  END TYPE tMultiGrid
  ! Primary multigrid controller for the coupled Q2/P1 system
  TYPE (tMultiGrid) :: myMG
  TYPE tMultiGrid_cc
    CHARACTER*10 :: cVariable
    CHARACTER*1  :: CycleType
    CHARACTER*10 :: MGProlongation
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
    INTEGER :: MinLev,MaxLev,MedLev,MinIterCycle,MaxIterCycle,nIterCoarse,nSmootherSteps,VANKA
    REAL*8  :: DefImprCoarse
    REAL*8  :: Criterion1,Criterion2,RLX
    REAL*8  :: RhoMG1,RhoMG2,DefInitial,DefFinal
    INTEGER :: UsedIterCycle
  END TYPE tMultiGrid_cc
  ! Multigrid descriptors for the coarse-correction (cc) subsystem
  TYPE (tMultiGrid_cc) :: myMG_cc
  ! Cell-wise material identifiers stored on each multigrid level
  TYPE(mg_kVector), ALLOCATABLE :: MaterialDistribution(:)
  ! Cell-centered density and diffusion coefficients per multigrid level
  TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: mgDensity(:),mgDiffCoeff(:)
  ! Normalised shear stress field cached per level for post-processing
  TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: mgNormShearStress(:)
  ! Hierarchical mesh description shared by multigrid routines
  type(tMultiMesh),save :: mg_mesh
  ! Marks DOFs that require MPI synchronisation during multigrid sweeps
  INTEGER, ALLOCATABLE :: ParKNPR(:)
  ! FictKNPR stores FBM particle ids; MixerKnpr marks mixer-influenced DOFs
  INTEGER, ALLOCATABLE :: FictKNPR(:),MixerKnpr(:)
  ! Packed 64-bit representation of FictKNPR entries for FBM lookups
  type(tUint64), allocatable :: FictKNPR_uint64(:)

  ! ============================================================================
  ! Vertex Cache for KVEL-based Force Acceleration
  ! ============================================================================
  TYPE tVertexCache
    INTEGER :: nVertices                    ! Number of cached DOF indices
    INTEGER, ALLOCATABLE :: dofIndices(:)   ! DOF indices (1..ndof)
    INTEGER :: particleID                   ! Particle index for debugging
  END TYPE tVertexCache

  TYPE(tVertexCache), ALLOCATABLE :: ParticleVertexCache(:)

  ! Runtime acceleration control flags (read from q2p1_param.dat)
  LOGICAL :: bUseHashGridAccel = .TRUE.
  LOGICAL :: bUseKVEL_Accel = .TRUE.

  ! Statistics for performance monitoring
  TYPE tKVEL_Stats
    INTEGER :: nCachedVertices
    INTEGER :: nCandidateElements
    INTEGER :: nBoundaryElements
    REAL*8  :: cacheTime, candidateTime, forceTime
  END TYPE tKVEL_Stats
  TYPE(tKVEL_Stats) :: myKVEL_Stats

  TYPE tElementMatrix
    REAL*8 :: A(85,85)
    INTEGER :: H(2)
  END TYPE
  TYPE tMGElementMatrix
    TYPE(tElementMatrix), ALLOCATABLE :: E(:)
  END TYPE
  ! Element-local coarse-correction matrices (UMFPACK factorizations)
  TYPE(tmgElementMatrix), ALLOCATABLE :: CC_EMat(:)
  ! CRS sparsity structure of the coarse-correction system matrix
  TYPE(TMatrix)          :: CC_crs_lMat
  ! Numeric values for the coarse-correction CRS matrix
  REAL*8 , ALLOCATABLE   :: CC_crs_AMat(:)
  ! UMFPACK handles for the coarse-correction matrix factorisation
  INTEGER :: CC_H(2)
  TYPE crs_e013_map
    INTEGER  :: ndof,cc_ndof
    INTEGER , allocatable :: ind(:)
    INTEGER , allocatable :: indE(:)
    REAL*8, allocatable   :: dBuffer(:)
  END TYPE crs_e013_map
  ! Mapping between E013 element blocks and global CRS numbering
  TYPE(crs_e013_map), ALLOCATABLE :: my_crs_e013_map(:)
  TYPE CC_Elem
    INTEGER, ALLOCATABLE :: pairE(:,:),pairV(:,:)
    INTEGER, ALLOCATABLE :: E_qq(:),E_lq(:),E_ql(:)
    REAL*8 , ALLOCATABLE :: a(:)
    INTEGER :: sym,num
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
  ! Cached coarse patches for assembling cc multigrid operators
  TYPE(mg_CCPiece),ALLOCATABLE :: my_mg_CCPiece(:)
  ! Global numbering used while assembling coarse-correction matrices
  INTEGER, ALLOCATABLE :: GlobalNumbering(:)
  TYPE tCoarseMat
    INTEGER :: na,nu
    INTEGER, dimension(:), allocatable   :: Row
    INTEGER, dimension(:), allocatable   :: Col
    REAL*8 , ALLOCATABLE   :: A(:)
    REAL*8 , ALLOCATABLE   :: D(:)
  END TYPE tCoarseMat
  ! Sparse coarse-grid matrix storage used by the coarse-correction solve
  TYPE(tCoarseMat), target :: myCrsMat

  ! Material & mesh
  ! Mesh descriptions, material fields, and geometry updates belong here.
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
  TYPE (tCGALObjects) :: myHeatObjects
  TYPE tParticleFBM
    CHARACTER cTYPE*10
    REAL*8 :: sizes(20),density
    REAL*8 :: ResistanceForce(3),TorqueForce(3)
    REAL*8 :: Position(3),Velocity(3),Angle(3),AngularVelocity(3)
    REAL*8 :: Acceleration(3),FrameVelocity(3)
  END TYPE tParticleFBM
  TYPE tFBM
    INTEGER :: nParticles
    REAL*8,ALLOCATABLE :: Force(:)
    TYPE (tParticleFBM), ALLOCATABLE :: ParticleOld(:),ParticleNew(:)
    integer, allocatable, dimension(:) :: iel_ug
    REAL*8, ALLOCATABLE :: ParticleRe(:)
    REAL*8, ALLOCATABLE :: ParticleSt(:)
  END TYPE tFBM
  CHARACTER cFBM_File*30
  TYPE (tFBM) :: myFBM
  LOGICAL, ALLOCATABLE :: BndrForce(:)
  REAL*8, ALLOCATABLE :: myQ2Coor(:,:), BoundaryNormal(:,:)
  TYPE tBoundary
    LOGICAL, ALLOCATABLE :: nWall(:)   ,nInflown(:),  nOutflow(:)   ,nSymmetry(:)
    LOGICAL, ALLOCATABLE :: bWall(:),bOutflow(:),bSymmetry(:,:),bSlip(:)
    INTEGER, ALLOCATABLE :: iInflow(:), iPhase(:),iTemperature(:)
    INTEGER, ALLOCATABLE :: LS_zero(:)
    LOGICAL, ALLOCATABLE :: bDisp_DBC(:)
  END TYPE tBoundary
  TYPE (tBoundary) :: myBoundary
  TYPE tDump
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Elements,Vertices
  END TYPE tDump
  TYPE(tDump) :: myDump
  REAL*8, ALLOCATABLE :: dPeriodicVector(:)
  REAL*8, allocatable :: mySegmentIndicator(:,:)
  TYPE tTriangle
    REAL*8 :: C(3,9)
  END TYPE tTriangle
  TYPE tTriangulation
    INTEGER :: nT
    TYPE(tTriangle), ALLOCATABLE :: T(:)
  END TYPE tTriangulation
  TYPE(tTriangulation), ALLOCATABLE :: myTSurf(:)
  LOGICAL :: bMeshAdaptation = .FALSE.
  CHARACTER*100 :: cAdaptedMeshFile
  INTEGER :: nUmbrellaSteps,nInitUmbrellaSteps
  integer :: nUmbrellaStepsLvl(9) = (/0, 0, 0, 0, 0, 0, 0, 0, 0/)
  integer :: nMainUmbrellaSteps = 0
  REAL*8 :: dIntegralHeat = 0d0, dNozzlePosition
  TYPE tALE
    REAL*8, ALLOCATABLE :: Monitor(:)
    REAL*8, ALLOCATABLE :: OldCoor(:,:), NewCoor(:,:),MeshVelo(:,:), OrigCoor(:,:)
    REAL*8, ALLOCATABLE :: Q2coor_old(:,:)
    LOGICAL :: bUseFrameVelocity = .false.
    REAL*8 :: dFrameVelocity(3),dFrameVelocityChange(3)
  END TYPE tALE
  TYPE(tALE),save :: myALE
  TYPE tTetInterface
    INTEGER :: nG,nT
    REAL*8  , ALLOCATABLE :: X(:,:),Y(:,:)
    INTEGER , ALLOCATABLE :: L(:)
  END TYPE tTetInterface
  TYPE(tTetInterface) :: lInterface,gInterface
  REAL*8, allocatable :: myFracField(:),myDistance(:)
  TYPE(tBoundingBox), dimension(:), allocatable :: mgBoundingBox
  INTEGER :: nSubCoarseMesh

  ! Output & diagnostics
  ! Post-processing data structures and reporting handles live in this region.
  TYPE tExport
    INTEGER :: Level,LevelMax
    CHARACTER*(3) :: Format
    CHARACTER*(40),ALLOCATABLE, DIMENSION(:) :: Fields
    CHARACTER*(4) :: cFileName(2)=['res ','main']
  END TYPE tExport
  TYPE(tExport) :: myExport
  INTEGER :: iOutput=0
  CHARACTER cProjectGridFile*220,cGridFileName*200,cProjectFile*200,cProjectFolder*200,cProjectNumber*4
  CHARACTER*13 :: outfile="OutFile  .txt"
  type(tPostprocessingParams) :: postParams
  TYPE t1DOutput
    REAL*8, ALLOCATABLE :: dMean(:),dMin(:),dMax(:),dLoc(:)
    CHARACTER cName*20
  END TYPE t1DOutput
  TYPE(t1DOutput), TARGET :: my1DOut(11)
  REAL*8, ALLOCATABLE :: my1DIntervals(:,:),my1DWeight(:),my1DTorque(:,:),my1DForceX(:,:),my1DForceY(:,:)
  INTEGER :: my1DOut_nol
  TYPE tViscFunc
    REAL*8 :: shear_rate(66),visc(66)
  END TYPE tViscFunc
  TYPE(tViscFunc) :: myViscFunc
  TYPE(tProperties),save :: Properties
  TYPE (tParticleParam) :: myParticleParam
  TYPE (tErrorCodes) ::  myErrorCode

  ! Communication
  ! Synchronisation lists and hostname bookkeeping should be added here.
  TYPE(mg_kVector),ALLOCATABLE :: GlobalParallelList1(:),GlobalParallelList2(:)
  REAL*8, ALLOCATABLE          :: GlobalParallelBufferOut(:),GlobalParallelBufferIn(:)
  INTEGER, ALLOCATABLE         :: GlobalNList(:),GlobalNBuffer(:)
  TYPE(mg_kVector),ALLOCATABLE :: HlobalParallelList1(:),HlobalParallelList2(:),HlobalParallelList3(:)
  REAL*8, ALLOCATABLE          :: HlobalParallelBufferOut(:),HlobalParallelBufferIn(:)
  INTEGER, ALLOCATABLE         :: HlobalNList(:),HlobalNBuffer(:)
  TYPE(mg_kVector),ALLOCATABLE :: CommOrder(:)
  TYPE tRecursiveCommunication
    character(len=256), allocatable :: all_hostnames(:)
    character(len=256), allocatable :: unique_hostnames(:)
    integer, allocatable :: hostleaders(:),groupIDs(:)
    integer, allocatable :: hostgroup(:)
    integer :: myid,numnodes,NumHosts,myNodeGroup
    character(len=256) :: HostName
  END TYPE tRecursiveCommunication
  TYPE(tRecursiveCommunication) :: myRecComm

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
