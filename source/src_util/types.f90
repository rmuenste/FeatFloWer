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
 REAL*8 coor(3)
 INTEGER indice
 INTEGER :: ID=1
END TYPE tParticle

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

TYPE tParticleParam
 REAL*8 dEps1,dEps2, D_Out,D_in, f, Z_seed,Epsilon,hSize,d_CorrDist
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

 LOGICAL :: bBacktrace=.false.
 
 INTEGER :: NumberOfInflowRegions=0
 TYPE(tParticleInflow), ALLOCATABLE :: InflowRegion(:)
 
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

end module types
