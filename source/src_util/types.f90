module types
!-------------------------------------------------------------------------------------------------
! A module that contains several variants of the Laplacian 
! smoother 'Umbrella'. Besides the standard version of this 
! smoother type, versions with user-defined weighting functions
! are available
!-------------------------------------------------------------------------------------------------
! No implicit variables in this module
implicit none

type tBoundingBox
  double precision, dimension(3,2)  :: vertices
end type tBoundingBox

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

TYPE tParticle
 REAL*8 time
 REAL*8 coor(3)
 INTEGER indice
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


TYPE tParticleParam
 REAL*8 dEps1,dEps2, D_Out,D_in, f, Z_seed,Epsilon,hSize
 REAL*8 :: minFrac
 INTEGER :: nTimeLevels,nParticles,nRotation,Raster
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
 real*8, allocatable :: cutplanePositions(:)

 ! Now many segments do we take for coloring?
 real*8 :: numberSegments
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

TYPE tMeshInfoParticle
  real*8 xmin, xmax, ymin, ymax, zmin, zmax
END TYPE tMeshInfoParticle
type(tMeshInfoParticle) :: myMeshInfo


end module types
