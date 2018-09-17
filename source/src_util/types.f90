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
 LOGICAL :: ParamTypes(4)=.FALSE.
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

end module types
