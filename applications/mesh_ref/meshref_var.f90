Module MeshRefVar

INTEGER lTriOutputLevel,lVTUOutputLevel

CHARACTER*(200) :: cOutputFolder,cShortProjectFile
LOGICAL :: bA_MD=.false.
LOGICAL :: bPDE_MD=.false.
Logical :: bDefTensor = .true.
integer :: nOfMarkers=1
integer, allocatable, dimension(:) :: markerE,markerV

type RefinerMesh
 integer :: nOfElem=0,nOfVert=0,nUniquePoints,nUniqueElems
 integer :: patchID=-2
 integer :: selection=0
 type(RefinerMesh), allocatable :: myRF(:)
 real*8, allocatable :: dcoor(:,:),dUniquedCoor(:,:)
 integer, allocatable :: knpr(:),kUniqueElem(:,:),PointerToMerged(:)
 integer, allocatable :: kvert(:,:)
 real*8 :: dTol = 1d8
end type RefinerMesh

type(RefinerMesh), allocatable :: myRF(:)

real*8, allocatable :: MergedMeshCoor(:,:)
integer, allocatable :: MergedMeshElem(:,:)
real*8, allocatable :: MergedMeshDist(:)

integer nUniquePoints,nUniqueElems

integer RROT(8,8),R2ROT(8),R3ROT(8)
data RROT/1,2,3,4,5,6,7,8,&
       2,3,4,1,6,7,8,5,&
       3,4,1,2,7,8,5,6,&
       4,1,2,3,8,5,6,7,&
       5,8,7,6,1,4,3,2,&
       6,5,8,7,2,1,4,3,&
       7,6,5,8,3,2,1,4,&
       8,7,6,5,4,3,2,1/
data R2ROT/1,5,6,2,4,8,7,3/
data R3ROT/1,4,8,5,2,3,7,6/

! character cPatches(0:16)*256,cSPatches(0:8)*256
integer :: nRefScheme = 6, nNonRefScheme = 2

integer myTemplate
logical templates(8,22)
data templates /1,0,0,0,0,0,0,0,&
                1,1,0,0,0,0,0,0,&
                1,0,1,0,0,0,0,0,&
                1,0,0,0,0,0,1,0,&
                1,1,1,0,0,0,0,0,&
                1,0,1,0,1,0,0,0,&
                1,0,1,0,0,1,0,0,&
                1,1,1,1,0,0,0,0,&
                1,1,1,0,1,0,0,0,&
                1,1,0,1,1,0,0,0,&
                1,0,1,1,1,0,0,0,&
                1,0,1,1,0,1,0,0,&
                1,0,1,0,1,0,1,0,&
                1,0,1,0,0,1,0,1,&
                1,1,1,1,1,0,0,0,&
                1,0,1,1,1,1,0,0,&
                1,1,0,1,1,0,1,0,&
                1,1,1,1,1,1,0,0,&
                1,1,1,1,1,0,1,0,&
                1,0,1,1,1,1,1,0,&
                1,1,1,1,1,1,1,0,&
                1,1,1,1,1,1,1,1/
                
character cTemplates(0:26)*256
data cTemplates /'_adc/PATCHES/MIXED/M-2.tri3d'    ,&  !Template:   1, #MarkedVerts:   1, Scheme:  F F F F F F F F  
                 '_adc/PATCHES/MIXED/M_Std.tri3d'  ,&  !Template:   1, #MarkedVerts:   1, Scheme:  T F F F F F F F  
                 '_adc/PATCHES/MIXED/M_Edge.tri3d' ,&  !Template:   2, #MarkedVerts:   2, Scheme:  T T F F F F F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:   3, #MarkedVerts:   2, Scheme:  T F T F F F F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:   4, #MarkedVerts:   2, Scheme:  T F F F F F T F  
                 '_adc/PATCHES/MIXED/M0_P0.tri3d'  ,&  !Template:   5, #MarkedVerts:   3, Scheme:  T T T F F F F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:   6, #MarkedVerts:   3, Scheme:  T F T F T F F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:   7, #MarkedVerts:   3, Scheme:  T F T F F T F F  
                 '_adc/PATCHES/MIXED/M_Face.tri3d' ,&  !Template:   8, #MarkedVerts:   4, Scheme:  T T T T F F F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:   9, #MarkedVerts:   4, Scheme:  T T T F T F F F  
                 '_adc/PATCHES/MIXED/M0_P3B.tri3d' ,&  !Template:  10, #MarkedVerts:   4, Scheme:  T T F T T F F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  11, #MarkedVerts:   4, Scheme:  T F T T T F F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  12, #MarkedVerts:   4, Scheme:  T F T T F T F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  13, #MarkedVerts:   4, Scheme:  T F T F T F T F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  14, #MarkedVerts:   4, Scheme:  T F T F F T F T  
                 '_adc/PATCHES/MIXED/M0_P2.tri3d'  ,&  !Template:  15, #MarkedVerts:   5, Scheme:  T T T T T F F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  16, #MarkedVerts:   5, Scheme:  T F T T T T F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  17, #MarkedVerts:   5, Scheme:  T T F T T F T F  
                 '_adc/PATCHES/MIXED/M0_LL.tri3d'  ,&  !Template:  18, #MarkedVerts:   6, Scheme:  T T T T T T F F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  19, #MarkedVerts:   6, Scheme:  T T T T T F T F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  20, #MarkedVerts:   6, Scheme:  T F T T T T T F  
                 '_adc/PATCHES/MIXED/M0_P3.tri3d'  ,&  !Template:  21, #MarkedVerts:   7, Scheme:  T T T T T T T F  
                 '_adc/PATCHES/MIXED/M0.tri3d'     ,&  !Template:  22, #MarkedVerts:   8, Scheme:  T T T T T T T T           
                 '_adc/PATCHES/MIXED/M0_Std.tri3d',&      
                 '_adc/PATCHES/MIXED/M0_Edge.tri3d',&      
                 '_adc/PATCHES/MIXED/M0.tri3d',&
                 '_adc/PATCHES/MIXED/M0_Face.tri3d'/

END Module MeshRefVar