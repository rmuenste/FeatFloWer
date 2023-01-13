Module MeshRefVar

INTEGER lTriOutputLevel,lVTUOutputLevel
REAL*8 RefinementThickness,MeshOutputScaleFactor

CHARACTER*(200) :: cIntputFolder=" ",cOutputFolder=" ",cShortProjectFile=" ",cReducedGridFile="ReducedMesh.tri", cReducedMeshdFolder="NEW_meshDir"
LOGICAL :: bA_MD=.false.
LOGICAL :: bPDE_MD=.false.
Logical :: bDefTensor = .true.
integer :: nOfMarkers=1,initfield
integer, allocatable, dimension(:) :: markerE,markerV

type RefinerMesh
 integer :: nOfElem=0,nOfVert=0,nUniquePoints,nUniqueElems
 integer :: patchID=-2
 integer :: selection=0
 type(RefinerMesh), allocatable :: myRF(:)
 real*8, allocatable :: dcoor(:,:),dUniquedCoor(:,:)
 integer, allocatable :: kUniqueKnpr(:)
 integer, allocatable :: knpr(:),kUniqueElem(:,:),PointerToMerged(:)
 integer, allocatable :: kvert(:,:)
 real*8 :: dTol = 1d8
end type RefinerMesh

type(RefinerMesh), allocatable :: myRF(:)

type tParList
 Logical, allocatable :: Wall(:),SideWall(:,:),Inflow(:,:),Outflow(:)
end type tParList
type (tParList) :: ParList,ParCleanList

real*8 OverallBoundingBox(3,2)
real*8 , allocatable :: AreaIntensity(:,:)
real*8, allocatable :: MergedMeshCoor(:,:),ReducedMeshCoor(:,:),ReducedCleanMeshCoor(:,:)
integer, allocatable :: MergedMeshElem(:,:),MergedMeshknpr(:),ReducedMeshElem(:,:),ReducedCleanMeshElem(:,:),ReducedMeshBC(:)
real*8, allocatable :: MergedMeshDist(:)
integer, allocatable :: level(:)

integer nUniquePoints,nUniqueElems,nReducedElems,nReducedPoints,nReducedCleanPoints,nReducedCleanElems

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
integer :: nRefScheme , nNonRefScheme 

integer myTemplate
logical templates(8,22)
data templates /.TRUE. ,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE. ,.FALSE.,&
                .TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.FALSE.,.TRUE. ,.FALSE.,.FALSE.,&
                .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,&
                .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,.FALSE.,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,.FALSE.,&
                .TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,&
                .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,.FALSE.,&
                .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,&
                .TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,&
                .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE.,&
                .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. /
                
character cTemplatesX(0:22)*256
data cTemplatesX /'PATCHES_C/X0.tri3d',&       !Template:   0, #MarkedVerts:   0, Scheme:  F F F F F F F F  
                 'PATCHES_C/X1.tri3d' ,&       !Template:   1, #MarkedVerts:   1, Scheme:  T F F F F F F F  
                 'PATCHES_C/x2.tri3d' ,&       !Template:   2, #MarkedVerts:   2, Scheme:  T T F F F F F F  
                 'PATCHES_C/X3.tri3d' ,&       !Template:   3, #MarkedVerts:   2, Scheme:  T F T F F F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   4, #MarkedVerts:   2, Scheme:  T F F F F F T F  
                 'PATCHES_C/X5.tri3d' ,&       !Template:   5, #MarkedVerts:   3, Scheme:  T T T F F F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   6, #MarkedVerts:   3, Scheme:  T F T F T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   7, #MarkedVerts:   3, Scheme:  T F T F F T F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   8, #MarkedVerts:   4, Scheme:  T T T T F F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   9, #MarkedVerts:   4, Scheme:  T T T F T F F F  
                 'PATCHES_C/X10.tri3d' ,&       !Template:  10, #MarkedVerts:   4, Scheme:  T T F T T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  11, #MarkedVerts:   4, Scheme:  T F T T T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  12, #MarkedVerts:   4, Scheme:  T F T T F T F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  13, #MarkedVerts:   4, Scheme:  T F T F T F T F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  14, #MarkedVerts:   4, Scheme:  T F T F F T F T  
                 'PATCHES_C/V.tri3d' ,&       !Template:  15, #MarkedVerts:   5, Scheme:  T T T T T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  16, #MarkedVerts:   5, Scheme:  T F T T T T F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  17, #MarkedVerts:   5, Scheme:  T T F T T F T F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  18, #MarkedVerts:   6, Scheme:  T T T T T T F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  19, #MarkedVerts:   6, Scheme:  T T T T T F T F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  20, #MarkedVerts:   6, Scheme:  T F T T T T T F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  21, #MarkedVerts:   7, Scheme:  T T T T T T T F  
                 'PATCHES_C/V.tri3d'/         !Template:  22, #MarkedVerts:   8, Scheme:  T T T T T T T T           

                 character cTemplates(0:26)*256
data cTemplates /'PATCHES_C/M_Null.tri3d'    ,&  !Template:   0, #MarkedVerts:   0, Scheme:  F F F F F F F F  
                 'PATCHES_C/M_Std.tri3d'  ,&  !Template:   1, #MarkedVerts:   1, Scheme:  T F F F F F F F  
                 'PATCHES_C/M_Edge.tri3d' ,&  !Template:   2, #MarkedVerts:   2, Scheme:  T T F F F F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   3, #MarkedVerts:   2, Scheme:  T F T F F F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   4, #MarkedVerts:   2, Scheme:  T F F F F F T F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   5, #MarkedVerts:   3, Scheme:  T T T F F F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   6, #MarkedVerts:   3, Scheme:  T F T F T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   7, #MarkedVerts:   3, Scheme:  T F T F F T F F  
                 'PATCHES_C/M_Face.tri3d' ,&  !Template:   8, #MarkedVerts:   4, Scheme:  T T T T F F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:   9, #MarkedVerts:   4, Scheme:  T T T F T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  10, #MarkedVerts:   4, Scheme:  T T F T T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  11, #MarkedVerts:   4, Scheme:  T F T T T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  12, #MarkedVerts:   4, Scheme:  T F T T F T F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  13, #MarkedVerts:   4, Scheme:  T F T F T F T F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  14, #MarkedVerts:   4, Scheme:  T F T F F T F T  
                 'PATCHES_C/V.tri3d' ,&       !Template:  15, #MarkedVerts:   5, Scheme:  T T T T T F F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  16, #MarkedVerts:   5, Scheme:  T F T T T T F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  17, #MarkedVerts:   5, Scheme:  T T F T T F T F  
                 'PATCHES_C/18.tri3d' ,&       !Template:  18, #MarkedVerts:   6, Scheme:  T T T T T T F F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  19, #MarkedVerts:   6, Scheme:  T T T T T F T F  
                 'PATCHES_C/V.tri3d' ,&       !Template:  20, #MarkedVerts:   6, Scheme:  T F T T T T T F  
                 'PATCHES_C/21.tri3d' ,&       !Template:  21, #MarkedVerts:   7, Scheme:  T T T T T T T F  
                 'PATCHES_C/Full.tri3d' ,&       !Template:  22, #MarkedVerts:   8, Scheme:  T T T T T T T T           
                 'PATCHES_C/M_Std.tri3d' ,&      
                 'PATCHES_C/M_Edge.tri3d' ,&      
                 'PATCHES_C/M_Std.tri3d' ,&               
                 'PATCHES_C/M_Face.tri3d' /

character cTemplatesF(0:26)*256
data cTemplatesF /'PATCHES/M-2.tri3d'    ,&  !Template:   0, #MarkedVerts:  0, Scheme:  F F F F F F F F  
                 'PATCHES/M_Std.tri3d'  ,&  !Template:   1, #MarkedVerts:   1, Scheme:  T F F F F F F F  
                 'PATCHES/M_Edge.tri3d' ,&  !Template:   2, #MarkedVerts:   2, Scheme:  T T F F F F F F  
                 'PATCHES/M0.tri3d'     ,&  !Template:   3, #MarkedVerts:   2, Scheme:  T F T F F F F F  
                 'PATCHES/M0.tri3d'     ,&  !Template:   4, #MarkedVerts:   2, Scheme:  T F F F F F T F  
                 'PATCHES/M0_P0.tri3d'  ,&  !Template:   5, #MarkedVerts:   3, Scheme:  T T T F F F F F  
                 'PATCHES/M0.tri3d'     ,&  !Template:   6, #MarkedVerts:   3, Scheme:  T F T F T F F F  
                 'PATCHES/M0.tri3d'     ,&  !Template:   7, #MarkedVerts:   3, Scheme:  T F T F F T F F  
                 'PATCHES/M_Face.tri3d' ,&  !Template:   8, #MarkedVerts:   4, Scheme:  T T T T F F F F  
                 'PATCHES/M0_P2B.tri3d' ,&  !Template:   9, #MarkedVerts:   4, Scheme:  T T T F T F F F  
                 'PATCHES/M0_P3B.tri3d' ,&  !Template:  10, #MarkedVerts:   4, Scheme:  T T F T T F F F  
                 'PATCHES/M0_P2C.tri3d' ,&  !Template:  11, #MarkedVerts:   4, Scheme:  T F T T T F F F  
                 'PATCHES/M0.tri3d'     ,&  !Template:  12, #MarkedVerts:   4, Scheme:  T F T T F T F F  
                 'PATCHES/M0_LL.tri3d'  ,&  !Template:  13, #MarkedVerts:   4, Scheme:  T F T F T F T F  
                 'PATCHES/M0.tri3d'     ,&  !Template:  14, #MarkedVerts:   4, Scheme:  T F T F F T F T  
                 'PATCHES/M0_P2.tri3d'  ,&  !Template:  15, #MarkedVerts:   5, Scheme:  T T T T T F F F  
                 'PATCHES/M0.tri3d'     ,&  !Template:  16, #MarkedVerts:   5, Scheme:  T F T T T T F F  
                 'PATCHES/M0.tri3d'     ,&  !Template:  17, #MarkedVerts:   5, Scheme:  T T F T T F T F  
                 'PATCHES/M0_LL.tri3d'  ,&  !Template:  18, #MarkedVerts:   6, Scheme:  T T T T T T F F  
                 'PATCHES/M0.tri3d'     ,&  !Template:  19, #MarkedVerts:   6, Scheme:  T T T T T F T F  
                 'PATCHES/M0_P5.tri3d'  ,&  !Template:  20, #MarkedVerts:   6, Scheme:  T F T T T T T F  
                 'PATCHES/M0_P3.tri3d'  ,&  !Template:  21, #MarkedVerts:   7, Scheme:  T T T T T T T F  
                 'PATCHES/M0.tri3d'     ,&  !Template:  22, #MarkedVerts:   8, Scheme:  T T T T T T T T           
                 'PATCHES/M0_Std.tri3d',&      
                 'PATCHES/M0_Edge.tri3d',&      
                 'PATCHES/M0.tri3d',&
                 'PATCHES/M0_Face.tri3d'/
END Module MeshRefVar