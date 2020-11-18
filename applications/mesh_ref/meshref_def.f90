MODULE MeshRefDef
USE types, ONLY: tMultiMesh,tMesh
USE Parametrization, ONLY: InitParametrization_STRCT,ParametrizeBndryPoints_STRCT,&
    DeterminePointParametrization_STRCT,ProlongateParametrization_STRCT,myParBndr,nBnds
USE var_QuadScalar
USE PP3D_MPI, ONLY:myid,master
USE MESH_Structures
use iniparser

implicit none

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

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CreateRefinedMesh
integer i,j,k,l,ii,ivt,jj
real*8 E(3,8)
logical bVert(8),bVertCopy(8)
integer iVert(8)
logical bFound,bUnique
logical, allocatable :: mytemp(:,:)
integer, allocatable :: mytempCount(:)

integer iCombination, nOfCombinations,nAdditionalElems

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nOfCombinations = 0
allocate(mytemp(8,22),mytempCount(22))
mytempCount = 0

11 continue

nAdditionalElems = 0
do i=1,nel
 if (markerE(i).eq.3) then
 
   bFound = .false.
   bVert = .FALSE.
   DO j=1,8
    ivt = mg_mesh%level(ilev)%kvert(j,i)
    do k=1,mg_mesh%level(ilev)%nvel
     l = mg_mesh%level(ilev)%kvel(k,ivt)
     if (l.ne.0) then
      if (MarkerE(l).eq.1.or.MarkerE(l).eq.2) then
       bVert(j) = .TRUE.
      end if
     end if
    end do
   END DO
  
   CALL RotatePatch(bVert,iVert)
   CALL DetermineTemplate(bVert,myTemplate)
   
   if (myTemplate.eq.1.or.&
       myTemplate.eq.2.or.&
       myTemplate.eq.5.or.&
       myTemplate.eq.8.or.&
       myTemplate.eq.10.or.&
       myTemplate.eq.15.or.&
       myTemplate.eq.18.or.&
       myTemplate.eq.21.or.&
       myTemplate.eq.22) then
    else
     markerE(i)= 1
     
     do j=1,8
      ivt = mg_mesh%level(ilev)%kvert(j,i)
      do k=1,mg_mesh%level(ilev)%nvel
       l = mg_mesh%level(ilev)%kvel(k,ivt)
       if (l.ne.0) then
        if (MarkerE(l).ne.1) MarkerE(l) = 3
       end if
      end do
     end do
     
     nAdditionalElems = nAdditionalElems + 1
     
    end if
    
   
 end if 
end do

write(*,*) nAdditionalElems, " additional element has been set to full refinement"

if (nAdditionalElems.gt.0) goto 11

allocate(myRF(nel))

do i=1,nel

 if (markerE(i).eq.0) then
 
   myRF(i)%nOfElem = nNonRefScheme**3
   myRF(i)%nOfVert = (nNonRefScheme+1)**3
   allocate(myRF(i)%kvert(8,myRF(i)%nOfElem))
   allocate(myRF(i)%dcoor(3,myRF(i)%nOfVert))
   allocate(myRF(i)%knpr(myRF(i)%nOfVert))
   myRF(i)%kvert=0
   do j=1,8
    E(:,j) = mg_mesh%level(ilev)%dcorvg(:,mg_mesh%level(ilev)%kvert(j,i))
   end do
   
   Call FillUpRefined(myRF(i),E,nNonRefScheme+1)
   myRF(i)%patchID = -1
 end if

 if (markerE(i).eq.1.or.markerE(i).eq.2) then
 
   myRF(i)%nOfElem = nRefScheme**3
   myRF(i)%nOfVert = (nRefScheme+1)**3
   allocate(myRF(i)%kvert(8,myRF(i)%nOfElem))
   allocate(myRF(i)%dcoor(3,myRF(i)%nOfVert))
   allocate(myRF(i)%knpr(myRF(i)%nOfVert))
   myRF(i)%knpr=0
   do j=1,8
    E(:,j) = mg_mesh%level(ilev)%dcorvg(:,mg_mesh%level(ilev)%kvert(j,i))
   end do
   
   Call FillUpRefined(myRF(i),E,nRefScheme+1)
   myRF(i)%patchID = 0

 end if

 if (markerE(i).eq.3) then
 
   bFound = .false.
   bVert = .FALSE.
   DO j=1,8
    ivt = mg_mesh%level(ilev)%kvert(j,i)
    do k=1,mg_mesh%level(ilev)%nvel
     l = mg_mesh%level(ilev)%kvel(k,ivt)
     if (l.ne.0) then
      if (MarkerE(l).eq.1.or.MarkerE(l).eq.2) then
       bVert(j) = .TRUE.
      end if
     end if
    end do
   END DO

   ii = 0
   do j=1,8
    if (bVert(j))  ii = ii + 1
   end  do
   
   bVertCopy = bVert
   CALL RotatePatch(bVert,iVert)
   CALL DetermineTemplate(bVert,myTemplate)
   bVert = bVertCopy
   
   if (myTemplate.eq.1.or.myTemplate.eq.2.or.myTemplate.eq.8) then
   
    iVert = mg_mesh%level(ilev)%kvert(:,i)

    CALL RotatePatch(bVert,iVert)
   
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    
    CALL FillUpRefinedElement(myRF(i),E,ii,bvert,.false.)
    myRF(i)%patchID = ii
    bFound=.true.
!     write(*,*) myTemplate,ii
    
   end if

   if (ii.eq.8) then
   
    myRF(i)%nOfElem = nRefScheme**3
    myRF(i)%nOfVert = (nRefScheme+1)**3
    allocate(myRF(i)%kvert(8,myRF(i)%nOfElem))
    allocate(myRF(i)%dcoor(3,myRF(i)%nOfVert))
    allocate(myRF(i)%knpr(myRF(i)%nOfVert))
    myRF(i)%knpr=0
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,mg_mesh%level(ilev)%kvert(j,i))
    end do
    
    Call FillUpRefined(myRF(i),E,nRefScheme+1)
    myRF(i)%patchID = ii
    bFound=.true.
    
   end if

   
   if (myTemplate.eq.5.or.myTemplate.eq.10.or.myTemplate.eq.15.or.myTemplate.eq.18.or.myTemplate.eq.21) then
!    if ((ii.eq.3.or.ii.eq.5.or.ii.eq.6.or.ii.eq.7)) then
   
    iVert = mg_mesh%level(ilev)%kvert(:,i)

!      write(*,'(8L,8((", "),I0))') bVert,iVert
     CALL RotatePatch(bVert,iVert)
!      write(*,'(8L,8((", "),I0))') bVert,iVert
   
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    
    CALL FillUpRefinedElement(myRF(i),E,ii,bvert,.true.)
    myRF(i)%patchID = ii
    bFound=.true.
    
    !!!!!!!!!!! RECURSIVE REFINEMENT !!!!!!!!!!!!!!!
    allocate(myRF(i)%myRF(myRF(i)%nOfElem))
    
    DO iel=1,myRF(i)%nOfElem
     
     DO j=1,8
      ivt = myRF(i)%kvert(j,iel)
      iVert(j) = ivt
!       write(*,*) ivt
      if (myRF(i)%knpr(ivt).eq.1) then
       bVert(j)=.true.
      else
       bVert(j)=.false.
      endif
     END DO
     
     jj = 0
     do j=1,8
      if (bVert(j))  jj = jj + 1
     end  do
     
     CALL RotatePatch(bVert,iVert)
!      write(*,'(4(I0,(": ")),8L,8((", "),I0))') jj,iel,myRF(i)%nOfElem,ii,bVert,iVert
     
     do j=1,8
      E(:,j) = myRF(i)%dcoor(:,iVert(j))
     end do
     
     if (jj.gt.0) then
      myTemplate=22+jj
      CALL FillUpRefinedElement(myRF(i)%myRF(iel),E,jj+8,bvert,.false.)
      myRF(i)%myRF(iel)%selection = jj
      myRF(i)%myRF(iel)%patchID = ii
!       write(*,*) 'jj=',jj, bVert,cPatches(jj+8),cSPatches(jj+8)
      
     else
!       write(*,*) 'jj=',jj
      myRF(i)%myRF(iel)%nOfElem = 1
      myRF(i)%myRF(iel)%nOfVert = 8
      allocate(myRF(i)%myRF(iel)%kvert(8,myRF(i)%myRF(iel)%nOfElem))
      allocate(myRF(i)%myRF(iel)%dcoor(3,myRF(i)%myRF(iel)%nOfVert))
      allocate(myRF(i)%myRF(iel)%knpr(myRF(i)%myRF(iel)%nOfVert))
      myRF(i)%myRF(iel)%kvert=0
      
      Call FillUpRefined(myRF(i)%myRF(iel),E,1+1)
      myRF(i)%myRF(iel)%patchID = ii
      myRF(i)%myRF(iel)%selection = jj
     
     end if
    end do
    
   end if
   
   if (.not.bFound) then
     
    iVert = mg_mesh%level(ilev)%kvert(:,i)
    CALL RotatePatch(bVert,iVert)
    
   
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    CALL FillUpRefinedElement(myRF(i),E,0,bVert,.false.)
    
    myRF(i)%patchID = -2

   end if

   bVert = bVertCopy
   iVert = 0
   CALL RotatePatch(bVert,iVert)

   bUnique = .true.
   DO iCombination = 1,nOfCombinations
    CALL ComparePatch(mytemp(:,iCombination),bVert,bUnique)
    if (.not.bUnique) exit
   END DO
   
   if (bUnique) then
    nOfCombinations = nOfCombinations+1
    mytemp(:,nOfCombinations) = bVert
    mytempCount(iCombination) = 1
   else
    mytempCount(iCombination) = mytempCount(iCombination) + 1
   end if
   
 end if

end do

 DO iCombination = 1,nOfCombinations
   j = 0
   do i=1,8
    if (mytemp(i,iCombination)) j = j + 1
   end do
   
   write(*,*) j,' :: ', mytemp(:,iCombination),mytempCount(iCombination)
   
 END DO

 END SUBROUTINE CreateRefinedMesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE InitMarker
logical bM
integer i,j,k,l,m,ivt
real*8 myRand 

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

allocate(markerE(nel))
allocate(markerV(nvt))

markerE = 0
markerV = 0

! j=0
! DO i=1,nel
!  CALL RANDOM_NUMBER(myRand)
!  if (myRand.gt.0.95d0) then
!   j = j +1
!   markerE(i) = 1
!  end if
! end do

CALL Initfield3(markerE,mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%dcorvg,nel)

END SUBROUTINE InitMarker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetUpMarker
logical bM
integer i,j,k,l,m,ivt
integer istep

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

!!! First the Marker Corner poionts of the 1st level have to be created
DO i=1,nel
 if (markerE(i).Eq.1) then
  do k=1,8
   j = mg_mesh%level(ilev)%kvert(k,i)
   markerV(j) = 1
  end do
 end if
end do

!!! Construction of the 2nd level Markers
DO i=1,nel
 if (markerE(i).Eq.1) then
  do j=1,8
   ivt = mg_mesh%level(ilev)%kvert(j,i)
   do k=1,mg_mesh%level(ilev)%nvel
    l = mg_mesh%level(ilev)%kvel(k,ivt)
    if (l.ne.0) then
     if (MarkerE(l).ne.1) MarkerE(l) = 3
    end if
   end do
  end do
 end if
end do


END SUBROUTINE SetUpMarker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GetParameters
CHARACTER*(10) :: cParam
integer ilong

!  CHARACTER*(200) :: cmd
 LOGICAL bExist

 OPEN(1,file='param.cfg')

 READ(1,*) cProjectFolder
 iLong = LEN(ADJUSTL(TRIM(cProjectFolder)))+1
 WRITE(cProjectFolder(iLong:),"(A)") "/"
 WRITE(*,*) adjustl(trim(cProjectFolder))
 
 READ(1,*) cShortProjectFile
 cProjectFile = adjustl(trim(cProjectFolder))//adjustl(trim(cShortProjectFile))
 WRITE(*,*) adjustl(trim(cProjectFile))
 
 CALL ExtractMeshfile()
 
! READ(1,*) mg_Mesh%nlmax
 mg_Mesh%nlmax = 1
 mg_Mesh%nlmin = 1
 WRITE(*,*) 'Min and Max levels: ', mg_Mesh%nlmin,mg_Mesh%nlmax
 mg_Mesh%nlmin = 1
 READ(1,*) dCGALtoRealFactor
 WRITE(*,'(A,ES12.4)') ' CGAL Scaling factor: ', dCGALtoRealFactor
 READ(1,*) cOutputFolder
 WRITE(*,*) 'Output Folder: "'//adjustl(trim(cOutputFolder))//'"'
 READ(1,*) lTriOutputLevel
 WRITE(*,*) 'Outputlevel for the ".tri" file: ', lTriOutputLevel
 READ(1,*) lVTUOutputLevel
 WRITE(*,*) 'Outputlevel for the ".vtu" file: ', lVTUOutputLevel

 CLOSE(1)

#ifdef __GFORTRAN__
 INQUIRE(file=adjustl(trim(cOutputFolder)),EXIST=bExist)
#elif defined __INTEL_COMPILER
 INQUIRE(directory=adjustl(trim(cOutputFolder)),EXIST=bExist)
#else
 INQUIRE(file=adjustl(trim(cOutputFolder)),EXIST=bExist)
#endif

 if (.not.bExist) then
!   cmd='mkdir '//adjustl(trim(cOutputFolder))
  CALL system('mkdir '//adjustl(trim(cOutputFolder)))
 end if

END SUBROUTINE GetParameters
!----------------------------------------------------------
SUBROUTINE Output_VTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER ive,ivt,ioffset
INTEGER :: iunit=123
CHARACTER*(100) filename

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

write(*,*) 'nel=',nel

filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/mesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nvt,""" NumberOfCells=""",nel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerV",""" format=""ascii"">"
 do ivt=1,nvt
  write(iunit, '(A,I10)')"        ",markerV(ivt)
 end do
 write(iunit, *)"        </DataArray>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
 do ivt=1,nvt
  write(iunit, '(A,E16.7)')"        ",REAL(ivt)
 end do
 write(iunit, *)"        </DataArray>"

 write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do ivt=1,nvt
 write(iunit,'(A10,3E16.7)')"          ",REAL(mg_mesh%level(ilev)%dcorvg(1,ivt)),REAL(mg_mesh%level(ilev)%dcorvg(2,ivt)),REAL(mg_mesh%level(ilev)%dcorvg(3,ivt))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"


write(iunit, '(A)')"    <CellData>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerE",""" format=""ascii"">"
do ivt=1,nel
 write(iunit, '(A,I10)')"        ",MarkerE(ivt)
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","OrigID",""" format=""ascii"">"
do ivt=1,nel
 write(iunit, '(A,I10)')"        ",ivt
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
do ive=1,nel   
 write(iunit, '(8I10)')mg_mesh%level(ilev)%kvert(1,ive)-1,mg_mesh%level(ilev)%kvert(2,ive)-1,mg_mesh%level(ilev)%kvert(3,ive)-1,mg_mesh%level(ilev)%kvert(4,ive)-1,&
                       mg_mesh%level(ilev)%kvert(5,ive)-1,mg_mesh%level(ilev)%kvert(6,ive)-1,mg_mesh%level(ilev)%kvert(7,ive)-1,mg_mesh%level(ilev)%kvert(8,ive)-1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_VTK
!----------------------------------------------------------
SUBROUTINE Output_RefVTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(100) filename


nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nnvt = 0
nnel = 0

DO iel=1,nel
 if (.not.allocated(myRF(iel)%myRF)) then
  nnel = nnel + myRF(iel)%nOfElem
  nnvt = nnvt + myRF(iel)%nOfVert
 else
  do jel =1,myRF(iel)%nOfElem
    nnel = nnel + myRF(iel)%myRF(jel)%nOfElem
    nnvt = nnvt + myRF(iel)%myRF(jel)%nOfVert
  end do
 end if
end do


filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/RefMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nnvt,""" NumberOfCells=""",nnel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

!  write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerV",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,I10)')"        ",markerV(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"
! 
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,E16.7)')"        ",REAL(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","knpr",""" format=""ascii"">"
do iel=1,nel
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfVert
    write(iunit,'(A10,I0)')"          ",myRF(iel)%knpr(i)
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfVert
     write(iunit,'(A10,I0)')"          ",myRF(iel)%myRF(jel)%knpr(i)
   end do
  end do
 end if  
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do iel=1,nel
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfVert
   write(iunit,'(A10,3E16.7)')"          ",REAL(myRF(iel)%dcoor(1,i)),REAL(myRF(iel)%dcoor(2,i)),REAL(myRF(iel)%dcoor(3,i))
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfVert
    write(iunit,'(A10,3E16.7)')"          ",REAL(myRF(iel)%myRF(jel)%dcoor(1,i)),REAL(myRF(iel)%myRF(jel)%dcoor(2,i)),REAL(myRF(iel)%myRF(jel)%dcoor(3,i))
   end do
  end do
 end if
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","NREFEL",""" format=""ascii"">"
do iel=1,nel
 n=myRF(iel)%nOfElem
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfElem
   write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfElem
    write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
   end do
  end do
 end if
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","PATCH_ID",""" format=""ascii"">"
do iel=1,nel
 n=myRF(iel)%nOfElem
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfElem
   write(iunit, '(A,I10)')"        ",myRF(iel)%patchID
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfElem
    write(iunit, '(A,I10)')"        ",myRF(iel)%myRF(jel)%patchID
   end do
  end do
 end if
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","Selection",""" format=""ascii"">"
do iel=1,nel
 n=myRF(iel)%nOfElem
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfElem
   write(iunit, '(A,I10)')"        ",myRF(iel)%selection
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfElem
    write(iunit, '(A,I10)')"        ",myRF(iel)%myRF(jel)%selection
   end do
  end do
 end if
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
nnvt = 0
do iel=1,nel
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfElem
   write(iunit, '(8I10)') nnvt + myRF(iel)%kvert(1,i) - 1, nnvt + myRF(iel)%kvert(2,i) - 1, nnvt + myRF(iel)%kvert(3,i) - 1, nnvt + myRF(iel)%kvert(4,i) - 1,&
                          nnvt + myRF(iel)%kvert(5,i) - 1, nnvt + myRF(iel)%kvert(6,i) - 1, nnvt + myRF(iel)%kvert(7,i) - 1, nnvt + myRF(iel)%kvert(8,i) - 1 
  end do
  nnvt = nnvt + myRF(iel)%nOfVert
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfElem
    write(iunit, '(8I10)') nnvt + myRF(iel)%myRF(jel)%kvert(1,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(2,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(3,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(4,i) - 1,&
                           nnvt + myRF(iel)%myRF(jel)%kvert(5,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(6,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(7,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(8,i) - 1 
   end do
   nnvt = nnvt + myRF(iel)%myRF(jel)%nOfVert
  end do
 end if
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nnel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nnel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_RefVTK
!----------------------------------------------------------
SUBROUTINE Output_UniqueRefVTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(100) filename


nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nnvt = 0
nnel = 0

DO iel=1,nel
  nnel = nnel + myRF(iel)%nUniqueElems
  nnvt = nnvt + myRF(iel)%nUniquePoints
end do


filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/UniqueRefMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nnvt,""" NumberOfCells=""",nnel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

!  write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerV",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,I10)')"        ",markerV(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"
! 
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,E16.7)')"        ",REAL(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","knpr",""" format=""ascii"">"
! do iel=1,nel
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfVert
!     write(iunit,'(A10,I0)')"          ",myRF(iel)%knpr(i)
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfVert
!      write(iunit,'(A10,I0)')"          ",myRF(iel)%myRF(jel)%knpr(i)
!    end do
!   end do
!  end if  
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do iel=1,nel
  do i=1,myRF(iel)%nUniquePoints
   write(iunit,'(A10,3E16.7)')"          ",REAL(myRF(iel)%dUniquedCoor(1,i)),REAL(myRF(iel)%dUniquedCoor(2,i)),REAL(myRF(iel)%dUniquedCoor(3,i))
  end do
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","NREFEL",""" format=""ascii"">"
! do iel=1,nel
!  n=myRF(iel)%nOfElem
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfElem
!    write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfElem
!     write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
!    end do
!   end do
!  end if
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","PATCH_ID",""" format=""ascii"">"
do iel=1,nel
  do i=1,myRF(iel)%nUniqueElems
   write(iunit, '(A,I10)')"        ",myRF(iel)%patchID
  end do
end do
write(iunit, *)"        </DataArray>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","Selection",""" format=""ascii"">"
! do iel=1,nel
!  n=myRF(iel)%nOfElem
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfElem
!    write(iunit, '(A,I10)')"        ",myRF(iel)%selection
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfElem
!     write(iunit, '(A,I10)')"        ",myRF(iel)%myRF(jel)%selection
!    end do
!   end do
!  end if
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
nnvt = 0
do iel=1,nel
  do i=1,myRF(iel)%nUniqueElems
   write(iunit, '(8I10)') nnvt + myRF(iel)%kUniqueElem(1,i) - 1, nnvt + myRF(iel)%kUniqueElem(2,i) - 1, nnvt + myRF(iel)%kUniqueElem(3,i) - 1, nnvt + myRF(iel)%kUniqueElem(4,i) - 1,&
                          nnvt + myRF(iel)%kUniqueElem(5,i) - 1, nnvt + myRF(iel)%kUniqueElem(6,i) - 1, nnvt + myRF(iel)%kUniqueElem(7,i) - 1, nnvt + myRF(iel)%kUniqueElem(8,i) - 1 
  end do
  nnvt = nnvt + myRF(iel)%nUniquePoints
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nnel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nnel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_UniqueRefVTK
!----------------------------------------------------------
SUBROUTINE Output_MergedRefVTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(100) filename

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nnel = nUniqueElems
nnvt = nUniquePoints

filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/MergedRefMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nnvt,""" NumberOfCells=""",nnel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

!  write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerV",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,I10)')"        ",markerV(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"
! 
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,E16.7)')"        ",REAL(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","knpr",""" format=""ascii"">"
! do iel=1,nel
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfVert
!     write(iunit,'(A10,I0)')"          ",myRF(iel)%knpr(i)
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfVert
!      write(iunit,'(A10,I0)')"          ",myRF(iel)%myRF(jel)%knpr(i)
!    end do
!   end do
!  end if  
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do i=1,nnvt
   write(iunit,'(A10,3E16.7)')"          ",REAL(MergedMeshCoor(1,i)),REAL(MergedMeshCoor(2,i)),REAL(MergedMeshCoor(3,i))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","NREFEL",""" format=""ascii"">"
! do iel=1,nel
!  n=myRF(iel)%nOfElem
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfElem
!    write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfElem
!     write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
!    end do
!   end do
!  end if
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","PATCH_ID",""" format=""ascii"">"
do iel=1,nel
  do i=1,myRF(iel)%nUniqueElems
   write(iunit, '(A,I10)')"        ",myRF(iel)%patchID 
  end do
end do
write(iunit, *)"        </DataArray>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","Selection",""" format=""ascii"">"
! do iel=1,nel
!  n=myRF(iel)%nOfElem
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfElem
!    write(iunit, '(A,I10)')"        ",myRF(iel)%selection
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfElem
!     write(iunit, '(A,I10)')"        ",myRF(iel)%myRF(jel)%selection
!    end do
!   end do
!  end if
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
do iel=1,nnvt
   write(iunit, '(8I10)') MergedMeshElem(1,iel) - 1,MergedMeshElem(2,iel) - 1,MergedMeshElem(3,iel) - 1,MergedMeshElem(4,iel) - 1,&
                          MergedMeshElem(5,iel) - 1,MergedMeshElem(6,iel) - 1,MergedMeshElem(7,iel) - 1,MergedMeshElem(8,iel) - 1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nnel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nnel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_MergedRefVTK
!----------------------------------------------------------
SUBROUTINE ExtractMeshfile()
IMPLICIT NONE
INTEGER LenStr,iEnd
CHARACTER(LEN=200) :: string,cFile
logical :: bFound=.false.

OPEN(unit=2,file=adjustl(trim(cProjectFile)))

DO
 READ(2,FMT='(200A)',IOSTAT=iEnd) string
 IF (iEnd.EQ.-1) EXIT
 LenStr = LEN(ADJUSTL(TRIM(string)))
 IF (LenStr.gt.4) THEN
  cFile = ADJUSTL(TRIM(string))
  IF (cFile(LenStr-3:LenStr).EQ.".tri") THEN
   cProjectGridFile = adjustl(trim(cFile))
!    cProjectGridFile = adjustl(trim(cProjectFolder))//'/'//adjustl(trim(cFile))
   Write(*,*) 'Mesh file: "'//ADJUSTL(TRIM(cProjectGridFile))//'"'
   bFound=.true.
   exit
  END IF
 END IF
END DO

CLOSE(2)
if (.not.bFound) then
 Write(*,*) 'Mesh file was NOT found in Project file! '
end if

END SUBROUTINE ExtractMeshfile
! ----------------------------------------------
SUBROUTINE Output_RefTriMesh()
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
INTEGER i,j,iBnds,nnvt,nnel
CHARACTER cf*(256)

nnvt = 0
nnel = 0

DO i=1,nel
 nnel = nnel + myRF(i)%nOfElem
 nnvt = nnvt + myRF(i)%nOfVert
end do

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/Ref'//adjustl(trim(cProjectGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: Ref'"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
WRITE(1,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(1,'(2I8,A)') nNEL,nNVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(1,'(A)') 'DCORVG'
do iel=1,nel
 do i=1,myRF(iel)%nOfVert
  write(1,'(A10,3E16.7)')"          ",REAL(myRF(iel)%dcoor(1,i)),REAL(myRF(iel)%dcoor(2,i)),REAL(myRF(iel)%dcoor(3,i))
 end do
end do

WRITE(1,'(A)') 'KVERT'
nnvt = 0
do iel=1,nel
 do i=1,myRF(iel)%nOfElem
  write(1, '(8I10)') nnvt + myRF(iel)%kvert(1,i), nnvt + myRF(iel)%kvert(2,i), nnvt + myRF(iel)%kvert(3,i), nnvt + myRF(iel)%kvert(4,i),&
                         nnvt + myRF(iel)%kvert(5,i), nnvt + myRF(iel)%kvert(6,i), nnvt + myRF(iel)%kvert(7,i), nnvt + myRF(iel)%kvert(8,i) 
 end do
 nnvt = nnvt + myRF(iel)%nOfVert
end do

WRITE(1,'(A)') 'KNPR'
do iel=1,nel
 do i=1,myRF(iel)%nOfVert
  WRITE(1,'(I8)') 0
 end do
end do

CLOSE(1)

END SUBROUTINE Output_RefTriMesh
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE Output_TriMesh()
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
INTEGER i,j,iBnds
CHARACTER cf*(256)

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/'//adjustl(trim(cProjectGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
WRITE(1,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(1,'(2I8,A)') mg_mesh%level(ilev)%NEL,mg_mesh%level(ilev)%NVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(1,'(A)') 'DCORVG'
DO i = 1,mg_mesh%level(ilev)%nvt
 WRITE(1,'(3ES13.5)') mg_mesh%level(ilev)%dcorvg(:,i)
END DO

WRITE(1,'(A)') 'KVERT'
DO i = 1,mg_mesh%level(ilev)%nel
 WRITE(1,'(8I8)') mg_mesh%level(ILEV)%kvert(:,i)
END DO

WRITE(1,'(A)') 'KNPR'
IF (allocated(myBoundary%bWall)) then
 DO i = 1,mg_mesh%level(ilev)%nvt
  IF (myBoundary%bWall(i)) THEN
   WRITE(1,'(I8)') 1
  ELSE
   WRITE(1,'(I8)') 0
  END IF
 END DO
ELSE
 DO i = 1,mg_mesh%level(ilev)%nvt
  WRITE(1,'(I8)') 0
 END DO
END IF

CLOSE(1)

OPEN(UNIT=2,FILE=ADJUSTL(TRIM(cOutputFolder))//'/'//adjustl(trim(cShortProjectFile)))
!WRITE(cf,'(A11,I3.3,A4)') 'cMESH_',iO, '.tri'
WRITE(2,'(A)') adjustl(trim(cProjectGridFile))
 
DO iBnds = 1, nBnds
 cf = ' '
 WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//"/"//ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
 WRITE(2,'(A)') ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
 WRITE(*,*) "Outputting actual parametrization into: '"//ADJUSTL(TRIM(cf))//"'"
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
 j=0
 DO i=1,mg_mesh%level(ilev)%nvt
  IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
   j = j + 1
  END IF
 END DO
 WRITE(1,'(I8,A)') j," "//myParBndr(iBnds)%Types
 WRITE(1,'(A)')    "'"//ADJUSTL(TRIM(myParBndr(iBnds)%Parameters))//"'"
 DO i=1,mg_mesh%level(ilev)%nvt
  IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
   WRITE(1,'(I8,A)') i
  END IF
 END DO
 CLOSE(1)
END DO
CLOSE(2)

END SUBROUTINE Output_TriMesh
! ----------------------------------------------
SUBROUTINE FillUpRefined(RF,E,nCBP)
integer nCBP
type(RefinerMesh) RF
real*8 E(3,8)
real*8 :: Q8=0.125d0
real*8 DJ(8,3),XX,YY,ZZ
real*8 :: xi1,xi2,xi3,DJAC(3,3),DETJ
integer i,j,k,ii

DJ(1,1)=( E(1,1)+E(1,2)+E(1,3)+E(1,4)+E(1,5)+E(1,6)+E(1,7)+E(1,8))*Q8
DJ(1,2)=( E(2,1)+E(2,2)+E(2,3)+E(2,4)+E(2,5)+E(2,6)+E(2,7)+E(2,8))*Q8
DJ(1,3)=( E(3,1)+E(3,2)+E(3,3)+E(3,4)+E(3,5)+E(3,6)+E(3,7)+E(3,8))*Q8
DJ(2,1)=(-E(1,1)+E(1,2)+E(1,3)-E(1,4)-E(1,5)+E(1,6)+E(1,7)-E(1,8))*Q8
DJ(2,2)=(-E(2,1)+E(2,2)+E(2,3)-E(2,4)-E(2,5)+E(2,6)+E(2,7)-E(2,8))*Q8
DJ(2,3)=(-E(3,1)+E(3,2)+E(3,3)-E(3,4)-E(3,5)+E(3,6)+E(3,7)-E(3,8))*Q8
DJ(3,1)=(-E(1,1)-E(1,2)+E(1,3)+E(1,4)-E(1,5)-E(1,6)+E(1,7)+E(1,8))*Q8
DJ(3,2)=(-E(2,1)-E(2,2)+E(2,3)+E(2,4)-E(2,5)-E(2,6)+E(2,7)+E(2,8))*Q8
DJ(3,3)=(-E(3,1)-E(3,2)+E(3,3)+E(3,4)-E(3,5)-E(3,6)+E(3,7)+E(3,8))*Q8
DJ(4,1)=(-E(1,1)-E(1,2)-E(1,3)-E(1,4)+E(1,5)+E(1,6)+E(1,7)+E(1,8))*Q8
DJ(4,2)=(-E(2,1)-E(2,2)-E(2,3)-E(2,4)+E(2,5)+E(2,6)+E(2,7)+E(2,8))*Q8
DJ(4,3)=(-E(3,1)-E(3,2)-E(3,3)-E(3,4)+E(3,5)+E(3,6)+E(3,7)+E(3,8))*Q8
DJ(5,1)=( E(1,1)-E(1,2)+E(1,3)-E(1,4)+E(1,5)-E(1,6)+E(1,7)-E(1,8))*Q8
DJ(5,2)=( E(2,1)-E(2,2)+E(2,3)-E(2,4)+E(2,5)-E(2,6)+E(2,7)-E(2,8))*Q8
DJ(5,3)=( E(3,1)-E(3,2)+E(3,3)-E(3,4)+E(3,5)-E(3,6)+E(3,7)-E(3,8))*Q8
DJ(6,1)=( E(1,1)-E(1,2)-E(1,3)+E(1,4)-E(1,5)+E(1,6)+E(1,7)-E(1,8))*Q8
DJ(6,2)=( E(2,1)-E(2,2)-E(2,3)+E(2,4)-E(2,5)+E(2,6)+E(2,7)-E(2,8))*Q8
DJ(6,3)=( E(3,1)-E(3,2)-E(3,3)+E(3,4)-E(3,5)+E(3,6)+E(3,7)-E(3,8))*Q8
DJ(7,1)=( E(1,1)+E(1,2)-E(1,3)-E(1,4)-E(1,5)-E(1,6)+E(1,7)+E(1,8))*Q8
DJ(7,2)=( E(2,1)+E(2,2)-E(2,3)-E(2,4)-E(2,5)-E(2,6)+E(2,7)+E(2,8))*Q8
DJ(7,3)=( E(3,1)+E(3,2)-E(3,3)-E(3,4)-E(3,5)-E(3,6)+E(3,7)+E(3,8))*Q8
DJ(8,1)=(-E(1,1)+E(1,2)-E(1,3)+E(1,4)+E(1,5)-E(1,6)+E(1,7)-E(1,8))*Q8
DJ(8,2)=(-E(2,1)+E(2,2)-E(2,3)+E(2,4)+E(2,5)-E(2,6)+E(2,7)-E(2,8))*Q8
DJ(8,3)=(-E(3,1)+E(3,2)-E(3,3)+E(3,4)+E(3,5)-E(3,6)+E(3,7)-E(3,8))*Q8

ii = 0
do i=1,nCBP
 do j=1,nCBP
  do k=1,nCBP
   
   xi1 = -1d0 + 2d0*DBLE(i-1)/DBLE(nCBP-1)
   xi2 = -1d0 + 2d0*DBLE(j-1)/DBLE(nCBP-1)
   xi3 = -1d0 + 2d0*DBLE(k-1)/DBLE(nCBP-1)
   ii = ii + 1

   DJAC(1,1)=DJ(2,1) + DJ(5,1)*XI2 + DJ(6,1)*XI3 + DJ(8,1)*XI2*XI3
   DJAC(1,2)=DJ(3,1) + DJ(5,1)*XI1 + DJ(7,1)*XI3 + DJ(8,1)*XI1*XI3
   DJAC(1,3)=DJ(4,1) + DJ(6,1)*XI1 + DJ(7,1)*XI2 + DJ(8,1)*XI1*XI2
   DJAC(2,1)=DJ(2,2) + DJ(5,2)*XI2 + DJ(6,2)*XI3 + DJ(8,2)*XI2*XI3
   DJAC(2,2)=DJ(3,2) + DJ(5,2)*XI1 + DJ(7,2)*XI3 + DJ(8,2)*XI1*XI3
   DJAC(2,3)=DJ(4,2) + DJ(6,2)*XI1 + DJ(7,2)*XI2 + DJ(8,2)*XI1*XI2
   DJAC(3,1)=DJ(2,3) + DJ(5,3)*XI2 + DJ(6,3)*XI3 + DJ(8,3)*XI2*XI3
   DJAC(3,2)=DJ(3,3) + DJ(5,3)*XI1 + DJ(7,3)*XI3 + DJ(8,3)*XI1*XI3
   DJAC(3,3)=DJ(4,3) + DJ(6,3)*XI1 + DJ(7,3)*XI2 + DJ(8,3)*XI1*XI2
   
   DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
        -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
        +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
     
   XX=DJ(1,1) + DJAC(1,1)*XI1 + DJ(3,1)*XI2 + DJ(4,1)*XI3 + DJ(7,1)*XI2*XI3
   YY=DJ(1,2) + DJ(2,2)*XI1 + DJAC(2,2)*XI2 + DJ(4,2)*XI3 + DJ(6,2)*XI1*XI3
   ZZ=DJ(1,3) + DJ(2,3)*XI1 + DJ(3,3)*XI2 + DJAC(3,3)*XI3 + DJ(5,3)*XI1*XI2
   
   RF%dcoor(:,ii) = [xx,yy,zz]
     
!    write(*,*) RF%dcoor(:,ii)
   
   end do
 end do
end do

ii = 0
do i=1,nCBP-1
 do j=1,nCBP-1
  do k=1,nCBP-1
   
   ii = ii + 1
   RF%kvert(:,ii) = [(k-1)*nCBP*nCBP+(j-1)*nCBP + i,&
                     (k-1)*nCBP*nCBP+(j-1)*nCBP + i + 1 ,&
                     (k-1)*nCBP*nCBP+(j)*nCBP + i + 1 ,&
                     (k-1)*nCBP*nCBP+(j)*nCBP + i ,&
                     (k)*nCBP*nCBP+(j-1)*nCBP + i ,&
                     (k)*nCBP*nCBP+(j-1)*nCBP + i + 1 ,&
                     (k)*nCBP*nCBP+(j)*nCBP + i + 1 ,&
                     (k)*nCBP*nCBP+(j)*nCBP + i ]
   
!    write(*,*) RF%kvert(:,ii)
   
   end do
 end do
end do

! pause

END SUBROUTINE FillUpRefined
! ----------------------------------------------
SUBROUTINE FillUpRefinedElement(RF,E,iP,bV,bRecursive)
logical bRecursive
integer :: iP
type(RefinerMesh) RF
logical bV(8)
real*8 E(3,8)
real*8 :: Q8=0.125d0
real*8 DJ(8,3),XX,YY,ZZ
real*8 :: CBP(3,64),xi1,xi2,xi3,DJAC(3,3),DETJ
integer i,j,k,ii
real*8 :: dPhi = 0.0d0
character*256 cF
integer, allocatable :: kvert(:,:)

 cF = cTemplates(myTemplate)

open(file=TRIM(ADJUSTL(cF)),unit=1)

read(1,*)
read(1,*) 
read(1,*) RF%nOfElem, RF%nOfVert
read(1,*)

allocate(RF%kvert(8,RF%nOfElem))
allocate(RF%dcoor(3,RF%nOfVert))
allocate(RF%knpr(RF%nOfVert))
RF%knpr = 0

DO i=1,RF%nOfVert
 read(1,*) CBP(:,i)
end do

read(1,*)
 
DO i=1,RF%nOfElem
 read(1,*) RF%kvert(:,i)
end do
 
close(1)

if (bRecursive) then

 if (myTemplate.eq.5) then
   RF%knpr(20) = 1
   RF%knpr(27) = 1
   RF%knpr(33) = 1
   RF%knpr(22) = 1
   RF%knpr(28) = 1
 end if

 if (myTemplate.eq.10) then
   RF%knpr(31) = 1
   RF%knpr(71) = 1
   RF%knpr(47) = 1
   RF%knpr(59) = 1
   RF%knpr(72) = 1
   RF%knpr(56) = 1
   RF%knpr(87) = 1
 end if
 
 if (myTemplate.eq.15) then
   RF%knpr(5)  = 1
   RF%knpr(12) = 1
   RF%knpr(17) = 1
   RF%knpr(19) = 1
   RF%knpr(29) = 1
   RF%knpr(49) = 1
   RF%knpr(54) = 1
   RF%knpr(68) = 1
   RF%knpr(70) = 1
   RF%knpr(65) = 1
   RF%knpr(72) = 1
 end if
 
 if (myTemplate.eq.18) then
  RF%knpr(8 ) = 1
  RF%knpr(10) = 1
  RF%knpr(13) = 1
  RF%knpr(21) = 1
  RF%knpr(25) = 1
  RF%knpr(32) = 1
  RF%knpr(38) = 1
  RF%knpr(42) = 1
  RF%knpr(49) = 1
  RF%knpr(4 ) = 1
  RF%knpr(9 ) = 1
  RF%knpr(19) = 1
  RF%knpr(29) = 1
  RF%knpr(36) = 1
  RF%knpr(46) = 1
 END IF
 
 if (myTemplate.eq.21) then
  RF%knpr(87 ) = 1
  RF%knpr(88 ) = 1
  RF%knpr(89 ) = 1
  RF%knpr(92 ) = 1
  RF%knpr(93 ) = 1
  RF%knpr(97 ) = 1
  RF%knpr(100) = 1
  RF%knpr(103) = 1
  RF%knpr(104) = 1
  RF%knpr(65 ) = 1
  RF%knpr(67 ) = 1
  RF%knpr(70 ) = 1
  RF%knpr(71 ) = 1
  RF%knpr(75 ) = 1
  RF%knpr(77 ) = 1
  RF%knpr(72 ) = 1
  RF%knpr(73 ) = 1
  RF%knpr(75 ) = 1
  RF%knpr(77 ) = 1
  RF%knpr(79 ) = 1
  RF%knpr(85 ) = 1
 END IF
  
end if

DJ(1,1)=( E(1,1)+E(1,2)+E(1,3)+E(1,4)+E(1,5)+E(1,6)+E(1,7)+E(1,8))*Q8
DJ(1,2)=( E(2,1)+E(2,2)+E(2,3)+E(2,4)+E(2,5)+E(2,6)+E(2,7)+E(2,8))*Q8
DJ(1,3)=( E(3,1)+E(3,2)+E(3,3)+E(3,4)+E(3,5)+E(3,6)+E(3,7)+E(3,8))*Q8
DJ(2,1)=(-E(1,1)+E(1,2)+E(1,3)-E(1,4)-E(1,5)+E(1,6)+E(1,7)-E(1,8))*Q8
DJ(2,2)=(-E(2,1)+E(2,2)+E(2,3)-E(2,4)-E(2,5)+E(2,6)+E(2,7)-E(2,8))*Q8
DJ(2,3)=(-E(3,1)+E(3,2)+E(3,3)-E(3,4)-E(3,5)+E(3,6)+E(3,7)-E(3,8))*Q8
DJ(3,1)=(-E(1,1)-E(1,2)+E(1,3)+E(1,4)-E(1,5)-E(1,6)+E(1,7)+E(1,8))*Q8
DJ(3,2)=(-E(2,1)-E(2,2)+E(2,3)+E(2,4)-E(2,5)-E(2,6)+E(2,7)+E(2,8))*Q8
DJ(3,3)=(-E(3,1)-E(3,2)+E(3,3)+E(3,4)-E(3,5)-E(3,6)+E(3,7)+E(3,8))*Q8
DJ(4,1)=(-E(1,1)-E(1,2)-E(1,3)-E(1,4)+E(1,5)+E(1,6)+E(1,7)+E(1,8))*Q8
DJ(4,2)=(-E(2,1)-E(2,2)-E(2,3)-E(2,4)+E(2,5)+E(2,6)+E(2,7)+E(2,8))*Q8
DJ(4,3)=(-E(3,1)-E(3,2)-E(3,3)-E(3,4)+E(3,5)+E(3,6)+E(3,7)+E(3,8))*Q8
DJ(5,1)=( E(1,1)-E(1,2)+E(1,3)-E(1,4)+E(1,5)-E(1,6)+E(1,7)-E(1,8))*Q8
DJ(5,2)=( E(2,1)-E(2,2)+E(2,3)-E(2,4)+E(2,5)-E(2,6)+E(2,7)-E(2,8))*Q8
DJ(5,3)=( E(3,1)-E(3,2)+E(3,3)-E(3,4)+E(3,5)-E(3,6)+E(3,7)-E(3,8))*Q8
DJ(6,1)=( E(1,1)-E(1,2)-E(1,3)+E(1,4)-E(1,5)+E(1,6)+E(1,7)-E(1,8))*Q8
DJ(6,2)=( E(2,1)-E(2,2)-E(2,3)+E(2,4)-E(2,5)+E(2,6)+E(2,7)-E(2,8))*Q8
DJ(6,3)=( E(3,1)-E(3,2)-E(3,3)+E(3,4)-E(3,5)+E(3,6)+E(3,7)-E(3,8))*Q8
DJ(7,1)=( E(1,1)+E(1,2)-E(1,3)-E(1,4)-E(1,5)-E(1,6)+E(1,7)+E(1,8))*Q8
DJ(7,2)=( E(2,1)+E(2,2)-E(2,3)-E(2,4)-E(2,5)-E(2,6)+E(2,7)+E(2,8))*Q8
DJ(7,3)=( E(3,1)+E(3,2)-E(3,3)-E(3,4)-E(3,5)-E(3,6)+E(3,7)+E(3,8))*Q8
DJ(8,1)=(-E(1,1)+E(1,2)-E(1,3)+E(1,4)+E(1,5)-E(1,6)+E(1,7)-E(1,8))*Q8
DJ(8,2)=(-E(2,1)+E(2,2)-E(2,3)+E(2,4)+E(2,5)-E(2,6)+E(2,7)-E(2,8))*Q8
DJ(8,3)=(-E(3,1)+E(3,2)-E(3,3)+E(3,4)+E(3,5)-E(3,6)+E(3,7)-E(3,8))*Q8

DO ii=1,RF%nOfVert

   xi1 = CBP(1,ii)
   xi2 = CBP(2,ii)
   xi3 = CBP(3,ii)

   DJAC(1,1)=DJ(2,1) + DJ(5,1)*XI2 + DJ(6,1)*XI3 + DJ(8,1)*XI2*XI3
   DJAC(1,2)=DJ(3,1) + DJ(5,1)*XI1 + DJ(7,1)*XI3 + DJ(8,1)*XI1*XI3
   DJAC(1,3)=DJ(4,1) + DJ(6,1)*XI1 + DJ(7,1)*XI2 + DJ(8,1)*XI1*XI2
   DJAC(2,1)=DJ(2,2) + DJ(5,2)*XI2 + DJ(6,2)*XI3 + DJ(8,2)*XI2*XI3
   DJAC(2,2)=DJ(3,2) + DJ(5,2)*XI1 + DJ(7,2)*XI3 + DJ(8,2)*XI1*XI3
   DJAC(2,3)=DJ(4,2) + DJ(6,2)*XI1 + DJ(7,2)*XI2 + DJ(8,2)*XI1*XI2
   DJAC(3,1)=DJ(2,3) + DJ(5,3)*XI2 + DJ(6,3)*XI3 + DJ(8,3)*XI2*XI3
   DJAC(3,2)=DJ(3,3) + DJ(5,3)*XI1 + DJ(7,3)*XI3 + DJ(8,3)*XI1*XI3
   DJAC(3,3)=DJ(4,3) + DJ(6,3)*XI1 + DJ(7,3)*XI2 + DJ(8,3)*XI1*XI2
   
   DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
        -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
        +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
     
   XX=DJ(1,1) + DJAC(1,1)*XI1 + DJ(3,1)*XI2 + DJ(4,1)*XI3 + DJ(7,1)*XI2*XI3
   YY=DJ(1,2) + DJ(2,2)*XI1 + DJAC(2,2)*XI2 + DJ(4,2)*XI3 + DJ(6,2)*XI1*XI3
   ZZ=DJ(1,3) + DJ(2,3)*XI1 + DJ(3,3)*XI2 + DJAC(3,3)*XI3 + DJ(5,3)*XI1*XI2
   
   RF%dcoor(:,ii) = [xx,yy,zz]
   
END DO


END SUBROUTINE FillUpRefinedElement
! ----------------------------------------------
SUBROUTINE RotatePatch(C,ID)
integer ID(8),IDR(8),IDR1(8),IDR2(8),IDR3(8)
logical :: C(8),R1(8),R2(8),R3(8),P(8)
integer i,iSum,jSum,j,jSum_init

jSum_init = 64*8
jSum = jSum_init
 
 DO i=1,8
  IF (C(i)) THEN
   
   R1 = C(RROT(:,i))
   IDR1 = ID(RROT(:,i))
   iSum = 0
   do j=1,8
    if (R1(j)) iSum = iSum + j*j
   end do
   if (iSum.lt.jSum) then
    P = R1
    IDR = IDR1
    jSum = iSum
   end if
   
   R2 = R1(R2ROT)
   IDR2 = IDR1(R2ROT)
   iSum = 0
   do j=1,8
    if (R2(j)) iSum = iSum + j*j
   end do
   if (iSum.lt.jSum) then
    P = R2
    IDR = IDR2
    jSum = iSum
   end if
   
   R3 = R1(R3ROT)
   IDR3 = IDR1(R3ROT)
   iSum = 0
   do j=1,8
    if (R3(j)) iSum = iSum + j*j
   end do
   if (iSum.lt.jSum) then
    P = R3
    IDR = IDR3
    jSum = iSum
   end if
   
  END IF
 END DO

 if (jSum.lt.jSum_init) then
  ID = IDR
  C=P
 end if
!   WRITE(*,'(8L,A,8L,A,8L,A,8L)') C , " :: ", R1,  " | ", R2,  " | ", R3

END SUBROUTINE RotatePatch
!
!------------------------------------------------------------------
!
SUBROUTINE DetermineTemplate(C,iT)
logical C(8),R(8)
integer iT
integer i
logical bSame

do i=1,22
 R = templates(:,i)
 bSame=.true.
 call ComparePatch(C,R,bSame)
!  write(*,*) C,":",R,":",bSame
 if (.not.bSame) then
  iT = i
  RETURN
 end if
end do

END SUBROUTINE DetermineTemplate
!
!------------------------------------------------------------------
!
SUBROUTINE ComparePatch(C,R,b)
logical :: C(8),R(8)
logical b
integer i,j

j = 0
DO i=1,8
 if ((C(i).and.R(i)).or.((.not.C(i)).and.(.not.R(i)))) j = j + 1
END DO
 
if (j.eq.8) b=.false.
!  write(*,*) C,":",R,":",j
 
END SUBROUTINE ComparePatch
!
!------------------------------------------------------------------
!
SUBROUTINE CleanUpPatches()

integer iel,jel,i,j,jj,iat,nUniquePoints,nAllPoints,nAllElements,il
integer nTotalPoints,nTotalElements
real*8 Pi(3),Pj(3),dist
logical bFound

integer I8(8)
real*8 P8(3,8),dVol,dTol

ilev = mg_Mesh%nlmax
nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nTotalPoints = 0
nTotalElements = 0
nAllPoints = 0
nAllElements = 0

do iel=1,nel

 nAllPoints = 0
 nAllElements = 0

 if (allocated(myRF(iel)%myRF)) then
  CALL CountSubPatch(myRF(iel))
 else
  nAllPoints = myRF(iel)%nOfVert
  nAllElements = myRF(iel)%nOfElem
 end if
 
 allocate(myRF(iel)%dUniquedCoor(3,nAllPoints))
 allocate(myRF(iel)%kUniqueElem(8,nAllElements))
 
 nUniquePoints = 0
 nAllElements  = 0
 if (allocated(myRF(iel)%myRF)) then
  il = 0
  CALL CleanUpSubPatch(myRF(iel))
 else
  nUniquePoints = myRF(iel)%nOfVert
  nAllElements = myRF(iel)%nOfElem
  myRF(iel)%dUniquedCoor(1:3,1:nUniquePoints) = myRF(iel)%dcoor(1:3,1:nUniquePoints)
  myRF(iel)%kUniqueElem(1:8,1:nAllElements) = myRF(iel)%kvert(1:8,1:nAllElements)
  
  I8 = mg_mesh%level(ilev)%kvert(:,iel)
  do jj=1,8
   P8(:,jj) = mg_mesh%level(ilev)%dcorvg(:,I8(jj))
  end do
  CALL getVol(P8,dVol)
  dTol = 0.01d0*dVol**(0.333333d0)
  myRF(iel)%dTol = min(myRF(iel)%dTol,dTol)
 end if

 myRF(iel)%nUniquePoints = nUniquePoints
 myRF(iel)%nUniqueElems  = nAllElements
 nTotalPoints   = nTotalPoints + nUniquePoints
 nTotalElements = nTotalElements + nAllElements
!  write(*,*) myRF(iel)%dTol

end do

write(*,*) "nTotalElements,nTotalPoints",nTotalElements,nTotalPoints

! nTotalPoints = 0
! nTotalElements = 0
! do iel=1,nel
!  nTotalPoints   = nTotalPoints + myRF(iel)%nUniquePoints
!  nTotalElements = nTotalElements + myRF(iel)%nUniqueElems
! enddo
! write(*,*) "nTotalElements,nTotalPoints",nTotalElements,nTotalPoints

 CONTAINS
 
 RECURSIVE SUBROUTINE CountSubPatch(M)
 type(RefinerMesh) :: M
 integer jj

 if (allocated(M%myRF)) then
  DO jj=1,M%nOfElem
   call CountSubPatch(M%myRF(jj))
  END DO
 else
  nAllPoints = nAllPoints + M%nOfVert
  nAllElements = nAllElements + M%nOfElem
 end if
 
 END SUBROUTINE CountSubPatch
! 
!---------------------------------------------------------------
!
 RECURSIVE SUBROUTINE CleanUpSubPatch(M)
 type(RefinerMesh) :: M
 integer ii,jj,kk,I8(8)
 real*8 Pi(3),Pj(3),dist,P8(3,8),dVol
 real*8 :: dTol=1d-4
 integer iFound
 integer, allocatable :: iBuff(:)
 
 allocate(iBuff(M%nOfVert))
 
 if (allocated(M%myRF)) then
  il = il + 1
  DO jj=1,M%nOfElem
   call CleanUpSubPatch(M%myRF(jj))
  END DO
 else
 
  I8 = mg_mesh%level(ilev)%kvert(:,iel)
  do jj=1,8
   P8(:,jj) = mg_mesh%level(ilev)%dcorvg(:,I8(jj))
  end do
  CALL getVol(P8,dVol)
  dTol = 0.01d0*dVol**(0.333333d0)
  myRF(iel)%dTol = min(myRF(iel)%dTol,dTol)
  
  DO jj = 1,M%nOfVert
   Pi = M%dcoor(:,jj)
   
   iFound = 0
   DO kk = 1, nUniquePoints
    Pj = myRF(iel)%dUniquedCoor(:,kk)
    dist =  sqrt((Pi(1)-Pj(1))**2d0 + (Pi(2)-Pj(2))**2d0 + (Pi(3)-Pj(3))**2d0)
    if (dist.lt.dTol) then
     iFound = kk
     exit
    END IF
   END DO
   
   if (iFound.gt.0) then
     iBuff(jj) = iFound
   else
     nUniquePoints = nUniquePoints + 1
     myRF(iel)%dUniquedCoor(:,nUniquePoints) = Pi
     iBuff(jj) = nUniquePoints
   end if
  END DO
  
  DO jj=1,M%nOfElem
   do kk=1,8
    myRF(iel)%kUniqueElem(kk,nAllElements+jj) =  iBuff(M%kvert(kk,jj))
   end do
  END DO
  
  nAllElements = nAllElements + M%nOfElem
  
 end if
 
 deallocate(iBuff)
 
!   WRITE(*,*) "<zsdsaddasdsa sadsa sa "
!   do jj=1,nUniquePoints
!     write(*,'(3ES12.3)') myRF(iel)%dUniquedCoor(:,jj)
!   end do
!   pause
!   do jj=1,nAllElements
!     write(*,'(8I8)') myRF(iel)%kUniqueElem(:,jj)
!   end do

 
 END SUBROUTINE CleanUpSubPatch
 
END SUBROUTINE CleanUpPatches
!
!------------------------------------------------------------------
!
SUBROUTINE CleanUpMesh()

logical, allocatable :: bDone(:)
integer iel,jel,i,j,iat,ivt,ivel
real*8 Pi(3),Pj(3),dist
integer iFound

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

allocate(bDone(nel))

bDone = .false.
nUniquePoints = 0
nUniqueElems = 0

do iel=1,nel

 do i=1,myRF(iel)%nUniquePoints
 
  Pi = myRF(iel)%dUniquedCoor(:,i)
  
  iFound = 0
  
  do iat=1,8
   ivt = mg_mesh%level(ilev)%kvert(iat,iel)
   do ivel=1,mg_mesh%level(ilev)%nvel
    jel =  mg_mesh%level(ilev)%kvel(ivel,ivt)
    if (jel.gt.0) then
     if (bDone(jel)) then
      do j=1,myRF(jel)%nUniquePoints
       Pj = myRF(jel)%dUniquedCoor(:,j)
       dist =  sqrt((Pi(1)-Pj(1))**2d0 + (Pi(2)-Pj(2))**2d0 + (Pi(3)-Pj(3))**2d0)
       if (dist.lt.myRF(jel)%dTol) then
        iFound = 1
        GOTO 1
       end if
      end do
     end if
    end if
   end do
  end do
1 continue
  if (iFound.ge.1) then
  else
   nUniquePoints = nUniquePoints + 1 
  end if
 end do
 bDone(iel) = .true.
 nUniqueElems = nUniqueElems + myRF(iel)%nUniqueElems
end do

write(*,*) "nUniquePoints=",nUniquePoints,"nUniqueElems=",nUniqueElems

allocate(MergedMeshCoor(3,nUniquePoints))
allocate(MergedMeshElem(8,nUniqueElems))

bDone = .false.
nUniquePoints = 0
nUniqueElems = 0

do iel=1,nel

 allocate(myRF(iel)%PointerToMerged(myRF(iel)%nUniquePoints))

 do i=1,myRF(iel)%nUniquePoints
 
  Pi = myRF(iel)%dUniquedCoor(:,i)
  
  iFound = 0
  
  do iat=1,8
   
   ivt = mg_mesh%level(ilev)%kvert(iat,iel)
   
   do ivel=1,mg_mesh%level(ilev)%nvel
   
      jel =  mg_mesh%level(ilev)%kvel(ivel,ivt)
      if (jel.gt.0) then
       if (bDone(jel)) then
        do j=1,myRF(jel)%nUniquePoints
         Pj = myRF(jel)%dUniquedCoor(:,j)
         dist =  sqrt((Pi(1)-Pj(1))**2d0 + (Pi(2)-Pj(2))**2d0 + (Pi(3)-Pj(3))**2d0)
         if (dist.lt.myRF(jel)%dTol) then
          iFound = myRF(jel)%PointerToMerged(j)
          if (iFound.lt.1) then
           write(*,*) 'o uuu ...'
          end if
          GOTO 2
         end if
        end do
       end if
      end if
   
   end do
  end do

2 continue

  if (iFound.ge.1) then
!    write(*,*) "iFound:", iFound
   myRF(iel)%PointerToMerged(i) = iFound
  else
   nUniquePoints = nUniquePoints + 1 
   MergedMeshCoor(:,nUniquePoints) = Pi
   myRF(iel)%PointerToMerged(i) = nUniquePoints
  end if
 end do
 
 do i=1,myRF(iel)%nUniqueElems
  do j=1,8
   ivt = myRF(iel)%kUniqueElem(j,i)
   MergedMeshElem(j,nUniqueElems+i) = myRF(iel)%PointerToMerged(ivt)
  end do
 end do
 
 bDone(iel) = .true.
 nUniqueElems = nUniqueElems + myRF(iel)%nUniqueElems
 
end do


END SUBROUTINE CleanUpMesh
!
!------------------------------------------------------------------
!
SUBROUTINE getVol(P,dVol)

REAL*8 P(3,8),dVol
REAL*8 :: A1 = 1d0/6d0

 dVol=A1*((DABS((P(1,4)-P(1,1))*(P(2,4)-P(2,3))*(P(3,4)-P(3,8))+(P(2,4)-P(2,1))*  &
      (P(3,4)-P(3,3))*(P(1,4)-P(1,8))+(P(3,4)-P(3,1))*(P(1,4)-P(1,3))*(P(2,4)-P(2,8))- &
      (P(1,4)-P(1,8))*(P(2,4)-P(2,3))*(P(3,4)-P(3,1))-(P(2,4)-P(2,8))*(P(3,4)-P(3,3))* &
      (P(1,4)-P(1,1))-(P(3,4)-P(3,8))*(P(1,4)-P(1,3))*(P(2,4)-P(2,1))))+       &
      (DABS((P(1,2)-P(1,3))*(P(2,2)-P(2,1))*(P(3,2)-P(3,6))+(P(2,2)-P(2,3))*   &
      (P(3,2)-P(3,1))*(P(1,2)-P(1,6))+(P(3,2)-P(3,3))*(P(1,2)-P(1,1))*(P(2,2)-P(2,6))- &
      (P(1,2)-P(1,6))*(P(2,2)-P(2,1))*(P(3,2)-P(3,3))-(P(2,2)-P(2,6))*(P(3,2)-P(3,1))* &
      (P(1,2)-P(1,3))-(P(3,2)-P(3,6))*(P(1,2)-P(1,1))*(P(2,2)-P(2,3))))+       &
      (DABS((P(1,5)-P(1,8))*(P(2,5)-P(2,6))*(P(3,5)-P(3,1))+(P(2,5)-P(2,8))*   &
      (P(3,5)-P(3,6))*(P(1,5)-P(1,1))+(P(3,5)-P(3,8))*(P(1,5)-P(1,6))*(P(2,5)-P(2,1))- &
      (P(1,5)-P(1,1))*(P(2,5)-P(2,6))*(P(3,5)-P(3,8))-(P(2,5)-P(2,1))*(P(3,5)-P(3,6))* &
      (P(1,5)-P(1,8))-(P(3,5)-P(3,1))*(P(1,5)-P(1,6))*(P(2,5)-P(2,8))))+       &
      (DABS((P(1,7)-P(1,6))*(P(2,7)-P(2,8))*(P(3,7)-P(3,3))+(P(2,7)-P(2,6))*   &
      (P(3,7)-P(3,8))*(P(1,7)-P(1,3))+(P(3,7)-P(3,6))*(P(1,7)-P(1,8))*(P(2,7)-P(2,3))- &
      (P(1,7)-P(1,3))*(P(2,7)-P(2,8))*(P(3,7)-P(3,6))-(P(2,7)-P(2,3))*(P(3,7)-P(3,8))* &
      (P(1,7)-P(1,6))-(P(3,7)-P(3,3))*(P(1,7)-P(1,8))*(P(2,7)-P(2,6))))+       &
      (DABS((P(1,1)-P(1,3))*(P(2,1)-P(2,8))*(P(3,1)-P(3,6))+(P(2,1)-P(2,3))*   &
      (P(3,1)-P(3,8))*(P(1,1)-P(1,6))+(P(3,1)-P(3,3))*(P(1,1)-P(1,8))*(P(2,1)-P(2,6))- &
      (P(1,1)-P(1,6))*(P(2,1)-P(2,8))*(P(3,1)-P(3,3))-(P(2,1)-P(2,6))*(P(3,1)-P(3,8))* &
      (P(1,1)-P(1,3))-(P(3,1)-P(3,6))*(P(1,1)-P(1,8))*(P(2,1)-P(2,3)))))
   dVol = abs(dVol)

END SUBROUTINE getVol

END MODULE MeshRefDef
