PROGRAM MeshRef
USE MeshRefDef
USE MeshRefRefine
USE MeshRefOutput
use f90getopt
USE Parametrization, ONLY: InitParametrization,ParametrizeBndr,&
    ParametrizeBndryPoints,InitBoundaryStructure

use cinterface

USE MESH_Structures

IMPLICIT NONE

character*(INIP_STRLEN) :: sfilename,sfilepath
integer :: unitProtfile = -1, unitTerminal = 6
integer :: MeshingScheme
integer :: sGENDIE=1,sPARTICLE=2,sDEFAULT=0

integer i,ivt
real*8 P8(3,8)
external e013

! character(len=*), parameter :: Folder = 'SimFolder'
type(option_s)              :: opts(3)

! MeshingScheme = sPARTICLE
! MeshingScheme = sDEFAULT
 MeshingScheme = sGENDIE

opts(1) = option_s('inputfolder', .true.,  'f')
opts(2) = option_s('outputfolder',  .false., 'o')
opts(3) = option_s('help',  .false., 'h')

! check the command line arguments
do
    select case (getopt('f:o:h', opts))
        case (char(0))
            exit
        case ('f')
            read(optarg,*) cIntputFolder
            write(*,*)'Input Folder: ', ADJUSTL(TRIM(cIntputFolder))
        case ('-f')
            read(optarg,*) cIntputFolder
            write(*,*)'Input Folder: ', ADJUSTL(TRIM(cIntputFolder))
        case ('o')
            read(optarg,*) cReducedMeshdFolder
            write(*,*)'Reduced Output Folder: ', ADJUSTL(TRIM(cReducedMeshdFolder))
        case ('-o')
            read(optarg,*) cReducedMeshdFolder
            write(*,*)'Reduced Output Folder: ', ADJUSTL(TRIM(cReducedMeshdFolder))
        case ('h')
          call print_help()
          call exit(0)
        case default
    end select
end do

MASTER = 0
bParallel = .false.
bBoundaryCheck = .true.

myid = 1
call init_fc_rigid_body(myid)      
myid = 0

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! READING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL GetParameters()
NLMAX = mg_Mesh%nlmax 
NLMIN = mg_Mesh%nlmin 
mg_Mesh%maxlevel = mg_Mesh%nlmax+1
allocate(mg_mesh%level(mg_Mesh%maxlevel))

write(*,*) 'Project Grid File: = "'//adjustl(trim(cProjectFolder))//'/'//adjustl(trim(cProjectGridFile))
call readTriCoarse(adjustl(trim(cProjectFolder))//'/'//adjustl(trim(cProjectGridFile)), mg_mesh)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! READING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! building up the mesh structures
call refineMesh(mg_mesh, mg_Mesh%maxlevel)  

CALL InitMarker()
call SetUpMarker()

! call CreateRefinedMesh()
if (MeshingScheme.eq.sDEFAULT) then
 call CreateRefinedMesh_Fine()
end if
if (MeshingScheme.eq.sGENDIE) then
 call CreateRefinedMesh_Fine()
end if
if (MeshingScheme.eq.sPARTICLE) then
 call CreateRefinedParticleMesh()
end if

call CleanUpPatches()

call CleanUpMesh()

call CutMesh()

call inip_output_init(1,1,unitProtfile,unitTerminal)
sfilepath=adjustl(trim(cOutputFolder))//'/meshDir_BU'
if (.not. inip_isDirectory(sfilepath)) then
  write(*,*) 'Creating Directory: '//adjustl(trim(sfilepath))
  call inip_makeDirectory(sfilepath)
end if
sfilepath=adjustl(trim(cOutputFolder))//'/'//adjustl(trim(cReducedMeshdFolder))
if (.not. inip_isDirectory(sfilepath)) then
  write(*,*) 'Creating Directory: '//adjustl(trim(sfilepath))
  call inip_makeDirectory(sfilepath)
end if
! sfilepath=adjustl(trim(cOutputFolder))//'/meshDir'
! call inip_makeDirectory(sfilepath)
ilev = lTriOutputLevel
CALL Output_TriMesh()
CALL Output_RefTriMesh()

if (MeshingScheme.eq.sGENDIE) then
 CALL Output_MergedRefTriMesh()
 CALL Output_MergedRefTriMeshPar()
end if
if (MeshingScheme.eq.sPARTICLE) then
 CALL Output_ParticleMergedRefTriMesh()
 CALL Output_ParticleMergedRefTriMeshPar()
end if
if (MeshingScheme.eq.sDEFAULT) then
 CALL Output_MergedRefTriMesh()
END IF

ilev = lVTUOutputLevel
CALL Output_VTK()
CALL Output_RefVTK()
CALL Output_UniqueRefVTK()
CALL Output_MergedRefVTK()

if (MeshingScheme.eq.sGENDIE) then
 CALL CreateReducedMesh()
 CALL Output_ReducedRefVTK()
 CALL Output_ReducedRefTRI()
 mg_ReducedMesh%nlmax = 1
 mg_ReducedMesh%nlmin = 1
 NLMAX = mg_ReducedMesh%nlmax 
 NLMIN = mg_ReducedMesh%nlmin 
 mg_ReducedMesh%maxlevel = mg_ReducedMesh%nlmax+1
 allocate(mg_ReducedMesh%level(mg_ReducedMesh%maxlevel))
 call readTriCoarse(adjustl(trim(cIntputFolder))//'/'//adjustl(trim(cReducedGridFile)), mg_ReducedMesh)
 call refineMesh(mg_ReducedMesh, mg_ReducedMesh%maxlevel)  
 CALL CreateCleanReducedMesh()
 CALL Output_ReducedCleanRefVTK()
 CALL Output_ReducedRefTRIandPAR()
end if

if (MeshingScheme.eq.sPARTICLE) then
 mg_ParticleMesh%nlmax = 2
 mg_ParticleMesh%nlmin = 1
 NLMAX = mg_ParticleMesh%nlmax 
 NLMIN = mg_ParticleMesh%nlmin 
 mg_ParticleMesh%maxlevel = mg_ParticleMesh%nlmax+1
 allocate(mg_ParticleMesh%level(mg_ParticleMesh%maxlevel))
 call readTriCoarse(adjustl(trim(cIntputFolder))//'/meshDir/Merged_Mesh.tri', mg_ParticleMesh)
 call refineMesh(mg_ParticleMesh, mg_ParticleMesh%maxlevel)
 cProjectFolder =  'myMesh/meshDir/'
 cProjectFile='myMesh/meshDir/file.prj'

 ILEV=NLMIN
 CALL InitParametrization(mg_ParticleMesh%level(ILEV),ILEV)
 
 DO ILEV=NLMIN,NLMAX

   CALL ParametrizeBndr(mg_ParticleMesh,ilev)

   CALL ProlongateCoordinates(mg_ParticleMesh%level(ILEV)%dcorvg,&
                              mg_ParticleMesh%level(ILEV+1)%dcorvg,&
                              mg_ParticleMesh%level(ILEV)%karea,&
                              mg_ParticleMesh%level(ILEV)%kvert,&
                              mg_ParticleMesh%level(ILEV)%kedge,&
                              mg_ParticleMesh%level(ILEV)%nel,&
                              mg_ParticleMesh%level(ILEV)%nvt,&
                              mg_ParticleMesh%level(ILEV)%net,&
                              mg_ParticleMesh%level(ILEV)%nat)
                                
 END DO
 
 ilev = nlmax+1
 CALL ParametrizeBndr(mg_ParticleMesh,ilev)
 
 DO i=1,8
  ivt = mg_ParticleMesh%level(ILEV)%kvert(i,1)
  P8(:,i) = mg_ParticleMesh%level(ILEV)%dcorvg(:,ivt)
 END DO
 P8(1,:) = [+1.0,-1.0,-1.0,+1.0,+1.0,-1.0,-1.0,+1.0]
 P8(2,:) = [+1.0,+1.0,-1.0,-1.0,+1.0,+1.0,-1.0,-1.0]
 P8(3,:) = [+1.0,+1.0,+1.0,+1.0,-1.0,-1.0,-1.0,-1.0]

 CALL DIFFQ2_Elem(P8,E013)

 ilev = nlmax+1
 CALL Output_ParticleFineVTK()
 
end if

 CONTAINS

 subroutine print_help()
     print '(a, /)', 'command-line options:'
     print '(a)',    '  -f      Input Folder'
     print '(a)',    '  -o      Output Folder'
     print '(a)',    '  -h      Print usage information and exit'
 end subroutine print_help  

END PROGRAM MeshRef
