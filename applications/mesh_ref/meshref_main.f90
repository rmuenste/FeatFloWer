PROGRAM MeshRef
USE MeshRefDef
USE MeshRefRefine
USE MeshRefOutput
use f90getopt

use cinterface

USE MESH_Structures

IMPLICIT NONE

character*(INIP_STRLEN) :: sfilename,sfilepath
integer :: unitProtfile = -1, unitTerminal = 6

! character(len=*), parameter :: Folder = 'SimFolder'
type(option_s)              :: opts(3)

opts(1) = option_s('inputfolder', .true.,  'i')
opts(2) = option_s('projectfolder', .true.,  'p')
opts(3) = option_s('help',  .false., 'h')

! check the command line arguments
do
    select case (getopt('i:p:h', opts))
        case (char(0))
            exit
        case ('i')
            read(optarg,*) cIntputFolder
            write(*,*)'Input Folder: ', ADJUSTL(TRIM(cIntputFolder))
        case ('-i')
            read(optarg,*) cIntputFolder
            write(*,*)'Input Folder: ', ADJUSTL(TRIM(cIntputFolder))
        case ('p')
            read(optarg,*) cProjectFolder
            write(*,*)'Project Folder: ', ADJUSTL(TRIM(cProjectFolder))
        case ('-p')
            read(optarg,*) cProjectFolder
            write(*,*)'Project Folder: ', ADJUSTL(TRIM(cProjectFolder))
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
call CreateRefinedMesh_Fine()

call CleanUpPatches()

call CleanUpMesh()

call CutMesh()

call inip_output_init(1,1,unitProtfile,unitTerminal)
sfilepath=adjustl(trim(cOutputFolder))//'/meshDir'
if (.not. inip_isDirectory(sfilepath)) then
  write(*,*) 'Creating Directory: '//adjustl(trim(sfilepath))
  call inip_makeDirectory(sfilepath)
end if
! sfilepath=adjustl(trim(cOutputFolder))//'/meshDir'
! call inip_makeDirectory(sfilepath)

ilev = lTriOutputLevel
CALL Output_TriMesh()
CALL Output_RefTriMesh()
CALL Output_MergedRefTriMesh()
CALL Output_MergedRefTriMeshPar()

ilev = lVTUOutputLevel
CALL Output_VTK()
CALL Output_RefVTK()
CALL Output_UniqueRefVTK()
CALL Output_MergedRefVTK()

 CONTAINS

 subroutine print_help()
     print '(a, /)', 'command-line options:'
     print '(a)',    '  -i      Input Folder'
     print '(a)',    '  -p      Project Folder'
     print '(a)',    '  -h      Print usage information and exit'
 end subroutine print_help  

END PROGRAM MeshRef
