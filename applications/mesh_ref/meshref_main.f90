PROGRAM MeshRef
USE MeshRefDef
USE MeshRefRefine
USE MeshRefOutput

use cinterface

USE MESH_Structures

IMPLICIT NONE

character*(INIP_STRLEN) :: sfilename,sfilepath
integer :: unitProtfile = -1, unitTerminal = 6

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

END PROGRAM MeshRef
