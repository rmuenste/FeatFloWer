PROGRAM AutoParam
USE MeshProcDef
USE MESH_Structures

!USE Parametrization, ONLY : InitParametrization
IMPLICIT NONE

dCGALtoRealFactor = 1d0
MASTER = 0
bParallel = .false.
bBoundaryCheck = .true.
mg_Mesh%nlmax = 1
mg_Mesh%nlmin = 1
NLMAX = mg_Mesh%nlmax 
NLMIN = mg_Mesh%nlmin 
mg_Mesh%maxlevel = mg_Mesh%nlmax+1
allocate(mg_mesh%level(mg_Mesh%maxlevel))

myid = 1
call init_fc_rigid_body(myid)      
myid = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! READING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL GetParameters()

call readTriCoarse(ADJUSTL(TRIM(cProjectGridFile)), mg_mesh)

call refineMesh(mg_mesh, mg_Mesh%maxlevel)  

ilev = 1
CALL InitParametrization_STRCT(mg_mesh%level(ilev),ilev)

DO ILEV=mg_Mesh%nlmin,mg_Mesh%nlmax
  CALL ProlongateParametrization_STRCT(mg_mesh,ilev)
END DO
CALL ProlongateParametrization_STRCT(mg_mesh,mg_Mesh%nlmax+1)
 
myid = 1
ilev = mg_Mesh%nlmax
CALL DeterminePointParametrization_STRCT(mg_mesh,ilev)
 
ilev = mg_Mesh%maxlevel
! CALL ParametrizeBndryPoints_STRCT(mg_mesh,ilev)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! READING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ilev = mg_Mesh%maxlevel
CALL Output_VTK()

END PROGRAM AutoParam
