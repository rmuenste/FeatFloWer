PROGRAM AutoParam
USE MeshProcDef
USE MeshProcPDE
USE MESH_Structures

!USE Parametrization, ONLY : InitParametrization
IMPLICIT NONE
! LOGICAL :: bA_MD=.true.
! LOGICAL :: bPDE_MD=.false.

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

write(*,*) 'Project Grid File: = "'//ADJUSTL(TRIM(cProjectGridFile))//'"'
call readTriCoarse(ADJUSTL(TRIM(cProjectGridFile)), mg_mesh)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! READING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! building up the mesh structures
call refineMesh(mg_mesh, mg_Mesh%maxlevel)  

!!! building up the parametrization structures
ilev = 1
CALL InitParametrization_STRCT(mg_mesh%level(ilev),ilev)

DO ILEV=mg_Mesh%nlmin,mg_Mesh%nlmax
  CALL ProlongateParametrization_STRCT(mg_mesh,ilev)
END DO
CALL ProlongateParametrization_STRCT(mg_mesh,mg_Mesh%nlmax+1)
 
myid = 1
ilev = mg_Mesh%nlmax
CALL DeterminePointParametrization_STRCT(mg_mesh,ilev)
 
IF (bPDE_MD) then
 CALL InitMeshDef(bDefTensor)
END IF

IF (bA_MD) then
 !!! performing parametrization
 ilev = mg_Mesh%maxlevel
 CALL ParametrizeBndryPoints_STRCT(mg_mesh,ilev)

 !!! smoothening of the mesh  + parametrization
 CALL SeqUmbrella(mg_Mesh%maxlevel,nUmbrellaSteps)
END IF

CALL InitBoundaryFields()

IF (bPDE_MD) then

 CALL MeshDefPDE(bDefTensor)

 CALL refineMeshLevel(mg_mesh%level(mg_Mesh%nlmax),mg_mesh%level(mg_Mesh%maxlevel))
 
 ilev = mg_Mesh%maxlevel
 CALL ParametrizeBndryPoints_STRCT(mg_mesh,ilev)
 
 CALL Correct_myQ2Coor()
 
END IF

CALL RecoverSurfaceNormals()

ilev = lTriOutputLevel
CALL Output_Mesh(1,cOutputFolder)

ilev = lVTUOutputLevel
CALL Output_VTK()

END PROGRAM AutoParam
