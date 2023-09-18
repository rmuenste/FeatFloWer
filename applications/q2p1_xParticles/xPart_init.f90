!========================================================================================
!                           Sub: init_q2p1_xParicles
!========================================================================================
subroutine init_q2p1_xParicles(log_unit)
    
  USE def_FEAT
  USE PP3D_MPI, ONLY : myid,master,showid
  USE var_QuadScalar, ONLY : mg_Mesh,QuadSc,myExport,Shell

  integer, intent(in) :: log_unit
  real*8, allocatable, Dimension(:) :: xField,yField,zField
  integer ndof
  Character cField*(128)

  !-------INIT PHASE-------

  call General_init_ext(log_unit)

  
  if (myid.ne.0) THEN
  
   call CreateDumpStructures(1)
  
   ILEV = NLMAX
   ndof = mg_Mesh%level(nlmax)%nvt + mg_Mesh%level(nlmax)%net + &
          mg_Mesh%level(nlmax)%nat + mg_Mesh%level(nlmax)%nel 
          
!    write(*,*) 'ndof:',ndof
   allocate(xField(ndof),yField(ndof),zField(ndof))
   xField = -100
   yField = -100
   zField = -100
   
   cField = 'coordinates'
   call subLoadMPIFieldQ2_NX(cField,0,3,xField,yField,zField)
   mg_Mesh%level(nlmax)%dcorvg(1,:) = xField
   mg_Mesh%level(nlmax)%dcorvg(2,:) = yField
   mg_Mesh%level(nlmax)%dcorvg(3,:) = zField
   deallocate(xField,yField,zField)
   
   cField = 'velocity'
   allocate(QuadSc%ValU(ndof),QuadSc%ValV(ndof),QuadSc%ValW(ndof))
   QuadSc%ValU = -100
   QuadSc%ValV = -100
   QuadSc%ValW = -100
   call subLoadMPIFieldQ2_NX(cField,0,3,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)
   ! lets invert the flowfield
   QuadSc%ValU = -QuadSc%ValU
   QuadSc%ValV = -QuadSc%ValV
   QuadSc%ValW = -QuadSc%ValW
   
   cField = 'shell'
   allocate(Shell(ndof))
   shell = -100
   call subLoadMPIFieldQ2_NX(cField,0,1,Shell)
   
  end if 

end subroutine init_q2p1_xParicles
!
!----------------------------------------------
!
SUBROUTINE OUTPUT_xParticleVTK(iOut)
 USE def_FEAT
 USE var_QuadScalar, ONLY : mg_Mesh,QuadSc,myExport
 USE PP3D_MPI, ONLY : myid,master,showid

 implicit none
 integer iOut

 if (myid.eq.1) THEN
  ILEV = myExport%Level
  myExport%Format = 'VTK'
  allocate(myExport%Fields(2))
  myExport%Fields(1) = 'Velocity'
  myExport%Fields(2) = 'Shell'
  CALL Output_VTK_piece(iOut,mg_mesh%level(ILEV)%dcorvg,mg_mesh%level(ILEV)%kvert)   
 end if

END SUBROUTINE OUTPUT_xParticleVTK
!
!----------------------------------------------
!
SUBROUTINE General_init_ext(MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : mg_mesh 
 USE cinterface 
 USE xPart_def

 IMPLICIT NONE
 INTEGER MFILE

 CALL INIT_MPI()
 
 CALL ReadParameters(cParamFile)
 
 mg_Mesh%nlmax = NLMAX
 mg_Mesh%nlmin = NLMIN
 mg_Mesh%maxlevel = mg_Mesh%nlmax+1
 allocate(mg_mesh%level(mg_Mesh%maxlevel))
 
 call readTriCoarse(cMeshFile, mg_mesh)
 
 call refineMesh(mg_mesh, mg_Mesh%maxlevel,.TRUE.)  
 
!     ----------------------------------------------------------            
 call init_fc_rigid_body(myid)      
 !     ----------------------------------------------------------        

END SUBROUTINE General_init_ext
!========================================================================================
!                           Sub: ReadParameters
!========================================================================================

SUBROUTINE ReadParameters()
use iniparser
USE PP3D_MPI, ONLY : myid,master,showid
USE def_FEAT
USE xPart_def
USE var_QuadScalar, ONLY : myExport

implicit none

type(t_parlist) :: parameterlist
integer :: unitProtfile = -1 ! I guess you use mfile here
integer :: unitTerminal = 6 ! I guess you use mterm here

if (myid.eq.1) WRITE(*,*) "Parameter File:",ADJUSTL(TRIM(cParamFile))

call inip_output_init(myid,showid,unitProtfile,unitTerminal)

! Init the parameterlist
call inip_init(parameterlist)


call inip_readfromfile(parameterlist,ADJUSTL(TRIM(cParamFile)))

call INIP_getvalue_Int(parameterlist,"MAIN","nlmax",NLMAX,1)
if (myid.eq.1) WRITE(*,*) "nlmax = ",NLMAX

call INIP_getvalue_Int(parameterlist,"MAIN","nlmin",NLMIN,1)
if (myid.eq.1) WRITE(*,*) "nlmin = ",NLMIN

call INIP_getvalue_Int(parameterlist,"MAIN","nIter",nIter,1000)
if (myid.eq.1) WRITE(*,*) "nIter = ",nIter

call INIP_getvalue_double(parameterlist,"MAIN","xFactor",xFactor,1d0)
if (myid.eq.1) WRITE(*,*) "xFactor = ",xFactor

call INIP_getvalue_double(parameterlist,"MAIN","d_CorrDist",d_CorrDist,0.25d0)
if (myid.eq.1) WRITE(*,*) "d_CorrDist = ",d_CorrDist

call INIP_getvalue_double(parameterlist,"MAIN","dTimeStep",dTimeStep,1.00d1)
if (myid.eq.1) WRITE(*,*) "dTimeStep = ",dTimeStep

call INIP_getvalue_double(parameterlist,"MAIN","minDist",minDist,-0.01d0)
if (myid.eq.1) WRITE(*,*) "minDist = ",minDist

call INIP_getvalue_double(parameterlist,"MAIN","hSize",myParticleParam%hSize,0.2d0)
if (myid.eq.1) WRITE(*,*) "hSize = ",myParticleParam%hSize

call INIP_getvalue_Int(parameterlist,"MAIN","nTime",nTime,5)
if (myid.eq.1) WRITE(*,*) "nTime = ",nTime

call INIP_getvalue_Int(parameterlist,"MAIN","xChunks",xChunks,5)
if (myid.eq.1) WRITE(*,*) "xChunks = ",xChunks

call INIP_getvalue_Int(parameterlist,"MAIN","ExportLevel",myExport%Level,NLMAX)
if (myid.eq.1) WRITE(*,*) "ExportLevel = ",myExport%Level
myExport%LevelMax = myExport%Level

call INIP_getvalue_string(parameterlist,"MAIN","MeshFile",cMeshFile,"_data/meshDir/ReducedMesh.tri")
if (myid.eq.1) WRITE(*,*) "MeshFile = '"//ADJUSTL(TRIM(cMeshFile))//"'"

! call inip_toupper_replace(cRotation)




! call inip_info(parameterlist)

! Clean up the parameterlist
call inip_done(parameterlist)

!pause

END SUBROUTINE ReadParameters