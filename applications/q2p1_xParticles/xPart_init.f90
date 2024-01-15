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
   ndof = mg_Mesh%level(ilev)%nvt + mg_Mesh%level(ilev)%net + &
          mg_Mesh%level(ilev)%nat + mg_Mesh%level(ilev)%nel 
          
!    write(*,*) 'ndof:',ndof
   allocate(xField(ndof),yField(ndof),zField(ndof))
   xField = -100
   yField = -100
   zField = -100
   
   cField = 'coordinates'
   call subLoadMPIFieldQ2_NX(cField,0,3,xField,yField,zField)
   mg_Mesh%level(ilev)%dcorvg(1,:) = xField
   mg_Mesh%level(ilev)%dcorvg(2,:) = yField
   mg_Mesh%level(ilev)%dcorvg(3,:) = zField
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
!========================================================================================
!                           Sub: init_q2p1_xParicles
!========================================================================================
subroutine OUTPUT_xParticleDUMP(iOut)
 USE def_FEAT
 USE var_QuadScalar, ONLY : mg_Mesh,QuadSc,myExport,GenLinScalar,MaterialDistribution,Shell
 USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI,mpi_comm_subs,subnodes
 USE Sigma_User, ONLY: myMultiMat,myProcess,mySigma

 implicit none
 include 'mpif.h'
 integer, intent(in) :: iOut
 
 integer iFld
 real*8, allocatable, Dimension(:) :: xField,yField,zField
 integer ndof,ierr,i,NumOfTooFarElements,NumOfTooFarElementsOld
 Character cField*(128)
 integer, allocatable :: E(:),F(:)

if (myid.ne.0) THEN
  ILEV = NLMAX
  ndof = mg_Mesh%level(ilev)%nvt + mg_Mesh%level(ilev)%net + &
         mg_Mesh%level(ilev)%nat + mg_Mesh%level(ilev)%nel 
end if

if (myid.ne.0) THEN
  ILEV = NLMAX
  ndof = mg_Mesh%level(ilev)%nvt + mg_Mesh%level(ilev)%net + &
         mg_Mesh%level(ilev)%nat + mg_Mesh%level(ilev)%nel 
          
  allocate(E(mg_mesh%level(ILEV)%nel))
  allocate(F(mg_mesh%level(ILEV)%nel))
  NumOfTooFarElementsOld = 100000000
  
1 CONTINUE

  E = MaterialDistribution(ILEV)%x
  MaterialDistribution(ILEV)%x = 0
  
  CALL ExtractMatrialProperties(E,MaterialDistribution(ILEV)%x,Shell,&
                                mg_mesh%level(ILEV)%dcorvg,&
                                mg_mesh%level(ILEV)%kvert,&
                                mg_mesh%level(ILEV)%kedge,&
                                mg_mesh%level(ILEV)%karea,&
                                mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                                mg_mesh%level(ILEV)%elementsAtVertex,&
                                mg_mesh%level(ILEV)%nvt,&
                                mg_mesh%level(ILEV)%net,&
                                mg_mesh%level(ILEV)%nat,&
                                mg_mesh%level(ILEV)%nel,&
                                NumOfTooFarElements)
  
!   write(*,*) "NumOfTooFarElements: ",NumOfTooFarElements
  call MPI_AllReduce(NumOfTooFarElements,iFld, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_SUBS, ierr)
  NumOfTooFarElements = iFld

  if (myid.eq.1) THEN
   write(*,*) "Num Of Elements being too far: ", NumOfTooFarElements
  end if

  F = MaterialDistribution(ILEV)%x
  MaterialDistribution(ILEV)%x = 0
  call MPI_AllReduce(F,MaterialDistribution(ILEV)%x, mg_mesh%level(ILEV)%nel, MPI_INTEGER, MPI_SUM, MPI_COMM_SUBS, ierr)
  
  do i=1, mg_mesh%level(ILEV)%nel
   if (MaterialDistribution(ILEV)%x(i).gt.0) then
    MaterialDistribution(ILEV)%x(i) = MaterialDistribution(ILEV)%x(i) / subnodes
   end if
   if (MaterialDistribution(ILEV)%x(i).lt.0) then
    MaterialDistribution(ILEV)%x(i) = abs(MaterialDistribution(ILEV)%x(i))
   end if
  end do
  
  do i=1, mg_mesh%level(ILEV)%nel
   IF ((E(i).lt.0).and.(MaterialDistribution(ILEV)%x(i).eq.0)) then
    MaterialDistribution(ILEV)%x(i) = E(i)
   end if
  end do
  
  if (NumOfTooFarElements.eq.NumOfTooFarElementsOld) GOTO 2
  NumOfTooFarElementsOld = NumOfTooFarElements
  
  if (NumOfTooFarElements.gt.0) GOTO 1
 
2 DeALLOCATE(E,F)
 
end if

 
 if (myid.eq.1) THEN
  myMultiMat%nOfMaterials = mySigma%NumberOfMat
  GenLinScalar%cName = "Temper"
  GenLinScalar%prm%nOfFields = myMultiMat%nOfMaterials
  GenLinScalar%nOfFields = myMultiMat%nOfMaterials
  ALLOCATE(GenLinScalar%prm%cField(GenLinScalar%prm%nOfFields))
  ALLOCATE(GenLinScalar%Fld(GenLinScalar%prm%nOfFields))
  
  DO iFld=1,myMultiMat%nOfMaterials
   write(GenLinScalar%prm%cField(iFld),'(A,I0)') 'alpha',iFld
  end do

  ALLOCATE(mg_mesh%level(ILEV)%dvol(mg_Mesh%level(ilev)%nel))
  mg_mesh%level(ILEV)%dvol = 1d0
  
  DO iFld = 1,myMultiMat%nOfMaterials
   ILEV = NLMAX
   GenLinScalar%Fld(iFld)%cName = TRIM(GenLinScalar%prm%cField(iFld))
   ALLOCATE(GenLinScalar%Fld(iFld)%aux(ndof))
   ALLOCATE(GenLinScalar%Fld(iFld)%val(ndof))
   GenLinScalar%Fld(iFld)%aux=0d0
   GenLinScalar%Fld(iFld)%val=0d0
   
  call INT_P0toQ2(MaterialDistribution(ILEV)%x,&
                      GenLinScalar%Fld(iFld)%val,&
                      GenLinScalar%Fld(iFld)%aux,&
                      mg_mesh%level(ILEV)%dvol,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%nvt,&
                      mg_mesh%level(ILEV)%net,&
                      mg_mesh%level(ILEV)%nat,&
                      mg_mesh%level(ILEV)%nel,iFld)
                      
   cField = adjustl(trim(GenLinScalar%Fld(iFld)%cName))
   call subWRITEMPIFieldQ2_NX(cField,iOut,1,GenLinScalar%Fld(iFld)%val)
  END DO  
 end if 

end subroutine OUTPUT_xParticleDUMP
!========================================================================================
!                           Sub: init_q2p1_xParicles
!========================================================================================
SUBROUTINE OUTPUT_xParticleVTK(iOut)
 USE def_FEAT
 USE var_QuadScalar, ONLY : mg_Mesh,QuadSc,myExport
 USE PP3D_MPI, ONLY : myid,master,showid

 implicit none
 integer iOut

 if (myid.eq.1) THEN
  NLMAX = NLMAX + 1
  ILEV = myExport%Level
  myExport%Format = 'VTK'
  allocate(myExport%Fields(5))
  myExport%Fields(1) = 'Velocity'
  myExport%Fields(2) = 'Shell'
  myExport%Fields(3) = 'Material_E'
  myExport%Fields(4) = 'GenScalar'
  myExport%Fields(5) = 'ParticleTime'
  CALL Output_VTK_piece(iOut,mg_mesh%level(ILEV)%dcorvg,mg_mesh%level(ILEV)%kvert)   
  NLMAX = NLMAX - 1
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
 character*(256) :: cTRIA
 character*(256) :: cArray(7)

 CALL INIT_MPI()
 
 CALL ReadParameters(cParamFile)
 
 mg_Mesh%nlmax = NLMAX
 mg_Mesh%nlmin = NLMIN
 mg_Mesh%maxlevel = mg_Mesh%nlmax+1
 allocate(mg_mesh%level(mg_Mesh%maxlevel))
 
 call readTriCoarse(cMeshFile, mg_mesh)
 
if (ADJUSTL(TRIM(xProcess)).eq.'TRACE_TRIA') THEN
 cTRIA='_data/tria'
 cArray = ["DCORVG             ",&
           "KVERT              ",&
           "KEDGE              ",&
           "KAREA              ",&
           "KADJ               ",&
           "elementsAtVertexIdx",&
           "elementsAtVertex   "]
 call read_tria(mg_mesh, mg_Mesh%maxlevel,cTRIA,cArray,7)  
END IF
  
if (ADJUSTL(TRIM(xProcess)).eq.'TRACE'.or.&
    ADJUSTL(TRIM(xProcess)).eq.'WRITE_TRIA') THEN
 call refineMesh(mg_mesh, mg_Mesh%maxlevel,.TRUE.)  
end if
 
if (ADJUSTL(TRIM(xProcess)).eq.'WRITE_TRIA') THEN
 if (myid.eq.1) then
  cTRIA='_data/tria'
  cArray = ["DCORVG             ",&
            "KVERT              ",&
            "KEDGE              ",&
            "KAREA              ",&
            "KADJ               ",&
            "elementsAtVertexIdx",&
            "elementsAtVertex   "]
  call write_tria(mg_mesh, mg_Mesh%maxlevel,cTRIA,cArray,7)
 end if
 CALL Barrier_myMPI()
 CALL MPI_Finalize(ierr)
 stop
end if

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
use Sigma_User, only : myProcess,mySigma

implicit none

type(t_parlist) :: parameterlist
integer :: unitProtfile = -1 ! I guess you use mfile here
integer :: unitTerminal = 6 ! I guess you use mterm here
character(len=INIP_STRLEN) cParserString

if (myid.eq.1) WRITE(*,*) "Parameter File:",ADJUSTL(TRIM(cParamFile))

call inip_output_init(myid,showid,unitProtfile,unitTerminal)

! Init the parameterlist
call inip_init(parameterlist)


call inip_readfromfile(parameterlist,ADJUSTL(TRIM(cParamFile)))

call INIP_getvalue_string(parameterlist,"MAIN","xProcess",xProcess,'TRACE')
call inip_toupper_replace(xProcess)
if (myid.eq.1) WRITE(*,*) "xProcess = '"//ADJUSTL(TRIM(xProcess))//"'"
if (.not.(ADJUSTL(TRIM(xProcess)).eq.'TRACE'.or.&
          ADJUSTL(TRIM(xProcess)).eq.'WRITE_TRIA'.or.&
          ADJUSTL(TRIM(xProcess)).eq.'TRACE_TRIA')) THEN
 write(*,*) "unknown process is being defined // STOP"
END IF

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

call INIP_getvalue_int(parameterlist,"E3DMaterialParameters",   "NoOfMaterials"      ,mySigma%NumberOfMat,1)
if (myid.eq.1) WRITE(*,*) "NumberOfMat = '",mySigma%NumberOfMat

 call INIP_getvalue_int(parameterlist,"E3DProcessParameters",   "nOfInflows"      ,myProcess%nOfInflows,0)
 if (myid.eq.1) WRITE(*,*) "nOfInflows = '",myProcess%nOfInflows
 ALLOCATE(myProcess%myInflow(myProcess%nOfInflows))
 cParserString = "E3DProcessParameters"
 CALL FillUpInflows(myProcess%nOfInflows,cParserString)

! call inip_toupper_replace(cRotation)




! call inip_info(parameterlist)

! Clean up the parameterlist
call inip_done(parameterlist)

contains

    SUBROUTINE FillUpInflows(nT,cINI)
    use iniparser
    USE PP3D_MPI, ONLY : myid,master,showid
    use Sigma_User, only : myProcess
    use, intrinsic :: ieee_arithmetic
    implicit none
    INTEGER :: nT 
    integer :: iSubInflow,iInflow
    character(len=INIP_STRLEN) cINI,cUnit,cCenter,cNormal,cInflow_i,cBCtype
    character(len=INIP_STRLEN) cMidpointA, cMidpointB 
    real*8 daux
    real*8 :: myInf
!     type(t_parlist) :: parameterlist
    
    if(ieee_support_inf(myInf))then
      myInf = ieee_value(myInf, ieee_negative_inf)
    endif
    
    DO iInflow=1,myProcess%nOfInflows
!      if (iInflow.gt.0.and.iInflow.le.9) WRITE(cInflow_i,'(A,I1.1)') 'E3DProcessParameters/Inflow_',iInflow
!      if (iInflow.gt.9.and.iInflow.le.99) WRITE(cInflow_i,'(A,I2.2)') 'E3DProcessParameters/Inflow_',iInflow
!      if (iInflow.gt.99.and.iInflow.le.999) WRITE(cInflow_i,'(A,I3.3)') 'E3DProcessParameters/Inflow_',iInflow
     
     WRITE(cInflow_i,'(A,A,I0)') ADJUSTL(TRIM(cINI)),"/Inflow_",iInflow
     
     if (myid.eq.1) write(*,*) "|",ADJUSTL(TRIM(cInflow_i)),"|"
     
     call INIP_getvalue_string(parameterlist,cInflow_i,"Unit",cUnit,'cm')
     call inip_toupper_replace(cUnit)
     IF (.NOT.(TRIM(cUnit).eq.'MM'.OR.TRIM(cUnit).eq.'CM'.OR.TRIM(cUnit).eq.'DM'.OR.TRIM(cUnit).eq.'M')) THEN
       WRITE(*,*) "Unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
       cUnit = 'cm'
     END IF
     if (TRIM(cUnit).eq.'MM') daux = 0.100d0
     if (TRIM(cUnit).eq.'CM') daux = 1.000d0
     if (TRIM(cUnit).eq.'DM') daux = 10.00d0
     if (TRIM(cUnit).eq.'M')  daux = 100.0d0

     call INIP_getvalue_int(parameterlist,cInflow_i,"nSubInflows",myProcess%myInflow(iInflow)%nSubInflows,0)
     
     IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
     
      call INIP_getvalue_int(parameterlist,cInflow_i,"Material",myProcess%myInflow(iInflow)%Material,0)
      if (myProcess%myInflow(iInflow)%Material.eq.0) write(*,*) 'UNDEFINED material from Inflow ',iInflow,' !!'

      call INIP_getvalue_string(parameterlist,cInflow_i,"Type",cBCtype,'unknown')
      call inip_toupper_replace(cBCtype)
      myProcess%myInflow(iInflow)%iBCtype = 0
      IF (ADJUSTL(TRIM(cBCtype)).eq."ROTATEDPARABOLA1") THEN
       myProcess%myInflow(iInflow)%iBCtype = 1
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."ROTATEDPARABOLA2") THEN
       myProcess%myInflow(iInflow)%iBCtype = 2
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."FLAT") THEN
       myProcess%myInflow(iInflow)%iBCtype = 3
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."CURVEDFLAT") THEN
       myProcess%myInflow(iInflow)%iBCtype = 4
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."RECTANGLE") THEN
       myProcess%myInflow(iInflow)%iBCtype = 5
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."CURVEDRECTANGLE") THEN
       myProcess%myInflow(iInflow)%iBCtype = 6
      END IF
      if (myProcess%myInflow(iInflow)%iBCtype.eq.0) then
       write(*,*) 'UNDEFINED Inflow type!!'
      end if
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"massflowrate",myProcess%myInflow(iInflow)%massflowrate,myInf)
      if (myProcess%myInflow(iInflow)%massflowrate.eq.myInf) write(*,*) 'UNDEFINED massflowrate through Inflow',iInflow,' !!'
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"temperature",myProcess%myInflow(iInflow)%temperature,myInf)
      if (myProcess%myInflow(iInflow)%temperature.eq.myInf) write(*,*) 'UNDEFINED temperature for Inflow',iInflow,' !!'
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"TemperatureRange",myProcess%myInflow(iInflow)%TemperatureRange,0d0)

      call INIP_getvalue_string(parameterlist,cInflow_i,"TemperatureType",cBCtype,'CNST') ! Default is cnst temperature (0)
      call inip_toupper_replace(cBCtype)
      myProcess%myInflow(iInflow)%temperatureType = 0
      IF (ADJUSTL(TRIM(cBCtype)).eq."LINEAR") THEN
       myProcess%myInflow(iInflow)%temperatureType = 1
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."QUADRATIC") THEN
       myProcess%myInflow(iInflow)%temperatureType = 2
      END IF
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"innerradius",myProcess%myInflow(iInflow)%innerradius,myInf)
      myProcess%myInflow(iInflow)%innerradius = daux*myProcess%myInflow(iInflow)%innerradius
      if (myProcess%myInflow(iInflow)%innerradius.eq.myInf) write(*,*) 'UNDEFINED inner radius for Inflow',iInflow,' !!'
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"outerradius",myProcess%myInflow(iInflow)%outerradius,myInf)
      myProcess%myInflow(iInflow)%outerradius = daux*myProcess%myInflow(iInflow)%outerradius
      if (myProcess%myInflow(iInflow)%outerradius.eq.myInf) write(*,*) 'UNDEFINED outer radius for Inflow',iInflow,' !!'
      
      call INIP_getvalue_string(parameterlist,cInflow_i,"center",cCenter,'unknown')
      call INIP_getvalue_string(parameterlist,cInflow_i,"normal",cNormal,'unknown')
      read(cCenter,*,err=55) myProcess%myInflow(iInflow)%Center
      myProcess%myInflow(iInflow)%Center = daux*myProcess%myInflow(iInflow)%Center
      read(cNormal,*,err=56) myProcess%myInflow(iInflow)%Normal

      call INIP_getvalue_string(parameterlist,cInflow_i,"midpointA",cMidpointA,'unknown')
      call INIP_getvalue_string(parameterlist,cInflow_i,"midpointB",cMidpointB,'unknown')
      read(cMidpointA,*,err=57) myProcess%myInflow(iInflow)%MidpointA
      myProcess%myInflow(iInflow)%MidpointA = daux*myProcess%myInflow(iInflow)%MidpointA
      read(cMidpointB,*,err=57) myProcess%myInflow(iInflow)%MidpointB
      myProcess%myInflow(iInflow)%MidpointB = daux*myProcess%myInflow(iInflow)%MidpointB
     ELSE
      
      ALLOCATE(myProcess%myInflow(iInflow)%mySubInflow(myProcess%myInflow(iInflow)%nSubInflows))
      
      DO iSubInflow = 1,myProcess%myInflow(iInflow)%nSubInflows
       WRITE(cInflow_i,'(A,A,I0,A,I0)') ADJUSTL(TRIM(cINI)),"/Inflow_",iInflow,"/Sub_",iSubInflow
!       WRITE(*,*)"'",ADJUSTL(TRIM(cInflow_i)),"'"
      
       call INIP_getvalue_int(parameterlist,cInflow_i,"Material",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Material,0)
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Material.eq.0) write(*,*) 'UNDEFINED material from Inflow ',iInflow,' !!'

       call INIP_getvalue_string(parameterlist,cInflow_i,"Type",cBCtype,'unknown')
       call inip_toupper_replace(cBCtype)
       myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 0
       IF (ADJUSTL(TRIM(cBCtype)).eq."ROTATEDPARABOLA1") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 1
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."ROTATEDPARABOLA2") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 2
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."FLAT") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 3
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."CURVEDFLAT") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 4
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."RECTANGLE") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 5
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."CURVEDRECTANGLE") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 6
       END IF
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype.eq.0) then
        write(*,*) 'UNDEFINED Inflow type!!'
       end if
       
       call INIP_getvalue_double(parameterlist,cInflow_i,"massflowrate",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate,myInf)
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate.eq.myInf) write(*,*) 'UNDEFINED massflowrate through Inflow',iInflow,' !!'
       
       call INIP_getvalue_double(parameterlist,cInflow_i,"temperature",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%temperature,myInf)
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%temperature.eq.myInf) write(*,*) 'UNDEFINED temperature for Inflow',iInflow,' !!'

       call INIP_getvalue_double(parameterlist,cInflow_i,"innerradius",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius,myInf)
       myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius = daux*myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius.eq.myInf) write(*,*) 'UNDEFINED inner radius for Inflow',iInflow,' !!'
       
       call INIP_getvalue_double(parameterlist,cInflow_i,"outerradius",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius,myInf)
       myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius = daux*myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius.eq.myInf) write(*,*) 'UNDEFINED outer radius for Inflow',iInflow,' !!'
       
       call INIP_getvalue_string(parameterlist,cInflow_i,"center",cCenter,'unknown')
       call INIP_getvalue_string(parameterlist,cInflow_i,"normal",cNormal,'unknown')
       read(cCenter,*,err=55) myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Center
       myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Center = daux*myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Center
       read(cNormal,*,err=56) myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Normal
       
      END DO
      
     END IF
     GOTO 57     
55   write(*,*) 'WRONGLY DEFINED center for Inflow',iInflow,' !!'//"|",ADJUSTL(TRIM(cCenter)),"|"
     GOTO 57
56   write(*,*) 'WRONGLY DEFINED normal for Inflow',iInflow,' !!'//"|",ADJUSTL(TRIM(cNormal)),"|"  
     GOTO 57
57   CONTINUE
     
    END DO
    
    if (myid.eq.1) THEN
     write(*,*) "myProcess%nOfInflows",'=',myProcess%nOfInflows
     DO iInflow=1,myProcess%nOfInflows
      if (iInflow.gt.0.and.iInflow.le.9) WRITE(cInflow_i,'(A,I1.1)') 'Inflow_',iInflow
      if (iInflow.gt.9.and.iInflow.le.99) WRITE(cInflow_i,'(A,I2.2)') 'Inflow_',iInflow
      if (iInflow.gt.99.and.iInflow.le.999) WRITE(cInflow_i,'(A,I3.3)') 'Inflow_',iInflow
      
      IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Type','=',myProcess%myInflow(iInflow)%iBCtype
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Material','=',myProcess%myInflow(iInflow)%Material
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Massflowrate','=',myProcess%myInflow(iInflow)%massflowrate
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Temperature','=',myProcess%myInflow(iInflow)%temperature
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_InnerRadius','=',myProcess%myInflow(iInflow)%InnerRadius
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_OuterRadius','=',myProcess%myInflow(iInflow)%OuterRadius
       write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_center'//'=',myProcess%myInflow(iInflow)%center
       write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_normal'//'=',myProcess%myInflow(iInflow)%normal
      ELSE
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'nSubInflows = ',myProcess%myInflow(iInflow)%nSubInflows
       DO iSubInflow=1,myProcess%myInflow(iInflow)%nSubInflows
        write(*,'(A,I0,A)') "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Sub_',iSubinflow
        write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Type','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype
        write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Material','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Material
        write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Massflowrate','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate
        write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Temperature','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%temperature
        write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_InnerRadius','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%InnerRadius
        write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_OuterRadius','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%OuterRadius
        write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_center'//'=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%center
        write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_normal'//'=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%normal
       END DO
      END IF
     END DO 
    END IF
    
    END SUBROUTINE FillUpInflows    

END SUBROUTINE ReadParameters

!================================================================================================
!                                 Sub: writeArrays  
!================================================================================================
subroutine write_tria(mgMesh, maxlevel, cF, cFields, nFields)
  USE PP3D_MPI, ONLY:myid,showid,master,coarse,MPI_SEEK_SET,MPI_REAL8,MPI_MODE_RDONLY,MPI_MAX
  USE PP3D_MPI, ONLY:MPI_COMM_WORLD,MPI_MODE_CREATE,MPI_MODE_WRONLY,MPI_INFO_NULL,mpi_integer,mpi_status_ignore,MPI_COMM_subs,MPI_Offset_kind,mpi_seek_cur,MPI_DOUBLE_PRECISION,subnodes,comm_summn
  use var_QuadScalar
  implicit none

  type(tMultiMesh), intent(in) :: mgMesh
  integer, intent(in) :: maxlevel
  CHARACTER (len = 256), intent(in) :: cF 
  integer :: nfine
  integer :: ilevel,i,nmax,ld00,ld11,iField
  integer :: cunit=1111
  CHARACTER (len = 256) :: cfile, cFMT
  CHARACTER (len = *), intent(in) :: cFields(*)
  integer, intent(in) :: nFields

  nfine = maxlevel-1 

  write(cFile,'(A)') adjustl(trim(cF))//'_header.dat'
  write(*,*) "cFile: ",adjustl(trim(cFile))
  open (unit=cunit,file=cFile)
  write(cunit, '(I0)') maxlevel
  do ilevel =1,maxlevel
   write(cunit, '(8(I0,(" ")))') mgMesh%level(ilevel)%nvt,&
                                 mgMesh%level(ilevel)%net,&
                                 mgMesh%level(ilevel)%nat,&
                                 mgMesh%level(ilevel)%nel,&
                                 mgMesh%level(ilevel)%nvel
  end do
  close(cunit)

  DO iField = 1, nFields
  
   IF (ADJUSTL(TRIM(cFields(iField))).eq.'DCORVG') THEN
    do ilevel =maxlevel,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_dcorvg_',ilevel,'.dat'
     CALL WRITE_MPI_dDATA(mgMesh%level(ilevel)%dcorvg,3*mgMesh%level(ilevel)%nvt,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KVERT') THEN
    do ilevel =1,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_kvert_',ilevel,'.dat'
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%kvert,8*mgMesh%level(ilevel)%nel,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KEDGE') THEN
    do ilevel =1,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_kedge_',ilevel,'.dat'
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%kedge,12*mgMesh%level(ilevel)%nel,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KAREA') THEN
    do ilevel =1,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_karea_',ilevel,'.dat'
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%karea,6*mgMesh%level(ilevel)%nel,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KADJ') THEN
    do ilevel =1,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_kadj_',ilevel,'.dat'
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%kadj,6*mgMesh%level(ilevel)%nel,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'elementsAtVertexIdx') THEN
    do ilevel =1,nfine
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_elementsAtVertexIdx_',ilevel,'.dat'
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%elementsAtVertexIdx,mgMesh%level(ilevel)%nvt+1,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'elementsAtVertex') THEN
    do ilevel =1,nfine
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_elementsAtVertex_',ilevel,'.dat'
     nmax = mgMesh%level(ilevel)%elementsAtVertexIdx(mgMesh%level(ilevel)%nvt+1)-1
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%elementsAtVertex,nmax,cFile)
    end do
   END IF
   
  END DO
 
 CONTAINS  
 
 SUBROUTINE WRITE_MPI_dDATA(dData,nData,cData)
  CHARACTER (len = 256) :: cData
  integer nData
  real*8 dData(*)
  integer ierr,file_handle
  integer(kind=MPI_Offset_kind) :: offset
  
  if (myid.ne.1) return
  write(*,*) "File: '",adjustl(trim(cData))//"' is being released ..."
  call MPI_File_open(1, trim(cData), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)
  offset = 0
  call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, ierr)
  ! Read data from the file collectively
  call MPI_File_write(file_handle, dData, nData, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  ! Close the file
  call MPI_File_close(file_handle, ierr)
 end  SUBROUTINE WRITE_MPI_dDATA
 
 SUBROUTINE WRITE_MPI_kDATA(dData,nData,cData)
  CHARACTER (len = 256) :: cData
  integer nData
  integer dData(*)
  integer ierr,file_handle
  integer(kind=MPI_Offset_kind) :: offset
  
  if (myid.ne.1) return
  write(*,*) "File: '",adjustl(trim(cData))//"' is being released ..."
  call MPI_File_open(1, trim(cData), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)
  offset = 0
  call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, ierr)
  ! Read data from the file collectively
  call MPI_File_write(file_handle, dData, nData, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  ! Close the file
  call MPI_File_close(file_handle, ierr)
 end  SUBROUTINE WRITE_MPI_kDATA
 
end SUBROUTINE WRITE_tria

!================================================================================================
!                                 Sub: writeArrays  
!================================================================================================
subroutine read_tria(mgMesh, maxlevel, cF, cFields, nFields)
  USE PP3D_MPI, ONLY:myid,showid,master,coarse,MPI_SEEK_SET,MPI_REAL8,MPI_MODE_RDONLY,MPI_MAX
  USE PP3D_MPI, ONLY:MPI_COMM_WORLD,MPI_MODE_CREATE,MPI_MODE_WRONLY,MPI_INFO_NULL,mpi_integer,mpi_status_ignore,MPI_COMM_subs,MPI_Offset_kind,mpi_seek_cur,MPI_DOUBLE_PRECISION,subnodes,comm_summn
  use var_QuadScalar
  implicit none

  type(tMultiMesh), intent(inout) :: mgMesh
  integer, intent(in) :: maxlevel
  CHARACTER (len = 256), intent(in) :: cF
  CHARACTER (len = *), intent(in) :: cFields(*)
  integer, intent(in) :: nFields
  
  integer :: nfine
  integer :: ilevel,i,nmax,ld00,ld11,iField
  integer :: cunit=1111
  CHARACTER (len = 256) :: cfile, cFMT

  nfine = maxlevel-1

  write(cFile,'(A)') adjustl(trim(cF))//'_header.dat'
  if (myid.eq.1) write(*,*) "'"//adjustl(trim(cFile))//"' is being read"
  open (unit=cunit,file=cFile)
  read(cunit, *) nmax
  if (maxlevel.gt.nmax) THEN
   WRITE(*,*) 'incompatible TRIA data  // data is available on too coarse level only'
  END IF
  
  do ilevel =1,maxlevel
   read(cunit, *) mgMesh%level(ilevel)%nvt,&
                  mgMesh%level(ilevel)%net,&
                  mgMesh%level(ilevel)%nat,&
                  mgMesh%level(ilevel)%nel,&
                  mgMesh%level(ilevel)%nvel
                  
   mg_mesh%level(ilevel)%nve = 8
   mg_mesh%level(ilevel)%nee = 12
   mg_mesh%level(ilevel)%nae = 6
  end do
  close(cunit)

  DO iField = 1, nFields
  
   IF (ADJUSTL(TRIM(cFields(iField))).eq.'DCORVG') THEN
    do ilevel =maxlevel,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_dcorvg_',ilevel,'.dat'
     if (associated(mgMesh%level(ilevel)%dcorvg)) deallocate(mgMesh%level(ilevel)%dcorvg)
     if (.not.associated(mgMesh%level(ilevel)%dcorvg)) allocate(mgMesh%level(ilevel)%dcorvg(3,mgMesh%level(ilevel)%nvt))
     CALL READ_MPI_dDATA(mgMesh%level(ilevel)%dcorvg,3*mgMesh%level(ilevel)%nvt,cFile)
    end do
    do ilevel=1,maxlevel-1
      if (associated(mgMesh%level(ilevel)%dcorvg)) deallocate(mgMesh%level(ilevel)%dcorvg)
      mgMesh%level(ilevel)%dcorvg => mgMesh%level(maxlevel)%dcorvg
    end do
   END IF
    
    
   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KVERT') THEN
    do ilevel =1,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_kvert_',ilevel,'.dat'
     if (allocated(mgMesh%level(ilevel)%kvert)) deallocate(mgMesh%level(ilevel)%kvert)
     if (.not.allocated(mgMesh%level(ilevel)%kvert)) allocate(mgMesh%level(ilevel)%kvert(8,mgMesh%level(ilevel)%nel))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%kvert,8*mgMesh%level(ilevel)%nel,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KEDGE') THEN
    do ilevel =1,nfine
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_kedge_',ilevel,'.dat'
     if (allocated(mgMesh%level(ilevel)%kedge)) deallocate(mgMesh%level(ilevel)%kedge)
     if (.not.allocated(mgMesh%level(ilevel)%kedge)) allocate(mgMesh%level(ilevel)%kedge(12,mgMesh%level(ilevel)%nel))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%kedge,12*mgMesh%level(ilevel)%nel,cFile)
    end do
   END IF
   
   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KAREA') THEN
    do ilevel =1,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_karea_',ilevel,'.dat'
     if (allocated(mgMesh%level(ilevel)%karea)) deallocate(mgMesh%level(ilevel)%karea)
     if (.not.allocated(mgMesh%level(ilevel)%karea)) allocate(mgMesh%level(ilevel)%karea(6,mgMesh%level(ilevel)%nel))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%karea,6*mgMesh%level(ilevel)%nel,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KADJ') THEN
    do ilevel =1,maxlevel
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_kadj_',ilevel,'.dat'
     if (allocated(mgMesh%level(ilevel)%kadj)) deallocate(mgMesh%level(ilevel)%kadj)
     if (.not.allocated(mgMesh%level(ilevel)%kadj)) allocate(mgMesh%level(ilevel)%kadj(6,mgMesh%level(ilevel)%nel))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%kadj,6*mgMesh%level(ilevel)%nel,cFile)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'elementsAtVertexIdx') THEN
    do ilevel =1,nfine
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_elementsAtVertexIdx_',ilevel,'.dat'
     if (allocated(mgMesh%level(ilevel)%elementsAtVertexIdx)) deallocate(mgMesh%level(ilevel)%elementsAtVertexIdx)
     if (.not.allocated(mgMesh%level(ilevel)%elementsAtVertexIdx)) allocate(mgMesh%level(ilevel)%elementsAtVertexIdx(mgMesh%level(ilevel)%nvt+1))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%elementsAtVertexIdx,mgMesh%level(ilevel)%nvt+1,cFile)
    end do

   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'elementsAtVertex') THEN
    do ilevel =1,nfine
     write(cFile,'(A,I0,A)') adjustl(trim(cF))//'_elementsAtVertex_',ilevel,'.dat'
     nmax = mgMesh%level(ilevel)%elementsAtVertexIdx(mgMesh%level(ilevel)%nvt+1)-1
     if (allocated(mgMesh%level(ilevel)%elementsAtVertex)) deallocate(mgMesh%level(ilevel)%elementsAtVertex)
     if (.not.allocated(mgMesh%level(ilevel)%elementsAtVertex)) allocate(mgMesh%level(ilevel)%elementsAtVertex(nmax))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%elementsAtVertex,nmax,cFile)
    end do
    
   END IF
   
  END DO
  
 CONTAINS  
 
SUBROUTINE READ_MPI_dDATA(dData,nData,cData)
  CHARACTER (len = 256) :: cData
  integer nData
  real*8 dData(*)
  integer ierr,file_handle
  integer(kind=MPI_Offset_kind) :: offset
    
  if (myid.eq.1) write(*,*) "File: '",adjustl(trim(cData))//"' is being loaded ..."
  call MPI_File_open(MPI_COMM_WORLD, trim(cData), MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, ierr)
  offset = 0
  call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, ierr)
  ! Read data from the file collectively
  call MPI_File_read_all(file_handle, dData, nData, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  ! Close the file
  call MPI_File_close(file_handle, ierr)
 end  SUBROUTINE READ_MPI_dDATA
 
 SUBROUTINE READ_MPI_kDATA(dData,nData,cData)
  CHARACTER (len = 256) :: cData
  integer nData
  integer dData(*)
  integer ierr,file_handle
  integer(kind=MPI_Offset_kind) :: offset
    
  if (myid.eq.1) write(*,*) "File: '",adjustl(trim(cData))//"' is being loaded ..."
  call MPI_File_open(MPI_COMM_WORLD, trim(cData), MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, ierr)
  offset = 0
  call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, ierr)
  ! Read data from the file collectively
  call MPI_File_read_all(file_handle, dData, nData, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  ! Close the file
  call MPI_File_close(file_handle, ierr)
 end  SUBROUTINE READ_MPI_kDATA
 
end subroutine read_tria

