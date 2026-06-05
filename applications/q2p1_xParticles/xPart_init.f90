!========================================================================================
!                           Sub: init_q2p1_xParicles
!========================================================================================
subroutine init_q2p1_xParicles(log_unit)
    
  USE def_FEAT
  USE PP3D_MPI, ONLY : myid,master,showid
  USE var_QuadScalar, ONLY : mg_Mesh,QuadSc,myExport,Shell

  integer, intent(in) :: log_unit
  real*8, allocatable, Dimension(:) :: xField,yField,zField
  integer ndof,ierr
  Character cField*(128)
  include 'mpif.h'

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
 USE PP3D_MPI, ONLY : myid
 implicit none
 integer, intent(in) :: iOut

 if (myid.eq.1) then
  call WriteFinalVTU('FINAL.vtu')
 end if

END SUBROUTINE OUTPUT_xParticleVTK
!
!----------------------------------------------
!
SUBROUTINE WriteFinalVTU(filename)
 USE def_FEAT, ONLY : NLMAX
 USE var_QuadScalar, ONLY : mg_Mesh,QuadSc,Shell,GenLinScalar,MaterialDistribution,myExport
 USE types, ONLY : myLostSet
 USE xPart_def, ONLY : DistanceToInflow
 implicit none
 character(*), intent(in) :: filename
 integer :: unitVTU,ILEV,nPoints,nCells,ivt,elemIdx,iFld
 logical :: hasVelocity,hasShell,hasMaterial,hasTime,hasDistance
 character(len=64) :: fieldName

 if (.not.allocated(mg_mesh%level)) return

 ILEV = myExport%Level
 if (ILEV.lt.1 .or. ILEV.gt.NLMAX) ILEV = NLMAX
 if (ILEV.gt.size(mg_mesh%level)) ILEV = size(mg_mesh%level)

 nPoints = mg_mesh%level(ILEV)%nvt
 nCells  = mg_mesh%level(ILEV)%nel
 if (nPoints.le.0 .or. nCells.le.0) return

 hasVelocity = allocated(QuadSc%ValU) .and. allocated(QuadSc%ValV) .and. allocated(QuadSc%ValW)
 if (hasVelocity) then
  hasVelocity = size(QuadSc%ValU).ge.nPoints .and. &
                size(QuadSc%ValV).ge.nPoints .and. &
                size(QuadSc%ValW).ge.nPoints
 end if

 hasShell = allocated(Shell)
 if (hasShell) hasShell = size(Shell).ge.nPoints

 hasMaterial = allocated(MaterialDistribution)
 if (hasMaterial) then
  if (ILEV.gt.size(MaterialDistribution)) then
   hasMaterial = .false.
  else if (.not.allocated(MaterialDistribution(ILEV)%x)) then
   hasMaterial = .false.
  else
   hasMaterial = size(MaterialDistribution(ILEV)%x).ge.nCells
  end if
 end if

 hasTime = allocated(myLostSet)
 if (hasTime) hasTime = size(myLostSet).ge.nCells

 hasDistance = allocated(DistanceToInflow)
 if (hasDistance) hasDistance = size(DistanceToInflow).ge.nCells

 open(newunit=unitVTU,file=adjustl(trim(filename)),status='replace',action='write')
 write(unitVTU,'(A)') '<?xml version="1.0"?>'
 write(unitVTU,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
 write(unitVTU,'(A)') '  <UnstructuredGrid>'
 write(unitVTU,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',nPoints,'" NumberOfCells="',nCells,'">'

 write(unitVTU,'(A)') '      <Points>'
 write(unitVTU,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
 do ivt=1,nPoints
  write(unitVTU,'(3(1X,ES23.16))') mg_mesh%level(ILEV)%dcorvg(:,ivt)
 end do
 write(unitVTU,'(A)') '        </DataArray>'
 write(unitVTU,'(A)') '      </Points>'

 write(unitVTU,'(A)') '      <Cells>'
 write(unitVTU,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
 do elemIdx=1,nCells
  write(unitVTU,'(8(1X,I0))') (mg_mesh%level(ILEV)%kvert(ivt,elemIdx)-1,ivt=1,8)
 end do
 write(unitVTU,'(A)') '        </DataArray>'
 write(unitVTU,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
 do elemIdx=1,nCells
  write(unitVTU,'(1X,I0)') elemIdx*8
 end do
 write(unitVTU,'(A)') '        </DataArray>'
 write(unitVTU,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
 do elemIdx=1,nCells
  write(unitVTU,'(1X,I0)') 12
 end do
 write(unitVTU,'(A)') '        </DataArray>'
 write(unitVTU,'(A)') '      </Cells>'

 write(unitVTU,'(A)') '      <PointData>'
 if (hasVelocity) then
  write(unitVTU,'(A)') '        <DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">'
  do ivt=1,nPoints
   write(unitVTU,'(3(1X,ES23.16))') QuadSc%ValU(ivt),QuadSc%ValV(ivt),QuadSc%ValW(ivt)
  end do
  write(unitVTU,'(A)') '        </DataArray>'
 end if

 if (hasShell) then
  write(unitVTU,'(A)') '        <DataArray type="Float64" Name="Shell" format="ascii">'
  do ivt=1,nPoints
   write(unitVTU,'(1X,ES23.16)') Shell(ivt)
  end do
  write(unitVTU,'(A)') '        </DataArray>'
 end if

 if (allocated(GenLinScalar%Fld)) then
  do iFld=1,size(GenLinScalar%Fld)
   if (.not.allocated(GenLinScalar%Fld(iFld)%val)) cycle
   if (size(GenLinScalar%Fld(iFld)%val).lt.nPoints) cycle
   fieldName = adjustl(trim(GenLinScalar%Fld(iFld)%cName))
   if (len_trim(fieldName).eq.0) cycle
   write(unitVTU,'(A,A,A)') '        <DataArray type="Float64" Name="',trim(fieldName),'" format="ascii">'
   do ivt=1,nPoints
    write(unitVTU,'(1X,ES23.16)') GenLinScalar%Fld(iFld)%val(ivt)
   end do
   write(unitVTU,'(A)') '        </DataArray>'
  end do
 end if
 write(unitVTU,'(A)') '      </PointData>'

 write(unitVTU,'(A)') '      <CellData>'
 if (hasMaterial) then
  write(unitVTU,'(A)') '        <DataArray type="Float64" Name="Material_E" format="ascii">'
  do elemIdx=1,nCells
   write(unitVTU,'(1X,ES23.16)') real(MaterialDistribution(ILEV)%x(elemIdx),8)
  end do
  write(unitVTU,'(A)') '        </DataArray>'
 end if

 if (hasDistance) then
  write(unitVTU,'(A)') '        <DataArray type="Float64" Name="InflowDistance_E" format="ascii">'
  do elemIdx=1,nCells
   write(unitVTU,'(1X,ES23.16)') DistanceToInflow(elemIdx)
  end do
  write(unitVTU,'(A)') '        </DataArray>'
 end if

 if (hasTime) then
  write(unitVTU,'(A)') '        <DataArray type="Float64" Name="ParticleTime" format="ascii">'
  do elemIdx=1,nCells
   write(unitVTU,'(1X,ES23.16)') myLostSet(elemIdx)%time
  end do
  write(unitVTU,'(A)') '        </DataArray>'
 end if

 write(unitVTU,'(A)') '      </CellData>'

 write(unitVTU,'(A)') '    </Piece>'
 write(unitVTU,'(A)') '  </UnstructuredGrid>'
 write(unitVTU,'(A)') '</VTKFile>'
 close(unitVTU)

END SUBROUTINE WriteFinalVTU
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
 integer :: istatus

 CALL INIT_MPI()
 
 CALL ReadParameters()
 CALL ResolveProjectTriFile()
 
 mg_Mesh%nlmax = NLMAX
 mg_Mesh%nlmin = NLMIN
 mg_Mesh%maxlevel = mg_Mesh%nlmax+1
 allocate(mg_mesh%level(mg_Mesh%maxlevel))
 
 call readTriCoarse(cMeshFile, mg_mesh)

 cTRIA='_tria/'
 call EnsureTriaDirectory(cTRIA)
 
if (ADJUSTL(TRIM(xProcess)).eq.'TRACE_TRIA') THEN
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

 CONTAINS

 SUBROUTINE EnsureTriaDirectory(cDir)
  character(len=*), intent(in) :: cDir
  integer :: iisdir,ierrLocal
  character(len=256) :: cDirNoSlash
  external isdirectory
  external mkdir_recursive

  cDirNoSlash = trim(cDir)
  if (len_trim(cDirNoSlash).gt.0) then
   if (cDirNoSlash(len_trim(cDirNoSlash):len_trim(cDirNoSlash)).eq.'/') then
    cDirNoSlash = cDirNoSlash(1:len_trim(cDirNoSlash)-1)
   end if
  end if

  if (myid.eq.1) then
   call isdirectory(trim(cDirNoSlash)//achar(0), iisdir, ierrLocal)
   if (ierrLocal.ne.0 .or. iisdir.eq.0) then
    call mkdir_recursive(trim(cDirNoSlash)//achar(0), istatus)
    if (istatus.ne.0) then
     write(*,*) 'Unable to create TRIA cache directory: ',trim(cDirNoSlash)
     call MPI_Abort(MPI_COMM_WORLD,1,ierrLocal)
    end if
   end if
  end if

  CALL Barrier_myMPI()
 END SUBROUTINE EnsureTriaDirectory

 SUBROUTINE ResolveProjectTriFile()
  character(len=256) :: cProjectFile,cProjectFolder,cLine,cTriFile
  integer :: iunit,ios,iLen
  logical :: bFound

  cProjectFile = '_data/meshDir/file.prj'
  cProjectFolder = '_data/meshDir/'
  bFound = .false.

  if (myid.eq.1) then
   write(*,*) "Project file = '"//trim(cProjectFile)//"'"
  end if

  open(newunit=iunit,file=trim(cProjectFile),status='old',action='read',iostat=ios)
  if (ios.ne.0) then
   write(*,*) 'Unable to open project file: ',trim(cProjectFile)
   call MPI_Abort(MPI_COMM_WORLD,1,ios)
  end if

  do
   read(iunit,'(A)',iostat=ios) cLine
   if (ios.ne.0) exit
   cLine = adjustl(trim(cLine))
   iLen = len_trim(cLine)
   if (iLen.gt.4) then
    if (cLine(iLen-3:iLen).eq.'.tri') then
     cTriFile = cLine
     bFound = .true.
     exit
    end if
   end if
  end do
  close(iunit)

  if (.not.bFound) then
   write(*,*) 'Mesh file was NOT found in project file: ',trim(cProjectFile)
   call MPI_Abort(MPI_COMM_WORLD,1,ios)
  end if

  cMeshFile = trim(cProjectFolder)//trim(cTriFile)
  if (myid.eq.1) then
   write(*,*) "Resolved MeshFile = '"//trim(cMeshFile)//"'"
  end if
 END SUBROUTINE ResolveProjectTriFile

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
USE mpi

implicit none

type(t_parlist) :: parameterlist
integer :: unitProtfile = -1 ! I guess you use mfile here
integer :: unitTerminal = 6 ! I guess you use mterm here
character(len=INIP_STRLEN) cParserString
integer :: ierr

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
if (abs(xFactor-1d0).gt.1d-12) then
 if (myid.eq.1) write(*,*) "ERROR: xFactor must be 1.0 for q2p1_xParticles."
 call MPI_Abort(MPI_COMM_WORLD,1,ierr)
end if

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
  USE PP3D_MPI, ONLY:MPI_COMM_WORLD,MPI_MODE_CREATE,MPI_MODE_WRONLY,MPI_INFO_NULL, &
                     mpi_integer,mpi_status_ignore,MPI_COMM_subs,MPI_Offset_kind, &
                     mpi_seek_cur,MPI_DOUBLE_PRECISION,subnodes,comm_summn
  use var_QuadScalar
  implicit none

  type(tMultiMesh), intent(in) :: mgMesh
  integer, intent(in) :: maxlevel
  CHARACTER (len = 256), intent(in) :: cF 
  integer :: nfine
  integer :: ilevel,nmax,iField
  integer :: n1,n2,totalItems,bytesPerItem,nChunks
  integer, allocatable :: chunkSizes(:)
  CHARACTER (len = 256) :: cfile
  CHARACTER (len = *), intent(in) :: cFields(*)
  integer, intent(in) :: nFields

  nfine = maxlevel-1 

  DO iField = 1, nFields
  
   IF (ADJUSTL(TRIM(cFields(iField))).eq.'DCORVG') THEN
    do ilevel =maxlevel,maxlevel
     n1 = 3
     n2 = mgMesh%level(ilevel)%nvt
     totalItems = n1*n2
     bytesPerItem = 8
     call ComputeChunkLayout(totalItems,bytesPerItem,chunkSizes,nChunks)
     call BuildMetaFileName(cF,'dcorvg',ilevel,cFile)
     call DeleteFileIfExists(cFile)
     CALL WRITE_MPI_dDATA(mgMesh%level(ilevel)%dcorvg,totalItems,cF,'dcorvg',ilevel,chunkSizes,nChunks)
     call WriteMetaInfo(cF,'dcorvg',ilevel,n1,n2,totalItems,mgMesh%level(ilevel)%nvt,&
                        mgMesh%level(ilevel)%net,mgMesh%level(ilevel)%nat,mgMesh%level(ilevel)%nel,&
                        mgMesh%level(ilevel)%nvel,chunkSizes,nChunks)
     deallocate(chunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KVERT') THEN
    do ilevel =1,maxlevel
     n1 = 8
     n2 = mgMesh%level(ilevel)%nel
     totalItems = n1*n2
     bytesPerItem = 4
     call ComputeChunkLayout(totalItems,bytesPerItem,chunkSizes,nChunks)
     call BuildMetaFileName(cF,'kvert',ilevel,cFile)
     call DeleteFileIfExists(cFile)
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%kvert,totalItems,cF,'kvert',ilevel,chunkSizes,nChunks)
     call WriteMetaInfo(cF,'kvert',ilevel,n1,n2,totalItems,mgMesh%level(ilevel)%nvt,&
                        mgMesh%level(ilevel)%net,mgMesh%level(ilevel)%nat,mgMesh%level(ilevel)%nel,&
                        mgMesh%level(ilevel)%nvel,chunkSizes,nChunks)
     deallocate(chunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KEDGE') THEN
    do ilevel =1,maxlevel
     n1 = 12
     n2 = mgMesh%level(ilevel)%nel
     totalItems = n1*n2
     bytesPerItem = 4
     call ComputeChunkLayout(totalItems,bytesPerItem,chunkSizes,nChunks)
     call BuildMetaFileName(cF,'kedge',ilevel,cFile)
     call DeleteFileIfExists(cFile)
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%kedge,totalItems,cF,'kedge',ilevel,chunkSizes,nChunks)
     call WriteMetaInfo(cF,'kedge',ilevel,n1,n2,totalItems,mgMesh%level(ilevel)%nvt,&
                        mgMesh%level(ilevel)%net,mgMesh%level(ilevel)%nat,mgMesh%level(ilevel)%nel,&
                        mgMesh%level(ilevel)%nvel,chunkSizes,nChunks)
     deallocate(chunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KAREA') THEN
    do ilevel =1,maxlevel
     n1 = 6
     n2 = mgMesh%level(ilevel)%nel
     totalItems = n1*n2
     bytesPerItem = 4
     call ComputeChunkLayout(totalItems,bytesPerItem,chunkSizes,nChunks)
     call BuildMetaFileName(cF,'karea',ilevel,cFile)
     call DeleteFileIfExists(cFile)
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%karea,totalItems,cF,'karea',ilevel,chunkSizes,nChunks)
     call WriteMetaInfo(cF,'karea',ilevel,n1,n2,totalItems,mgMesh%level(ilevel)%nvt,&
                        mgMesh%level(ilevel)%net,mgMesh%level(ilevel)%nat,mgMesh%level(ilevel)%nel,&
                        mgMesh%level(ilevel)%nvel,chunkSizes,nChunks)
     deallocate(chunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KADJ') THEN
    do ilevel =1,maxlevel
     n1 = 6
     n2 = mgMesh%level(ilevel)%nel
     totalItems = n1*n2
     bytesPerItem = 4
     call ComputeChunkLayout(totalItems,bytesPerItem,chunkSizes,nChunks)
     call BuildMetaFileName(cF,'kadj',ilevel,cFile)
     call DeleteFileIfExists(cFile)
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%kadj,totalItems,cF,'kadj',ilevel,chunkSizes,nChunks)
     call WriteMetaInfo(cF,'kadj',ilevel,n1,n2,totalItems,mgMesh%level(ilevel)%nvt,&
                        mgMesh%level(ilevel)%net,mgMesh%level(ilevel)%nat,mgMesh%level(ilevel)%nel,&
                        mgMesh%level(ilevel)%nvel,chunkSizes,nChunks)
     deallocate(chunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'elementsAtVertexIdx') THEN
    do ilevel =1,nfine
     n1 = mgMesh%level(ilevel)%nvt + 1
     n2 = 1
     totalItems = n1
     bytesPerItem = 4
     call ComputeChunkLayout(totalItems,bytesPerItem,chunkSizes,nChunks)
     call BuildMetaFileName(cF,'elementsAtVertexIdx',ilevel,cFile)
     call DeleteFileIfExists(cFile)
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%elementsAtVertexIdx,totalItems,cF,'elementsAtVertexIdx',ilevel,chunkSizes,nChunks)
     call WriteMetaInfo(cF,'elementsAtVertexIdx',ilevel,n1,n2,totalItems,mgMesh%level(ilevel)%nvt,&
                        mgMesh%level(ilevel)%net,mgMesh%level(ilevel)%nat,mgMesh%level(ilevel)%nel,&
                        mgMesh%level(ilevel)%nvel,chunkSizes,nChunks)
     deallocate(chunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'elementsAtVertex') THEN
    do ilevel =1,nfine
     nmax = mgMesh%level(ilevel)%elementsAtVertexIdx(mgMesh%level(ilevel)%nvt+1)-1
     n1 = nmax
     n2 = 1
     totalItems = nmax
     bytesPerItem = 4
     call ComputeChunkLayout(totalItems,bytesPerItem,chunkSizes,nChunks)
     call BuildMetaFileName(cF,'elementsAtVertex',ilevel,cFile)
     call DeleteFileIfExists(cFile)
     CALL WRITE_MPI_kDATA(mgMesh%level(ilevel)%elementsAtVertex,totalItems,cF,'elementsAtVertex',ilevel,chunkSizes,nChunks)
     call WriteMetaInfo(cF,'elementsAtVertex',ilevel,n1,n2,totalItems,mgMesh%level(ilevel)%nvt,&
                        mgMesh%level(ilevel)%net,mgMesh%level(ilevel)%nat,mgMesh%level(ilevel)%nel,&
                        mgMesh%level(ilevel)%nvel,chunkSizes,nChunks)
     deallocate(chunkSizes)
    end do
   END IF
   
  END DO
 
 CONTAINS  

 SUBROUTINE ComputeChunkLayout(totalItems,bytesPerItem,chunkSizes,nChunks)
  integer, intent(in) :: totalItems,bytesPerItem
  integer, allocatable, intent(out) :: chunkSizes(:)
  integer, intent(out) :: nChunks
  integer :: itemsPerChunk,remaining,idx

  itemsPerChunk = int(DataSizeThresholdMPI / dble(bytesPerItem))
  if (itemsPerChunk.le.0) itemsPerChunk = 1

  nChunks = totalItems / itemsPerChunk
  if (mod(totalItems,itemsPerChunk).ne.0) nChunks = nChunks + 1
  if (nChunks.le.0) nChunks = 1

  allocate(chunkSizes(nChunks))
  remaining = totalItems
  do idx = 1,nChunks
   chunkSizes(idx) = min(itemsPerChunk,remaining)
   remaining = remaining - chunkSizes(idx)
  end do
 END SUBROUTINE ComputeChunkLayout

 SUBROUTINE BuildMetaFileName(baseDir,fieldName,ilevel,fileName)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: ilevel
  CHARACTER(len=256), intent(out) :: fileName
  write(fileName,'(A,A,"_",I0,".meta")') adjustl(trim(baseDir)),trim(fieldName),ilevel
 END SUBROUTINE BuildMetaFileName

 SUBROUTINE BuildChunkFileName(baseDir,fieldName,ilevel,iChunk,fileName)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: ilevel,iChunk
  CHARACTER(len=256), intent(out) :: fileName
  write(fileName,'(A,A,"_",I0,"_chunk",I0,".dat")') adjustl(trim(baseDir)),trim(fieldName),ilevel,iChunk
 END SUBROUTINE BuildChunkFileName

 SUBROUTINE FailTriaIO(message)
  CHARACTER(len=*), intent(in) :: message
  integer :: ierrAbort
  write(*,*) 'TRIA cache I/O error: ',trim(message)
  call MPI_Abort(MPI_COMM_WORLD,1,ierrAbort)
 END SUBROUTINE FailTriaIO

 SUBROUTINE DeleteFileIfExists(fileName)
  CHARACTER(len=*), intent(in) :: fileName
  logical :: exists
  integer :: unitLocal
  inquire(file=trim(fileName),exist=exists)
  if (.not.exists) return
  open(newunit=unitLocal,file=trim(fileName),status='old')
  close(unitLocal,status='delete')
 END SUBROUTINE DeleteFileIfExists

 SUBROUTINE WriteMetaInfo(baseDir,fieldName,ilevel,n1,n2,totalItems,nvt,net,nat,nel,nvel,chunkSizes,nChunks)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: ilevel,n1,n2,totalItems,nvt,net,nat,nel,nvel,nChunks
  integer, intent(in) :: chunkSizes(nChunks)
  integer :: unitMeta,idx

  if (myid.ne.1) return

  call BuildMetaFileName(baseDir,fieldName,ilevel,cFile)
  call DeleteFileIfExists(cFile)
  open(newunit=unitMeta,file=trim(cFile),status='replace',action='write')
  write(unitMeta,'(A)') 'TRIA_CHUNKED_V1'
  write(unitMeta,'(A)') trim(fieldName)
  write(unitMeta,'(I0)') ilevel
  write(unitMeta,'(3(I0,1X))') n1,n2,totalItems
  write(unitMeta,'(5(I0,1X))') nvt,net,nat,nel,nvel
  write(unitMeta,'(I0)') nChunks
  do idx = 1,nChunks
   write(unitMeta,'(2(I0,1X))') idx-1,chunkSizes(idx)
  end do
  close(unitMeta)
 end SUBROUTINE WriteMetaInfo
 
 SUBROUTINE WRITE_MPI_dDATA(dData,nData,baseDir,fieldName,ilevel,chunkSizes,nChunks)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: nData,ilevel,nChunks
  integer, intent(in) :: chunkSizes(nChunks)
  real*8 dData(*)
  integer ierr,file_handle
  integer(kind=MPI_Offset_kind) :: offset
  integer :: iChunk,startIdx,chunkItems
  
  if (myid.ne.1) return
  startIdx = 1
  do iChunk = 1,nChunks
   chunkItems = chunkSizes(iChunk)
   call BuildChunkFileName(baseDir,fieldName,ilevel,iChunk-1,cFile)
   write(*,*) "File: '",adjustl(trim(cFile))//"' is being released ..."
   call MPI_File_open(1, trim(cFile), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to open chunk file '"//trim(cFile)//"' for writing")
   offset = 0
   call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to seek chunk file '"//trim(cFile)//"'")
   call MPI_File_write(file_handle, dData(startIdx), chunkItems, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to write chunk file '"//trim(cFile)//"'")
   call MPI_File_close(file_handle, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to close chunk file '"//trim(cFile)//"'")
   startIdx = startIdx + chunkItems
  end do
 end  SUBROUTINE WRITE_MPI_dDATA
 
 SUBROUTINE WRITE_MPI_kDATA(dData,nData,baseDir,fieldName,ilevel,chunkSizes,nChunks)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: nData,ilevel,nChunks
  integer, intent(in) :: chunkSizes(nChunks)
  integer dData(*)
  integer ierr,file_handle
  integer(kind=MPI_Offset_kind) :: offset
  integer :: iChunk,startIdx,chunkItems
  
  if (myid.ne.1) return
  startIdx = 1
  do iChunk = 1,nChunks
   chunkItems = chunkSizes(iChunk)
   call BuildChunkFileName(baseDir,fieldName,ilevel,iChunk-1,cFile)
   write(*,*) "File: '",adjustl(trim(cFile))//"' is being released ..."
   call MPI_File_open(1, trim(cFile), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to open chunk file '"//trim(cFile)//"' for writing")
   offset = 0
   call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to seek chunk file '"//trim(cFile)//"'")
   call MPI_File_write(file_handle, dData(startIdx), chunkItems, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to write chunk file '"//trim(cFile)//"'")
   call MPI_File_close(file_handle, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to close chunk file '"//trim(cFile)//"'")
   startIdx = startIdx + chunkItems
  end do
 end  SUBROUTINE WRITE_MPI_kDATA
 
end SUBROUTINE WRITE_tria

!================================================================================================
!                                 Sub: writeArrays  
!================================================================================================
subroutine read_tria(mgMesh, maxlevel, cF, cFields, nFields)
  USE PP3D_MPI, ONLY:myid,showid,master,coarse,MPI_SEEK_SET,MPI_REAL8,MPI_MODE_RDONLY,MPI_MAX
  USE PP3D_MPI, ONLY:MPI_COMM_WORLD,MPI_MODE_CREATE,MPI_MODE_WRONLY,MPI_INFO_NULL, &
                     mpi_integer,mpi_status_ignore,MPI_COMM_subs,MPI_Offset_kind, &
                     mpi_seek_cur,MPI_DOUBLE_PRECISION,subnodes,comm_summn
  use var_QuadScalar
  implicit none

  type(tMultiMesh), intent(inout) :: mgMesh
  integer, intent(in) :: maxlevel
  CHARACTER (len = 256), intent(in) :: cF
  CHARACTER (len = *), intent(in) :: cFields(*)
  integer, intent(in) :: nFields
  
  integer :: nfine
  integer :: ilevel,nmax,iField
  integer :: metaN1,metaN2,metaTotalItems,metaNvt,metaNet,metaNat,metaNel,metaNvel
  integer :: metaNChunks
  integer, allocatable :: metaChunkSizes(:)
  CHARACTER (len = 256) :: cfile

  nfine = maxlevel-1

  DO iField = 1, nFields
  
   IF (ADJUSTL(TRIM(cFields(iField))).eq.'DCORVG') THEN
    do ilevel =maxlevel,maxlevel
     call ReadMetaInfo(cF,'dcorvg',ilevel,metaN1,metaN2,metaTotalItems,metaNvt,metaNet, &
                       metaNat,metaNel,metaNvel,metaChunkSizes,metaNChunks)
     call ApplyLevelInfo(mgMesh%level(ilevel),metaNvt,metaNet,metaNat,metaNel,metaNvel)
     if (associated(mgMesh%level(ilevel)%dcorvg)) deallocate(mgMesh%level(ilevel)%dcorvg)
     if (.not.associated(mgMesh%level(ilevel)%dcorvg)) allocate(mgMesh%level(ilevel)%dcorvg(metaN1,metaN2))
     CALL READ_MPI_dDATA(mgMesh%level(ilevel)%dcorvg,metaTotalItems,cF,'dcorvg',ilevel,metaChunkSizes,metaNChunks)
     deallocate(metaChunkSizes)
    end do
    do ilevel=1,maxlevel-1
      if (associated(mgMesh%level(ilevel)%dcorvg)) deallocate(mgMesh%level(ilevel)%dcorvg)
      mgMesh%level(ilevel)%dcorvg => mgMesh%level(maxlevel)%dcorvg
    end do
   END IF
    
   
   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KVERT') THEN
    do ilevel =1,maxlevel
     call ReadMetaInfo(cF,'kvert',ilevel,metaN1,metaN2,metaTotalItems,metaNvt,metaNet, &
                       metaNat,metaNel,metaNvel,metaChunkSizes,metaNChunks)
     call ApplyLevelInfo(mgMesh%level(ilevel),metaNvt,metaNet,metaNat,metaNel,metaNvel)
     if (allocated(mgMesh%level(ilevel)%kvert)) deallocate(mgMesh%level(ilevel)%kvert)
     if (.not.allocated(mgMesh%level(ilevel)%kvert)) allocate(mgMesh%level(ilevel)%kvert(metaN1,metaN2))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%kvert,metaTotalItems,cF,'kvert',ilevel,metaChunkSizes,metaNChunks)
     deallocate(metaChunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KEDGE') THEN
    do ilevel =1,nfine
     call ReadMetaInfo(cF,'kedge',ilevel,metaN1,metaN2,metaTotalItems,metaNvt,metaNet, &
                       metaNat,metaNel,metaNvel,metaChunkSizes,metaNChunks)
     call ApplyLevelInfo(mgMesh%level(ilevel),metaNvt,metaNet,metaNat,metaNel,metaNvel)
     if (allocated(mgMesh%level(ilevel)%kedge)) deallocate(mgMesh%level(ilevel)%kedge)
     if (.not.allocated(mgMesh%level(ilevel)%kedge)) allocate(mgMesh%level(ilevel)%kedge(metaN1,metaN2))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%kedge,metaTotalItems,cF,'kedge',ilevel,metaChunkSizes,metaNChunks)
     deallocate(metaChunkSizes)
    end do
   END IF
   
   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KAREA') THEN
    do ilevel =1,maxlevel
     call ReadMetaInfo(cF,'karea',ilevel,metaN1,metaN2,metaTotalItems,metaNvt,metaNet, &
                       metaNat,metaNel,metaNvel,metaChunkSizes,metaNChunks)
     call ApplyLevelInfo(mgMesh%level(ilevel),metaNvt,metaNet,metaNat,metaNel,metaNvel)
     if (allocated(mgMesh%level(ilevel)%karea)) deallocate(mgMesh%level(ilevel)%karea)
     if (.not.allocated(mgMesh%level(ilevel)%karea)) allocate(mgMesh%level(ilevel)%karea(metaN1,metaN2))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%karea,metaTotalItems,cF,'karea',ilevel,metaChunkSizes,metaNChunks)
     deallocate(metaChunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'KADJ') THEN
    do ilevel =1,maxlevel
     call ReadMetaInfo(cF,'kadj',ilevel,metaN1,metaN2,metaTotalItems,metaNvt,metaNet, &
                       metaNat,metaNel,metaNvel,metaChunkSizes,metaNChunks)
     call ApplyLevelInfo(mgMesh%level(ilevel),metaNvt,metaNet,metaNat,metaNel,metaNvel)
     if (allocated(mgMesh%level(ilevel)%kadj)) deallocate(mgMesh%level(ilevel)%kadj)
     if (.not.allocated(mgMesh%level(ilevel)%kadj)) allocate(mgMesh%level(ilevel)%kadj(metaN1,metaN2))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%kadj,metaTotalItems,cF,'kadj',ilevel,metaChunkSizes,metaNChunks)
     deallocate(metaChunkSizes)
    end do
   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'elementsAtVertexIdx') THEN
    do ilevel =1,nfine
     call ReadMetaInfo(cF,'elementsAtVertexIdx',ilevel,metaN1,metaN2,metaTotalItems,metaNvt,metaNet, &
                       metaNat,metaNel,metaNvel,metaChunkSizes,metaNChunks)
     call ApplyLevelInfo(mgMesh%level(ilevel),metaNvt,metaNet,metaNat,metaNel,metaNvel)
     if (allocated(mgMesh%level(ilevel)%elementsAtVertexIdx)) deallocate(mgMesh%level(ilevel)%elementsAtVertexIdx)
     if (.not.allocated(mgMesh%level(ilevel)%elementsAtVertexIdx)) allocate(mgMesh%level(ilevel)%elementsAtVertexIdx(metaN1))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%elementsAtVertexIdx,metaTotalItems,cF,'elementsAtVertexIdx',ilevel,metaChunkSizes,metaNChunks)
     deallocate(metaChunkSizes)
    end do

   END IF

   IF (ADJUSTL(TRIM(cFields(iField))).eq.'elementsAtVertex') THEN
    do ilevel =1,nfine
     call ReadMetaInfo(cF,'elementsAtVertex',ilevel,metaN1,metaN2,metaTotalItems,metaNvt,metaNet, &
                       metaNat,metaNel,metaNvel,metaChunkSizes,metaNChunks)
     call ApplyLevelInfo(mgMesh%level(ilevel),metaNvt,metaNet,metaNat,metaNel,metaNvel)
     nmax = metaTotalItems
     if (allocated(mgMesh%level(ilevel)%elementsAtVertex)) deallocate(mgMesh%level(ilevel)%elementsAtVertex)
     if (.not.allocated(mgMesh%level(ilevel)%elementsAtVertex)) allocate(mgMesh%level(ilevel)%elementsAtVertex(nmax))
     CALL READ_MPI_kDATA(mgMesh%level(ilevel)%elementsAtVertex,nmax,cF,'elementsAtVertex',ilevel,metaChunkSizes,metaNChunks)
     deallocate(metaChunkSizes)
    end do
    
   END IF
   
  END DO
  
 CONTAINS  

 SUBROUTINE BuildMetaFileName(baseDir,fieldName,ilevel,fileName)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: ilevel
  CHARACTER(len=256), intent(out) :: fileName
  write(fileName,'(A,A,"_",I0,".meta")') adjustl(trim(baseDir)),trim(fieldName),ilevel
 END SUBROUTINE BuildMetaFileName

 SUBROUTINE BuildChunkFileName(baseDir,fieldName,ilevel,iChunk,fileName)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: ilevel,iChunk
  CHARACTER(len=256), intent(out) :: fileName
  write(fileName,'(A,A,"_",I0,"_chunk",I0,".dat")') adjustl(trim(baseDir)),trim(fieldName),ilevel,iChunk
 END SUBROUTINE BuildChunkFileName

 SUBROUTINE FailTriaIO(message)
  CHARACTER(len=*), intent(in) :: message
  integer :: ierrAbort
  write(*,*) 'TRIA cache I/O error: ',trim(message)
  call MPI_Abort(MPI_COMM_WORLD,1,ierrAbort)
 END SUBROUTINE FailTriaIO

 SUBROUTINE ApplyLevelInfo(level,metaNvt,metaNet,metaNat,metaNel,metaNvel)
  type(tMesh), intent(inout) :: level
  integer, intent(in) :: metaNvt,metaNet,metaNat,metaNel,metaNvel
  level%nvt = metaNvt
  level%net = metaNet
  level%nat = metaNat
  level%nel = metaNel
  level%nvel = metaNvel
  level%nve = 8
  level%nee = 12
  level%nae = 6
 END SUBROUTINE ApplyLevelInfo

 SUBROUTINE ReadMetaInfo(baseDir,fieldName,ilevel,metaN1,metaN2,metaTotalItems,metaNvt,metaNet, &
                         metaNat,metaNel,metaNvel,chunkSizes,metaNChunks)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: ilevel
  integer, intent(out) :: metaN1,metaN2,metaTotalItems,metaNvt,metaNet,metaNat,metaNel
  integer, intent(out) :: metaNvel,metaNChunks
  integer, allocatable, intent(out) :: chunkSizes(:)
  integer :: unitMeta,idx,chunkId
  logical :: exists
  CHARACTER(len=256) :: magic,fieldInMeta

  call BuildMetaFileName(baseDir,fieldName,ilevel,cFile)
  inquire(file=trim(cFile),exist=exists)
  if (.not.exists) call FailTriaIO("missing meta file '"//trim(cFile)//"'")

  open(newunit=unitMeta,file=trim(cFile),status='old',action='read')
  read(unitMeta,'(A)') magic
  if (trim(magic).ne.'TRIA_CHUNKED_V1') then
   close(unitMeta)
   call FailTriaIO("unsupported meta format in '"//trim(cFile)//"'")
  end if

  read(unitMeta,'(A)') fieldInMeta
  if (trim(fieldInMeta).ne.trim(fieldName)) then
   close(unitMeta)
   call FailTriaIO("field mismatch in meta file '"//trim(cFile)//"'")
  end if

  read(unitMeta,*) chunkId
  if (chunkId.ne.ilevel) then
   close(unitMeta)
   call FailTriaIO("level mismatch in meta file '"//trim(cFile)//"'")
  end if

  read(unitMeta,*) metaN1,metaN2,metaTotalItems
  read(unitMeta,*) metaNvt,metaNet,metaNat,metaNel,metaNvel
  read(unitMeta,*) metaNChunks
  if (metaNChunks.le.0) then
   close(unitMeta)
   call FailTriaIO("invalid chunk count in '"//trim(cFile)//"'")
  end if

  allocate(chunkSizes(metaNChunks))
  do idx = 1,metaNChunks
   read(unitMeta,*) chunkId,chunkSizes(idx)
   if (chunkId.ne.idx-1) then
    close(unitMeta)
    call FailTriaIO("unexpected chunk index in '"//trim(cFile)//"'")
   end if
   if (chunkSizes(idx).le.0) then
    close(unitMeta)
    call FailTriaIO("non-positive chunk size in '"//trim(cFile)//"'")
   end if
  end do
  close(unitMeta)

  if (sum(chunkSizes).ne.metaTotalItems) then
   call FailTriaIO("chunk size sum mismatch in '"//trim(cFile)//"'")
  end if
 END SUBROUTINE ReadMetaInfo
 
SUBROUTINE READ_MPI_dDATA(dData,nData,baseDir,fieldName,ilevel,chunkSizes,nChunks)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: nData,ilevel,nChunks
  integer, intent(in) :: chunkSizes(nChunks)
  real*8 dData(*)
  integer ierr,file_handle
  integer(kind=MPI_Offset_kind) :: offset
  integer :: iChunk,startIdx,chunkItems
  logical :: exists
    
  startIdx = 1
  do iChunk = 1,nChunks
   chunkItems = chunkSizes(iChunk)
   call BuildChunkFileName(baseDir,fieldName,ilevel,iChunk-1,cFile)
   inquire(file=trim(cFile),exist=exists)
   if (.not.exists) call FailTriaIO("missing chunk file '"//trim(cFile)//"'")
   if (myid.eq.1) write(*,*) "File: '",adjustl(trim(cFile))//"' is being loaded ..."
   call MPI_File_open(MPI_COMM_WORLD, trim(cFile), MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to open chunk file '"//trim(cFile)//"' for reading")
   offset = 0
   call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to seek chunk file '"//trim(cFile)//"'")
   call MPI_File_read_all(file_handle, dData(startIdx), chunkItems, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to read chunk file '"//trim(cFile)//"'")
   call MPI_File_close(file_handle, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to close chunk file '"//trim(cFile)//"'")
   startIdx = startIdx + chunkItems
  end do
 end  SUBROUTINE READ_MPI_dDATA
 
 SUBROUTINE READ_MPI_kDATA(dData,nData,baseDir,fieldName,ilevel,chunkSizes,nChunks)
  CHARACTER(len=*), intent(in) :: baseDir,fieldName
  integer, intent(in) :: nData,ilevel,nChunks
  integer, intent(in) :: chunkSizes(nChunks)
  integer dData(*)
  integer ierr,file_handle
  integer(kind=MPI_Offset_kind) :: offset
  integer :: iChunk,startIdx,chunkItems
  logical :: exists
    
  startIdx = 1
  do iChunk = 1,nChunks
   chunkItems = chunkSizes(iChunk)
   call BuildChunkFileName(baseDir,fieldName,ilevel,iChunk-1,cFile)
   inquire(file=trim(cFile),exist=exists)
   if (.not.exists) call FailTriaIO("missing chunk file '"//trim(cFile)//"'")
   if (myid.eq.1) write(*,*) "File: '",adjustl(trim(cFile))//"' is being loaded ..."
   call MPI_File_open(MPI_COMM_WORLD, trim(cFile), MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to open chunk file '"//trim(cFile)//"' for reading")
   offset = 0
   call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to seek chunk file '"//trim(cFile)//"'")
   call MPI_File_read_all(file_handle, dData(startIdx), chunkItems, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to read chunk file '"//trim(cFile)//"'")
   call MPI_File_close(file_handle, ierr)
   if (ierr.ne.0) call FailTriaIO("unable to close chunk file '"//trim(cFile)//"'")
   startIdx = startIdx + chunkItems
  end do
 end  SUBROUTINE READ_MPI_kDATA
 
end subroutine read_tria
