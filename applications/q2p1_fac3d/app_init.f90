subroutine init_q2p1_ext(log_unit)
    
  USE def_FEAT
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar,ProlongateSolution, &
    bTracer,bViscoElastic,StaticMeshAdaptation,&
    LinScalar_InitCond
  USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
    Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File,mg_Mesh
  USE prov_dump_config, ONLY : use_prov_dump_io

  integer, intent(in) :: log_unit

  !-------INIT PHASE-------

  ! Initialization for FEATFLOW
  call General_init_ext(79,log_unit)

  call Init_QuadScalar_Stuctures(log_unit)

  IF(bViscoElastic)call Init_ViscoScalar_Stuctures(log_unit)

  call Init_LinScalar(log_unit)

  call InitCond_LinScalar()

  ! Normal start from inital configuration
  if (istart.eq.0) then
    if (myid.ne.0) call CreateDumpStructures(1)
    call InitCond_QuadScalar()
    IF(bViscoElastic)call IniProf_ViscoScalar()

  ! Start from a solution on the same lvl
  ! with the same number of partitions
  elseif (istart.eq.1) then
    if (myid.ne.0) call CreateDumpStructures(1)
    if (use_prov_dump_io) then
      call SolFromFileProv(CSTART,1)
    else
      call SolFromFile(CSTART,1)
    end if

  ! Start from a solution on a lower lvl
  ! with the same number of partitions
  elseif (istart.eq.2)then
    ! In order to read in from a lower level
    ! the lower level structures are needed
    if (myid.ne.0) call CreateDumpStructures(0)
    if (use_prov_dump_io) then
      call SolFromFileProv(CSTART,0)
    else
      call SolFromFile(CSTART,0)
    end if
    call ProlongateSolution()

    ! Now generate the structures for the actual level 
    if (myid.ne.0) call CreateDumpStructures(1)

  ! Start from a solution on the same lvl
  ! with a different number of partitions
  elseif (istart.eq.3) then
    IF (myid.ne.0) CALL CreateDumpStructures(1)
    if (use_prov_dump_io) then
      call SolFromFileRepartProv(CSTART,1)
    else
      call SolFromFileRepart(CSTART,1)
    end if
  end if

end subroutine init_q2p1_ext
!
!----------------------------------------------
!
SUBROUTINE General_init_ext(MDATA,MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
   cProjectFolder,cProjectNumber,mg_mesh,nInitUmbrellaSteps,&
   bApplyFAC3DMeshDeformation,bFAC3D_CylUmbrellaWeight,&
   dFAC3D_CylCenter,dFAC3D_CylRadius,dFAC3D_CylLength,myErrorCode
 USE Transport_Q2P1, ONLY : Init_FAC3D_Handlers, Init_QuadScalar,LinSc,QuadSc
 USE Parametrization, ONLY: InitParametrization,ParametrizeBndr,&
     ProlongateParametrization_STRCT,InitParametrization_STRCT,ParametrizeBndryPoints,&
     DeterminePointParametrization_STRCT,ParametrizeBndryPoints_STRCT
 USE Parametrization, ONLY: ParametrizeQ2Nodes
 USE param_parser, ONLY: GDATNEW
 USE cinterface 

 IMPLICIT NONE
 ! -------------- workspace -------------------
 INTEGER  NNWORK
 PARAMETER (NNWORK=1)
 INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

 INTEGER            :: KWORK(1)
 REAL               :: VWORK(1)
 DOUBLE PRECISION   :: DWORK(NNWORK)

 COMMON       NWORK,IWORK,IWMAX,L,DWORK
 EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
 ! -------------- workspace -------------------

 INTEGER MDATA,MFILE
 INTEGER abortError
 INTEGER ISE,ISA,ISVEL,ISEEL,ISAEL,ISVED,ISAED,ISVAR
 INTEGER ISEAR,ISEVE,ISAVE,ISVBD,ISEBD,ISABD,IDISP
 INTEGER NEL0,NEL1,NEL2
 REAL ttt0,ttt1
 INTEGER II,NDOF
 LOGICAL BLIN
 CHARACTER CFILE*60 !CFILE1*60,
 INTEGER kSubPart,iSubPart,iPart,LenFile
 CHARACTER command*100,CSimPar*7
 CHARACTER (len = 60) :: afile,cXX
 integer iXX
 CHARACTER (len = 60) :: bfile 

 INTEGER nLengthV,nLengthE,LevDif
 REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
 REAL*8 , ALLOCATABLE :: dCylDist(:), dCylWeight(:)
 logical :: bwait = .true.
 
 character(MPI_MAX_PROCESSOR_NAME) :: processor_name
 integer :: name_len

 CHARACTER (len = 60) :: ctemp 
 integer, dimension(1) :: processRanks
 integer :: MPI_W0, MPI_EX0
 integer :: MPI_Comm_EX0
 integer :: error_indicator
 integer :: numParticles


 CALL ZTIME(TTT0)

 !=======================================================================
 !     Data input
 !=======================================================================
 !
 CALL INIT_MPI()                                 ! PARALLEL
 
 CALL FindNodes()
 
 CSimPar = "SimPar"
 CALL  GDATNEW (CSimPar,0)

 IF (bApplyFAC3DMeshDeformation .AND. ISTART.NE.0) THEN
   IF (myid.EQ.master) THEN
     WRITE(*,'(A)') 'FATAL: FAC3D startup mesh deformation cannot be applied on a restart.'
     WRITE(*,'(A)') 'Set SimPar@ApplyFAC3DMeshDeformation = No; dump coordinates are'
     WRITE(*,'(A)') 'already processed.'
     FLUSH(6)
   END IF
   CALL MPI_Abort(MPI_COMM_WORLD,myErrorCode%PARAMETER_CONFLICT,abortError)
 END IF

 CFILE=CFILE1
 MFILE=MFILE1

 !=======================================================================
 !     Grid generation
 !=======================================================================

 CALL CommBarrier()
 
 ! partitioning
#ifdef HAVE_PE
#ifdef PE_SERIAL_MODE
  include 'PartitionReader2.f90'
#else
  include 'PartitionReader.f90'
#endif
#else
  include 'PartitionReader2.f90'
#endif

 call MPI_Get_processor_name(processor_name, name_len, ierr)
 write(*,'(A,I0,A)') 'Hello from MPI process ', myid, ' on processor "'//trim(processor_name)//'" with mesh '//TRIM(ADJUSTL(CMESH1))
 
 CALL Init_QuadScalar(mfile)
 call Init_FAC3D_Handlers()

 bFAC3D_CylUmbrellaWeight = bApplyFAC3DMeshDeformation
 dFAC3D_CylCenter = (/0.5d0, 0.2d0, 0.205d0/)
 dFAC3D_CylRadius = 0.05d0
 dFAC3D_CylLength = 0.4105d0
 if ((myid.eq.0).or.(myid.eq.1)) then
   IF (bApplyFAC3DMeshDeformation) THEN
     WRITE(*,'(A)') 'FAC3D startup mesh deformation: enabled'
   ELSE
     WRITE(*,'(A)') 'FAC3D startup mesh deformation: disabled'
   END IF
 end if

 IF (myid.EQ.0) THEN
  NLMAX = LinSc%prm%MGprmIn%MedLev
 END IF

 if(NLMAX.eq.0)then
   write(*,*)'NLMAX=0 is invalid, exiting...'
 end if

 IF (IER.NE.0) return

 CLOSE(MMESH1)

 ISE=0
 ISA=0
 ISVEL=0
 ISEEL=0
 ISAEL=0
 ISVED=0
 ISAED=0
 ISVAR=0
 ISEAR=0
 ISEVE=0
 ISAVE=0
 ISVBD=0
 ISEBD=0
 ISABD=0
 IDISP=1

 IF (myid.NE.0) NLMAX = NLMAX + 1
 
 IF (myid.EQ.0) then
   mg_Mesh%nlmax = LinSc%prm%MGprmIn%MedLev
   mg_Mesh%nlmin = 1
   mg_Mesh%maxlevel = LinSc%prm%MGprmIn%MedLev+1
   allocate(mg_mesh%level(LinSc%prm%MGprmIn%MedLev+1))
 else
   allocate(mg_mesh%level(NLMAX))
   mg_Mesh%maxlevel = nlmax
   mg_Mesh%nlmax = nlmax-1
   mg_Mesh%nlmin = 1
 end if

 call readTriCoarse(CMESH1, mg_mesh)

 call refineMesh(mg_mesh, mg_Mesh%maxlevel,.true.)  

 
 II=NLMIN
 IF (myid.eq.1) WRITE(*,*) 'setting up general parallel structures on level : ',II

 CALL PARENTCOMM(mg_mesh%level(II)%nat,&
                 mg_mesh%level(II)%nel,&
                 mg_mesh%level(II)%nvt,&
                 mg_mesh%level(II)%dcorvg,&
                 mg_mesh%level(II)%dcorag,&
                 mg_mesh%level(II)%karea,&
                 mg_mesh%level(II)%kvert)


 IF (myid.EQ.0) NLMAX = NLMAX - 1

 DO II=NLMIN+1,NLMAX
 IF (myid.eq.1) WRITE(*,*) 'setting up general parallel structures on level : ',II
 BLIN = .FALSE.

 CALL CREATECOMM(II,&
                 mg_mesh%level(II)%nat,&
                 mg_mesh%level(II)%nel,&
                 mg_mesh%level(II)%nvt,&
                 mg_mesh%level(II)%dcorag,&
                 mg_mesh%level(II)%dcorvg,&
                 mg_mesh%level(II)%kadj,&
                 mg_mesh%level(II)%karea,&
                 mg_mesh%level(II)%kvert,&
                 BLIN)

 END DO

 IF (myid.eq.1) WRITE(*,*) 'setting up general parallel structures : done!'
 IF (myid.EQ.0) THEN
  NLMAX = LinSc%prm%MGprmIn%MedLev
 END IF
 !     THIS PART WILL BUILD THE REQUIRED COMMUNICATION STRUCTURES
 !     ----------------------------------------------------------

 ! Set up the communication structures for the Quadratic element
 ILEV=NLMIN

 IF (myid.eq.1) write(*,*) 'setting up parallel structures for Q2  on level : ',ILEV

 CALL E013_CreateComm_coarse(mg_mesh%level(ILEV)%dcorvg,&
                             mg_mesh%level(ILEV)%dcorag,&
                             mg_mesh%level(ILEV)%kvert,&
                             mg_mesh%level(ILEV)%kedge,&
                             mg_mesh%level(ILEV)%karea,&
                             mg_mesh%level(ILEV)%nvt,&
                             mg_mesh%level(ILEV)%net,&
                             mg_mesh%level(ILEV)%nat,&
                             mg_mesh%level(ILEV)%nel,&
                             LinSc%prm%MGprmIn%MedLev)

 ILEV = LinSc%prm%MGprmIn%MedLev

 CALL Create_GlobalNumbering(mg_mesh%level(ILEV)%dcorvg,&
                             mg_mesh%level(ILEV)%kvert,&
                             mg_mesh%level(ILEV)%kedge,&
                             mg_mesh%level(ILEV)%karea,&
                             mg_mesh%level(ILEV)%nvt,&
                             mg_mesh%level(ILEV)%net,&
                             mg_mesh%level(ILEV)%nat,&
                             mg_mesh%level(ILEV)%nel)

IF (myid.NE.0) NLMAX = NLMAX - 1
DO ILEV=NLMIN+1,NLMAX

 IF (myid.eq.1) write(*,*) '||->>>setting up parallel structures for Q2  on level : ',ILEV

 CALL E013_CreateComm(mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%nvt,&
                      mg_mesh%level(ILEV)%net,&
                      mg_mesh%level(ILEV)%nat,&
                      mg_mesh%level(ILEV)%nel,&
                      LinSc%prm%MGprmIn%MedLev)

 END DO
 IF (myid.eq.1) write(*,*) 'setting up parallel structures for Q2 :  done!'

 CALL ExtraxtParallelPattern()

 NDOF = mg_mesh%level(NLMAX)%nvt + mg_mesh%level(NLMAX)%nat + &
        mg_mesh%level(NLMAX)%nel + mg_mesh%level(NLMAX)%net

 CALL E011_CreateComm(NDOF)
 

!     ----------------------------------------------------------            
call init_fc_rigid_body(myid)      
call FBM_GetParticles()
CALL FBM_ScatterParticles()
!     ----------------------------------------------------------        
 
ILEV=NLMIN
CALL InitParametrization_STRCT(mg_mesh%level(ILEV),ILEV)
 
DO ILEV=NLMIN,NLMAX
   CALL ProlongateParametrization_STRCT(mg_mesh,ilev)
END DO
if(myid.ne.0)then
   CALL ProlongateParametrization_STRCT(mg_mesh,nlmax+1)
endif
 
CALL DeterminePointParametrization_STRCT(mg_mesh,nlmax)
DO ILEV=NLMIN,NLMAX
   CALL ParametrizeBndryPoints_STRCT(mg_mesh,ilev)
!    CALL ProjectPointToSTL(ilev)
   IF (.not.(myid.eq.0.AND.ilev.gt.LinSc%prm%MGprmIn%MedLev)) THEN
    CALL ProlongateCoordinates(mg_mesh%level(ILEV)%dcorvg,&
                               mg_mesh%level(ILEV+1)%dcorvg,&
                               mg_mesh%level(ILEV)%karea,&
                               mg_mesh%level(ILEV)%kvert,&
                               mg_mesh%level(ILEV)%kedge,&
                               mg_mesh%level(ILEV)%nel,&
                               mg_mesh%level(ILEV)%nvt,&
                               mg_mesh%level(ILEV)%net,&
                               mg_mesh%level(ILEV)%nat)
   END IF
END DO

IF (bApplyFAC3DMeshDeformation) THEN
  CALL ApplyFAC3DInitialMeshDeformation(NLMAX,nInitUmbrellaSteps,MFILE,MTERM)
END IF
 
IF (myid.ne.0) THEN
 
  CALL ProlongateCoordinates(mg_mesh%level(nlmax)%dcorvg,&
                             mg_mesh%level(nlmax+1)%dcorvg,&
                             mg_mesh%level(nlmax)%karea,&
                             mg_mesh%level(nlmax)%kvert,&
                             mg_mesh%level(nlmax)%kedge,&
                             mg_mesh%level(nlmax)%nel,&
                             mg_mesh%level(nlmax)%nvt,&
                             mg_mesh%level(nlmax)%net,&
                             mg_mesh%level(nlmax)%nat)
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial mesh smoothening !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! FINAL Projection to NLMAX +1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (myid.ne.0) THEN
   CALL ParametrizeBndryPoints_STRCT(mg_mesh,nlmax+1)
!    CALL ProjectPointToSTL(ilev+1)
END IF
!!!!!!!!!!!!!!!!!!!!!!!!! FINAL Projection to NLMAX +1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!     Debug output: Write parametrized mesh + cylinder distance to VTK
!=======================================================================
IF (myid.eq.1) THEN
  WRITE(MTERM,*) 'Writing parametrized mesh with cylinder distance to VTK...'
  WRITE(MFILE,*) 'Writing parametrized mesh with cylinder distance to VTK...'
END IF

IF (myid.eq.0) THEN
  call execute_command_line('mkdir -p _vtk', wait=.true.)
END IF
CALL CommBarrier()

IF (myid.NE.0) THEN
  ALLOCATE(dCylDist(mg_mesh%level(mg_mesh%maxlevel)%nvt))
  ALLOCATE(dCylWeight(mg_mesh%level(mg_mesh%maxlevel)%nvt))
  CALL ComputeCylDist(mg_mesh%level(mg_mesh%maxlevel)%nvt, &
    mg_mesh%level(mg_mesh%maxlevel)%dcorvg, dCylDist)
  CALL ComputeCylAttractionWeight(mg_mesh%level(mg_mesh%maxlevel)%nvt, &
    dCylDist, dCylWeight)
  CALL Output_VTK_mesh_piece_with_2fields('initial_mesh', &
    mg_mesh%level(mg_mesh%maxlevel)%dcorvg, &
    mg_mesh%level(mg_mesh%maxlevel)%kvert, &
    dCylDist, 'CylinderDistance', &
    dCylWeight, 'AttractionWeight')
  DEALLOCATE(dCylDist, dCylWeight)
ELSE
  CALL Output_VTK_mesh_main_with_2fields('initial_mesh', &
    'CylinderDistance', 'AttractionWeight')
END IF

IF (myid.eq.1) THEN
  WRITE(MTERM,*) 'Mesh + distance written to _vtk/initial_mesh.pvtu'
  WRITE(MFILE,*) 'Mesh + distance written to _vtk/initial_mesh.pvtu'
END IF

! This part here is responsible for creation of structures enabling the mesh coordinate 
! transfer to the master node so that it can create the corresponding matrices
IF (myid.EQ.0) THEN
   CALL CreateDumpStructures(0)
ELSE
   LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
   CALL CreateDumpStructures(LevDif)
END IF
 !     ----------------------------------------------------------            
 !     ----------------------------------------------------------            
 !     ----------------------------------------------------------            
 !     ----------------------------------------------------------            
 !     ----------------------------------------------------------            
 !     ----------------------------------------------------------            

 ILEV = LinSc%prm%MGprmIn%MedLev

 nLengthV = (2**(ILEV-1)+1)**3
 nLengthE = mg_mesh%level(NLMIN)%nel

 ALLOCATE(SendVect(3,nLengthV,nLengthE))

 CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
                             mg_mesh%level(ILEV)%kvert,&
                             nLengthV,&
                             nLengthE,&
                             mg_mesh%level(ILEV)%nel,&
                             mg_mesh%level(ILEV)%nvt)
 DEALLOCATE(SendVect)

 showid = 1

 IF (myid.eq.showid) THEN
   WRITE(MTERM,'(10(2XA8))') 'ILEV','NVT','NAT','NEL','NET','NDOF'
   WRITE(MFILE,'(10(2XA8))') 'ILEV','NVT','NAT','NEL','NET','NDOF'
 END IF

 DO II=NLMIN,NLMAX

 ILEV=II

 NVT=mg_mesh%level(II)%nvt
 NAT=mg_mesh%level(II)%nat
 NET=mg_mesh%level(II)%net
 NEL=mg_mesh%level(II)%nel

 IF (myid.eq.showid) THEN
   WRITE(MTERM,'(10(2XI8))')ILEV,NVT,NAT,NEL,NET,NVT+NAT+NEL+NET
   WRITE(MFILE,'(10(2XI8))')ILEV,NVT,NAT,NEL,NET,NVT+NAT+NEL+NET
 END IF

 if(.not.allocated(mg_mesh%level(II)%dvol))then
   allocate(mg_mesh%level(II)%dvol(NEL+1))
 end if

 CALL  SETARE(mg_mesh%level(II)%dvol,&
              NEL,&
              mg_mesh%level(II)%kvert,&
              mg_mesh%level(II)%dcorvg)

 END DO

 IF (myid.ne.0) THEN
   ILEV=NLMAX +1 

   if(.not.allocated(mg_mesh%level(ILEV)%dvol))then
     allocate(mg_mesh%level(ILEV)%dvol(NEL+1))
   end if

   CALL  SETARE(mg_mesh%level(ILEV)%dvol,&
                NEL,&
                mg_mesh%level(ILEV)%kvert,&
                mg_mesh%level(ILEV)%dcorvg)
                
 END IF

 CALL ZTIME(TTT1)
 TTGRID=TTT1-TTT0

 IF (myid.eq.showid) THEN
   WRITE(MTERM,*)
   WRITE(MFILE,*)
   WRITE(MTERM,*) 'time for grid initialization : ',TTGRID
   WRITE(MFILE,*) 'time for grid initialization : ',TTGRID
   WRITE(MTERM,*)
   WRITE(MFILE,*)
 END IF

 !=======================================================================
 !     Set up the rigid body C++ library
 !=======================================================================
 processRanks(1) = 0
 CALL MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_W0, error_indicator)
 CALL MPI_GROUP_EXCL(MPI_W0, 1, processRanks, MPI_EX0, error_indicator)
 CALL MPI_COMM_CREATE(MPI_COMM_WORLD, MPI_EX0, MPI_Comm_EX0, error_indicator)


#ifdef HAVE_PE 
  if (myid .ne. 0) then
    call commf2c_fsi(MPI_COMM_WORLD, MPI_Comm_Ex0, myid)
  end if
#endif

 RETURN

END SUBROUTINE General_init_ext

!-----------------------------------------------------------------------
! Apply the optional FAC3D startup deformation in the historical order.
! Runtime SimPar@Umbrella smoothing remains an independent setting.
!-----------------------------------------------------------------------
SUBROUTINE ApplyFAC3DInitialMeshDeformation(nLevel,nInitialSteps,mfile,mterm)
  USE PP3D_MPI, ONLY: myid
  USE var_QuadScalar, ONLY: mg_mesh
  USE Parametrization, ONLY: ParametrizeBndryPoints_STRCT
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nLevel,nInitialSteps,mfile,mterm
  INTEGER :: iStep
  INTEGER, PARAMETER :: FAC3D_ATTRACTION_STEPS = 8
  INTEGER, PARAMETER :: FAC3D_POST_SMOOTHING_STEPS = 2
  REAL*8, PARAMETER :: FAC3D_ATTRACTION_OMEGA = 0.2d0

  IF (myid.EQ.1) THEN
    WRITE(mterm,'(A,I0,A)') 'FAC3D initial umbrella smoothing: ',nInitialSteps,' step(s)'
    WRITE(mfile,'(A,I0,A)') 'FAC3D initial umbrella smoothing: ',nInitialSteps,' step(s)'
  END IF
  DO iStep=1,nInitialSteps
    CALL UmbrellaSmoother_STRCT(0d0,1)
  END DO

  IF (myid.NE.0) THEN
    IF (myid.EQ.1) THEN
      WRITE(mterm,'(A,I0,A)') 'FAC3D cylinder attraction: ',FAC3D_ATTRACTION_STEPS,' step(s)'
      WRITE(mfile,'(A,I0,A)') 'FAC3D cylinder attraction: ',FAC3D_ATTRACTION_STEPS,' step(s)'
    END IF
    DO iStep=1,FAC3D_ATTRACTION_STEPS
      CALL CylinderAttraction(mg_mesh%level(nLevel)%nvt, &
        mg_mesh%level(nLevel)%dcorvg,FAC3D_ATTRACTION_OMEGA)
      CALL ParametrizeBndryPoints_STRCT(mg_mesh,nLevel)
    END DO
  END IF

  IF (myid.EQ.1) THEN
    WRITE(mterm,'(A,I0,A)') 'FAC3D post-attraction smoothing: ', &
      FAC3D_POST_SMOOTHING_STEPS,' step(s)'
    WRITE(mfile,'(A,I0,A)') 'FAC3D post-attraction smoothing: ', &
      FAC3D_POST_SMOOTHING_STEPS,' step(s)'
  END IF
  DO iStep=1,FAC3D_POST_SMOOTHING_STEPS
    CALL UmbrellaSmoother_STRCT(0d0,1)
  END DO
END SUBROUTINE ApplyFAC3DInitialMeshDeformation
 !
 !-----------------------------------------------------------------------
 !

 SUBROUTINE ComputeCylDist(nVtx, dcoor, dist)
   USE var_QuadScalar, ONLY : dFAC3D_CylCenter, dFAC3D_CylRadius, dFAC3D_CylLength
   IMPLICIT NONE
   INTEGER, intent(in) :: nVtx
   REAL*8, intent(in)  :: dcoor(3,*)
   REAL*8, intent(out) :: dist(*)
   INTEGER :: i
   REAL*8  :: dx, dy, dzAbs, rxy, zGap, halfLen

   halfLen = 0.5d0 * dFAC3D_CylLength

   DO i = 1, nVtx
     dx    = dcoor(1,i) - dFAC3D_CylCenter(1)
     dy    = dcoor(2,i) - dFAC3D_CylCenter(2)
     dzAbs = ABS(dcoor(3,i) - dFAC3D_CylCenter(3))
     rxy   = SQRT(dx*dx + dy*dy)
     zGap  = dzAbs - halfLen
     dist(i) = MAX(rxy - dFAC3D_CylRadius, zGap)
   END DO

 END SUBROUTINE ComputeCylDist

 SUBROUTINE CylinderAttraction(nVtx, dcorvg, dOmega)
 !
 ! Attracts mesh nodes toward the cylinder surface using distance bands:
 !   Exterior:
 !     [0.0, dBand1):    strong attraction (alpha = 0.6 to 1.0)
 !     [dBand1, dBand2]:  quadratic fade-out
 !     > dBand2:          no movement
 !   Interior:
 !     [0.0, dBandIn):    gentle outward push (alpha linearly fading)
 !     > dBandIn:          no movement
 !
   USE var_QuadScalar, ONLY : dFAC3D_CylCenter, dFAC3D_CylRadius, dFAC3D_CylLength
   IMPLICIT NONE
   INTEGER, intent(in) :: nVtx
   REAL*8, intent(inout) :: dcorvg(3,*)
   REAL*8, intent(in) :: dOmega

   INTEGER :: i
   REAL*8  :: dx, dy, dz, rxy, dist, alpha
   REAL*8  :: nearX, nearY, halfLen, dzAbs, zGap, absDist
   ! Exterior bands
   REAL*8, parameter :: dBand1 = 0.04d0
   REAL*8, parameter :: dBand2 = 0.15d0
   ! Interior band: gentle push, max alpha = 0.3, linear fade over dBandIn
   REAL*8, parameter :: dBandIn   = 0.04d0
   REAL*8, parameter :: dAlphaIn  = 0.45d0

   halfLen = 0.5d0 * dFAC3D_CylLength

   DO i = 1, nVtx
     dx = dcorvg(1,i) - dFAC3D_CylCenter(1)
     dy = dcorvg(2,i) - dFAC3D_CylCenter(2)
     dz = dcorvg(3,i) - dFAC3D_CylCenter(3)
     dzAbs = ABS(dz)
     rxy = SQRT(dx*dx + dy*dy)

     ! Signed distance to cylinder surface (negative = inside)
     zGap = dzAbs - halfLen
     dist = MAX(rxy - dFAC3D_CylRadius, zGap)

     ! Skip nodes whose closest surface point is on the cap, not the barrel
     IF (zGap .GT. rxy - dFAC3D_CylRadius) CYCLE

     ! Nearest point on cylinder barrel (radial projection)
     IF (rxy .GT. 1d-14) THEN
       nearX = dFAC3D_CylCenter(1) + dFAC3D_CylRadius * dx / rxy
       nearY = dFAC3D_CylCenter(2) + dFAC3D_CylRadius * dy / rxy
     ELSE
       CYCLE
     END IF

     IF (dist .GT. 0d0 .AND. dist .LE. dBand2) THEN
       ! --- Exterior attraction ---
       IF (dist .LT. dBand1) THEN
         alpha = 0.6d0 + 0.4d0 * dist / dBand1
       ELSE
         alpha = ((dBand2 - dist) / (dBand2 - dBand1))**2
       END IF
     ELSE IF (dist .LT. 0d0) THEN
       ! --- Interior: gentle outward push ---
       absDist = ABS(dist)
       IF (absDist .GT. dBandIn) CYCLE
       ! Linear fade: strongest near surface, zero at dBandIn
       alpha = dAlphaIn * (1d0 - absDist / dBandIn)
     ELSE
       CYCLE
     END IF

     ! Move node radially toward the barrel surface
     dcorvg(1,i) = dcorvg(1,i) + dOmega * alpha * (nearX - dcorvg(1,i))
     dcorvg(2,i) = dcorvg(2,i) + dOmega * alpha * (nearY - dcorvg(2,i))
   END DO

 END SUBROUTINE CylinderAttraction

 SUBROUTINE ComputeCylAttractionWeight(nVtx, dist, weight)
 !
 ! Computes the attraction weight (alpha * omega) for visualization.
 ! Must match the alpha logic in CylinderAttraction exactly.
 !
   IMPLICIT NONE
   INTEGER, intent(in) :: nVtx
   REAL*8, intent(in)  :: dist(*)
   REAL*8, intent(out) :: weight(*)
   INTEGER :: i
   REAL*8  :: d, alpha, absDist
   REAL*8, parameter :: dBand1 = 0.04d0
   REAL*8, parameter :: dBand2 = 0.15d0
   REAL*8, parameter :: dBandIn  = 0.04d0
   REAL*8, parameter :: dAlphaIn = 0.45d0

   DO i = 1, nVtx
     d = dist(i)
     IF (d .GT. 0d0 .AND. d .LT. dBand1) THEN
       weight(i) = 0.6d0 + 0.4d0 * d / dBand1
     ELSE IF (d .GE. dBand1 .AND. d .LE. dBand2) THEN
       weight(i) = ((dBand2 - d) / (dBand2 - dBand1))**2
     ELSE IF (d .LT. 0d0) THEN
       absDist = ABS(d)
       IF (absDist .LE. dBandIn) THEN
         weight(i) = -dAlphaIn * (1d0 - absDist / dBandIn)
       ELSE
         weight(i) = 0d0
       END IF
     ELSE
       weight(i) = 0d0
     END IF
   END DO

 END SUBROUTINE ComputeCylAttractionWeight
