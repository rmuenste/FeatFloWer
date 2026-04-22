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
  use el_frozen_driver, only: el_force_kernel, el_write_diagnostics, el_apply_forces, &
                              el_enable_buoyancy, el_fluid_density, el_kinematic_viscosity, &
                              el_particle_density, el_print_configuration, el_parse_yes_no

  integer, intent(in) :: log_unit
  integer :: ierr_abort

  include 'mpif.h'


  !-------INIT PHASE-------

  ! Initialization for FEATFLOW
  call General_init_ext(79,log_unit)

  call Init_QuadScalar_Stuctures(log_unit)

  IF(bViscoElastic)call Init_ViscoScalar_Stuctures(log_unit)

  call Init_LinScalar(log_unit)

  call InitCond_LinScalar()

  if (istart.eq.0 .or. len_trim(adjustl(CSTART)).eq.0) then
    if (myid.eq.1) then
      write(mterm,*) 'q2p1_el_frozen_trace requires SimPar@StartingProc /= 0 and a valid SimPar@StartFile.'
      write(log_unit,*) 'q2p1_el_frozen_trace requires SimPar@StartingProc /= 0 and a valid SimPar@StartFile.'
    end if
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr_abort)
  end if

  ! Normal start from inital configuration
  if (istart.eq.0) then
    if (myid.ne.0) call CreateDumpStructures(1)
    call InitCond_QuadScalar()
    IF(bViscoElastic)call IniProf_ViscoScalar()

  ! Start from a solution on the same lvl
  ! with the same number of partitions
  elseif (istart.eq.1) then
    if (myid.ne.0) call CreateDumpStructures(1)
    call SolFromFile(CSTART,1)

  ! Start from a solution on a lower lvl
  ! with the same number of partitions
  elseif (istart.eq.2)then
    ! In order to read in from a lower level
    ! the lower level structures are needed
    if (myid.ne.0) call CreateDumpStructures(0)
    call SolFromFile(CSTART,0)
    call ProlongateSolution()

    ! Now generate the structures for the actual level 
    if (myid.ne.0) call CreateDumpStructures(1)

  ! Start from a solution on the same lvl
  ! with a different number of partitions
  elseif (istart.eq.3) then
    IF (myid.ne.0) CALL CreateDumpStructures(1)
    call SolFromFileRepart(CSTART,1)
  end if

  call report_loaded_frozen_field(log_unit)
  call el_print_configuration(log_unit)

end subroutine init_q2p1_ext

subroutine report_loaded_frozen_field(log_unit)

  use def_FEAT
  use PP3D_MPI, only: myid, showid
  use Transport_Q2P1, only: LinSc, QuadSc
  use var_QuadScalar, only: istep_ns

  implicit none

  integer, intent(in) :: log_unit

  real*8 :: umin_local, umax_local, vmin_local, vmax_local
  real*8 :: wmin_local, wmax_local, pmin_local, pmax_local
  real*8 :: umin_global, umax_global, vmin_global, vmax_global
  real*8 :: wmin_global, wmax_global, pmin_global, pmax_global
  real*8 :: unorm_local, vnorm_local, wnorm_local, pnorm_local
  real*8 :: unorm_global, vnorm_global, wnorm_global, pnorm_global
  real*8 :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8 :: lx, ly, lz
  integer :: ierr

  include 'mpif.h'

  umin_local = minval(QuadSc%ValU)
  umax_local = maxval(QuadSc%ValU)
  vmin_local = minval(QuadSc%ValV)
  vmax_local = maxval(QuadSc%ValV)
  wmin_local = minval(QuadSc%ValW)
  wmax_local = maxval(QuadSc%ValW)
  pmin_local = minval(LinSc%ValP(NLMAX)%x)
  pmax_local = maxval(LinSc%ValP(NLMAX)%x)

  unorm_local = sum(QuadSc%ValU*QuadSc%ValU)
  vnorm_local = sum(QuadSc%ValV*QuadSc%ValV)
  wnorm_local = sum(QuadSc%ValW*QuadSc%ValW)
  pnorm_local = sum(LinSc%ValP(NLMAX)%x*LinSc%ValP(NLMAX)%x)

  call MPI_Allreduce(umin_local, umin_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(umax_local, umax_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(vmin_local, vmin_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(vmax_local, vmax_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(wmin_local, wmin_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(wmax_local, wmax_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(pmin_local, pmin_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(pmax_local, pmax_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  call MPI_Allreduce(unorm_local, unorm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(vnorm_local, vnorm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(wnorm_local, wnorm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(pnorm_local, pnorm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  unorm_global = sqrt(unorm_global)
  vnorm_global = sqrt(vnorm_global)
  wnorm_global = sqrt(wnorm_global)
  pnorm_global = sqrt(pnorm_global)

  call get_global_domain_extents(xmin, xmax, ymin, ymax, zmin, zmax)

  lx = xmax - xmin
  ly = ymax - ymin
  lz = zmax - zmin

  IF (myid.eq.showid) THEN
    write(mterm,'(A)') 'Frozen-field load summary:'
    write(log_unit,'(A)') 'Frozen-field load summary:'
    write(mterm,'(A,ES14.6,A,I0,A,I0)') '  time=', timens, ' | istep_ns=', istep_ns, ' | nlmax=', NLMAX
    write(log_unit,'(A,ES14.6,A,I0,A,I0)') '  time=', timens, ' | istep_ns=', istep_ns, ' | nlmax=', NLMAX
    write(mterm,'(A,2ES14.6,A,ES14.6)') '  X[min,max] = ', xmin, xmax, ' | Lx = ', lx
    write(log_unit,'(A,2ES14.6,A,ES14.6)') '  X[min,max] = ', xmin, xmax, ' | Lx = ', lx
    write(mterm,'(A,2ES14.6,A,ES14.6)') '  Y[min,max] = ', ymin, ymax, ' | Ly = ', ly
    write(log_unit,'(A,2ES14.6,A,ES14.6)') '  Y[min,max] = ', ymin, ymax, ' | Ly = ', ly
    write(mterm,'(A,2ES14.6,A,ES14.6)') '  Z[min,max] = ', zmin, zmax, ' | Lz = ', lz
    write(log_unit,'(A,2ES14.6,A,ES14.6)') '  Z[min,max] = ', zmin, zmax, ' | Lz = ', lz
    write(mterm,'(A,2ES14.6,A,ES14.6)') '  U[min,max] = ', umin_global, umax_global, ' | l2-sum = ', unorm_global
    write(log_unit,'(A,2ES14.6,A,ES14.6)') '  U[min,max] = ', umin_global, umax_global, ' | l2-sum = ', unorm_global
    write(mterm,'(A,2ES14.6,A,ES14.6)') '  V[min,max] = ', vmin_global, vmax_global, ' | l2-sum = ', vnorm_global
    write(log_unit,'(A,2ES14.6,A,ES14.6)') '  V[min,max] = ', vmin_global, vmax_global, ' | l2-sum = ', vnorm_global
    write(mterm,'(A,2ES14.6,A,ES14.6)') '  W[min,max] = ', wmin_global, wmax_global, ' | l2-sum = ', wnorm_global
    write(log_unit,'(A,2ES14.6,A,ES14.6)') '  W[min,max] = ', wmin_global, wmax_global, ' | l2-sum = ', wnorm_global
    write(mterm,'(A,2ES14.6,A,ES14.6)') '  P[min,max] = ', pmin_global, pmax_global, ' | l2-sum = ', pnorm_global
    write(log_unit,'(A,2ES14.6,A,ES14.6)') '  P[min,max] = ', pmin_global, pmax_global, ' | l2-sum = ', pnorm_global
  END IF

end subroutine report_loaded_frozen_field
!
!----------------------------------------------
!
SUBROUTINE General_init_ext(MDATA,MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
   cProjectFolder,cProjectNumber,nUmbrellaSteps,mg_mesh,nInitUmbrellaSteps,&
   bConstForce, ConstForce
 USE Transport_Q2P1, ONLY : Init_QuadScalar,LinSc,QuadSc
 USE Parametrization, ONLY: InitParametrization,ParametrizeBndr,&
     ProlongateParametrization_STRCT,InitParametrization_STRCT,ParametrizeBndryPoints,&
     DeterminePointParametrization_STRCT,ParametrizeBndryPoints_STRCT
 USE Parametrization, ONLY: ParametrizeQ2Nodes
 USE cinterface 
 use dem_query
 use post_utils,  only: sim_finalize

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
 INTEGER ISE,ISA,ISVEL,ISEEL,ISAEL,ISVED,ISAED,ISVAR
 INTEGER ISEAR,ISEVE,ISAVE,ISVBD,ISEBD,ISABD,IDISP
 INTEGER NEL0,NEL1,NEL2
 REAL ttt0,ttt1
 INTEGER II,NDOF,iUmbrella
 LOGICAL BLIN
 CHARACTER CFILE*60 !CFILE1*60,
 INTEGER kSubPart,iSubPart,iPart,LenFile
 CHARACTER command*100,CSimPar*7
 CHARACTER (len = 60) :: afile 
 CHARACTER (len = 60) :: bfile 
 CHARACTER (len = 60) :: ctemp 

 character(50) :: filename
 character(10) :: int_as_string

 INTEGER nLengthV,nLengthE,LevDif
 REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
 logical :: bwait = .true.
 real*8 :: localMax, totalMax

 integer, dimension(1) :: processRanks
 integer :: MPI_W0, MPI_EX0
 integer :: MPI_Comm_EX0, new_comm
 integer :: error_indicator
 integer :: numParticles, ierror, ndims, reorder
 integer, dimension(3) :: dim_size
 logical periods(0:2)
 real*8 :: domain_xmin, domain_xmax, domain_ymin, domain_ymax
 real*8 :: domain_zmin, domain_zmax

 ndims = 3
 dim_size = (/2, 2, 2/)
 periods = (/.false., .false., .false./)
 reorder = 0

 CALL ZTIME(TTT0)


 !=======================================================================
 !     Data input
 !=======================================================================
 !
 CALL INIT_MPI()                                 ! PARALLEL
 CSimPar = "SimPar"
 CALL  myGDATNEW (CSimPar,0)

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

 CALL Init_QuadScalar(mfile)

 IF (myid.EQ.0) NLMAX = LinSc%prm%MGprmIn%MedLev

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

 call refineMesh(mg_mesh, mg_Mesh%maxlevel)  

 filename = "meshL4"

 ! Convert the integer int_as_string to a string
 write(int_as_string, '(I0)') myid

 filename = trim(filename) // "_" // trim(int_as_string) // ".tri"

 II=NLMIN
 !call writeTriFile(mg_mesh%level(NLMAX), filename)
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
 IF (myid.EQ.0) NLMAX = LinSc%prm%MGprmIn%MedLev
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

 IF (myid.eq.1) write(*,*) 'setting up parallel structures for Q2  on level : ',ILEV

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


 subnodes=numnodes-1


 !=======================================================================
 !     Boundary parametrization
 !=======================================================================
 ILEV=NLMIN
 CALL InitParametrization(mg_mesh%level(ILEV),ILEV)
 
 DO ILEV=NLMIN,NLMAX

   CALL ParametrizeBndr(mg_mesh,ilev)

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

 filename = "meshLP"

 filename = trim(filename) // "_" // trim(int_as_string) // ".tri"

 ! Parametrize the highest level
 if(myid.ne.0)then
   CALL ParametrizeBndr(mg_mesh,nlmax+1)
 endif

 ! This part here is responsible for creation of structures enabling the mesh coordinate 
 ! transfer to the master node so that it can create the corresponding matrices
 IF (myid.EQ.0) THEN
   CALL CreateDumpStructures(0)
 ELSE
   LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
   CALL CreateDumpStructures(LevDif)
 END IF

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

call MPI_Barrier(MPI_COMM_WORLD, error_indicator)
IF (myid.eq.1) write(*,*) 'done!'

 !=======================================================================
 !     Set up the rigid body C++ library
 !=======================================================================
 processRanks(1) = 0
 CALL MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_W0, error_indicator)
 CALL MPI_GROUP_EXCL(MPI_W0, 1, processRanks, MPI_EX0, error_indicator)
 CALL MPI_COMM_CREATE(MPI_COMM_WORLD, MPI_EX0, MPI_Comm_EX0, error_indicator)

#ifdef HAVE_PE 
 call get_global_domain_extents(domain_xmin, domain_xmax, domain_ymin, domain_ymax, &
                                domain_zmin, domain_zmax)
 if (myid .ne. 0) then
   call commf2c_el_frozen_trace(MPI_COMM_WORLD, MPI_Comm_Ex0, myid, &
                                domain_xmin, domain_xmax, domain_ymin, domain_ymax, &
                                domain_zmin, domain_zmax)
 end if
#endif

 call init_fc_rigid_body(myid)      

 call MPI_Barrier(MPI_COMM_WORLD, error_indicator)

 !=======================================================================
 !     Debug output: Write parametrized mesh to VTK for visual reference
 !=======================================================================
 IF (myid.eq.1) THEN
   WRITE(MTERM,*) 'Writing parametrized mesh to VTK for debugging...'
   WRITE(MFILE,*) 'Writing parametrized mesh to VTK for debugging...'
 END IF

 ! Create output directory
 IF (myid.eq.0) THEN
   call execute_command_line('mkdir -p _vtk', wait=.true.)
 END IF
 CALL CommBarrier()

 ! Output the parametrized mesh (mesh-only, no field data)
 IF (myid.NE.0) THEN
   CALL Output_VTK_mesh_piece('initial_mesh', &
     mg_mesh%level(mg_Mesh%maxlevel)%dcorvg, &
     mg_mesh%level(mg_Mesh%maxlevel)%kvert)
 ELSE
   CALL Output_VTK_mesh_main('initial_mesh')
 END IF

 IF (myid.eq.1) THEN
   WRITE(MTERM,*) 'Parametrized mesh written to _vtk/initial_mesh.pvtu'
   WRITE(MFILE,*) 'Parametrized mesh written to _vtk/initial_mesh.pvtu'
 END IF


RETURN

END SUBROUTINE General_init_ext
 !
 !-----------------------------------------------------------------------
 !
SUBROUTINE get_global_domain_extents(xmin, xmax, ymin, ymax, zmin, zmax)
  USE def_FEAT
  USE PP3D_MPI, ONLY: myid, showid
  USE var_QuadScalar, ONLY: mg_mesh

  IMPLICIT NONE

  REAL*8, INTENT(OUT) :: xmin, xmax, ymin, ymax, zmin, zmax
  REAL*8 :: xmin_local, xmax_local, ymin_local, ymax_local
  REAL*8 :: zmin_local, zmax_local
  INTEGER :: ierr

  include 'mpif.h'

  IF (myid.ne.0) THEN
    xmin_local = MINVAL(mg_mesh%level(NLMAX)%dcorvg(1,:))
    xmax_local = MAXVAL(mg_mesh%level(NLMAX)%dcorvg(1,:))
    ymin_local = MINVAL(mg_mesh%level(NLMAX)%dcorvg(2,:))
    ymax_local = MAXVAL(mg_mesh%level(NLMAX)%dcorvg(2,:))
    zmin_local = MINVAL(mg_mesh%level(NLMAX)%dcorvg(3,:))
    zmax_local = MAXVAL(mg_mesh%level(NLMAX)%dcorvg(3,:))
  ELSE
    xmin_local = HUGE(1d0)
    ymin_local = HUGE(1d0)
    zmin_local = HUGE(1d0)
    xmax_local = -HUGE(1d0)
    ymax_local = -HUGE(1d0)
    zmax_local = -HUGE(1d0)
  END IF

  CALL MPI_Allreduce(xmin_local, xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  CALL MPI_Allreduce(xmax_local, xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  CALL MPI_Allreduce(ymin_local, ymin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  CALL MPI_Allreduce(ymax_local, ymax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  CALL MPI_Allreduce(zmin_local, zmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  CALL MPI_Allreduce(zmax_local, zmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  IF (myid.eq.showid) THEN
    WRITE(mterm,'(A,6ES14.6)') 'Global CFD domain extents: ', xmin, xmax, ymin, ymax, zmin, zmax
    WRITE(mfile1,'(A,6ES14.6)') 'Global CFD domain extents: ', xmin, xmax, ymin, ymax, zmin, zmax
  END IF
END SUBROUTINE get_global_domain_extents
 !
 !-----------------------------------------------------------------------
 !
SUBROUTINE myGDATNEW (cName,iCurrentStatus)
  USE PP3D_MPI
  use el_frozen_driver, only: el_force_kernel, el_write_diagnostics, el_apply_forces, &
                              el_enable_buoyancy, el_fluid_density, el_kinematic_viscosity, &
                              el_particle_density, el_parse_yes_no
  USE var_QuadScalar, ONLY : myMatrixRenewal,bNonNewtonian,cGridFileName,&
     nSubCoarseMesh,cFBM_File,bTracer,cProjectFile,bMeshAdaptation,&
     myExport,cAdaptedMeshFile,nUmbrellaSteps,bNoOutflow,myDataFile,&
     bViscoElastic,bRefFrame,bConstForce,ConstForce

   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   PARAMETER (NNLEV=9)
   CHARACTER*7 cName
   CHARACTER letter
   INTEGER :: myFile=888
   INTEGER iEnd,iAt,iEq,iLen,iCurrentStatus,istat
   INTEGER iOutShift
   CHARACTER string*500,cVar*7,cPar*25,cLongString*400
   CHARACTER cParam*32,cParam2*20
   LOGICAL bOK,bOutNMAX
   INTEGER, ALLOCATABLE :: iPos(:)

   !-----------------------------------------------------------------------
   !     C O M M O N S 
   !-----------------------------------------------------------------------
   ! *** Standard COMMON blocks
   COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8

   ! *** COMMON blocks for multigrid data management
   COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,&
     ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
   ! *** COMMON blocks for time discretization
   COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
   COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,&
     EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
   COMMON /NSSAV/  INSAV,INSAVN
   COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,&
     IFUSAV,IFPSAV,IFXSAV,IGID,IGMV,IFINIT

   CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
   COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,&
     ISTART,MSTART,CSTART,ISOL,MSOL,CSOL

   SAVE 

   OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action="read",iostat=istat)
   if(istat .ne. 0)then
     write(*,*)'Could not open data file: ',myDataFile
     stop
   end if

   bOutNMAX = .FALSE.

   DO
   READ (UNIT=myFile,FMT='(A500)',IOSTAT=iEnd) string
   IF (iEnd.EQ.-1) EXIT
   CALL StrStuct()
   IF (bOK) THEN

     READ(string(1:iAt-1),*) cVar

     !   IF (myid.eq.1) WRITE(*,*) TRIM(ADJUSTL(cName))//"#"//TRIM(ADJUSTL(cVar)),TRIM(ADJUSTL(cVar)).EQ.TRIM(ADJUSTL(cName))
     IF (TRIM(ADJUSTL(cVar)).EQ.TRIM(ADJUSTL(cName))) THEN

       READ(string(iAt+1:iEq-1),*) cPar
       !     IF (myid.eq.1) WRITE(*,*) TRIM(ADJUSTL(cName))//"#"//TRIM(ADJUSTL(cPar))
       SELECT CASE (TRIM(ADJUSTL(cPar)))

       CASE ("MeshFolder")
         READ(string(iEq+1:),*) cGridFileName
       CASE ("SubMeshNumber")
         READ(string(iEq+1:),*) nSubCoarseMesh
       CASE ("ParticleFile")
         READ(string(iEq+1:),*) cFBM_File
       CASE ("ProjectFile")
         READ(string(iEq+1:),*) cProjectFile
         MMESH1=61
       CASE ("ProtocolFile")
         READ(string(iEq+1:),*) CFILE1
         MFILE1=62
         MFILE=MFILE1
       CASE ("StartingProc")
         READ(string(iEq+1:),*) ISTART
       CASE ("Umbrella")
         READ(string(iEq+1:),*) nUmbrellaSteps
       CASE ("StartFile")
         READ(string(iEq+1:),*) CSTART
         !      iLen = LEN(TRIM(ADJUSTL(CSTART)))
         !      IF     (myid.lt.10 ) THEN
         !       WRITE(CSTART(iLen+1:),'(A,I1)') "00",myid
         !      ELSEIF (myid.lt.100) THEN
         !       WRITE(CSTART(iLen+1:),'(A,I2)') "0",myid
         !      ELSE 
         !       WRITE(CSTART(iLen+1:),'(I3)') myid
         !      END IF
         MSTART=63
       CASE ("LoadAdaptedMesh")
         bMeshAdaptation = .TRUE. 
         READ(string(iEq+1:),*) cAdaptedMeshFile
       CASE ("SolFile")
         READ(string(iEq+1:),*) CSOL
         iLen = LEN(TRIM(ADJUSTL(CSOL)))
         IF     (myid.lt.10 ) THEN
           WRITE(CSOL(iLen+1:),'(A,I1)') "00",myid
         ELSEIF (myid.lt.100) THEN
           WRITE(CSOL(iLen+1:),'(A,I2)') "0",myid
         ELSE 
           WRITE(CSOL(iLen+1:),'(I3)') myid
         END IF
         MSOL=64
         ISOL = 1
       CASE ("MinMeshLevel")
         READ(string(iEq+1:),*) NLMIN
       CASE ("MaxMeshLevel")
         IF (myid.ne.master) THEN
           READ(string(iEq+1:),*) NLMAX
           myExport%LevelMax = NLMAX
         END IF
         !     IF (myid.eq.MASTER) NLMAX=2
       CASE ("TimeScheme")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         THETA = 0.5d0
         IF (TRIM(ADJUSTL(cParam)).EQ."BE") THETA = 1.0d0
         IF (TRIM(ADJUSTL(cParam)).EQ."FE") THETA = 0.0d0
       CASE ("TimeStep")
         READ(string(iEq+1:),*) TSTEP
       CASE ("TimeAdaptivity")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         IADTIM = 0
         IF (TRIM(ADJUSTL(cParam)).EQ."Yes") IADTIM = 1
       CASE ("ELForceKernel")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         el_force_kernel = trim(adjustl(cParam))
       CASE ("ELWriteDiagnostics")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         el_write_diagnostics = el_parse_yes_no(cParam)
       CASE ("ELApplyForces")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         el_apply_forces = el_parse_yes_no(cParam)
       CASE ("ELEnableBuoyancy")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         el_enable_buoyancy = el_parse_yes_no(cParam)
       CASE ("ELFluidDensity")
         READ(string(iEq+1:),*) el_fluid_density
       CASE ("ELKinematicViscosity")
         READ(string(iEq+1:),*) el_kinematic_viscosity
       CASE ("ELParticleDensity")
         READ(string(iEq+1:),*) el_particle_density
       CASE ("Tracer")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         bTracer = .FALSE.
         IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bTracer = .TRUE.
       CASE ("ViscoElastic")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         bViscoElastic = .FALSE.
         IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bViscoElastic = .TRUE.
       CASE ("ReferenceFrame")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         bRefFrame = .FALSE.
         IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bRefFrame = .TRUE.
       CASE ("NoOutflow")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         bNoOutflow = .FALSE.
         IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bNoOutflow = .TRUE.
       CASE ("UseConstantForcing")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         bConstForce = .FALSE.
         IF (TRIM(ADJUSTL(cParam)).EQ."Yes" .OR. &
             TRIM(ADJUSTL(cParam)).EQ."YES") bConstForce = .TRUE.
       CASE ("ConstantForcing")
         READ(string(iEq+1:),*) ConstForce
       CASE ("MinTimeAdapt")
         READ(string(iEq+1:),*) DTMIN
       CASE ("MaxTimeAdapt")
         READ(string(iEq+1:),*) DTMAX
       CASE ("StartSimTime")
         READ(string(iEq+1:),*) TIMENS
       CASE ("MaxSimTime")
         READ(string(iEq+1:),*) TIMEMX
       CASE ("MaxNumStep")
         READ(string(iEq+1:),*) NITNS
       CASE ("BackUpFreq")
         READ(string(iEq+1:),*) INSAV
       CASE ("BackUpNum")
         READ(string(iEq+1:),*) INSAVN
       CASE ("OutputFreq")
         READ(string(iEq+1:),*) DTGMV
       CASE ("MatrixRenewal")
         READ(string(iEq+1:),*) cParam2
         iLoc=INDEX (cParam2, 'M', .TRUE.)+1
         READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%M
         iLoc=INDEX (cParam2, 'D', .TRUE.)+1
         READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%D
         iLoc=INDEX (cParam2, 'K', .TRUE.)+1
         READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%K
         iLoc=INDEX (cParam2, 'C', .TRUE.)+1
         READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%C
         iLoc=INDEX (cParam2, 'S', .TRUE.)+1
         READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%S
       CASE ("FlowType")
         READ(string(iEq+1:),*) cParam2
         bNonNewtonian=.TRUE.
         IF (TRIM(ADJUSTL(cParam2)).EQ."Newtonian") bNonNewtonian=.FALSE.
       CASE ("OutputLevel")
         READ(string(iEq+1:),*) cParam2
         IF (TRIM(ADJUSTL(cParam2)).EQ."MAX") THEN
           bOutNMAX = .TRUE.
           iOutShift = 0
         END IF
         IF (TRIM(ADJUSTL(cParam2)).EQ."MAX+1") THEN
           bOutNMAX = .TRUE.
           iOutShift = 1
         END IF
         IF (TRIM(ADJUSTL(cParam2)).EQ."MAX-1") THEN
           bOutNMAX = .TRUE.
           iOutShift = -1
         END IF
         IF (.NOT.bOutNMAX) THEN
           READ(string(iEq+1:),*) myExport%Level
           myExport%Level = MAX(MIN(myExport%Level,myExport%LevelMax+1),1)
         END IF

       CASE ("OutputFormat")
         READ(string(iEq+1:),*) myExport%Format
       CASE ("OutputFields")
         READ(string(iEq+1:),*) cLongString
         cLongString = ADJUSTL(TRIM(cLongString))
         mylength = LEN(ADJUSTL(TRIM(cLongString)))
         nFields = 0
         DO i=1,mylength
         READ(cLongString(i:i),'(A)') letter
         IF (letter.eq.',') nFields = nFields + 1
         END DO
         IF (ALLOCATED(myExport%Fields)) DEALLOCATE(myExport%Fields)
         IF (ALLOCATED(iPos)) DEALLOCATE(iPos)
         ALLOCATE(myExport%Fields(nFields+1),iPos(nFields+2))
         iPos(1) = 0
         nFields = 0
         DO i=1,mylength
         READ(cLongString(i:i),'(A)') letter
         IF (letter.eq.',') THEN
           nFields = nFields + 1
           iPos(nFields+1) = i
         END IF
         END DO
         iPos(nFields+2) = mylength+1
         DO i=1,nFields+1
         READ(cLongString(iPos(i)+1:iPos(i+1)-1),'(A)') myExport%Fields(i)
         myExport%Fields(i) = ADJUSTL(TRIM(myExport%Fields(i)))
         END DO

       END SELECT

     END IF
   END IF
   END DO

   CLOSE (myFile)

   M     = 1
   MT    = 1
   IGMV   = NLMAX
   THSTEP=TSTEP*THETA
   IF (bOutNMAX) myExport%Level = NLMAX + iOutShift

   IF (myid.eq.showid) THEN
     OPEN (UNIT=mfile1,FILE=cfile1,action="write",status="replace",iostat=istat)
     if(istat .ne. 0)then
       write(*,*)'Could not open protocal file for writing.'
       stop
     end if
   end if

   IF (iCurrentStatus.EQ.0) THEN
     IF (myid.eq.showid) WRITE(UNIT=mterm,FMT=101)
     IF (myid.eq.showid) WRITE(UNIT=mfile,FMT=101)
   END IF

   ! Printout of all loaded parameters

   IF (myid.eq.showid) THEN
     WRITE(mfile,'(A,A)') "Meshfolder = ",cGridFileName
     WRITE(mterm,'(A,A)') "Meshfolder = ",cGridFileName

     WRITE(mfile,'(A,I10)') "nSubCoarseMesh = ",nSubCoarseMesh
     WRITE(mterm,'(A,I10)') "nSubCoarseMesh = ",nSubCoarseMesh

     WRITE(mfile,'(A,A)') "ParticleFile = ",cFBM_File
     WRITE(mterm,'(A,A)') "ParticleFile = ",cFBM_File

     WRITE(mfile,'(A,A)') "ProjectFile = ",CProjectFile
     WRITE(mterm,'(A,A)') "ProjectFile = ",CProjectFile

     WRITE(mfile,'(A,I1)') "StartingProc = ", ISTART
     WRITE(mterm,'(A,I1)') "StartingProc = ", ISTART

     WRITE(mfile,'(A,A)') "StartFile = ",CSTART
     WRITE(mterm,'(A,A)') "StartFile = ",CSTART

     WRITE(mfile,'(A,A)') "SolFile = ",CSOL
     WRITE(mterm,'(A,A)') "SolFile = ",CSOL

     WRITE(mfile,'(A,I1)') "MinMeshLevel = ",NLMIN
     WRITE(mterm,'(A,I1)') "MinMeshLevel = ",NLMIN

     WRITE(mfile,'(A,I1)') "MaxMeshLevel = ",NLMAX
     WRITE(mterm,'(A,I1)') "MaxMeshLevel = ",NLMAX

     WRITE(mfile,'(A,D12.4)') "TimeScheme = ",THETA
     WRITE(mterm,'(A,D12.4)') "TimeScheme = ",THETA

     WRITE(mfile,'(A,D12.4)') "TimeStep = ",TSTEP
     WRITE(mterm,'(A,D12.4)') "TimeStep = ",TSTEP

     WRITE(mfile,'(A,I1)') "TimeAdaptivity = ", IADTIM
     WRITE(mterm,'(A,I1)') "TimeAdaptivity = ", IADTIM

     WRITE(mfile,'(A,A)') "ELForceKernel = ", trim(el_force_kernel)
     WRITE(mterm,'(A,A)') "ELForceKernel = ", trim(el_force_kernel)

     WRITE(mfile,'(A,L1)') "ELWriteDiagnostics = ", el_write_diagnostics
     WRITE(mterm,'(A,L1)') "ELWriteDiagnostics = ", el_write_diagnostics

     WRITE(mfile,'(A,L1)') "ELApplyForces = ", el_apply_forces
     WRITE(mterm,'(A,L1)') "ELApplyForces = ", el_apply_forces

     WRITE(mfile,'(A,L1)') "ELEnableBuoyancy = ", el_enable_buoyancy
     WRITE(mterm,'(A,L1)') "ELEnableBuoyancy = ", el_enable_buoyancy

     WRITE(mfile,'(A,ES12.4)') "ELFluidDensity = ", el_fluid_density
     WRITE(mterm,'(A,ES12.4)') "ELFluidDensity = ", el_fluid_density

     WRITE(mfile,'(A,ES12.4)') "ELKinematicViscosity = ", el_kinematic_viscosity
     WRITE(mterm,'(A,ES12.4)') "ELKinematicViscosity = ", el_kinematic_viscosity

     WRITE(mfile,'(A,ES12.4)') "ELParticleDensity = ", el_particle_density
     WRITE(mterm,'(A,ES12.4)') "ELParticleDensity = ", el_particle_density

     WRITE(mfile,'(A,D12.4)') "MinTimeAdapt = ",DTMIN
     WRITE(mterm,'(A,D12.4)') "MinTimeAdapt = ",DTMIN

     WRITE(mfile,'(A,D12.4)') "MaxTimeAdapt = ",DTMAX
     WRITE(mterm,'(A,D12.4)') "MaxTimeAdapt = ",DTMAX

     WRITE(mfile,'(A,D12.4)') "StartSimTime = ", TIMENS
     WRITE(mterm,'(A,D12.4)') "StartSimTime = ", TIMENS

     WRITE(mfile,'(A,D12.4)') "MaxSimTime = ", TIMEMX
     WRITE(mterm,'(A,D12.4)') "MaxSimTime = ", TIMEMX

     WRITE(mfile,'(A,I10)') "MaxNumStep = ", NITNS
     WRITE(mterm,'(A,I10)') "MaxNumStep = ", NITNS

     IF (bConstForce) THEN
       WRITE(mfile,'(A)') "UseConstantForcing = ON"
       WRITE(mterm,'(A)') "UseConstantForcing = ON"
       WRITE(mfile,'(A,3ES14.4)') "ConstantForcing = ", ConstForce
       WRITE(mterm,'(A,3ES14.4)') "ConstantForcing = ", ConstForce
     ELSE
       WRITE(mfile,'(A)') "UseConstantForcing = OFF"
       WRITE(mterm,'(A)') "UseConstantForcing = OFF"
     END IF

     WRITE(mfile,'(A,I10)') "BackUpFreq = ", INSAV
     WRITE(mterm,'(A,I10)') "BackUpFreq = ", INSAV

     WRITE(mfile,'(A,I10)') "BackUpNum = ", INSAVN
     WRITE(mterm,'(A,I10)') "BackUpNum = ", INSAVN

     WRITE(mfile,'(A,D12.4)') "OutputFreq = ", DTGMV
     WRITE(mterm,'(A,D12.4)') "OutputFreq = ", DTGMV

     WRITE(mfile,'(A,5(A6I1))') "Matrix Renewal scheme : ","  M = ",myMatrixRenewal%M,&
       ", D = ",myMatrixRenewal%D,", K = ",myMatrixRenewal%K,", S = ",myMatrixRenewal%S,&
       ", C = ",myMatrixRenewal%C
     WRITE(mterm,'(A,5(A6I1))') "Matrix Renewal scheme : ","  M = ",myMatrixRenewal%M,&
       ", D = ",myMatrixRenewal%D,", K = ",myMatrixRenewal%K,", S = ",myMatrixRenewal%S,&
       ", C = ",myMatrixRenewal%C

     IF (bMeshAdaptation) THEN 
       WRITE(mfile,'(A,A)') "Use initial Mesh Adaptation file: ",ADJUSTL(TRIM(cAdaptedMeshFile))
       WRITE(mterm,'(A,A)') "Use initial Mesh Adaptation file: ",ADJUSTL(TRIM(cAdaptedMeshFile))
     ELSE
       WRITE(mfile,'(A)') "No Initial Mesh Adaptation"
       WRITE(mterm,'(A)') "No Initial Mesh Adaptation"
     END IF

     WRITE(mfile,'(A,I10)') "Number of Umbrella smoothening steps",nUmbrellaSteps
     WRITE(mterm,'(A,I10)') "Number of Umbrella smoothening steps",nUmbrellaSteps

     IF (bNoOutflow) THEN 
       WRITE(mfile,'(A)') "Matrix modification is to be performed due to the NoOuflow Condition"
       WRITE(mterm,'(A)') "Matrix modification is to be performed due to the NoOuflow Condition"
     END IF

     IF (bTracer) THEN 
       WRITE(mfile,'(A)') "Tracer equation is included"
       WRITE(mterm,'(A)') "Tracer equation is included"
     END IF

     IF (bViscoElastic) THEN 
       WRITE(mfile,'(A)') "Visco-elastic equation is included"
       WRITE(mterm,'(A)') "Visco-elastic equation is included"
     END IF

     IF (bNonNewtonian) THEN 
       WRITE(mfile,'(A)') "FlowType = non-Newtonian"
       WRITE(mterm,'(A)') "FlowType = non-Newtonian"
     ELSE
       WRITE(mfile,'(A)') "FlowType = Newtonian"
       WRITE(mterm,'(A)') "FlowType = Newtonian"
     END IF

     WRITE(mfile,'(4A,I3,3A,100A)') "Exporting into ", "'",myExport%Format,"' on level: ",&
       myExport%Level," Fields: '",("["//TRIM(ADJUSTL(myExport%Fields(i)))//"]",i=1,nFields+1),"'"!,mylength,iPos
     WRITE(mterm,'(4A,I3,3A,100A)') "Exporting into ", "'",myExport%Format,"' on level: ",&
       myExport%Level," Fields: '",("["//TRIM(ADJUSTL(myExport%Fields(i)))//"]",i=1,nFields+1),"'"!,mylength,iPos
   END IF


   !-----------------------------------------------------------------------
   101 FORMAT(/2X,100('=')/&
     2X,"|",10X,"                                                "&
     40X,"|"/&
     2X,"|",10X,"Parellel Q2/P1 FEM Fluid Dynamics code          "&
     40X,"|"/&
     2X,"|",10X,"Developed by:                                   "&
     40X,"|"/&
     2X,"|",10X,"Otto Mierka, Dmitri Kuzmin and Stefan Turek     "&
     40X,"|"/&
     2X,"|",10X,"Developed at:                                   "&
     40X,"|"/&
     2X,"|",10X,"                                                "&
     40X,"|"/&
     2X,"|",10X,"##########  ##      ##      ",&
     "                                                     ",7X,"|"/&
     2X,"|",10X,"#   ##   #  ##      ##      ",&
     "###     ###   ###   ##### #   #  #    #  #   #  ###  ",7X,"|"/&
     2X,"|",10X,"    ##      ##      ##      ",&
     "#  #   #   #  #  #    #   ## ##  #    #  ##  #  #  # ",7X,"|"/&
     2X,"|",10X,"    ##      ##      ##  ####",&
     "#   #  #   #  ###     #   # # #  #    #  # # #  #   #",7X,"|"/&
     2X,"|",10X,"    ##      ##      ##      ",&
     "#  #   #   #  #  #    #   #   #  #    #  #  ##  #  # ",7X,"|"/&
     2X,"|",10X,"    ##      ##      ##      ",&
     "###     ###   #   #   #   #   #   ####   #   #  ###  ",7X,"|"/&
     2X,"|",10X,"    ##        ######        ",&
     "                                                     ",7X,"|"/&
     2X,"|",57X,"            Chair of Mathematics III",5X,"|"/&
     2X,"|",57X,"    Applied Mathematics and Numerics",5X,"|"/&
     2X,"|",57X,"                    Vogelopthsweg 87",5X,"|"/&
     2X,"|",57X,"                      Dortmund 44225",5X,"|"/&
     2X,"|",57X,"                                    ",5X,"|"/&
     2X,"|",10X,"Based on FeatFlow (c)     ",&
     "see also: http://www.featflow.de",30X,"|"/&
     2X,"|",10X,"Correspondance:",73X,"|"/&
     2X,"|",10X," otto.mierka@math.tu-dortmund.de, ",&
     "stefan.turek@math.tu-dortmund.de",22X,"|"/&
     2X,"|",98X,"|"/&
     2X,100('=')/)

 CONTAINS

 SUBROUTINE StrStuct()
   IMPLICIT NONE
   INTEGER i,n

   n = len(string)
   iAt = 0
   iEq = 0
   DO i=1,n
   IF (string(i:i).EQ. '@') iAt = i
   IF (string(i:i).EQ. '=') iEq = i
   END DO

   bOk=.FALSE.
   IF (iAt.ne.0.AND.iEq.ne.0) bOk=.TRUE.

 END SUBROUTINE StrStuct

 END SUBROUTINE myGDATNEW
