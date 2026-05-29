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

end subroutine init_q2p1_ext
!
!----------------------------------------------
!
SUBROUTINE set_lubrication_threshold_from_mesh(myid)
  USE PP3D_MPI, ONLY: master
  USE var_QuadScalar, ONLY: mg_Mesh
  USE cinterface, ONLY: set_lubrication_threshold
  USE def_FEAT, ONLY: NLMAX

  IMPLICIT NONE

  include 'mpif.h'

  INTEGER, INTENT(IN) :: myid
  REAL*8 :: h_min, h_elem, dVol, lub_threshold, dLubricationFactor
  REAL*8 :: x(8), y(8), z(8)
  INTEGER :: iel, ive, ierr
  INTEGER :: NEL, NVT
  INTEGER, DIMENSION(:,:), POINTER :: KVERT
  REAL*8, DIMENSION(:,:), POINTER :: DCORVG

#ifdef HAVE_PE
  ! Configuration parameter: lubrication threshold = factor * h_min
  ! User can modify this value based on their simulation requirements
  dLubricationFactor = 1.0d0  ! Default: 3x minimum mesh size

  ! Get mesh data structures
  NEL = mg_Mesh%level(NLMAX)%nel
  NVT = mg_Mesh%level(NLMAX)%nvt
  !KVERT => mg_Mesh%level(NLMAX)%kvert
  !DCORVG => mg_Mesh%level(NLMAX)%dcorvg

  ! Compute minimum element size on this process
  h_min = HUGE(1.0d0)

  DO iel = 1, NEL
    ! Extract vertex coordinates for this element
    DO ive = 1, 8
      x(ive) = mg_Mesh%level(NLMAX)%dcorvg(1, mg_Mesh%level(NLMAX)%kvert(ive, iel))
      y(ive) = mg_Mesh%level(NLMAX)%dcorvg(2, mg_Mesh%level(NLMAX)%kvert(ive, iel))
      z(ive) = mg_Mesh%level(NLMAX)%dcorvg(3, mg_Mesh%level(NLMAX)%kvert(ive, iel))
    END DO

    ! Compute element volume
    CALL GetElemVol(x, y, z, dVol)

    ! Characteristic element size: h = V^(1/3)
    h_elem = dVol**(1.0d0/3.0d0)

    ! Track minimum
    h_min = MIN(h_min, h_elem)
  END DO

  ! MPI reduction to get global minimum across all processes
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, h_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

  ! Compute lubrication threshold
  lub_threshold = dLubricationFactor * h_min

  ! Set the threshold in PE library
  CALL set_lubrication_threshold(lub_threshold)

  ! Output information (master process only)
  IF (myid .EQ. master) THEN
    WRITE(*,'(A)') '=========================================='
    WRITE(*,'(A)') 'Lubrication Threshold Configuration'
    WRITE(*,'(A)') '=========================================='
    WRITE(*,'(A,ES12.5)') '  Minimum mesh size (h_min):    ', h_min
    WRITE(*,'(A,F6.2)') '  Lubrication factor:           ', dLubricationFactor
    WRITE(*,'(A,ES12.5)') '  Lubrication threshold:        ', lub_threshold
    WRITE(*,'(A)') '=========================================='
  END IF

#endif 
END SUBROUTINE set_lubrication_threshold_from_mesh
!
!----------------------------------------------
!
SUBROUTINE General_init_ext(MDATA,MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
   cProjectFolder,cProjectNumber,nUmbrellaSteps,mg_mesh,nInitUmbrellaSteps
 USE Transport_Q2P1, ONLY : Init_QuadScalar,LinSc,QuadSc
 USE Parametrization, ONLY: InitParametrization,ParametrizeBndr,&
     ProlongateParametrization_STRCT,InitParametrization_STRCT,ParametrizeBndryPoints,&
     DeterminePointParametrization_STRCT,ParametrizeBndryPoints_STRCT
 USE Parametrization, ONLY: ParametrizeQ2Nodes
 USE param_parser, ONLY: GDATNEW
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
 CALL  GDATNEW (CSimPar,0)

 CFILE=CFILE1
 MFILE=MFILE1

 dPeriodicity(1)= 0.1d0
 dPeriodicity(2)= 0.1d0

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

 !!call writeTriFile(mg_mesh%level(NLMAX+1), filename)


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
 if (myid .ne. 0) then
   call commf2c_lubrication_lab(MPI_COMM_WORLD, MPI_Comm_Ex0, myid)
 end if

#ifdef LUBRICATION_PIPELINE
 ! Set lubrication threshold based on mesh size
 call set_lubrication_threshold_from_mesh(myid)
#endif
#endif

 call init_fc_rigid_body(myid)

 call MPI_Barrier(MPI_COMM_WORLD, error_indicator)

RETURN

END SUBROUTINE General_init_ext
 !
 !-----------------------------------------------------------------------
 !
 SUBROUTINE myGDATNEW (cName,iCurrentStatus)
   USE PP3D_MPI
   USE var_QuadScalar, ONLY : myMatrixRenewal,bNonNewtonian,cGridFileName,&
     nSubCoarseMesh,cFBM_File,bTracer,cProjectFile,bMeshAdaptation,&
     myExport,cAdaptedMeshFile,nUmbrellaSteps,bNoOutflow,myDataFile,&
     bViscoElastic,bRefFrame

   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   PARAMETER (NNLEV=9)
   CHARACTER*7 cName
   CHARACTER letter
   INTEGER :: myFile=888
   INTEGER iEnd,iAt,iEq,iLen,iCurrentStatus,istat
   INTEGER iOutShift
   CHARACTER string*500,cVar*7,cPar*25,cLongString*400
   CHARACTER cParam*8,cParam2*20
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

  SUBROUTINE myGDATNEW2 (cName,iCurrentStatus)
   USE PP3D_MPI
   USE var_QuadScalar, ONLY : myMatrixRenewal,bNonNewtonian,cGridFileName,&
     nSubCoarseMesh,cFBM_File,bTracer,cProjectFile,bMeshAdaptation,&
     myExport,cAdaptedMeshFile,nUmbrellaSteps,bNoOutflow,myDataFile,&
     bViscoElastic,bRefFrame, GammaDot, AlphaRelax, RadParticle ! Added GammaDot, AlphaRelax

   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   PARAMETER (NNLEV=9)
   CHARACTER*7 cName
   CHARACTER letter
   INTEGER :: myFile=888
   INTEGER iEnd,iAt,iEq,iLen,iCurrentStatus,istat
   INTEGER iOutShift
   CHARACTER string*500,cVar*7,cPar*25,cLongString*400
   CHARACTER cParam*8,cParam2*20
   LOGICAL bOK,bOutNMAX
   INTEGER, ALLOCATABLE :: iPos(:)

   INTEGER :: iLoc ! Declaration for iLoc was missing

   ! Local loop variable for OutputFields parsing (using a specific name for clarity)
   INTEGER :: idx_loop_output_fields 
   ! Local variable for string length in OutputFields (using a specific name for clarity)
   INTEGER :: len_output_fields_str 
   ! Local variable for number of fields in OutputFields (using specific name for clarity)
   INTEGER :: num_output_fields 


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
   
   ! MFILE is used for writing protocol. It gets its value from MFILE1.
   ! It is implicitly INTEGER due to the IMPLICIT statement.
   INTEGER MFILE 

   SAVE 

   OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action="read",iostat=istat)
   IF (istat .NE. 0) THEN
     WRITE(*,*) 'Could not open data file: ', TRIM(ADJUSTL(myDataFile))
     STOP 'Error opening data file in myGDATNEW'
   END IF

   bOutNMAX = .FALSE.

   DO
     READ (UNIT=myFile,FMT='(A500)',IOSTAT=iEnd) string
     IF (iEnd.EQ.-1) EXIT ! End of file
     IF (iEnd.NE.0) THEN  ! Other read error
        WRITE(*,*) 'Error reading from data file: ', TRIM(ADJUSTL(myDataFile))
        CLOSE(myFile)
        STOP 'Read error in myGDATNEW'
     END IF

     CALL StrStuct() ! Parses 'string' to find '@' and '=', sets iAt, iEq, bOK

     IF (bOK) THEN
       READ(string(1:iAt-1),*,IOSTAT=istat) cVar
       IF (istat .NE. 0) CYCLE ! Skip malformed line (variable name part)

       IF (TRIM(ADJUSTL(cVar)).EQ.TRIM(ADJUSTL(cName))) THEN
         READ(string(iAt+1:iEq-1),*,IOSTAT=istat) cPar
         IF (istat .NE. 0) CYCLE ! Skip malformed line (parameter name part)

         ! The value part of the string is string(iEq+1:)
         ! Most READ statements below will parse from this.
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
           MFILE=MFILE1 ! MFILE (local or common) gets its unit number from MFILE1
         CASE ("StartingProc")
           READ(string(iEq+1:),*) ISTART
         CASE ("Umbrella")
           READ(string(iEq+1:),*) nUmbrellaSteps
         CASE ("StartFile")
           READ(string(iEq+1:),*) CSTART
           ! The commented out section for appending myid is MPI rank specific.
           ! iLen = LEN(TRIM(ADJUSTL(CSTART)))
           ! IF     (myid.lt.10 ) THEN
           !  WRITE(CSTART(iLen+1:),'(A,I1)') "00",myid
           ! ELSEIF (myid.lt.100) THEN
           !  WRITE(CSTART(iLen+1:),'(A,I2)') "0",myid
           ! ELSE 
           !  WRITE(CSTART(iLen+1:),'(I3)') myid
           ! END IF
           MSTART=63
         CASE ("LoadAdaptedMesh")
           bMeshAdaptation = .TRUE. 
           READ(string(iEq+1:),*) cAdaptedMeshFile
         CASE ("SolFile")
           READ(string(iEq+1:),*) CSOL
           iLen = LEN(TRIM(ADJUSTL(CSOL)))
           IF     (myid.LT.10 ) THEN
             WRITE(CSOL(iLen+1:),'(A,I1)') "00",myid
           ELSEIF (myid.LT.100) THEN
             WRITE(CSOL(iLen+1:),'(A,I2)') "0",myid
           ELSE 
             WRITE(CSOL(iLen+1:),'(I3)') myid
           END IF
           MSOL=64
           ISOL = 1
         CASE ("MinMeshLevel")
           READ(string(iEq+1:),*) NLMIN
         CASE ("MaxMeshLevel")
           IF (myid.NE.master) THEN ! master is likely from PP3D_MPI
             READ(string(iEq+1:),*) NLMAX
             myExport%LevelMax = NLMAX
           END IF
           ! IF (myid.eq.MASTER) NLMAX=2 ! Original commented out code
         CASE ("TimeScheme")
           cParam = " " ! Initialize/clear
           READ(string(iEq+1:),*) cParam
           THETA = 0.5d0 ! Default to Crank-Nicolson (CN)
           IF (TRIM(ADJUSTL(cParam)).EQ."BE") THETA = 1.0d0 ! Backward Euler
           IF (TRIM(ADJUSTL(cParam)).EQ."FE") THETA = 0.0d0 ! Forward Euler
         CASE ("TimeStep")
           READ(string(iEq+1:),*) TSTEP
         CASE ("TimeAdaptivity")
           cParam = " "
           READ(string(iEq+1:),*) cParam
           IADTIM = 0 ! Default to No
           IF (TRIM(ADJUSTL(cParam)).EQ."Yes") IADTIM = 1
         CASE ("Tracer")
           cParam = " "
           READ(string(iEq+1:),*) cParam
           bTracer = .FALSE. ! Default to No
           IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bTracer = .TRUE.
         CASE ("ViscoElastic")
           cParam = " "
           READ(string(iEq+1:),*) cParam
           bViscoElastic = .FALSE. ! Default to No
           IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bViscoElastic = .TRUE.
         CASE ("ReferenceFrame")
           cParam = " "
           READ(string(iEq+1:),*) cParam
           bRefFrame = .FALSE. ! Default to No
           IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bRefFrame = .TRUE.
         CASE ("NoOutflow")
           cParam = " "
           READ(string(iEq+1:),*) cParam
           bNoOutflow = .FALSE. ! Default to No
           IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bNoOutflow = .TRUE.
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
           ! Expected format: M<d>D<d>K<d>S<d>C<d> where <d> is a single digit
           READ(string(iEq+1:),*) cParam2 ! e.g., "M1D1K3S0C1"
           iLoc=INDEX (cParam2, 'M', .TRUE.)+1
           READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%M
           iLoc=INDEX (cParam2, 'D', .TRUE.)+1
           READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%D
           iLoc=INDEX (cParam2, 'K', .TRUE.)+1
           READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%K
           iLoc=INDEX (cParam2, 'C', .TRUE.)+1 ! Original code had C after S
           READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%C
           iLoc=INDEX (cParam2, 'S', .TRUE.)+1
           READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%S
         CASE ("FlowType")
           READ(string(iEq+1:),*) cParam2
           bNonNewtonian=.TRUE. ! Default to non-Newtonian
           IF (TRIM(ADJUSTL(cParam2)).EQ."Newtonian") bNonNewtonian=.FALSE.
         CASE ("OutputLevel")
           READ(string(iEq+1:),*) cParam2
           IF (TRIM(ADJUSTL(cParam2)).EQ."MAX") THEN
             bOutNMAX = .TRUE.
             iOutShift = 0
           ELSE IF (TRIM(ADJUSTL(cParam2)).EQ."MAX+1") THEN ! Added ELSE IF for clarity
             bOutNMAX = .TRUE.
             iOutShift = 1
           ELSE IF (TRIM(ADJUSTL(cParam2)).EQ."MAX-1") THEN ! Added ELSE IF for clarity
             bOutNMAX = .TRUE.
             iOutShift = -1
           END IF
           IF (.NOT.bOutNMAX) THEN
             READ(string(iEq+1:),*) myExport%Level ! Assumes value is integer if not MAX,MAX+1,MAX-1
             myExport%Level = MAX(MIN(myExport%Level,myExport%LevelMax+1),1)
           END IF

         CASE ("OutputFormat")
           READ(string(iEq+1:),*) myExport%Format
         CASE ("OutputFields")
           ! Parses a comma-separated list of field names
           READ(string(iEq+1:),*) cLongString
           cLongString = ADJUSTL(TRIM(cLongString))
           len_output_fields_str = LEN(ADJUSTL(TRIM(cLongString))) ! Use specific variable name

           num_output_fields = 0 ! This will count commas
           DO idx_loop_output_fields=1,len_output_fields_str
             READ(cLongString(idx_loop_output_fields:idx_loop_output_fields),'(A)') letter
             IF (letter.EQ.',') num_output_fields = num_output_fields + 1
           END DO
           
           IF (ALLOCATED(myExport%Fields)) DEALLOCATE(myExport%Fields)
           IF (ALLOCATED(iPos)) DEALLOCATE(iPos)
           ! num_output_fields is comma count, so actual fields = num_output_fields + 1
           ! iPos needs N_fields + 1 positions (0, comma_pos_1, ..., comma_pos_N, end_string_pos)
           ALLOCATE(myExport%Fields(num_output_fields+1), iPos(num_output_fields+2))
           
           iPos(1) = 0
           Dim_nfields_idx_tracker = 0 ! Renaming 'nFields' from original when used as tracker
           DO idx_loop_output_fields=1,len_output_fields_str
             READ(cLongString(idx_loop_output_fields:idx_loop_output_fields),'(A)') letter
             IF (letter.EQ.',') THEN
               nfields_idx_tracker = nfields_idx_tracker + 1
               iPos(nfields_idx_tracker+1) = idx_loop_output_fields
             END IF
           END DO
           iPos(nfields_idx_tracker+2) = len_output_fields_str+1
           
           DO idx_loop_output_fields=1,num_output_fields+1 ! Loop for each field string
             READ(cLongString(iPos(idx_loop_output_fields)+1:iPos(idx_loop_output_fields+1)-1),'(A)') &
               myExport%Fields(idx_loop_output_fields)
             myExport%Fields(idx_loop_output_fields) = &
               ADJUSTL(TRIM(myExport%Fields(idx_loop_output_fields)))
           END DO
         
         ! --- Add new parameters here ---
         CASE ("GammaDot")
           READ(string(iEq+1:),*) GammaDot
         CASE ("AlphaRelax")
           READ(string(iEq+1:),*) AlphaRelax
         CASE ("RadParticle")
           READ(string(iEq+1:),*) RadParticle
         
         CASE DEFAULT
           ! Optionally, handle unknown parameters for cName
           ! IF (myid.eq.showid) THEN
           !   WRITE(MTERM,*) "Warning: Unknown parameter '", TRIM(ADJUSTL(cPar)), &
           !                  "' for section '", TRIM(ADJUSTL(cName)), "'"
           ! ENDIF
         END SELECT
       END IF
     END IF
   END DO

   CLOSE (myFile)

   ! Set derived parameters or defaults after reading all data
   M     = 1 ! These seem to be output control flags/units, MTERM is used for terminal.
   MT    = 1 ! M is perhaps also for terminal or general messages.
   IGMV   = NLMAX
   THSTEP=TSTEP*THETA
   IF (bOutNMAX) myExport%Level = NLMAX + iOutShift

   ! Open protocol file for writing if this MPI rank is responsible
   IF (myid.eq.showid) THEN ! showid is likely from PP3D_MPI (e.g. master rank)
     OPEN (UNIT=mfile1,FILE=TRIM(ADJUSTL(cfile1)),action="write",status="replace",iostat=istat)
     IF (istat .NE. 0) THEN
       WRITE(*,*) 'Could not open protocol file for writing: ', TRIM(ADJUSTL(cfile1))
       STOP 'Error opening protocol file in myGDATNEW'
     END IF
   END IF

   ! Print banner if this is the initial call (iCurrentStatus == 0)
   IF (iCurrentStatus.EQ.0) THEN
     IF (myid.eq.showid) WRITE(UNIT=MTERM,FMT=101)
     IF (myid.eq.showid) WRITE(UNIT=MFILE,FMT=101) ! MFILE got unit from MFILE1
   END IF

   ! Printout of all loaded parameters by the 'showid' rank
   IF (myid.eq.showid) THEN
     WRITE(MFILE,'(A,A)') "Meshfolder = ",TRIM(ADJUSTL(cGridFileName))
     WRITE(MTERM,'(A,A)') "Meshfolder = ",TRIM(ADJUSTL(cGridFileName))

     WRITE(MFILE,'(A,I10)') "nSubCoarseMesh = ",nSubCoarseMesh
     WRITE(MTERM,'(A,I10)') "nSubCoarseMesh = ",nSubCoarseMesh

     WRITE(MFILE,'(A,A)') "ParticleFile = ",TRIM(ADJUSTL(cFBM_File))
     WRITE(MTERM,'(A,A)') "ParticleFile = ",TRIM(ADJUSTL(cFBM_File))

     WRITE(MFILE,'(A,A)') "ProjectFile = ",TRIM(ADJUSTL(CProjectFile))
     WRITE(MTERM,'(A,A)') "ProjectFile = ",TRIM(ADJUSTL(CProjectFile))

     WRITE(MFILE,'(A,I10)') "StartingProc = ", ISTART ! Changed format to I10 for consistency
     WRITE(MTERM,'(A,I10)') "StartingProc = ", ISTART

     WRITE(MFILE,'(A,A)') "StartFile = ",TRIM(ADJUSTL(CSTART))
     WRITE(MTERM,'(A,A)') "StartFile = ",TRIM(ADJUSTL(CSTART))

     WRITE(MFILE,'(A,A)') "SolFile = ",TRIM(ADJUSTL(CSOL))
     WRITE(MTERM,'(A,A)') "SolFile = ",TRIM(ADJUSTL(CSOL))

     WRITE(MFILE,'(A,I10)') "MinMeshLevel = ",NLMIN
     WRITE(MTERM,'(A,I10)') "MinMeshLevel = ",NLMIN

     WRITE(MFILE,'(A,I10)') "MaxMeshLevel = ",NLMAX
     WRITE(MTERM,'(A,I10)') "MaxMeshLevel = ",NLMAX

     WRITE(MFILE,'(A,F12.4)') "TimeScheme (Theta) = ",THETA ! Use F format for THETA
     WRITE(MTERM,'(A,F12.4)') "TimeScheme (Theta) = ",THETA

     WRITE(MFILE,'(A,D12.4)') "TimeStep = ",TSTEP
     WRITE(MTERM,'(A,D12.4)') "TimeStep = ",TSTEP

     WRITE(MFILE,'(A,I10)') "TimeAdaptivity (1=Yes,0=No) = ", IADTIM
     WRITE(MTERM,'(A,I10)') "TimeAdaptivity (1=Yes,0=No) = ", IADTIM

     WRITE(MFILE,'(A,D12.4)') "MinTimeAdapt = ",DTMIN
     WRITE(MTERM,'(A,D12.4)') "MinTimeAdapt = ",DTMIN

     WRITE(MFILE,'(A,D12.4)') "MaxTimeAdapt = ",DTMAX
     WRITE(MTERM,'(A,D12.4)') "MaxTimeAdapt = ",DTMAX

     WRITE(MFILE,'(A,D12.4)') "StartSimTime = ", TIMENS
     WRITE(MTERM,'(A,D12.4)') "StartSimTime = ", TIMENS

     WRITE(MFILE,'(A,D12.4)') "MaxSimTime = ", TIMEMX
     WRITE(MTERM,'(A,D12.4)') "MaxSimTime = ", TIMEMX

     WRITE(MFILE,'(A,I10)') "MaxNumStep = ", NITNS
     WRITE(MTERM,'(A,I10)') "MaxNumStep = ", NITNS

     WRITE(MFILE,'(A,I10)') "BackUpFreq = ", INSAV
     WRITE(MTERM,'(A,I10)') "BackUpFreq = ", INSAV

     WRITE(MFILE,'(A,I10)') "BackUpNum = ", INSAVN
     WRITE(MTERM,'(A,I10)') "BackUpNum = ", INSAVN

     WRITE(MFILE,'(A,D12.4)') "OutputFreq = ", DTGMV
     WRITE(MTERM,'(A,D12.4)') "OutputFreq = ", DTGMV

     WRITE(MFILE,'(A,5(A4,I1))') "Matrix Renewal scheme : ",&
       " M=",myMatrixRenewal%M,", D=",myMatrixRenewal%D,&
       ", K=",myMatrixRenewal%K,", S=",myMatrixRenewal%S,", C=",myMatrixRenewal%C
     WRITE(MTERM,'(A,5(A4,I1))') "Matrix Renewal scheme : ",&
       " M=",myMatrixRenewal%M,", D=",myMatrixRenewal%D,&
       ", K=",myMatrixRenewal%K,", S=",myMatrixRenewal%S,", C=",myMatrixRenewal%C

     IF (bMeshAdaptation) THEN 
       WRITE(MFILE,'(A,A)') "Use initial Mesh Adaptation file: ",ADJUSTL(TRIM(cAdaptedMeshFile))
       WRITE(MTERM,'(A,A)') "Use initial Mesh Adaptation file: ",ADJUSTL(TRIM(cAdaptedMeshFile))
     ELSE
       WRITE(MFILE,'(A)') "No Initial Mesh Adaptation"
       WRITE(MTERM,'(A)') "No Initial Mesh Adaptation"
     END IF

     WRITE(MFILE,'(A,I10)') "Number of Umbrella smoothening steps = ",nUmbrellaSteps
     WRITE(MTERM,'(A,I10)') "Number of Umbrella smoothening steps = ",nUmbrellaSteps

     IF (bNoOutflow) THEN 
       WRITE(MFILE,'(A)') "Matrix modification for NoOutflow Condition: Yes"
       WRITE(MTERM,'(A)') "Matrix modification for NoOutflow Condition: Yes"
     ELSE
       WRITE(MFILE,'(A)') "Matrix modification for NoOutflow Condition: No"
       WRITE(MTERM,'(A)') "Matrix modification for NoOutflow Condition: No"
     END IF

     IF (bTracer) THEN 
       WRITE(MFILE,'(A)') "Tracer equation included: Yes"
       WRITE(MTERM,'(A)') "Tracer equation included: Yes"
     ELSE
       WRITE(MFILE,'(A)') "Tracer equation included: No"
       WRITE(MTERM,'(A)') "Tracer equation included: No"
     END IF

     IF (bViscoElastic) THEN 
       WRITE(MFILE,'(A)') "Visco-elastic equation included: Yes"
       WRITE(MTERM,'(A)') "Visco-elastic equation included: Yes"
     ELSE
       WRITE(MFILE,'(A)') "Visco-elastic equation included: No"
       WRITE(MTERM,'(A)') "Visco-elastic equation included: No"
     END IF

     IF (bNonNewtonian) THEN 
       WRITE(MFILE,'(A)') "FlowType = non-Newtonian"
       WRITE(MTERM,'(A)') "FlowType = non-Newtonian"
     ELSE
       WRITE(MFILE,'(A)') "FlowType = Newtonian"
       WRITE(MTERM,'(A)') "FlowType = Newtonian"
     END IF
     
     ! Print new parameters
     WRITE(MFILE,'(A,D12.4)') "GammaDot = ", GammaDot
     WRITE(MTERM,'(A,D12.4)') "GammaDot = ", GammaDot

     WRITE(MFILE,'(A,D12.4)') "AlphaRelax = ", AlphaRelax
     WRITE(MTERM,'(A,D12.4)') "AlphaRelax = ", AlphaRelax
     
     WRITE(MFILE,'(A,D12.4)') "RadParticle = ", RadParticle
     WRITE(MTERM,'(A,D12.4)') "RadParticle = ", RadParticle

     WRITE(MFILE,'(A,A,A,I3,A,A,100A)') "Exporting format: ",TRIM(ADJUSTL(myExport%Format)), &
        ", Level: ", myExport%Level, ", Fields: '", &
        (TRIM("["//TRIM(ADJUSTL(myExport%Fields(idx_loop_output_fields)))//"]"), &
         idx_loop_output_fields=1,num_output_fields+1),"'"
     WRITE(MTERM,'(A,A,A,I3,A,A,100A)') "Exporting format: ",TRIM(ADJUSTL(myExport%Format)), &
        ", Level: ", myExport%Level, ", Fields: '", &
        (TRIM("["//TRIM(ADJUSTL(myExport%Fields(idx_loop_output_fields)))//"]"), &
         idx_loop_output_fields=1,num_output_fields+1),"'"
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
   ! This subroutine parses 'string' (from host) to find positions of '@' and '='.
   ! It sets iAt, iEq, and bOK (variables from host).
   IMPLICIT NONE
   INTEGER :: i, n_len ! Declare local loop counter and length variable

   n_len = LEN(string) ! Use the declared length of string from host
   iAt = 0
   iEq = 0
   DO i=1,n_len
     IF (string(i:i).EQ.'@') iAt = i
     IF (string(i:i).EQ.'=') iEq = i
   END DO

   bOK = .FALSE.
   ! A valid structure requires '@' to appear before '=' and both must be present.
   IF (iAt.GT.0 .AND. iEq.GT.iAt) bOK = .TRUE. 
 END SUBROUTINE StrStuct

 END SUBROUTINE myGDATNEW2
