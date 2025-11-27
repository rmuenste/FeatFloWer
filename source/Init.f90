SUBROUTINE General_init(MDATA,MFILE)
  USE def_FEAT
  USE PP3D_MPI
  USE MESH_Structures
  USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
    cProjectFolder,cProjectNumber,nUmbrellaSteps,mg_mesh
  USE Transport_Q2P1, ONLY : Init_QuadScalar,LinSc,QuadSc
  USE Parametrization, ONLY: InitParametrization,ParametrizeBndr
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
  REAL  :: ttt0,ttt1
  INTEGER II,NDOF
  LOGICAL BLIN
  CHARACTER CFILE*60 !CFILE1*60,
  INTEGER kSubPart,iSubPart,iPart,LenFile
  CHARACTER command*100,CSimPar*7
  CHARACTER (len = 60) :: afile 
  CHARACTER (len = 60) :: bfile 

  INTEGER nLengthV,nLengthE,LevDif
  REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
  logical :: bwait = .true.


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

  !=======================================================================
  !     Grid generation
  !=======================================================================

  IF (myid.eq.1) THEN
    WRITE(MTERM,'(37("-"),A30,37("-"))') " Grid Partititioning starts  "
    WRITE(MFILE,'(37("-"),A30,37("-"))') " Grid Partititioning starts  "
    WRITE(MFILE,*)  "Partitioning command:"
    WRITE(command,'(A,I4,A3,I3,5A)') "./partitioner",subnodes," 1 ",nSubCoarseMesh,&
      ' "',TRIM(ADJUSTL(cGridFileName)),'" "',TRIM(ADJUSTL(cProjectFile)),'" 1> /dev/null'
    WRITE(MTERM,*) "Partitioning command: "
    WRITE(MFILE,*)  "Partitioning command:"
    WRITE(MTERM,*) command
    WRITE(MFILE,*) command
    CALL execute_command_line(command,wait=bwait)
    WRITE(MTERM,'(37("-"),A30,37("-"))') " Grid Partititioning finished "
    WRITE(MFILE,'(37("-"),A30,37("-"))') " Grid Partititioning finished "
  END IF
  CALL CommBarrier()
  CMESH1="_mesh/                 "                     ! PARALLEL
  LenFile = LEN((TRIM(ADJUSTL(cGridFileName))))
  WRITE(CMESH1(7:7+LenFile),'(A,A1)') TRIM(ADJUSTL(cGridFileName)),"/"
  IF (myid.ne.0) THEN                                  ! PARALLEL
    kSubPart = FLOOR(DBLE(subnodes)/DBLE(nSubCoarseMesh)-1d-10)+1
    iSubPart = FLOOR(DBLE(myid)/DBLE(kSubPart)-1d-10)+1
    iPart    = myid - (iSubPart-1)*kSubPart
    WRITE(CMESH1(7+LenFile+1:13+LenFile+1),'(A3,I4.4,A1)') "sub",iSubpart,"/"  ! PARALLEL
    
!     IF     (iSubpart.lt.10 ) THEN
!       WRITE(CMESH1(7+LenFile+1:13+LenFile+1),'(A5,I1,A1)') "sub00",iSubpart,"/"  ! PARALLEL
!     ELSEIF (iSubpart.lt.100) THEN
!       WRITE(CMESH1(7+LenFile+1:13+LenFile+1),'(A4,I2,A1)') "sub0",iSubpart,"/"  ! PARALLEL
!     ELSE
!       WRITE(CMESH1(7+LenFile+1:13+LenFile+1),'(A3,I3,A1)') "sub",iSubpart,"/"  ! PARALLEL
!     END IF

    cProjectFolder = CMESH1

    WRITE(CMESH1(14+LenFile+1:24+LenFile+1),'(A4,I4.4,A4)') "GRID",iPart,".tri"  ! PARALLEL
    WRITE(cProjectNumber(1:3),'(I4.4)') iPart
      
!     IF      (iPart.lt.10) THEN
!       WRITE(CMESH1(14+LenFile+1:24+LenFile+1),'(A6,I1,A4)') "GRID00",iPart,".tri"  ! PARALLEL
!       WRITE(cProjectNumber(1:3),'(A2,I1)') "00",iPart
!     ELSE IF (iPart.lt.100) THEN
!       WRITE(CMESH1(14+LenFile+1:24+LenFile+1),'(A5,I2,A4)') "GRID0",iPart,".tri"  ! PARALLEL
!       WRITE(cProjectNumber(1:3),'(A1,I2)') "0",iPart
!     ELSE
!       WRITE(CMESH1(14+LenFile+1:24+LenFile+1),'(A4,I3,A4)') "GRID",iPart,".tri"  ! PARALLEL
!       WRITE(cProjectNumber(1:3),'(I3)') iPart
!     END IF
  ELSE                                                 ! PARALLEL
    cProjectFolder = CMESH1
    WRITE(CMESH1(7+LenFile+1:14+LenFile+1),'(A8)') "GRID.tri"  ! PARALLEL
  END IF                                               ! PARALLEL

  CALL Init_QuadScalar(mfile)

  IF (myid.EQ.0) NLMAX = LinSc%prm%MGprmIn%MedLev

  if(NLMAX.eq.0)then
    write(*,*)'NLMAX=0 is invalid, exiting...'
  end if

  IF (IER.NE.0) GOTO 99999
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
    mg_Mesh%maxlevel = LinSc%prm%MGprmIn%MedLev
    allocate(mg_mesh%level(LinSc%prm%MGprmIn%MedLev))
  else
    allocate(mg_mesh%level(NLMAX))
    mg_Mesh%maxlevel = nlmax
    mg_Mesh%nlmax = nlmax-1
    mg_Mesh%nlmin = 1
  end if

  call readTriCoarse(CMESH1, mg_mesh)
  IF (myid.EQ.0) then
    call refineMesh(mg_mesh, LinSc%prm%MGprmIn%MedLev)  
  else
    call refineMesh(mg_mesh, NLMAX)  
  end if

  write(*,*)'Refinement finished: ',myid

  II=NLMIN
  IF (myid.eq.1) WRITE(*,*) 'setting up general parallel structures on level : ',II

  CALL PARENTCOMM(mg_mesh%level(II)%nat,&
                  mg_mesh%level(II)%nel,&
                  mg_mesh%level(II)%nvt,&
                  mg_mesh%level(II)%dcorvg,&
                  mg_mesh%level(II)%dcorag,&
                  mg_mesh%level(II)%karea,&
                  mg_mesh%level(II)%kvert)

  IF (myid.EQ.0) NLMAX = 1

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

  IF (myid.NE.0) NLMAX = NLMAX - 1

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

  NDOF = mg_mesh%level(NLMAX)%nvt + mg_mesh%level(NLMAX)%nat + &
         mg_mesh%level(NLMAX)%nel + mg_mesh%level(NLMAX)%net

  CALL E011_CreateComm(NDOF)

  DO ILEV=NLMIN,NLMAX
  if(myid.eq.1)then

  write(*,*)"new:",mg_mesh%level(ILEV)%nvt,&
                   mg_mesh%level(ILEV)%nat,&
                   mg_mesh%level(ILEV)%nel,&
                   mg_mesh%level(ILEV)%net
  
  end if
  end do

  DO ILEV=NLMIN,NLMAX
    !CALL SETLEV(2)
    IF (ILEV.EQ.1) THEN
      CALL InitParametrization(mg_mesh%level(ILEV),ilev)
    ELSE

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
    END IF
  END DO

  ! This part here is responsible for creation of structures enabling the mesh coordinate 
  ! transfer to the master node so that it can create the corresponding matrices
  IF (myid.EQ.0) THEN
    CALL CreateDumpStructures(0)
  ELSE
    LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
    CALL CreateDumpStructures(LevDif)
  END IF
  !       IF (nUmbrellaSteps.GT.0) THEN
  !        CALL UmbrellaSmoother(timens,nUmbrellaSteps)
  !       END IF


  ILEV = LinSc%prm%MGprmIn%MedLev
  !CALL SETLEV(2)

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

  IF (myid.eq.showid) THEN
    WRITE(MTERM,'(10(2XA8))') 'ILEV','NVT','NAT','NEL','NET','NDOF'
    WRITE(MFILE,'(10(2XA8))') 'ILEV','NVT','NAT','NEL','NET','NDOF'
  END IF

  DO II=NLMIN,NLMAX

!  mg_mesh%level(ILEV)%dcorvg,&
!  mg_mesh%level(ILEV)%karea,&
!  mg_mesh%level(ILEV)%kvert,&
!  mg_mesh%level(ILEV)%kedge,&
!  mg_mesh%level(ILEV)%nel,&
!  mg_mesh%level(ILEV)%nvt,&
!  mg_mesh%level(ILEV)%net,&
!  mg_mesh%level(ILEV)%nat

  ILEV=II

!  NVT=KNVT(II)
!  NAT=KNAT(II)
!  NET=KNET(II)
!  NEL=KNEL(II)

  NVT=mg_mesh%level(II)%nvt
  NAT=mg_mesh%level(II)%nat
  NET=mg_mesh%level(II)%net
  NEL=mg_mesh%level(II)%nel

  IF (myid.eq.showid) THEN
    WRITE(MTERM,'(10(2XI8))')ILEV,NVT,NAT,NEL,NET,NVT+NAT+NEL+NET
    WRITE(MFILE,'(10(2XI8))')ILEV,NVT,NAT,NEL,NET,NVT+NAT+NEL+NET
  END IF

  !CALL SETLEV(2)

  IF (myid.eq.showid) THEN
    WRITE(MTERM,'(10(2XI8))')ILEV,NVT,NAT,NEL,NET,NVT+NAT+NEL+NET
    WRITE(MFILE,'(10(2XI8))')ILEV,NVT,NAT,NEL,NET,NVT+NAT+NEL+NET
  END IF

  if(.not.allocated(mg_mesh%level(II)%dvol))then
    allocate(mg_mesh%level(II)%dvol(NEL))
  end if

  CALL  SETARE(mg_mesh%level(II)%dvol,&
               NEL,&
               mg_mesh%level(II)%kvert,&
               mg_mesh%level(II)%dcorvg)

  END DO

  IF (myid.ne.0) THEN
    ILEV=NLMAX +1 

    if(.not.allocated(mg_mesh%level(ILEV)%dvol))then
      allocate(mg_mesh%level(ILEV)%dvol(NEL))
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

  RETURN

  99999 CONTINUE
  END
  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE GDATNEW (cName,iCurrentStatus)
    USE PP3D_MPI
    use iniparser
    USE var_QuadScalar, ONLY : iCommSwitch,BaSynch
    USE var_QuadScalar, ONLY : myMatrixRenewal,bNonNewtonian,cGridFileName,&
      nSubCoarseMesh,cFBM_File,bTracer,cProjectFile,bMeshAdaptation,&
      myExport,cAdaptedMeshFile,nUmbrellaSteps,nInitUmbrellaSteps,bNoOutflow,myDataFile,&
      bViscoElastic,bRefFrame,bSteadyState,Properties,dCGALtoRealFactor,&
      nUmbrellaStepsLvl, nMainUmbrellaSteps,bBoundaryCheck,Transform,postParams,&
      ProlongationDirection,bNS_Stabilization,b2DViscoBench,b3DViscoBench,&
      SSE_HAS_ANGLE, extruder_angle, ApplicationString,VersionString,MaxLevelKnownToMaster, GammaDot, AlphaRelax, RadParticle

    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: NNLEV=9

    ! File unit constants
    INTEGER, PARAMETER :: PARAM_FILE_UNIT = 888
    INTEGER, PARAMETER :: MESH_FILE_UNIT = 61
    INTEGER, PARAMETER :: PROTOCOL_FILE_UNIT = 62
    INTEGER, PARAMETER :: START_FILE_UNIT = 63
    INTEGER, PARAMETER :: SOL_FILE_UNIT = 64

    ! Subroutine arguments
    CHARACTER(len=7) :: cName
    INTEGER :: iCurrentStatus

    ! Local variables - integers
    INTEGER :: myFile = PARAM_FILE_UNIT
    INTEGER :: iEnd, iAt, iEq, iLen, istat
    INTEGER :: iOutShift
    INTEGER :: iVisco, iLoc, iangle, mylength, nFields, i
    INTEGER :: MFILE
    INTEGER, ALLOCATABLE :: iPos(:)

    ! Local variables - characters
    CHARACTER :: letter
    CHARACTER(len=500) :: string
    CHARACTER(len=7) :: cVar
    CHARACTER(len=25) :: cPar
    CHARACTER(len=400) :: cLongString
    CHARACTER(len=8) :: cParam
    CHARACTER(len=20) :: cParam2

    ! Local variables - logicals
    LOGICAL :: bOK, bOutNMAX
    LOGICAL :: is_open

    ! COMMON block variables - /OUTPUT/
    INTEGER :: M, MT, MKEYB, MTERM, MERR, MPROT, MSYS, MTRC, IRECL8

    ! COMMON block variables - /MGPAR/
    INTEGER :: ILEV, NLEV, NLMIN, NLMAX, ICYCLE
    INTEGER :: KPRSM, KPOSM
    ! Note: Array dimensions for KPRSM and KPOSM are specified in COMMON statement

    ! COMMON block variables - /NSPAR/
    DOUBLE PRECISION :: TSTEP, THETA, THSTEP, TIMENS, EPSNS
    INTEGER :: NITNS, ITNS

    ! COMMON block variables - /NSADAT/
    DOUBLE PRECISION :: TIMEMX, DTMIN, DTMAX, DTFACT, TIMEIN
    DOUBLE PRECISION :: EPSADI, EPSADL, EPSADU
    DOUBLE PRECISION :: PRDIF1, PRDIF2
    INTEGER :: IEPSAD, IADIN, IREPIT, IADTIM

    ! COMMON block variables - /NSSAV/
    INTEGER :: INSAV, INSAVN

    ! COMMON block variables - /NSSAVF/
    DOUBLE PRECISION :: DTFILM, DTFILO, DTAVS, DTAVSO, DTGMV, DTGMVO
    INTEGER :: IFUSAV, IFPSAV, IFXSAV, IGID, IGMV, IFINIT

    ! COMMON block variables - /FILES/
    INTEGER :: IMESH1, MMESH1, MFILE1
    INTEGER :: ISTART, MSTART
    INTEGER :: ISOL, MSOL
    CHARACTER(len=60) :: CPARM1, CMESH1, CFILE1, CSTART, CSOL

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

    ! Note: CHARACTER variables for /FILES/ declared above
    COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,&
      ISTART,MSTART,CSTART,ISOL,MSOL,CSOL

    SAVE 

    inquire(unit=myFile, OPENED=is_open)
    if (.not. is_open) then
      OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action="read",iostat=istat)
      if (istat /= 0) then
        write(*,'(A)') 'ERROR: Could not open parameter data file'
        write(*,'(A,A)') '  File: ', TRIM(myDataFile)
        write(*,'(A,I0)') '  iostat = ', istat
        write(*,'(A)') 'Please check that the file exists and is readable.'
        stop 1
      end if
    end if

    bOutNMAX = .false.
    DO
    READ (UNIT=myFile,FMT='(A500)',IOSTAT=iEnd) string
    IF (iEnd == -1) EXIT
    CALL StrStuct()
    IF (bOK) THEN

      READ(string(1:iAt-1),*) cVar

      IF (TRIM(ADJUSTL(cVar)) == TRIM(ADJUSTL(cName))) THEN

        READ(string(iAt+1:iEq-1),*) cPar
        SELECT CASE (TRIM(ADJUSTL(cPar)))

        CASE ("MeshFolder")
          READ(string(iEq+1:),*) cGridFileName
        CASE ("SubMeshNumber")
          READ(string(iEq+1:),*) nSubCoarseMesh
        CASE ("ParticleFile")
          READ(string(iEq+1:),*) cFBM_File
        CASE ("ProjectFile")
          READ(string(iEq+1:),*) cProjectFile
          MMESH1 = MESH_FILE_UNIT
        CASE ("ProtocolFile")
          READ(string(iEq+1:),*) CFILE1
          if (SSE_HAS_ANGLE)then
            CFILE1=''
            iangle = int(extruder_angle)
            write(cfile1,'(a, I4.4,a)') '_data/prot.',iangle,'.txt'
          end if
          MFILE1 = PROTOCOL_FILE_UNIT
          MFILE = MFILE1
        CASE ("StartingProc")
          READ(string(iEq+1:),*) ISTART
        CASE ("Umbrella")
          READ(string(iEq+1:),*) nUmbrellaSteps
        CASE ("InitUmbrella")
          READ(string(iEq+1:),*) nInitUmbrellaSteps
        CASE ("UmbrellaStepM")
          READ(string(iEq+1:),*) nMainUmbrellaSteps
        CASE ("UmbrellaStepL")
          READ(string(iEq+1:),*) nUmbrellaStepsLvl
        CASE ("StartFile")
          READ(string(iEq+1:),*) CSTART
          MSTART = START_FILE_UNIT
        CASE ("LoadAdaptedMesh")
          bMeshAdaptation = .true.
          READ(string(iEq+1:),*) cAdaptedMeshFile
        CASE ("SolFile")
          READ(string(iEq+1:),*) CSOL
          iLen = LEN(TRIM(ADJUSTL(CSOL)))
          IF (myid < 10) THEN
            WRITE(CSOL(iLen+1:),'(A,I1)') "00",myid
          ELSEIF (myid < 100) THEN
            WRITE(CSOL(iLen+1:),'(A,I2)') "0",myid
          ELSE
            WRITE(CSOL(iLen+1:),'(I3)') myid
          END IF
          MSOL = SOL_FILE_UNIT
          ISOL = 1
        CASE ("MinMeshLevel")
          READ(string(iEq+1:),*) NLMIN
        CASE ("MaxMeshLevel")
          IF (myid /= master) THEN
            READ(string(iEq+1:),*) NLMAX
            myExport%LevelMax = NLMAX
            MaxLevelKnownToMaster = NLMAX
          else
            READ(string(iEq+1:),*) MaxLevelKnownToMaster
            myExport%LevelMax = MaxLevelKnownToMaster
          END IF
        CASE ("TimeScheme")
          cParam = " "
          READ(string(iEq+1:),*) cParam
          THETA = 0.5d0
          IF (TRIM(ADJUSTL(cParam)) == "BE") THETA = 1.0d0
          IF (TRIM(ADJUSTL(cParam)) == "FE") THETA = 0.0d0
        CASE ("TimeStep")
          READ(string(iEq+1:),*) TSTEP
        CASE ("Bench_U_mean")
          READ(string(iEq+1:),*) postParams%U_mean
          write(*,*)'umean:',postParams%U_mean
        CASE ("Bench_H")
          READ(string(iEq+1:),*) postParams%H
          write(*,*)'h:',postParams%h
        CASE ("Bench_D")
          READ(string(iEq+1:),*) postParams%D
          write(*,*)'d:',postParams%d
        CASE ("Bench_Sc_U")
          READ(string(iEq+1:),*) postParams%Sc_U
        CASE ("Bench_Sc_Mu")
          READ(string(iEq+1:),*) postParams%Sc_Mu
        CASE ("Bench_Sc_a")
          READ(string(iEq+1:),*) postParams%Sc_a
        CASE ("TimeAdaptivity")
          IADTIM = 0
          IF (read_yes_no_param(string, iEq)) IADTIM = 1
         CASE ("ProlongationDirection")
          READ(string(iEq+1:),*) ProlongationDirection
         CASE ("Tracer")
          bTracer = read_yes_no_param(string, iEq)
        CASE ("ViscoElastic")
          bViscoElastic = read_yes_no_param(string, iEq)
        CASE ("ViscoElasticBench")
          READ(string(iEq+1:),*) iVisco
          b2DViscoBench = .false.
          b3DViscoBench = .false.
          if (iVisco == 2) b2DViscoBench = .true.
          if (iVisco == 3) b3DViscoBench = .true.
        CASE ("SteadyState")
          bSteadyState = read_yes_no_param(string, iEq)
        CASE ("ReferenceFrame")
          bRefFrame = read_yes_no_param(string, iEq)
        CASE ("NoOutflow")
          bNoOutflow = read_yes_no_param(string, iEq)
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
        CASE ("ElemTrans")
          READ(string(iEq+1:),*) Transform%ilint
          If (Transform%ilint.lt.1.or.Transform%ilint.gt.2) Transform%ilint = 2
        CASE ("BackUpFreq")
          READ(string(iEq+1:),*) INSAV
        CASE ("BackUpNum")
          READ(string(iEq+1:),*) INSAVN
        CASE ("BoundaryCheck")
          bBoundaryCheck = read_yes_no_param(string, iEq)
        CASE ("NS_Stabilization")
          bNS_Stabilization = read_yes_no_param(string, iEq)
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
          bNonNewtonian = .true.
          IF (TRIM(ADJUSTL(cParam2)) == "Newtonian") bNonNewtonian = .false.
        CASE ("OutputLevel")
          READ(string(iEq+1:),*) cParam2
          IF (TRIM(ADJUSTL(cParam2)) == "MAX") THEN
            bOutNMAX = .true.
            iOutShift = 0
          END IF
          IF (TRIM(ADJUSTL(cParam2)) == "MAX+1") THEN
            bOutNMAX = .true.
            iOutShift = 1
          END IF
          IF (TRIM(ADJUSTL(cParam2)) == "MAX-1") THEN
            bOutNMAX = .true.
            iOutShift = -1
          END IF
          IF (.not. bOutNMAX) THEN
            READ(string(iEq+1:),*) myExport%Level
            myExport%Level = MAX(MIN(myExport%Level,myExport%LevelMax+1),1)
          END IF

        CASE ("CGALtoRealFactor")
         READ(string(iEq+1:),*)dCGALtoRealFactor

         ! --- Add new parameters here ---
        CASE ("GammaDot")
         READ(string(iEq+1:),*) GammaDot
        CASE ("AlphaRelax")
         READ(string(iEq+1:),*) AlphaRelax
        CASE ("RadParticle")
         READ(string(iEq+1:),*) RadParticle

        CASE ("aSynchComm")
          baSynch = read_yes_no_param(string, iEq)
        CASE ("CommSwitch")
          READ(string(iEq+1:),*) iCommSwitch
        CASE ("OutputFormat")
          READ(string(iEq+1:),*) myExport%Format
        CASE ("OutputFields")
          READ(string(iEq+1:),*) cLongString
          cLongString = ADJUSTL(TRIM(cLongString))
          mylength = LEN(ADJUSTL(TRIM(cLongString)))
          nFields = 0
          DO i=1,mylength
          READ(cLongString(i:i),'(A)') letter
          IF (letter == ',') nFields = nFields + 1
          END DO
          IF (ALLOCATED(myExport%Fields)) DEALLOCATE(myExport%Fields)
          IF (ALLOCATED(iPos)) DEALLOCATE(iPos)
          ALLOCATE(myExport%Fields(nFields+1),iPos(nFields+2))
          iPos(1) = 0
          nFields = 0
          DO i=1,mylength
          READ(cLongString(i:i),'(A)') letter
          IF (letter == ',') THEN
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
    IGMV   = MaxLevelKnownToMaster
    THSTEP=TSTEP*THETA
    IF (bOutNMAX) myExport%Level = MaxLevelKnownToMaster + iOutShift



    inquire(unit=mfile1, OPENED=is_open)
    if (.not. is_open) then
      IF (myid == showid) THEN
        OPEN (UNIT=mfile1,FILE=cfile1,action="write",status="replace",iostat=istat)
        if (istat /= 0) then
          write(*,'(A)') 'ERROR: Could not open protocol file for writing'
          write(*,'(A,A)') '  File: ', TRIM(cfile1)
          write(*,'(A,I0)') '  iostat = ', istat
          write(*,'(A)') 'Please check write permissions and disk space.'
          stop 1
        end if
      end if
    end if

    IF (iCurrentStatus == 0) THEN
      IF (myid == showid) WRITE(UNIT=mterm,FMT=101) ApplicationString,VersionString
      IF (myid == showid) WRITE(UNIT=mfile,FMT=101) ApplicationString,VersionString
    END IF

    ! Printout of all loaded parameters

    IF (myid == showid) THEN
      CALL write_param_int(mfile, mterm, "CommSwitch = ", iCommSwitch)

      IF (baSynch) THEN
       CALL write_param_str(mfile, mterm, "aSynchComm = ", "YES")
      ELSE
       CALL write_param_str(mfile, mterm, "aSynchComm = ", "NO")
      END IF

      CALL write_param_int(mfile, mterm, "nSubCoarseMesh = ", nSubCoarseMesh)
      CALL write_param_str(mfile, mterm, "ParticleFile = ", cFBM_File)
      CALL write_param_str(mfile, mterm, "ProjectFile = ", TRIM(CProjectFile))
      CALL write_param_int(mfile, mterm, "StartingProc = ", ISTART)
      CALL write_param_str(mfile, mterm, "StartFile = ", CSTART)
      CALL write_param_str(mfile, mterm, "SolFile = ", CSOL)
      CALL write_param_int(mfile, mterm, "MinMeshLevel = ", NLMIN)
      CALL write_param_int(mfile, mterm, "MaxMeshLevel = ", NLMAX)
      CALL write_param_real(mfile, mterm, "TimeScheme = ", THETA)
      CALL write_param_real(mfile, mterm, "TimeStep = ", TSTEP)
      CALL write_param_int(mfile, mterm, "TimeAdaptivity = ", IADTIM)
      CALL write_param_real(mfile, mterm, "MinTimeAdapt = ", DTMIN)
      CALL write_param_real(mfile, mterm, "MaxTimeAdapt = ", DTMAX)
      CALL write_param_real(mfile, mterm, "StartSimTime = ", TIMENS)

      CALL write_param_real(mfile, mterm, "MaxSimTime = ", TIMEMX)
      CALL write_param_int(mfile, mterm, "MaxNumStep = ", NITNS)
      CALL write_param_int(mfile, mterm, "BackUpFreq = ", INSAV)
      CALL write_param_int(mfile, mterm, "BackUpNum = ", INSAVN)
      CALL write_param_real(mfile, mterm, "OutputFreq = ", DTGMV)

      WRITE(mfile,'(A,I1)') "ElemTransform = Q", Transform%ilint
      WRITE(mterm,'(A,I1)') "ElemTransform = Q", Transform%ilint
      
      WRITE(mfile,'(A,5(A6I1))') "Matrix Renewal scheme : ","  M = ",myMatrixRenewal%M,&
        ", D = ",myMatrixRenewal%D,", K = ",myMatrixRenewal%K,", S = ",myMatrixRenewal%S,&
        ", C = ",myMatrixRenewal%C
      WRITE(mterm,'(A,5(A6I1))') "Matrix Renewal scheme : ","  M = ",myMatrixRenewal%M,&
        ", D = ",myMatrixRenewal%D,", K = ",myMatrixRenewal%K,", S = ",myMatrixRenewal%S,&
        ", C = ",myMatrixRenewal%C

      IF (bMeshAdaptation) THEN
        CALL write_param_str(mfile, mterm, "Use initial Mesh Adaptation file: ", &
                            ADJUSTL(TRIM(cAdaptedMeshFile)))
      ELSE
        CALL write_param_str(mfile, mterm, "No Initial Mesh Adaptation", "")
      END IF

      IF (bNS_Stabilization) THEN
        CALL write_param_str(mfile, mterm, "Stabilization of Navier-Stokes is ::  ", "ON")
      ELSE
        CALL write_param_str(mfile, mterm, "Stabilization of Navier-Stokes is ::  ", "OFF")
      END IF

      CALL write_param_int(mfile, mterm, "Number of Initial Umbrella smoothening steps", &
                          nInitUmbrellaSteps)
      CALL write_param_int(mfile, mterm, "Number of Umbrella smoothening steps: ", nUmbrellaSteps)
      CALL write_param_int(mfile, mterm, "Number of Umbrella Loops: ", nMainUmbrellaSteps)
      
      WRITE(mfile,'(A,10I10)') "Number of Umbrella steps per levels: ", nUmbrellaStepsLvl
      WRITE(mterm,'(A,10I10)') "Number of Umbrella steps per levels: ", nUmbrellaStepsLvl

      IF (bNoOutflow) THEN
        CALL write_param_str(mfile, mterm, &
          "Matrix modification is to be performed due to the NoOuflow Condition", "")
      END IF

      IF (bTracer) THEN
        CALL write_param_str(mfile, mterm, "Tracer equation is included", "")
      END IF

      IF (bViscoElastic) THEN
        CALL write_param_str(mfile, mterm, "Visco-elastic equation is included", "")
      END IF

      IF (b2DViscoBench .or. b3DViscoBench) THEN
        WRITE(mfile,'(A,I1,A)') "Visco-elastic benchmark computation for ",iVisco,"D"
        WRITE(mterm,'(A,I1,A)') "Visco-elastic benchmark computation for ",iVisco,"D"
      END IF

      IF (bBoundaryCheck) THEN
        CALL write_param_str(mfile, mterm, "BoundaryCheck is ", "ON")
      ELSE
        CALL write_param_str(mfile, mterm, "BoundaryCheck is ", "OFF")
      END IF

      WRITE(mfile,'(A,3ES14.4)') "Newtonian FAC Benchamrk params (U,H,D) : ",postParams%U_mean,postParams%H,postParams%D
      WRITE(mterm,'(A,3ES14.4)') "Newtonian FAC Benchamrk params (U,H,D) : ",postParams%U_mean,postParams%H,postParams%D
      
      WRITE(mfile,'(A,3ES14.4)') "Viscoelastic Benchamrk params (Sc_U,Sc_Mu,Sc_a) : ",postParams%Sc_U,postParams%Sc_Mu,postParams%Sc_a
      WRITE(mterm,'(A,3ES14.4)') "Viscoelastic Benchamrk params (Sc_U,Sc_Mu,Sc_a) : ",postParams%Sc_U,postParams%Sc_Mu,postParams%Sc_a
          
      IF (bNonNewtonian) THEN
        CALL write_param_str(mfile, mterm, "FlowType = ", "non-Newtonian")
      ELSE
        CALL write_param_str(mfile, mterm, "FlowType = ", "Newtonian")
      END IF

      ! Print new parameters
      CALL write_param_real(mfile, mterm, "GammaDot = ", GammaDot)
      CALL write_param_real(mfile, mterm, "AlphaRelax = ", AlphaRelax)
      CALL write_param_real(mfile, mterm, "RadParticle = ", RadParticle)

      
      IF (ProlongationDirection == 0) THEN
       WRITE(mfile,'(A,D12.4)') "Mesh Prolongation is set to  = STANDARD"
       WRITE(mterm,'(A,D12.4)') "Mesh Prolongation is set to  = STANDARD"
      ELSE
       IF (ProlongationDirection == 1) THEN
        WRITE(mfile,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in X axis"
        WRITE(mterm,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in X axis"
       END IF
       IF (ProlongationDirection == 2) THEN
        WRITE(mfile,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in Y axis"
        WRITE(mterm,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in Y axis"
       END IF
       IF (ProlongationDirection == 3) THEN
        WRITE(mfile,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in Z axis"
        WRITE(mterm,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in Z axis"
       END IF
      END IF

      CALL write_param_real(mfile, mterm, "CGALtoRealFactor = ", dCGALtoRealFactor)
      

      WRITE(mfile,'(4A,I3,3A,100A)') "Exporting into ", "'",myExport%Format,"' on level: ",&
        myExport%Level," Fields: '",("["//TRIM(ADJUSTL(myExport%Fields(i)))//"]",i=1,nFields+1),"'"
      WRITE(mterm,'(4A,I3,3A,100A)') "Exporting into ", "'",myExport%Format,"' on level: ",&
        myExport%Level," Fields: '",("["//TRIM(ADJUSTL(myExport%Fields(i)))//"]",i=1,nFields+1),"'"
    END IF


    !-----------------------------------------------------------------------
    101 FORMAT(/2X,100('=')/&
      2X,"|",10X,"                                                                                        |"/&
      2X,"|",10X,"Parellel Q2/P1 FEM Fluid Dynamics code          FeatFloWer                              |"/&
      A104,/&
      A104,/&
      2X,"|",10X,"                                                                                        |"/&
      2X,"|",10X,"Developed by:                                   Otto Mierka, Raphael Münster            |"/&
      2X,"|",10X,"under the supervision of:                       Stefan Turek and Dmitri Kuzmin          |"/&
      2X,"|",10X,"additional contributions from:                  Robert Jendrny, Christoph Lohmann       |"/&
      2X,"|",10X,"                                                Tatiana Theis, Malte Schuh              |"/&
      2X,"|",10X,"                                                Michael Fast                            |"/&
      2X,"|",10X,"                                                                                        |"/&
      2X,"|",10X,"Developed at:                                                                           |"/&
      2X,"|",10X,"                                                                                        |"/&
      2X,"|",10X,"                                                                                        |"/&
      2X,"|",10X,"     #####   #    #         ###     ###   ###   ##### #   #  #    #  #   #  ###         |"/&
      2X,"|",10X,"       #     #    #         #  #   #   #  #  #    #   ## ##  #    #  ##  #  #  #        |"/&
      2X,"|",10X,"       #     #    #   ###   #   #  #   #  ###     #   # # #  #    #  # # #  #   #       |"/&
      2X,"|",10X,"       #     #    #         #  #   #   #  #  #    #   #   #  #    #  #  ##  #  #        |"/&
      2X,"|",10X,"       #      ####          ###     ###   #   #   #   #   #   ####   #   #  ###         |"/&
      2X,"|",10X,"                                                                                        |"/&
      2X,"|",10X,"                                                                                        |"/&
      2X,"|",57X,"            Chair of Mathematics III",5X,"|"/&
      2X,"|",57X,"    Applied Mathematics and Numerics",5X,"|"/&
      2X,"|",57X,"                    Vogelopthsweg 87",5X,"|"/&
      2X,"|",57X,"                      Dortmund 44225",5X,"|"/&
      2X,"|",57X,"                                    ",5X,"|"/&
      2X,"|",10X,"Based on FeatFlow (c)     ",&
      "see also: http://www.featflow.de",30X,"|"/&
      2X,"|",10X,"Correspondance:",73X,"|"/&
      2X,"|",10X,"otto.mierka@math.tu-dortmund.de, ",&
      "stefan.turek@math.tu-dortmund.de",23X,"|"/&
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
      IF (string(i:i) == '@') iAt = i
      IF (string(i:i) == '=') iEq = i
      END DO

      bOk = .false.
      IF (iAt /= 0 .and. iEq /= 0) bOk = .true.

    END SUBROUTINE StrStuct

    !-----------------------------------------------------------------------
    ! Helper subroutine to write parameter to both file and terminal
    !-----------------------------------------------------------------------
    SUBROUTINE write_param_int(mf, mt, label, value)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mf, mt, value
      CHARACTER(*), INTENT(IN) :: label

      WRITE(mf,'(A,I10)') label, value
      WRITE(mt,'(A,I10)') label, value
    END SUBROUTINE write_param_int

    SUBROUTINE write_param_real(mf, mt, label, value)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mf, mt
      DOUBLE PRECISION, INTENT(IN) :: value
      CHARACTER(*), INTENT(IN) :: label

      WRITE(mf,'(A,D12.4)') label, value
      WRITE(mt,'(A,D12.4)') label, value
    END SUBROUTINE write_param_real

    SUBROUTINE write_param_str(mf, mt, label, value)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mf, mt
      CHARACTER(*), INTENT(IN) :: label, value

      WRITE(mf,'(A,A)') label, value
      WRITE(mt,'(A,A)') label, value
    END SUBROUTINE write_param_str

    !-----------------------------------------------------------------------
    ! Helper function to parse Yes/No string to logical value
    !-----------------------------------------------------------------------
    FUNCTION read_yes_no_param(input_string, iEq_pos) RESULT(bool_value)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: input_string
      INTEGER, INTENT(IN) :: iEq_pos
      LOGICAL :: bool_value
      CHARACTER(len=8) :: cParam_local

      cParam_local = " "
      READ(input_string(iEq_pos+1:),*) cParam_local
      bool_value = .false.
      IF (TRIM(ADJUSTL(cParam_local)) == "Yes") bool_value = .true.
    END FUNCTION read_yes_no_param

  END SUBROUTINE GDATNEW
