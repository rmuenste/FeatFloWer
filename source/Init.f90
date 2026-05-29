SUBROUTINE General_init(MDATA,MFILE)
  USE def_FEAT
  USE PP3D_MPI
  USE MESH_Structures
  USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
    cProjectFolder,cProjectNumber,nUmbrellaSteps,mg_mesh
  USE Transport_Q2P1, ONLY : Init_QuadScalar,LinSc,QuadSc
  USE Parametrization, ONLY: InitParametrization,ParametrizeBndr
  USE param_parser, ONLY: GDATNEW
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