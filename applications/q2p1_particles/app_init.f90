subroutine init_q2p1_ext(log_unit)
    
  USE def_FEAT
  USE PLinScalar, ONLY : Init_PLinScalar,InitCond_PLinLS, &
    UpdateAuxVariables,Transport_PLinLS,Reinitialize_PLinLS, &
    Reinit_Interphase,dMaxSTF
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar,ProlongateSolution, &
    ResetTimer,bTracer,bViscoElastic,StaticMeshAdaptation
  USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
    Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File

  integer, intent(in) :: log_unit

  !-------INIT PHASE-------

  ! Initialization for FEATFLOW
  CALL General_init_ext(79,log_unit)

  CALL Init_QuadScalar_Stuctures(log_unit)

  IF(bViscoElastic)CALL Init_ViscoScalar_Stuctures(log_unit)

  CALL Init_LinScalar(log_unit)

  CALL InitCond_LinScalar()

  IF (ISTART.EQ.0) THEN
    IF (myid.ne.0) CALL CreateDumpStructures(1)
    CALL InitCond_QuadScalar()
    IF(bViscoElastic)CALL IniProf_ViscoScalar()
  ELSE
    IF (ISTART.EQ.1) THEN
      IF (myid.ne.0) CALL CreateDumpStructures(1)
      CALL SolFromFile(CSTART,1)
    ELSE
      IF (myid.ne.0) CALL CreateDumpStructures(0)
      CALL SolFromFile(CSTART,0)
      CALL ProlongateSolution()
      IF (myid.ne.0) CALL CreateDumpStructures(1)
    END IF
  END IF

end subroutine init_q2p1_ext
!
!----------------------------------------------
!
SUBROUTINE General_init_ext(MDATA,MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
   cProjectFolder,cProjectNumber,nUmbrellaSteps,mg_mesh
 USE Transport_Q2P1, ONLY : Init_QuadScalar,LinSc,QuadSc
 USE Parametrization, ONLY: InitParametrization,ParametrizeBndr
 USE Parametrization, ONLY: ParametrizeQ2Nodes
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
 INTEGER ISE,ISA,ISVEL,ISEEL,ISAEL,ISVED,ISAED,ISVAR
 INTEGER ISEAR,ISEVE,ISAVE,ISVBD,ISEBD,ISABD,IDISP
 INTEGER NEL0,NEL1,NEL2
 REAL ttt0,ttt1
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

 CALL CommBarrier()
 CMESH1="_mesh/                 "                     ! PARALLEL
 LenFile = LEN((TRIM(ADJUSTL(cGridFileName))))
 WRITE(CMESH1(7:7+LenFile),'(A,A1)') TRIM(ADJUSTL(cGridFileName)),"/"
 IF (myid.ne.0) THEN                                  ! PARALLEL
   kSubPart = FLOOR(DBLE(subnodes)/DBLE(nSubCoarseMesh)-1d-10)+1
   iSubPart = FLOOR(DBLE(myid)/DBLE(kSubPart)-1d-10)+1
   iPart    = myid - (iSubPart-1)*kSubPart
   IF     (iSubpart.lt.10 ) THEN
     WRITE(CMESH1(7+LenFile+1:13+LenFile+1),'(A5,I1,A1)') "sub00",iSubpart,"/"  ! PARALLEL
   ELSEIF (iSubpart.lt.100) THEN
     WRITE(CMESH1(7+LenFile+1:13+LenFile+1),'(A4,I2,A1)') "sub0",iSubpart,"/"  ! PARALLEL
   ELSE
     WRITE(CMESH1(7+LenFile+1:13+LenFile+1),'(A3,I3,A1)') "sub",iSubpart,"/"  ! PARALLEL
   END IF

   cProjectFolder = CMESH1

   IF      (iPart.lt.10) THEN
     WRITE(CMESH1(14+LenFile+1:24+LenFile+1),'(A6,I1,A4)') "GRID00",iPart,".tri"  ! PARALLEL
     WRITE(cProjectNumber(1:3),'(A2,I1)') "00",iPart
   ELSE IF (iPart.lt.100) THEN
     WRITE(CMESH1(14+LenFile+1:24+LenFile+1),'(A5,I2,A4)') "GRID0",iPart,".tri"  ! PARALLEL
     WRITE(cProjectNumber(1:3),'(A1,I2)') "0",iPart
   ELSE
     WRITE(CMESH1(14+LenFile+1:24+LenFile+1),'(A4,I3,A4)') "GRID",iPart,".tri"  ! PARALLEL
     WRITE(cProjectNumber(1:3),'(I3)') iPart
   END IF
 ELSE                                                 ! PARALLEL
   cProjectFolder = CMESH1
   WRITE(CMESH1(7+LenFile+1:14+LenFile+1),'(A8)') "GRID.tri"  ! PARALLEL
 END IF                                               ! PARALLEL

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

 ! --------------------------------------------------------------            

#ifdef WITH_ODE 
 call init_fc_rigid_body(myid)      
 call FBM_GetParticles(1)
#else
 call init_fc_rigid_body(myid)      
 call FBM_GetParticles()
#endif
 CALL FBM_ScatterParticles()

 ! --------------------------------------------------------------        

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
     bViscoElastic,bViscoElasticFAC,bRefFrame

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
       CASE ("ViscoElasticFAC")
         cParam = " "
         READ(string(iEq+1:),*) cParam
         bViscoElasticFAC = .FALSE.
         IF (TRIM(ADJUSTL(cParam)).EQ."Yes") bViscoElasticFAC = .TRUE.
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

     IF (bViscoElasticFAC) THEN 
       WRITE(mfile,'(A)') "Visco-elastic flow around cylinder"
       WRITE(mterm,'(A)') "Visco-elastic flow around cylinder"
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
