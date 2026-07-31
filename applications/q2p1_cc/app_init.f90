subroutine init_q2p1_cc(log_unit)
    
  USE def_FEAT
!   USE PLinScalar, ONLY : Init_PLinScalar,InitCond_PLinLS, &
!     UpdateAuxVariables,Transport_PLinLS,Reinitialize_PLinLS, &
!     Reinit_Interphase,dMaxSTF
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar,ProlongateSolution, &
    bTracer,StaticMeshAdaptation
!   USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
!     Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File
  use Transport_cc
  use transport_linear_cc

  integer, intent(in) :: log_unit

  !-------INIT PHASE-------

  ! Initialization for FEATFLOW
  CALL General_init_cc(79,log_unit)

  call Init_Q2P1_Structures_cc(log_unit)

  call Init_LinScalar_cc()

  call InitCond_LinScalar_cc()

  IF (ISTART.EQ.0) THEN
    IF (myid.ne.0) CALL myCreateDumpStructures(1)
    CALL InitCond_QuadScalar_cc()
!    IF(bViscoElastic)CALL IniProf_ViscoScalar()
  ELSE
    IF (ISTART.EQ.1) THEN
      IF (myid.ne.0) CALL myCreateDumpStructures(1)
      CALL mySolFromFile(CSTART,1)
!      CALL myLoadSmartDumpFiles(CSTART,1)
    ELSE
      IF (myid.ne.0) CALL myCreateDumpStructures(0)
      CALL mySolFromFile(CSTART,0)
!      CALL myLoadSmartDumpFiles(CSTART,0)
      CALL ProlongateSolution()
      IF (myid.ne.0) CALL myCreateDumpStructures(1)
    END IF
  END IF

end subroutine init_q2p1_cc
!
!----------------------------------------------
!
SUBROUTINE General_init_cc(MDATA,MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
   cProjectFolder,cProjectNumber,nUmbrellaSteps,mg_mesh
 USE Transport_Q2P1, ONLY : Init_QuadScalar,LinSc,QuadSc
 USE Parametrization, ONLY: InitParametrization,ParametrizeBndr
 USE Parametrization, ONLY: ParametrizeQ2Nodes
 USE param_parser, ONLY: GDATNEW
 USE Transport_CC, ONLY: Init_CCParam
 USE var_QuadScalar_newton, ONLY: ccParams
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
 
#ifdef HAVE_PE
  include 'PartitionReader.f90'
#else
  include 'PartitionReader2.f90'
#endif

 CALL Init_CCParam(mfile)

 IF (myid.EQ.0) NLMAX = ccParams%MedLev

 if(NLMAX.eq.0)then
   write(*,*)'NLMAX=0 is invalid, exiting...'
   stop
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
   mg_Mesh%nlmax = ccParams%MedLev
   mg_Mesh%nlmin = 1
   mg_Mesh%maxlevel = ccParams%MedLev+1
   allocate(mg_mesh%level(ccParams%MedLev+1))
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

 ! This line of code I guess is wrong if you want to use multi-grid
 ! as the coarse grid solver (so beware):
 IF (myid.EQ.0) NLMAX = ccParams%MedLev
 !     THIS PART WILL BUILD THE REQUIRED COMMUNICATION STRUCTURES
 !     ----------------------------------------------------------

 ! Set up the communication structures for the Quadratic element
 ! Computation minimum level
 ILEV=ccParams%MinLev

 IF (myid.eq.1) write(*,*) 'setting up mapping CC structures for Q2  on level : ',ILEV
 CALL GetQ2CoarseMapper(mg_mesh%level(ILEV)%dcorvg,&
                        mg_mesh%level(ILEV)%kvert,&
                        mg_mesh%level(ILEV)%kedge,&
                        mg_mesh%level(ILEV)%karea,&
                        mg_mesh%level(ILEV)%nvt,&
                        mg_mesh%level(ILEV)%net,&
                        mg_mesh%level(ILEV)%nat,&
                        mg_mesh%level(ILEV)%nel)


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
                             ccParams%MedLev)

 CALL E013_Comm_Master(mg_mesh%level(ILEV)%dcorvg,&
                       mg_mesh%level(ILEV)%kvert,&
                       mg_mesh%level(ILEV)%kedge,&
                       mg_mesh%level(ILEV)%karea,&
                       mg_mesh%level(ILEV)%nvt,&
                       mg_mesh%level(ILEV)%net,&
                       mg_mesh%level(ILEV)%nat,&
                       mg_mesh%level(ILEV)%nel)


 ILEV = ccParams%MedLev

IF (myid.NE.0) NLMAX = NLMAX - 1
 
 ! Structures above the Computation minimum level
 DO ILEV=ccParams%MinLev+1,NLMAX

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
                      ccParams%MedLev)

END DO

 IF (myid.eq.1) write(*,*) 'setting up parallel structures for Q2 :  done!'

 CALL ExtraxtParallelPattern()

 NDOF = mg_mesh%level(NLMAX)%nvt + mg_mesh%level(NLMAX)%nat + &
        mg_mesh%level(NLMAX)%nel + mg_mesh%level(NLMAX)%net

 !    ----------------------------------------------------------            
 call init_fc_rigid_body(myid)      
 call FBM_GetParticles()
 CALL FBM_ScatterParticles()
 !    ----------------------------------------------------------        

 ILEV=NLMIN
 CALL InitParametrization(mg_mesh%level(ILEV),ILEV)
 
 DO ILEV=NLMIN,NLMAX

   CALL ParametrizeBndr(mg_mesh,ilev)

   IF (.not.(myid.eq.0.AND.ilev.gt.ccParams%MedLev)) THEN

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
   CALL myCreateDumpStructures(0)
 ELSE
   LevDif = ccParams%MedLev - NLMAX
   CALL myCreateDumpStructures(LevDif)
 END IF

 ILEV = ccParams%MedLev

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

 ! Set NLMIN to computational NLMIN 
 NLMIN = ccParams%MinLev
 RETURN

END SUBROUTINE General_init_cc
 !
 !-----------------------------------------------------------------------
 !
 SUBROUTINE General_init_ext(MDATA,MFILE)
 END SUBROUTINE General_init_ext
