subroutine init_q2p1_ext(log_unit)
    
  USE def_FEAT
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar,ProlongateSolution, &
    bTracer,bViscoElastic,StaticMeshAdaptation, QuadSc,&
    OperatorRegenaration,UpdateMaterialProperties
  USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
    Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar, &
    LinSc_InitCond_Weber,Boundary_LinSc_Val_Weber
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File, mg_Mesh
  USE Parametrization, ONLY: InitParametrization,ParametrizeBndr,&
      ProlongateParametrization_STRCT,InitParametrization_STRCT,ParametrizeBndryPoints,&
      DeterminePointParametrization_STRCT,ParametrizeBndryPoints_STRCT
  USE app_initialization, only:init_sol_same_level,init_sol_lower_level,init_sol_repart

  integer, intent(in) :: log_unit
  logical            :: I_EXIST

  !-------INIT PHASE-------

  ! Initialization for FEATFLOW
  CALL General_init_ext(79,log_unit)

  INQUIRE (FILE='_data/rheo.s3d', EXIST=I_EXIST)
  if (I_EXIST) then
   CALL ReadS3Dfile('_data/rheo.s3d')
  end if

  CALL Init_QuadScalar_Stuctures(log_unit)

  IF(bViscoElastic)CALL Init_ViscoScalar_Stuctures(log_unit)

  CALL Init_LinScalar(log_unit)

  ! Normal start from inital configuration
  if (istart.eq.0) then
  
    if (myid.ne.0) call CreateDumpStructures(1)
    call InitCond_QuadScalar()

  ! Start from a solution on the same lvl
  ! with the same number of partitions
  elseif (istart.eq.1) then

  call init_sol_same_level(CSTART)

  ! Start from a solution on a lower lvl
  ! with the same number of partitions
  elseif (istart.eq.2)then

    call init_sol_lower_level(CSTART)

  ! Start from a solution on the same lvl
  ! with a different number of partitions
  elseif (istart.eq.3) then

    call init_sol_repart(CSTART)

  end if

  CALL UpdateMaterialProperties()
  
end subroutine init_q2p1_ext
!
!----------------------------------------------
!
SUBROUTINE General_init_ext(MDATA,MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
   cProjectFolder,cProjectNumber,nInitUmbrellaSteps,mg_mesh
 USE Transport_Q2P1, ONLY : Init_QuadScalar,Init_Die_Handlers,LinSc,QuadSc
 USE Parametrization, ONLY: InitParametrization,ParametrizeBndr,&
     ProlongateParametrization_STRCT,InitParametrization_STRCT,ParametrizeBndryPoints,&
     DeterminePointParametrization_STRCT,ParametrizeBndryPoints_STRCT
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
 INTEGER iUmbrella
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
 
 include 'PartitionReader.f90'

 CALL Init_QuadScalar(mfile)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CALL Init_Die_Handlers()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial mesh smoothening !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 DO iUmbrella=1,nInitUmbrellaSteps
  CALL UmbrellaSmoother_STRCT(0d0,1)
!   CALL ProjectPointToSTL(nlmax)
 END DO
 
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
END IF
!!!!!!!!!!!!!!!!!!!!!!!!! FINAL Projection to NLMAX +1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
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

 !CALL SETLEV(2)

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
! -----------------------------------------------------------------------
!
#if !defined WIN32
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE xSEND_START()
USE Sigma_User, ONLY: myProcess
USE PP3D_MPI, ONLY : myid

CHARACTER command*(200),CaseFile*(22)
character(8)  :: cdate
character(10) :: ctime
character(5)  :: czone
integer,dimension(8) :: values

call date_and_time(cdate,ctime,czone,values)
WRITE(CaseFile,'(15A)') 'Case_',cdate(7:8),".",cdate(5:6),".",cdate(3:4)," ",ctime(1:2),":",ctime(3:4),":",ctime(5:6)

IF (myid.eq.1) THEN
 command = " "
 WRITE(command,'(3A)') 'echo `pwd` | cat - _data/rheo.s3d | mail -s "',CaseFile,'_start" omierka@mathematik.uni-dortmund.de'
 CALL system(TRIM(ADJUSTL(command)))
END IF

END SUBROUTINE xSEND_START
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE xSEND_FINISH()
USE Sigma_User, ONLY: myProcess
USE PP3D_MPI, ONLY : myid

CHARACTER command*(200),CaseFile*(22)
character(8)  :: cdate
character(10) :: ctime
character(5)  :: czone
integer,dimension(8) :: values

call date_and_time(cdate,ctime,czone,values)
WRITE(CaseFile,'(15A)') 'Case_',cdate(7:8),".",cdate(5:6),".",cdate(3:4)," ",ctime(1:2),":",ctime(3:4),":",ctime(5:6)

IF (myid.eq.1) THEN
 command = " "
 WRITE(command,'(3A)') 'echo `pwd` | cat - _data/rheo.s3d | mail -s "',CaseFile,'_finish" omierka@mathematik.uni-dortmund.de'
! WRITE(*,'(3A)') "[",TRIM(ADJUSTL(command)),"]"
 CALL system(TRIM(ADJUSTL(command)))
END IF

END SUBROUTINE xSEND_FINISH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
!
! -----------------------------------------------------------------------
!
SUBROUTINE AdjustTimeStepping(dt)
USE PP3D_MPI

IMPLICIT DOUBLE PRECISION(A-H,O-Z)
REAL*8 dt
integer nn0

!-----------------------------------------------------------------------
!     C O M M O N S 
!-----------------------------------------------------------------------
! *** Standard COMMON blocks
! *** COMMON blocks for time discretization
COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,&
                      IFUSAV,IFPSAV,IFXSAV,IGID,IGMV,IFINIT
COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,&
                      EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
SAVE 

nn = NINT(DTGMV/TSTEP)

TSTEP  = dt
DTGMV  = DBLE(nn)*dt
TIMEMX = 1d8 

END SUBROUTINE AdjustTimeStepping
