subroutine init_q2p1_ext(log_unit)
    
  USE def_FEAT
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Structures_sse, &
    InitCond_QuadScalar,ProlongateSolution,updateFBMGeometry, &
    bTracer,bViscoElastic,StaticMeshAdaptation,&
    LinScalar_InitCond, QuadSc, InitMeshDeform, InitOperators 
    
  USE Transport_Q1, ONLY : Init_GenLinSc_HEATALPHA_Q1
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File,mg_Mesh,tQuadScalar,nUmbrellaStepsLvl,&
      ApplicationString,bMultiMat
  use solution_io, only: read_sol_from_file
  use Sigma_User, only: myProcess
  USE iniparser, ONLY : inip_output_init


  integer, intent(in) :: log_unit

  !-------INIT PHASE-------
  ApplicationString = &
"  |                                                          SSE-FluidDynamics module                |"
  
  ! Initialization for FEATFLOW
  call General_init_ext(79,log_unit)

  call Init_QuadScalar_Structures_sse(log_unit)

  if (bMultiMat) then
   call Init_GenLinSc_HEATALPHA_Q1(log_unit)
  else
   call Init_LinScalar(log_unit)

   call InitCond_LinScalar()
  end if
  
  CALL inip_output_init(myid,showid,log_unit,mterm)
  
  ! Normal start from inital configuration
  if (istart.eq.0) then
    if (myid.ne.0) call CreateDumpStructures(1)
    call InitMeshDeform(log_unit, mg_mesh)
    call InitOperators(log_unit, mg_mesh,.true.)
    call InitCond_QuadScalar()

  ! Start from a solution on the same lvl
  ! with the same number of partitions
  elseif (istart.eq.1) then
    if (myid.ne.0) call CreateDumpStructures(1)
    call LoadMPIDumpFiles(int(myProcess%Angle),'p,v,d,x,t,q')
!    call Load_ListFiles_General(int(myProcess%Angle),'p,v,d,x,t,q')
!     call Load_ListFiles_SSE(int(myProcess%Angle))
!    call read_sol_from_file(CSTART,1,timens)
    if (myid.ne.0) call CreateDumpStructures(1)
    call InitOperators(log_unit, mg_mesh,.true.)

  ! Start from a solution on a lower lvl
  ! with the same number of partitions
  elseif (istart.eq.2)then
    ! In order to read in from a lower level
    ! the lower level structures are needed
    if (myid.ne.0) call CreateDumpStructures(0)
    call read_sol_from_file(CSTART,0,timens)
    call ProlongateSolution()
    ! Now generate the structures for the actual level 
    if (myid.ne.0) call CreateDumpStructures(1)
    call InitOperators(log_unit, mg_mesh,.true.)
    
  ! Start from a solution on the same lvl
  ! with a different number of partitions
  elseif (istart.eq.3) then
    IF (myid.ne.0) CALL CreateDumpStructures(1)
    call SolFromFileRepart(CSTART,1)
    if (myid.ne.0) call CreateDumpStructures(1)
    call InitOperators(log_unit, mg_mesh,.true.)
  end if

 

end subroutine init_q2p1_ext
!
!----------------------------------------------
!
SUBROUTINE General_init_ext(MDATA,MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,mySSE_covergence,&
   cProjectFolder,cProjectNumber,nInitUmbrellaSteps,mg_mesh,MaxLevelKnownToMaster
 USE Transport_Q2P1, ONLY : Init_QuadScalar,LinSc,QuadSc
 USE Parametrization, ONLY: InitParametrization,ParametrizeBndr,&
     ProlongateParametrization_STRCT,InitParametrization_STRCT,ParametrizeBndryPoints,&
     DeterminePointParametrization_STRCT,ParametrizeBndryPoints_STRCT
! USE Parametrization, ONLY: ParametrizeQ2Nodes
 USE Sigma_User, ONLY: mySigma,myProcess,mySetup
 USE cinterface 
 use iniparser

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
 INTEGER iUmbrella
 CHARACTER CFILE*60 !CFILE1*60,
 INTEGER kSubPart,iSubPart,iPart,LenFile
 CHARACTER command*100,CSimPar*7
 CHARACTER (len = 60) :: afile 
 CHARACTER (len = 60) :: bfile 
 CHARACTER (len = 120) :: cExtrud3DFile

 REAL*8 AlphaViscosityMatModel
 REAL*8 dCharVisco,dCharSize,dCharVelo,dCharShear,TimeStep
 CHARACTER sTimeStep*(9)
 
 INTEGER nLengthV,nLengthE,LevDif
 REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
 logical :: bwait = .true.
 logical :: bexist = .false.

 REAL*8 :: dTemp(3),dVisco(3),dShear
 INTEGER i,j

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
 
 IF (MaxLevelKnownToMaster.gt.2) THEN
  QuadSc%prm%MGprmIn%MinLev = max(2,QuadSc%prm%MGprmIn%MinLev)
  QuadSc%prm%MGprmIn%MedLev = max(2,QuadSc%prm%MGprmIn%MedLev)
  if (myid.eq.1) WRITE(MTERM,*) 'Min/Med MG level increasement for velocity solver!'
  if (myid.eq.1) WRITE(MFILE,*) 'Min/Med MG level increasement for velocity solver!'
 END IF

 !------------------------------------------------------------------
 cExtrud3DFile = '_data/Extrud3D.dat'
 inquire(file=cExtrud3DFile,Exist=bExist)
 if (bExist) then
  call ReadS3Dfile(cExtrud3DFile)
  call Setup_STL_Segments()
 else
  write(*,*) 'file: "',adjustl(trim(cExtrud3DFile)),'" does not exist!'
 end if
 !------------------------------------------------------------------

!=====================================================================================
!=====================================================================================
!=====================================================================================
  IF (myid.eq.1) THEN
   WRITE(MTERM,'(A)') 
   WRITE(MTERM,'(A)') " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
   WRITE(MTERM,'(A)') " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
  END IF
  
  !Release the flow-curve to json database for 3 temperatures and a range of shearrates
  if (myid.eq.1) then
   WRITE(mterm,'(A)') '========================== Flow Curve ==================='
   WRITE(mfile,'(A)') '========================== Flow Curve ==================='
   do j=1,3
    dTemp(j)  = myProcess%T0 + dble(j-2)*10.0
   end do
   WRITE(mterm,'(5(A14))') 'Gamma_[1/s]', 'mu_[Pa.s]', 'mu_[Pa.s]', 'mu_[Pa.s]'
   WRITE(mterm,'(A14,4(A4,ES10.3))') '- ', '  T=',dTemp(1), '  T=',dTemp(2), '  T=',dTemp(3)
   WRITE(mfile,'(5(A14))') 'Gamma_[1/s]', 'mu_[Pa.s]', 'mu_[Pa.s]', 'mu_[Pa.s]'
   WRITE(mfile,'(A14,4(A4,ES10.3))') '- ', '  T=',dTemp(1), '  T=',dTemp(2), '  T=',dTemp(3)
   WRITE(mterm,'(A)') '---------------------------------------------------------'
   WRITE(mfile,'(A)') '---------------------------------------------------------'
   do i= -2,5
    dShear = 10**dble(i)
    do j=-1,1
     dTemp(j+2)  = myProcess%T0 + dble(j)*10.0
     dVisco(j+2)     = AlphaViscosityMatModel((dShear**2d0)/2d0,1,dTemp(j+2))
    end do
    WRITE(mterm,'(5(A1,ES13.5))') ' ',dShear,' ',0.1d0*dVisco(1),' ',0.1d0*dVisco(2),' ',0.1d0*dVisco(3)
    WRITE(mfile,'(5(A1,ES13.5))') ' ',dShear,' ',0.1d0*dVisco(1),' ',0.1d0*dVisco(2),' ',0.1d0*dVisco(3)
   end do
   WRITE(mterm,'(A)') '========================================================='
   WRITE(mfile,'(A)') '========================================================='
  end if
  

   
  IF (mySetup%bAutomaticTimeStepControl) THEN
   IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."XSE") THEN
    ! get the characteristic viscosity for characteristic shear rate (10.0[1/s])
     dCharSize      = 0.5d0*(mySigma%Dz_out-mySigma%Dz_in)
     dCharVelo      = 3.14d0*mySigma%Dz_out*(myProcess%Umdr/60d0)
     dCharShear     = dCharVelo/dCharSize
     dCharVisco     = AlphaViscosityMatModel(dCharShear,1,myProcess%T0)
  !    dCharVisco     = ViscosityMatModel(mySetup%CharacteristicShearRate,1,myProcess%T0)
     TimeStep       = 1d-2 * (dCharSize/dCharVisco)
     WRITE(sTimeStep,'(ES9.1)') TimeStep
     READ(sTimeStep,*) TimeStep
     mySSE_covergence%dCharVisco = 0.1d0*dCharVisco

     IF (myid.eq.1) THEN
      WRITE(MTERM,'(A,5ES12.4,A)') " Characteristic size[cm],velo[cm/s],shear[1/s]_E/U: ",dCharSize,dCharVelo,dCharShear,mySetup%CharacteristicShearRate
      WRITE(MFILE,'(A,5ES12.4,A)') " Characteristic size[cm],velo[cm/s],shear[1/s]_E/U: ",dCharSize,dCharVelo,dCharShear,mySetup%CharacteristicShearRate
      WRITE(MTERM,'(A,2ES12.4,A)') " Characteristic viscosity [Pa.s] and corresponding Timestep [s]: ",0.1d0*dCharVisco,TimeStep
      WRITE(MFILE,'(A,2ES12.4,A)') " Characteristic viscosity [Pa.s] and corresponding Timestep [s]: ",0.1d0*dCharVisco,TimeStep
     END IF
   END IF
   
   IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
      ! get the characteristic viscosity for characteristic shear rate (10.0[1/s])
     dCharSize      = 1d-1*myProcess%ExtrusionGapSize
     dCharVelo      = myProcess%ExtrusionSpeed
     dCharShear     = dCharVelo/dCharSize
     dCharVisco     = AlphaViscosityMatModel(dCharShear,1,myProcess%T0)
  !    dCharVisco     = ViscosityMatModel(mySetup%CharacteristicShearRate,1,myProcess%T0)
     TimeStep       = 5d-3 * (dCharSize/dCharVisco)
     WRITE(sTimeStep,'(ES9.1)') TimeStep
     READ(sTimeStep,*) TimeStep

     IF (myid.eq.1) THEN
      WRITE(MTERM,'(A,5ES12.4,A)') " Characteristic size[cm],velo[cm/s],shear[1/s]: ",dCharSize,dCharVelo,dCharShear
      WRITE(MFILE,'(A,5ES12.4,A)') " Characteristic size[cm],velo[cm/s],shear[1/s]: ",dCharSize,dCharVelo,dCharShear
      WRITE(MTERM,'(A,2ES12.4,A)') " Characteristic viscosity [Pa.s] and corresponding Timestep [s]: ",0.1d0*dCharVisco,TimeStep
      WRITE(MFILE,'(A,2ES12.4,A)') " Characteristic viscosity [Pa.s] and corresponding Timestep [s]: ",0.1d0*dCharVisco,TimeStep
     END IF
   END IF
   
   CALL AdjustTimeStepping(TimeStep)
   
 END IF
     
  IF (myid.eq.1) THEN
    WRITE(MTERM,'(A,3ES12.4,I10)') " TSTEP,DTGMV,TIMEMX,NITNS ",TSTEP,DTGMV, TIMEMX, NITNS
    WRITE(MFILE,'(A,3ES12.4,I10)') " TSTEP,DTGMV,TIMEMX,NITNS ",TSTEP,DTGMV, TIMEMX, NITNS
    WRITE(MTERM,'(A)') " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
    WRITE(MTERM,'(A)') " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
    WRITE(MTERM,'(A)') 
  END IF
!=====================================================================================
!=====================================================================================
!=====================================================================================
 
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

!  ILEV = LinSc%prm%MGprmIn%MedLev
! 
!  CALL Create_GlobalNumbering(mg_mesh%level(ILEV)%dcorvg,&
!                              mg_mesh%level(ILEV)%kvert,&
!                              mg_mesh%level(ILEV)%kedge,&
!                              mg_mesh%level(ILEV)%karea,&
!                              mg_mesh%level(ILEV)%nvt,&
!                              mg_mesh%level(ILEV)%net,&
!                              mg_mesh%level(ILEV)%nat,&
!                              mg_mesh%level(ILEV)%nel)

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
!    CALL ProjectPointToSTL(ilev+1)
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
SUBROUTINE Setup_STL_Segments()
USE PP3D_MPI
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,DistTolerance,myOutput
use geometry_processing, only: dEpsDist
implicit none
INTEGER iSeg,iFile,NumberOfSTLDescription

    NumberOfSTLDescription = 0
    DO iSeg=1,mySigma%NumberOfSeg
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL".OR.&
         ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_L".OR.&
         ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_R") THEN
      ALLOCATE(mySigma%mySegment(iSeg)%idxCgal(mySigma%mySegment(iSeg)%nOFFfiles))
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       NumberOfSTLDescription = NumberOfSTLDescription + 1
       mySigma%mySegment(iSeg)%idxCgal(iFile) = NumberOfSTLDescription
      END DO
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_LR") THEN
      ALLOCATE(mySigma%mySegment(iSeg)%idxCgalL(mySigma%mySegment(iSeg)%nOFFfilesL))
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfilesL
       NumberOfSTLDescription = NumberOfSTLDescription + 1
       mySigma%mySegment(iSeg)%idxCgalL(iFile) = NumberOfSTLDescription
      END DO
      ALLOCATE(mySigma%mySegment(iSeg)%idxCgalR(mySigma%mySegment(iSeg)%nOFFfilesR))
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfilesR
       NumberOfSTLDescription = NumberOfSTLDescription + 1
       mySigma%mySegment(iSeg)%idxCgalR(iFile) = NumberOfSTLDescription
      END DO
     END IF
    END DO

    IF (myid.eq.1) OPEN(UNIT=633,FILE='mesh_names.offs')
    IF (myid.eq.1) write(633,'(I0)') NumberOfSTLDescription
    DO iSeg=1,mySigma%NumberOfSeg
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL".OR.&
         ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_L".OR.&
         ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_R") THEN
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       IF (myid.eq.1) write(633,'(A)') adjustl(trim(mySigma%mySegment(iSeg)%OFFfiles(iFile)))
      END DO
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_LR") THEN
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfilesL
       IF (myid.eq.1) write(633,'(A)') adjustl(trim(mySigma%mySegment(iSeg)%OFFfilesL(iFile)))
      END DO
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfilesR
       IF (myid.eq.1) write(633,'(A)') adjustl(trim(mySigma%mySegment(iSeg)%OFFfilesR(iFile)))
      END DO
     END IF
    END DO
    IF (myid.eq.1) CLOSE(633)
    
    dEpsDist = 0.20d0*mySigma%Dz_Out
    
END SUBROUTINE Setup_STL_Segments

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

IF (myProcess%Phase.eq.0.or.myProcess%Phase.eq.3) THEN

call date_and_time(cdate,ctime,czone,values)
WRITE(CaseFile,'(15A)') 'Case_',cdate(7:8),".",cdate(5:6),".",cdate(3:4)," ",ctime(1:2),":",ctime(3:4),":",ctime(5:6)

IF (myid.eq.1) THEN
 command = " "
 WRITE(command,'(3A)') 'echo `pwd` | cat - _data/Extrud3D.dat | mail -s "',CaseFile,'_start" omierka@mathematik.uni-dortmund.de'
 CALL system(TRIM(ADJUSTL(command)))
END IF

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

IF (myProcess%Phase.eq.2.or.myProcess%Phase.eq.3) THEN

call date_and_time(cdate,ctime,czone,values)
WRITE(CaseFile,'(15A)') 'Case_',cdate(7:8),".",cdate(5:6),".",cdate(3:4)," ",ctime(1:2),":",ctime(3:4),":",ctime(5:6)

IF (myid.eq.1) THEN
 command = " "
 WRITE(command,'(3A)') 'echo `pwd` | cat - _data/Extrud3D.dat | mail -s "',CaseFile,'_finish" omierka@mathematik.uni-dortmund.de'
! WRITE(*,'(3A)') "[",TRIM(ADJUSTL(command)),"]"
 CALL system(TRIM(ADJUSTL(command)))
END IF

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
DTGMV  = (DBLE(nn)-0.5d0)*dt
TIMEMX = 1d8 

! TSTEP  = dt
! DTGMV  = DBLE(100)*dt
! TIMEMX = DBLE(NITNS+1)*dt

END SUBROUTINE AdjustTimeStepping
