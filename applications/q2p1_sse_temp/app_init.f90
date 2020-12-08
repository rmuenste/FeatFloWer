subroutine init_q2p1_ext(log_unit)
    
  USE def_FEAT
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar,ProlongateSolution, &
    bTracer,bViscoElastic,StaticMeshAdaptation, QuadSc,&
    OperatorRegenaration
!   USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
!     Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar_XSE, &
    LinSc_InitCond_XSE,Boundary_LinSc_Val_XSE,Tracer
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File, mg_Mesh, myVelo,Screw,Shell,Shearrate,&
      Viscosity,Temperature,myExport,ApplicationString,mySegmentIndicator
!   USE Parametrization, ONLY: InitParametrization,ParametrizeBndr,&
!       ProlongateParametrization_STRCT,InitParametrization_STRCT,ParametrizeBndryPoints,&
!       DeterminePointParametrization_STRCT,ParametrizeBndryPoints_STRCT
  USE app_initialization, only:init_sol_same_level,init_sol_lower_level,init_sol_repart
  USE Sigma_User, only : myProcess, myTransientSolution

  integer, intent(in) :: log_unit
  integer iPeriodicityShift,jFile,iAngle,dump_in_file
  real*8 dTimeStep,dPeriod

  !-------INIT PHASE-------
  ApplicationString = &
"  |                                                          SSE-Temperature module                  |"

  ! Initialization for FEATFLOW
  CALL General_init_ext(79,log_unit)

  CALL Init_QuadScalar_Stuctures(log_unit)

  CALL Init_LinScalar(log_unit)
  
  myExport%cFileName(1)='tmp'
  myExport%cFileName(2)='temp'

!!!!!!!!!!!!!!!!!!!  ---- Velocity Fields are to be loaded -----   !!!!!!!!!!!!!!!!!!!
  !myProcess%nTimeLevels = 36
  !myProcess%Periodicity = 2
  !myTransientSolution%DumpFormat = 2
  !myTransientSolution%nTimeSubStep = 6

  dPeriod = 6d1/myProcess%Umdr
  dTimeStep = dPeriod/DBLE(myProcess%nTimeLevels)
  iAngle = 360/myProcess%nTimeLevels
  iPeriodicityShift = INT(myProcess%nTimeLevels/myProcess%Periodicity)

  IF (.NOT.ALLOCATED(Screw)) ALLOCATE(Screw(QuadSc%ndof))
  IF (.NOT.ALLOCATED(Shell)) ALLOCATE(Shell(QuadSc%ndof))
  IF (.NOT.ALLOCATED(Shearrate)) ALLOCATE(Shearrate(QuadSc%ndof))
  IF (.NOT.ALLOCATED(Viscosity)) ALLOCATE(Viscosity(QuadSc%ndof))
  IF (myProcess%SegmentThermoPhysProps) THEN
   IF (.NOT.ALLOCATED(mySegmentIndicator)) ALLOCATE(mySegmentIndicator(2,QuadSc%ndof))
  END IF
  
  ALLOCATE(myTransientSolution%Velo(3,0:myProcess%nTimeLevels-1))
  ALLOCATE(myTransientSolution%Coor(3,0:myProcess%nTimeLevels-1))
  ALLOCATE(myTransientSolution%Dist(0:myProcess%nTimeLevels-1))
  ALLOCATE(myTransientSolution%Temp(0:myProcess%nTimeLevels-1))
  IF (myProcess%SegmentThermoPhysProps) THEN
   ALLOCATE(myTransientSolution%iSeg(0:myProcess%nTimeLevels-1))
  END IF

  DO iFile=0,myProcess%nTimeLevels/myProcess%Periodicity-1 !myProcess%nTimeLevels
   dump_in_file = iFile*iAngle
   if (myTransientSolution%DumpFormat.eq.2) CALL Load_ListFiles_SSE_temp(dump_in_file)
   ALLOCATE(myTransientSolution%Velo(1,iFile)%x(QuadSc%ndof))
   ALLOCATE(myTransientSolution%Velo(2,iFile)%x(QuadSc%ndof))
   ALLOCATE(myTransientSolution%Velo(3,iFile)%x(QuadSc%ndof))
   ALLOCATE(myTransientSolution%Coor(1,iFile)%x(QuadSc%ndof))
   ALLOCATE(myTransientSolution%Coor(2,iFile)%x(QuadSc%ndof))
   ALLOCATE(myTransientSolution%Coor(3,iFile)%x(QuadSc%ndof))
   ALLOCATE(myTransientSolution%Dist(iFile)%x(QuadSc%ndof))
   IF (myProcess%SegmentThermoPhysProps) THEN
    ALLOCATE(myTransientSolution%iSeg(iFile)%x(QuadSc%ndof))
   END IF
   ALLOCATE(myTransientSolution%Temp(iFile)%x(QuadSc%ndof))
   myTransientSolution%Velo(1,iFile)%x = QuadSc%ValU
   myTransientSolution%Velo(2,iFile)%x = QuadSc%ValV
   myTransientSolution%Velo(3,iFile)%x = QuadSc%ValW
   myTransientSolution%Coor(1,iFile)%x = mg_mesh%level(NLMAX+1)%dcorvg(1,:)
   myTransientSolution%Coor(2,iFile)%x = mg_mesh%level(NLMAX+1)%dcorvg(2,:)
   myTransientSolution%Coor(3,iFile)%x = mg_mesh%level(NLMAX+1)%dcorvg(3,:)
   myTransientSolution%Dist(  iFile)%x = Screw
   myTransientSolution%Temp(  iFile)%x = Temperature
   IF (myProcess%SegmentThermoPhysProps) THEN
    myTransientSolution%iSeg(  iFile)%x = mySegmentIndicator(2,:)
   END IF
   
   IF (myid.eq.1) WRITE(*,*) 'File ',iFile, ' is loaded... ==> angle :', iFile*iAngle
   
   if (myProcess%Periodicity.gt.1) then
    DO iPerio = 1,myProcess%Periodicity-1
     jFile = iFile + iPeriodicityShift*(iPerio) 
     ALLOCATE(myTransientSolution%Velo(1,jFile)%x(QuadSc%ndof))
     ALLOCATE(myTransientSolution%Velo(2,jFile)%x(QuadSc%ndof))
     ALLOCATE(myTransientSolution%Velo(3,jFile)%x(QuadSc%ndof))
     ALLOCATE(myTransientSolution%Coor(1,jFile)%x(QuadSc%ndof))
     ALLOCATE(myTransientSolution%Coor(2,jFile)%x(QuadSc%ndof))
     ALLOCATE(myTransientSolution%Coor(3,jFile)%x(QuadSc%ndof))
     ALLOCATE(myTransientSolution%Dist(jFile)%x(QuadSc%ndof))
     ALLOCATE(myTransientSolution%Temp(jFile)%x(QuadSc%ndof))
     myTransientSolution%Velo(1,jFile)%x = QuadSc%ValU
     myTransientSolution%Velo(2,jFile)%x = QuadSc%ValV
     myTransientSolution%Velo(3,jFile)%x = QuadSc%ValW
     myTransientSolution%Coor(1,jFile)%x = mg_mesh%level(NLMAX+1)%dcorvg(1,:)
     myTransientSolution%Coor(2,jFile)%x = mg_mesh%level(NLMAX+1)%dcorvg(2,:)
     myTransientSolution%Coor(3,jFile)%x = mg_mesh%level(NLMAX+1)%dcorvg(3,:)
     myTransientSolution%Dist(  jFile)%x = Screw
     myTransientSolution%Temp(  jFile)%x = Temperature
     IF (myProcess%SegmentThermoPhysProps) THEN
      myTransientSolution%iSeg(  jFile)%x = mySegmentIndicator(2,:)
     END IF
     IF (myid.eq.1) WRITE(*,*) 'File ',jFile, ' is loaded... ==> angle :', jFile*iAngle
    END DO
   end if
  END DO
  
  call InitCond_LinScalar_XSE(LinSc_InitCond_XSE,Boundary_LinSc_Val_XSE)
  
  !!!! SetThe True Initial Condition for angle 0 !!!!!
  if (myid.ne.0) THEN
   if (istart.ne.0) then
    Tracer%Val(NLMAX+1)%x = myTransientSolution%Temp(0)%x
   else
    if (myid.eq.1) write(*,*) "Setting initial tempeature to :",myProcess%T0
    Tracer%Val(NLMAX+1)%x = myProcess%T0
   end if
  end if
  !!! SetThe True Initial Condition for angle 0 !!!!!
  
end subroutine init_q2p1_ext
!
!----------------------------------------------
!
SUBROUTINE General_init_ext(MDATA,MFILE)
 USE def_FEAT
 USE PP3D_MPI
 USE MESH_Structures
 USE var_QuadScalar, ONLY : cGridFileName,nSubCoarseMesh,cProjectFile,&
   cProjectFolder,cProjectNumber,nInitUmbrellaSteps,mg_mesh,myDataFile
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
 CHARACTER (len = 120) :: cExtrud3DFile
 logical :: bexist = .false.


 CALL ZTIME(TTT0)


 !=======================================================================
 !     Data input
 !=======================================================================
 !
 CALL INIT_MPI()                                 ! PARALLEL
 
 myDataFile='_data/q2p1_paramT.dat'
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

 !------------------------------------------------------------------
 cExtrud3DFile = '_data/Extrud3D_0.dat'
 inquire(file=cExtrud3DFile,Exist=bExist)
 if (bExist) then
  call ReadS3Dfile(cExtrud3DFile)
 else
  write(*,*) 'file: "',adjustl(trim(cExtrud3DFile)),'" does not exist!'
 end if
 !------------------------------------------------------------------

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
!-----------------------------------------------------------------------
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
