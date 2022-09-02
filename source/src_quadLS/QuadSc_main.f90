! TODO: new name Navier-Stokes + Q2 or similar
MODULE Transport_Q2P1  

USE def_QuadScalar
! USE PP3D_MPI
USE PP3D_MPI, ONLY:myid,master,E011Sum,COMM_Maximum,COMM_Minimum,&
                   COMM_NLComplete,Comm_Summ,Comm_SummN,&
                   myMPI_Barrier,coarse
USE Parametrization,ONLY : InitBoundaryStructure,ReviseWallBC,myParBndr,&
ParametrizeQ2Nodes

USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,mySetup,myMultiMat,BKTPRELEASE
! USE PP3D_MPI, ONLY:E011Sum,E011True_False,Comm_NLComplete,&
!               Comm_Maximum,Comm_Summ,knprmpi,myid,master
! USE LinScalar, ONLY: AddSurfaceTension
use fbm

use var_QuadScalar, only: QuadSc, LinSc, ViscoSc, PLinSc

use, intrinsic :: ieee_arithmetic

IMPLICIT NONE


REAL*8, ALLOCATABLE :: ST_force(:)
REAL*8 :: Density_Secondary=1d0,Density_Primary=1d0
REAL*8 :: myPowerLawFluid(3),ViscoElasticForce(3)
REAL*8 :: Sigma=0.034D0,DiracEps=0.00625d0
INTEGER, ALLOCATABLE :: QuadScBoundary(:)
INTEGER PressureSample(2)
REAL tttt0,tttt1

!!!!!!!!!!   Artificial - TimeStepControl !!!!!!!!!!!
REAL*8  :: xTimeStepMultiplier=1d0,old_TSTEP
INTEGER :: TimeStepIncrease = 0,TimeStepIncreaseCrit = 3, MaxSmootheningSteps = 32
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! interfaces for the fbm_update and fbm_geom function
! handlers that process the dynamics update and
! the geometric computations for fbm objects
include 'fbm_geom_include.h'
include 'fbm_up_include.h'
include 'fbm_vel_bc_include.h'

! The handler function for the dynamics update
procedure(update_fbm_handler), pointer :: fbm_up_handler_ptr => null()

! The handler function for the geometry update
procedure(fbm_geom_handler), pointer :: fbm_geom_handler_ptr => null()

! The handler function for the velocity boundary condition update
! for the fictitious boundary object
procedure(fbm_velBC_handler), pointer :: fbm_vel_bc_handler_ptr => null()

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE Transport_Q2P1_UxyzP(mfile,inl_u,itns)

use cinterface, only: calculateDynamics,calculateFBM
use fbm, only: fbm_updateFBM
INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit
INTEGER INLComplete,I,J,IERR,iOuter,iITER

CALL updateFBMGeometry()

thstep = tstep*(1d0-theta)

CALL OperatorRegenaration(2)

CALL OperatorRegenaration(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------

! GOTO 1
IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)

 ! Assemble the right hand side
 CALL Matdef_General_QuadScalar(QuadSc,1)

 ! Add the pressure gradient to the rhs
 CALL AddPressureGradient()
END IF

 ! Add the viscoelastic stress to the rhs
 IF(bViscoElastic)THEN
   CALL AddViscoStress()
 END IF

IF (myid.ne.master) THEN
 ! Add the gravity force to the rhs
 CALL AddGravForce()

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

END IF

thstep = tstep*theta

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

CALL COMM_Maximum(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

DO INL=1,QuadSc%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

!!!!          Checking the quality of the result           !!!!
!!!! ----------------------------------------------------- !!!!

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN
! Restore the constant right hand side
 CALL ZTIME(tttt0)
 QuadSc%defU = QuadSc%rhsU
 QuadSc%defV = QuadSc%rhsV
 QuadSc%defW = QuadSc%rhsW
END IF

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

! Checking convergence rates against criterions
RhsUVW=DefUVW
CALL COMM_Maximum(RhsUVW)
CALL Protocol_QuadScalar(mfile,QuadSc,INL,&
     ResU,ResV,ResW,DefUVW,RhsUVW)
IF (ISNAN(RhsUVW)) stop

IF ((DefUVW.LE.DefUVWCrit).AND.&
    (INL.GE.QuadSc%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

IF (INLComplete.eq.1) GOTO 1
!IF (timens.lt.tstep+1d-8) GOTO 1

END DO

1 CONTINUE

! return
myStat%iNonLin = myStat%iNonLin + INL
inl_u = INL

! -------------------------------------------------
! Compute the pressure correction
! -------------------------------------------------
IF (myid.ne.0) THEN

 CALL ZTIME(tttt0)
 ! Save the old solution
 LinSc%valP_old = LinSc%valP(NLMAX)%x
 LinSc%valP(NLMAX)%x = 0d0

 ! Assemble the right hand side (RHS=1/k B^T U~)
 CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)

!  ! Assemble the right hand side (RHS:=RHS-C*Q)
!  CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,2)

 ! Save the right hand side
 LinSc%rhsP(NLMAX)%x = LinSc%defP(NLMAX)%x

 CALL ZTIME(tttt1)
 myStat%tDefP = myStat%tDefP + (tttt1-tttt0)
END IF

! Calling the solver
CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,Boundary_LinScalar_Mat,Boundary_LinScalar_Def,mfile)

CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson equation")

2 CONTINUE

IF (myid.ne.0) THEN
 CALL ZTIME(tttt0)
 CALL Velocity_Correction()
 CALL Pressure_Correction()
 CALL ZTIME(tttt1)
 myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
END IF

CALL QuadScP1toQ2(LinSc,QuadSc)

CALL FAC_GetForces(mfile)

CALL GetNonNewtViscosity()

call fbm_updateFBM(Properties%Density(1),tstep,timens,&
                   Properties%Gravity,mfile,myid,&
                   QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                   LinSc%valP(NLMAX)%x,fbm_up_handler_ptr) 

IF (myid.ne.0) THEN
 CALL STORE_OLD_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
END IF
 
 CALL UmbrellaSmoother_ext(0d0,nUmbrellaSteps)
 
IF (myid.ne.0) THEN
 CALL STORE_NEW_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
END IF
 
 CALL GET_MESH_VELO()
 
 ILEV=NLMAX
 CALL SETLEV(2)
 CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%kedge)

 CALL updateFBMGeometry()

RETURN

END SUBROUTINE Transport_q2p1_UxyzP
!
! ----------------------------------------------
!
! Include custom implementations of the Q2 transport equation
include 'QuadSc_transport_extensions.f90'
!
! ----------------------------------------------
!
SUBROUTINE Init_QuadScalar(mfile)
INTEGER I,J,ndof,mfile

 call Init_Default_Handlers()

 QuadSc%cName = "Velo"
 LinSc%cName = "Pres"

 CALL GetVeloParameters(QuadSc%prm,QuadSc%cName,mfile)
 CALL GetPresParameters(LinSc%prm,LinSc%cName,mfile)

END SUBROUTINE Init_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE Init_Die_Handlers()
! In this function we set the function handlers
! for Die FBM simulations --> as normal but return surface distances as well
implicit none

 fbm_geom_handler_ptr => GetFictKnpr_Die
 fbm_vel_bc_handler_ptr => FictKnpr_velBC

END SUBROUTINE Init_Die_Handlers
!
! ----------------------------------------------
!
SUBROUTINE Init_Wangen_Handlers()
! In this function we set the function handlers
! for Die FBM simulations --> as normal but return surface distances as well
implicit none

 fbm_geom_handler_ptr => GetFictKnpr_Wangen
 fbm_vel_bc_handler_ptr => FictKnpr_velBC_Wangen

END SUBROUTINE Init_Wangen_Handlers
!
! ----------------------------------------------
!
SUBROUTINE Init_Laser_Handlers()
implicit none

 fbm_up_handler_ptr => fbm_updateLaser

END SUBROUTINE Init_Laser_Handlers
!
! ----------------------------------------------
!
SUBROUTINE Init_Default_Handlers()
! In this function we set the function handlers
! for FBM, etc. to their default values
implicit none

 fbm_up_handler_ptr => fbm_updateDefault
 fbm_geom_handler_ptr => fbm_getFictKnpr
 fbm_vel_bc_handler_ptr => fbm_velBC

END SUBROUTINE Init_Default_Handlers
!
! ----------------------------------------------
!
SUBROUTINE Init_QuadScalar_Structures_sse(mfile)
use, intrinsic :: ieee_arithmetic
implicit none
LOGICAL bExist
INTEGER I,J,ndof,mfile,LevDif
integer :: mydof
integer :: maxlevel
Real*8 :: dabl
real*8 :: myInf

 IF (myid.eq.0) then
  IF (LinSc%prm%MGprmIn%CrsSolverType.EQ.7.or.LinSc%prm%MGprmIn%CrsSolverType.EQ.8) THEN
   NLMAX = NLMIN
  END IF
 end if
 
 ILEV=NLMAX
 CALL SETLEV(2)

 ! Initialize the scalar quantity
 CALL InitializeQuadScalar(QuadSc)

 ! Initialize the scalar quantity
 CALL InitializeLinScalar(LinSc)

 ! Initialize the boundary list (QuadScBoundary)
 ALLOCATE (QuadScBoundary(mg_mesh%level(ilev)%nvt+&
                          mg_mesh%level(ilev)%net+&
                          mg_mesh%level(ilev)%nat+&
                          mg_mesh%level(ilev)%nel))

 CALL InitBoundaryList(mg_mesh%level(ILEV)%knpr,&
                       mg_mesh%level(ILEV)%kvert,&
                       mg_mesh%level(ILEV)%kedge,&
                       mg_mesh%level(ILEV)%karea)

 ILEV=NLMAX
 CALL SETLEV(2)

 ! Set up the Coordinate Vector
 ALLOCATE (myQ2Coor(3,mg_mesh%level(ilev)%nvt+&
                      mg_mesh%level(ilev)%net+&
                      mg_mesh%level(ilev)%nat+&
                      mg_mesh%level(ilev)%nel))

 CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%kedge)

 !
 !IF (myid.ne.0) CALL ParametrizeQ2Nodes(myQ2Coor)
 !

 ALLOCATE(myALE%Q2coor_old(3,&
 mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))

 myALE%Q2coor_old = myQ2Coor

 ALLOCATE(myALE%MeshVelo(3,&
 mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))
 myALE%MeshVelo = 0d0
 
! IF (myProcess%SegmentThermoPhysProps) THEN
  allocate(mySegmentIndicator(2,mg_mesh%level(ilev)%nvt+&
  mg_mesh%level(ilev)%net+&
  mg_mesh%level(ilev)%nat+&
  mg_mesh%level(ilev)%nel))
! END IF

 CALL InitBoundaryStructure(mg_mesh%level(ILEV)%kvert,&
                            mg_mesh%level(ILEV)%kedge)

 ILEV=NLMAX
 CALL SETLEV(2)
 CALL ReviseWallBC(mg_mesh,ilev)
!  DO ILEV=NLMIN,NLMAX
!  END DO
 ILEV=NLMAX
 CALL SETLEV(2)
                            
 Properties%cName = "Prop"
 CALL GetPhysiclaParameters(Properties,Properties%cName,mfile)

 if (.not.ALLOCATED(myMultiMat%Mat)) then
  myMultiMat%nOfMaterials = 1
  ALLOCATE(myMultiMat%Mat(myMultiMat%nOfMaterials))
!   myMultiMat%Mat(1)%Rheology%Equation = 5
!   myMultiMat%Mat(1)%Rheology%AtFunc = 1
 end if
 IF (.not.ALLOCATED(MaterialDistribution)) ALLOCATE(MaterialDistribution(1:NLMAX))
 DO ilev=NLMIN,NLMAX
  IF (.not.ALLOCATED(MaterialDistribution(ilev)%x)) ALLOCATE(MaterialDistribution(ilev)%x(mg_mesh%level(ilev)%nel))
  MaterialDistribution(ilev)%x = myMultiMat%initMaterial
 END DO
 
 myPowerLawFluid(2) = 0.001d0
 myPowerLawFluid(3) = 0.75d0

 Properties%Density(1:2) = myThermodyn%Density
 Properties%Gravity(1:3) = 0d0

 ! Initialize the arrays and the distribution of physical properties
 ALLOCATE (mgDensity(NLMIN:NLMAX))
 ALLOCATE (mgNormShearStress(NLMIN:NLMAX))
 DO ILEV=NLMIN,NLMAX

  ALLOCATE (mgDensity(ILEV)%x(mg_mesh%level(ilev)%nel))
  ALLOCATE (mgNormShearStress(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDensity(ILEV)%x          = Properties%Density(1)
  mgNormShearStress(ILEV)%x  = 0d0

 END DO

 if(myid.ne.0)then
!---------------------                          
 ALLOCATE (mgDiffCoeff(NLMIN:NLMAX+1))
 DO ILEV=NLMIN,NLMAX+1
  ALLOCATE (mgDiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDiffCoeff(ILEV)%x = Properties%DiffCoeff(1)
 END DO
!---------------------                          
else
 maxlevel = mg_Mesh%nlmax
 ALLOCATE (mgDiffCoeff(NLMIN:maxlevel))
 DO ILEV=NLMIN,maxlevel
  ALLOCATE (mgDiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDiffCoeff(ILEV)%x = Properties%DiffCoeff(1)
 END DO
end if

 ILEV = NLMAX
 ALLOCATE (Viscosity(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 ALLOCATE (Shearrate(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 ALLOCATE (MaxShearRate(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))
 MaxShearRate = 0d0
                     
 ALLOCATE (Temperature(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 CALL ExtractElemSizeDistribution()
 
 Shearrate = 1d0

 Viscosity = Properties%Viscosity(1)


 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel
         
! Temperature = myProcess%T0
 CALL SetInitialTemperature(Temperature,myQ2Coor,mydof)

 ALLOCATE (myALE%Monitor(mydof))
 ALLOCATE (myALE%NewCoor(3,mydof))
 ALLOCATE (myALE%OldCoor(3,mydof))
! ALLOCATE (myALE%MeshVelo(3,mydof))
 ALLOCATE (myALE%OrigCoor(3,mydof))

 myALE%Monitor   = 1d0
 myALE%MeshVelo  = 0d0

 ! Building up the E013/E013 matrix strucrures
 CALL Create_QuadMatStruct()

 ! Iteration matrix (only allocation)
 CALL Create_AMat() !(A)

 ! Building up the E012/E013 E013/E012 and matrix structures
 CALL Create_QuadLinMatStruct() 

 ! Building up the E012/E012 matrix strucrures
 CALL Create_LinMatStruct ()

 ! Pressure gradient matrix
 CALL Create_BMat() !(B,BT)

 IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "

 IF (myid.ne.master) THEN
  ! Parallel E012/E013 matrix structure
  CALL Create_QuadLinParMatStruct(PLinSc) !(pB)
  
  ! Building up the Parallel E012/E012 matrix strucrures
  CALL Create_ParLinMatStruct ()

  CALL BuildUpPressureCommunicator(LinSc,PLinSc)
END IF

! Set up the boundary condition types (knpr)
 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)
  CALL QuadScalar_Knpr()
 END DO

 ILEV=NLMAX
 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

 ALLOCATE (FictKNPR(mydof))
 FictKNPR=0
 ALLOCATE (Distance(mydof))
 Distance = 0d0

 ALLOCATE (MixerKNPR(mydof))
 MixerKNPR=0
 ALLOCATE (Distamce(mydof))
 Distamce = 0d0

 IF (.NOT.ALLOCATED(Screw)) ALLOCATE(Screw(QuadSc%ndof))
 IF (.NOT.ALLOCATED(ScrewDist)) ALLOCATE(ScrewDist(2,QuadSc%ndof))
 IF (.NOT.ALLOCATED(Shell)) ALLOCATE(Shell(QuadSc%ndof))

 ! SEt up the knpr vector showing dofs with parallel property ...
 IF (myid.ne.0) THEN
  ALLOCATE (ParKNPR(mydof))
  QuadSc%auxU = 1d0
  CALL E013Sum(QuadSc%auxU)
  DO I=1,mydof
   IF (QuadSc%auxU(I).EQ.1d0) THEN
    ParKNPR(I) = 0
   ELSE
    ParKNPR(I) = 1
   END IF
  END DO
 END IF

 if(ieee_support_inf(myInf))then
   myInf = ieee_value(myInf, ieee_negative_inf)
 endif
 
 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."XSE") THEN
  mySSE_covergence%iC = 1
  mySSE_covergence%nC = NINT(DTGMV/TSTEP)
  mySSE_covergence%start = min(NINT(DBLE(NITNS)/2d0),3*mySSE_covergence%nC)
  mySSE_covergence%average = myInf
  mySSE_covergence%std_dev = myInf
  IF (myid.eq.1) WRITE(MTERM,*) "SSE convergence buffer size & start of control:",mySSE_covergence%nC,mySSE_covergence%start
  IF (myid.eq.1) WRITE(MFILE,*) "SSE convergence buffer size & start of control:",mySSE_covergence%nC,mySSE_covergence%start
  allocate(mySSE_covergence%Monitor(mySSE_covergence%nC))
  mySSE_covergence%Monitor = myInf
 ELSE
  mySSE_covergence%start = 10*NITNS
 END IF
 
 IF (myid.eq.showID) THEN
  INQUIRE (FILE="_data/BenchValues.txt", EXIST=bExist)
  IF (ISTART.EQ.0.OR.(.NOT.bExist)) THEN
   OPEN(666,FILE="_data/BenchValues.txt")
   WRITE(666,'(4A16)') "Time","Drag","Lift","ZForce"
  ELSE
   OPEN(666,FILE="_data/BenchValues.txt",ACCESS='APPEND')
  END IF
 END IF

 CALL InitializeProlRest(QuadSc,LinSc)

 CALL OperatorRegenaration(1)
 
 CALL SetUp_HYPRE_Solver(LinSc,PLinSc,mfile)
 
END SUBROUTINE Init_QuadScalar_Structures_sse
!
! ----------------------------------------------
!
SUBROUTINE Init_QuadScalar_Structures_sse_mesh(mfile)
use, intrinsic :: ieee_arithmetic
implicit none
LOGICAL bExist
INTEGER I,J,ndof,mfile,LevDif
integer :: mydof
integer :: maxlevel
Real*8 :: dabl
real*8 :: myInf

 ILEV=NLMAX
 CALL SETLEV(2)

 ! Initialize the scalar quantity
 CALL InitializeQuadScalar(QuadSc)

 ! Initialize the scalar quantity
 CALL InitializeLinScalar(LinSc)

 ILEV=NLMAX
 CALL SETLEV(2)

 ! Set up the Coordinate Vector
 ALLOCATE (myQ2Coor(3,mg_mesh%level(ilev)%nvt+&
                      mg_mesh%level(ilev)%net+&
                      mg_mesh%level(ilev)%nat+&
                      mg_mesh%level(ilev)%nel))

 CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%kedge)
 
 allocate(mySegmentIndicator(2,mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))

 CALL InitBoundaryStructure(mg_mesh%level(ILEV)%kvert,&
                            mg_mesh%level(ILEV)%kedge)

 ILEV=NLMAX
 CALL SETLEV(2)
 CALL ReviseWallBC(mg_mesh,ilev)
 ILEV=NLMAX
 CALL SETLEV(2)
                     
 ALLOCATE (Viscosity(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 ALLOCATE (Shearrate(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))
                     
 ALLOCATE (Temperature(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel
         
 CALL SetInitialTemperature(Temperature,myQ2Coor,mydof)

 ALLOCATE (myALE%Monitor(mydof))
 ALLOCATE (myALE%NewCoor(3,mydof))
 ALLOCATE (myALE%OldCoor(3,mydof))
 ALLOCATE (myALE%OrigCoor(3,mydof))

 myALE%Monitor   = 1d0

 ALLOCATE (FictKNPR(mydof))
 FictKNPR=0
 ALLOCATE (Distance(mydof))
 Distance = 0d0

 ALLOCATE (MixerKNPR(mydof))
 MixerKNPR=0
 ALLOCATE (Distamce(mydof))
 Distamce = 0d0

 IF (.NOT.ALLOCATED(Screw)) ALLOCATE(Screw(QuadSc%ndof))
 IF (.NOT.ALLOCATED(ScrewDist)) ALLOCATE(ScrewDist(2,QuadSc%ndof))
 IF (.NOT.ALLOCATED(Shell)) ALLOCATE(Shell(QuadSc%ndof))
  
END SUBROUTINE Init_QuadScalar_Structures_sse_mesh
!
! ----------------------------------------------
!
SUBROUTINE Init_QuadScalar_Structures_sse_PF(mfile)
implicit none
LOGICAL bExist
INTEGER I,J,ndof,mfile,LevDif
integer :: mydof
integer :: maxlevel
Real*8 :: dabl,X

 ILEV=NLMAX
 CALL SETLEV(2)

 ! Initialize the scalar quantity
 CALL InitializeQuadScalar(QuadSc)

 ! Initialize the scalar quantity
 CALL InitializeLinScalar(LinSc)

 ! Initialize the boundary list (QuadScBoundary)
 ALLOCATE (QuadScBoundary(mg_mesh%level(ilev)%nvt+&
                          mg_mesh%level(ilev)%net+&
                          mg_mesh%level(ilev)%nat+&
                          mg_mesh%level(ilev)%nel))

 CALL InitBoundaryList(mg_mesh%level(ILEV)%knpr,&
                       mg_mesh%level(ILEV)%kvert,&
                       mg_mesh%level(ILEV)%kedge,&
                       mg_mesh%level(ILEV)%karea)

 ILEV=NLMAX
 CALL SETLEV(2)

 ! Set up the Coordinate Vector
 ALLOCATE (myQ2Coor(3,mg_mesh%level(ilev)%nvt+&
                      mg_mesh%level(ilev)%net+&
                      mg_mesh%level(ilev)%nat+&
                      mg_mesh%level(ilev)%nel))

 CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%kedge)

 ALLOCATE(myALE%Q2coor_old(3,&
 mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))

 myALE%Q2coor_old = myQ2Coor

 ALLOCATE(myALE%MeshVelo(3,&
 mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))
 myALE%MeshVelo = 0d0
 
  allocate(mySegmentIndicator(2,mg_mesh%level(ilev)%nvt+&
  mg_mesh%level(ilev)%net+&
  mg_mesh%level(ilev)%nat+&
  mg_mesh%level(ilev)%nel))

 CALL InitBoundaryStructure(mg_mesh%level(ILEV)%kvert,&
                            mg_mesh%level(ILEV)%kedge)

 Properties%cName = "Prop"
 CALL GetPhysiclaParameters(Properties,Properties%cName,mfile)

 if (.not.ALLOCATED(myMultiMat%Mat)) then
  myMultiMat%nOfMaterials = 1
  ALLOCATE(myMultiMat%Mat(myMultiMat%nOfMaterials))
 end if
 IF (.not.ALLOCATED(MaterialDistribution)) ALLOCATE(MaterialDistribution(1:NLMAX))
 DO ilev=NLMIN,NLMAX
  IF (.not.ALLOCATED(MaterialDistribution(ilev)%x)) ALLOCATE(MaterialDistribution(ilev)%x(mg_mesh%level(ilev)%nel))
  MaterialDistribution(ilev)%x = myMultiMat%initMaterial
 END DO
 
 myPowerLawFluid(2) = 0.001d0
 myPowerLawFluid(3) = 0.75d0

 Properties%Density(1:2) = myThermodyn%Density

 ! Initialize the arrays and the distribution of physical properties
 ALLOCATE (mgDensity(NLMIN:NLMAX))
 ALLOCATE (mgNormShearStress(NLMIN:NLMAX))
 DO ILEV=NLMIN,NLMAX

  ALLOCATE (mgDensity(ILEV)%x(mg_mesh%level(ilev)%nel))
  ALLOCATE (mgNormShearStress(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDensity(ILEV)%x          = Properties%Density(1)
  mgNormShearStress(ILEV)%x  = 0d0
  
 END DO

 if(myid.ne.0)then
!---------------------                          
 ALLOCATE (mgDiffCoeff(NLMIN:NLMAX+1))
 DO ILEV=NLMIN,NLMAX+1
  ALLOCATE (mgDiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDiffCoeff(ILEV)%x = Properties%DiffCoeff(1)
 END DO
!---------------------                          
else
 maxlevel = mg_Mesh%nlmax
 ALLOCATE (mgDiffCoeff(NLMIN:maxlevel))
 DO ILEV=NLMIN,maxlevel
  ALLOCATE (mgDiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDiffCoeff(ILEV)%x = Properties%DiffCoeff(1)
 END DO
end if

 ILEV = NLMAX
 ALLOCATE (Viscosity(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 ALLOCATE (Shearrate(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 ALLOCATE (MaxShearRate(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))
 MaxShearRate = 0d0
                     
 ALLOCATE (Temperature(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 CALL ExtractElemSizeDistribution()
 
 Shearrate = 1d0

 Viscosity = Properties%Viscosity(1)


 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel
         
! Temperature = myProcess%T0
 CALL SetInitialTemperature(Temperature,myQ2Coor,mydof)

 ALLOCATE (myALE%Monitor(mydof))
 ALLOCATE (myALE%NewCoor(3,mydof))
 ALLOCATE (myALE%OldCoor(3,mydof))
! ALLOCATE (myALE%MeshVelo(3,mydof))
 ALLOCATE (myALE%OrigCoor(3,mydof))

 myALE%Monitor   = 1d0
 myALE%MeshVelo  = 0d0

 ! Building up the E013/E013 matrix strucrures
 CALL Create_QuadMatStruct()

 ! Iteration matrix (only allocation)
 CALL Create_AMat() !(A)

 ! Building up the E012/E013 E013/E012 and matrix structures
 CALL Create_QuadLinMatStruct() 

 ! Building up the E012/E012 matrix strucrures
 CALL Create_LinMatStruct ()

 ! Pressure gradient matrix
 CALL Create_BMat() !(B,BT)

 IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "

 IF (myid.ne.master) THEN
  ! Parallel E012/E013 matrix structure
  CALL Create_QuadLinParMatStruct(PLinSc) !(pB)
  
  ! Building up the Parallel E012/E012 matrix strucrures
  CALL Create_ParLinMatStruct ()

  CALL BuildUpPressureCommunicator(LinSc,PLinSc)
END IF

! Set up the boundary condition types (knpr)
 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)
  CALL QuadScalar_Knpr()
 END DO

 ILEV=NLMAX
 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

 ALLOCATE (FictKNPR(mydof))
 FictKNPR=0
 ALLOCATE (Distance(mydof))
 Distance = 0d0

 ALLOCATE (MixerKNPR(mydof))
 MixerKNPR=0
 ALLOCATE (Distamce(mydof))
 Distamce = 0d0

 IF (.NOT.ALLOCATED(Screw)) ALLOCATE(Screw(QuadSc%ndof))
 IF (.NOT.ALLOCATED(ScrewDist)) ALLOCATE(ScrewDist(2,QuadSc%ndof))
 IF (.NOT.ALLOCATED(Shell)) ALLOCATE(Shell(QuadSc%ndof))

 ! SEt up the knpr vector showing dofs with parallel property ...
 IF (myid.ne.0) THEN
  ALLOCATE (ParKNPR(mydof))
  QuadSc%auxU = 1d0
  CALL E013Sum(QuadSc%auxU)
  DO I=1,mydof
   IF (QuadSc%auxU(I).EQ.1d0) THEN
    ParKNPR(I) = 0
   ELSE
    ParKNPR(I) = 1
   END IF
  END DO
 END IF

 IF (myid.eq.showID) THEN
  INQUIRE (FILE="_data/BenchValues.txt", EXIST=bExist)
  IF (ISTART.EQ.0.OR.(.NOT.bExist)) THEN
   OPEN(666,FILE="_data/BenchValues.txt")
   WRITE(666,'(4A16)') "Time","Drag","Lift","ZForce"
  ELSE
   OPEN(666,FILE="_data/BenchValues.txt",ACCESS='APPEND')
  END IF
 END IF

 CALL InitializeProlRest(QuadSc,LinSc)

END SUBROUTINE Init_QuadScalar_Structures_sse_PF
!
! ----------------------------------------------
!
SUBROUTINE Init_operators_sse_PF()
implicit none
INTEGER I,J
integer :: maxlevel
Real*8 :: X,daux

if (myid.ne.master) then
!  DO ILEV = NLMIN,mg_Mesh%maxlevel-1
  ILEV=mg_Mesh%maxlevel-1
  do i=1,mg_mesh%level(ilev)%nel
   j = mg_mesh%level(ilev)%nvt+mg_mesh%level(ilev)%net+mg_mesh%level(ilev)%nat+i
   if (GenLinScalar%Fld(2)%val(j).gt.GenLinScalar%Fld(3)%val(j)) THEN
    mgDensity(ILEV)%x(i)  = myMultiMat%Mat(1)%Thermodyn%density
   else
    mgDensity(ILEV)%x(i)  = myMultiMat%Mat(2)%Thermodyn%density
   END IF
  end do
!  end do
 
 DO ILEV = mg_Mesh%maxlevel-2,NLMIN,-1
  DO i=1,mg_mesh%level(ilev)%nel
   DAUX = mgDensity(ILEV+1)%x(i)
   DO j=1,7
    DAUX = DAUX + mgDensity(ILEV+1)%x(mg_mesh%level(ilev)%nel + (i-1)*7 + j)
   END DO
   mgDensity(ILEV)%x(i) = 0.125d0*DAUX
  END DO
 END DO
!  
!  write(*,'(I,100000ES12.4)') myid,mgDensity(NLMIN)%x(:)
 
end if

CALL E010_CollectCoarseVector(mgDensity(NLMIN)%x,mg_mesh%level(NLMIN)%nel)

CALL OperatorRegenaration(1)

END SUBROUTINE Init_operators_sse_PF
!
! ----------------------------------------------
!
SUBROUTINE InitCond_Velocity_PF()

QuadSc%ValU = QuadSc%ValU
QuadSc%ValV = QuadSc%ValV
QuadSc%ValW = QuadSc%ValW + myProcess%umdr/6d1*mySigma%mySegment(1)%t

END SUBROUTINE InitCond_Velocity_PF
!
! ----------------------------------------------
!
SUBROUTINE Init_QuadScalar_Stuctures(mfile)
implicit none
LOGICAL bExist
INTEGER I,J,ndof,mfile,LevDif
integer :: mydof
integer :: maxlevel
Real*8 :: dabl

 bMasterTurnedON = .TRUE.
 IF (myid.eq.0) then
  IF (LinSc%prm%MGprmIn%CrsSolverType.EQ.7.or.LinSc%prm%MGprmIn%CrsSolverType.EQ.8) THEN
   NLMAX = NLMIN
   bMasterTurnedON = .FALSE.
  END IF
 end if
 
 ILEV=NLMAX
 CALL SETLEV(2)

 ! Initialize the scalar quantity
 CALL InitializeQuadScalar(QuadSc)

 ! Initialize the scalar quantity
 CALL InitializeLinScalar(LinSc)

 ! Initialize the boundary list (QuadScBoundary)
 ALLOCATE (QuadScBoundary(mg_mesh%level(ilev)%nvt+&
                          mg_mesh%level(ilev)%net+&
                          mg_mesh%level(ilev)%nat+&
                          mg_mesh%level(ilev)%nel))

 CALL InitBoundaryList(mg_mesh%level(ILEV)%knpr,&
                       mg_mesh%level(ILEV)%kvert,&
                       mg_mesh%level(ILEV)%kedge,&
                       mg_mesh%level(ILEV)%karea)

 ILEV=NLMAX
 CALL SETLEV(2)

 ! Set up the Coordinate Vector
 ALLOCATE (myQ2Coor(3,mg_mesh%level(ilev)%nvt+&
                      mg_mesh%level(ilev)%net+&
                      mg_mesh%level(ilev)%nat+&
                      mg_mesh%level(ilev)%nel))

 CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%kedge)

 !
 !IF (myid.ne.0) CALL ParametrizeQ2Nodes(myQ2Coor)
 !
 ALLOCATE(myALE%Q2coor_old(3,&
 mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))

 myALE%Q2coor_old = myQ2Coor

 ALLOCATE(myALE%MeshVelo(3,&
 mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))
 myALE%MeshVelo = 0d0

 CALL InitBoundaryStructure(mg_mesh%level(ILEV)%kvert,&
                            mg_mesh%level(ILEV)%kedge)

 Properties%cName = "Prop"
 CALL GetPhysiclaParameters(Properties,Properties%cName,mfile)

 if (.not.ALLOCATED(myMultiMat%Mat)) then
  myMultiMat%nOfMaterials = 1
  ALLOCATE(myMultiMat%Mat(myMultiMat%nOfMaterials))
!   myMultiMat%Mat(1)%Rheology%Equation = 5
!   myMultiMat%Mat(1)%Rheology%AtFunc = 1
 end if
 IF (.not.ALLOCATED(MaterialDistribution)) ALLOCATE(MaterialDistribution(1:NLMAX))
 DO ilev=NLMIN,NLMAX
  IF (.not.ALLOCATED(MaterialDistribution(ilev)%x)) ALLOCATE(MaterialDistribution(ilev)%x(mg_mesh%level(ilev)%nel))
  MaterialDistribution(ilev)%x = myMultiMat%initMaterial
 END DO

 myPowerLawFluid(2) = 0.001d0
 myPowerLawFluid(3) = 0.75d0

 !!!!!!!!!!!!!!!!!!!!! otherwise the code takes the density from the q2p1_param.dat file  !!!!!!!!!!!!!!
 IF (TRIM(ADJUSTL(myThermodyn%DensityModel)).ne.'NO') THEN
  IF (myThermodyn%density.gt.0d0) Properties%Density(1) = myThermodyn%density
 END IF
 !!!!!!!!!!!!!!!!!!!!! otherwise the code takes the density from the q2p1_param.dat file  !!!!!!!!!!!!!!

 ! Initialize the arrays and the distribution of physical properties
 ALLOCATE (mgDensity(NLMIN:NLMAX))
 ALLOCATE (mgNormShearStress(NLMIN:NLMAX))
 DO ILEV=NLMIN,NLMAX

  ALLOCATE (mgDensity(ILEV)%x(mg_mesh%level(ilev)%nel))
  ALLOCATE (mgNormShearStress(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDensity(ILEV)%x          = Properties%Density(1)
  mgNormShearStress(ILEV)%x  = 0d0

 END DO

 if(myid.ne.0) then
!---------------------                          
 ALLOCATE (mgDiffCoeff(NLMIN:NLMAX+1))
 DO ILEV=NLMIN,NLMAX+1
  ALLOCATE (mgDiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDiffCoeff(ILEV)%x = Properties%DiffCoeff(1)
 END DO
!---------------------                          
else
 maxlevel = mg_Mesh%nlmax
 ALLOCATE (mgDiffCoeff(NLMIN:maxlevel))
 DO ILEV=NLMIN,maxlevel
  ALLOCATE (mgDiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDiffCoeff(ILEV)%x = Properties%DiffCoeff(1)
 END DO
end if

 ILEV = NLMAX
 ALLOCATE (Viscosity(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 ALLOCATE (Temperature(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))
                     
 CALL ExtractElemSizeDistribution()
 
 CALL ExtractBoundaryNormals(QuadSc)

 Viscosity = Properties%Viscosity(1)
 
 Temperature = myProcess%T0

 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

 ALLOCATE (myALE%Monitor(mydof))
 ALLOCATE (myALE%NewCoor(3,mydof))
 ALLOCATE (myALE%OldCoor(3,mydof))
! ALLOCATE (myALE%MeshVelo(3,mydof))
 ALLOCATE (myALE%OrigCoor(3,mydof))

 myALE%Monitor   = 1d0
 myALE%MeshVelo  = 0d0

 ! Building up the E013/E013 matrix strucrures
 CALL Create_QuadMatStruct()

 ! Iteration matrix (only allocation)
 CALL Create_AMat() !(A)

 ! Building up the E012/E013 E013/E012 and matrix structures
 CALL Create_QuadLinMatStruct() 

 ! Building up the E012/E012 matrix strucrures
 CALL Create_LinMatStruct ()

 ! Pressure gradient matrix
 CALL Create_BMat() !(B,BT)

 IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "

 IF (myid.ne.master) THEN
  ! Parallel E012/E013 matrix structure
  CALL Create_QuadLinParMatStruct(PLinSc) !(pB)
  
  ! Building up the Parallel E012/E012 matrix strucrures
  CALL Create_ParLinMatStruct ()

  CALL BuildUpPressureCommunicator(LinSc,PLinSc)
END IF

 ! Correct the wall BCs
 IF (allocated(mg_mesh%BndryNodes))  then
  ilev = nlmax
  CALL RestrictWallBC()
 END IF
 
! Set up the boundary condition types (knpr)
 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)
  CALL QuadScalar_Knpr()
 END DO

 ILEV=NLMAX
 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

 ALLOCATE (FictKNPR(mydof))
 FictKNPR=0
 ALLOCATE (Distance(mydof))
 Distance = 0d0

 ALLOCATE (MixerKNPR(mydof))
 MixerKNPR=0
 ALLOCATE (Distamce(mydof))
 Distamce = 0d0

 ! SEt up the knpr vector showing dofs with parallel property ...
 IF (myid.ne.0) THEN
  ALLOCATE (ParKNPR(mydof))
  QuadSc%auxU = 1d0
  CALL E013Sum(QuadSc%auxU)
  DO I=1,mydof
   IF (QuadSc%auxU(I).EQ.1d0) THEN
    ParKNPR(I) = 0
   ELSE
    ParKNPR(I) = 1
   END IF
  END DO
 END IF

 IF (myid.eq.showID) THEN
  INQUIRE (FILE="_data/BenchValues.txt", EXIST=bExist)
  IF (ISTART.EQ.0.OR.(.NOT.bExist)) THEN
   OPEN(666,FILE="_data/BenchValues.txt")
   WRITE(666,'(10A16)') "Time","Drag","Lift","ZForce","ForceVx","ForceVy","ForceVz","ForcePx","ForcePy","ForcePz"
  ELSE
   OPEN(666,FILE="_data/BenchValues.txt",ACCESS='APPEND')
  END IF
 END IF
 
 CALL InitializeProlRest(QuadSc,LinSc)

 CALL OperatorRegenaration(1)
 
 CALL Create_MMat()
 
 CALL SetUp_HYPRE_Solver(LinSc,PLinSc,mfile)

END SUBROUTINE Init_QuadScalar_Stuctures
!
! ----------------------------------------------
!
SUBROUTINE InitCond_QuadScalar()
INTEGER i

ILEV=NLMAX
CALL SETLEV(2)

! Set initial conditions
IF (myid.ne.0)then
  IF (myFBM%nParticles.GT.0) THEN
    !    CALL updateFBMGeometry()
  END IF
end if

CALL UmbrellaSmoother_ext(0d0,nUmbrellaSteps)

ILEV=NLMAX
CALL SETLEV(2)
CALL SetUp_myQ2Coor(mg_mesh%level(ilev)%dcorvg,&
  mg_mesh%level(ilev)%dcorag,&
  mg_mesh%level(ilev)%kvert,&
  mg_mesh%level(ilev)%karea,&
  mg_mesh%level(ilev)%kedge)

CALL StoreOrigCoor(mg_mesh%level(NLMAX)%dcorvg)

IF (myid.ne.0) THEN

  CALL QuadScalar_InitCond()

  ! Set dirichlet boundary conditions on the solution
  CALL Boundary_QuadScalar_Val()

  ! Set initial conditions
  CALL LinScalar_InitCond(mg_mesh%level(ilev)%dcorvg,&
    mg_mesh%level(ilev)%kvert)

END IF

END SUBROUTINE InitCond_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE LinScalar_InitCond(dcorvg,kvert)
  REAL*8 dcorvg(3,*)
  INTEGER kvert (8,*)
  REAL*8 PX,PY,PZ,dd
  INTEGER i,j,ivt,iPos

  DO i=1,nel
   PX = 0d0
   PY = 0d0
   PZ = 0d0

   DO j =1,8
   ivt = kvert(j,i)
   PX = PX + 0.125d0*dcorvg(1,ivt)
   PY = PY + 0.125d0*dcorvg(2,ivt)
   PZ = PZ + 0.125d0*dcorvg(3,ivt)
   END DO
   iPos = 4*(i-1)
   CALL GetPresInitVal(PX,PY,PZ,LinSc%valP(NLMAX)%x(iPos+1:iPos+4))
  END DO
  
  LinSc%valP_old(:) = LinSc%valP(NLMAX)%x(:)

END SUBROUTINE LinScalar_InitCond
!
! ----------------------------------------------
!
SUBROUTINE QuadScalar_Knpr()
  INTEGER i,ndof

  ndof  = mg_mesh%level(ilev)%nvt+&
  mg_mesh%level(ilev)%net+&
  mg_mesh%level(ilev)%nat+&
  mg_mesh%level(ilev)%nel

  QuadSc%knprU(ILEV)%x = 0
  QuadSc%knprV(ILEV)%x = 0
  QuadSc%knprW(ILEV)%x = 0

  DO i=1,ndof

  IF (myBoundary%bWall(i).OR.myBoundary%iInflow(i).NE.0) THEN
    QuadSc%knprU(ILEV)%x(i) = 1
    QuadSc%knprV(ILEV)%x(i) = 1
    QuadSc%knprW(ILEV)%x(i) = 1
  END IF

  IF (myBoundary%bSymmetry(1,i)) THEN
    QuadSc%knprU(ILEV)%x(i) = 1
    !  WRITE(*,*) myid,"Symmetry u"
  END IF
  IF (myBoundary%bSymmetry(2,i)) THEN
    QuadSc%knprV(ILEV)%x(i) = 1
    !  WRITE(*,*) myid,"Symmetry v"
  END IF
  IF (myBoundary%bSymmetry(3,i)) THEN
    QuadSc%knprW(ILEV)%x(i) = 1
    !  WRITE(*,*) myid,"Symmetry w"
  END IF

  END DO

  IF (myid.eq.0) THEN
    ! bNoOutFlow = .TRUE.
    ! bNoOutFlow = .FALSE.
    ! DO i=1,nvt+net+nat+nel
    !  iBndr = QuadScBoundary(i)
    !  inprU  = QuadSc%knprU(ILEV)%x(i)
    !  inprV  = QuadSc%knprV(ILEV)%x(i)
    !  inprW  = QuadSc%knprW(ILEV)%x(i)
    !  IF (iBndr.EQ.1.AND.inprU.EQ.0) bNoOutFlow = .FALSE.
    ! END DO
  END IF

END SUBROUTINE QuadScalar_Knpr
!
! ----------------------------------------------
!
SUBROUTINE QuadScalar_FictKnpr(dcorvg,dcorag,kvert,kedge,karea)
  use fbm, only: fbm_updateFBMGeom
  REAL*8  dcorvg(3,*),dcorag(3,*)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
  REAL*8 PX,PY,PZ,DIST
  REAL tttx0,tttx1
  INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

  if (myid.eq.0) return
  
  CALL myMPI_Barrier()
  call ztime(tttx0)

  DO i=1,nvt
  PX = dcorvg(1,I)
  PY = dcorvg(2,I)
  PZ = dcorvg(3,I)
  call fbm_updateFBMGeom(PX,PY,PZ,QuadScBoundary(i),FictKNPR(i),Distance(i),fbm_geom_handler_ptr)
  END DO

  k=1
  DO i=1,nel
  DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
    ivt1 = kvert(NeighE(1,j),i)
    ivt2 = kvert(NeighE(2,j),i)
    PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
    PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
    PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
    call fbm_updateFBMGeom(PX,PY,PZ,QuadScBoundary(nvt+k),FictKNPR(nvt+k),Distance(nvt+k),fbm_geom_handler_ptr)
    k = k + 1
  END IF
  END DO
  END DO

  k=1
  DO i=1,nel
  DO j=1,6
  IF (k.eq.karea(j,i)) THEN
    ivt1 = kvert(NeighA(1,j),i)
    ivt2 = kvert(NeighA(2,j),i)
    ivt3 = kvert(NeighA(3,j),i)
    ivt4 = kvert(NeighA(4,j),i)
    PX = 0.25d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2)+dcorvg(1,ivt3)+dcorvg(1,ivt4))
    PY = 0.25d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2)+dcorvg(2,ivt3)+dcorvg(2,ivt4))
    PZ = 0.25d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2)+dcorvg(3,ivt3)+dcorvg(3,ivt4))
    call fbm_updateFBMGeom(PX,PY,PZ,QuadScBoundary(nvt+net+k),FictKNPR(nvt+net+k),Distance(nvt+net+k),fbm_geom_handler_ptr)
    k = k + 1
  END IF
  END DO
  END DO

  ! DO i=1,nat
  !  PX = dcorag(1,I)
  !  PY = dcorag(2,I)
  !  PZ = dcorag(3,I)
  !  CALL GetFictKnpr(PX,PY,PZ,QuadScBoundary(nvt+net+i),FictKNPR(nvt+net+i),Distance(nvt+net+i))
  ! END DO

  DO i=1,nel
  PX = 0d0
  PY = 0d0
  PZ = 0d0
  DO j=1,8
  PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
  PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
  PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
  END DO
  call fbm_updateFBMGeom(PX,PY,PZ,QuadScBoundary(nvt+net+nat+i),FictKNPR(nvt+net+nat+i),Distance(nvt+net+nat+i),fbm_geom_handler_ptr)
  END DO

  CALL myMPI_Barrier()
  call ztime(tttx1)
  if (myid.eq.1) WRITE(*,*) 'FBM time : ', tttx1-tttt0, ' s'

  do i=1,nvt+net+nat+nel
  myALE%Monitor(i)=distance(i)
  end do


END SUBROUTINE QuadScalar_FictKnpr
!
! ----------------------------------------------
!
SUBROUTINE QuadScalar_FictKnpr_Wangen(dcorvg,dcorag,kvert,kedge,karea)
  use fbm, only: fbm_updateFBMGeom,myFBM
  REAL*8  dcorvg(3,*),dcorag(3,*)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
  REAL*8 PX,PY,PZ,DIST,P(3),dR,distX,distY
  REAL tttx0,tttx1
  INTEGER i,j,k,ivt,ndof,IP,minpr,IPP
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
  LOGICAL, allocatable :: bMark(:)
  REAL*8 :: iCount(3)
  
  iCount = 0d0
  ndof = nvt+nat+net+nel
  
  if (myid.eq.0) goto 1
  
  CALL myMPI_Barrier()
  call ztime(tttx0)
  
  ALLOCATE(bMark(ndof))
  bMark = .false.
  Distance = 1d8
  
  DO IPP = 1,myFBM%nParticles
  IP = IPP
  
  DO i=1,nel
   PX = 0d0
   PY = 0d0
   PZ = 0d0
   DO j=1,8
    PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
    PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
    PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
   END DO
   call fbm_updateFBMGeom(PX,PY,PZ,IP,minpr,distX,fbm_geom_handler_ptr)
   
   if (distX.lt.Distance(nvt+net+nat+i)) then
    FictKNPR(nvt+net+nat+i) = minpr
    Distance(nvt+net+nat+i) = distX
    bMark(nvt+net+nat+i) = .true.
   end if
   iCount(1) = iCount(1) + 1d0
   
   dR = 0d0
   do j = 1,8
    ivt = kvert(j,i)
    P = dcorvg(:,ivt)
    dR = max(dR, SQRT((P(1)-PX)**2d0 + (P(2)-PY)**2d0 + (P(3)-PZ)**2d0))
   end do
   
   IF (1.4d0*dR.gt.distX) then
    ! nvt
    do j = 1,8
     ivt = kvert(j,i)
     P = dcorvg(:,ivt)
     if (.not.bMark(ivt).or.DistX.lt.Distance(ivt)) then
      call fbm_updateFBMGeom(P(1),P(2),P(3),IP,minpr,distY,fbm_geom_handler_ptr)
      if (distY.lt.Distance(ivt)) then
       FictKNPR(ivt) = minpr
       Distance(ivt) = distY
       bMark(ivt) = .True.
      end if
      iCount(2) = iCount(2) + 1d0
     end if
    end do
   
    ! net
    do j = 1,12
     ivt = nvt + kedge(j,i)
     P = dcorvg(:,ivt)
     if (.not.bMark(ivt).or.DistX.lt.Distance(ivt)) then
      call fbm_updateFBMGeom(P(1),P(2),P(3),IP,minpr,distY,fbm_geom_handler_ptr)
      if (distY.lt.Distance(ivt)) then
       FictKNPR(ivt) = minpr
       Distance(ivt) = distY
       bMark(ivt) = .True.
      end if
      iCount(2) = iCount(2) + 1d0
     end if
    end do
    
    ! nat
    do j = 1,6
     ivt = nvt + net + karea(j,i)
     P = dcorvg(:,ivt)
     if (.not.bMark(ivt).or.DistX.lt.Distance(ivt)) then
      call fbm_updateFBMGeom(P(1),P(2),P(3),IP,minpr,distY,fbm_geom_handler_ptr)
      if (distY.lt.Distance(ivt)) then
       FictKNPR(ivt) = minpr
       Distance(ivt) = distY
       bMark(ivt) = .True.
      end if
      iCount(3) = iCount(3) + 1d0
     end if
    end do
    
   ELSE
   
    ! nvt
    do j = 1,8
     ivt = kvert(j,i)
     if ((.not.bMark(ivt)).and.(distX.lt.Distance(ivt))) then
      Distance(ivt) = distX
      FictKNPR(ivt) = minpr
     end if
    end do
   
    ! net
    do j = 1,12
     ivt = nvt + kedge(j,i)
     if ((.not.bMark(ivt)).and.(distX.lt.Distance(ivt))) then
      Distance(ivt) = distX
      FictKNPR(ivt) = minpr
     end if
    end do
    
    ! nat
    do j = 1,6
     ivt = nvt + net + karea(j,i)
     if ((.not.bMark(ivt)).and.(distX.lt.Distance(ivt))) then
      Distance(ivt) = distX
      FictKNPR(ivt) = minpr
     end if
    end do
   
   END IF
   
  END DO
  
  END DO ! IP
  
1 continue

  iCount(3) = 4*ndof
  call Comm_SummN(iCount,3)
  call ztime(tttx1)
  if (myid.eq.1) WRITE(*,'(A,ES12.4,A,I0,A,I0,A,2(ES12.4,A))') &
  'FBM time : ', tttx1-tttt0, ' s | ',int(iCount(1))+int(iCount(2)),' / ',int(iCount(3)) ,'=> ',1d2*iCount(1)/iCount(3),'%, ',1d2*iCount(2)/iCount(3),'%'
!  if (myid.eq.1) WRITE(*,'(3I10,2ES12.4)') int(iCount(1:3)),1d2*iCount(1)/iCount(3),1d2*iCount(2)/iCount(3)

  do i=1,nvt+net+nat+nel
   myALE%Monitor(i)=distance(i)
  end do
  if (myid.ne.0) then
!    do i=1,nvt+net+nat+nel
!     myALE%Monitor(i)=0d0
!     if (bMark(i))  myALE%Monitor(i)=1d0
!    end do

   DEALLOCATE(bMark)
  end if

END SUBROUTINE QuadScalar_FictKnpr_Wangen
!
! ----------------------------------------------
!
SUBROUTINE QuadScalar_InitCond()
  REAL*8 PX,PY,PZ
  INTEGER i,ndof

  ndof = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net +&
    mg_mesh%level(ilev)%nat + mg_mesh%level(ilev)%nel


  DO i=1,ndof
  PX = myQ2Coor(1,I)
  PY = myQ2Coor(2,I)
  PZ = myQ2Coor(3,I)
  CALL GetVeloInitVal(PX,PY,PZ,QuadSc%valU(i),QuadSc%valV(i),QuadSc%valW(i))
  END DO

END SUBROUTINE QuadScalar_InitCond
!
! ----------------------------------------------
!
SUBROUTINE Boundary_QuadScalar_Def()
  INTEGER i
  REAL*8 daux

  DO i=1,QuadSc%ndof

    IF (QuadSc%knprU(ILEV)%x(i).eq.1) QuadSc%defU(i) = 0d0
    IF (QuadSc%knprV(ILEV)%x(i).eq.1) QuadSc%defV(i) = 0d0
    IF (QuadSc%knprW(ILEV)%x(i).eq.1) QuadSc%defW(i) = 0d0

    IF (FictKNPR(i).ne.0) THEN
      QuadSc%defU(i) = 0d0
      QuadSc%defV(i) = 0d0
      QuadSc%defW(i) = 0d0
    END IF
    IF (MixerKNPR(i).ne.0) THEN
      QuadSc%defU(i) = 0d0
      QuadSc%defV(i) = 0d0
      QuadSc%defW(i) = 0d0
    END IF
    
    IF (myBoundary%bSlip(i).and.(.not.(myBoundary%bWall(i).or.myBoundary%iInflow(i).gt.0))) then
    
     DAUX = QuadSc%defU(i) * BoundaryNormal(1,i) + &
            QuadSc%defV(i) * BoundaryNormal(2,i) + &
            QuadSc%defW(i) * BoundaryNormal(3,i)
           
     QuadSc%defU(i) = QuadSc%defU(i) - DAUX*BoundaryNormal(1,i)
     QuadSc%defV(i) = QuadSc%defV(i) - DAUX*BoundaryNormal(2,i)
     QuadSc%defW(i) = QuadSc%defW(i) - DAUX*BoundaryNormal(3,i)
     
    END IF
       

  END DO

END SUBROUTINE Boundary_QuadScalar_Def
!
! ----------------------------------------------
!
SUBROUTINE Boundary_QuadScalar_Val()
  use fbm, only: fbm_velBCUpdate
  implicit none
  REAL*8 PX,PY,PZ,DAUX
  INTEGER i,inpr,finpr,minpr,inprU,inprV,inprW,ndof,iType

  ilev = NLMAX
  ndof = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net +&
    mg_mesh%level(ilev)%nat + mg_mesh%level(ilev)%nel


  DO i=1,ndof
    PX = myQ2Coor(1,i);  PY = myQ2Coor(2,i);  PZ = myQ2Coor(3,i)
    inpr = 0

    IF (QuadSc%knprU(ilev)%x(i).EQ.1) QuadSc%valU(i)=0d0
    IF (QuadSc%knprV(ilev)%x(i).EQ.1) QuadSc%valV(i)=0d0
    IF (QuadSc%knprW(ilev)%x(i).EQ.1) QuadSc%valW(i)=0d0

    IF (myBoundary%iInflow(i).NE.0) THEN 
      inpr = 1
      iType = myBoundary%iInflow(i)
      CALL GetVeloBCVal(PX,PY,PZ,QuadSc%valU(i),QuadSc%valV(i),QuadSc%valW(i),iType,timens)
    END IF
    finpr = FictKNPR(i)
    minpr = MixerKNPR(i)
    IF (finpr.ne.0.and.inpr.eq.0) THEN
      CALL fbm_velBCUpdate(PX,PY,PZ,&
                           QuadSc%valU(i),QuadSc%valV(i),&
                           QuadSc%valW(i),finpr,timens,&
                           fbm_vel_bc_handler_ptr)
    END IF
    IF (minpr.ne.0.and.inpr.eq.0) THEN
      CALL GetVeloMixerVal(PX,PY,PZ,QuadSc%valU(i),QuadSc%valV(i),QuadSc%valW(i),minpr,timens)
    END IF
  END DO

  DO i=1,ndof
    IF (myBoundary%bSlip(i).and.(.not.(myBoundary%bWall(i).or.myBoundary%iInflow(i).gt.0))) then
     DAUX = QuadSc%ValU(i) * BoundaryNormal(1,i) + &
            QuadSc%ValV(i) * BoundaryNormal(2,i) + &
            QuadSc%ValW(i) * BoundaryNormal(3,i)
           
     QuadSc%ValU(i) = QuadSc%ValU(i) - DAUX*BoundaryNormal(1,i)
     QuadSc%ValV(i) = QuadSc%ValV(i) - DAUX*BoundaryNormal(2,i)
     QuadSc%ValW(i) = QuadSc%ValW(i) - DAUX*BoundaryNormal(3,i)
    END IF
  END DO
  
END SUBROUTINE Boundary_QuadScalar_Val
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinScalar_Mat(DA,KLD,KNPRP,NEL,iS)
  REAL*8  DA(*)
  INTEGER KLD(*),KNPRP(*),ICOL,I,NEL,J,JJ,iS
 
  DO I=1,NEL
   IF (KNPRP(I).EQ.1) THEN
    DO J=1,4
     JJ = 4*(I-1) + J
     IF (iS.eq.1) DA(KLD(JJ)) = 1d-8
     DO ICOL=KLD(JJ)+iS,KLD(JJ+1)-1
      DA(ICOL) = 0d0
     END DO
    END DO
   END IF
  END DO

END SUBROUTINE Boundary_LinScalar_Mat
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinScalar_Def(DD,KNPRP,NEL)
  REAL*8  DD(*)
  INTEGER KNPRP(*),ICOL,I,NEL,J,JJ
 
  DO I=1,NEL
   IF (KNPRP(I).EQ.1) THEN
    DO J=1,4
     JJ = 4*(I-1) + J
     DD(JJ) = 0d0
    END DO
   END IF
  END DO

END SUBROUTINE Boundary_LinScalar_Def
!
! ----------------------------------------------
!
SUBROUTINE Boundary_QuadScalar_Mat(DA11,DA22,DA33,KLD,&
    KNPRU,KNPRV,KNPRW,NDOF)
  REAL*8  DA11(*),DA22(*),DA33(*)
  INTEGER KLD(*),KNPRU(*),KNPRV(*),KNPRW(*),ICOL,I,NDOF
  REAL*8 DAUX

  DO I=1,NDOF
  IF (KNPRU(I).EQ.1) THEN
    ICOL = KLD(I)
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = 0d0
    END DO
  END IF
  IF (KNPRV(I).EQ.1) THEN
    ICOL = KLD(I)
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA22(ICOL) = 0d0
    END DO
  END IF
  IF (KNPRW(I).EQ.1) THEN
    ICOL = KLD(I)
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA33(ICOL) = 0d0
    END DO
  END IF
  END DO

  DO I=1,NDOF
  IF (FictKNPR(I).NE.0.OR.MixerKNPR(I).NE.0) THEN
    !    DA(KLD(I))=1d-8
    DO ICOL=KLD(I)+1,KLD(I+1)-1
      DA11(ICOL)=0E0
      DA22(ICOL)=0E0
      DA33(ICOL)=0E0
    END DO
  END IF
  END DO

  DO I=1,NDOF
    IF (myBoundary%bSlip(i).and.(.not.(myBoundary%bWall(i).or.myBoundary%iInflow(i).gt.0))) then
    ICOL = KLD(I)
    DO ICOL=KLD(I)+1,KLD(I+1)-1
      DA11(ICOL)=0E0
      DA22(ICOL)=0E0
      DA33(ICOL)=0E0
    END DO
     
    END IF
  END DO

  
  
END SUBROUTINE Boundary_QuadScalar_Mat
!
! ----------------------------------------------
!
SUBROUTINE Boundary_QuadScalar_Mat_9(DA11,DA22,DA33,DA12,DA13,DA23,DA21,DA31,DA32,&
    KLD,KNPRU,KNPRV,KNPRW,NDOF)
  REAL*8 DA11(*),DA22(*),DA33(*),DA12(*),DA13(*),DA23(*),DA21(*),DA31(*),DA32(*)
  INTEGER KLD(*),KNPRU(*),KNPRV(*),KNPRW(*),ICOL,I,NDOF

  DO I=1,NDOF
  IF (KNPRU(I).EQ.1) THEN
    ICOL = KLD(I)
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = 0d0
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    END DO
  END IF
  IF (KNPRV(I).EQ.1) THEN
    ICOL = KLD(I)
    DA23(ICOL) = 0d0
    DA21(ICOL) = 0d0
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA22(ICOL) = 0d0
    DA23(ICOL) = 0d0
    DA21(ICOL) = 0d0
    END DO
  END IF
  IF (KNPRW(I).EQ.1) THEN
    ICOL = KLD(I)
    DA31(ICOL) = 0d0
    DA32(ICOL) = 0d0
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA33(ICOL) = 0d0
    DA31(ICOL) = 0d0
    DA32(ICOL) = 0d0
    END DO
  END IF
  END DO

  DO I=1,NDOF
  IF (FictKNPR(I).NE.0.OR.MixerKNPR(I).NE.0) THEN
    ICOL = KLD(I)
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    DA23(ICOL) = 0d0
    DA21(ICOL) = 0d0
    DA31(ICOL) = 0d0
    DA32(ICOL) = 0d0
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = 0d0
    DA22(ICOL) = 0d0
    DA33(ICOL) = 0d0
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    DA23(ICOL) = 0d0
    DA21(ICOL) = 0d0
    DA31(ICOL) = 0d0
    DA32(ICOL) = 0d0
    END DO
  END IF
  END DO

END SUBROUTINE Boundary_QuadScalar_Mat_9
!
! ----------------------------------------------
!
SUBROUTINE Velocity_Correction()
  INTEGER i

  ! *** Update of U = U~ - k M^-1 B P 

  ILEV = NLMAX
  qlMat     => mg_qlMat(ILEV)
  BXMat     => mg_BXMat(ILEV)%a
  BYMat     => mg_BYMat(ILEV)%a
  BZMat     => mg_BZMat(ILEV)%a
  MlRhoPmat => mg_MlRhoPmat(ILEV)%a
  
  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%ValP(NLMAX)%x,&
    QuadSc%defU,QuadSc%defV,QuadSc%defW,QuadSc%ndof,-TSTEP,0d0)

  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)

  CALL Boundary_QuadScalar_Def()

  DO I=1,QuadSc%ndof
  QuadSc%valU(i) = QuadSc%valU(i) - QuadSc%defU(i)/MlRhoPmat(i)
  QuadSc%valV(i) = QuadSc%valV(i) - QuadSc%defV(i)/MlRhoPmat(i)
  QuadSc%valW(i) = QuadSc%valW(i) - QuadSc%defW(i)/MlRhoPmat(i)
  END DO

END SUBROUTINE Velocity_Correction
!
! ----------------------------------------------
!
SUBROUTINE Pressure_Correction()
  INTEGER I
  real*8 dR

  DO I=1,lMat%nu
   LinSc%valP(NLMAX)%x(i) = LinSc%valP(NLMAX)%x(i) + LinSc%valP_old(i)
   LinSc%P_new(i) = 1.5d0*LinSc%valP(NLMAX)%x(i) - 0.5d0*LinSc%valP_old(i)
  END DO

  ! ! Set dirichlet boundary conditions on the solution
  ! CALL Boundary_LinScalar_Val(DWORK(L(LCORVG)))

END SUBROUTINE Pressure_Correction
!
! ----------------------------------------------
!
SUBROUTINE AddPressureGradient()
  INTEGER I,J,IEL
  REAL*8 ddx,ddy,ddz,ddp

  ILEV = NLMAX
  qlMat    => mg_qlMat(ILEV)
  BXMat    => mg_BXMat(ILEV)%a
  BYMat    => mg_BYMat(ILEV)%a
  BZMat    => mg_BZMat(ILEV)%a
  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%valP(NLMAX)%x,&
    QuadSc%defU,QuadSc%defV,QuadSc%defW,QuadSc%ndof,TSTEP,1d0)

END SUBROUTINE AddPressureGradient
!
! ----------------------------------------------
!
SUBROUTINE AddPeriodicPressureGradient()
  INTEGER I,J,IEL
  REAL*8 ddx,ddy,ddz,ddp

  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%valP(NLMAX)%x,&
  QuadSc%defU,QuadSc%defV,QuadSc%defW,QuadSc%ndof,TSTEP,1d0)

  DO I=1,QuadSc%ndof
   QuadSc%defW(i) = QuadSc%defW(i) + TSTEP*dPeriodicVector(i)
  END DO
  
END SUBROUTINE AddPeriodicPressureGradient
!
! ----------------------------------------------
!
SUBROUTINE AddGravForce
  INTEGER I,LINT
  REAL*8 DAUX
  EXTERNAL E013

  ILEV=NLMAX
  CALL SETLEV(2)

  DAUX=SQRT(Properties%Gravity(1)**2d0+Properties%Gravity(2)**2d0+Properties%Gravity(3)**2d0)

  IF (DAUX.GT.0d0) THEN
    CALL Grav_QuadSc(QuadSc%defU,QuadSc%defV,QuadSc%defW,mgDensity(ILEV)%x,&
      Properties%Gravity,QuadSc%ndof,&
      mg_mesh%level(ilev)%kvert,&
      mg_mesh%level(ilev)%karea,&
      mg_mesh%level(ilev)%kedge,&
      KWORK(L(KLINT(NLMAX))),&
      mg_mesh%level(ilev)%dcorvg,&
      tstep,E013)

  END IF

END SUBROUTINE AddGravForce
!
! ----------------------------------------------
!
SUBROUTINE OperatorRegenaration(iType)
  INTEGER iType
  LOGICAL bHit

  bHit = .FALSE.

  IF (iType.EQ.myMatrixRenewal%D) THEN
    CALL Create_DiffMat(QuadSc)
    bHit = .TRUE.
  END IF

  IF (iType.EQ.myMatrixRenewal%K) THEN
    CALL Create_KMat(QuadSc)
    bHit = .TRUE.
  END IF

  IF (iType.EQ.myMatrixRenewal%M) THEN
    CALL Create_MRhoMat()
    bHit = .TRUE.
  END IF

  IF (iType.EQ.myMatrixRenewal%S) THEN
    CALL Create_SMat(QuadSc)
    bHit = .TRUE.
  END IF

  IF (iType.EQ.myMatrixRenewal%C) THEN

    CALL Create_BMat()

    IF (myid.ne.master) THEN
      CALL Fill_QuadLinParMat()
    END IF
    
    CALL SetSlipOnBandBT()
    
    CALL Create_CMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW,LinSc%knprP,LinSc%prm%MGprmIn%MinLev,LinSc%prm%MGprmIn%CrsSolverType)
    IF (myid.ne.master) THEN
      CALL Create_ParCMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW)
    END IF
    bHit = .TRUE.
  END IF

  IF (iType.EQ.1.and.bNS_Stabilization) then
   CALL Create_hDiffMat()
   bHit = .TRUE.
  END IF

  IF (iType.EQ.3.and.(NewtonForBurgers.ne.0d0)) then
   CALL Create_barMMat_iso(QuadSc)
   bHit = .TRUE.
  END IF

  IF (myid.EQ.ShowID.AND.bHit) WRITE(MTERM,'(A)', advance='yes') " "

END SUBROUTINE OperatorRegenaration
!
! ----------------------------------------------
!
SUBROUTINE ProlongateSolution()

 IF (allocated(Temperature)) then
  CALL ProlongateSolutionSub(QuadSc,LinSc,Boundary_QuadScalar_Val,Temperature)
 else
  CALL ProlongateSolutionSub(QuadSc,LinSc,Boundary_QuadScalar_Val)
 end if
 CALL QuadScP1toQ2(LinSc,QuadSc)

END SUBROUTINE ProlongateSolution
!
! ----------------------------------------------
!
SUBROUTINE InitBoundaryList(knpr,kvert,kedge,karea)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*),knpr(*)
  INTEGER i,j,k,ivt1,ivt2
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

  DO i=1,nvt
  QuadScBoundary(i) = mg_mesh%level(ilev)%knpr(i)
  !  IF (QuadScBoundary(i).eq.1) write(*,*) "type 1"
  END DO

  k=1
  DO i=1,nel
  DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
    ivt1 = kvert(NeighE(1,j),i)
    ivt2 = kvert(NeighE(2,j),i)
    !IF (knpr(ivt1).EQ.1.AND.knpr(ivt2).EQ.1) THEN
    QuadScBoundary(nvt+k) = 0
    !END IF
    !    IF (QuadScBoundary(nvt+k).eq.1) write(*,*) "type 2"
    k = k + 1
  END IF
  END DO
  END DO

  DO i=1,nat
  !QuadScBoundary(nvt+net+i) = knpr(nvt+i)
  QuadScBoundary(nvt+net+i) = 0! knpr(nvt+i)
  !  IF (QuadScBoundary(nvt+net+i).eq.1) write(*,*) "type 3"
  END DO

  DO i=1,nel
  QuadScBoundary(nvt+net+nat+i) = 0
  END DO

END SUBROUTINE InitBoundaryList
!
! ----------------------------------------------
!
SUBROUTINE TESTER(DX,DXP)
  REAL*8 DX(*),DXP(*)

  CALL GetParPressure(DX,DXP)

END SUBROUTINE TESTER
!
! ----------------------------------------------
!
SUBROUTINE FBM_GetForces()

  CALL EvaluateDragLift(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
    LinSc%valP(NLMAX)%x,Viscosity,FictKNPR)

END SUBROUTINE FBM_GetForces
!
! ----------------------------------------------
!
SUBROUTINE FAC_GetForcesParT(mfile,iT)
  INTEGER mfile,iT
  REAL*8 :: Force2(3),ForceV(3),ForceP(3),Force(3),Factor,PI = 3.141592654d0
  REAL*8 :: Scale
  REAL*8 :: U_mean=1.0d0,H=0.20d0,D=1d0
  INTEGER i,nn
  EXTERNAL E013

  ILEV=NLMAX
  CALL SETLEV(2)

  IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN 
    CALL EvaluateDragLift9_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      LinSc%P_new,BndrForce,ForceV,ForceP)
  ELSE
    CALL EvaluateDragLift_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      LinSc%P_new,BndrForce,ForceV,ForceP)
  END IF

  Factor = 2d0/(postParams%U_mean*postParams%U_mean*postParams%D*postParams%H)
  ForceP = Factor*ForceP
  ForceV = Factor*ForceV

  IF (myid.eq.showID) THEN
   write(mfile,'(A16,3ES15.7E2,I3)') "BenchForceParT: ",timens,ForceV(1:2)+forceP(1:2),iT
   write(MTERM,'(A16,3ES15.7E2,I3)') "BenchForceParT: ",timens,ForceV(1:2)+forceP(1:2),iT
   WRITE(MTERM,5)
   WRITE(MFILE,5)

!    WRITE(666,'(10ES16.8,I3)') Timens,ForceV+forceP,ForceV,forceP,iT
  END IF

  5  FORMAT(104('-'))

END SUBROUTINE FAC_GetForcesParT
!
! ----------------------------------------------
!
SUBROUTINE FAC_GetForces(mfile)
  INTEGER mfile
  REAL*8 :: Force2(3),ForceV(3),ForceP(3),Force(3),Factor,PI = 3.141592654d0
  REAL*8 :: Scale
  REAL*8 :: U_mean=1.0d0,H=0.20d0,D=1d0
  INTEGER i,nn
  EXTERNAL E013

  ILEV=NLMAX
  CALL SETLEV(2)

  IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN 
    CALL EvaluateDragLift9_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      LinSc%P_new,BndrForce,ForceV,ForceP)
  ELSE
    CALL EvaluateDragLift_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      LinSc%P_new,BndrForce,ForceV,ForceP)
  END IF

  if(bViscoElastic)then
    Force = ForceV + ForceP + ViscoElasticForce
    IF (b2DViscoBench) THEN
     Scale = 2d0/(U_mean*U_mean*D*H)
     Force = Scale*Force
     ViscoElasticForce = (ViscoElasticForce)*Scale
    END IF
    IF (b3DViscoBench) THEN
     Scale = 6d0*PI*postParams%Sc_Mu*postParams%Sc_U*postParams%Sc_a
     Force = (4d0*Force)/Scale
     ViscoElasticForce = (4d0*ViscoElasticForce)/Scale
    END IF
  else
    Factor = 2d0/(postParams%U_mean*postParams%U_mean*postParams%D*postParams%H)
    ForceP = Factor*ForceP
    ForceV = Factor*ForceV
  end if

  IF (myid.eq.showID) THEN
   
    if(bViscoElastic)then
      WRITE(MTERM,5)
      WRITE(MFILE,5)
      IF (b2DViscoBench) THEN
       write(mfile,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(1),(Force(1)-ViscoElasticForce(1)),Force(1)
       write(mterm,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(1),(Force(1)-ViscoElasticForce(1)),Force(1)
       WRITE(666,'(10ES13.5)')timens,ViscoElasticForce,&
         (Force-ViscoElasticForce),Force
      END IF
      IF (b3DViscoBench) THEN
       write(mfile,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(3),(Force(3)-ViscoElasticForce(3)),Force(3)
       write(mterm,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(3),(Force(3)-ViscoElasticForce(3)),Force(3)
       WRITE(666,'(10ES13.5)')timens,ViscoElasticForce,&
         (Force-ViscoElasticForce),Force
      END IF
    else
      WRITE(MTERM,5)
      WRITE(MFILE,5)
       write(mfile,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
       write(*,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
    end if
      write(mfile,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
      write(*,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
      write(mfile,'(A12,3ES15.7E2)') "BenchForce: ",timens,ForceV(1:2)+forceP(1:2)
      write(*,'(A12,3ES15.7E2)') "BenchForce: ",timens,ForceV(1:2)+forceP(1:2)
!      write(*,'(A12,3ES15.7E2)') "FBMForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2)

!      write(mfile,'(A12,3ES15.7E2)') "TotalForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2) + &
!      ForceV(1:2)+forceP(1:2)
!      write(*,'(A12,3ES15.7E2)') "TotalForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2) + &
!      ForceV(1:2)+forceP(1:2)

      WRITE(666,'(10ES16.8)') Timens,ForceV+forceP,ForceV,forceP
      WRITE(MTERM,5)
      WRITE(MFILE,5)
  END IF

  5  FORMAT(104('-'))

END SUBROUTINE FAC_GetForces
!
! ----------------------------------------------
!
SUBROUTINE GetPressureSample(dcorvg,nvt)
  INTEGER NVT
  REAL*8 dcorvg(3,*)
  INTEGER I,J1,J2
  REAL*8 :: P1X=0.55d0,P1Y=0.20d0,P1Z=0.205d0
  REAL*8 :: P2X=0.45d0,P2Y=0.20d0,P2Z=0.205d0
  REAL*8 DIST1,DIST2,MINDIST1,MINDIST2,PX,PY,PZ

  IF (myid.ne.0) THEN
    MINDIST1 = 1d30
    MINDIST2 = 1d30
    DO i=1,nvt
    PX = dcorvg(1,I)
    PY = dcorvg(2,I)
    PZ = dcorvg(3,I)
    DIST1 = SQRT((PX-P1X)**2d0 + (PY-P1Y)**2d0 + (PZ-P1Z)**2d0)
    DIST2 = SQRT((PX-P2X)**2d0 + (PY-P2Y)**2d0 + (PZ-P2Z)**2d0)
    IF (DIST1.LT.MINDIST1) THEN
      MINDIST1 = DIST1
      J1 = I
    END IF
    IF (DIST2.LT.MINDIST2) THEN
      MINDIST2 = DIST2
      J2 = I
    END IF
    END DO
    MINDIST1 = -MINDIST1
    MINDIST2 = -MINDIST2
    DIST1 = MINDIST1
    DIST2 = MINDIST2
    !  WRITE(*,*) myid,J1,J2,MINDIST1,MINDIST2!PressureSample
  END IF

  CALL COMM_Maximum(MINDIST1)
  CALL COMM_Maximum(MINDIST2)

  IF (myid.ne.0) THEN
    !  WRITE(*,*) myid,J1,J2,MINDIST1,MINDIST2!PressureSample
    PressureSample = 0
    IF (DIST1.EQ.MINDIST1) THEN
      PressureSample(1) = J1
    END IF
    IF (DIST2.EQ.MINDIST2) THEN
      PressureSample(2) = J2
    END IF
    !  WRITE(*,*) myid,MINDIST1,MINDIST2,PressureSample
  END IF

END SUBROUTINE GetPressureSample


SUBROUTINE Analyzer
  INTEGER I,J

  J=0
  DO I=1,QuadSc%ndof
  IF (ABS(QuadSC%auxU(I)).GT.1d-10) THEN
    J = J + 1
    !   WRITE(*,*) I,J,QuadSC%auxU(I)
  END IF
  END DO

END SUBROUTINE Analyzer
!
! ----------------------------------------------
!
SUBROUTINE updateFBMGeometry_Wangen()
use cinterface, only: calculateFBM

 if (calculateFBM()) then
  if (myid.eq.showid) write(*,*) '> FBM computation step'

  ILEV=NLMAX
  CALL SETLEV(2)

  CALL QuadScalar_FictKnpr_Wangen(mg_mesh%level(ilev)%dcorvg,&
    mg_mesh%level(ilev)%dcorag,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%karea)
    
  CALL E013Max_SUPER(FictKNPR)
  
 end if

END SUBROUTINE  updateFBMGeometry_Wangen
!
! ----------------------------------------------
!
SUBROUTINE updateFBMGeometry()
use cinterface, only: calculateFBM

 if (calculateFBM()) then
  if (myid.eq.showid) write(*,*) '> FBM computation step'

  ILEV=NLMAX
  CALL SETLEV(2)

  CALL QuadScalar_FictKnpr(mg_mesh%level(ilev)%dcorvg,&
    mg_mesh%level(ilev)%dcorag,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%karea)
    
  CALL E013Max_SUPER(FictKNPR)
  
 end if

END SUBROUTINE  updateFBMGeometry
!
! ----------------------------------------------
!
SUBROUTINE updateMixerGeometry(mfile)
use geometry_processing, only : calcDistanceFunction_sse, QuadScalar_MixerKnpr,&
    calcDistanceFunction_netzsch,dEpsDist

integer, intent(in) :: mfile

integer :: i

REAL :: tttt0,tttt1

CALL myMPI_Barrier()
CALL ZTIME(tttt0)

ILEV=NLMAX
CALL SETLEV(2)
QuadSc%AuxU = dEpsDist
QuadSc%AuxV = dEpsDist

MixerKNPR(:) = 0

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 CALL calcDistanceFunction_sse(mg_mesh%level(ilev)%dcorvg,&
                           mg_mesh%level(ilev)%kvert,&
                           mg_mesh%level(ilev)%kedge,&
                           mg_mesh%level(ilev)%karea,&
                           mg_mesh%level(ilev)%nel,&
                           mg_mesh%level(ilev)%nvt,&
                           mg_mesh%level(ilev)%nat,&
                           mg_mesh%level(ilev)%net,&
                           QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."NETZSCH") THEN
 CALL calcDistanceFunction_netzsch(mg_mesh%level(ilev)%dcorvg,&
                           mg_mesh%level(ilev)%kvert,&
                           mg_mesh%level(ilev)%kedge,&
                           mg_mesh%level(ilev)%karea,&
                           mg_mesh%level(ilev)%nel,&
                           mg_mesh%level(ilev)%nvt,&
                           mg_mesh%level(ilev)%nat,&
                           mg_mesh%level(ilev)%net,&
                           QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
 CALL QuadScalar_MixerKnpr(mg_mesh%level(ilev)%dcorvg,&
                           mg_mesh%level(ilev)%kvert,&
                           mg_mesh%level(ilev)%kedge,&
                           mg_mesh%level(ilev)%karea,&
                           mg_mesh%level(ilev)%nel,&
                           mg_mesh%level(ilev)%nvt,&
                           mg_mesh%level(ilev)%nat,&
                           mg_mesh%level(ilev)%net,&
                           QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW)
END IF


CALL myMPI_Barrier()
CALL ZTIME(tttt1)
IF (myid.eq.1) WRITE(mterm,"(A,F6.1,A)") "Time used for FINE mesh distance estimation was: ", tttt1-tttt0, "s!"
IF (myid.eq.1) WRITE(mfile,"(A,F6.1,A)") "Time used for FINE mesh distance estimation was: ", tttt1-tttt0, "s!"

END SUBROUTINE  updateMixerGeometry
!
! ----------------------------------------------
!
SUBROUTINE StaticMeshAdaptation()
  INTEGER iAdaptMeshLevel,idL

  IF (.NOT.bMeshAdaptation) RETURN

  OPEN(474,FILE=ADJUSTL(TRIM(cAdaptedMeshFile))//"/level.prf")
  READ(474,*) iAdaptMeshLevel
  CLOSE(474)

  idL = iAdaptMeshLevel - NLMAX

  CALL CreateDumpStructures(idL)
  CALL LoadSmartAdaptedMeshFile(DWORK(L(KLCVG(1))),cAdaptedMeshFile,idL)

  ! ---------------------------- -  -    - -- -- - - - - - - - -   ----------------

  DO ILEV = iAdaptMeshLevel,NLMAX
  CALL SETLEV(2)
  WRITE(*,*) 'mesh levels', ilev,iAdaptMeshLevel,NLMAX
  CALL RefreshCoordinates(DWORK(L(KLCVG(ILEV+1))),DWORK(L(KLCAG(ILEV))),&
    KWORK(L(KLVERT(ILEV))),KWORK(L(KLEDGE(ILEV))),KWORK(L(KLAREA(ILEV))))
  END DO

  ! ---------------------------- -  -    - -- -- - - - - - - - -   ----------------

END SUBROUTINE StaticMeshAdaptation
!
! ----------------------------------------------
!
SUBROUTINE CoorWriter(dcorvg,nvt,cF)
  CHARACTER*(*) cF
  REAL*8 dcorvg(3,*)
  INTEGER i,nvt

  OPEN(547,FILE=TRIM(ADJUSTL(cF)))
  DO i=1,nvt
  WRITE(547,*) dcorvg(:,i)
  END DO
  CLOSE(547)

END SUBROUTINE CoorWriter
!
! ----------------------------------------------
!
SUBROUTINE RefreshCoordinates(dcorvg,dcorag,kvert,kedge,karea)
  REAL*8  dcorvg(3,*),dcorag(3,*)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
  REAL*8 PX,PY,PZ,DIST
  INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

  k=1
  DO i=1,nel
  DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
    ivt1 = kvert(NeighE(1,j),i)
    ivt2 = kvert(NeighE(2,j),i)
    PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
    PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
    PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
    dcorvg(1,nvt+k) = PX
    dcorvg(2,nvt+k) = PY
    dcorvg(3,nvt+k) = PZ
    k = k + 1
  END IF
  END DO
  END DO

  k=1
  DO i=1,nel
  DO j=1,6
  IF (k.eq.karea(j,i)) THEN
    ivt1 = kvert(NeighA(1,j),i)
    ivt2 = kvert(NeighA(2,j),i)
    ivt3 = kvert(NeighA(3,j),i)
    ivt4 = kvert(NeighA(4,j),i)
    PX = 0.25d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2)+dcorvg(1,ivt3)+dcorvg(1,ivt4))
    PY = 0.25d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2)+dcorvg(2,ivt3)+dcorvg(2,ivt4))
    PZ = 0.25d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2)+dcorvg(3,ivt3)+dcorvg(3,ivt4))
    dcorag(1,k) = PX
    dcorag(2,k) = PY
    dcorag(3,k) = PZ
    dcorvg(1,nvt+net+k) = PX
    dcorvg(2,nvt+net+k) = PY
    dcorvg(3,nvt+net+k) = PZ
    k = k + 1
  END IF
  END DO
  END DO

  DO i=1,nel
  PX = 0d0
  PY = 0d0
  PZ = 0d0
  DO j=1,8
  PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
  PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
  PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
  END DO
  dcorvg(1,nvt+net+nat+i) = PX
  dcorvg(2,nvt+net+nat+i) = PY
  dcorvg(3,nvt+net+nat+i) = PZ
  END DO

END SUBROUTINE RefreshCoordinates
!
! ----------------------------------------------
!
SUBROUTINE  GetMonitor()
  implicit none
  INTEGER i
  REAL*8 daux,px,py,pz
  return
  ILEV = NLMAX
  CALL SETLEV(2)

  !  
  ! 
  !  DO i=1,nvt
  !  
  !    PX = dcorvg(1,i)
  !    PY = dcorvg(2,i)
  !    PZ = dcorvg(3,i)
  !    getdistanceid(px,py,pz,daux,i);
  !   
  !   myALE%Monitor(i) = sqrt(daux) !HogenPowerlaw(daux)
  ! 
  !  END DO

END SUBROUTINE  GetMonitor
!
! ----------------------------------------------
!
SUBROUTINE  GetNonNewtViscosity()
INTEGER i
REAL*8 ViscosityModel
integer :: ilevel

EXTERNAL E013

ilev   = mg_mesh%nlmax
ilevel = mg_mesh%nlmax

 IF (myid.ne.0) then
  QuadSc%defU = 0d0
  QuadSc%defV = 0d0
  QuadSc%defW = 0d0
  CALL L2ProjVisco(QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                   Temperature,MaterialDistribution(ilev)%x,&
                   QuadSc%defU,QuadSc%defV,QuadSc%defW,&
                   mg_mesh%level(ilevel)%kvert,&
                   mg_mesh%level(ilevel)%karea,&
                   mg_mesh%level(ilevel)%kedge,&
                   mg_mesh%level(ilevel)%dcorvg,E013)

  CALL E013Sum(QuadSc%defU)
  CALL E013Sum(QuadSc%defV)
  CALL E013Sum(QuadSc%defW)

  if(.not.allocated(Shearrate)) allocate(Shearrate(QuadSc%ndof))
  
  DO i=1,QuadSc%ndof
   Shearrate(i) = QuadSc%defV(i)/QuadSc%defW(i)
   Viscosity(i) = QuadSc%defU(i)/QuadSc%defW(i)
  END DO

 END IF

END SUBROUTINE  GetNonNewtViscosity
!
! ----------------------------------------------
!
SUBROUTINE ExtractElemSizeDistribution()
real*8 x(8),y(8),z(8),dVol
real*8, allocatable :: daux(:),daux2(:)
integer iel,i,nQ2_dof

if (myid.ne.0) then

 ILEV = NLMAX
 CALL SETLEV(2)
 
 nQ2_dof = mg_mesh%level(ilev)%nvt + &
           mg_mesh%level(ilev)%net + &
           mg_mesh%level(ilev)%nat + &
           mg_mesh%level(ilev)%nel 
           
 allocate(daux(nQ2_dof))
 allocate(daux2(nQ2_dof))
 IF (.not.allocated(ElemSizeDist)) allocate(ElemSizeDist(mg_mesh%level(ilev)%nvt))
 
 daux2 = 0d0
 daux  = 0d0
 
 do iel=1,mg_mesh%level(ilev)%nel
  
  x(:) = mg_mesh%level(ilev)%dcorvg(1,mg_mesh%level(ilev)%kvert(:,iel))
  y(:) = mg_mesh%level(ilev)%dcorvg(2,mg_mesh%level(ilev)%kvert(:,iel))
  z(:) = mg_mesh%level(ilev)%dcorvg(3,mg_mesh%level(ilev)%kvert(:,iel))
  CALL GetElemVol(x,y,z,dVol)

  daux2(mg_mesh%level(ilev)%kvert(:,iel)) = daux2(mg_mesh%level(ilev)%kvert(:,iel)) + 0.125d0*(dVol**0.333d0)
  daux(mg_mesh%level(ilev)%kvert(:,iel))  = daux(mg_mesh%level(ilev)%kvert(:,iel)) + 0.125d0
 
 end do
  
 CALL E013Sum(daux2)
 CALL E013Sum(daux)

 do i=1,mg_mesh%level(ilev)%nvt
  ElemSizeDist(i) = daux2(i)/daux(i)
 end do
  
 deallocate(daux,daux2)

end if
 
END SUBROUTINE ExtractElemSizeDistribution
!
! ----------------------------------------------
!
SUBROUTINE ExtractVeloGradients()

  ILEV = NLMAX
  CALL SETLEV(2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValU)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,1)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValV)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValW)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,3)

END SUBROUTINE ExtractVeloGradients
!
! ----------------------------------------------
!
SUBROUTINE  GetNonNewtViscosity_sse()
  INTEGER i
  REAL*8 daux,taux
  REAL*8 HogenPowerlaw
  REAL*8 ViscosityMatModel

  ILEV = NLMAX
  CALL SETLEV(2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValU)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,1)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValV)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValW)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,3)

  DO i=1,SIZE(QuadSc%ValU)
   daux = QuadSc%ValUx(i)**2d0 + QuadSc%ValVy(i)**2d0 + QuadSc%ValWz(i)**2d0 + &
     0.5d0*(QuadSc%ValUy(i)+QuadSc%ValVx(i))**2d0 + &
     0.5d0*(QuadSc%ValUz(i)+QuadSc%ValWx(i))**2d0 + &
     0.5d0*(QuadSc%ValVz(i)+QuadSc%ValWy(i))**2d0
   taux = Temperature(i)

   Shearrate(i) = sqrt(2d0 * daux)
   Viscosity(i) = ViscosityMatModel(daux,1,taux)

  END DO

END SUBROUTINE  GetNonNewtViscosity_sse
!
! ----------------------------------------------
!
SUBROUTINE  GetAlphaNonNewtViscosity_sse()
  INTEGER i
  REAL*8 daux,taux,dAlpha
  REAL*8 AlphaViscosityMatModel,dMaxMat
  integer ifld,iMat

  ILEV = NLMAX
  CALL SETLEV(2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValU)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,1)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValV)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValW)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,3)

  DO i=1,SIZE(QuadSc%ValU)
   daux = QuadSc%ValUx(i)**2d0 + QuadSc%ValVy(i)**2d0 + QuadSc%ValWz(i)**2d0 + &
     0.5d0*(QuadSc%ValUy(i)+QuadSc%ValVx(i))**2d0 + &
     0.5d0*(QuadSc%ValUz(i)+QuadSc%ValWx(i))**2d0 + &
     0.5d0*(QuadSc%ValVz(i)+QuadSc%ValWy(i))**2d0

   taux   = GenLinScalar%Fld(1)%val(i)
   
   iMat = myMultiMat%InitMaterial
   dMaxMat = 1d-5
   do iFld=2,GenLinScalar%nOfFields
    if (GenLinScalar%Fld(iFld)%val(i).gt.dMAxMat) then
     iMat = iFld-1
     dMaxMat = GenLinScalar%Fld(iFld)%val(i)
    end if
   end do
!   dalpha = GenLinScalar%Fld(2)%val(i)
     
   Shearrate(i) = sqrt(2d0 * daux)
   Viscosity(i) = AlphaViscosityMatModel(daux,iMat,taux)

  END DO

END SUBROUTINE  GetAlphaNonNewtViscosity_sse
!
! ----------------------------------------------
!
SUBROUTINE Calculate_Torque(mfile)
implicit none
INTEGER mfile,i
REAL*8 Torque1(3), Torque2(3),dVolFlow1,dVolFlow2,myPI,daux
REAL*8 dHeat,Ml_i,Shear,Visco,dVol,dArea1,dArea2
REAL*8 dIntPres1,dIntPres2,dPressureDifference,zMin, zMax


integer :: ilevel

EXTERNAL E013

ilevel = mg_mesh%nlmax

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque1, E013,101)
                     
 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque2, E013,102)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE") THEN
 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque1, E013,103)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque1, E013,103)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 IF (myid.ne.0) then
  call Integrate_DIE_Flowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow1,0)

  call Integrate_DIE_Flowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow2,1)
 END IF
ELSE
 IF (myid.ne.0) then
  call IntegrateFlowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow1,mySigma%L0)

  call IntegrateFlowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow2,mySigma%L0+mySigma%L)
 END IF
END IF

IF (myid.ne.0) then
 call IntegratePressure(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dIntPres1,dArea1,mySigma%L0)

 call IntegratePressure(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dIntPres2,dArea2,mySigma%L0+mySigma%L)
END IF

dHeat = 0d0
dVol  = 0d0

DO i=1,QuadSc%ndof
 IF (MixerKNPR(i).eq.0) THEN
  Shear = Shearrate(i)
  Visco = 0.1d0*Viscosity(i)
  Ml_i = mg_MlRhoMat(NLMAX)%a(i)*1e-6
  dHeat = dHeat + Ml_i*Shear*Shear*Visco
  dVol = dVol + mg_MlRhoMat(NLMAX)%a(i)*1e-3
 END IF
END DO

CALL COMM_SUMM(dVolFlow1)
CALL COMM_SUMM(dVolFlow2)
CALL COMM_SUMM(dArea1)
CALL COMM_SUMM(dArea2)
CALL COMM_SUMM(dIntPres1)
CALL COMM_SUMM(dIntPres2)
CALL COMM_SUMM(dHeat)
CALL COMM_SUMM(dVol)
dPressureDifference = ((dIntPres2/dArea2) - (dIntPres1/dArea1))*1e-6 !!! [Bar]

myPI = dATAN(1d0)*4d0
daux = 1D0*1e-7*myPI*(myProcess%umdr/3d1)

dHeat = dHeat / myThermodyn%density
dVol = dVol / myThermodyn%density

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."XSE") THEN

 mySSE_covergence%Monitor(mySSE_covergence%iC) = dPressureDifference
 mySSE_covergence%iC = MOD(mySSE_covergence%iC,mySSE_covergence%nC) + 1

 IF (itns.gt.mySSE_covergence%nC) THEN
  mySSE_covergence%average = 0d0
  DO i=1,mySSE_covergence%nC
   mySSE_covergence%average = mySSE_covergence%average +  mySSE_covergence%Monitor(i)
  END DO
  mySSE_covergence%average = mySSE_covergence%average/DBLE(mySSE_covergence%nC)
  if (mySSE_covergence%average.lt.0d0) then
   mySSE_covergence%average = MIN(mySSE_covergence%average,(1d1*mySSE_covergence%dCharVisco*1d-5))
  else
   mySSE_covergence%average = MAX(mySSE_covergence%average,(1d1*mySSE_covergence%dCharVisco*1d-5))
  end if

  mySSE_covergence%std_dev = 0d0
  DO i=1,mySSE_covergence%nC
   mySSE_covergence%std_dev = mySSE_covergence%std_dev +  (mySSE_covergence%Monitor(i)-mySSE_covergence%average)**2d0
  END DO
  mySSE_covergence%std_dev = (mySSE_covergence%std_dev/DBLE(mySSE_covergence%nC))**0.5d0
  mySSE_covergence%std_dev = ABS(1d2*mySSE_covergence%std_dev/mySSE_covergence%average)
 END IF

 mySetup%bPressureConvergence = .false.

 if (itns.gt.mySSE_covergence%start) then
  if (mySSE_covergence%std_dev.lt.2d-1) then
   mySetup%bPressureConvergence = .true.
  end if
 end if

END IF
 
IF (myid.eq.showID) THEN
  WRITE(MTERM,5)
  WRITE(MFILE,5)
  write(mfile,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_Volume_[l]_&_RT_[s]:",timens,dVolFlow1*3.6d0,dVolFlow1*3.6d0*myThermodyn%density,dVol,dVol/dVolFlow1*1000d0
  write(mterm,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_Volume_[l]_&_RT_[s]:",timens,dVolFlow1*3.6d0,dVolFlow1*3.6d0*myThermodyn%density,dVol,dVol/dVolFlow1*1000d0
  write(mfile,'(A,7ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_PressureDiff[bar]_stdERR[%]: ",timens,dVolFlow2*3.6d0,dVolFlow2*3.6d0*myThermodyn%density,dPressureDifference,mySSE_covergence%std_dev
  write(mterm,'(A,7ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_PressureDiff[bar]_stdERR[%]: ",timens,dVolFlow2*3.6d0,dVolFlow2*3.6d0*myThermodyn%density,dPressureDifference,mySSE_covergence%std_dev
IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
  write(mfile,'(A,3ES14.4,A,ES14.4)') "Power_acting_on_the_screws_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),1e-3*daux*Torque2(3),' & ',1e-3*dHeat
  write(mterm,'(A,3ES14.4,A,ES14.4)') "Power_acting_on_the_screws_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),1e-3*daux*Torque2(3),' & ',1e-3*dHeat
END IF
IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE") THEN
  write(mfile,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
  write(mterm,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
END IF
IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
  write(mfile,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
  write(mterm,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
END IF
!  WRITE(666,'(7G16.8)') Timens,Torque1,Torque2 
END IF

5  FORMAT(100('-'))

END SUBROUTINE Calculate_Torque
!
! ----------------------------------------------
!
SUBROUTINE Integrate_DIE_Flowrate(dcorvg,karea,kvert,nel,dVolFlow,iPar)
REAL*8 dcorvg(3,*),dVolFlow
INTEGER karea(6,*),kvert(8,*),nel,iPar
!---------------------------------
INTEGER NeighA(4,6)
REAL*8 P(3),dAN(3),dV
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4

dVolFlow = 0d0

if (iPar.eq.0) then
 k=1
 DO i=1,nel
  DO j=1,6
   IF (k.eq.karea(j,i)) THEN
    IF (myBoundary%iInflow(nvt+net+k).ne.0) THEN
     ivt1 = kvert(NeighA(1,j),i)
     ivt2 = kvert(NeighA(2,j),i)
     ivt3 = kvert(NeighA(3,j),i)
     ivt4 = kvert(NeighA(4,j),i)
     CALL GET_NormalArea(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dcorvg(1:3,nvt+net+nat+i),dAN)
     dV = dAN(1)*QuadSc%ValU(nvt+net+k) + dAN(2)*QuadSc%ValV(nvt+net+k) + dAN(3)*QuadSc%ValW(nvt+net+k)
!      write(*,'(4ES12.4)') dAN
     dVolFlow = dVolFlow + dV
    END IF
    k = k + 1
   END IF
  END DO
 END DO
end if

if (iPar.eq.1) then
 k=1
 DO i=1,nel
  DO j=1,6
   IF (k.eq.karea(j,i)) THEN
    IF (myBoundary%bOutflow(nvt+net+k)) THEN
     ivt1 = kvert(NeighA(1,j),i)
     ivt2 = kvert(NeighA(2,j),i)
     ivt3 = kvert(NeighA(3,j),i)
     ivt4 = kvert(NeighA(4,j),i)
     CALL GET_NormalArea(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dcorvg(1:3,nvt+net+nat+i),dAN)
     dAN = - dAN
     dV = dAN(1)*QuadSc%ValU(nvt+net+k) + dAN(2)*QuadSc%ValV(nvt+net+k) + dAN(3)*QuadSc%ValW(nvt+net+k)
     dVolFlow = dVolFlow + dV
    END IF
    k = k + 1
   END IF
  END DO
 END DO
end if

END SUBROUTINE Integrate_DIE_Flowrate
!
! ----------------------------------------------
!
SUBROUTINE IntegrateFlowrate(dcorvg,karea,kvert,nel,dVolFlow,dPar)
REAL*8 dcorvg(3,*),dVolFlow
INTEGER karea(6,*),kvert(8,*),nel
REAL*8 dPar
!---------------------------------
INTEGER NeighA(4,6)
REAL*8 P(3),dA
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4

dVolFlow = 0d0

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   IF (abs(dcorvg(3,ivt1)-dPar).lt.1d-4.and.abs(dcorvg(3,ivt2)-dPar).lt.1d-4.and. &
       abs(dcorvg(3,ivt3)-dPar).lt.1d-4.and.abs(dcorvg(3,ivt4)-dPar).lt.1d-4) then
       CALL GET_area(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dA)
       dVolFlow = dVolFlow + dA*(QuadSc%ValW(nvt+net+k))
   END IF
   k = k + 1
  END IF
 END DO
END DO

END SUBROUTINE IntegrateFlowrate
!
! ----------------------------------------------
!
SUBROUTINE IntegratePressure(dcorvg,karea,kvert,nel,dIntPres,dArea,dPar)
REAL*8 dcorvg(3,*),dIntPres
INTEGER karea(6,*),kvert(8,*),nel
REAL*8 dPar,dArea
!---------------------------------
INTEGER NeighA(4,6)
REAL*8 P(3),dA
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4

dIntPres = 0d0
dArea   = 0d0

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   IF (abs(dcorvg(3,ivt1)-dPar).lt.1d-4.and.abs(dcorvg(3,ivt2)-dPar).lt.1d-4.and. &
       abs(dcorvg(3,ivt3)-dPar).lt.1d-4.and.abs(dcorvg(3,ivt4)-dPar).lt.1d-4) then
       CALL GET_area(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dA)
       dIntPres = dIntPres + dA*(LinSc%ValP(NLMAX)%x(4*(i-1)+1))
       dArea = dArea + dA
   END IF
   k = k + 1
  END IF
 END DO
END DO

END SUBROUTINE IntegratePressure
!
! ----------------------------------------------
!
SUBROUTINE  AddViscoStress()
  EXTERNAL E013

  ILEV=NLMAX
  CALL SETLEV(2)

  CALL AssembleViscoStress(QuadSc%defU,QuadSc%defV,QuadSc%defW,BndrForce,&
    ViscoSc%Val11,ViscoSc%Val22,ViscoSc%Val33,ViscoSc%Val12,ViscoSc%Val13,ViscoSc%Val23,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%karea,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%dcorvg,&
    Properties%Viscosity(2),tstep,Properties%ViscoLambda,ViscoElasticForce,E013)

  CALL Comm_SummN(ViscoElasticForce,3)

END SUBROUTINE  AddViscoStress
!
! ----------------------------------------------
!
SUBROUTINE  STORE_OLD_MESH(dcoor)
  REAL*8 dcoor(3,*) 
  INTEGER i

  !write(*,*)'ndof:',QuadSc%ndof,mg_mesh%level(NLMAX+1)%nvt
  DO i=1,QuadSc%ndof
  myALE%OldCoor(:,i) = dcoor(:,i)
  END DO

END SUBROUTINE  STORE_OLD_MESH
!
! ----------------------------------------------
!
SUBROUTINE  STORE_NEW_MESH(dcoor)
  REAL*8 dcoor(3,*) 
  INTEGER i

  DO i=1,QuadSc%ndof
  myALE%NewCoor(:,i) = dcoor(:,i)
  END DO

END SUBROUTINE  STORE_NEW_MESH
!
! ----------------------------------------------
!
SUBROUTINE  GET_MESH_VELO()
  INTEGER i
  REAL*8 dmax,daux

  dmax = 0d0
  DO i=1,QuadSc%ndof
  myALE%MeshVelo(:,i) = (myALE%NewCoor(:,i)-myALE%OldCoor(:,i))/tstep
  ! myALE%MeshVelo(:,i) = (myALE%NewCoor(:,i)-myALE%OldCoor(:,i))/tstep
  daux = SQRT(myALE%MeshVelo(1,i)**2d0 + myALE%MeshVelo(2,i)**2d0 +&
    myALE%MeshVelo(3,i)**2d0)
  IF (daux.gt.dmax) dmax = daux
  END DO

  CALL COMM_Maximum(dmax)
  IF (myid.eq.showid) WRITE(*,*) 'max mesh velocity: ', dmax

END SUBROUTINE  GET_MESH_VELO
!
! ----------------------------------------------
!
SUBROUTINE  RotateMyMesh(dcoor)
  REAL*8 dcoor(3,*) 
  REAL*8 X,Y,Z,R,A,dAlpha
  REAL*8 :: PI = 3.141592654d0
  INTEGER i,nnn

  dAlpha = 1d0*(250d0/60d0)*2d0*timens*PI
  ! dAlpha = 2d0*(250d0/60d0)*2d0*tstep*PI

  IF (myid.eq.0) THEN
    nnn = KNVT(NLMAX)
  ELSE
    nnn = KNVT(NLMAX+1)
  END IF

  DO i=1,nnn
  !  X = dcoor(1,i)
  !  Y = dcoor(2,i)
  X = myALE%OrigCoor(1,i)
  Y = myALE%OrigCoor(2,i)
  R = DSQRT(X*X + Y*Y)
  A = DATAN(Y/X)
  IF (X.LT.0d0) A = A + PI
  A = A + dAlpha
  dcoor(1,i) = R*DCOS(A)
  dcoor(2,i) = R*DSIN(A)
  !IF (myid.eq.1) WRITE(*,'(2(2E12.4,A))') X,Y,' : ', dcoor(1:2,i)
  END DO

  ILEV = NLMIN
  CALL SETLEV(2)

  CALL ExchangeNodeValuesOnCoarseLevel(DWORK(L(LCORVG)),KWORK(L(LVERT)),NVT,NEL)

  ILEV = NLMAX
  CALL SETLEV(2)

END SUBROUTINE  RotateMyMesh
!
! ----------------------------------------------
!
SUBROUTINE  StoreOrigCoor(dcoor)
  INTEGER i,nnn
  REAL*8 dcoor(3,*) 

  IF (myid.eq.0) THEN
    nnn = KNVT(NLMAX)
  ELSE
    nnn = KNVT(NLMAX+1)
  END IF

  DO i=1,nnn
  myALE%OrigCoor(:,i) = dcoor(:,i)
  END DO

END SUBROUTINE  StoreOrigCoor
!
! ----------------------------------------------
!
SUBROUTINE GetMeshVelocity2(mfile)
  integer mfile
  REAL*8 dMaxVelo,daux
  INTEGER i

  dMaxVelo = 0d0
  IF (myid.ne.0) then
    DO i=1,QuadSc%ndof
    myALE%MeshVelo(:,i) = (myQ2Coor(:,i) -  myALE%Q2Coor_old(:,i))/tstep
    daux = myALE%MeshVelo(1,i)**2d0+myALE%MeshVelo(2,i)**2d0+myALE%MeshVelo(3,i)**2d0
    IF (dMaxVelo.lt.daux) dMaxVelo = daux
    END DO
  END IF

  CALL COMM_Maximum(dMaxVelo)

  IF (myid.eq.1) THEN
    WRITE(mfile,*)  "Maximum Mesh Velocity: ", SQRT(dMaxVelo)
    WRITE(mterm,*)  "Maximum Mesh Velocity: ", SQRT(dMaxVelo)
  END IF

END SUBROUTINE GetMeshVelocity2
!
! ----------------------------------------------
!
SUBROUTINE MoveInterfacePoints(dcoor,MFILE)
  INTEGER mfile
  REAL*8 dcoor(3,*)
  REAL*8 Velo(3),Displacement(3),dMaxVelo,daux,dArea
  INTEGER i,iInterface

  write(*,*)'New untested subroutine'
  stop

  IF (myid.ne.0) then
    DO i=1,QuadSc%ndof
    Velo = [QuadSc%ValU(i),QuadSc%ValV(i),QuadSc%ValW(i)]
    Displacement = 1d0*tstep*Velo
    myQ2Coor(:,i) = myQ2Coor(:,i) + Displacement
    END DO
  END IF

  IF (.NOT.ALLOCATED(myTSurf)) ALLOCATE(myTSurf(Properties%nInterface))

  IF (myid.ne.0) then
    ILEV=NLMAX
    CALL SETLEV(2)
    DO i=1,QuadSc%ndof
    dcoor(:,i) = myQ2Coor(:,i)
    END DO

    DO iInterface=1,Properties%nInterface
    CALL BuildUpTriangulation(KWORK(L(LVERT)),KWORK(L(LEDGE)),KWORK(L(LAREA)),myQ2Coor,iInterface)
    END DO
  END IF

  CALL CommunicateSurface()

  IF (myid.ne.0) then
    ILEV=NLMAX
    CALL SETLEV(2)
    DO i=1,QuadSc%ndof
    dcoor(:,i) = myALE%Q2coor_old(:,i)
    myQ2Coor(:,i) = myALE%Q2coor_old(:,i)
    END DO
  END IF

  ! IF (myid.eq.1) THEN
  !  CALL GetCompleteArea(dArea,1)
  !  WRITE(MFILE,'(A,3ES14.6)') "CompleteSurfaceAreaAndCircularity: ", TIMENS,dArea, (0.05*(2*(4d0*ATAN(1d0))*0.25d0))/dArea
  !  WRITE(MTERM,'(A,3ES14.6)') "CompleteSurfaceAreaAndCircularity: ", TIMENS,dArea, (0.05*(2*(4d0*ATAN(1d0))*0.25d0))/dArea
  ! END IF

END SUBROUTINE MoveInterfacePoints
!
! ----------------------------------------------
!
SUBROUTINE GetCompleteArea(DCompleteArea,iIF)
  REAL*8 DCompleteArea
  REAL*8 DA(3,3),dArea
  INTEGER i,j,IP1,IP2,IP3,iIF

  write(*,*)'New untested subroutine'
  stop

  DCompleteArea = 0d0
  DO i=1,myTSurf(iIF)%nT
  DO j=1,8
  IP1 = j
  IP2 = MOD(j,8)+1
  IP3 = 9
  DA(1,2)=myTSurf(iIF)%T(i)%C(1,IP2) - myTSurf(iIF)%T(i)%C(1,IP1) !P2X-P1X
  DA(2,2)=myTSurf(iIF)%T(i)%C(2,IP2) - myTSurf(iIF)%T(i)%C(2,IP1) !P2Y-P1Y
  DA(3,2)=myTSurf(iIF)%T(i)%C(3,IP2) - myTSurf(iIF)%T(i)%C(3,IP1) !P2Z-P1Z
  DA(1,3)=myTSurf(iIF)%T(i)%C(1,IP3) - myTSurf(iIF)%T(i)%C(1,IP1) !P3X-P1X
  DA(2,3)=myTSurf(iIF)%T(i)%C(2,IP3) - myTSurf(iIF)%T(i)%C(2,IP1) !P3Y-P1Y
  DA(3,3)=myTSurf(iIF)%T(i)%C(3,IP3) - myTSurf(iIF)%T(i)%C(3,IP1) !P3Z-P1Z

  DA(1,1) = DA(2,3)*DA(3,2) - DA(3,3)*DA(2,2)
  DA(2,1) = DA(3,3)*DA(1,2) - DA(1,3)*DA(3,2)
  DA(3,1) = DA(1,3)*DA(2,2) - DA(2,3)*DA(1,2)
  dArea = 0.5*SQRT(DA(1,1)*DA(1,1) + DA(2,1)*DA(2,1) + DA(3,1)*DA(3,1))
  DCompleteArea = DCompleteArea + dArea
  END DO
  END DO

END SUBROUTINE GetCompleteArea
!
! ----------------------------------------------
!
subroutine BuildUpTriangulation(kvert,kedge,karea,dcorvg,iIF)
  REAL*8 dcorvg(3,*)
  INTEGER karea(6,*),kvert(8,*),kedge(12,*),iIF
  INTEGER iel,i,j,k,ivt1,ivt2,ivt3,ivt4,ivt5,iT
  INTEGER NeighA(4,6),NeighU(4,6)
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
  DATA NeighU/1,2,3,4, 1,6,9,5, 2,7,10,6, 3,8,11,7, 4,5,12,8, 9,10,11,12/

  write(*,*)'New untested subroutine'
  stop

  iT = 0
  k=1
  DO i=1,nel
  DO j=1,6
  IF (k.eq.karea(j,i)) THEN
    IF (myBoundary%LS_zero(nvt+net+k).eq.iIF) THEN
      iT = iT + 1
    END IF
    k = k + 1
  END IF
  END DO
  END DO

  IF (ALLOCATED(myTSurf(iIF)%T)) DEALLOCATE(myTSurf(iIF)%T)

  myTSurf(iIF)%nT = iT
  ALLOCATE(myTSurf(iIF)%T(myTSurf(iIF)%nT))

  iT = 0
  k=1
  DO i=1,nel
  DO j=1,6
  IF (k.eq.karea(j,i)) THEN
    IF (myBoundary%LS_zero(nvt+net+k).eq.iIF) THEN
      iT = iT + 1
      ivt1 = kvert(NeighA(1,j),i)
      ivt2 = kvert(NeighA(2,j),i)
      ivt3 = kvert(NeighA(3,j),i)
      ivt4 = kvert(NeighA(4,j),i)
      myTSurf(iIF)%T(iT)%C(:,1) = dcorvg(:,ivt1)
      myTSurf(iIF)%T(iT)%C(:,3) = dcorvg(:,ivt2)
      myTSurf(iIF)%T(iT)%C(:,5) = dcorvg(:,ivt3)
      myTSurf(iIF)%T(iT)%C(:,7) = dcorvg(:,ivt4)
      ivt1 = kedge(NeighU(1,j),i)
      ivt2 = kedge(NeighU(2,j),i)
      ivt3 = kedge(NeighU(3,j),i)
      ivt4 = kedge(NeighU(4,j),i)
      myTSurf(iIF)%T(iT)%C(:,2) = dcorvg(:,nvt+ivt1)
      myTSurf(iIF)%T(iT)%C(:,4) = dcorvg(:,nvt+ivt2)
      myTSurf(iIF)%T(iT)%C(:,6) = dcorvg(:,nvt+ivt3)
      myTSurf(iIF)%T(iT)%C(:,8) = dcorvg(:,nvt+ivt4)
      myTSurf(iIF)%T(iT)%C(:,9) = dcorvg(:,nvt+net+k)
    END IF
    k = k + 1
  END IF
  END DO
  END DO

end subroutine BuildUpTriangulation
!
! ----------------------------------------------
!
subroutine IntegrateQuantities(mfile)
  INTEGER mfile
  REAL*8 dArray(8)
  REAL*8 :: dR=0.25d0, myPI=4d0*ATAN(1d0), dWidth=0.05d0
  EXTERNAL E013

  write(*,*)'New untested subroutine'
  stop

  IF (myid.ne.0) THEN
    dArray = 0d0
    ILEV = NLMAX
    CALL SETLEV(2)
    CALL GetSurface(KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),&
      myQ2coor,E013,dArray(1))
    CALL GetVolume(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),&
      myQ2coor,E013,dArray(2),dArray(3:5),dArray(6:8))
  END IF

  CALL Comm_SummN(dArray,8)

  IF (myid.eq.1) THEN
    WRITE(MTERM,'(A)') "Time Circularity Mass Center RiseVelo RefFrame"
    WRITE(MTERM,'(A,ES12.4,2ES14.6,10ES12.4)') "Stats: ", timens,&
      (dWidth*(2*myPI*dR))/dArray(1),(dWidth*(myPI*dR*dR))/dArray(2),dArray(3:5)/dArray(2),dArray(6:8)/dArray(2),myALE%dFrameVelocity(2)
    WRITE(MFILE,'(A)') "Time Circularity Mass Center RiseVelo RefFrame"
    WRITE(MFILE,'(A,ES12.4,2ES14.6,10ES12.4)') "Stats: ", timens,&
      (dWidth*(2*myPI*dR))/dArray(1),(dWidth*(myPI*dR*dR))/dArray(2),dArray(3:5)/dArray(2),dArray(6:8)/dArray(2),myALE%dFrameVelocity(2)
  END IF

  IF (myALE%bUseFrameVelocity) THEN
    myALE%dFrameVelocityChange = dArray(6:8)/dArray(2)
    myALE%dFrameVelocity = myALE%dFrameVelocity + 1d0*myALE%dFrameVelocityChange
  ELSE
    myALE%dFrameVelocityChange = 0d0
    myALE%dFrameVelocity       = 0d0
  END IF

  IF (myid.eq.1) THEN
    WRITE(MTERM,'(A,3ES12.4)') "ReferenceFrame: ", myALE%dFrameVelocity(2), myALE%dFrameVelocityChange(2)/TSTEP
    WRITE(MFILE,'(A,3ES12.4)') "ReferenceFrame: ", myALE%dFrameVelocity(2), myALE%dFrameVelocityChange(2)/TSTEP
  END IF
end subroutine IntegrateQuantities
!
! ----------------------------------------------
!
subroutine InitMeshDeform(mfile, mgMesh)
USE var_QuadScalar, ONLY : nMainUmbrellaSteps,tMultiMesh
USE var_QuadScalar, ONLY : nUmbrellaStepsLvl
use umbrella_smoother, only : us_UmbrellaSmoother
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
USE def_FEAT

implicit none

integer, intent(in) :: mfile

type(tMultiMesh), intent(inout) :: mgMesh

! local variables
real :: dttt0,dttt1
integer :: i

CALL myMPI_Barrier()
CALL ztime(dttt0)

ilev=nlmax
call setlev(2)

if (nMainUmbrellaSteps.ne.0) then
 nMainUmbrellaSteps =nMainUmbrellaSteps + (mySetup%MeshResolution-1)
end if

do i=1,nMainUmbrellaSteps

 call us_UmbrellaSmoother(0d0, nUmbrellaStepsLvl, mgMesh, QuadSc)

 ilev=nlmax
 call setlev(2)

 call SetUp_myQ2Coor( mgMesh%level(ilev)%dcorvg,&
                      mgMesh%level(ilev)%dcorag,&
                      mgMesh%level(ilev)%kvert,&
                      mgMesh%level(ilev)%karea,&
                      mgMesh%level(ilev)%kedge)

END DO

call myMPI_Barrier()
call ztime(dttt1)
if (myid.eq.1) write(mfile,"(A,F6.1,A)") "Time used for mesh smoothening was: ", dttt1-dttt0, "s!"
if (myid.eq.1) write(*,"(A,F6.1,A)") "Time used for mesh smoothening was: ", dttt1-dttt0, "s!"

end subroutine InitMeshDeform
!
! ----------------------------------------------
!
subroutine InitOperators(mfile, mgMesh,bCreate)
use PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
use var_QuadScalar, only : tMultiMesh
implicit none

integer, intent(in) :: mfile
logical :: bCreate
type(tMultiMesh), intent(inout) :: mgMesh

! local variables
integer :: i

ilev = mgMesh%nlmax
call setlev(2)

call SetUp_myQ2Coor( mgMesh%level(ilev)%dcorvg,&
                     mgMesh%level(ilev)%dcorag,&
                     mgMesh%level(ilev)%kvert,&
                     mgMesh%level(ilev)%karea,&
                     mgMesh%level(ilev)%kedge)

call StoreOrigCoor(mgMesh%level(mgMesh%nlmax)%dcorvg)
call store_old_mesh(mgMesh%level(mgMesh%nlmax)%dcorvg)

ilev = mgMesh%nlmax
call setlev(2)

if (myid.ne.0) call updateMixerGeometry(mfile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       PRESS BC        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (mySetup%bPressureFBM) THEN
 ilev=nlmin
 CALL SETLEV(2)
 CALL SetPressBC_NewGen(mgMesh)
 ! send them to the master
 ilev=nlmin
 CALL SETLEV(2)
 CALL SendPressBCElemsToCoarse(LinSc%knprP(ilev)%x,nel)
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SET BC !!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if (myid.ne.0) then
!   ilev=nlmin+1
!   CALL SETLEV(2)
!  END IF
!  CALL SetPressBC(mgMesh)
! 
 do ilev=nlmin+1,nlmax
  CALL SETLEV(2)
  CALL GetMG_KNPRP(mgMesh)
 end do
! 
 ! Set up the boundary condition types (knpr)
 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)
  CALL IncludeFBM_BCs(mgMesh)
  CALL QuadScalar_Knpr()
 END DO
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       PRESS BC        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (bCreate) THEN
 CALL Release_cgal_structures()
 call OperatorRegenaration(1)
 call OperatorRegenaration(2)
 call OperatorRegenaration(3)

 CALL Create_MMat()
END IF

end subroutine InitOperators
!
! ----------------------------------------------
!
SUBROUTINE SetPressBC_NewGen(mgMesh)
type(tMultiMesh), intent(inout) :: mgMesh
integer iel,jlev
real*8 dnn
logical bKick

ilev=nlmin
dnn=0d0

if (myid.ne.0) then
 DO iel=1,nel

  jlev=nlmin
  
  bKick = .true.
  CALL FindPressBC_NewGenREC(iel,jlev)
  
  if (bKick) then
   dnn = dnn + 1d0 
   LinSc%knprP(jlev)%x(iel) = 1
   CALL SetPressBC_NewGenREC(iel,jlev)
  end if
  
 END DO
 
END IF

call Comm_Summ(dnn)

if (Myid.eq.showid) write(*,*) 'Number of Pressure BC Elements:',int(dnn)

 CONTAINS
 
RECURSIVE SUBROUTINE FindPressBC_NewGenREC(iiel,iilev)
integer iiel,iilev

integer JEL(8)
integer jjlev,i,ivt
real*8 dSize

 if (iilev.eq.mgMesh%nlmax+1) then
  dSize = 0.1d0*0.1d0*(mgMesh%level(iilev)%dvol(iiel)**(1d0/3d0)) ! 0.1 because of the cm --> mm scaling
 else
  dSize = 0d0
 end if
 
 do i=1,8
  ivt = mgMesh%level(iilev)%kvert(i,iiel)
  if (screw(ivt).gt.-dSize.and.shell(ivt).gt.-dSize) THEN
   bKick = .false.
   RETURN
  end if
 end do

 if (iilev+1.gt.mgMesh%nlmax) RETURN
! Possible canditate found
 jjlev = iilev + 1

 JEL(1)  = iiel
 JEL(2)  = mgMesh%level(jjlev)%kadj(3,JEL(1))
 JEL(3)  = mgMesh%level(jjlev)%kadj(3,JEL(2))
 JEL(4)  = mgMesh%level(jjlev)%kadj(3,JEL(3))
 JEL(5)  = mgMesh%level(jjlev)%kadj(6,JEL(1))
 JEL(6)  = mgMesh%level(jjlev)%kadj(6,JEL(2))
 JEL(7)  = mgMesh%level(jjlev)%kadj(6,JEL(3))
 JEL(8)  = mgMesh%level(jjlev)%kadj(6,JEL(4))
 
 DO i=1,8
  CALL FindPressBC_NewGenREC(JEL(i),jjlev)
  if (.not.bKick) RETURN
 end do

END SUBROUTINE FindPressBC_NewGenREC

RECURSIVE SUBROUTINE SetPressBC_NewGenREC(iiel,iilev)
integer iiel,iilev

integer JEL(8)
integer jjlev,i,ivt

LinSc%knprP(iilev)%x(iiel) = 1
! write(*,*) iilev,iiel

if ((iilev + 1).gt.mgMesh%nlmax) RETURN
! Possible canditate found
 jjlev = iilev + 1

 JEL(1)  = iiel
 JEL(2)  = mgMesh%level(jjlev)%kadj(3,JEL(1))
 JEL(3)  = mgMesh%level(jjlev)%kadj(3,JEL(2))
 JEL(4)  = mgMesh%level(jjlev)%kadj(3,JEL(3))
 JEL(5)  = mgMesh%level(jjlev)%kadj(6,JEL(1))
 JEL(6)  = mgMesh%level(jjlev)%kadj(6,JEL(2))
 JEL(7)  = mgMesh%level(jjlev)%kadj(6,JEL(3))
 JEL(8)  = mgMesh%level(jjlev)%kadj(6,JEL(4))

 DO i=1,8
  LinSc%knprP(jjlev)%x(JEL(i)) = 1
  CALL SetPressBC_NewGenREC(JEL(i),jjlev)
 end do

END SUBROUTINE SetPressBC_NewGenREC

END SUBROUTINE SetPressBC_NewGen
!
! ----------------------------------------------
!
SUBROUTINE SetPressBC(mgMesh)
type(tMultiMesh), intent(inout) :: mgMesh
integer i,iel,ivt
real*8 dnn

dnn=0d0
if (myid.ne.0) then
 DO iel=1,nel
  do i=1,8
   ivt = mgMesh%level(ilev)%kvert(i,iel)
   if (screw(ivt).gt.0d0.and.shell(ivt).gt.0d0) goto 1
  end do
  dnn = dnn + 1d0 
!   write(*,*) 'Pressure BC !!!',myid,iel
  LinSc%knprP(ilev)%x(iel) = 1
 1 continue
 END DO
END IF

call Comm_Summ(dnn)

if (Myid.eq.showid) write(*,*) 'Number of Pressure BC Elements:',int(dnn)

END SUBROUTINE SetPressBC
!
! ----------------------------------------------
!
SUBROUTINE GetMG_KNPRP(mgMesh)
type(tMultiMesh), intent(inout) :: mgMesh
integer iel,jel(8),i,jjj

jjj = 0
do iel = 1,nel/8
 if (LinSc%knprP(ilev-1)%x(iel).eq.1) then
  JEL(1)  = iel
  JEL(2)  = mgMesh%level(ilev)%kadj(3,JEL(1))
  JEL(3)  = mgMesh%level(ilev)%kadj(3,JEL(2))
  JEL(4)  = mgMesh%level(ilev)%kadj(3,JEL(3))
  JEL(5)  = mgMesh%level(ilev)%kadj(6,JEL(1))
  JEL(6)  = mgMesh%level(ilev)%kadj(6,JEL(2))
  JEL(7)  = mgMesh%level(ilev)%kadj(6,JEL(3))
  JEL(8)  = mgMesh%level(ilev)%kadj(6,JEL(4))
  do i=1,8
    LinSc%knprP(ilev)%x(jel(i)) = 1
    jjj = jjj + 1
  end do
 end if
end do

if (myid.eq.1) write(*,'(A,I0,A,I0)') 'KNPRP nodes on level ',ilev, ' are :', jjj

END SUBROUTINE GetMG_KNPRP
!
! ----------------------------------------------
!
SUBROUTINE IncludeFBM_BCs(mgMesh)

type(tMultiMesh), intent(inout) :: mgMesh
integer iel,i

do iel = 1,nel
 if (LinSc%knprP(ilev)%x(iel).eq.1) then
  do i=1,8
   myBoundary%bWall(mgMesh%level(ilev)%kvert(i,iel)) = .true.
  end do

  do i=1,12
   myBoundary%bWall(nvt + mgMesh%level(ilev)%kedge(i,iel)) = .true.
  end do

  do i=1,6
   myBoundary%bWall(nvt + net + mgMesh%level(ilev)%karea(i,iel)) = .true.
  end do

  myBoundary%bWall(nvt + net + nat + iel) = .true.
 end if
end do

QuadSc%auxU = 0d0
DO i=1,QuadSc%ndof
 if (myBoundary%bWall(i)) QuadSc%auxU(i) = 1d0
END DO

CALL E013Sum(QuadSc%auxU)

DO i=1,QuadSc%ndof
 if (QuadSc%auxU(i).gt.0d0) myBoundary%bWall(i) = .true.
END DO

!!!! COMMUNICATION needed !!!!!


END SUBROUTINE IncludeFBM_BCs
!
! ----------------------------------------------
!
SUBROUTINE  Setup_PeriodicVelocityRHS()
INTEGER i

IF (.NOT.(ALLOCATED(dPeriodicVector))) ALLOCATE(dPeriodicVector(QuadSc%ndof))

IF (ieee_is_finite(myProcess%dPress)) THEN
!IF (.NOT.ISNAN(myProcess%dPress)) THEN

DO i=1,SIZE(LinSc%AuxP(NLMAX)%x)
 IF (MOD(i,4).EQ.1) then
  LinSc%AuxP(NLMAX)%x(i) = myProcess%dPress
 ELSE
  LinSc%AuxP(NLMAX)%x(i) = 0d0
 END IF
END DO

CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%auxP(NLMAX)%x,&
     QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof,+1d0,0d0)

DO I=1,QuadSc%ndof
 IF ((myQ2Coor(3,i)     .LT.+1e-3)) THEN
  dPeriodicVector(i) = QuadSc%auxW(i)
 ELSE
  dPeriodicVector(i) = 0d0
 END IF
END DO
ELSE
  dPeriodicVector = 0d0
END IF

END SUBROUTINE  Setup_PeriodicVelocityRHS
!
! ----------------------------------------------
!
SUBROUTINE RestrictWallBC()
INTEGER ndof
INTEGER i,j,k
integer iat,ivt1,ivt2,ivt3,ivt4
INTEGER NeighA(4,6),NeighU(4,6)
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
DATA NeighU/1,2,3,4, 1,6,9,5, 2,7,10,6, 3,8,11,7, 4,5,12,8, 9,10,11,12/

 ndof  = mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel
 
 nvt  = mg_mesh%level(ilev)%nvt
 net  = mg_mesh%level(ilev)%net
 nat  = mg_mesh%level(ilev)%nat
 nel  = mg_mesh%level(ilev)%nel
 
 ! overlap outper points with wall ==> inner 'wall' markers will disappear
 DO i=1,ndof
  IF (myBoundary%bWall(i)) THEN
   IF (.not.(mg_mesh%BndryNodes(i)%bOuterPoint)) myBoundary%bWall(i) = .FALSE.
  END IF
 END DO 
 
!  remove all outflow / inflow dofs
  k=1
  DO i=1,nel
   DO j=1,6
    IF (k.eq.mg_mesh%level(ilev)%karea(j,i)) THEN
     if (myBoundary%iInflow(nvt+net+k).ne.0.or.myBoundary%bOutflow(nvt+net+k)) then
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(1,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(2,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(3,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(4,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(1,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(2,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(3,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(4,j),i)) = .false.
      myBoundary%bWall(nvt+net+k) = .false.
     end if
     k = k + 1
    END IF
   END DO
  END DO
 
!  remove all symmetry dofs
  k=1
  DO i=1,nel
   DO j=1,6
    IF (k.eq.mg_mesh%level(ilev)%karea(j,i)) THEN
     if (myBoundary%bSymmetry(1,nvt+net+k).or.myBoundary%bSymmetry(2,nvt+net+k).or.myBoundary%bSymmetry(3,nvt+net+k)) then
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(1,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(2,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(3,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(4,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(1,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(2,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(3,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(4,j),i)) = .false.
      myBoundary%bWall(nvt+net+k) = .false.
     end if
     k = k + 1
    END IF
   END DO
  END DO

  ! refresh Wall
  k=1
  DO i=1,nel
   DO j=1,6
    IF (k.eq.mg_mesh%level(ilev)%karea(j,i)) THEN
     if (myBoundary%bWall(nvt+net+k)) then
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(1,j),i)) = .true.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(2,j),i)) = .true.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(3,j),i)) = .true.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(4,j),i)) = .true.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(1,j),i)) = .true.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(2,j),i)) = .true.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(3,j),i)) = .true.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(4,j),i)) = .true.
     end if
     k = k + 1
    END IF
   END DO
  END DO
   
END SUBROUTINE RestrictWallBC
!
! ----------------------------------------------
!
SUBROUTINE UpdateMaterialProperties
integer ii,iMat,jel(8),ndof_nel,iMaxFrac
REAL*8 dMaxFrac
REAL*8, allocatable :: dFrac(:)


if (myid.eq.1) WRITE(*,*) 'updating Material Properties Distribution... '
IF (myid.ne.0) then

!     do iel = 1,mg_mesh%level(nlmax)%nel
!      ii = mg_mesh%level(nlmax)%nvt + mg_mesh%level(nlmax)%net + mg_mesh%level(nlmax)%nat + iel
!      if (mg_mesh%level(nlmax)%dcorvg(3,ii).gt.92d0) then
!       MaterialDistribution(nlmax)%x(iel) = 2
!      end if
!     end do
     
!  write(*,*) MaterialDistribution(NLMAX+0-1)%x(1:knel(NLMAX+0-1))
! pause

 IF (allocated(MaterialDistribution)) then
 
  if (istart.eq.2) then
    do iel = 1,mg_mesh%level(nlmax-1)%nel
      jel(1) = iel
      jel(2) = mg_mesh%level(ilev+1)%kadj(3,jel(1))
      jel(3) = mg_mesh%level(ilev+1)%kadj(3,jel(2))
      jel(4) = mg_mesh%level(ilev+1)%kadj(3,jel(3))
      jel(5) = mg_mesh%level(ilev+1)%kadj(6,jel(1))
      jel(6) = mg_mesh%level(ilev+1)%kadj(6,jel(2))
      jel(7) = mg_mesh%level(ilev+1)%kadj(6,jel(3))
      jel(8) = mg_mesh%level(ilev+1)%kadj(6,jel(4))
   
      iMat = MaterialDistribution(nlmax-1)%x(jel(1))
!       write(*,*) iMAt
      do ii=1,8
       MaterialDistribution(nlmax)%x(jel(1)) = iMat
      end do
     end do
  end if
  
  ndof_nel = (knvt(NLMAX) + knat(NLMAX) + knet(NLMAX))
  allocate(dFrac(myMultiMat%nOfMaterials))

  DO ilev=NLMAX-1,NLMIN,-1
   ndof_nel = (knvt(ilev+1) + knat(ilev+1) + knet(ilev+1))
   do iel = 1,mg_mesh%level(ilev)%nel
     jel(1) = iel
     jel(2) = mg_mesh%level(ilev+1)%kadj(3,jel(1))
     jel(3) = mg_mesh%level(ilev+1)%kadj(3,jel(2))
     jel(4) = mg_mesh%level(ilev+1)%kadj(3,jel(3))
     jel(5) = mg_mesh%level(ilev+1)%kadj(6,jel(1))
     jel(6) = mg_mesh%level(ilev+1)%kadj(6,jel(2))
     jel(7) = mg_mesh%level(ilev+1)%kadj(6,jel(3))
     jel(8) = mg_mesh%level(ilev+1)%kadj(6,jel(4))

!      write(*,*) jel
     dFrac = 0d0
     do ii=1,8
      iMat = MaterialDistribution(ilev+1)%x(jel(ii))
      dFrac(iMat) = dFrac(iMat) + 1d0 !mg_mesh%level(ilev+1)%dvol(jel(ii))
     end do
     dMaxFrac = 0d0
!      write(*,*) dFrac
     do ii=1,myMultiMat%nOfMaterials
      IF (dFrac(ii).gt.dMaxFrac) THEN
       dMaxFrac = dFrac(ii)
       iMaxFrac = ii
      END IF
     end do
     
     MaterialDistribution(ilev)%x(jel(1)) = iMaxFrac
     
    end do
  END DO
  
  deallocate(dFrac)
 end if

else 
 DO ilev=NLMAX,NLMIN,-1
  MaterialDistribution(ilev)%x = myMultiMat%InitMaterial
 end do
! write(*,*) MaterialDistribution(nlmax)%x
end if

! write(*,*) MaterialDistribution(NLMAX)%x(1:knel(NLMAX))
! pause
END SUBROUTINE UpdateMaterialProperties
!
! ----------------------------------------------
!
SUBROUTINE DetermineIfGoalsWereReached(bGoalsReached)
use, intrinsic :: ieee_arithmetic
REAL*8 myinf
LOGICAL bGoalsReached

if(ieee_support_inf(myInf))then
  myInf = ieee_value(myInf, ieee_negative_inf)
endif

bGoalsReached = .true.

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."XSE") THEN
 if (myProcess%FillingDegree.eq.myInf .or. myProcess%FillingDegree .eq. 1d0) then
    if (itns.ge.nitns) bGoalsReached=.false.
 end if
END IF
   
IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 if (istart.eq.1) then
   if (itns.ge.nitns) bGoalsReached=.false.
 end if
END IF


if (.not.bGoalsReached) THEN
 if (myid.eq.1) write(*,*) 'max time steps have been reached // the simulation has - most probably - not converged! '
end if

END SUBROUTINE DetermineIfGoalsWereReached

END MODULE Transport_Q2P1
