! TODO: new name Navier-Stokes + Q2 or similar
MODULE Transport_Q2P1  

USE def_QuadScalar
! USE PP3D_MPI
USE PP3D_MPI, ONLY:myid,master,E011Sum,COMM_Maximum,&
                   COMM_NLComplete,Comm_Summ,Comm_SummN,&
                   myMPI_Barrier
USE Parametrization,ONLY : InitBoundaryStructure,myParBndr,&
ParametrizeQ2Nodes

USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,GetMixerKnpr
! USE PP3D_MPI, ONLY:E011Sum,E011True_False,Comm_NLComplete,&
!               Comm_Maximum,Comm_Summ,knprmpi,myid,master
! USE LinScalar, ONLY: AddSurfaceTension
IMPLICIT NONE

TYPE(TQuadScalar), target :: QuadSc
! TODO: Move these to var
TYPE(TLinScalar)    LinSc
TYPE(TViscoScalar)  ViscoSc
TYPE(TParLinScalar) PLinSc
REAL*8, ALLOCATABLE :: ST_force(:)
REAL*8 :: Density_Secondary=1d0,Density_Primary=1d0
REAL*8 :: myPowerLawFluid(3),ViscoElasticForce(3)
REAL*8 :: Sigma=0.034D0,DiracEps=0.00625d0
INTEGER, ALLOCATABLE :: QuadScBoundary(:)
INTEGER PressureSample(2)
REAL tttt0,tttt1

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE Transport_Q2P1_UxyzP(mfile,inl_u,itns)

use cinterface, only: calculateDynamics,calculateFBM
INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit
INTEGER INLComplete,I,J,IERR,iOuter,iITER

 IF (calculateFBM()) THEN
  CALL updateFBMGeometry()
 END IF

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
CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,mfile)

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

IF (calculateDynamics()) THEN
 CALL FBM_GetForces()
 CALL updateFBM(Properties%Density(1),tstep,timens,Properties%Gravity,mfile,myid)
END IF

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

 IF (calculateFBM()) THEN
  CALL updateFBMGeometry()
 END IF

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

 QuadSc%cName = "Velo"
 LinSc%cName = "Pres"

 CALL GetVeloParameters(QuadSc%prm,QuadSc%cName,mfile)
 CALL GetPresParameters(LinSc%prm,LinSc%cName,mfile)

END SUBROUTINE Init_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE Init_QuadScalar_Structures_sse(mfile)
implicit none
LOGICAL bExist
INTEGER I,J,ndof,mfile,LevDif
integer :: mydof
integer :: maxlevel
Real*8 :: dabl

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
 Shearrate = 1d0

 Viscosity = Properties%Viscosity(1)

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

 CALL OperatorRegenaration(1)

END SUBROUTINE Init_QuadScalar_Structures_sse
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

 myPowerLawFluid(2) = 0.001d0
 myPowerLawFluid(3) = 0.75d0

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

 Viscosity = Properties%Viscosity(1)

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

  IF (myBoundary%bWall(i).OR.myBoundary%iInflow(i).GT.0) THEN
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
  REAL*8  dcorvg(3,*),dcorag(3,*)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
  REAL*8 PX,PY,PZ,DIST
  INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

  DO i=1,nvt
  PX = dcorvg(1,I)
  PY = dcorvg(2,I)
  PZ = dcorvg(3,I)
  CALL GetFictKnpr(PX,PY,PZ,QuadScBoundary(i),FictKNPR(i),Distance(i))
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
    CALL GetFictKnpr(PX,PY,PZ,QuadScBoundary(nvt+k),FictKNPR(nvt+k),Distance(nvt+k))
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
    CALL GetFictKnpr(PX,PY,PZ,QuadScBoundary(nvt+net+k),FictKNPR(nvt+net+k),Distance(nvt+net+k))
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
  CALL GetFictKnpr(PX,PY,PZ,QuadScBoundary(nvt+net+i),FictKNPR(nvt+net+nat+i),Distance(nvt+net+nat+i))
  END DO

  do i=1,nvt+net+nat+nel
  myALE%Monitor(i)=distance(i)
  end do


END SUBROUTINE QuadScalar_FictKnpr
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

  END DO

END SUBROUTINE Boundary_QuadScalar_Def
!
! ----------------------------------------------
!
SUBROUTINE Boundary_QuadScalar_Val()
  implicit none
  REAL*8 PX,PY,PZ
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

    IF (myBoundary%iInflow(i).GT.0) THEN 
      inpr = 1
      iType = myBoundary%iInflow(i)
      CALL GetVeloBCVal(PX,PY,PZ,QuadSc%valU(i),QuadSc%valV(i),QuadSc%valW(i),iType,timens)
    END IF
    finpr = FictKNPR(i)
    minpr = MixerKNPR(i)
    IF (finpr.ne.0.and.inpr.eq.0) THEN
      CALL GetVeloFictBCVal(PX,PY,PZ,QuadSc%valU(i),QuadSc%valV(i),QuadSc%valW(i),finpr,timens)
    END IF
    IF (minpr.ne.0.and.inpr.eq.0) THEN
      CALL GetVeloMixerVal(PX,PY,PZ,QuadSc%valU(i),QuadSc%valV(i),QuadSc%valW(i),minpr,timens)
    END IF
  END DO

END SUBROUTINE Boundary_QuadScalar_Val
!
! ----------------------------------------------
!
SUBROUTINE Boundary_QuadScalar_Mat(DA11,DA22,DA33,KLD,&
    KNPRU,KNPRV,KNPRW,NDOF)
  REAL*8  DA11(*),DA22(*),DA33(*)
  INTEGER KLD(*),KNPRU(*),KNPRV(*),KNPRW(*),ICOL,I,NDOF

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

  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%valP(NLMAX)%x,&
    QuadSc%defU,QuadSc%defV,QuadSc%defW,QuadSc%ndof,TSTEP,1d0)

END SUBROUTINE AddPressureGradient
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

    CALL Create_CMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW,LinSc%prm%MGprmIn%MinLev,LinSc%prm%MGprmIn%CrsSolverType)
    IF (myid.ne.master) THEN
      CALL Create_ParCMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW)
    END IF
    bHit = .TRUE.
  END IF

  IF (myid.EQ.ShowID.AND.bHit) WRITE(MTERM,'(A)', advance='yes') " "

END SUBROUTINE OperatorRegenaration
!
! ----------------------------------------------
!
SUBROUTINE ProlongateSolution()

  CALL ProlongateSolutionSub(QuadSc,LinSc,Boundary_QuadScalar_Val)
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
SUBROUTINE FAC_GetForces(mfile)
  INTEGER mfile
  REAL*8 :: Force(3),U_mean=0.2d0,H=0.205d0,D=0.1d0,Factor
!   REAL*8 :: Force(3),U_mean=0.2d0,H=0.05d0,D=0.1d0,Factor
  REAL*8 :: Sc_U = 1d0, Sc_Mu = 1d0, Sc_a = 1d0, PI = 3.141592654d0
  REAL*8 :: Force2(3),ForceV(3),ForceP(3)
  REAL*8 :: Scale
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
    Scale = 6d0*PI*Sc_Mu*Sc_U*Sc_a
    Force = (4d0*Force)/Scale
    ViscoElasticForce = (4d0*ViscoElasticForce)/Scale
  else
    Factor = 2d0/(U_mean*U_mean*D*H)
    ForceP = Factor*ForceP
    ForceV = Factor*ForceV
  end if

  IF (myid.eq.showID) THEN
   
    if(bViscoElastic)then
      WRITE(MTERM,5)
      WRITE(MFILE,5)
      write(mfile,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
        timens,ViscoElasticForce(3),(Force(3)-ViscoElasticForce(3)),Force(3)
      write(mterm,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
        timens,ViscoElasticForce(3),(Force(3)-ViscoElasticForce(3)),Force(3)
      WRITE(666,'(10ES13.5)')timens,ViscoElasticForce,&
        (Force-ViscoElasticForce),Force
!       WRITE(666,'(10G16.8)') Timens,Force
    else
      WRITE(MTERM,5)
      WRITE(MFILE,5)
      if (itns.eq.1) then
       write(mfile,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
       write(mterm,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
      end if
      write(mfile,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
      write(mterm,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
      WRITE(666,'(10ES16.8)') Timens,ForceV+forceP,ForceV,forceP
    end if
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
SUBROUTINE updateFBMGeometry()

  IF (myid.eq.showid) WRITE(*,*) '> FBM computation step'

  ILEV=NLMAX
  CALL SETLEV(2)
  CALL QuadScalar_FictKnpr(mg_mesh%level(ilev)%dcorvg,&
    mg_mesh%level(ilev)%dcorag,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%karea)

  CALL E013Max_SUPER(FictKNPR)
  if (myid.eq.1) write(*,*) 'CALL E013Max_SUPER(FictKNPR)'

END SUBROUTINE  updateFBMGeometry
!
! ----------------------------------------------
!
SUBROUTINE updateMixerGeometry(mfile)
use geometry_processing, only : calcDistanceFunction, QuadScalar_MixerKnpr,dEpsDist

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
 CALL calcDistanceFunction(mg_mesh%level(ilev)%dcorvg,&
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
                           QuadSc%AuxU,QuadSc%AuxV)
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
ilev = ilevel

 IF (myid.ne.0) then
  QuadSc%defU = 0d0
  QuadSc%defV = 0d0
  QuadSc%defW = 0d0
  CALL L2ProjVisco(QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                   QuadSc%defU,QuadSc%defV,QuadSc%defW,&
                   mg_mesh%level(ilevel)%kvert,&
                   mg_mesh%level(ilevel)%karea,&
                   mg_mesh%level(ilevel)%kedge,&
                   mg_mesh%level(ilevel)%dcorvg,E013)

  CALL E013Sum(QuadSc%defU)
  CALL E013Sum(QuadSc%defV)
  CALL E013Sum(QuadSc%defW)

  DO i=1,QuadSc%ndof
!    Shearrate(i) = QuadSc%defV(i)/QuadSc%defW(i)
   Viscosity(i) = QuadSc%defU(i)/QuadSc%defW(i)
  END DO

 END IF

END SUBROUTINE  GetNonNewtViscosity
!
! ----------------------------------------------
!
SUBROUTINE  GetNonNewtViscosityOld()
  INTEGER i
  REAL*8 daux
  REAL*8 HogenPowerlaw
  LOGICAL bCondition
  REAL*8 ViscosityModel

  bCondition = .FALSE.

  IF (bNonNewtonian) THEN
    DO i=1,SIZE(myExport%Fields)
    IF (ADJUSTL(TRIM(myExport%Fields(i))).EQ.'Viscosity') bCondition=.TRUE.
    END DO
  END IF

  IF (bCondition) THEN
    ILEV = NLMAX
    CALL SETLEV(2)

    CALL GetGradVelo_rhs(QuadSc,QuadSc%ValU)
    CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
    CALL GetGradVelo_val(QuadSc,1,Properties%Density(1))

    CALL GetGradVelo_rhs(QuadSc,QuadSc%ValV)
    CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
    CALL GetGradVelo_val(QuadSc,2,Properties%Density(1))

    CALL GetGradVelo_rhs(QuadSc,QuadSc%ValW)
    CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
    CALL GetGradVelo_val(QuadSc,3,Properties%Density(1))

    DO i=1,SIZE(QuadSc%ValU)
    daux = QuadSc%ValUx(i)**2d0 + QuadSc%ValVy(i)**2d0 + QuadSc%ValWz(i)**2d0 + &
      0.5d0*(QuadSc%ValUy(i)+QuadSc%ValVx(i))**2d0 + &
      0.5d0*(QuadSc%ValUz(i)+QuadSc%ValWx(i))**2d0 + &
      0.5d0*(QuadSc%ValVz(i)+QuadSc%ValWy(i))**2d0

    if(allocated(Shearrate))then
      Shearrate(i) = sqrt(2d0 * daux)
    end if
    Viscosity(i) = ViscosityModel(daux)

    END DO

  END IF

END SUBROUTINE  GetNonNewtViscosityOld
!
! ----------------------------------------------
!
SUBROUTINE Calculate_Torque(mfile)
implicit none
INTEGER mfile,i
REAL*8 Torque1(3), Torque2(3),dVolFlow1,dVolFlow2,myPI,daux
REAL*8 dHeat,Ml_i,Shear,Visco

integer :: ilevel

EXTERNAL E013

ilevel = mg_mesh%nlmax

 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque1, E013,103)

IF (myid.ne.0) then
 call IntegrateFlowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow1,0.0d0)

 call IntegrateFlowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow2,mySigma%L)

END IF

dHeat = 0d0
DO i=1,QuadSc%ndof
 IF (MixerKNPR(i).eq.0) THEN
  Shear = Shearrate(i)
  Visco = 0.1d0*Viscosity(i)
  Ml_i = mg_MlRhoMat(NLMAX)%a(i)*1e-6
  dHeat = dHeat + Ml_i*Shear*Shear*Visco
 END IF
END DO

CALL COMM_SUMM(dVolFlow1)
CALL COMM_SUMM(dVolFlow2)
CALL COMM_SUMM(dHeat)

myPI = dATAN(1d0)*4d0
daux = 1D0*1e-7*myPI*(myProcess%umdr/3d1)

IF (myid.eq.showID) THEN
  WRITE(MTERM,5)
  WRITE(MFILE,5)
  write(mfile,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]:",timens,dVolFlow1*3.6d0,dVolFlow1*3.6d0*myThermodyn%density
  write(mterm,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]:",timens,dVolFlow1*3.6d0,dVolFlow1*3.6d0*myThermodyn%density
  write(mfile,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]:",timens,dVolFlow2*3.6d0,dVolFlow2*3.6d0*myThermodyn%density
  write(mterm,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]:",timens,dVolFlow2*3.6d0,dVolFlow2*3.6d0*myThermodyn%density
  write(mfile,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
  write(mterm,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
!  WRITE(666,'(7G16.8)') Timens,Torque1,Torque2 
END IF

5  FORMAT(100('-'))

END SUBROUTINE Calculate_Torque
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
SUBROUTINE GetMeshVelocity()
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
    WRITE(*,*)  "Maximum Mesh Velocity: ", SQRT(dMaxVelo)
  END IF

END SUBROUTINE GetMeshVelocity
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
subroutine InitOperators(mfile, mgMesh)
use PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
use var_QuadScalar, only : tMultiMesh

implicit none

integer, intent(in) :: mfile

type(tMultiMesh), intent(inout) :: mgMesh

! local variables
integer :: i

call StoreOrigCoor(mgMesh%level(mgMesh%nlmax)%dcorvg)
call store_old_mesh(mgMesh%level(mgMesh%nlmax)%dcorvg)

ilev = mgMesh%nlmax
call setlev(2)

if (myid.ne.0) call updateMixerGeometry(mfile)

call OperatorRegenaration(1)
call OperatorRegenaration(2)
call OperatorRegenaration(3)

end subroutine InitOperators
!
! ----------------------------------------------
!
END MODULE Transport_Q2P1
