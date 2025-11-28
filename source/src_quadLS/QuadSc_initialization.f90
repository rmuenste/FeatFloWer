!=========================================================================
! QuadSc_initialization.f90
!
! Initialization routines for Q2/P1 transport solver structures
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE Init_QuadScalar_Structures_sse(mfile)
use, intrinsic :: ieee_arithmetic
use param_parser, only: GetPhysiclaParameters
implicit none
LOGICAL bExist
INTEGER I,J,ndof,mfile,LevDif
integer :: mydof
integer :: maxlevel
Real*8 :: dabl
real*8 :: myInf
integer :: idx


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

 ALLOCATE (FictKNPR_uint64(mydof))
 ! loop and initialize
 do idx = 1,mydof
   FictKNPR_uint64(idx)%bytes(:) = -1
 end do

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

!  CALL InitializeProlRest(QuadSc,LinSc)
!
!  CALL OperatorRegenaration(1)
!
!  CALL SetUp_HYPRE_Solver(LinSc,PLinSc,mfile)

END SUBROUTINE Init_QuadScalar_Structures_sse
!=========================================================================
!
!=========================================================================
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
!=========================================================================
!
!=========================================================================
SUBROUTINE Init_QuadScalar_Structures_sse_PF(mfile)
use param_parser, only: GetPhysiclaParameters
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
!=========================================================================
!
!=========================================================================
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
!=========================================================================
!
!=========================================================================
SUBROUTINE InitCond_Velocity_PF()

QuadSc%ValU = QuadSc%ValU
QuadSc%ValV = QuadSc%ValV
QuadSc%ValW = QuadSc%ValW + myProcess%umdr/6d1*mySigma%mySegment(1)%t

END SUBROUTINE InitCond_Velocity_PF
!=========================================================================
!
!=========================================================================
SUBROUTINE Init_QuadScalar_ReducedStuctures(mfile)
use param_parser, only: GetPhysiclaParameters
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
!
!  ! Iteration matrix (only allocation)
!  CALL Create_AMat() !(A)
!
!  ! Building up the E012/E013 E013/E012 and matrix structures
!  CALL Create_QuadLinMatStruct()
!
!  ! Building up the E012/E012 matrix strucrures
!  CALL Create_LinMatStruct ()
!
!  ! Pressure gradient matrix
!  CALL Create_BMat() !(B,BT)
!
!  IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "
!
!  IF (myid.ne.master) THEN
!   ! Parallel E012/E013 matrix structure
!   CALL Create_QuadLinParMatStruct(PLinSc) !(pB)
!
!   ! Building up the Parallel E012/E012 matrix strucrures
!   CALL Create_ParLinMatStruct ()
!
!   CALL BuildUpPressureCommunicator(LinSc,PLinSc)
! END IF

!  ! Correct the wall BCs
!  IF (allocated(mg_mesh%BndryNodes))  then
!   ilev = nlmax
!   CALL RestrictWallBC()
!  END IF
!
! ! Set up the boundary condition types (knpr)
!  DO ILEV=NLMIN,NLMAX
!   CALL SETLEV(2)
!   CALL QuadScalar_Knpr()
!  END DO
!
!  ILEV=NLMAX
!  mydof = mg_mesh%level(ilev)%nvt+&
!          mg_mesh%level(ilev)%net+&
!          mg_mesh%level(ilev)%nat+&
!          mg_mesh%level(ilev)%nel
!
!  ALLOCATE (FictKNPR(mydof))
!  FictKNPR=0
!  ALLOCATE (Distance(mydof))
!  Distance = 0d0
!
!  ALLOCATE (MixerKNPR(mydof))
!  MixerKNPR=0
!  ALLOCATE (Distamce(mydof))
!  Distamce = 0d0
!
!  ! SEt up the knpr vector showing dofs with parallel property ...
!  IF (myid.ne.0) THEN
!   ALLOCATE (ParKNPR(mydof))
!   QuadSc%auxU = 1d0
!   CALL E013Sum(QuadSc%auxU)
!   DO I=1,mydof
!    IF (QuadSc%auxU(I).EQ.1d0) THEN
!     ParKNPR(I) = 0
!    ELSE
!     ParKNPR(I) = 1
!    END IF
!   END DO
!  END IF
!
!  IF (myid.eq.showID) THEN
!   INQUIRE (FILE="_data/BenchValues.txt", EXIST=bExist)
!   IF (ISTART.EQ.0.OR.(.NOT.bExist)) THEN
!    OPEN(666,FILE="_data/BenchValues.txt")
!    WRITE(666,'(10A16)') "Time","Drag","Lift","ZForce","ForceVx","ForceVy","ForceVz","ForcePx","ForcePy","ForcePz"
!   ELSE
!    OPEN(666,FILE="_data/BenchValues.txt",ACCESS='APPEND')
!   END IF
!  END IF
!
!  CALL InitializeProlRest(QuadSc,LinSc)
!
!  CALL OperatorRegenaration(1)
!
 CALL Create_MMat()
!
!  CALL SetUp_HYPRE_Solver(LinSc,PLinSc,mfile)

END SUBROUTINE Init_QuadScalar_ReducedStuctures
!=========================================================================
!
!=========================================================================
SUBROUTINE Init_QuadScalar_Stuctures(mfile)
use param_parser, only: GetPhysiclaParameters
implicit none
LOGICAL bExist
INTEGER I,J,ndof,mfile,LevDif
integer :: mydof
integer :: maxlevel, idx
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


 ALLOCATE (FictKNPR_uint64(mydof))
 ! loop and initialize
 do idx = 1,mydof
   FictKNPR_uint64(idx)%bytes(:) = -1
 end do

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

 CALL Create_GradDivMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW,LinSc%knprP)

 CALL OperatorRegenaration(1)

 CALL Create_MMat()

 CALL SetUp_HYPRE_Solver(LinSc,PLinSc,mfile)

END SUBROUTINE Init_QuadScalar_Stuctures
!=========================================================================
!
!=========================================================================
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
!=========================================================================
!
!=========================================================================
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
!=========================================================================
!
!=========================================================================
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
