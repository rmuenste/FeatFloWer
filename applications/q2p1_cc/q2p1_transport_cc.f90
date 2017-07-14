MODULE Transport_CC

USE Transport_Q2P1
! USE PP3D_MPI
USE PP3D_MPI, ONLY:myid,master,E011Sum,COMM_Maximum,&
                   COMM_NLComplete,Comm_Summ,myMPI_Barrier,Barrier_myMPI
USE Parametrization,ONLY : InitBoundaryStructure,myParBndr
! USE PP3D_MPI, ONLY:E011Sum,E011True_False,Comm_NLComplete,&
!               Comm_Maximum,Comm_Summ,knprmpi,myid,master
! USE LinScalar, ONLY: AddSurfaceTension
use def_cc
IMPLICIT NONE


CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE Init_Q2P1_Structures_cc(mfile)
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

 CALL InitBoundaryList(KWORK(L(LNPR)),&
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
 CALL Create_BMat_iso() !(B,BT)

 IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "

 IF (myid.ne.master) THEN
  ! Parallel E012/E013 matrix structure
  CALL Create_QuadLinParMatStruct(PLinSc) !(pB)

  ! Building up the Parallel E012/E012 matrix strucrures
  CALL Create_ParLinMatStruct ()
 END IF

 CALL Create_P1MMat() !(MP1,iMP1)

 IF (QuadSc%prm%MGprmIn%VANKA.eq.1) THEN
   CALL Create_Special_CCStructures()
 END IF 
 
 CALL Create_GlobalNumbering_CC()

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
  ALLOCATE (ParKNPR(NVT+NET+NAT+NEL))
  QuadSc%auxU = 1d0
  CALL E013Sum(QuadSc%auxU)
  DO I=1,NVT+NET+NAT+NEL
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

 CALL OperatorRegenaration_iso(1)


END SUBROUTINE Init_Q2P1_Structures_cc
!
! ----------------------------------------------
!
SUBROUTINE Transport_q2p1_UxyzP_cc(mfile,inl_u)
implicit none

INTEGER mfile,INL,inl_u
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit,iIter,iIterges
INTEGER INLComplete,I,J,IERR,iOuter
REAL*8  DefNormUVWP0(4),DefNormUVWP(4),DefNorm0,DefNorm,ni
REAL*8 :: FORCES_NEW(3),FORCES_OLD(3),MGnonLin,digitcriterion
REAL*8 stopOne,diffOne,diffTwo,myTolerance,alpha,DefNormOld

thstep = 0d0
FORCES_OLD = 0d0
iIterges = 0d0
digitcriterion = 1d-1
myTolerance = QuadSc%prm%MGprmIn%Criterion1
ni = 0d0
alpha = QuadSc%prm%Alpha


 CALL ExchangeVelocitySolutionCoarse()

 CALL OperatorRegenaration_iso(2)

 CALL OperatorRegenaration_iso(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------

 CALL ZTIME(tttt0)

! Assemble the right hand side
 CALL Matdef_general_QuadScalar_cc(QuadSc,1)

IF (myid.ne.master) THEN

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

END IF

thstep = tstep

 CALL ExchangeVelocitySolutionCoarse()

 CALL OperatorRegenaration_iso(3)

! Assemble the defect vector and fine level matrix
 CALL Matdef_general_QuadScalar_cc(QuadSc,-1)
 CALL OperatorDeallocation()

DO i=1,QuadSc%prm%NLmax
ni = ni + 1d0

IF (myid.eq.1) THEN
  write(mfile,66) 
  write(mterm,66)
  write(mfile,'(A,I6)') "Iteration number: ",i
  write(mterm,'(A,I6)') "Iteration number: ",i
  write(mfile,'(A,F6.3,A,ES12.4)')" Parameter for Fixpoint-Newton-Adaptivity:", alpha,", Criterion for MG:", digitcriterion 
  write(mterm,'(A,F6.3,A,ES12.4)')" Parameter for Fixpoint-Newton-Adaptivity:", alpha,", Criterion for MG:", digitcriterion
  write(mfile,5)
  write(mterm,5)
END IF

CALL myMPI_Barrier()
CALL ZTIME(tttt0)

IF (myid.eq.1) WRITE(*,'(A)') " Factorization of element matrices starts:"

 CALL Barrier_myMPI()

IF (QuadSc%prm%MGprmIn%VANKA.eq.0) THEN
  CALL CC_Extraction(QuadSc)
ELSE
  CALL Special_CCExtraction(QuadSc)
END IF

 CALL Barrier_myMPI()

 CALL ZTIME(tttt1)
IF (myid.eq.1) WRITE(*,'(A,ES12.4,A)') " was done in: ",tttt1-tttt0 ," s"

 CALL Special_CC_Coarse(QuadSc)

IF (myid.ne.master) THEN
! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()
 QuadSc%valU_old = QuadSc%valU
 QuadSc%valV_old = QuadSc%valV
 QuadSc%valW_old = QuadSc%valW

 LinSc%valP_old  = LinSc%valP(NLMAX)%x
 CALL CC_GetDefect(QuadSc,LinSc)
 CALL Boundary_QuadScalar_Def()
END IF

 CALL GetDefNorms(QuadSc,LinSc,DefNormUVWP0)
IF (i.eq.1) DefNorm0 = MAX(DefNormUVWP0(1),DefNormUVWP0(2),DefNormUVWP0(3),DefNormUVWP0(4))

IF (myid.ne.master) THEN
 QuadSc%valU          = 0d0 
 QuadSc%valV          = 0d0 
 QuadSc%valW          = 0d0 
 LinSc%valP(NLMAX)%x  = 0d0 
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Solver !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CALL myMPI_Barrier()
 CALL ZTIME(tttt0)

 CALL CC_mgSolve(QuadSc,LinSc,mfile,iIter,digitcriterion)
iIterges = iIterges + iIter

 CALL myMPI_Barrier()
 CALL ZTIME(tttt1)
IF (myid.eq.1) WRITE(*,'(A,ES12.4,A)') " Solution of linear system ... done! --> Time: ",tttt1-tttt0 ," s"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Solver !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

IF (QuadSc%prm%MGprmIn%VANKA.eq.0) THEN
  CALL CC_MemFree()
ELSE
  CALL Special_CCMemFree()
END IF

IF (myid.ne.master) THEN
 QuadSc%valU         = QuadSc%valU         + QuadSc%valU_old
 QuadSc%valV         = QuadSc%valV         + QuadSc%valV_old
 QuadSc%valW         = QuadSc%valW         + QuadSc%valW_old
 LinSc%valP(NLMAX)%x = LinSc%valP(NLMAX)%x + LinSc%valP_old 
! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()
END IF

 CALL ExchangeVelocitySolutionCoarse()

 CALL OperatorRegenaration_iso(3)

! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)
 CALL OperatorDeallocation()

! Set dirichlet boundary conditions on the solution
IF (myid.ne.master) THEN
 CALL Boundary_QuadScalar_Val()
 CALL CC_GetDefect(QuadSc,LinSc)
 CALL Boundary_QuadScalar_Def()
END IF



!!!!!!!!!!!!!!!!!!!! "ADAPTIVITY" !!!!!!!!!!!!!!!!!!!!
 if(i.eq.1) then
	DefNormOld = MAX(DefNormUVWP0(1),DefNormUVWP0(2),DefNormUVWP0(3),DefNormUVWP0(4))
 else
	DefNormOld = DefNorm
 end if

! PLEASE DO NOT COMMENT THESE TWO LINES
 CALL GetDefNorms(QuadSc,LinSc,DefNormUVWP)
 DefNorm = MAX(DefNormUVWP(1),DefNormUVWP(2),DefNormUVWP(3),DefNormUVWP(4))
! PLEASE DO NOT COMMENT THESE TWO LINES

! Fixpoint-Newton-adaptivity
 alpha = (0.2d0 + 1.46d0/(-0.48d0 + EXP(0.94d0*DefNorm/DefNormOld)))*alpha
 if(alpha.gt.1d0) alpha = 1d0
! digits to gain in MG
 digitcriterion = (DefNorm/DefNormOld)**(2d0**alpha)
 if(digitcriterion.gt.1d-1) digitcriterion = 1d-1
!!!!!!!!!!!!!!!!!!!! "ADAPTIVITY" !!!!!!!!!!!!!!!!!!!!



IF (myid.eq.showid) THEN
  write(mfile,5)
  write(mterm,5)
  write(mfile,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') "INITIAL-DEF:",DefNorm0,",  ACTUAL-DEF:",DefNorm,",  ACTUAL-criterion:",DefNorm/DefNorm0,",  needed:",myTolerance
  write(mterm,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') "INITIAL-DEF:",DefNorm0,",  ACTUAL-DEF:",DefNorm,",  ACTUAL-criterion:",DefNorm/DefNorm0,",  needed:",myTolerance
END IF

 CALL FAC_GetForces_CC(mfile,FORCES_NEW)
diffOne = ABS(FORCES_NEW(1)-FORCES_OLD(1))
diffTwo = ABS(FORCES_NEW(2)-FORCES_OLD(2))


!IF (FORCES_NEW(1).LT.1d-7) THEN
!  diffOne = diffOne/1d-7
!ELSE
!  diffOne = diffOne/ABS(FORCES_New(1))
!END IF
!IF (FORCES_NEW(2).LT.1d-7) THEN
!  diffTwo = diffTwo/1d-7
!ELSE
!  diffTwo= diffTwo/ABS(FORCES_New(2))
!END IF


stopOne = max(diffOne,diffTwo)

FORCES_OLD = FORCES_NEW

!!!!!!!!!!!!!!! STOPPING CRITERION !!!!!!!!!!!!!!!!!!!!
IF (DefNorm/DefNorm0.LT.myTolerance) exit
!IF (stopOne.LT.1d-5) THEN
!	IF (myid.eq.1) THEN
!	write(mfile,55) 
!  	write(mterm,55)
!	write(mfile,*) " !!!! FORCES REACHED CONVERGENCE CRITERION !!!!"
!  	write(mterm,*) " !!!! FORCES REACHED CONVERGENCE CRITERION !!!!"
!        END IF
!        exit
!END IF
END DO

IF (myid.eq.showid) THEN
MGnonLin = iIterges/ni
  write(mfile,55) 
  write(mterm,55)
  write(mfile,'(A,F6.3)')"     AVERAGE MG ITERATIONS:", MGnonLin
  write(mterm,'(A,F6.3)')"     AVERAGE MG ITERATIONS:", MGnonLin
END IF


 CALL GetNonNewtViscosity()

CALL QuadScP1toQ2_cc(LinSc,QuadSc)

IF (myid.ne.0) THEN
 CALL STORE_OLD_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
END IF

66  FORMAT(104('_'))
5  FORMAT(104('-'))
55  FORMAT(104('='))
END SUBROUTINE Transport_q2p1_UxyzP_cc
!
! ----------------------------------------------
!
SUBROUTINE ExchangeVelocitySolutionCoarse()
implicit none
integer :: ndof

ILEV=NLMIN
ndof = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net +&
       mg_mesh%level(ilev)%nat + mg_mesh%level(ilev)%nel

CALL ExchangeVelocitySolutionCoarseSub(QuadSc%ValU,&
                                       QuadSc%ValV,&
                                       QuadSc%ValW,&
                                       QuadSc%AuxU,&
                                       QuadSc%AuxV,&
                                       QuadSc%AuxW,&
                                       ndof)


END SUBROUTINE ExchangeVelocitySolutionCoarse
!
! ----------------------------------------------
!
SUBROUTINE InitCond_QuadScalar_cc()
INTEGER i

ILEV=NLMAX
CALL SETLEV(2)

 ! Set initial conditions
 IF (myid.ne.0)then
  IF (myFBM%nParticles.GT.0) THEN
    CALL updateFBMGeometry()
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

 CALL QuadScalar_InitCond_cc()

 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

 ! Set initial conditions
 CALL LinScalar_InitCond_cc(mg_mesh%level(ilev)%dcorvg,&
                            mg_mesh%level(ilev)%kvert)

END IF

END SUBROUTINE InitCond_QuadScalar_cc
!
! ----------------------------------------------
!
SUBROUTINE QuadScalar_InitCond_cc()
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

END SUBROUTINE QuadScalar_InitCond_cc
!
!----------------------------------------------
!
SUBROUTINE LinScalar_InitCond_cc(dcorvg,kvert)
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

END SUBROUTINE LinScalar_InitCond_cc
!
! ----------------------------------------------
!
SUBROUTINE OperatorDeallocation()

 !DO ILEV=NLMIN,NLMAX

 DO ILEV=1,1

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (ALLOCATED(mg_S11mat(ILEV)%a)) DEALLOCATE(mg_S11mat(ILEV)%a)
  IF (ALLOCATED(mg_S22mat(ILEV)%a)) DEALLOCATE(mg_S22mat(ILEV)%a)
  IF (ALLOCATED(mg_S33mat(ILEV)%a)) DEALLOCATE(mg_S33mat(ILEV)%a)
  IF (ALLOCATED(mg_S12mat(ILEV)%a)) DEALLOCATE(mg_S12mat(ILEV)%a)
  IF (ALLOCATED(mg_S13mat(ILEV)%a)) DEALLOCATE(mg_S13mat(ILEV)%a)
  IF (ALLOCATED(mg_S23mat(ILEV)%a)) DEALLOCATE(mg_S23mat(ILEV)%a)
  IF (ALLOCATED(mg_S21mat(ILEV)%a)) DEALLOCATE(mg_S21mat(ILEV)%a)
  IF (ALLOCATED(mg_S31mat(ILEV)%a)) DEALLOCATE(mg_S31mat(ILEV)%a)
  IF (ALLOCATED(mg_S32mat(ILEV)%a)) DEALLOCATE(mg_S32mat(ILEV)%a)

  IF (ALLOCATED(mg_Kmat(ILEV)%a)) DEALLOCATE(mg_KMat(ILEV)%a)

END DO

END SUBROUTINE OperatorDeallocation
!
! ----------------------------------------------
!
SUBROUTINE FAC_GetForces_CC(mfile,Force)
INTEGER mfile
!REAL*8 :: Force(3),U_mean=1.0d0,R=0.5d0,dens_const=1.0d0,Factor
REAL*8 :: Force(3),U_mean=0.2d0,H=0.05d0,D=0.1d0,dens_const=1.0d0,Factor
REAL*8 :: PI=dATAN(1d0)*4d0 
REAL*8 :: Force2(3)
INTEGER i,nn
EXTERNAL E013

 ILEV=NLMAX
 CALL SETLEV(2)
 IF (bNonNewtonian) THEN

 CALL GetForceCyl_cc(QuadSc%valU,QuadSc%valV,&
                     QuadSc%valW,LinSc%valP(NLMAX)%x,&
                     BndrForce,&
                     mg_mesh%level(ILEV)%kvert,&
                     mg_mesh%level(ILEV)%karea,&
                     mg_mesh%level(ILEV)%kedge,&
                     mg_mesh%level(ILEV)%dcorvg,&
                     Force, E013)

 ELSE
  CALL EvaluateDragLift_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
  LinSc%valP(NLMAX)%x,BndrForce,Force)
 END IF

!Pipe
! Factor = 2d0/(dens_const*U_mean*U_mean*PI*R*R)
!Halfpipe/Quarterpipe
! Factor = 2d0/(dens_const*U_mean*U_mean*PI*R*R/2d0)
!FAC
 Factor = 2d0/(dens_const*U_mean*U_mean*H*D)
 Force = Factor*Force

 IF (myid.eq.showID) THEN
  WRITE(MTERM,5)
  WRITE(MFILE,5)
  write(mfile,'(A30,4E16.8)') "Force acting on the cylinder:",timens,Force
  write(mterm,'(A30,4E16.8)') "Force acting on the cylinder:",timens,Force
 END IF

5  FORMAT(104('-'))

END SUBROUTINE FAC_GetForces_CC
!
!----------------------------------------------
!
SUBROUTINE myFAC_GetForces(mfile,Force)
INTEGER mfile
!REAL*8 :: Force(3),U_mean=1.0d0,R=0.5d0,dens_const=1.0d0,Factor
REAL*8 :: Force(3),U_mean=0.2d0,H=0.05d0,D=0.1d0,dens_const=1.0d0,Factor
REAL*8 :: PI=dATAN(1d0)*4d0 
REAL*8 :: Force2(3)
INTEGER i,nn
EXTERNAL E013

 ILEV=NLMAX
 CALL SETLEV(2)

 
 IF (bNonNewtonian) THEN
   ! At the moment we have this case 
   CALL EvaluateDragLift9_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
   LinSc%valP(NLMAX)%x,BndrForce,Force)
 ELSE
  CALL EvaluateDragLift_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
  LinSc%valP(NLMAX)%x,BndrForce,Force)
 END IF

!Pipe
! Factor = 2d0/(dens_const*U_mean*U_mean*PI*R*R)
!Halfpipe/Quarterpipe
! Factor = 2d0/(dens_const*U_mean*U_mean*PI*R*R/2d0)
!FAC
 Factor = 2d0/(dens_const*U_mean*U_mean*H*D)
 Force = Factor*Force


 IF (myid.eq.showID) THEN
  WRITE(MTERM,5)
  WRITE(MFILE,5)
  write(mfile,'(A30,4E16.8)') "Force acting on the sphere:",timens,Force
  write(mterm,'(A30,4E16.8)') "Force acting on the sphere:",timens,Force
  WRITE(666,'(7G16.8)') Timens,Force
 END IF

5  FORMAT(104('-'))

END SUBROUTINE myFAC_GetForces
!
!----------------------------------------------
!
SUBROUTINE OperatorRegenaration_iso(iType)
INTEGER iType
LOGICAL bHit

bHit = .FALSE.

IF (iType.EQ.myMatrixRenewal%D) THEN
 CALL Create_DiffMat_iso(QuadSc)
 bHit = .TRUE.
END IF

IF (iType.EQ.myMatrixRenewal%K) THEN
 CALL Create_KMat_iso(QuadSc)
 CALL Create_barMMat_iso(QuadSc)
 bHit = .TRUE.
END IF

IF (iType.EQ.myMatrixRenewal%M) THEN
 CALL Create_MRhoMat_iso()
 bHit = .TRUE.
END IF

IF (iType.EQ.myMatrixRenewal%S) THEN
 CALL Create_SMat_iso(QuadSc)
 bHit = .TRUE.
END IF

IF (iType.EQ.myMatrixRenewal%C) THEN

 CALL Create_BMat_iso()
 
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

END SUBROUTINE OperatorRegenaration_iso
END MODULE Transport_CC


