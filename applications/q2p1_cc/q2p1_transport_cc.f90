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
use def_QuadScalar
use var_QuadScalar_newton, ONLY:zeitstep
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
 IF (ccParams%BDF.ne.0) THEN
 	ALLOCATE(QuadSc%valU_old1(QuadSc%ndof))
 	ALLOCATE(QuadSc%valU_help(QuadSc%ndof))
 	ALLOCATE(QuadSc%valV_old1(QuadSc%ndof))
 	ALLOCATE(QuadSc%valV_help(QuadSc%ndof))
 	ALLOCATE(QuadSc%valW_old1(QuadSc%ndof))
 	ALLOCATE(QuadSc%valW_help(QuadSc%ndof))
 	ALLOCATE(QuadSc%valU_old2(QuadSc%ndof))
 	ALLOCATE(QuadSc%valV_old2(QuadSc%ndof))
 	ALLOCATE(QuadSc%valW_old2(QuadSc%ndof))
 END IF
 ! Initialize the scalar quantity
 CALL InitializeLinScalar(LinSc)

 if (.not.ALLOCATED(myMultiMat%Mat)) then
  myMultiMat%nOfMaterials = 1
  ALLOCATE(myMultiMat%Mat(myMultiMat%nOfMaterials))
  myMultiMat%Mat(1)%Rheology%Equation = 5
  myMultiMat%Mat(1)%Rheology%AtFunc = 1
 end if

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

                     
 ALLOCATE (Temperature(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))
 Temperature = 1d0                    

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
 CALL Create_AMat_new() !(A)

 ! Building up the E012/E013 E013/E012 and matrix structures
 CALL Create_QuadLinMatStruct() 

 ! Building up the E012/E012 matrix strucrures
 CALL Create_LinMatStruct ()

 ! Pressure gradient matrix
 CALL Create_BMat_mod() !(B,BT)

 IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "

 IF (myid.ne.master) THEN
  ! Parallel E012/E013 matrix structure
  CALL Create_QuadLinParMatStruct(PLinSc) !(pB)

  ! Building up the Parallel E012/E012 matrix strucrures
  CALL Create_ParLinMatStruct ()
 END IF

 CALL Create_P1MMat() !(MP1,iMP1)

 IF (ccParams%VANKA.eq.1) THEN
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

  INQUIRE (FILE="_data/BenchValues_time.txt", EXIST=bExist)
  IF (ISTART.EQ.0.OR.(.NOT.bExist)) THEN
   OPEN(777,FILE="_data/BenchValues_time.txt")
   WRITE(777,'(6A16)') "Time","Drag","Lift","ZForce","MGiter","nonLin"
  ELSE
   OPEN(777,FILE="_data/BenchValues_time.txt",ACCESS='APPEND')
  END IF
 END IF

 CALL InitializeProlRest_cc(ccParams)

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
REAL*8 :: FORCES_NEW(7),FORCES_OLD(7),MGnonLin,digitcriterion
REAL*8 stopOne,diffOne,diffTwo,myTolerance,alpha,DefNormOld
REAL*8 B,a,h1,h2,h3

thstep = 0d0
FORCES_OLD = 0d0
iIterges = 0d0
myTolerance = ccParams%StoppingCriterion
ni = 0d0
tsm = ccParams%BDF

!Adaptivity
! F(x) = a + h1/(h2+exp(h3 x)); F(0) = B, F(1) = 0.9, F(0.5) = B/4 + 3/4
alpha = ccParams%Alpha
B = ccParams%ValAdap(1)
a = ccParams%ValAdap(2)
h1 = 1d0/(9d0*(10d0*a-9d0)*(B-1d0)**2d0) * (a**3d0*(96d0-80d0*B)+12d0*a**2d0*(5d0*B**2d0+6d0*B-15d0)+3d0*a*(10d0*B**3d0-57d0*B**2d0+36d0*B+27d0)+B*(-10d0*B**3d0+3d0*B**2d0+72d0*B-81d0))
h2 = 2d0*(a-B)*(8d0*a*(5d0*B-6d0)+5d0*B**2d0-42d0*B+45d0)/(9d0*(10d0*a-9d0)*(B-1d0)**2d0)
h3 = 2d0*log(2d0*(5d0*B-3d0)*(a-B)/(3d0*(10d0*a-9d0)*(B-1d0)))

IF (myid.eq.1) THEN
  write(mfile,77) 
  write(mterm,77)
  write(mfile,'(A,F6.3,A,F6.3,A,F6.3,A,F6.3,A)')"	ADAPTIVE-FUNCTION: ",a," + ",h1," / ( ",h2," + exp( ",h3," x ))" 
  write(mterm,'(A,F6.3,A,F6.3,A,F6.3,A,F6.3,A)')"	ADAPTIVE-FUNCTION: ",a," + ",h1," / ( ",h2," + exp( ",h3," x ))"
  write(mfile,77)
  write(mterm,77)
END IF

digitcriterion = 10d0**(-1d0-alpha)



 CALL ExchangeVelocitySolutionCoarse()

 CALL OperatorRegenaration_iso(2)

 CALL OperatorRegenaration_iso(3)

IF (ccParams%BDF.ne.0) THEN
!#######################################
! Store the old solution for BDF to help
!#######################################
IF (myid.ne.master) THEN
 QuadSc%valU_help = QuadSc%valU
 QuadSc%valV_help = QuadSc%valV
 QuadSc%valW_help = QuadSc%valW
END IF
END IF

IF (myid.ne.master) THEN
 LinSc%P_old  = LinSc%valP(NLMAX)%x
END IF


 CALL ZTIME(tttt0)

! Assemble the right hand side
 zeitstep = tstep*(1d0-theta)
 CALL Matdef_general_QuadScalar_cc(QuadSc,1,alpha)

IF (myid.ne.master) THEN

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW
END IF

!##########################################
! Set zeitstep as the right scaling factor
!##########################################
IF (tsm.EQ.0 .OR. itns.EQ.1) THEN
	zeitstep = tstep*theta
	if (myid.eq.showid) then
		write(mfile,'(A)')'BE/CN'
		write(mterm,'(A)')'BE/CN'
	end if
ELSE IF (tsm.EQ.2 .OR. itns.EQ.2) THEN
	zeitstep = 2d0/3d0*tstep
	if (myid.eq.showid) then
		write(mfile,'(A)')'BDF(2)'
		write(mterm,'(A)')'BDF(2)'
	end if
ELSE IF (tsm.EQ.3) THEN
	zeitstep = 6d0/11d0*tstep
	if (myid.eq.showid) then
		write(mfile,'(A)')'BDF(3)'
		write(mterm,'(A)')'BDF(3)'
	end if
END IF


 CALL ExchangeVelocitySolutionCoarse()

! Assemble the defect vector and fine level matrix
 CALL Matdef_general_QuadScalar_cc(QuadSc,-1,alpha)


! OUTPUT at the beginning
! CALL myOutputMatrix("AA11",qMat,AA11mat)
! CALL myOutputMatrix("AA12",qMat,AA12mat)
! CALL myOutputMatrix("AA13",qMat,AA13mat)
! CALL myOutputMatrix("AA21",qMat,AA21mat)
! CALL myOutputMatrix("AA22",qMat,AA22mat)
! CALL myOutputMatrix("AA23",qMat,AA23mat)
! CALL myOutputMatrix("AA31",qMat,AA31mat)
! CALL myOutputMatrix("AA32",qMat,AA32mat)
! CALL myOutputMatrix("AA33",qMat,AA33mat)
!
! CALL myOutputMatrix("Bx",qlMat,BXMat)
! CALL myOutputMatrix("By",qlMat,BYMat)
! CALL myOutputMatrix("Bz",qlMat,BZMat)
! CALL myOutputMatrix("BTx",lqMat,BTXMat)
! CALL myOutputMatrix("BTy",lqMat,BTYMat)
! CALL myOutputMatrix("BTz",lqMat,BTZMat)


 CALL OperatorDeallocation()

DO i=1,ccParams%NLmax
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

IF (myid.eq.1) WRITE(*,'(A)') " Factorization of element matrices"

 CALL Barrier_myMPI()

IF (ccParams%VANKA.eq.0) THEN
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
 CALL Matdef_general_QuadScalar_cc(QuadSc,-1,alpha)


! CALL OperatorDeallocation() ! after force calculation

! Set dirichlet boundary conditions on the solution
IF (myid.ne.master) THEN
 CALL Boundary_QuadScalar_Val()
 CALL CC_GetDefect(QuadSc,LinSc)
 CALL Boundary_QuadScalar_Def()
END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!---------------------------
! Fixpoint-Newton-adaptivity
!---------------------------
!	continous
 alpha = (a+h1/(h2+EXP(h3 * DefNorm/DefNormOld)))*alpha
 if(alpha.gt.1d0) alpha = 1d0
!	discontinous
! if(DefNorm/DefNormOld.lt.1d0) then
!	alpha = 1.25d0*alpha
! else
!	alpha = 0.75d0*alpha
! end if
! if(alpha.gt.1d0) alpha = 1d0

!---------------------
! digits to gain in MG
!---------------------
 digitcriterion = (DefNorm/DefNormOld)**(2d0**alpha)
 if(digitcriterion.gt.10d0**(-1d0-alpha)) digitcriterion = 10d0**(-1d0-alpha)
!!!!!!!!!!!!!!!!!!!! "ADAPTIVITY" !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF (myid.eq.showid) THEN
  write(mfile,5)
  write(mterm,5)
  IF (bSteadyState) THEN
  write(mfile,'(A,ES12.4,A,ES12.4,A,ES12.4)') "INITIAL-DEF:",DefNorm0,",  ACTUAL-DEF:",DefNorm,",  needed:",myTolerance
  write(mterm,'(A,ES12.4,A,ES12.4,A,ES12.4)') "INITIAL-DEF:",DefNorm0,",  ACTUAL-DEF:",DefNorm,",  needed:",myTolerance
  ELSE
  write(mfile,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') "INITIAL-DEF:",DefNorm0,",  ACTUAL-DEF:",DefNorm,",  ACTUAL-criterion:",DefNorm/DefNorm0,",  needed:",myTolerance
  write(mterm,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') "INITIAL-DEF:",DefNorm0,",  ACTUAL-DEF:",DefNorm,",  ACTUAL-criterion:",DefNorm/DefNorm0,",  needed:",myTolerance
  END IF
END IF

 CALL FAC_GetForces_CC(mfile,FORCES_NEW)
 CALL OperatorDeallocation()

diffOne = ABS(FORCES_NEW(1)-FORCES_OLD(1))
diffTwo = ABS(FORCES_NEW(2)-FORCES_OLD(2))


IF (FORCES_NEW(1).LT.1d-7) THEN
  diffOne = diffOne/1d-7
ELSE
  diffOne = diffOne/ABS(FORCES_New(1))
END IF
IF (FORCES_NEW(2).LT.1d-7) THEN
  diffTwo = diffTwo/1d-7
ELSE
  diffTwo= diffTwo/ABS(FORCES_New(2))
END IF


stopOne = max(diffOne,diffTwo)

FORCES_OLD = FORCES_NEW

!!!!!!!!!!!!!!! STOPPING CRITERION !!!!!!!!!!!!!!!!!!!!

IF (bSteadyState) THEN
	IF(DefNorm.LT.myTolerance) exit
ELSE 
	IF (DefNorm/DefNorm0.LT.myTolerance) exit
END IF

IF (DefNorm.EQ.0d0) exit

!IF (stopOne.LT.1d-4) THEN
!	IF (myid.eq.showid) THEN
!	write(mfile,55) 
!  	write(mterm,55)
!	write(mfile,*) " !!!! FORCES REACHED CONVERGENCE CRITERION !!!!"
!  	write(mterm,*) " !!!! FORCES REACHED CONVERGENCE CRITERION !!!!"
!        END IF
!        exit
!END IF
END DO

IF (ccParams%BDF.ne.0) THEN
!#########################################
! Store the old solution for BDF to old1/2
!#########################################
IF (myid.ne.master) THEN 
 QuadSc%valU_old2 = QuadSc%valU_old1
 QuadSc%valV_old2 = QuadSc%valV_old1
 QuadSc%valW_old2 = QuadSc%valW_old1
 QuadSc%valU_old1 = QuadSc%valU_help
 QuadSc%valV_old1 = QuadSc%valV_help
 QuadSc%valW_old1 = QuadSc%valW_help
END IF
END IF

! OUTPUT at the end
! CALL Matdef_general_QuadScalar_cc(QuadSc,-1,alpha)
! CALL myOutputMatrix("EA11",qMat,AA11mat)
! CALL myOutputMatrix("EA12",qMat,AA12mat)
! CALL myOutputMatrix("EA13",qMat,AA13mat)
! CALL myOutputMatrix("EA21",qMat,AA21mat)
! CALL myOutputMatrix("EA22",qMat,AA22mat)
! CALL myOutputMatrix("EA23",qMat,AA23mat)
! CALL myOutputMatrix("EA31",qMat,AA31mat)
! CALL myOutputMatrix("EA32",qMat,AA32mat)
! CALL myOutputMatrix("EA33",qMat,AA33mat)
! CALL OperatorDeallocation()

IF (myid.eq.showid) THEN
  WRITE(777,'(10G16.8)') Timens,FORCES_NEW,iIterges,ni

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
77  FORMAT(104('*'))
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

DO ILEV=NLMIN,NLMAX

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

  IF (myMatrixRenewal%K.ne.0) THEN  
    IF (ALLOCATED(mg_Kmat(ILEV)%a)) DEALLOCATE(mg_KMat(ILEV)%a)

    IF (ALLOCATED(mg_barM11mat(ILEV)%a)) DEALLOCATE(mg_barM11mat(ILEV)%a)
    IF (ALLOCATED(mg_barM12mat(ILEV)%a)) DEALLOCATE(mg_barM12mat(ILEV)%a)
    IF (ALLOCATED(mg_barM13mat(ILEV)%a)) DEALLOCATE(mg_barM13mat(ILEV)%a)
    IF (ALLOCATED(mg_barM21mat(ILEV)%a)) DEALLOCATE(mg_barM21mat(ILEV)%a)
    IF (ALLOCATED(mg_barM22mat(ILEV)%a)) DEALLOCATE(mg_barM22mat(ILEV)%a)
    IF (ALLOCATED(mg_barM23mat(ILEV)%a)) DEALLOCATE(mg_barM23mat(ILEV)%a)
    IF (ALLOCATED(mg_barM31mat(ILEV)%a)) DEALLOCATE(mg_barM31mat(ILEV)%a)
    IF (ALLOCATED(mg_barM32mat(ILEV)%a)) DEALLOCATE(mg_barM32mat(ILEV)%a)
    IF (ALLOCATED(mg_barM33mat(ILEV)%a)) DEALLOCATE(mg_barM33mat(ILEV)%a)
  END IF


END DO

END SUBROUTINE OperatorDeallocation
!
! ----------------------------------------------
!
SUBROUTINE FAC_GetForces_CC(mfile,Force)
INTEGER mfile
REAL*8 :: Force(7),dens_const=1.0d0,Factor
REAL*8 :: PI=dATAN(1d0)*4d0 
REAL*8 :: Force2(3)
INTEGER i,nn
EXTERNAL E013


 ILEV=NLMAX
 CALL SETLEV(2)
 IF (bNonNewtonian) THEN

! CALL GetForceCyl_cc(QuadSc%valU,QuadSc%valV,&
!                     QuadSc%valW,LinSc%valP(NLMAX)%x,&
!                     BndrForce,&
!                     mg_mesh%level(ILEV)%kvert,&
!                     mg_mesh%level(ILEV)%karea,&
!                     mg_mesh%level(ILEV)%kedge,&
!                     mg_mesh%level(ILEV)%dcorvg,&
!                     Force, E013)

 CALL GetForceCyl_cc_iso(QuadSc%valU,QuadSc%valV,&
                     QuadSc%valW,LinSc%ValP(NLMAX)%x,&
                     LinSc%P_old,BndrForce,&
                     mg_mesh%level(ILEV)%kvert,&
                     mg_mesh%level(ILEV)%karea,&
                     mg_mesh%level(ILEV)%kedge,&
                     mg_mesh%level(ILEV)%dcorvg,&
                     Force, E013,bNonNewtonian)

 END IF

 if (postParams%D.eq.0d0) then
   Factor = 2d0/(dens_const*postParams%U_mean*postParams%U_mean*PI*postParams%H*postParams%H)
 else
   Factor = 2d0/(dens_const*postParams%U_mean*postParams%U_mean*postParams%H*postParams%D)
 end if

 Force = Factor*Force

 IF (myid.eq.showID) THEN
  WRITE(MTERM,5)
  WRITE(MFILE,5)
  write(mfile,'(A,8E16.8)') "Force acting:",timens,Force
  write(mterm,'(A,8E16.8)') "Force acting:",timens,Force
  WRITE(mfile,'(8G16.8)') Timens,Force
  write(mfile,'(A,3E16.8)') "BenchForce:",timens,Force(1),Force(4)
  write(mterm,'(A,3E16.8)') "BenchForce:",timens,Force(1),Force(4)
 END IF

5  FORMAT(104('-'))

END SUBROUTINE FAC_GetForces_CC
!
!----------------------------------------------
!
SUBROUTINE myFAC_GetForces(mfile,Force)
INTEGER mfile
!REAL*8 :: Force(3),U_mean=1.0d0,R=0.5d0,dens_const=1.0d0,Factor
REAL*8 :: Force(3),U_mean=1.0d0,H=0.205d0,D=0.1d0,dens_const=1.0d0,Factor
REAL*8 :: PI=dATAN(1d0)*4d0 
REAL*8 :: Force2(3)
INTEGER i,nn
EXTERNAL E013

 ILEV=NLMAX
 CALL SETLEV(2)

 
 IF (bNonNewtonian) THEN
   ! At the moment we have this case 
   CALL EvaluateDragLift9_mod(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
   LinSc%valP(NLMAX)%x,BndrForce,Force)
 END IF


!Pipe
! Factor = 2d0/(dens_const*U_mean*U_mean*PI*R*R)
! Factor = 2d0/(dens_const*U_mean*U_mean*PI*R*R*1d0/4d0)
!FAC
 Factor = 2d0/(dens_const*U_mean*U_mean*H*D)
 Force = Factor*Force

 IF (myid.eq.showID) THEN
  WRITE(MTERM,5)
  WRITE(MFILE,5)
  write(mfile,'(A30,4E16.8)') "Force acting:",timens,Force
  write(mterm,'(A30,4E16.8)') "Force acting:",timens,Force
  WRITE(666,'(4G16.8)') Timens,Force
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
 CALL Create_DiffMat(QuadSc)
 bHit = .TRUE.
END IF

IF (iType.EQ.myMatrixRenewal%K) THEN
 CALL Create_KMat(QuadSc)
 CALL Create_barMMat_iso(QuadSc)
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

 CALL Create_CMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW,LinSc%knprP,LinSc%prm%MGprmIn%MinLev,LinSc%prm%MGprmIn%CrsSolverType)
 IF (myid.ne.master) THEN
  CALL Create_ParCMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW)
 END IF
 bHit = .TRUE.
END IF

IF (myid.EQ.ShowID.AND.bHit) WRITE(MTERM,'(A)', advance='yes') " "

END SUBROUTINE OperatorRegenaration_iso
!
!----------------------------------------------
!
SUBROUTINE EvaluateDragLift9_mod(U,V,W,P,CYL,df)
REAL*8 U(*),V(*),W(*),P(*)
LOGICAL CYL(*)
REAL*8 df(3)
INTEGER I,J,K,JJ

IF (myid.ne.0) THEN

ILEV=NLMAX
CALL SETLEV(2)

S11Mat   => mg_S11Mat(ILEV)%a
S22Mat   => mg_S22Mat(ILEV)%a
S33Mat   => mg_S33Mat(ILEV)%a
S12Mat   => mg_S12Mat(ILEV)%a
S13Mat   => mg_S13Mat(ILEV)%a
S23Mat   => mg_S23Mat(ILEV)%a
S21Mat   => mg_S21Mat(ILEV)%a
S31Mat   => mg_S31Mat(ILEV)%a
S32Mat   => mg_S32Mat(ILEV)%a
BXMat_new    => mg_BXMat_new(ILEV)%a
BYMat_new    => mg_BYMat_new(ILEV)%a
BZMat_new    => mg_BZMat_new(ILEV)%a
qMat     => mg_qMat(ILEV)

df = 0d0
DO I=1,qMat%nu
 IF (CYL(I)) THEN
  DO J=qMat%LdA(I),qMat%LdA(I+1)-1
   K = qMat%ColA(J)
   df(1) = df(1) - (S11Mat(J)*U(K) + S12Mat(J)*V(K) + S13Mat(J)*W(K))
   df(2) = df(2) - (S21Mat(J)*U(K) + S22Mat(J)*V(K) + S23Mat(J)*W(K))
   df(3) = df(3) - (S31Mat(J)*U(K) + S32Mat(J)*V(K) + S33Mat(J)*W(K))
  END DO
  DO J=qlMat%LdA(I),qlMat%LdA(I+1)-1
   K = qlMat%ColA(J)
   df(1) = df(1) + BXMat_new(J)*P(K)
   df(2) = df(2) + BYMat_new(J)*P(K)
   df(3) = df(3) + BZMat_new(J)*P(K)
  END DO
 END IF
END DO

END IF

CALL COMM_SUMM(df(1))
CALL COMM_SUMM(df(2))
CALL COMM_SUMM(df(3))


END SUBROUTINE EvaluateDragLift9_mod
!
! ----------------------------------------------
!
SUBROUTINE Init_CCParam(mfile)
IMPLICIT NONE
CHARACTER*9 cName
INTEGER mfile
INTEGER iEnd,iAt,iEq,istat
CHARACTER string*100,param*50,cVar*7,cPar*50
LOGICAL bOK
INTEGER :: myFile=90909090
CHARACTER(20) :: out_string

cName = "CCuvwp"
  
OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action='read',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open data file: ",myDataFile  
  stop          
end if

! IF (myid.eq.showid) WRITE(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
! IF (myid.eq.showid) WRITE(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

DO
 READ (UNIT=myFile,FMT='(A100)',IOSTAT=iEnd) string
 IF (iEnd.EQ.-1) EXIT
 CALL StrStuct()
  IF (bOK) THEN

  READ(string(1:iAt-1),*) cVar

  IF (TRIM(ADJUSTL(cVar)).EQ.TRIM(ADJUSTL(cName))) THEN

   READ(string(iAt+1:iEq-1),*) cPar

   SELECT CASE (TRIM(ADJUSTL(cPar)))

    ! Problems with GCC: 
    ! The gfortran compiler does not permit a statement like write(*,'(A,I)') where
    ! the length of the integer I is not specified.
    ! A way to solve this problem is to write the integer to a string:
    ! write(string,'(i20)')ccParams%iMass
    ! that is long enough to hold all reasonable integer values
    ! then adjust the length of the string write it out:
    ! write(*,'(A)')adjustl(string)

    CASE ("NLmin")
    READ(string(iEq+1:),*) ccParams%NLmin
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%NLmin)
    CASE ("NLmax")
    READ(string(iEq+1:),*) ccParams%NLmax
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%NLmax)
    CASE ("Alpha")
    READ(string(iEq+1:),*) ccParams%Alpha
    call write_param_real(mfile,cVar,cPar,out_string,ccParams%Alpha)
    CASE ("ValAdap")
    READ(string(iEq+1:),*) ccParams%ValAdap
    IF (myid.eq.showid) write(mterm,'(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",ccParams%ValAdap
    IF (myid.eq.showid) write(mfile,'(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",ccParams%ValAdap
    CASE ("Stopping")
    READ(string(iEq+1:),*) ccParams%StoppingCriterion
    call write_param_real(mfile,cVar,cPar,out_string,ccParams%StoppingCriterion)
    CASE ("MGMinLev")
    READ(string(iEq+1:),*) ccParams%MinLev
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%MinLev)
    CASE ("MGMedLev")
    READ(string(iEq+1:),*) ccParams%MedLev
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%MedLev)
    CASE ("MGMinIterCyc")
    READ(string(iEq+1:),*) ccParams%MinIterCycle
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%MinIterCycle)
    CASE ("MGMaxIterCyc")
    READ(string(iEq+1:),*) ccParams%MaxIterCycle
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%MaxIterCycle)
    CASE ("MGSmoothSteps")
    READ(string(iEq+1:),*) ccParams%nSmootherSteps
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%nSmootherSteps)
    CASE ("MGCriterion")
    READ(string(iEq+1:),*) ccParams%Criterion
    call write_param_real(mfile,cVar,cPar,out_string,ccParams%Criterion)
    CASE ("MGRelaxPrm")
    READ(string(iEq+1:),*) ccParams%RLX
    call write_param_real(mfile,cVar,cPar,out_string,ccParams%RLX)
    CASE ("MGCycType")
    READ(string(iEq+1:),*) param
    ccParams%CycleType = TRIM(ADJUSTL(param))
    IF (myid.eq.showid) write(mterm,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",ccParams%CycleType
    IF (myid.eq.showid) write(mfile,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",ccParams%CycleType
    CASE ("MG_VANKA")
    READ(string(iEq+1:),*) ccParams%VANKA
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%VANKA)
    CASE ("BDF")
    READ(string(iEq+1:),*) ccParams%BDF
    call write_param_int(mfile,cVar,cPar,out_string,ccParams%BDF)


  END SELECT

  END IF
 END IF
END DO

! IF (myid.eq.showid) write(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
! IF (myid.eq.showid) write(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

CLOSE (myFile)

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

END SUBROUTINE Init_CCParam

END MODULE Transport_CC


