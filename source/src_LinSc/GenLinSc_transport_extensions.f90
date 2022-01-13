SUBROUTINE SetupTransportOperators_Q1_Melt()
USE PP3D_MPI, ONLY : myid,master,showid

IF (myid.ne.0) NLMAX = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev

! Generate the necessary operators
IF (myid.ne.0) THEN

  ! Recover the physical parameters
  CALL Create_GenLinSc_Q1_PhysParamCoeffs(GenLinScalar,mg_RhoCoeff,'Rho')
  CALL Create_GenLinSc_Q1_PhysParamCoeffs(GenLinScalar,mg_CpCoeff,'Cp')
  CALL Create_GenLinSc_Q1_PhysParamCoeffs(GenLinScalar,mg_LambdaCoeff,'Lambda')
  
  ! Convection + stabilization
  CALL Create_GenLinSc_Q1_RhoCpConvection()
  CALL InitAFC_GenLinSc_Q1()

  CALL Create_GenLinSc_Q1_RhoCpMass()

  CALL Create_GenLinSc_Q1_Heat_LambdaDiffusion()

END IF

IF (myid.ne.0) NLMAX = NLMAX - GenLinScalar%prm%MGprmIn%MaxDifLev

Temperature = GenLinScalar%Fld(1)%val

END SUBROUTINE SetupTransportOperators_Q1_Melt
!
! ----------------------------------------------
!
SUBROUTINE Transport_GenLinSc_Q1_Melt(mfile,INL)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
INTEGER mfile,INL
REAL*8, allocatable :: ResTemp(:),DefTemp(:),DefTempCrit(:),RhsTemp(:)
REAL*8 tstep_old,thstep_old,defXXX
INTEGER INLComplete,I,J,iFld,iEnd,iStart
logical bDefTemp
logical :: bInit=.true.

if (.not.allocated(ResTemp)) allocate(ResTemp(GenLinScalar%nOfFields))
if (.not.allocated(DefTemp)) allocate(DefTemp(GenLinScalar%nOfFields))
if (.not.allocated(DefTempCrit)) allocate(DefTempCrit(GenLinScalar%nOfFields))
if (.not.allocated(RhsTemp)) allocate(RhsTemp(GenLinScalar%nOfFields))

IF (myid.ne.0) NLMAX = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev

! pause
thstep = tstep*(1d0-theta)

IF (myid.ne.0) THEN
 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_GenLinSc_HEATALPHA_Q1_Val(mg_mesh%level(nlmax)%dcorvg)

 DO iFld = 1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%def = 0d0
 END DO

 CALL Matdef_HEATALPHA_GenLinSc_Q1(GenLinScalar,1,0)
END IF

CALL Add_DissipativeEnergy(mfile)

IF (myid.ne.0) THEN
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%rhs = GenLinScalar%Fld(iFld)%def
 END DO
END IF

thstep = tstep*(theta)
IF (myid.ne.0) THEN

 CALL Matdef_HEATALPHA_GenLinSc_Q1(GenLinScalar,0,1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_GenLinSc_Q1_Def()

 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%aux = GenLinScalar%Fld(iFld)%def
  CALL E011Sum(GenLinScalar%fld(iFld)%def)
 END DO

 ! Save the old solution
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%val_old = GenLinScalar%Fld(iFld)%val
 END DO

 !! Compute the defect
 CALL Resdfk_GenLinSc_Q1(GenLinScalar,ResTemp,DefTemp,RhsTemp)

end if

DO iFld=1,GenLinScalar%nOfFields
 CALL COMM_Maximum(RhsTemp(iFld))
 CALL COMM_Maximum(DefTemp(iFld))
END DO

DO iFld=1,GenLinScalar%nOfFields
 DefTempCrit(iFld)=MAX((DefTemp(iFld))*GenLinScalar%prm%defCrit,GenLinScalar%prm%MinDef)
END DO

CALL Protocol_GenLinSc_Q1(mfile,GenLinScalar,0,&
       ResTemp,DefTemp,DefTempCrit," Heat&Alpha equation ")

do INL=1,GenLinScalar%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_GenLinSc_Q1_MGLinScalar(GenLinScalar,&
                                   Boundary_GenLinSc_Q1_Mat,&
                                   mfile)

IF (myid.ne.0) THEN

 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_GenLinSc_HEATALPHA_Q1_Val(mg_mesh%level(nlmax)%dcorvg)
 
 ! Restore the constant right hand side
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%fld(iFld)%def = GenLinScalar%fld(iFld)%rhs
 END DO

! Assemble the defect vector and fine level matrix
 CALL Matdef_HEATALPHA_GenLinSc_Q1(GenLinScalar,-1,0)
 
 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_GenLinSc_Q1_Def()

 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%aux = GenLinScalar%Fld(iFld)%def
  CALL E011Sum(GenLinScalar%fld(iFld)%def)
 END DO

! Save the old solution
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%val_old = GenLinScalar%Fld(iFld)%val
 END DO

! Compute the defect
 CALL Resdfk_GenLinSc_Q1(GenLinScalar,ResTemp,DefTemp,RhsTemp)

END IF

! Checking convergence rates against criterions
RhsTemp=DefTemp
DO iFld=1,GenLinScalar%nOfFields
 CALL COMM_Maximum(RhsTemp(iFld))
 CALL COMM_Maximum(DefTemp(iFld))
END DO
CALL Protocol_GenLinSc_Q1(mfile,GenLinScalar,INL,&
       ResTemp,DefTemp,DefTempCrit," GenLinScalar equation ")

bDefTemp = .True.
DO iFld=1,GenLinScalar%nOfFields
 bDefTemp = (bDefTemp.and.DefTemp(iFld).LE.DefTempCrit(iFld))
END DO

IF (bDefTemp.and.(INL.GE.GenLinScalar%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
IF (INLComplete.eq.1) GOTO 1

END DO

1 CONTINUE

2 CONTINUE

IF (myid.ne.0) NLMAX = NLMAX - GenLinScalar%prm%MGprmIn%MaxDifLev

IF (myid.ne.0) THEN
 iStart = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev
 iEnd   = NLMAX
 do NLMAX = iStart,iEnd
  CALL ProlongateSolution_GenLinSc_Q1()
 END DO
 NLMAX = iEnd
END IF

Temperature = GenLinScalar%Fld(1)%val

END SUBROUTINE Transport_GenLinSc_Q1_Melt
!
! ----------------------------------------------
!
SUBROUTINE Transport_GenLinSc_Q1_MultiMat(mfile,INL)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
INTEGER mfile,INL
REAL*8, allocatable :: ResTemp(:),DefTemp(:),DefTempCrit(:),RhsTemp(:)
REAL*8 tstep_old,thstep_old,defXXX
INTEGER INLComplete,I,J,iFld,iEnd,iStart
logical bDefTemp
logical :: bInit=.true.

if (.not.allocated(ResTemp)) allocate(ResTemp(GenLinScalar%nOfFields))
if (.not.allocated(DefTemp)) allocate(DefTemp(GenLinScalar%nOfFields))
if (.not.allocated(DefTempCrit)) allocate(DefTempCrit(GenLinScalar%nOfFields))
if (.not.allocated(RhsTemp)) allocate(RhsTemp(GenLinScalar%nOfFields))

IF (myid.ne.0) NLMAX = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev

! Generate the necessary operators
IF (myid.ne.0) THEN

 if (bInit) then
  ! Convection + stabilization
  CALL Create_GenLinSc_Q1_Convection()
  CALL InitAFC_GenLinSc_Q1()

  ! Mass Matrix scales with factor 1.0 because the equation is divided by Rho*Cp
  CALL Create_GenLinSc_Q1_Mass(1d0)

  ! Diffusion Matrix for the alpha fields
  CALL Create_GenLinSc_Q1_Alpha_Diffusion(3d-4)
  bInit = .false.
 end if
 
 ! Diffusion Matrix scales with Lambda/Rho*Cp ==> which is material specific
 CALL Create_GenLinSc_Q1_DiffCoeff(GenLinScalar)
 CALL Create_GenLinSc_Q1_Heat_Diffusion()

END IF

thstep = tstep*(1d0-theta)

IF (myid.ne.0) THEN
 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_GenLinSc_HEATALPHA_Q1_Val(mg_mesh%level(nlmax)%dcorvg)

 DO iFld = 1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%def = 0d0
 END DO

 CALL Matdef_HEATALPHA_GenLinSc_Q1(GenLinScalar,1,0)

 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%rhs = GenLinScalar%Fld(iFld)%def
 END DO
END IF

thstep = tstep*(theta)
IF (myid.ne.0) THEN

 CALL Matdef_HEATALPHA_GenLinSc_Q1(GenLinScalar,0,1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_GenLinSc_Q1_Def()

 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%aux = GenLinScalar%Fld(iFld)%def
  CALL E011Sum(GenLinScalar%fld(iFld)%def)
 END DO

 ! Save the old solution
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%val_old = GenLinScalar%Fld(iFld)%val
 END DO

 !! Compute the defect
 CALL Resdfk_GenLinSc_Q1(GenLinScalar,ResTemp,DefTemp,RhsTemp)

end if

DO iFld=1,GenLinScalar%nOfFields
 CALL COMM_Maximum(RhsTemp(iFld))
 CALL COMM_Maximum(DefTemp(iFld))
END DO

DO iFld=1,GenLinScalar%nOfFields
 DefTempCrit(iFld)=MAX((DefTemp(iFld))*GenLinScalar%prm%defCrit,GenLinScalar%prm%MinDef)
END DO

CALL Protocol_GenLinSc_Q1(mfile,GenLinScalar,0,&
       ResTemp,DefTemp,DefTempCrit," Heat&Alpha equation ")

do INL=1,GenLinScalar%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_GenLinSc_Q1_MGLinScalar(GenLinScalar,&
                                   Boundary_GenLinSc_Q1_Mat,&
                                   mfile)

IF (myid.ne.0) THEN

 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_GenLinSc_HEATALPHA_Q1_Val(mg_mesh%level(nlmax)%dcorvg)
 
 ! Restore the constant right hand side
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%fld(iFld)%def = GenLinScalar%fld(iFld)%rhs
 END DO

! Assemble the defect vector and fine level matrix
 CALL Matdef_HEATALPHA_GenLinSc_Q1(GenLinScalar,-1,0)
 
 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_GenLinSc_Q1_Def()

 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%aux = GenLinScalar%Fld(iFld)%def
  CALL E011Sum(GenLinScalar%fld(iFld)%def)
 END DO

! Save the old solution
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%val_old = GenLinScalar%Fld(iFld)%val
 END DO

! Compute the defect
 CALL Resdfk_GenLinSc_Q1(GenLinScalar,ResTemp,DefTemp,RhsTemp)

END IF

! Checking convergence rates against criterions
RhsTemp=DefTemp
DO iFld=1,GenLinScalar%nOfFields
 CALL COMM_Maximum(RhsTemp(iFld))
 CALL COMM_Maximum(DefTemp(iFld))
END DO
CALL Protocol_GenLinSc_Q1(mfile,GenLinScalar,INL,&
       ResTemp,DefTemp,DefTempCrit," GenLinScalar equation ")

bDefTemp = .True.
DO iFld=1,GenLinScalar%nOfFields
 bDefTemp = (bDefTemp.and.DefTemp(iFld).LE.DefTempCrit(iFld))
END DO

IF (bDefTemp.and.(INL.GE.GenLinScalar%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
IF (INLComplete.eq.1) GOTO 1

END DO

1 CONTINUE

2 CONTINUE

CALL CheckAlphaConvergence(mfile)

IF (myid.ne.0) NLMAX = NLMAX - GenLinScalar%prm%MGprmIn%MaxDifLev

IF (myid.ne.0) THEN
 iStart = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev
 iEnd   = NLMAX
 do NLMAX = iStart,iEnd
  CALL ProlongateSolution_GenLinSc_Q1()
 END DO
 NLMAX = iEnd
END IF


END SUBROUTINE Transport_GenLinSc_Q1_MultiMat
!
! ----------------------------------------------
!
SUBROUTINE Transport_GenLinSc_Q1(mfile,INL)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
INTEGER mfile,INL
REAL*8, allocatable :: ResTemp(:),DefTemp(:),DefTempCrit(:),RhsTemp(:)
REAL*8 tstep_old,thstep_old,defXXX
INTEGER INLComplete,I,J,iFld,iEnd,iStart

if (.not.allocated(ResTemp)) allocate(ResTemp(GenLinScalar%nOfFields))
if (.not.allocated(DefTemp)) allocate(DefTemp(GenLinScalar%nOfFields))
if (.not.allocated(DefTempCrit)) allocate(DefTempCrit(GenLinScalar%nOfFields))
if (.not.allocated(RhsTemp)) allocate(RhsTemp(GenLinScalar%nOfFields))

IF (myid.ne.0) NLMAX = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev

! Generate the necessary operators
IF (myid.ne.0) THEN

 CALL Create_GenLinSc_Q1_Convection()
 CALL InitAFC_GenLinSc_Q1()

 CALL Create_GenLinSc_Q1_Mass(1d0)

 CALL Create_GenLinSc_Q1_Diffusion(Properties%DiffCoeff(1))

END IF

thstep = tstep*(1d0-theta)

IF (myid.ne.0) THEN
 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_GenLinSc_Q1_Val(mg_mesh%level(nlmax)%dcorvg)

 DO iFld = 1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%def = 0d0
 END DO

 IF (TRIM(GenLinScalar%prm%cEquation).eq."INSTATIONARY") then
  CALL Matdef_INSTATIONARY_GenLinSc_Q1(GenLinScalar,1,0)
 END IF

 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%rhs = GenLinScalar%Fld(iFld)%def
 END DO
END IF

thstep = tstep*(theta)
IF (myid.ne.0) THEN

 IF (TRIM(GenLinScalar%prm%cEquation).eq."STATIONARY") then
  CALL Matdef_STATIONARY_GenLinSc_Q1(GenLinScalar,0,1)
 END IF
 IF (TRIM(GenLinScalar%prm%cEquation).eq."INSTATIONARY") then
  CALL Matdef_INSTATIONARY_GenLinSc_Q1(GenLinScalar,0,1)
 END IF

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_GenLinSc_Q1_Def()

 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%aux = GenLinScalar%Fld(iFld)%def
  CALL E011Sum(GenLinScalar%fld(iFld)%def)
 END DO

 ! Save the old solution
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%val_old = GenLinScalar%Fld(iFld)%val
 END DO

 !! Compute the defect
 CALL Resdfk_GenLinSc_Q1(GenLinScalar,ResTemp,DefTemp,RhsTemp)

end if

DO iFld=1,GenLinScalar%nOfFields
 CALL COMM_Maximum(RhsTemp(iFld))
END DO

DO iFld=1,GenLinScalar%nOfFields
 DefTempCrit(iFld)=MAX((RhsTemp(iFld))*GenLinScalar%prm%defCrit,GenLinScalar%prm%MinDef)
END DO

CALL Protocol_GenLinSc_Q1(mfile,GenLinScalar,0,&
       ResTemp,DefTemp,DefTempCrit," Laplace equation ")

do INL=1,GenLinScalar%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_GenLinSc_Q1_MGLinScalar(GenLinScalar,&
                                   Boundary_GenLinSc_Q1_Mat,&
                                   mfile)

IF (myid.ne.0) THEN

 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_GenLinSc_Q1_Val(mg_mesh%level(nlmax)%dcorvg)
 
 ! Restore the constant right hand side
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%fld(iFld)%def = GenLinScalar%fld(iFld)%rhs
 END DO

! Assemble the defect vector and fine level matrix
 IF (TRIM(GenLinScalar%prm%cEquation).eq."STATIONARY") THEN
  CALL Matdef_STATIONARY_GenLinSc_Q1(GenLinScalar,-1,0)
 END IF
 IF (TRIM(GenLinScalar%prm%cEquation).eq."INSTATIONARY") THEN
  CALL Matdef_INSTATIONARY_GenLinSc_Q1(GenLinScalar,-1,0)
 END IF
 
 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_GenLinSc_Q1_Def()

 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%aux = GenLinScalar%Fld(iFld)%def
  CALL E011Sum(GenLinScalar%fld(iFld)%def)
 END DO

! Save the old solution
 DO iFld=1,GenLinScalar%nOfFields
  GenLinScalar%Fld(iFld)%val_old = GenLinScalar%Fld(iFld)%val
 END DO

! Compute the defect
 CALL Resdfk_GenLinSc_Q1(GenLinScalar,ResTemp,DefTemp,RhsTemp)

END IF

! Checking convergence rates against criterions
RhsTemp=DefTemp
DO iFld=1,GenLinScalar%nOfFields
 CALL COMM_Maximum(RhsTemp(iFld))
END DO
CALL Protocol_GenLinSc_Q1(mfile,GenLinScalar,INL,&
       ResTemp,DefTemp,DefTempCrit," GenLinScalar equation ")

IF ((DefTemp(1).LE.DefTempCrit(1)).AND.&
    (DefTemp(2).LE.DefTempCrit(2)).AND.&
    (DefTemp(3).LE.DefTempCrit(3)).AND.&
    (INL.GE.GenLinScalar%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
IF (INLComplete.eq.1) GOTO 1

END DO

1 CONTINUE

2 CONTINUE

IF (myid.ne.0) NLMAX = NLMAX - GenLinScalar%prm%MGprmIn%MaxDifLev

IF (myid.ne.0) THEN
 iStart = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev
 iEnd   = NLMAX
 do NLMAX = iStart,iEnd
  CALL ProlongateSolution_GenLinSc_Q1()
 END DO
 NLMAX = iEnd
END IF

END SUBROUTINE Transport_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Init_GenLinSc_Q1(log_unit)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
USE var_QuadScalar
implicit none
integer, intent(in) :: log_unit
integer :: n, ndof, iFld

 GenLinScalar%cName = "Temper"
 CALL Get_GenLinSc_Q1_Parameters(GenLinScalar%prm,GenLinScalar%cName,log_unit)
 GenLinScalar%nOfFields = GenLinScalar%prm%nOfFields
 
 IF (myid.ne.0) NLMAX = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev
 
 ! Building up the matrix strucrures
 CALL Create_MatStruct_Q1()

! Iteration matrix (only allocation)
 CALL Create_AMat_GenLinSc_Q1(GenLinScalar)

 CALL Initialize_GenLinSc_Q1(GenLinScalar)
 
 ILEV=NLMAX
 CALL SETLEV(2)

 ! Set the types of boundary conditions (set up knpr)
 call Knpr_GenLinSc_Q1(mg_mesh%level(ilev)%dcorvg)

 DO iFld = 1,GenLinScalar%nOfFields
  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'FAC') GenLinScalar%Fld(iFld)%val = 0d0
 END DO

 !  ! Building up the matrix strucrures
 CALL Create_GenLinSc_Q1_AFCStruct()

 IF (myid.ne.0) NLMAX = NLMAX - GenLinScalar%prm%MGprmIn%MaxDifLev
 
!!=========================================================================
 IF (myid.ne.0) NLMAX = NLMAX + 1
 ! Prolongation Restriction Matrices
 ALLOCATE(mg_E011Prol(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011Rest(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011ProlM(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011RestM(NLMIN:NLMAX-1))

 DO ILEV=NLMIN,NLMAX-1
  N = KNVT(ILEV+1) + 2*KNET(ILEV+1) + 4*KNAT(ILEV+1) + 8*KNEL(ILEV+1)
  ALLOCATE(mg_E011Prol(ILEV)%a(N))
  ALLOCATE(mg_E011Rest(ILEV)%a(N))
  ALLOCATE(mg_E011ProlM(ILEV)%ColA(N))
  ALLOCATE(mg_E011RestM(ILEV)%ColA(N))
  mg_E011ProlM(ILEV)%na = N
  mg_E011RestM(ILEV)%na = N
  NDOF = KNVT(ILEV+1)
  mg_E011ProlM(ILEV)%nu = NDOF
  ALLOCATE(mg_E011ProlM(ILEV)%LdA(NDOF+1))
  NDOF = KNVT(ILEV)
  mg_E011RestM(ILEV)%nu = NDOF
  ALLOCATE(mg_E011RestM(ILEV)%LdA(NDOF+1))
 END DO

 CALL InitializeProlRest_GenLinSc_Q1(GenLinScalar,GenLinScalar%cName)
 IF (myid.ne.0) NLMAX = NLMAX - 1
!!=========================================================================

END SUBROUTINE Init_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Init_GenLinSc_HEATALPHA_Q1(log_unit)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
USE Sigma_User, ONLY: myMultiMat
USE var_QuadScalar
implicit none
integer, intent(in) :: log_unit
integer :: n, ndof, iFld,nMaterials

 GenLinScalar%cName = "Temper"
 nMaterials = myMultiMat%nOfMaterials+1
 GenLinScalar%prm%nOfFields = nMaterials
 ALLOCATE(GenLinScalar%prm%cField(GenLinScalar%prm%nOfFields))
 DO iFld=1,nMaterials
  if (iFld.eq.1) GenLinScalar%prm%cField(iFld) = 'temp'
  if (iFld.gt.1) then
   write(GenLinScalar%prm%cField(iFld),'(A,I0)') 'alpha',iFld-1
  end if
 end do
 
 CALL Get_GenLinSc_Q1_Parameters(GenLinScalar%prm,GenLinScalar%cName,log_unit)
 GenLinScalar%nOfFields = GenLinScalar%prm%nOfFields
 
 IF (myid.ne.0) NLMAX = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev
 
 ! Building up the matrix strucrures
 CALL Create_MatStruct_Q1()

! Iteration matrix (only allocation)
 CALL Create_AMat_GenLinSc_Q1(GenLinScalar)

 CALL Initialize_GenLinSc_Q1(GenLinScalar)
 DO iFld=1,nMaterials
  if (iFld.gt.1) then
   GenLinScalar%Fld(iFld)%ID = iFld-1
  end if
 end do
 
 ILEV=NLMAX
 CALL SETLEV(2)

 ! Set the types of boundary conditions (set up knpr)
 call Knpr_GenLinSc_HEATALPHA_Q1(mg_mesh%level(ilev)%dcorvg)

 DO iFld = 1,GenLinScalar%nOfFields
  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'temp')  THEN
   GenLinScalar%Fld(iFld)%val = myProcess%T0
  ELSE
   GenLinScalar%Fld(iFld)%val = 0d0
  END IF
 END DO

 !  ! Building up the matrix strucrures
 CALL Create_GenLinSc_Q1_AFCStruct()

 IF (myid.ne.0) NLMAX = NLMAX - GenLinScalar%prm%MGprmIn%MaxDifLev
 
!!=========================================================================
 IF (myid.ne.0) NLMAX = NLMAX + 1
 ! Prolongation Restriction Matrices
 ALLOCATE(mg_E011Prol(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011Rest(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011ProlM(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011RestM(NLMIN:NLMAX-1))

 DO ILEV=NLMIN,NLMAX-1
  N = KNVT(ILEV+1) + 2*KNET(ILEV+1) + 4*KNAT(ILEV+1) + 8*KNEL(ILEV+1)
  ALLOCATE(mg_E011Prol(ILEV)%a(N))
  ALLOCATE(mg_E011Rest(ILEV)%a(N))
  ALLOCATE(mg_E011ProlM(ILEV)%ColA(N))
  ALLOCATE(mg_E011RestM(ILEV)%ColA(N))
  mg_E011ProlM(ILEV)%na = N
  mg_E011RestM(ILEV)%na = N
  NDOF = KNVT(ILEV+1)
  mg_E011ProlM(ILEV)%nu = NDOF
  ALLOCATE(mg_E011ProlM(ILEV)%LdA(NDOF+1))
  NDOF = KNVT(ILEV)
  mg_E011RestM(ILEV)%nu = NDOF
  ALLOCATE(mg_E011RestM(ILEV)%LdA(NDOF+1))
 END DO

 CALL InitializeProlRest_GenLinSc_Q1(GenLinScalar,GenLinScalar%cName)
 IF (myid.ne.0) NLMAX = NLMAX - 1
!!=========================================================================

END SUBROUTINE Init_GenLinSc_HEATALPHA_Q1
!
! ----------------------------------------------
!
SUBROUTINE InitializeProlRest_GenLinSc_Q1(lSc,cF)
implicit none
TYPE(lScalarGen), INTENT(INOUT), TARGET :: lSc
CHARACTER cF*(*)

 MyMG%MinLev  = NLMIN
 myMG%MedLev  = NLMIN
 myMG%MaxLev  = NLMAX

 IF(myid.eq.showid) WRITE(*,*) "Initialization of "//ADJUSTL(TRIM(cF))//" prolongation matrix"
 myMG%bProlRest => lSc%bProlRest
 MyMG%cVariable = ADJUSTL(TRIM(cF))
 CALL mgProlRestInit()

END SUBROUTINE InitializeProlRest_GenLinSc_Q1  
!
! ----------------------------------------------
!
SUBROUTINE Create_GenLinSc_Q1_RhoCpConvection()
INTEGER I,J,jLEV
EXTERNAL E011
integer invt, inet, inat, inel

 IF (.not.allocated(mg_ConvMat)) allocate(mg_ConvMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX-1
   JLEV = ILEV-1

   CALL SETLEV(2)
   IF (.not.allocated(mg_ConvMat(ILEV)%a)) allocate(mg_ConvMat(ILEV)%a(mg_lMat(ILEV)%na))

   IF (myid.eq.showID) THEN
    IF (ILEV.EQ.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') " [K]: [", ILEV,"]"
    ELSE
     WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
    END IF
   END IF
   
   ConvectionMat=>mg_ConvMat(ILEV)%a
   plMat=>mg_lMat(ILEV)

   ConvectionMat = 0d0

   CALL RhoCpConv_LinSc1(QuadSc%valU,QuadSc%valV,QuadSc%valW,ConvectionMat,&
   mg_RhoCoeff(ILEV)%x,mg_CpCoeff(ILEV)%x,&
   plMat%nu,plMat%ColA,plMat%LdA,&
   mg_mesh%level(ilev)%kvert,&
   mg_mesh%level(ilev)%karea,&
   mg_mesh%level(ilev)%kedge,&
   mg_mesh%level(ilev)%kedge,&
   mg_mesh%level(ilev)%dcorvg,&
   E011)
   
 END DO

 ILEV=NLMAX
 JLEV = ILEV-1
 CALL SETLEV(2)
 IF (.not.allocated(mg_ConvMat(ILEV)%a)) allocate(mg_ConvMat(ILEV)%a(mg_lMat(ILEV)%na))
 ConvectionMat=>mg_ConvMat(ILEV)%a
 plMat=>mg_lMat(ILEV)
 ConvectionMat = 0d0
 IF (myid.eq.showID) THEN
  IF (ILEV.EQ.NLMIN) THEN
   WRITE(MTERM,'(A,I1,A)', advance='no') " [K]: [", ILEV,"]"
  ELSE
   WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
  END IF
 END IF
 CALL RhoCpConv_LinSc2(QuadSc%valU,QuadSc%valV,QuadSc%valW,ConvectionMat,&
 mg_RhoCoeff(ILEV)%x,mg_CpCoeff(ILEV)%x,&
 plMat%nu,plMat%ColA,plMat%LdA,&
 mg_mesh%level(ilev)%kvert,&
 mg_mesh%level(ilev)%karea,&
 mg_mesh%level(ilev)%kedge,&
 mg_mesh%level(ilev)%dcorvg,&
 mg_mesh%level(ilev)%kadj,&
 mg_mesh%level(jlev)%kvert,&
 mg_mesh%level(jlev)%karea,&
 mg_mesh%level(jlev)%kedge,&
 mg_mesh%level(jlev)%nel,&
 mg_mesh%level(jlev)%nvt,&
 mg_mesh%level(jlev)%net,&
 mg_mesh%level(jlev)%nat,E011)

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') " |"
 
END SUBROUTINE Create_GenLinSc_Q1_RhoCpConvection
!
! ----------------------------------------------
!
SUBROUTINE Create_GenLinSc_Q1_Convection()
INTEGER I,J,jLEV
EXTERNAL E011
integer invt, inet, inat, inel

 IF (.not.allocated(mg_ConvMat)) allocate(mg_ConvMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX-1
   JLEV = ILEV-1

   CALL SETLEV(2)
   IF (.not.allocated(mg_ConvMat(ILEV)%a)) allocate(mg_ConvMat(ILEV)%a(mg_lMat(ILEV)%na))

   IF (myid.eq.showID) THEN
    IF (ILEV.EQ.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') " [K]: [", ILEV,"]"
    ELSE
     WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
    END IF
   END IF
   
   ConvectionMat=>mg_ConvMat(ILEV)%a
   plMat=>mg_lMat(ILEV)

   ConvectionMat = 0d0

   CALL Conv_LinSc1(QuadSc%valU,QuadSc%valV,QuadSc%valW,ConvectionMat,&
   plMat%nu,plMat%ColA,plMat%LdA,&
   mg_mesh%level(ilev)%kvert,&
   mg_mesh%level(ilev)%karea,&
   mg_mesh%level(ilev)%kedge,&
   mg_mesh%level(ilev)%kedge,&
   mg_mesh%level(ilev)%dcorvg,&
   E011)
   
 END DO

 ILEV=NLMAX
 JLEV = ILEV-1
 CALL SETLEV(2)
 IF (.not.allocated(mg_ConvMat(ILEV)%a)) allocate(mg_ConvMat(ILEV)%a(mg_lMat(ILEV)%na))
 ConvectionMat=>mg_ConvMat(ILEV)%a
 plMat=>mg_lMat(ILEV)
 ConvectionMat = 0d0
 IF (myid.eq.showID) THEN
  IF (ILEV.EQ.NLMIN) THEN
   WRITE(MTERM,'(A,I1,A)', advance='no') " [K]: [", ILEV,"]"
  ELSE
   WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
  END IF
 END IF
 CALL Conv_LinSc2(QuadSc%valU,QuadSc%valV,QuadSc%valW,ConvectionMat,&
 plMat%nu,plMat%ColA,plMat%LdA,&
 mg_mesh%level(ilev)%kvert,&
 mg_mesh%level(ilev)%karea,&
 mg_mesh%level(ilev)%kedge,&
 mg_mesh%level(ilev)%dcorvg,&
 mg_mesh%level(ilev)%kadj,&
 mg_mesh%level(jlev)%kvert,&
 mg_mesh%level(jlev)%karea,&
 mg_mesh%level(jlev)%kedge,&
 mg_mesh%level(jlev)%nel,&
 mg_mesh%level(jlev)%nvt,&
 mg_mesh%level(jlev)%net,&
 mg_mesh%level(jlev)%nat,E011)

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') " |"
 
END SUBROUTINE Create_GenLinSc_Q1_Convection
!
! ----------------------------------------------
!
SUBROUTINE ProlongateSolution_GenLinSc_Q1()

 CALL ProlongateSolution_GenLinSc_Q1_sub(GenLinScalar) 
 
END SUBROUTINE ProlongateSolution_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Reinit_GenLinSc_Q1()
integer iFld,i

do iFld=1,GenLinScalar%nOfFields

 IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'temp') then
!   GenLinScalar%Fld(iFld)%val = 0.5d0
  DO i=1,GenLinScalar%ndof
   IF (GenLinScalar%Fld(iFld)%val(i).gt.0.5d0) then
    GenLinScalar%Fld(iFld)%val(i) = 1d0
   else
    GenLinScalar%Fld(iFld)%val(i) = 0d0
   end if
  END DO
 END IF
END DO

END SUBROUTINE Reinit_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Correct_GenLinSc_Q1_ALPHA(mfile)
use Sigma_User, only: mySegmentIndicator
INTEGER mfile
integer iFld,i,iMat,iFldMax,iSeg
integer iStart,iEnd
real*8 dFldMax,dMaxValue

IF (myid.ne.0) THEN

if (myid.eq.1) then
 write(mterm,*) 'Reinitialization of the fields takes place now!'
 write(mfile,*) 'Reinitialization of the fields takes place now!'
end if

NLMAX = NLMAX + 1! GenLinScalar%prm%MGprmIn%MaxDifLev

DO i=1,mg_mesh%level(nlmax)%nvt
 dMaxValue = 1d0
!  iSeg = INT(mySegmentIndicator(2,i))
!  IF (iSeg.eq.0) dMaxValue = 0d0
  
 iMat =0
 dFldMax = -1d0
 do iFld=2,GenLinScalar%nOfFields
  if (GenLinScalar%Fld(iFld)%val(i).gt.0.5d0)  then
   iMat = iFld - 1 
  else
   if (GenLinScalar%Fld(iFld)%val(i).gt.dFldMax)  then
    iFldMax = iFld-1
    dFldMax = GenLinScalar%Fld(iFld)%val(i)
   end if
  end if
 end do
 
 if (iMat.ne.0) then ! Material with clear majority is known

  do iFld=2,GenLinScalar%nOfFields
   if (iMat.eq.iFld-1) then
    GenLinScalar%Fld(iFld)%val(i) = dMaxValue
   else
    GenLinScalar%Fld(iFld)%val(i) = 0d0
   end if
  end do
  
 else               ! Material with largest fraction is used
 
  do iFld=2,GenLinScalar%nOfFields
   if (iFldMax.eq.iFld-1) then
    GenLinScalar%Fld(iFld)%val(i) = dMaxValue
   else
    GenLinScalar%Fld(iFld)%val(i) = 0d0
   end if
  end do
  
 end if

END DO !i

NLMAX = NLMAX - 1 !GenLinScalar%prm%MGprmIn%MaxDifLev

! iStart = NLMAX + GenLinScalar%prm%MGprmIn%MaxDifLev
! iEnd   = NLMAX
! do NLMAX = iStart,iEnd
!  CALL ProlongateSolution_GenLinSc_Q1()
! END DO
! NLMAX = iEnd

END IF

END SUBROUTINE Correct_GenLinSc_Q1_ALPHA
!
! ----------------------------------------------
!
SUBROUTINE CheckAlphaConvergence(mfile)
use Sigma_User, only: mySegmentIndicator
USE var_QuadScalar, ONLY: bAlphaConverged
integer mfile
integer i,iPhase,ifld,iSeg
real*8 dMaxValue
real*8, allocatable :: dVolPhase(:)
real*8 :: dAlphaThreshold = 0.3333d0

allocate(dVolPhase(GenLinScalar%nOfFields))
dVolPhase = 0d0

bAlphaConverged=.false.

IF (myid.ne.0) THEN

  LMassMat   => mg_LMassMat(nlmax)%a
  
  DO i=1,mg_mesh%level(nlmax)%nvt

  iSeg = INT(mySegmentIndicator(2,i))

  IF (iSeg.ne.0) THEN
   iPhase = 0
   dMaxValue = 0d0
   
   do iFld=2,GenLinScalar%nOfFields
    if (GenLinScalar%Fld(iFld)%val(i).gt.dAlphaThreshold.and.GenLinScalar%Fld(iFld)%val(i).gt.dMaxValue) then
     dMaxValue = GenLinScalar%Fld(iFld)%val(i)
     iPhase = iFld
    end if
   end do
   
   if (iPhase.ne.0) then
    dVolPhase(iPhase) = dVolPhase(iPhase) + LMassMat(i)
   end if
   dVolPhase(1) = dVolPhase(1) + lMassMat(i)
   
  END IF
  
 END DO

END IF

CALL COMM_SUMMN(dVolPhase,GenLinScalar%nOfFields)

dMaxValue=0d0
DO iFld=2,GenLinScalar%nOfFields
 dMaxValue = dMaxValue + dVolPhase(iFld)
end do

if (dMaxValue/dVolPhase(1).gt.0.9d0) bAlphaConverged=.true.

if (myid.eq.1) then
 write(mterm,'(A,100ES12.4)') 'VolumeFractions[%]: ',1d2*dMaxValue/dVolPhase(1),1d2*dVolPhase(2:GenLinScalar%nOfFields)/dVolPhase(1)
 write(mfile,'(A,100ES12.4)') 'VolumeFractions[%]: ',1d2*dMaxValue/dVolPhase(1),1d2*dVolPhase(2:GenLinScalar%nOfFields)/dVolPhase(1)
end if

deallocate(dVolPhase)

END SUBROUTINE CheckAlphaConvergence
!
! ----------------------------------------------
!
SUBROUTINE EstimateAlphaTimeStepSize(mfile)
use Sigma_User, only: mySegmentIndicator,myProcess,myThermodyn
USE var_QuadScalar, ONLY: bAlphaConverged
integer mfile
integer iseg,iInflow,i,iSubInflow
real*8 dVolume,dFlow,ResidenceTime

dVolume = 0d0

IF (myid.ne.0) THEN

 ! Mass Matrix scales with factor 1.0 because the equation is divided by Rho*Cp
 CALL Create_GenLinSc_Q1_Mass(1d0)

 LMassMat   => mg_LMassMat(nlmax)%a
  
 DO i=1,mg_mesh%level(nlmax)%nvt
  iSeg = INT(mySegmentIndicator(2,i))
  IF (iSeg.ne.0) THEN
   dVolume = dVolume + lMassMat(i) ! cm3
  end if
 END DO
end if

CALL COMM_SUMM(dVolume)

dFlow = 0d0
do iInflow=1,myProcess%nOfInflows
 IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
  dFlow = dFlow + myProcess%myInflow(iInflow)%massflowrate  
 ELSE
  DO iSubInflow=1,myProcess%myInflow(iInflow)%nSubInflows
   dFlow = dFlow + myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate  
  END DO
 END IF
enddo
dFlow = (dFlow/3600d0)/(1e3*myThermodyn%density) ! (kg/s)/(kg/m3) = m3/s
dFlow = 1e6*dFlow ! cm3/s

ResidenceTime = dVolume/dFlow
tstep = ResidenceTime/600d0
timemx = 5000d0*tstep
dtgmv = 250d0*tstep

if (myid.eq.1) then
 write(mterm,'(A,100ES12.4)') 'TimeEstimations[Vol,Volrate,RT,dt]: ',dVolume,dFlow,ResidenceTime,tstep
 write(mfile,'(A,100ES12.4)') 'TimeEstimations[Vol,Volrate,RT,dt]: ',dVolume,dFlow,ResidenceTime,tstep
end if

END SUBROUTINE EstimateAlphaTimeStepSize
!
! ----------------------------------------------
!
SUBROUTINE Add_DissipativeEnergy(mfile)
USE var_QuadScalar, ONLY : Viscosity,Shearrate
integer mfile
integer j
real*8 daux,dauxKW,dDissipEnergy,Rho,Cp,Temp,MF,dMasMat,dFactor,dShearSquare
real*8 :: T_Fin = 407.00d0-273.15d0
real*8 :: AlphaViscosityMatModel

dDissipEnergy = 0d0

if (myid.ne.0) then 
 ILEV = NLMAX
 LMassMat      => mg_LMassMat(ILEV)%a

 DO j=1,plMat%nu
   Temp = GenLinScalar%Fld(1)%val(j)
   CALL MeltFunction_MF(MF,Temp)
   CALL MeltFunction_Rho(Rho,Temp,MF)
   CALL MeltFunction_Cp(Cp,Temp,MF)
   Rho = 1d-6*1d-3*Rho ! kg/m3 = 1d-3T/1d6cm3 = 1d-9 T/cm3
   Cp  = 1d4*Cp        !J/kg/K = m2/(s2*K) = 1d+4cm2/(s2*K)
   
   IF (Temp.lt.T_Fin) then
    Temp = T_Fin
    dShearSquare = (0.5d0*Shearrate(j))**2d0
    viscosity(j) = AlphaViscosityMatModel(dShearSquare,1,Temp)
    dFactor = MF*Viscosity(j)*Shearrate(j)**2d0
   else
    dFactor = 1d0*Viscosity(j)*Shearrate(j)**2d0
   end if
   
   dMasMat = LMassMat(j)/(Rho*Cp) ! cm3
   
   dauxKW = 1d-7*dFactor*dMasMat ! g/(cm*s)*1/(s2)*cm3 ==> 1d-3*kg*1e-4m2/(s3) = 1d-7 kg/(m2*s3) = 1d-7*W 
   daux   = 1d-6*dFactor*dMasMat ! g/(cm*s)*1/(s2)*cm3 ==> 1d-6*T*cm2/(s3) = 1d-6 T/(cm2*s3) = 
   
   
   GenLinScalar%Fld(1)%def(j) = GenLinScalar%fld(1)%def(j) + 1d0*daux*tstep
   dDissipEnergy = dDissipEnergy + 1d-3*dauxKW ! W = 1d-3 kW
 END DO
END IF

CALL COMM_SUMM(dDissipEnergy)

IF (myid.eq.1) THEN
 write(mterm,'(A,100ES12.4)') 'DissipationEnergy[kW]: ',dDissipEnergy
 write(mfile,'(A,100ES12.4)') 'DissipationEnergy[kW]: ',dDissipEnergy
END IF

END SUBROUTINE Add_DissipativeEnergy

