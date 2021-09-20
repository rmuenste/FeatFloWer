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
 CALL Boundary_GenLinSc_Q1_Val(mg_mesh%level(ilev)%dcorvg)

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
 CALL Boundary_GenLinSc_Q1_Val(mg_mesh%level(ilev)%dcorvg)
 
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

! DO iFld = 1,GenLinScalar%nOfFields
  GenLinScalar%Fld(1)%val = 250d0
  GenLinScalar%Fld(2)%val = 1d0
! END DO

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
