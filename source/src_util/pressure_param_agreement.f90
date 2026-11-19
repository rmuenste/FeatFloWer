!=========================================================================
! pressure_param_agreement.f90 - collective agreement on reloaded pressure parameters
!
! A runtime reload (ProcessControl Reload_Pres) parses the parameter file on every rank. If the
! ranks read different contents (the file edited while it is read, a stale copy on one node),
! they would run the pressure multigrid with different settings: a different number of smoothing
! steps or cycles changes the number of exchanges (GetParPressure) of each rank and the run
! mismatches or deadlocks. PRESSURE_PARAMETERS_AGREE compares, after the parse, every setting
! of the pressure parameters that the parse produces (GetPresParameters: the multigrid input
! and the two keys of the periodic coarse operator; tMGParamOut holds results) over all ranks:
! integers and characters (as character codes) by an allreduce of minimum and maximum, reals by
! the same on their values (equal only when bitwise equal; a NaN makes them unequal). Equal on
! every rank: the parsed values stay, nothing else happens (a consistent reload behaves as
! before). Otherwise every rank restores the previous values (prm_old) and agreed = .FALSE.;
! the caller logs it. Collective over MPI_COMM_WORLD; called only for a reload.
!=========================================================================
MODULE pressure_param_agreement

  USE PP3D_MPI, ONLY: MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_MAX
  USE types, ONLY: tParamP

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: NI = 10, NR = 5

  PUBLIC :: PRESSURE_PARAMETERS_AGREE, PRESSURE_PARAMETERS_PACK

CONTAINS

  ! the settings GetPresParameters parses (integers, the cycle type as its character code, and
  ! reals). Fields it does not set are left out: MaxLev and MaxDifLev are derived from NLMAX,
  ! which differs between the master and the workers by design.
  SUBROUTINE PRESSURE_PARAMETERS_PACK(prm, iv, rv)
    TYPE(tParamP), INTENT(IN) :: prm
    INTEGER, INTENT(OUT) :: iv(NI)
    REAL*8, INTENT(OUT) :: rv(NR)
    iv(1) = prm%MGprmIn%MinLev
    iv(2) = prm%MGprmIn%MedLev
    iv(3) = prm%MGprmIn%MinIterCycle
    iv(4) = prm%MGprmIn%MaxIterCycle
    iv(5) = prm%MGprmIn%nSmootherSteps
    iv(6) = prm%MGprmIn%nIterCoarse
    iv(7) = prm%MGprmIn%CrsSolverType
    iv(8) = prm%PeriodicCoarseOperator
    iv(9) = prm%PeriodicCoarseCheck
    iv(10) = ICHAR(prm%MGprmIn%CycleType)
    rv(1) = prm%MGprmIn%DefImprCoarse
    rv(2) = prm%MGprmIn%CrsRelaxPrm
    rv(3) = prm%MGprmIn%Criterion1
    rv(4) = prm%MGprmIn%Criterion2
    rv(5) = prm%MGprmIn%RLX
  END SUBROUTINE PRESSURE_PARAMETERS_PACK

  ! see the module header
  SUBROUTINE PRESSURE_PARAMETERS_AGREE(prm, prm_old, agreed)
    TYPE(tParamP), INTENT(INOUT) :: prm
    TYPE(tParamP), INTENT(IN) :: prm_old
    LOGICAL, INTENT(OUT) :: agreed
    INTEGER :: iv(NI), imin(NI), imax(NI), ierr, k
    REAL*8 :: rv(NR), rmin(NR), rmax(NR)
    CALL PRESSURE_PARAMETERS_PACK(prm, iv, rv)
    CALL MPI_Allreduce(iv, imin, NI, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(iv, imax, NI, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(rv, rmin, NR, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(rv, rmax, NR, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    agreed = ALL(imin.EQ.imax)
    DO k = 1, NR
      ! equal on every rank: min = max = own value (a NaN anywhere fails this on its rank,
      ! and the rank-local results are combined below)
      IF (.NOT.(rmin(k).EQ.rmax(k) .AND. rv(k).EQ.rmin(k))) agreed = .FALSE.
    END DO
    CALL agree_logical(agreed)
    IF (.NOT.agreed) THEN
      prm%MGprmIn = prm_old%MGprmIn
      prm%PeriodicCoarseOperator = prm_old%PeriodicCoarseOperator
      prm%PeriodicCoarseCheck = prm_old%PeriodicCoarseCheck
    END IF
  END SUBROUTINE PRESSURE_PARAMETERS_AGREE

  ! logical AND over all ranks (as integers: no logical reduction type assumptions)
  SUBROUTINE agree_logical(flag)
    LOGICAL, INTENT(INOUT) :: flag
    INTEGER :: i, imin, ierr
    i = MERGE(1, 0, flag)
    CALL MPI_Allreduce(i, imin, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    flag = (imin.EQ.1)
  END SUBROUTINE agree_logical

END MODULE pressure_param_agreement
