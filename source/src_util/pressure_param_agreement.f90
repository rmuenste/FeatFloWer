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
! the same on their bit patterns (equal only when bitwise equal, +0 and -0 differ; a NaN is
! rejected). Equal on
! every rank: the parsed values stay, nothing else happens (a consistent reload behaves as
! before). Otherwise every rank restores the previous values (prm_old) and agreed = .FALSE.;
! the caller logs it. Collective over MPI_COMM_WORLD; called only for a reload.
!=========================================================================
MODULE pressure_param_agreement

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
  USE PP3D_MPI, ONLY: MPI_COMM_WORLD, MPI_INTEGER, MPI_INTEGER8, MPI_MIN, MPI_MAX
  USE types, ONLY: tParamP

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: NI = 10, NR = 5

  PUBLIC :: PRESSURE_PARAMETERS_AGREE, PRESSURE_PARAMETERS_PACK, PCO_KEYS_AGREE

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
    INTEGER :: iv(NI), imin(NI), imax(NI), ierr
    REAL*8 :: rv(NR)
    INTEGER(8) :: bv(NR), bmin(NR), bmax(NR)
    CALL PRESSURE_PARAMETERS_PACK(prm, iv, rv)
    ! reals as bit patterns: equal only when bitwise equal (+0 and -0 differ); a NaN is
    ! rejected explicitly on the rank that holds it
    bv = TRANSFER(rv, bv)
    CALL MPI_Allreduce(iv, imin, NI, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(iv, imax, NI, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(bv, bmin, NR, MPI_INTEGER8, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(bv, bmax, NR, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, ierr)
    agreed = ALL(imin.EQ.imax) .AND. ALL(bmin.EQ.bmax) .AND. .NOT.ANY(IEEE_IS_NAN(rv))
    CALL agree_logical(agreed)
    IF (.NOT.agreed) THEN
      prm%MGprmIn = prm_old%MGprmIn
      prm%PeriodicCoarseOperator = prm_old%PeriodicCoarseOperator
      prm%PeriodicCoarseCheck = prm_old%PeriodicCoarseCheck
    END IF
  END SUBROUTINE PRESSURE_PARAMETERS_AGREE

  ! the two keys of the periodic coarse operator (Pres@PeriodicCoarseOperator,
  ! Pres@PeriodicCoarseCheck) equal on every rank? Collective; the caller stops every rank on
  ! .FALSE. (they decide the mode and the schedule of collective checks).
  LOGICAL FUNCTION PCO_KEYS_AGREE(op, chk) RESULT(ok)
    INTEGER, INTENT(IN) :: op, chk
    INTEGER :: v(2), vmin(2), vmax(2), ierr
    v = (/op, chk/)
    CALL MPI_Allreduce(v, vmin, 2, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(v, vmax, 2, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    ok = ALL(vmin.EQ.vmax)
  END FUNCTION PCO_KEYS_AGREE

  ! logical AND over all ranks (as integers: no logical reduction type assumptions)
  SUBROUTINE agree_logical(flag)
    LOGICAL, INTENT(INOUT) :: flag
    INTEGER :: i, imin, ierr
    i = MERGE(1, 0, flag)
    CALL MPI_Allreduce(i, imin, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    flag = (imin.EQ.1)
  END SUBROUTINE agree_logical

END MODULE pressure_param_agreement
