PROGRAM TEST_PCO_RELOAD

  ! Collective agreement on reloaded pressure parameters (src_util/pressure_param_agreement.f90),
  ! run on several ranks (ctest pco-reload-mpi-3). A file that differs between the ranks cannot
  ! be staged on a shared file system, so the test hook is the parsed parameter set itself: every
  ! rank holds the same previous values and the same "parsed" values, then the parsed
  ! value is perturbed (an integer: the number of smoothing steps; a real: the coarse
  ! relaxation; the cycle type character).
  !   1 equal parsed values on every rank: kept, agreed
  !   2 nSmootherSteps 2 on the last rank, 1 elsewhere: previous values on EVERY rank
  !   3 CrsRelaxPrm differing in the last bit on one rank: previous values everywhere
  !   4 CycleType 'W' on one rank, 'V' elsewhere: previous values everywhere
  !   5 a NaN in a real setting on every rank: not agreed, previous values

  USE MPI
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
  USE types, ONLY: tParamP
  USE pressure_param_agreement, ONLY: PRESSURE_PARAMETERS_AGREE
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, EL_TEST_CHECK, EL_TEST_REPORT

  IMPLICIT NONE

  INTEGER :: rank, nproc, ierr, nf, last
  TYPE(tParamP) :: old, new, ref
  LOGICAL :: agreed

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  CALL EL_TEST_INIT()
  nf = 0
  last = nproc - 1

  CALL base(old)
  CALL base(ref)
  ref%MGprmIn%MaxIterCycle = 40
  ref%MGprmIn%CrsSolverType = 2

  ! 1
  new = ref
  CALL PRESSURE_PARAMETERS_AGREE(new, old, agreed)
  CALL EL_TEST_CHECK(agreed .AND. same(new, ref), '1 equal parsed values: kept on every rank', nf, rank)

  ! 2
  new = ref
  IF (rank.EQ.last) new%MGprmIn%nSmootherSteps = 2
  CALL PRESSURE_PARAMETERS_AGREE(new, old, agreed)
  CALL EL_TEST_CHECK(.NOT.agreed .AND. same(new, old), &
    '2 smoothing steps differ on one rank: previous values on every rank', nf, rank)

  ! 3
  new = ref
  IF (rank.EQ.0) new%MGprmIn%CrsRelaxPrm = NEAREST(new%MGprmIn%CrsRelaxPrm, 1d0)
  CALL PRESSURE_PARAMETERS_AGREE(new, old, agreed)
  CALL EL_TEST_CHECK(.NOT.agreed .AND. same(new, old), &
    '3 a real setting differs in the last bit on one rank: previous values', nf, rank)

  ! 4
  new = ref
  IF (rank.EQ.last) new%MGprmIn%CycleType = 'W'
  CALL PRESSURE_PARAMETERS_AGREE(new, old, agreed)
  CALL EL_TEST_CHECK(.NOT.agreed .AND. same(new, old), &
    '4 the cycle type differs on one rank: previous values', nf, rank)

  ! 5
  new = ref
  new%MGprmIn%Criterion1 = IEEE_VALUE(1d0, IEEE_QUIET_NAN)
  CALL PRESSURE_PARAMETERS_AGREE(new, old, agreed)
  CALL EL_TEST_CHECK(.NOT.agreed .AND. same(new, old), '5 NaN in a setting: not agreed, previous values', &
    nf, rank)

  CALL EL_TEST_REPORT(nf, rank, nproc, 'pco reload agreement')
  CALL MPI_Finalize(ierr)
  IF (nf.NE.0) STOP 1

CONTAINS

  SUBROUTINE base(p)
    TYPE(tParamP), INTENT(OUT) :: p
    p%MGprmIn%MinLev = 1
    p%MGprmIn%MedLev = 1
    p%MGprmIn%MaxLev = 3 + rank          ! rank-dependent by design (NLMAX): not compared
    p%MGprmIn%MaxDifLev = 1
    p%MGprmIn%MinIterCycle = 1
    p%MGprmIn%MaxIterCycle = 20
    p%MGprmIn%nSmootherSteps = 1
    p%MGprmIn%nIterCoarse = 50
    p%MGprmIn%CrsSolverType = 4
    p%MGprmIn%SmootherType = 1
    p%MGprmIn%vanka = 0
    p%MGprmIn%DefImprCoarse = 1d-3
    p%MGprmIn%Criterion1 = 1d-8
    p%MGprmIn%Criterion2 = 1d-1
    p%MGprmIn%RLX = 0.66d0
    p%MGprmIn%CrsRelaxPrm = 2d0/3d0
    p%MGprmIn%CrsRelaxParPrm = 1d0/3d0
    p%MGprmIn%CycleType = 'V'
    p%MGprmIn%MGProlongation = 'standard'
    p%PeriodicCoarseOperator = 0
    p%PeriodicCoarseCheck = 0
  END SUBROUTINE base

  LOGICAL FUNCTION same(a, b)
    TYPE(tParamP), INTENT(IN) :: a, b
    same = a%MGprmIn%MinLev.EQ.b%MGprmIn%MinLev .AND. a%MGprmIn%MedLev.EQ.b%MGprmIn%MedLev .AND. &
           a%MGprmIn%MaxLev.EQ.b%MGprmIn%MaxLev .AND. a%MGprmIn%MaxIterCycle.EQ.b%MGprmIn%MaxIterCycle .AND. &
           a%MGprmIn%nSmootherSteps.EQ.b%MGprmIn%nSmootherSteps .AND. &
           a%MGprmIn%CrsSolverType.EQ.b%MGprmIn%CrsSolverType .AND. &
           a%MGprmIn%CrsRelaxPrm.EQ.b%MGprmIn%CrsRelaxPrm .AND. a%MGprmIn%CycleType.EQ.b%MGprmIn%CycleType .AND. &
           a%MGprmIn%Criterion1.EQ.b%MGprmIn%Criterion1
  END FUNCTION same

END PROGRAM TEST_PCO_RELOAD
