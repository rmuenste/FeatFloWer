!=========================================================================
! CHI_SPARSE_DIRECT - instance-based sparse direct solver of the Chimera
! component (design: chimera-integration-design.md v3, Layer M; v2-review
! blocker 5).
!
! Wraps the umf4* Fortran interface of UMFPACK.  Unlike the legacy
! module UMFPackSolver (module-global symbolic/numeric handles, in-place
! 1->0 index mutation of the caller's CSR), every tSparseDirectSolver
! instance owns:
!   - its 0-based, column-sorted copy of the CSR structure (the caller's
!     arrays are never touched),
!   - a value-permutation from the caller's ordering to the sorted one,
!   - its own symbolic/numeric handles.  The umf4 F77 wrapper
!     (extern/libraries/umfpack4/src/umf4_f77wrapper_port.c) stores the
!     factorization pointers in a handle table, so concurrent instances
!     are supported by design.
!
! Matrix convention: the caller passes a CSR matrix A (row pointers LdA,
! column indices ColA).  UMFPACK interprets the arrays as CSC, i.e. it
! factorizes A^T; SOLVE therefore requests the transpose system
! (sys = 1, UMFPACK_At; identical to UMFPACK_Aat for real matrices) so
! that A x = b is solved.
!=========================================================================
MODULE CHI_SPARSE_DIRECT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tSparseDirectSolver
  PUBLIC :: SD_INIT
  PUBLIC :: SD_FACTORIZE
  PUBLIC :: SD_SOLVE
  PUBLIC :: SD_FREE

  TYPE tSparseDirectSolver
    LOGICAL :: initialized = .FALSE.
    LOGICAL :: factorized = .FALSE.
    INTEGER :: n = 0, nz = 0
    INTEGER :: symbolic = -1, numeric = -1
    ! 0-based, per-row column-sorted structure copy
    INTEGER, ALLOCATABLE :: Ap0(:), Ai0(:)
    ! caller position -> sorted position of each value
    INTEGER, ALLOCATABLE :: perm(:)
    ! sorted value copy of the last factorized matrix (needed by the
    ! iterative refinement of umf4solr)
    REAL*8, ALLOCATABLE :: Ax(:)
    REAL*8 :: control(20) = 0d0
    REAL*8 :: info(90) = 0d0
  END TYPE tSparseDirectSolver

CONTAINS

  !-----------------------------------------------------------------------
  ! Store the structure (1-based CSR: LdA(n+1), ColA(nz)) as a 0-based,
  ! per-row sorted copy plus the value permutation.  No factorization
  ! happens here.
  !-----------------------------------------------------------------------
  SUBROUTINE SD_INIT(slv, n, LdA, ColA, ok)
    TYPE(tSparseDirectSolver), INTENT(INOUT) :: slv
    INTEGER, INTENT(IN) :: n, LdA(n+1), ColA(*)
    LOGICAL, INTENT(OUT) :: ok

    INTEGER :: i, j, j1, j2, m, k, key, keyp

    ok = .FALSE.
    CALL SD_FREE(slv)
    IF (n .LE. 0) RETURN

    slv%n = n
    slv%nz = LdA(n+1) - 1
    ALLOCATE(slv%Ap0(n+1), slv%Ai0(slv%nz), slv%perm(slv%nz), slv%Ax(slv%nz))

    DO i = 1, n+1
      slv%Ap0(i) = LdA(i) - 1
    END DO
    DO j = 1, slv%nz
      slv%Ai0(j) = ColA(j) - 1
      slv%perm(j) = j
    END DO

    ! insertion sort of each row's columns, carrying the permutation
    DO i = 1, n
      j1 = LdA(i)
      j2 = LdA(i+1) - 1
      DO m = j1+1, j2
        key  = slv%Ai0(m)
        keyp = slv%perm(m)
        k = m - 1
        DO WHILE (k .GE. j1)
          IF (slv%Ai0(k) .LE. key) EXIT
          slv%Ai0(k+1)  = slv%Ai0(k)
          slv%perm(k+1) = slv%perm(k)
          k = k - 1
        END DO
        slv%Ai0(k+1)  = key
        slv%perm(k+1) = keyp
      END DO
    END DO

    CALL umf4def(slv%control)
    slv%initialized = .TRUE.
    ok = .TRUE.
  END SUBROUTINE SD_INIT

  !-----------------------------------------------------------------------
  ! (Re)factorize with the values Avals given in the CALLER's CSR
  ! ordering.  The symbolic factorization is computed once per structure;
  ! the numeric factorization is renewed on every call.
  !-----------------------------------------------------------------------
  SUBROUTINE SD_FACTORIZE(slv, Avals, ok)
    TYPE(tSparseDirectSolver), INTENT(INOUT) :: slv
    REAL*8, INTENT(IN) :: Avals(*)
    LOGICAL, INTENT(OUT) :: ok

    INTEGER :: j

    ok = .FALSE.
    IF (.NOT. slv%initialized) RETURN

    DO j = 1, slv%nz
      slv%Ax(j) = Avals(slv%perm(j))
    END DO

    IF (slv%symbolic .LT. 0) THEN
      CALL umf4sym(slv%n, slv%n, slv%Ap0, slv%Ai0, slv%Ax, slv%symbolic, &
                   slv%control, slv%info)
      IF (slv%symbolic .LT. 0 .OR. slv%info(1) .LT. 0d0) THEN
        WRITE(*,'(A,I0,A,I0,A,F6.0)') 'SD_FACTORIZE: UMFPACK symbolic failed, n = ', &
          slv%n, ', nnz = ', slv%nz, ', info(1) = ', slv%info(1)
        RETURN
      END IF
    END IF

    IF (slv%numeric .GE. 0) THEN
      CALL umf4fnum(slv%numeric)
      slv%numeric = -1
    END IF
    CALL umf4num(slv%Ap0, slv%Ai0, slv%Ax, slv%symbolic, slv%numeric, &
                 slv%control, slv%info)
    IF (slv%numeric .LT. 0 .OR. slv%info(1) .LT. 0d0) THEN
      WRITE(*,'(A,I0,A,I0,A,F6.0)') 'SD_FACTORIZE: UMFPACK numeric failed, n = ', &
        slv%n, ', nnz = ', slv%nz, ', info(1) = ', slv%info(1)
      RETURN
    END IF

    slv%factorized = .TRUE.
    ok = .TRUE.
  END SUBROUTINE SD_FACTORIZE

  !-----------------------------------------------------------------------
  ! Solve A x = b (with iterative refinement).
  !-----------------------------------------------------------------------
  SUBROUTINE SD_SOLVE(slv, x, b, ok)
    TYPE(tSparseDirectSolver), INTENT(INOUT) :: slv
    REAL*8, INTENT(OUT) :: x(*)
    REAL*8, INTENT(IN)  :: b(*)
    LOGICAL, INTENT(OUT) :: ok

    ! sys = 1: solve the transpose system of the factorized CSC matrix,
    ! i.e. A x = b for our CSR input (see module header).
    INTEGER, PARAMETER :: sys_transpose = 1

    ok = .FALSE.
    IF (.NOT. slv%factorized) RETURN

    CALL umf4solr(sys_transpose, slv%Ap0, slv%Ai0, slv%Ax, x, b, &
                  slv%numeric, slv%control, slv%info)
    IF (slv%info(1) .LT. 0d0) RETURN
    ok = .TRUE.
  END SUBROUTINE SD_SOLVE

  !-----------------------------------------------------------------------
  ! Release handles and storage.  Idempotent; safe on partial init.
  !-----------------------------------------------------------------------
  SUBROUTINE SD_FREE(slv)
    TYPE(tSparseDirectSolver), INTENT(INOUT) :: slv

    IF (slv%numeric .GE. 0) THEN
      CALL umf4fnum(slv%numeric)
      slv%numeric = -1
    END IF
    IF (slv%symbolic .GE. 0) THEN
      CALL umf4fsym(slv%symbolic)
      slv%symbolic = -1
    END IF
    IF (ALLOCATED(slv%Ap0))  DEALLOCATE(slv%Ap0)
    IF (ALLOCATED(slv%Ai0))  DEALLOCATE(slv%Ai0)
    IF (ALLOCATED(slv%perm)) DEALLOCATE(slv%perm)
    IF (ALLOCATED(slv%Ax))   DEALLOCATE(slv%Ax)
    slv%initialized = .FALSE.
    slv%factorized = .FALSE.
    slv%n = 0
    slv%nz = 0
  END SUBROUTINE SD_FREE

END MODULE CHI_SPARSE_DIRECT
