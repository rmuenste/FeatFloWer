!=========================================================================
! test_chi_sparse_direct - Phase 1 unit test of CHI_SPARSE_DIRECT
! (chimera-integration-design.md v3, section 9).
!
! The finding-5 regression: two solver INSTANCES with different matrices
! are initialized, factorized, solved and refactorized in interleaved
! order - with the legacy module-global handles this scenario corrupts
! one factorization; with per-instance handles both must stay correct.
! Also checks: unsorted input columns (exercises the value permutation),
! nonsymmetric matrices (exercises the CSR/CSC transpose convention),
! caller arrays untouched, refactorization with new values, idempotent
! SD_FREE.
!
! Serial; prints PASS / STOP 1.
!=========================================================================
PROGRAM test_chi_sparse_direct

  USE CHI_SPARSE_DIRECT, ONLY: tSparseDirectSolver, SD_INIT, SD_FACTORIZE, &
    SD_SOLVE, SD_FREE

  IMPLICIT NONE

  ! System 1: 5x5 nonsymmetric, diagonally dominant, with deliberately
  ! UNSORTED column order inside rows.
  INTEGER, PARAMETER :: n1 = 5
  INTEGER :: Ld1(n1+1), Co1(12), Co1_ref(12)
  REAL*8  :: Av1(12), Av1_ref(12), x1(n1), b1(n1), xe1(n1)

  ! System 2: 7x7 tridiagonal-ish nonsymmetric.
  INTEGER, PARAMETER :: n2 = 7
  INTEGER :: Ld2(n2+1), Co2(19)
  REAL*8  :: Av2(19), x2(n2), b2(n2), xe2(n2)

  TYPE(tSparseDirectSolver) :: s1, s2
  LOGICAL :: ok
  INTEGER :: i, nfail, pos

  nfail = 0

  !--- define system 1 (rows with unsorted columns) ----------------------
  ! row1: cols (3,1,5)   row2: (2,1)   row3: (3)   row4: (5,4,2)  row5: (1,5,4)
  Ld1 = (/ 1, 4, 6, 7, 10, 13 /)
  Co1 = (/ 3, 1, 5,   2, 1,   3,   5, 4, 2,   1, 5, 4 /)
  Av1 = (/ 1.0d0, 10.0d0, 0.5d0, &
           12.0d0, 2.0d0, &
           9.0d0, &
           0.7d0, 11.0d0, 1.5d0, &
           0.3d0, 8.0d0, 1.1d0 /)
  Co1_ref = Co1
  Av1_ref = Av1
  xe1 = (/ 1d0, 2d0, 3d0, 4d0, 5d0 /)
  CALL csr_matvec(n1, Ld1, Co1, Av1, xe1, b1)

  !--- define system 2 ---------------------------------------------------
  Ld2(1) = 1
  DO i = 1, n2
    IF (i .EQ. 1 .OR. i .EQ. n2) THEN
      Ld2(i+1) = Ld2(i) + 2
    ELSE
      Ld2(i+1) = Ld2(i) + 3
    END IF
  END DO
  ! rows: (i-1, i, i+1) with values (-1.0, 6+i*0.1, -2.0)
  pos = 0
  DO i = 1, n2
    IF (i .GT. 1) THEN
      pos = pos + 1
      Co2(pos) = i-1
      Av2(pos) = -1.0d0
    END IF
    pos = pos + 1
    Co2(pos) = i
    Av2(pos) = 6d0 + 0.1d0*DBLE(i)
    IF (i .LT. n2) THEN
      pos = pos + 1
      Co2(pos) = i+1
      Av2(pos) = -2.0d0
    END IF
  END DO
  DO i = 1, n2
    xe2(i) = DBLE(n2 + 1 - i)
  END DO
  CALL csr_matvec(n2, Ld2, Co2, Av2, xe2, b2)

  !--- interleaved lifecycle --------------------------------------------
  CALL SD_INIT(s1, n1, Ld1, Co1, ok)
  IF (.NOT. ok) CALL fail('SD_INIT s1', nfail)
  CALL SD_INIT(s2, n2, Ld2, Co2, ok)
  IF (.NOT. ok) CALL fail('SD_INIT s2', nfail)

  CALL SD_FACTORIZE(s2, Av2, ok)
  IF (.NOT. ok) CALL fail('SD_FACTORIZE s2', nfail)
  CALL SD_FACTORIZE(s1, Av1, ok)
  IF (.NOT. ok) CALL fail('SD_FACTORIZE s1', nfail)

  CALL SD_SOLVE(s1, x1, b1, ok)
  IF (.NOT. ok) CALL fail('SD_SOLVE s1', nfail)
  IF (MAXVAL(ABS(x1 - xe1)) .GT. 1d-10) THEN
    WRITE(*,*) '  s1 error:', MAXVAL(ABS(x1 - xe1))
    CALL fail('s1 solution wrong (interleaved)', nfail)
  END IF

  CALL SD_SOLVE(s2, x2, b2, ok)
  IF (.NOT. ok) CALL fail('SD_SOLVE s2', nfail)
  IF (MAXVAL(ABS(x2 - xe2)) .GT. 1d-10) THEN
    WRITE(*,*) '  s2 error:', MAXVAL(ABS(x2 - xe2))
    CALL fail('s2 solution wrong (interleaved)', nfail)
  END IF

  !--- caller arrays untouched ------------------------------------------
  DO i = 1, 12
    IF (Co1(i) .NE. Co1_ref(i) .OR. ABS(Av1(i)-Av1_ref(i)) .GT. 0d0) THEN
      CALL fail('caller CSR arrays mutated', nfail)
      EXIT
    END IF
  END DO

  !--- refactorize s1 with scaled values; s2 must stay intact ------------
  Av1 = 2d0*Av1
  CALL SD_FACTORIZE(s1, Av1, ok)
  IF (.NOT. ok) CALL fail('SD_FACTORIZE s1 (rescaled)', nfail)
  CALL csr_matvec(n1, Ld1, Co1, Av1, xe1, b1)
  CALL SD_SOLVE(s1, x1, b1, ok)
  IF (.NOT. ok .OR. MAXVAL(ABS(x1 - xe1)) .GT. 1d-10) &
    CALL fail('s1 solve after refactorization', nfail)

  CALL SD_SOLVE(s2, x2, b2, ok)
  IF (.NOT. ok .OR. MAXVAL(ABS(x2 - xe2)) .GT. 1d-10) &
    CALL fail('s2 solve after s1 refactorization', nfail)

  !--- free, idempotent --------------------------------------------------
  CALL SD_FREE(s1)
  CALL SD_FREE(s1)
  CALL SD_FREE(s2)

  IF (nfail .GT. 0) THEN
    WRITE(*,*) 'test_chi_sparse_direct: ', nfail, ' failure(s)'
    STOP 1
  END IF
  WRITE(*,*) 'test_chi_sparse_direct: PASS'

CONTAINS

  SUBROUTINE csr_matvec(n, Ld, Co, Av, x, y)
    INTEGER, INTENT(IN) :: n, Ld(n+1), Co(*)
    REAL*8, INTENT(IN) :: Av(*), x(n)
    REAL*8, INTENT(OUT) :: y(n)
    INTEGER :: ii, jj
    DO ii = 1, n
      y(ii) = 0d0
      DO jj = Ld(ii), Ld(ii+1)-1
        y(ii) = y(ii) + Av(jj)*x(Co(jj))
      END DO
    END DO
  END SUBROUTINE csr_matvec

  SUBROUTINE fail(msg, cnt)
    CHARACTER(*), INTENT(IN) :: msg
    INTEGER, INTENT(INOUT) :: cnt
    WRITE(*,*) 'FAIL: ', msg
    cnt = cnt + 1
  END SUBROUTINE fail

END PROGRAM test_chi_sparse_direct
