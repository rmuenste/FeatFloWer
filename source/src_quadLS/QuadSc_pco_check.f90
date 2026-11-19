!=========================================================================
! QuadSc_pco_check.f90 - diagnostics of the periodic coarse pressure operator
!
! Part of the periodic coarse pressure operator (overview: QuadSc_pco_operator.f90). The
! check of Pres@PeriodicCoarseCheck and the comparison of two operators; nothing here changes
! the operator, the solver copy or the production statistics of the coarse solve.
!
! PCO_COMPARE_OPERATORS(a, ref): entry-by-entry difference, scaled by the largest modulus of
!   the row of ref; a row whose ref entries are all zero is compared against the largest
!   modulus of ref overall (a nonzero value there is a difference, not a skipped row); a
!   non-finite entry of a counts as an infinite difference.
! PCO_RUN_CHECK: given the gathered (ag), the wrapped (aw) and the per-image-mass wrapped (ai)
!   operators on the master (the coordinator builds them, collectively):
!   - entry-by-entry gathered - wrapped (with the location, the values and the legacy master
!     value at the largest difference) and gathered - per-image variant;
!   - for three vectors (mean modes = 1; a four-mode impulse on the element at the largest
!     vertex class; a reproducible vector exciting every mode) the worker product with the
!     E012_DAX arithmetic (C x + CP xp on the workers, distributed and collected with the
!     coarse vector transfers of the multigrid) against the master products, every row,
!     scaled by sum_j |A_ij| max|v|;
!   - closure (max |A e| / (|A| |e|)) and the left measure e^T A on closed components;
!   - a compatible solve with the ungauged residual, the deviation from v up to the kernel and
!     the value of the reference unknowns (PCO_SOLVE_CHECK, diagnostic path);
!   - an incompatible right-hand side (a constant added to the mean modes) must show an
!     imbalance near 1 (PCO_PROBE_IMBALANCE, no solve).
!   Output: PCO_CHECK lines on the terminal of the master. Collective (the worker product).
!=========================================================================
MODULE pco_diagnostics

  USE PP3D_MPI, ONLY: myid, master
  USE def_FEAT, ONLY: MTERM, ILEV
  USE var_QuadScalar, ONLY: mg_mesh, mg_lMat, mg_CMat, mg_lPMat, mg_CPMat
  USE pco_periodic_topology, ONLY: tPcoClasses, PCO_GHOST_BUFFER_SIZE
  USE pco_coarse_solver, ONLY: pco, PCO_SOLVE_CHECK, PCO_PROBE_IMBALANCE, PCO_RAW_MATVEC

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: PCO_COMPARE_OPERATORS, PCO_RUN_CHECK

CONTAINS

  !-----------------------------------------------------------------------
  ! max_k |a(k) - ref(k)| / scale(row of k), see the module header; kmax: its slot (0: none)
  !-----------------------------------------------------------------------
  REAL*8 FUNCTION PCO_COMPARE_OPERATORS(a, ref, kmax) RESULT(d)
    REAL*8, INTENT(IN) :: a(:), ref(:)
    INTEGER, INTENT(OUT) :: kmax
    INTEGER :: row, k
    REAL*8 :: rmax, gmax, x
    d = 0d0
    kmax = 0
    gmax = MAXVAL(ABS(ref))
    DO row = 1, pco%n
      rmax = MAXVAL(ABS(ref(pco%ldA(row):pco%ldA(row+1)-1)))
      IF (rmax.LE.0d0) rmax = gmax
      IF (rmax.LE.0d0) rmax = 1d0
      DO k = pco%ldA(row), pco%ldA(row+1)-1
        IF (a(k).NE.a(k) .OR. ABS(a(k)).GT.HUGE(1d0)) THEN
          x = HUGE(1d0)
        ELSE
          x = ABS(a(k) - ref(k))/rmax
        END IF
        IF (x.GT.d) THEN
          d = x
          kmax = k
        END IF
      END DO
    END DO
  END FUNCTION PCO_COMPARE_OPERATORS

  !-----------------------------------------------------------------------
  ! The check (module header). ag/aw/ai: master operators in the pco layout (dummies on the
  ! workers); weighted: ai is the per-image-mass variant (else equal to aw); crs: coarse
  ! solver type; nrenew: renewal number for the output.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_RUN_CHECK(lev, crs, nrenew, cls, ag, aw, ai, weighted)
    INTEGER, INTENT(IN) :: lev, crs, nrenew
    TYPE(tPcoClasses), INTENT(IN) :: cls
    REAL*8, INTENT(IN) :: ag(:), aw(:), ai(:)
    LOGICAL, INTENT(IN) :: weighted
    REAL*8, ALLOCATABLE :: v(:), yd(:), yg(:), yw(:), yi(:), sc(:), b(:), cs(:), ca(:)
    INTEGER :: n, nnz, k, row, iv, e, kmax, idum
    REAL*8 :: dgw, dgi, eg, ew, ei, resid, devi, cnorm, tol_res, tol_dev, cleft, xref, imb
    CHARACTER(len=8) :: vname(3) = (/'mean    ', 'impulse ', 'random  '/)

    IF (myid.EQ.master) THEN
      n = pco%n
      nnz = SIZE(pco%colA)
      ALLOCATE(v(n), yd(n), yg(n), yw(n), yi(n), sc(n))
      dgw = PCO_COMPARE_OPERATORS(aw, ag, kmax)
      dgi = PCO_COMPARE_OPERATORS(ai, ag, idum)
      WRITE(MTERM,'(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A,L1,A,A)') 'PCO_CHECK renewal ', nrenew, ' level ', lev, &
        ' entries ', nnz, ' |gathered-wrapped|/rowmax ', dgw, ' per-image-mass variant ', dgi, &
        ' (weighted ', weighted, ') ', MERGE('PASS', 'FAIL', dgw.LE.1d-12)
      IF (kmax.GT.0) THEN
        row = row_of(kmax)
        WRITE(MTERM,'(A,I0,A,2(I0,A,I0,A,3ES11.3,A),3(A,ES12.4))') 'PCO_CHECK renewal ', nrenew, &
          ' largest difference at ', (row-1)/4+1, ' mode ', MOD(row-1,4)+1, ' (', centroid((row-1)/4+1), ')', &
          (pco%colA(kmax)-1)/4+1, ' mode ', MOD(pco%colA(kmax)-1,4)+1, ' (', centroid((pco%colA(kmax)-1)/4+1), &
          ')', ' gathered ', ag(kmax), ' wrapped ', aw(kmax), ' legacy master C ', legacy_entry(row, pco%colA(kmax))
      END IF
    ELSE
      n = 4*mg_mesh%level(lev)%nel
      ALLOCATE(v(n), yd(n), yg(1), yw(1), yi(1), sc(1))
    END IF

    DO iv = 1, 3
      IF (myid.EQ.master) CALL test_vector(iv, v)
      CALL worker_product(lev, v, yd)
      IF (myid.EQ.master) THEN
        CALL matvec(ag, v, yg)
        CALL matvec(aw, v, yw)
        CALL matvec(ai, v, yi)
        DO row = 1, n
          sc(row) = SUM(ABS(ag(pco%ldA(row):pco%ldA(row+1)-1)))*MAXVAL(ABS(v))
        END DO
        eg = 0d0
        ew = 0d0
        ei = 0d0
        DO row = 1, n
          IF (sc(row).LE.0d0) CYCLE
          eg = MAX(eg, ABS(yd(row)-yg(row))/sc(row))
          ew = MAX(ew, ABS(yd(row)-yw(row))/sc(row))
          ei = MAX(ei, ABS(yd(row)-yi(row))/sc(row))
        END DO
        cnorm = MAXVAL(ABS(yd)/MAX(sc, 1d-300))
        WRITE(MTERM,'(A,I0,A,A,A,I0,3(A,ES10.3),A,ES10.3,A,A)') 'PCO_CHECK renewal ', nrenew, ' vector ', &
          TRIM(vname(iv)), ' rows ', n, ' dax-gathered ', eg, ' dax-wrapped ', ew, ' dax-per-image ', ei, &
          ' max |Av|/(|A| |v|) ', cnorm, ' ', MERGE('PASS', 'FAIL', eg.LE.1d-12 .AND. ew.LE.1d-12)
      END IF
    END DO

    IF (myid.EQ.master) THEN
      ! direct solvers: rounding level; BiCGStab (types 1, 3): its tolerance 1e-11 in
      ! PCO_SOLVE_CHECK, the solution error amplified by the conditioning
      tol_res = MERGE(1d-10, 1d-8, crs.EQ.2 .OR. crs.EQ.4)
      tol_dev = MERGE(1d-8, 1d-3, crs.EQ.2 .OR. crs.EQ.4)
      ALLOCATE(cs(n), ca(n))
      cs = 0d0
      ca = 0d0
      DO row = 1, n, 4
        DO k = pco%ldA(row), pco%ldA(row+1)-1
          cs(pco%colA(k)) = cs(pco%colA(k)) + ag(k)
          ca(pco%colA(k)) = ca(pco%colA(k)) + ABS(ag(k))
        END DO
      END DO
      cleft = 0d0
      DO k = 1, n
        IF (pco%comp_closed(pco%comp_of((k-1)/4+1)) .AND. ca(k).GT.0d0) cleft = MAX(cleft, ABS(cs(k))/ca(k))
      END DO
      DEALLOCATE(cs, ca)
      CALL test_vector(3, v)
      CALL PCO_SOLVE_CHECK(v, resid, devi, xref)
      WRITE(MTERM,'(A,I0,A,I0,4(A,ES10.3),A,A)') 'PCO_CHECK renewal ', nrenew, ' solver type ', &
        crs, ' left |e^T A|/(e^T |A|) ', cleft, ' compatible rhs: ungauged residual ', resid, &
        ' solution - v off the kernel ', devi, ' reference value ', xref, ' ', &
        MERGE('PASS', 'FAIL', resid.LE.tol_res .AND. devi.LE.tol_dev .AND. cleft.LE.1d-12 .AND. xref.LE.tol_res)
      IF (pco%nclosed.GT.0) THEN
        ! raw v + a constant on the mean modes: not in the range of a closed operator
        ALLOCATE(b(n))
        CALL PCO_RAW_MATVEC(v, b)
        imb = MAXVAL(ABS(b))
        DO e = 1, pco%nel
          b(4*(e-1)+1) = b(4*(e-1)+1) + imb
        END DO
        imb = PCO_PROBE_IMBALANCE(b)
        WRITE(MTERM,'(A,I0,A,ES10.3,A,A)') 'PCO_CHECK renewal ', nrenew, ' incompatible rhs: imbalance ', &
          imb, ' ', MERGE('PASS', 'FAIL', imb.GT.0.5d0)
        DEALLOCATE(b)
      END IF
    END IF
    DEALLOCATE(v, yd, yg, yw, yi, sc)

  CONTAINS

    INTEGER FUNCTION row_of(slot)
      INTEGER, INTENT(IN) :: slot
      INTEGER :: l, r, mid
      l = 1
      r = pco%n
      DO WHILE (l.LT.r)
        mid = (l + r + 1)/2
        IF (pco%ldA(mid).LE.slot) THEN
          l = mid
        ELSE
          r = mid - 1
        END IF
      END DO
      row_of = l
    END FUNCTION row_of

    ! value of the own (legacy, non-periodic) matrix of the master at (r, c), 0 if not in its pattern
    REAL*8 FUNCTION legacy_entry(r, c) RESULT(val)
      INTEGER, INTENT(IN) :: r, c
      INTEGER :: kk
      val = 0d0
      DO kk = mg_lMat(lev)%LdA(r), mg_lMat(lev)%LdA(r+1)-1
        IF (mg_lMat(lev)%ColA(kk).EQ.c) val = mg_CMat(lev)%a(kk)
      END DO
    END FUNCTION legacy_entry

    FUNCTION centroid(e0) RESULT(c)
      INTEGER, INTENT(IN) :: e0
      REAL*8 :: c(3)
      INTEGER :: j
      c = 0d0
      DO j = 1, 8
        c = c + 0.125d0*mg_mesh%level(lev)%dcorvg(:, mg_mesh%level(lev)%kvert(j,e0))
      END DO
    END FUNCTION centroid

    SUBROUTINE matvec(a, x, y)
      REAL*8, INTENT(IN) :: a(:), x(:)
      REAL*8, INTENT(OUT) :: y(:)
      INTEGER :: i, kk
      DO i = 1, pco%n
        y(i) = 0d0
        DO kk = pco%ldA(i), pco%ldA(i+1)-1
          y(i) = y(i) + a(kk)*x(pco%colA(kk))
        END DO
      END DO
    END SUBROUTINE matvec

    ! the three test vectors in the master numbering
    SUBROUTINE test_vector(iv0, w)
      INTEGER, INTENT(IN) :: iv0
      REAL*8, INTENT(OUT) :: w(:)
      INTEGER :: ee, kk, best, msz, i
      INTEGER(8) :: sd
      w = 0d0
      SELECT CASE (iv0)
      CASE (1)
        DO ee = 1, pco%nel
          w(4*(ee-1)+1) = 1d0
        END DO
      CASE (2)
        best = 1
        msz = 0
        DO kk = 1, cls%ncls
          IF (cls%mem(cls%ptr(kk)).GT.cls%nvt) CYCLE
          IF (cls%ptr(kk+1)-cls%ptr(kk).GT.msz) THEN
            msz = cls%ptr(kk+1)-cls%ptr(kk)
            best = cls%mem(cls%ptr(kk))
          END IF
        END DO
        ee = 1
        DO i = 1, mg_mesh%level(lev)%nel
          IF (ANY(mg_mesh%level(lev)%kvert(:,i).EQ.best)) THEN
            ee = i
            EXIT
          END IF
        END DO
        w(4*(ee-1)+1:4*ee) = (/1d0, 0.5d0, -0.25d0, 0.75d0/)
      CASE DEFAULT
        sd = 20261119_8
        DO i = 1, SIZE(w)
          sd = MOD(1103515245_8*sd + 12345_8, 2147483648_8)
          w(i) = 2d0*DBLE(sd)/2147483648d0 - 1d0
        END DO
      END SELECT
    END SUBROUTINE test_vector

  END SUBROUTINE PCO_RUN_CHECK

  ! y = worker operator applied to v (master numbering on the master; collective):
  ! distribute v to the workers, E012_DAX arithmetic there, collect the result
  SUBROUTINE worker_product(lev, v, y)
    INTEGER, INTENT(IN) :: lev
    REAL*8, INTENT(INOUT) :: v(:), y(:)
    INTEGER :: nel, i, k, ilev_save
    REAL*8, ALLOCATABLE :: xp(:)

    ilev_save = ILEV
    ILEV = lev
    nel = mg_mesh%level(lev)%nel
    CALL coarse_to_workers(lev, v, nel)
    IF (myid.NE.master) THEN
      ALLOCATE(xp(PCO_GHOST_BUFFER_SIZE(lev)))
      xp = 0d0
    ELSE
      ALLOCATE(xp(4))
    END IF
    CALL GetParPressure(v, xp)
    IF (myid.NE.master) THEN
      DO i = 1, 4*nel
        y(i) = 0d0
        DO k = mg_lMat(lev)%LdA(i), mg_lMat(lev)%LdA(i+1)-1
          y(i) = y(i) + mg_CMat(lev)%a(k)*v(mg_lMat(lev)%ColA(k))
        END DO
        DO k = mg_lPMat(lev)%LdA(i), mg_lPMat(lev)%LdA(i+1)-1
          y(i) = y(i) + mg_CPMat(lev)%a(k)*xp(mg_lPMat(lev)%ColA(k))
        END DO
      END DO
    END IF
    DEALLOCATE(xp)
    CALL workers_to_coarse(lev, y, nel)
    ILEV = ilev_save
  END SUBROUTINE worker_product

  ! the coarse vector transfers of the multigrid (QuadSc_mpi.f90), by level
  SUBROUTINE coarse_to_workers(lev, d, nel)
    INTEGER, INTENT(IN) :: lev, nel
    REAL*8, INTENT(INOUT) :: d(:)
    EXTERNAL E012GATHR_L1, E012GATHR_L2, E012GATHR_L3
    IF (lev.EQ.1) CALL E012GATHR_L1(d, nel)
    IF (lev.EQ.2) CALL E012GATHR_L2(d, nel)
    IF (lev.EQ.3) CALL E012GATHR_L3(d, nel)
  END SUBROUTINE coarse_to_workers

  SUBROUTINE workers_to_coarse(lev, d, nel)
    INTEGER, INTENT(IN) :: lev, nel
    REAL*8, INTENT(INOUT) :: d(:)
    EXTERNAL E012DISTR_L1, E012DISTR_L2, E012DISTR_L3
    IF (lev.EQ.1) CALL E012DISTR_L1(d, nel)
    IF (lev.EQ.2) CALL E012DISTR_L2(d, nel)
    IF (lev.EQ.3) CALL E012DISTR_L3(d, nel)
  END SUBROUTINE workers_to_coarse

END MODULE pco_diagnostics
