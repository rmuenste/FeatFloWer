!=========================================================================
! test_chi_algebra - Phase 4 penalty-algebra test (chimera-integration-
! design.md v3, section 9; paper eqs. (7), (10), (12)).
!
! A uniform box mesh [0,1]^3 (8x8x8, loaded through the legacy mesh
! adapter path so the Q2 numbering is the production one) carries one
! sphere body with an atmosphere.  The interior-penalty operator D and
! vector g are tabulated/assembled with CHI_PENALTY on a Q2 CSR pattern
! built here (FEAT layout, diagonal first).  Checks:
!   1. D is symmetric;
!   2. consistency: for uhat = const the quadrature identity D*(c 1) = g
!      holds to round-off (sum of the Q2 basis = 1);
!   3. sum_ij d_ij = gamma * int beta dx against the analytic volume of
!      the damping function (kinked integrand: quadrature tolerance 3 %);
!   4. D is positive semi-definite on random vectors;
!   5. the correction CG (paper eq. (12)) reduces the residual of
!      [ML + coef D] x = rhs below 1e-9 with an active Dirichlet filter
!      (filtered rows stay exactly zero);
!   6. with D = 0 the CG returns rhs/ML bitwise in one iteration (the
!      algebraic counterpart of the API fast path).
!=========================================================================
PROGRAM test_chi_algebra
  USE types, ONLY: tMultiMesh
  USE mesh_structures, ONLY: readTriCoarse, genMeshStructures, &
    getNumberOfEdgesOnVerts
  USE CHI_MARKERS, ONLY: tChiBody, CHI_BODY_SPHERE
  USE CHI_FEM_EVAL, ONLY: CHI_Q2_DOFMAP
  USE CHI_PENALTY, ONLY: tChiPenaltyTab, CHI_PENALTY_NQ, CHI_PENALTY_TABULATE, &
    CHI_PENALTY_RELEASE, CHI_PENALTY_ASSEMBLE_D, CHI_PENALTY_ASSEMBLE_G, &
    CHI_PENALTY_DIAG, CHI_PENALTY_MATVEC, CHI_PCG3, CHI_PENALTY_LUMP, CHI_PENALTY_NODAL
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: NX = 8
  REAL*8, PARAMETER :: PI = 3.14159265358979323846d0
  TYPE(tMultiMesh) :: m
  TYPE(tChiBody) :: body(1)
  TYPE(tChiPenaltyTab) :: tab
  INTEGER :: ierr, nfail, ndof, nel, nvt, net, nat, na, i, j, p, q, a, it, itrial
  INTEGER, ALLOCATABLE :: LdA(:), ColA(:)
  REAL*8, ALLOCATABLE :: d(:), g1(:), g2(:), g3(:), t1(:), t2(:), t3(:), one(:)
  REAL*8, ALLOCATABLE :: ml(:), dglob(:), wts(:), diag(:), rhs(:,:), x(:,:), r(:,:)
  REAL*8, ALLOCATABLE :: xr(:), dl(:), xnd(:,:)
  INTEGER, ALLOCATABLE :: nbo(:)
  LOGICAL, ALLOCATABLE :: inb(:)
  REAL*8 :: gamma, hw(1), dmax, err, vol, va, vb, aa, bb, s, coef, resid, rn, xn
  LOGICAL, ALLOCATABLE :: masked(:)

  nfail = 0
  CALL MPI_INIT(ierr)

  CALL write_box('tca_box.tri', NX)
  CALL load('tca_box.tri', m)
  nel = m%level(1)%nel; nvt = m%level(1)%nvt
  net = m%level(1)%net; nat = m%level(1)%nat
  ndof = nvt + net + nat + nel
  CALL build_q2_csr(m, ndof, LdA, ColA)
  na = LdA(ndof+1) - 1
  WRITE(*,'(A,I0,A,I0,A,I0)') 'test_chi_algebra: nel = ', nel, ', Q2 dofs = ', ndof, &
    ', nnz = ', na

  body(1)%shape = CHI_BODY_SPHERE
  body(1)%center = (/ 0.5d0, 0.5d0, 0.5d0 /)
  body(1)%radius = 0.1d0
  hw(1) = 0.6d0
  gamma = 10d0

  CALL CHI_PENALTY_TABULATE(nel, m%level(1)%kvert, m%level(1)%dcorvg, 1, body, hw, tab, &
    0.5d0, 0.75d0)
  WRITE(*,'(A,I0)') 'test_chi_algebra: active elements = ', tab%nact
  IF (tab%nact .LE. 0) CALL fail('no active elements', nfail)

  ALLOCATE(d(na), g1(ndof), g2(ndof), g3(ndof), t1(ndof), t2(ndof), t3(ndof), one(ndof))
  CALL CHI_PENALTY_ASSEMBLE_D(tab, gamma, m%level(1)%kvert, m%level(1)%kedge, &
    m%level(1)%karea, nvt, net, nat, LdA, ColA, ndof, d)
  dmax = MAXVAL(ABS(d(1:na)))

  ! ---- 1. symmetry ------------------------------------------------------
  err = 0d0
  DO i = 1, ndof
    DO p = LdA(i), LdA(i+1)-1
      j = ColA(p)
      DO q = LdA(j), LdA(j+1)-1
        IF (ColA(q) .EQ. i) THEN
          err = MAX(err, ABS(d(p) - d(q)))
          EXIT
        END IF
      END DO
    END DO
  END DO
  WRITE(*,'(A,ES10.2)') 'test_chi_algebra: max |d_ij - d_ji| / max|d| = ', err/dmax
  IF (err .GT. 1d-13*dmax) CALL fail('D not symmetric', nfail)

  ! ---- 2. consistency D*(c 1) = g(uhat = c) ------------------------------
  DO a = 1, tab%nact
    DO q = 1, CHI_PENALTY_NQ
      tab%uhat(:,q,a) = (/ 1d0, 2d0, 3d0 /)
    END DO
  END DO
  CALL CHI_PENALTY_ASSEMBLE_G(tab, gamma, m%level(1)%kvert, m%level(1)%kedge, &
    m%level(1)%karea, nvt, net, nat, ndof, g1, g2, g3)
  one = 1d0
  t1 = 0d0; t2 = 0d0; t3 = 0d0
  CALL CHI_PENALTY_MATVEC(LdA, ColA, ndof, d, 1d0, one, t1)
  CALL CHI_PENALTY_MATVEC(LdA, ColA, ndof, d, 2d0, one, t2)
  CALL CHI_PENALTY_MATVEC(LdA, ColA, ndof, d, 3d0, one, t3)
  err = MAX(MAXVAL(ABS(t1-g1)), MAXVAL(ABS(t2-g2)), MAXVAL(ABS(t3-g3)))
  s = MAX(MAXVAL(ABS(g1)), MAXVAL(ABS(g2)), MAXVAL(ABS(g3)))
  WRITE(*,'(A,ES10.2)') 'test_chi_algebra: max |D c1 - g(c)| / max|g| = ', err/s
  IF (err .GT. 1d-12*s) CALL fail('D*1 differs from g(uhat = const)', nfail)

  ! ---- 3. penalised volume ---------------------------------------------
  vol = SUM(t1)/gamma
  aa = body(1)%radius + 0.5d0*hw(1)
  bb = body(1)%radius + 0.75d0*hw(1)
  va = 4d0/3d0*PI*aa**3
  vb = 4d0*PI/(0.25d0*hw(1)) * (bb*(bb**3 - aa**3)/3d0 - (bb**4 - aa**4)/4d0)
  err = ABS(vol - (va+vb))/(va+vb)
  WRITE(*,'(A,ES12.5,A,ES12.5,A,ES10.2)') 'test_chi_algebra: int beta = ', vol, &
    ' analytic = ', va+vb, ' rel. err = ', err
  IF (err .GT. 3d-2) CALL fail('penalised volume off', nfail)

  ! ---- 3b. lumping: row sums equal D*1 (diagnostic, may be negative);
  !          the nodal Lobatto lumping is positive and integrates beta ----
  ALLOCATE(dl(ndof), nbo(ndof), inb(ndof), xnd(3,ndof))
  CALL CHI_PENALTY_LUMP(LdA, ndof, d, dl)
  IF (MAXVAL(ABS(dl - t1)) .GT. 1d-12*MAXVAL(ABS(t1))) &
    CALL fail('row-sum lumping differs from D*1', nfail)
  WRITE(*,'(A,ES10.2)') 'test_chi_algebra: row-sum lumping min weight = ', MINVAL(dl)
  CALL CHI_PENALTY_NODAL(nel, nvt, net, nat, m%level(1)%kvert, m%level(1)%kedge, &
    m%level(1)%karea, m%level(1)%dcorvg, 1, body, hw, gamma, ndof, dl, nbo, inb, xnd, &
    0.5d0, 0.75d0)
  IF (MINVAL(dl) .LT. 0d0) CALL fail('nodal lumped penalty has a negative weight', nfail)
  err = ABS(SUM(dl)/gamma - (va+vb))/(va+vb)
  WRITE(*,'(A,ES12.5,A,ES10.2)') 'test_chi_algebra: nodal int beta = ', SUM(dl)/gamma, &
    ' rel. err = ', err
  IF (err .GT. 5d-2) CALL fail('nodal lumped volume off', nfail)
  DO i = 1, ndof
    IF ((dl(i) .GT. 0d0) .NEQV. (nbo(i) .GT. 0)) THEN
      CALL fail('nodal body index inconsistent with the weight', nfail)
      EXIT
    END IF
  END DO

  ! ---- 4. positive semi-definiteness -------------------------------------
  ALLOCATE(xr(ndof))
  DO itrial = 1, 5
    DO i = 1, ndof
      xr(i) = pseudo_random(itrial*7919 + i) - 0.5d0
    END DO
    t1 = 0d0
    CALL CHI_PENALTY_MATVEC(LdA, ColA, ndof, d, 1d0, xr, t1)
    s = SUM(xr*t1)
    IF (s .LT. -1d-12*dmax*SUM(xr*xr)) CALL fail('x^T D x < 0', nfail)
  END DO

  ! ---- 5. correction CG with a Dirichlet filter -------------------------
  ALLOCATE(ml(ndof), dglob(ndof), wts(ndof), diag(ndof), rhs(ndof,3), x(ndof,3), &
           r(ndof,3), masked(ndof))
  DO i = 1, ndof
    ml(i) = 1d0 + DBLE(MOD(i,7))/7d0
    wts(i) = 1d0
  END DO
  CALL CHI_PENALTY_DIAG(LdA, ndof, d, diag)
  coef = 0.05d0
  dglob = ml + coef*diag
  ! mask: every 11th dof plus the first 20 (a fake Dirichlet set)
  DO i = 1, ndof
    masked(i) = (MOD(i,11) .EQ. 0) .OR. (i .LE. 20)
    rhs(i,1) = pseudo_random(31*i) - 0.5d0
    rhs(i,2) = pseudo_random(37*i+1) - 0.5d0
    rhs(i,3) = pseudo_random(41*i+2) - 0.5d0
  END DO
  CALL mask_filter(rhs(:,1), rhs(:,2), rhs(:,3), ndof)
  CALL CHI_PCG3(ndof, LdA, ColA, d, coef, ml, dglob, wts, rhs(:,1), rhs(:,2), rhs(:,3), &
    x(:,1), x(:,2), x(:,3), 1d-12, 500, nosum3, mask_filter, noallsum, it, resid)
  WRITE(*,'(A,I0,A,ES10.2)') 'test_chi_algebra: CG iterations = ', it, &
    ', reported rel. residual = ', resid
  ! explicit residual
  err = 0d0
  DO j = 1, 3
    r(:,j) = 0d0
    CALL CHI_PENALTY_MATVEC(LdA, ColA, ndof, d, coef, x(:,j), r(:,j))
    r(:,j) = r(:,j) + ml*x(:,j)
  END DO
  CALL mask_filter(r(:,1), r(:,2), r(:,3), ndof)
  r = r - rhs
  rn = SQRT(SUM(r*r)); xn = SQRT(SUM(rhs*rhs))
  WRITE(*,'(A,ES10.2)') 'test_chi_algebra: explicit |A x - rhs| / |rhs| = ', rn/xn
  IF (rn .GT. 1d-9*xn) CALL fail('CG residual too large', nfail)
  IF (it .GE. 500) CALL fail('CG did not converge', nfail)
  DO i = 1, ndof
    IF (masked(i)) THEN
      IF (x(i,1) .NE. 0d0 .OR. x(i,2) .NE. 0d0 .OR. x(i,3) .NE. 0d0) THEN
        CALL fail('filtered dof moved', nfail)
        EXIT
      END IF
    END IF
  END DO

  ! ---- 6. D = 0: x = rhs/ml bitwise ---------------------------------------
  d = 0d0
  dglob = ml
  DO i = 1, ndof
    rhs(i,1) = pseudo_random(53*i) - 0.5d0
    rhs(i,2) = pseudo_random(59*i+1) - 0.5d0
    rhs(i,3) = pseudo_random(61*i+2) - 0.5d0
  END DO
  CALL CHI_PCG3(ndof, LdA, ColA, d, coef, ml, dglob, wts, rhs(:,1), rhs(:,2), rhs(:,3), &
    x(:,1), x(:,2), x(:,3), 1d-12, 500, nosum3, identity_filter, noallsum, it, resid)
  WRITE(*,'(A,I0)') 'test_chi_algebra: D = 0 CG iterations = ', it
  DO j = 1, 3
    DO i = 1, ndof
      IF (x(i,j) .NE. rhs(i,j)/ml(i)) THEN
        CALL fail('D = 0 correction is not rhs/ml bitwise', nfail)
        EXIT
      END IF
    END DO
  END DO

  CALL CHI_PENALTY_RELEASE(tab)
  CALL MPI_FINALIZE(ierr)
  IF (nfail .EQ. 0) THEN
    WRITE(*,'(A)') 'test_chi_algebra: PASSED'
  ELSE
    WRITE(*,'(A,I0,A)') 'test_chi_algebra: FAILED (', nfail, ' check(s))'
    STOP 1
  END IF

CONTAINS

  SUBROUTINE fail(msg, n)
    CHARACTER(*), INTENT(IN) :: msg
    INTEGER, INTENT(INOUT) :: n
    n = n + 1
    WRITE(*,'(A)') 'test_chi_algebra: FAIL - ' // msg
  END SUBROUTINE fail

  REAL*8 FUNCTION pseudo_random(k)
    INTEGER, INTENT(IN) :: k
    INTEGER*8 :: z
    z = INT(k, 8)*2654435761_8 + 12345_8
    z = MOD(z, 4294967296_8)
    IF (z .LT. 0) z = z + 4294967296_8
    pseudo_random = DBLE(z)/4294967296d0
  END FUNCTION pseudo_random

  SUBROUTINE nosum3(y1, y2, y3)
    REAL*8, INTENT(INOUT) :: y1(*), y2(*), y3(*)
  END SUBROUTINE nosum3

  SUBROUTINE noallsum(vals, n)
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(INOUT) :: vals(n)
  END SUBROUTINE noallsum

  SUBROUTINE identity_filter(y1, y2, y3, n)
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(INOUT) :: y1(n), y2(n), y3(n)
  END SUBROUTINE identity_filter

  SUBROUTINE mask_filter(y1, y2, y3, n)
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(INOUT) :: y1(n), y2(n), y3(n)
    INTEGER :: ii
    DO ii = 1, n
      IF ((MOD(ii,11) .EQ. 0) .OR. (ii .LE. 20)) THEN
        y1(ii) = 0d0; y2(ii) = 0d0; y3(ii) = 0d0
      END IF
    END DO
  END SUBROUTINE mask_filter

  SUBROUTINE write_box(fname, n)
    CHARACTER(*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: n
    INTEGER :: iu, ix, iy, iz, i
    iu = 91
    OPEN(iu, FILE=fname, STATUS='REPLACE')
    WRITE(iu,'(A)') 'test_chi_algebra box'
    WRITE(iu,'(A)') 'Parametrisierung PARXC, PARYC, TMAXC'
    WRITE(iu,'(I0,1X,I0,A)') n*n*n, (n+1)**3, ' 1 8 12 6     NEL NVT NBCT NVE NEE NAE'
    WRITE(iu,'(A)') 'DCORVG'
    DO iz = 0, n
      DO iy = 0, n
        DO ix = 0, n
          WRITE(iu,'(3ES24.16)') DBLE(ix)/n, DBLE(iy)/n, DBLE(iz)/n
        END DO
      END DO
    END DO
    WRITE(iu,'(A)') 'KVERT'
    DO iz = 0, n-1
      DO iy = 0, n-1
        DO ix = 0, n-1
          WRITE(iu,'(8(I0,1X))') vid(n,ix,iy,iz), vid(n,ix+1,iy,iz), vid(n,ix+1,iy+1,iz), &
            vid(n,ix,iy+1,iz), vid(n,ix,iy,iz+1), vid(n,ix+1,iy,iz+1), &
            vid(n,ix+1,iy+1,iz+1), vid(n,ix,iy+1,iz+1)
        END DO
      END DO
    END DO
    WRITE(iu,'(A)') 'KNPR'
    DO i = 1, (n+1)**3
      WRITE(iu,'(I0)') 1
    END DO
    CLOSE(iu)
  END SUBROUTINE write_box

  INTEGER FUNCTION vid(n, ix, iy, iz)
    INTEGER, INTENT(IN) :: n, ix, iy, iz
    vid = 1 + ix + (n+1)*(iy + (n+1)*iz)
  END FUNCTION vid

  SUBROUTINE load(fname, mm)
    CHARACTER(*), INTENT(IN) :: fname
    TYPE(tMultiMesh), INTENT(INOUT) :: mm
    INTEGER :: noe
    ALLOCATE(mm%level(1))
    mm%nlmin = 1; mm%nlmax = 1; mm%maxlevel = 1
    mm%level(1)%nel = 0
    CALL readTriCoarse(fname, mm)
    CALL getNumberOfEdgesOnVerts(mm%level(1), noe)
    CALL genMeshStructures(mm, .FALSE., 1, noe)
  END SUBROUTINE load

  ! Q2 CSR pattern in FEAT layout (diagonal first, then ascending columns).
  SUBROUTINE build_q2_csr(mm, n, kld, kcol)
    TYPE(tMultiMesh), INTENT(IN) :: mm
    INTEGER, INTENT(IN) :: n
    INTEGER, ALLOCATABLE, INTENT(OUT) :: kld(:), kcol(:)
    INTEGER, PARAMETER :: MAXN = 125
    INTEGER, ALLOCATABLE :: nbr(:,:), cnt(:)
    INTEGER :: e, idx(27), ii, jj, k, pos, tmp, l
    LOGICAL :: have
    ALLOCATE(nbr(MAXN, n), cnt(n))
    cnt = 0
    DO e = 1, mm%level(1)%nel
      CALL CHI_Q2_DOFMAP(e, mm%level(1)%kvert, mm%level(1)%kedge, mm%level(1)%karea, &
        mm%level(1)%nvt, mm%level(1)%net, mm%level(1)%nat, idx)
      DO ii = 1, 27
        DO jj = 1, 27
          have = .FALSE.
          DO k = 1, cnt(idx(ii))
            IF (nbr(k, idx(ii)) .EQ. idx(jj)) THEN
              have = .TRUE.
              EXIT
            END IF
          END DO
          IF (.NOT. have) THEN
            cnt(idx(ii)) = cnt(idx(ii)) + 1
            nbr(cnt(idx(ii)), idx(ii)) = idx(jj)
          END IF
        END DO
      END DO
    END DO
    ALLOCATE(kld(n+1))
    kld(1) = 1
    DO ii = 1, n
      kld(ii+1) = kld(ii) + cnt(ii)
    END DO
    ALLOCATE(kcol(kld(n+1)-1))
    DO ii = 1, n
      ! insertion sort of the neighbours, then rotate the diagonal to the front
      DO k = 2, cnt(ii)
        tmp = nbr(k, ii)
        l = k - 1
        DO WHILE (l .GE. 1)
          IF (nbr(l, ii) .LE. tmp) EXIT
          nbr(l+1, ii) = nbr(l, ii)
          l = l - 1
        END DO
        nbr(l+1, ii) = tmp
      END DO
      pos = kld(ii)
      kcol(pos) = ii
      pos = pos + 1
      DO k = 1, cnt(ii)
        IF (nbr(k, ii) .EQ. ii) CYCLE
        kcol(pos) = nbr(k, ii)
        pos = pos + 1
      END DO
    END DO
    DEALLOCATE(nbr, cnt)
  END SUBROUTINE build_q2_csr

END PROGRAM test_chi_algebra
