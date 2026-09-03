!=========================================================================
! CHI_PENALTY - distributed interior-penalty operator of the weak
! Chimera variant (Chimera-W; paper arXiv:2506.22831 eqs. (7)-(12);
! chimera-integration-design.md v3, sections 3 and 5).
!
! Layer M (ff_chimera): reentrant, COMMON-free, no solver state.  The
! operator is
!
!   d_ij = gamma_max * sum_k [ int_{Omega_k}  beta_k phi_j phi_i dx
!                            + int_{B_k}      phi_j phi_i dx ],
!   g_i  = gamma_max * sum_k [ int_{Omega_k}  beta_k uhat . phi_i dx
!                            + int_{B_k}      U    . phi_i dx ],
!
! with the damping function of the paper,
!
!   beta_k(x) = min(1, max(0, (R_k + 0.75 H_k - |x - X_k|) / (0.25 H_k))),
!
! i.e. beta = 1 inside the body and up to R + H/2, a linear ramp to 0 at
! R + 3H/4, and 0 in the outer quarter of the atmosphere (so that the
! penalty never interferes with the Robin condition on Gamma_k).
!
! Both integrals are evaluated with the 3x3x3 Gauss rule on every
! background element whose points carry beta > 0 ("active" elements).
! The point data are tabulated once (static geometry): quadrature
! weight x |J| x beta, the body index, an inside-body flag and the
! physical point.  The Dirichlet datum uhat at the points is supplied by
! the caller (Layer H evaluates the replicated submesh solutions there)
! before every g assembly; inside a body it is the rigid velocity.
!
! The matrix shares the Q2 CSR pattern of the background level (all
! contributions are element-local, hence inside the Q2 stencil); FEAT
! convention: the first entry of every row is the diagonal.
!
! CHI_PCG3 solves the correction system of paper eq. (12),
!   [M_L + coef D] delta = rhs        (coef = dt in FeatFloWer scaling)
! for three components at once with a Jacobi-preconditioned CG whose
! parallel contract is the design's: the caller passes the assembly sum
! (E013Sum3) for partial vectors, the Dirichlet defect filter, the
! communicator-wide scalar sum, and partition-of-unity weights for the
! scalar products.
!=========================================================================
MODULE CHI_PENALTY
  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP, CHI_GAUSS3
  USE CHI_FEM_EVAL, ONLY: CHI_Q2_BASIS, CHI_Q2_DOFMAP, CHI_Q2_REFNODES
  USE CHI_MARKERS, ONLY: tChiBody, CHI_BODY_CYLINDER_Z, CHI_BODY_SPHERE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tChiPenaltyTab
  PUBLIC :: CHI_PENALTY_NQ
  PUBLIC :: CHI_PENALTY_BETA
  PUBLIC :: CHI_PENALTY_TABULATE
  PUBLIC :: CHI_PENALTY_RELEASE
  PUBLIC :: CHI_PENALTY_ASSEMBLE_D
  PUBLIC :: CHI_PENALTY_ASSEMBLE_G
  PUBLIC :: CHI_PENALTY_DIAG
  PUBLIC :: CHI_PENALTY_LUMP
  PUBLIC :: CHI_PENALTY_NODAL
  PUBLIC :: CHI_PENALTY_MATVEC
  PUBLIC :: CHI_PCG3
  PUBLIC :: chi_sum3_iface, chi_filter3_iface, chi_allsum_iface

  INTEGER, PARAMETER :: CHI_PENALTY_NQ = 27

  TYPE tChiPenaltyTab
    INTEGER :: nact = 0                       ! active background elements
    INTEGER, ALLOCATABLE :: iel(:)            ! (nact) element numbers
    REAL*8,  ALLOCATABLE :: w(:,:)            ! (NQ,nact) weight*|J|*beta (0: beta = 0)
    INTEGER, ALLOCATABLE :: body(:,:)         ! (NQ,nact) body index, 0 where beta = 0
    LOGICAL, ALLOCATABLE :: inbody(:,:)       ! (NQ,nact) point inside B_k
    REAL*8,  ALLOCATABLE :: x(:,:,:)          ! (3,NQ,nact) physical points
    REAL*8,  ALLOCATABLE :: uhat(:,:,:)       ! (3,NQ,nact) Dirichlet data (caller)
    REAL*8  :: phi(27, CHI_PENALTY_NQ) = 0d0  ! Q2 basis at the Gauss points
  END TYPE tChiPenaltyTab

  ABSTRACT INTERFACE
    SUBROUTINE chi_sum3_iface(y1, y2, y3)
      REAL*8, INTENT(INOUT) :: y1(*), y2(*), y3(*)
    END SUBROUTINE chi_sum3_iface
    SUBROUTINE chi_filter3_iface(y1, y2, y3, n)
      INTEGER, INTENT(IN) :: n
      REAL*8, INTENT(INOUT) :: y1(n), y2(n), y3(n)
    END SUBROUTINE chi_filter3_iface
    SUBROUTINE chi_allsum_iface(vals, n)
      INTEGER, INTENT(IN) :: n
      REAL*8, INTENT(INOUT) :: vals(n)
    END SUBROUTINE chi_allsum_iface
  END INTERFACE

CONTAINS

  !-----------------------------------------------------------------------
  ! Damping function beta_k of the paper at x for a body with atmosphere
  ! width hwidth; inbody = .TRUE. when x lies in the closed body.
  ! Cylinders are treated as infinite along z (milestone convention).
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_PENALTY_BETA(body, hwidth, x, beta, inbody, rfull, rzero)
    TYPE(tChiBody), INTENT(IN) :: body
    REAL*8, INTENT(IN)  :: hwidth, x(3)
    REAL*8, INTENT(OUT) :: beta
    LOGICAL, INTENT(OUT) :: inbody
    REAL*8, INTENT(IN), OPTIONAL :: rfull, rzero   ! ramp bounds as fractions of H
    REAL*8 :: d(3), r, f1, f0

    d = x - body%center
    IF (body%shape .EQ. CHI_BODY_CYLINDER_Z) THEN
      r = SQRT(d(1)*d(1) + d(2)*d(2))
    ELSE
      r = SQRT(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
    END IF
    f1 = 0.5d0
    f0 = 0.75d0
    IF (PRESENT(rfull)) f1 = rfull
    IF (PRESENT(rzero)) f0 = rzero
    inbody = (r .LE. body%radius)
    IF (inbody) THEN
      beta = 1d0
    ELSE IF (hwidth .LE. 0d0 .OR. f0 .LE. f1) THEN
      beta = 0d0
    ELSE
      beta = (body%radius + f0*hwidth - r) / ((f0 - f1)*hwidth)
      beta = MIN(1d0, MAX(0d0, beta))
    END IF
  END SUBROUTINE CHI_PENALTY_BETA

  !-----------------------------------------------------------------------
  ! Tabulate the quadrature data of all active elements of one mesh level.
  ! The caller fills tab%uhat (finest level) before assembling g.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_PENALTY_TABULATE(nel, kvert, dcorvg, nbody, bodies, hwidth, tab, &
                                  rfull, rzero)
    INTEGER, INTENT(IN) :: nel, nbody
    INTEGER, INTENT(IN) :: kvert(8,*)
    REAL*8,  INTENT(IN) :: dcorvg(3,*), hwidth(nbody)
    TYPE(tChiBody), INTENT(IN) :: bodies(nbody)
    TYPE(tChiPenaltyTab), INTENT(INOUT) :: tab
    REAL*8,  INTENT(IN) :: rfull, rzero

    REAL*8 :: gp(3, CHI_PENALTY_NQ), gw(CHI_PENALTY_NQ), dphi(3,27)
    REAL*8 :: nodes(3,8), xq(3), jac(3,3), detj, beta, bmax
    REAL*8 :: lo(3), hi(3), reach
    REAL*8 :: wq(CHI_PENALTY_NQ), xtmp(3, CHI_PENALTY_NQ)
    INTEGER :: bq(CHI_PENALTY_NQ)
    LOGICAL :: inq(CHI_PENALTY_NQ), inb, near, anyact
    INTEGER :: e, q, k, i, nact, pass
    INTEGER, ALLOCATABLE :: actlist(:)

    CALL CHI_PENALTY_RELEASE(tab)
    CALL CHI_GAUSS3(gp, gw)
    DO q = 1, CHI_PENALTY_NQ
      CALL CHI_Q2_BASIS(gp(:,q), tab%phi(:,q), dphi)
    END DO

    ALLOCATE(actlist(MAX(nel,1)))
    nact = 0
    ! pass 1: count/collect, pass 2: fill
    DO pass = 1, 2
      IF (pass .EQ. 2) THEN
        tab%nact = nact
        ALLOCATE(tab%iel(MAX(nact,1)), tab%w(CHI_PENALTY_NQ, MAX(nact,1)), &
                 tab%body(CHI_PENALTY_NQ, MAX(nact,1)), &
                 tab%inbody(CHI_PENALTY_NQ, MAX(nact,1)), &
                 tab%x(3, CHI_PENALTY_NQ, MAX(nact,1)), &
                 tab%uhat(3, CHI_PENALTY_NQ, MAX(nact,1)))
        tab%uhat = 0d0
      END IF
      DO i = 1, MERGE(nel, nact, pass .EQ. 1)
        e = MERGE(i, actlist(MIN(i, MAX(nact,1))), pass .EQ. 1)
        DO k = 1, 8
          nodes(:,k) = dcorvg(:, kvert(k,e))
        END DO
        IF (pass .EQ. 1) THEN
          ! bounding-box reject against every body's atmosphere
          lo = MINVAL(nodes, DIM=2)
          hi = MAXVAL(nodes, DIM=2)
          near = .FALSE.
          DO k = 1, nbody
            reach = bodies(k)%radius + hwidth(k)
            IF (bodies(k)%shape .EQ. CHI_BODY_CYLINDER_Z) THEN
              near = near .OR. (lo(1) .LE. bodies(k)%center(1)+reach .AND. &
                                hi(1) .GE. bodies(k)%center(1)-reach .AND. &
                                lo(2) .LE. bodies(k)%center(2)+reach .AND. &
                                hi(2) .GE. bodies(k)%center(2)-reach)
            ELSE
              near = near .OR. ALL(lo .LE. bodies(k)%center+reach) .AND. &
                               ALL(hi .GE. bodies(k)%center-reach)
            END IF
          END DO
          IF (.NOT. near) CYCLE
        END IF
        anyact = .FALSE.
        DO q = 1, CHI_PENALTY_NQ
          CALL CHI_Q1_MAP(nodes, gp(:,q), xq, jac, detj)
          bmax = 0d0
          bq(q) = 0
          inq(q) = .FALSE.
          DO k = 1, nbody
            CALL CHI_PENALTY_BETA(bodies(k), hwidth(k), xq, beta, inb, rfull, rzero)
            IF (beta .GT. bmax) THEN
              bmax = beta
              bq(q) = k
              inq(q) = inb
            END IF
          END DO
          wq(q) = gw(q)*ABS(detj)*bmax
          xtmp(:,q) = xq
          IF (bmax .GT. 0d0) anyact = .TRUE.
        END DO
        IF (pass .EQ. 1) THEN
          IF (anyact) THEN
            nact = nact + 1
            actlist(nact) = e
          END IF
        ELSE
          tab%iel(i) = e
          tab%w(:,i) = wq
          tab%body(:,i) = bq
          tab%inbody(:,i) = inq
          tab%x(:,:,i) = xtmp
        END IF
      END DO
    END DO
    DEALLOCATE(actlist)
  END SUBROUTINE CHI_PENALTY_TABULATE

  SUBROUTINE CHI_PENALTY_RELEASE(tab)
    TYPE(tChiPenaltyTab), INTENT(INOUT) :: tab
    tab%nact = 0
    IF (ALLOCATED(tab%iel)) DEALLOCATE(tab%iel)
    IF (ALLOCATED(tab%w)) DEALLOCATE(tab%w)
    IF (ALLOCATED(tab%body)) DEALLOCATE(tab%body)
    IF (ALLOCATED(tab%inbody)) DEALLOCATE(tab%inbody)
    IF (ALLOCATED(tab%x)) DEALLOCATE(tab%x)
    IF (ALLOCATED(tab%uhat)) DEALLOCATE(tab%uhat)
  END SUBROUTINE CHI_PENALTY_RELEASE

  !-----------------------------------------------------------------------
  ! D on the CSR pattern (LdA, ColA) of the level; dvals is overwritten.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_PENALTY_ASSEMBLE_D(tab, gamma, kvert, kedge, karea, nvt, net, nat, &
                                    LdA, ColA, nu, dvals)
    TYPE(tChiPenaltyTab), INTENT(IN) :: tab
    REAL*8,  INTENT(IN) :: gamma
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*), nvt, net, nat
    INTEGER, INTENT(IN) :: LdA(*), ColA(*), nu
    REAL*8,  INTENT(OUT) :: dvals(*)

    INTEGER :: a, q, i, j, idx(27), pos(27,27), row, col, p, na
    REAL*8  :: wq, pi

    na = LdA(nu+1) - 1
    dvals(1:na) = 0d0
    DO a = 1, tab%nact
      CALL CHI_Q2_DOFMAP(tab%iel(a), kvert, kedge, karea, nvt, net, nat, idx)
      ! element-local position map into the CSR rows
      DO i = 1, 27
        row = idx(i)
        DO j = 1, 27
          col = idx(j)
          pos(i,j) = 0
          DO p = LdA(row), LdA(row+1)-1
            IF (ColA(p) .EQ. col) THEN
              pos(i,j) = p
              EXIT
            END IF
          END DO
          IF (pos(i,j) .EQ. 0) THEN
            WRITE(*,'(A,I0,A,I0,A,I0)') 'CHI_PENALTY error: pattern miss, row ', &
              row, ' col ', col, ' element ', tab%iel(a)
            STOP 1
          END IF
        END DO
      END DO
      DO q = 1, CHI_PENALTY_NQ
        wq = gamma*tab%w(q,a)
        IF (wq .EQ. 0d0) CYCLE
        DO i = 1, 27
          pi = wq*tab%phi(i,q)
          DO j = 1, 27
            dvals(pos(i,j)) = dvals(pos(i,j)) + pi*tab%phi(j,q)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE CHI_PENALTY_ASSEMBLE_D

  !-----------------------------------------------------------------------
  ! g from the tabulated Dirichlet data; g1..g3 are overwritten.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_PENALTY_ASSEMBLE_G(tab, gamma, kvert, kedge, karea, nvt, net, nat, &
                                    nu, g1, g2, g3)
    TYPE(tChiPenaltyTab), INTENT(IN) :: tab
    REAL*8,  INTENT(IN) :: gamma
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*), nvt, net, nat, nu
    REAL*8,  INTENT(OUT) :: g1(nu), g2(nu), g3(nu)

    INTEGER :: a, q, i, idx(27)
    REAL*8  :: wq, pi

    g1 = 0d0; g2 = 0d0; g3 = 0d0
    DO a = 1, tab%nact
      CALL CHI_Q2_DOFMAP(tab%iel(a), kvert, kedge, karea, nvt, net, nat, idx)
      DO q = 1, CHI_PENALTY_NQ
        wq = gamma*tab%w(q,a)
        IF (wq .EQ. 0d0) CYCLE
        DO i = 1, 27
          pi = wq*tab%phi(i,q)
          g1(idx(i)) = g1(idx(i)) + pi*tab%uhat(1,q,a)
          g2(idx(i)) = g2(idx(i)) + pi*tab%uhat(2,q,a)
          g3(idx(i)) = g3(idx(i)) + pi*tab%uhat(3,q,a)
        END DO
      END DO
    END DO
  END SUBROUTINE CHI_PENALTY_ASSEMBLE_G

  ! Diagonal of a CSR matrix in FEAT layout (first entry of the row).
  SUBROUTINE CHI_PENALTY_DIAG(LdA, nu, dvals, diag)
    INTEGER, INTENT(IN) :: LdA(*), nu
    REAL*8,  INTENT(IN) :: dvals(*)
    REAL*8,  INTENT(OUT) :: diag(nu)
    INTEGER :: i
    DO i = 1, nu
      diag(i) = dvals(LdA(i))
    END DO
  END SUBROUTINE CHI_PENALTY_DIAG

  ! Row-sum lumping D_L(i) = sum_j d_ij.  NOTE: with a spatially varying
  ! beta the row sums of a Q2 mass matrix can be slightly negative (the Q2
  ! basis changes sign), so this is a diagnostic; the production lumped
  ! operator is the nodal-quadrature one below (CHI_PENALTY_NODAL).
  SUBROUTINE CHI_PENALTY_LUMP(LdA, nu, dvals, dl)
    INTEGER, INTENT(IN) :: LdA(*), nu
    REAL*8,  INTENT(IN) :: dvals(*)
    REAL*8,  INTENT(OUT) :: dl(nu)
    INTEGER :: i, p
    REAL*8  :: s
    DO i = 1, nu
      s = 0d0
      DO p = LdA(i), LdA(i+1)-1
        s = s + dvals(p)
      END DO
      dl(i) = s
    END DO
  END SUBROUTINE CHI_PENALTY_LUMP

  !-----------------------------------------------------------------------
  ! Nodal-quadrature (Lobatto 3x3x3 = the Q2 nodes) lumped penalty:
  !   D_L(i) = gamma * sum_{e owning i} w_i |J_e(xi_i)| beta(x_i),
  ! positive by construction and exact for constants when paired with
  ! g_i = D_L(i) * uhat(x_i).  Returns per node the body index (0 where
  ! beta = 0), the inside-body flag and the node coordinates (for the
  ! donor search of the finest level).  dl is accumulated (caller zeroes).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_PENALTY_NODAL(nel, nvt, net, nat, kvert, kedge, karea, dcorvg, &
                               nbody, bodies, hwidth, gamma, nu, dl, nbody_of, &
                               inbody, xnode, rfull, rzero)
    INTEGER, INTENT(IN) :: nel, nvt, net, nat, nbody, nu
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*)
    REAL*8,  INTENT(IN) :: dcorvg(3,*), hwidth(nbody), gamma
    REAL*8,  INTENT(IN) :: rfull, rzero
    TYPE(tChiBody), INTENT(IN) :: bodies(nbody)
    REAL*8,  INTENT(OUT) :: dl(nu), xnode(3,nu)
    INTEGER, INTENT(OUT) :: nbody_of(nu)
    LOGICAL, INTENT(OUT) :: inbody(nu)

    REAL*8 :: refxi(3,27), w1d(-1:1), wn(27), nodes(3,8), x(3), jac(3,3), detj
    REAL*8 :: beta, bmax, lo(3), hi(3), reach
    INTEGER :: e, i, k, idx(27), ib
    LOGICAL :: inb, near, inmax

    w1d(-1) = 1d0/3d0; w1d(0) = 4d0/3d0; w1d(1) = 1d0/3d0
    CALL CHI_Q2_REFNODES(refxi)
    DO i = 1, 27
      wn(i) = w1d(NINT(refxi(1,i)))*w1d(NINT(refxi(2,i)))*w1d(NINT(refxi(3,i)))
    END DO
    dl = 0d0
    nbody_of = 0
    inbody = .FALSE.
    xnode = 0d0
    DO e = 1, nel
      DO k = 1, 8
        nodes(:,k) = dcorvg(:, kvert(k,e))
      END DO
      lo = MINVAL(nodes, DIM=2)
      hi = MAXVAL(nodes, DIM=2)
      near = .FALSE.
      DO k = 1, nbody
        reach = bodies(k)%radius + hwidth(k)
        IF (bodies(k)%shape .EQ. CHI_BODY_CYLINDER_Z) THEN
          near = near .OR. (lo(1) .LE. bodies(k)%center(1)+reach .AND. &
                            hi(1) .GE. bodies(k)%center(1)-reach .AND. &
                            lo(2) .LE. bodies(k)%center(2)+reach .AND. &
                            hi(2) .GE. bodies(k)%center(2)-reach)
        ELSE
          near = near .OR. ALL(lo .LE. bodies(k)%center+reach) .AND. &
                           ALL(hi .GE. bodies(k)%center-reach)
        END IF
      END DO
      IF (.NOT. near) CYCLE
      CALL CHI_Q2_DOFMAP(e, kvert, kedge, karea, nvt, net, nat, idx)
      DO i = 1, 27
        CALL CHI_Q1_MAP(nodes, refxi(:,i), x, jac, detj)
        bmax = 0d0
        ib = 0
        inmax = .FALSE.
        DO k = 1, nbody
          CALL CHI_PENALTY_BETA(bodies(k), hwidth(k), x, beta, inb, rfull, rzero)
          IF (beta .GT. bmax) THEN
            bmax = beta
            ib = k
            inmax = inb
          END IF
        END DO
        IF (bmax .LE. 0d0) CYCLE
        dl(idx(i)) = dl(idx(i)) + gamma*wn(i)*ABS(detj)*bmax
        nbody_of(idx(i)) = ib
        inbody(idx(i)) = inmax
        xnode(:,idx(i)) = x
      END DO
    END DO
  END SUBROUTINE CHI_PENALTY_NODAL

  ! y = y + coef * D x
  SUBROUTINE CHI_PENALTY_MATVEC(LdA, ColA, nu, dvals, coef, x, y)
    INTEGER, INTENT(IN) :: LdA(*), ColA(*), nu
    REAL*8,  INTENT(IN) :: dvals(*), coef, x(nu)
    REAL*8,  INTENT(INOUT) :: y(nu)
    INTEGER :: i, p
    REAL*8  :: s
    DO i = 1, nu
      s = 0d0
      DO p = LdA(i), LdA(i+1)-1
        s = s + dvals(p)*x(ColA(p))
      END DO
      y(i) = y(i) + coef*s
    END DO
  END SUBROUTINE CHI_PENALTY_MATVEC

  !-----------------------------------------------------------------------
  ! Jacobi-preconditioned CG for [ML + coef D] x = rhs, three components.
  !   ml     global (assembly-summed) lumped mass, dglob = ml + coef*diag(D)
  !          assembly-summed; dvals is the rank-partial D
  !   wts    partition-of-unity weights (1/share count) for scalar products
  !   sum3   assembly sum of partial vectors (E013Sum3); filter3 the
  !          Dirichlet defect filter; allsum the communicator-wide sum
  ! rhs must already be assembly-summed and filtered; x starts at 0.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_PCG3(nu, LdA, ColA, dvals, coef, ml, dglob, wts, &
                      rhs1, rhs2, rhs3, x1, x2, x3, tol, maxit, &
                      sum3, filter3, allsum, iters, resid)
    INTEGER, INTENT(IN) :: nu, LdA(*), ColA(*), maxit
    REAL*8,  INTENT(IN) :: dvals(*), coef, ml(nu), dglob(nu), wts(nu)
    REAL*8,  INTENT(IN) :: rhs1(nu), rhs2(nu), rhs3(nu), tol
    REAL*8,  INTENT(OUT) :: x1(nu), x2(nu), x3(nu)
    PROCEDURE(chi_sum3_iface)    :: sum3
    PROCEDURE(chi_filter3_iface) :: filter3
    PROCEDURE(chi_allsum_iface)  :: allsum
    INTEGER, INTENT(OUT) :: iters
    REAL*8,  INTENT(OUT) :: resid

    REAL*8, ALLOCATABLE :: r(:,:), z(:,:), p(:,:), q(:,:), t(:,:)
    REAL*8  :: rz(3), rznew(3), pq(3), alpha(3), beta(3), rnorm(3), rhsn(3), s(6)
    LOGICAL :: done(3)
    INTEGER :: c, it, i

    ALLOCATE(r(nu,3), z(nu,3), p(nu,3), q(nu,3), t(nu,3))
    x1 = 0d0; x2 = 0d0; x3 = 0d0
    r(:,1) = rhs1; r(:,2) = rhs2; r(:,3) = rhs3
    DO c = 1, 3
      s(c) = wdot(r(:,c), r(:,c))
    END DO
    CALL allsum(s(1:3), 3)
    rhsn = SQRT(MAX(s(1:3), 0d0))
    rnorm = rhsn
    iters = 0
    resid = MAXVAL(rnorm)
    done = (rnorm .LE. tol*MAX(rhsn, TINY(1d0)))
    IF (ALL(done)) THEN
      DEALLOCATE(r, z, p, q, t)
      RETURN
    END IF
    DO c = 1, 3
      z(:,c) = r(:,c)/dglob
    END DO
    CALL filter3(z(:,1), z(:,2), z(:,3), nu)
    p = z
    DO c = 1, 3
      s(c) = wdot(r(:,c), z(:,c))
    END DO
    CALL allsum(s(1:3), 3)
    rz = s(1:3)

    DO it = 1, maxit
      ! q = [ML + coef D] p  (D partial -> assembly sum), filtered
      t = 0d0
      DO c = 1, 3
        CALL CHI_PENALTY_MATVEC(LdA, ColA, nu, dvals, coef, p(:,c), t(:,c))
      END DO
      CALL sum3(t(:,1), t(:,2), t(:,3))
      DO c = 1, 3
        q(:,c) = ml*p(:,c) + t(:,c)
      END DO
      CALL filter3(q(:,1), q(:,2), q(:,3), nu)
      DO c = 1, 3
        s(c) = wdot(p(:,c), q(:,c))
      END DO
      CALL allsum(s(1:3), 3)
      pq = s(1:3)
      DO c = 1, 3
        IF (done(c) .OR. pq(c) .LE. 0d0) THEN
          alpha(c) = 0d0
        ELSE
          alpha(c) = rz(c)/pq(c)
        END IF
      END DO
      x1 = x1 + alpha(1)*p(:,1)
      x2 = x2 + alpha(2)*p(:,2)
      x3 = x3 + alpha(3)*p(:,3)
      DO c = 1, 3
        r(:,c) = r(:,c) - alpha(c)*q(:,c)
        z(:,c) = r(:,c)/dglob
      END DO
      CALL filter3(z(:,1), z(:,2), z(:,3), nu)
      DO c = 1, 3
        s(c)   = wdot(r(:,c), r(:,c))
        s(3+c) = wdot(r(:,c), z(:,c))
      END DO
      CALL allsum(s, 6)
      rnorm = SQRT(MAX(s(1:3), 0d0))
      rznew = s(4:6)
      iters = it
      DO c = 1, 3
        IF (.NOT. done(c)) done(c) = (rnorm(c) .LE. tol*MAX(rhsn(c), TINY(1d0)))
      END DO
      resid = 0d0
      DO c = 1, 3
        resid = MAX(resid, rnorm(c)/MAX(rhsn(c), TINY(1d0)))
      END DO
      IF (ALL(done)) EXIT
      DO c = 1, 3
        IF (done(c) .OR. rz(c) .EQ. 0d0) THEN
          beta(c) = 0d0
        ELSE
          beta(c) = rznew(c)/rz(c)
        END IF
        p(:,c) = z(:,c) + beta(c)*p(:,c)
      END DO
      rz = rznew
    END DO
    DEALLOCATE(r, z, p, q, t)

  CONTAINS
    REAL*8 FUNCTION wdot(a, b)
      REAL*8, INTENT(IN) :: a(nu), b(nu)
      INTEGER :: ii
      wdot = 0d0
      DO ii = 1, nu
        wdot = wdot + wts(ii)*a(ii)*b(ii)
      END DO
    END FUNCTION wdot
  END SUBROUTINE CHI_PCG3

END MODULE CHI_PENALTY
