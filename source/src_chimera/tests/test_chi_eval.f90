!=========================================================================
! test_chi_eval - Phase 1 unit test of CHI_FEM_EVAL
! (chimera-integration-design.md v3, section 9).
!
! 1. Delta property: phi_j(xi_i) = delta_ij at the 27 reference nodes
!    (this pins the local ordering to the RETURN_Velo/E013 convention).
! 2. Partition of unity and zero derivative-sum at random points.
! 3. Exact reproduction of an arbitrary quadratic field (value AND
!    physical gradient) on an AFFINE-mapped element - physical
!    quadratics lie in the Q2 space only for affine geometry, which is
!    what makes the check exact to round-off.
! 4. P1 pressure evaluation.
! 5. CHI_Q2_DOFMAP index layout.
!
! Serial; prints PASS / STOP 1.
!=========================================================================
PROGRAM test_chi_eval

  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP
  USE CHI_FEM_EVAL, ONLY: CHI_Q2_BASIS, CHI_Q2_REFNODES, CHI_Q2_DOFMAP, &
    CHI_EVAL_Q2_SCALAR, CHI_EVAL_Q2_GRADIENT, CHI_EVAL_P1

  IMPLICIT NONE

  REAL*8 :: refxi(3,27), phi(27), dphi(3,27)
  REAL*8 :: nodes(3,8), amat(3,3), b(3)
  REAL*8 :: xi(3), xnode(3,27), vals(27)
  REAL*8 :: point(3), jac(3,3), detj, grad(3), gradref(3)
  REAL*8 :: val, valref, qc(10)
  REAL*8 :: pdofs(4), x(3), xc(3), p, pref
  INTEGER :: i, j, k, itest, nfail
  INTEGER :: kvert(8,2), kedge(12,2), karea(6,2), idx(27)
  LOGICAL :: ok
  INTEGER*8 :: seed

  REAL*8, PARAMETER :: signs(3,8) = RESHAPE((/ &
    -1d0,-1d0,-1d0,  1d0,-1d0,-1d0,  1d0, 1d0,-1d0, -1d0, 1d0,-1d0, &
    -1d0,-1d0, 1d0,  1d0,-1d0, 1d0,  1d0, 1d0, 1d0, -1d0, 1d0, 1d0 /), (/3,8/))

  nfail = 0
  seed = 987654321_8

  !--- 1. delta property -------------------------------------------------
  CALL CHI_Q2_REFNODES(refxi)
  DO i = 1, 27
    CALL CHI_Q2_BASIS(refxi(:,i), phi, dphi)
    DO j = 1, 27
      IF (i .EQ. j) THEN
        IF (ABS(phi(j) - 1d0) .GT. 1d-13) THEN
          WRITE(*,*) 'FAIL: delta property phi_i(xi_i), i=', i, phi(j)
          nfail = nfail + 1
        END IF
      ELSE
        IF (ABS(phi(j)) .GT. 1d-13) THEN
          WRITE(*,*) 'FAIL: delta property phi_j(xi_i), i,j=', i, j, phi(j)
          nfail = nfail + 1
        END IF
      END IF
    END DO
  END DO

  !--- 2. partition of unity --------------------------------------------
  DO itest = 1, 25
    xi = (/ 2d0*lcg(seed)-1d0, 2d0*lcg(seed)-1d0, 2d0*lcg(seed)-1d0 /)
    CALL CHI_Q2_BASIS(xi, phi, dphi)
    IF (ABS(SUM(phi) - 1d0) .GT. 1d-13) THEN
      WRITE(*,*) 'FAIL: partition of unity at', xi, SUM(phi)
      nfail = nfail + 1
    END IF
    DO k = 1, 3
      IF (ABS(SUM(dphi(k,:))) .GT. 1d-12) THEN
        WRITE(*,*) 'FAIL: derivative sum nonzero, dir', k
        nfail = nfail + 1
      END IF
    END DO
  END DO

  !--- 3. quadratic reproduction on an affine element -------------------
  ! affine map x = A*xihat + b with a nonsingular, nonsymmetric A
  amat = RESHAPE((/ 0.8d0, 0.1d0, -0.05d0, &
                    0.2d0, 0.9d0,  0.15d0, &
                   -0.1d0, 0.05d0, 0.7d0 /), (/3,3/))
  b = (/ 0.3d0, -0.2d0, 0.5d0 /)
  DO i = 1, 8
    nodes(:,i) = MATMUL(amat, signs(:,i)) + b
  END DO
  ! physical coordinates of the 27 Q2 nodes (affine: map of ref nodes)
  DO i = 1, 27
    xnode(:,i) = MATMUL(amat, refxi(:,i)) + b
  END DO
  ! arbitrary quadratic f(x) = qc . (1,x,y,z,x^2,y^2,z^2,xy,xz,yz)
  DO i = 1, 10
    qc(i) = 2d0*lcg(seed) - 1d0
  END DO
  DO i = 1, 27
    vals(i) = quad_f(qc, xnode(:,i))
  END DO
  DO itest = 1, 25
    xi = (/ 2d0*lcg(seed)-1d0, 2d0*lcg(seed)-1d0, 2d0*lcg(seed)-1d0 /)
    CALL CHI_Q2_BASIS(xi, phi, dphi)
    CALL CHI_Q1_MAP(nodes, xi, point, jac, detj)
    val = CHI_EVAL_Q2_SCALAR(vals, phi)
    valref = quad_f(qc, point)
    IF (ABS(val - valref) .GT. 1d-12*MAX(1d0, ABS(valref))) THEN
      WRITE(*,*) 'FAIL: quadratic value reproduction', val, valref
      nfail = nfail + 1
    END IF
    CALL CHI_EVAL_Q2_GRADIENT(vals, dphi, jac, grad, ok)
    gradref = quad_grad(qc, point)
    IF (.NOT. ok .OR. MAXVAL(ABS(grad - gradref)) .GT. 1d-11) THEN
      WRITE(*,*) 'FAIL: quadratic gradient reproduction', grad, gradref
      nfail = nfail + 1
    END IF
  END DO

  !--- 4. P1 evaluation --------------------------------------------------
  pdofs = (/ 1.5d0, -0.4d0, 0.7d0, 2.1d0 /)
  xc = (/ 0.2d0, 0.3d0, -0.1d0 /)
  DO itest = 1, 10
    x = (/ lcg(seed), lcg(seed), lcg(seed) /)
    p = CHI_EVAL_P1(pdofs, x, xc)
    pref = pdofs(1) + pdofs(2)*(x(1)-xc(1)) + pdofs(3)*(x(2)-xc(2)) &
                    + pdofs(4)*(x(3)-xc(3))
    IF (ABS(p - pref) .GT. 1d-14) THEN
      WRITE(*,*) 'FAIL: P1 evaluation'
      nfail = nfail + 1
    END IF
  END DO

  !--- 5. DOF map layout -------------------------------------------------
  ! synthetic connectivity of "element 2" in a mesh with nvt=100, net=200,
  ! nat=50
  DO i = 1, 8
    kvert(i,2) = 10 + i
  END DO
  DO i = 1, 12
    kedge(i,2) = 20 + i
  END DO
  DO i = 1, 6
    karea(i,2) = 30 + i
  END DO
  CALL CHI_Q2_DOFMAP(2, kvert, kedge, karea, 100, 200, 50, idx)
  DO i = 1, 8
    IF (idx(i) .NE. 10+i) nfail = nfail + 1
  END DO
  DO i = 1, 12
    IF (idx(8+i) .NE. 100 + 20+i) nfail = nfail + 1
  END DO
  DO i = 1, 6
    IF (idx(20+i) .NE. 100 + 200 + 30+i) nfail = nfail + 1
  END DO
  IF (idx(27) .NE. 100 + 200 + 50 + 2) nfail = nfail + 1

  IF (nfail .GT. 0) THEN
    WRITE(*,*) 'test_chi_eval: ', nfail, ' failure(s)'
    STOP 1
  END IF
  WRITE(*,*) 'test_chi_eval: PASS'

CONTAINS

  FUNCTION lcg(s) RESULT(r)
    INTEGER*8, INTENT(INOUT) :: s
    REAL*8 :: r
    ! Lehmer/minstd: product stays < 2^47, no signed-overflow UB
    s = MOD(s*48271_8, 2147483647_8)
    r = DBLE(s)/2147483647d0
  END FUNCTION lcg

  FUNCTION quad_f(c, pt) RESULT(f)
    REAL*8, INTENT(IN) :: c(10), pt(3)
    REAL*8 :: f
    f = c(1) + c(2)*pt(1) + c(3)*pt(2) + c(4)*pt(3) &
      + c(5)*pt(1)*pt(1) + c(6)*pt(2)*pt(2) + c(7)*pt(3)*pt(3) &
      + c(8)*pt(1)*pt(2) + c(9)*pt(1)*pt(3) + c(10)*pt(2)*pt(3)
  END FUNCTION quad_f

  FUNCTION quad_grad(c, pt) RESULT(g)
    REAL*8, INTENT(IN) :: c(10), pt(3)
    REAL*8 :: g(3)
    g(1) = c(2) + 2d0*c(5)*pt(1) + c(8)*pt(2) + c(9)*pt(3)
    g(2) = c(3) + 2d0*c(6)*pt(2) + c(8)*pt(1) + c(10)*pt(3)
    g(3) = c(4) + 2d0*c(7)*pt(3) + c(9)*pt(1) + c(10)*pt(2)
  END FUNCTION quad_grad

END PROGRAM test_chi_eval
