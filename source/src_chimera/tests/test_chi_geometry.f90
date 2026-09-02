!=========================================================================
! test_chi_geometry - Phase 1 unit test of CHI_GEOMETRY
! (chimera-integration-design.md v3, section 9).
!
! 1. Pins CHI_Q1_MAP against the established EL_Q1_MAP (EL_QUADRATURE)
!    on distorted hexes and random reference points: any drift between
!    the twin implementations fails the test.
! 2. Inverse-map round trip: xi -> forward -> CHI_INVERSE_MAP -> xi.
! 3. Cubature sanity: weights sum to 8; integrating 1 over the unit cube
!    via detj gives volume 1; a trilinear integrand is integrated
!    exactly on a distorted hex.
!
! Serial; exits with status 0 and prints PASS on success, STOP 1 on any
! failure.
!=========================================================================
PROGRAM test_chi_geometry

  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP, CHI_INVERSE_MAP, CHI_GAUSS3, CHI_M33INV
  USE EL_QUADRATURE, ONLY: EL_Q1_MAP

  IMPLICIT NONE

  REAL*8 :: nodes(3,8), x8(8), y8(8), z8(8)
  REAL*8 :: xi(3), xi2(3), p1(3), p2(3), j1(3,3), j2(3,3), d1, d2
  REAL*8 :: pts(3,27), wts(27), vol, fint, fref
  REAL*8 :: point(3), jac(3,3), detj
  INTEGER :: itest, i, k
  LOGICAL :: conv
  INTEGER :: nfail
  INTEGER*8 :: seed

  REAL*8, PARAMETER :: signs(3,8) = RESHAPE((/ &
    -1d0,-1d0,-1d0,  1d0,-1d0,-1d0,  1d0, 1d0,-1d0, -1d0, 1d0,-1d0, &
    -1d0,-1d0, 1d0,  1d0,-1d0, 1d0,  1d0, 1d0, 1d0, -1d0, 1d0, 1d0 /), (/3,8/))

  nfail = 0
  seed = 20260902_8

  !--- 1. pin against EL_Q1_MAP -----------------------------------------
  DO itest = 1, 20
    CALL make_distorted_hex(seed, nodes)
    DO i = 1, 8
      x8(i) = nodes(1,i)
      y8(i) = nodes(2,i)
      z8(i) = nodes(3,i)
    END DO
    DO k = 1, 10
      xi = (/ 2d0*lcg(seed)-1d0, 2d0*lcg(seed)-1d0, 2d0*lcg(seed)-1d0 /)
      CALL CHI_Q1_MAP(nodes, xi, p1, j1, d1)
      CALL EL_Q1_MAP(x8, y8, z8, xi, p2, j2, d2)
      IF (MAXVAL(ABS(p1-p2)) .GT. 1d-14 .OR. &
          MAXVAL(ABS(j1-j2)) .GT. 1d-14 .OR. &
          ABS(d1-d2) .GT. 1d-14) THEN
        WRITE(*,*) 'FAIL: CHI_Q1_MAP deviates from EL_Q1_MAP, test', itest
        nfail = nfail + 1
      END IF
    END DO
  END DO

  !--- 2. inverse-map round trip ----------------------------------------
  DO itest = 1, 20
    CALL make_distorted_hex(seed, nodes)
    DO k = 1, 10
      xi = (/ 1.9d0*lcg(seed)-0.95d0, 1.9d0*lcg(seed)-0.95d0, &
              1.9d0*lcg(seed)-0.95d0 /)
      CALL CHI_Q1_MAP(nodes, xi, point, jac, detj)
      CALL CHI_INVERSE_MAP(nodes, point, xi2, conv)
      IF (.NOT. conv .OR. MAXVAL(ABS(xi-xi2)) .GT. 1d-10) THEN
        WRITE(*,*) 'FAIL: inverse-map round trip, test', itest, &
          ' err=', MAXVAL(ABS(xi-xi2)), ' conv=', conv
        nfail = nfail + 1
      END IF
    END DO
    ! a clearly outside point must not report converged-inside
    point = nodes(:,1) + 5d0*(nodes(:,1) - nodes(:,7))
    CALL CHI_INVERSE_MAP(nodes, point, xi2, conv)
    IF (conv) THEN
      WRITE(*,*) 'FAIL: outside point classified inside, test', itest
      nfail = nfail + 1
    END IF
  END DO

  !--- 3. cubature ------------------------------------------------------
  CALL CHI_GAUSS3(pts, wts)
  IF (ABS(SUM(wts) - 8d0) .GT. 1d-13) THEN
    WRITE(*,*) 'FAIL: Gauss weights sum to', SUM(wts)
    nfail = nfail + 1
  END IF

  ! unit cube [0,1]^3: volume 1
  DO i = 1, 8
    nodes(:,i) = 0.5d0*(signs(:,i) + 1d0)
  END DO
  vol = 0d0
  DO k = 1, 27
    CALL CHI_Q1_MAP(nodes, pts(:,k), point, jac, detj)
    vol = vol + wts(k)*ABS(detj)
  END DO
  IF (ABS(vol - 1d0) .GT. 1d-13) THEN
    WRITE(*,*) 'FAIL: unit-cube volume =', vol
    nfail = nfail + 1
  END IF

  ! trilinear integrand on a distorted hex: integral of f(x)=x*y*z over
  ! the image equals the cubature result computed at high exactness -
  ! compare degree-3-in-each-variable integrand (within Gauss-3 range)
  ! against itself under a refined 2x2x2 subdivision of the reference
  ! cube (consistency check of map + rule).
  CALL make_distorted_hex(seed, nodes)
  fint = 0d0
  DO k = 1, 27
    CALL CHI_Q1_MAP(nodes, pts(:,k), point, jac, detj)
    fint = fint + wts(k)*ABS(detj)*point(1)*point(2)*point(3)
  END DO
  fref = subdivided_integral(nodes, pts, wts)
  IF (ABS(fint - fref) .GT. 1d-11*MAX(1d0, ABS(fref))) THEN
    WRITE(*,*) 'FAIL: trilinear integrand mismatch', fint, fref
    nfail = nfail + 1
  END IF

  IF (nfail .GT. 0) THEN
    WRITE(*,*) 'test_chi_geometry: ', nfail, ' failure(s)'
    STOP 1
  END IF
  WRITE(*,*) 'test_chi_geometry: PASS'

CONTAINS

  ! Deterministic uniform(0,1) LCG (Numerical Recipes constants).
  FUNCTION lcg(s) RESULT(r)
    INTEGER*8, INTENT(INOUT) :: s
    REAL*8 :: r
    ! Lehmer/minstd: product stays < 2^47, no signed-overflow UB
    s = MOD(s*48271_8, 2147483647_8)
    r = DBLE(s)/2147483647d0
  END FUNCTION lcg

  ! Unit-cube-based hex with deterministic vertex perturbations that keep
  ! it convex and well-shaped.
  SUBROUTINE make_distorted_hex(s, nds)
    INTEGER*8, INTENT(INOUT) :: s
    REAL*8, INTENT(OUT) :: nds(3,8)
    INTEGER :: iv
    DO iv = 1, 8
      nds(:,iv) = 0.5d0*(signs(:,iv) + 1d0)
      nds(1,iv) = nds(1,iv) + 0.15d0*(lcg(s)-0.5d0)
      nds(2,iv) = nds(2,iv) + 0.15d0*(lcg(s)-0.5d0)
      nds(3,iv) = nds(3,iv) + 0.15d0*(lcg(s)-0.5d0)
    END DO
  END SUBROUTINE make_distorted_hex

  ! Integral of x*y*z over the mapped hex via 8 reference-subcube
  ! applications of the same Gauss rule (independent path through the
  ! same map; agrees when map+rule are consistent).
  FUNCTION subdivided_integral(nds, gp, gw) RESULT(total)
    REAL*8, INTENT(IN) :: nds(3,8), gp(3,27), gw(27)
    REAL*8 :: total
    REAL*8 :: xi(3), pnt(3), jc(3,3), dj, c(3)
    INTEGER :: sx, sy, sz, kk
    total = 0d0
    DO sz = 0, 1
      DO sy = 0, 1
        DO sx = 0, 1
          c = (/ -0.5d0 + DBLE(sx), -0.5d0 + DBLE(sy), -0.5d0 + DBLE(sz) /)
          DO kk = 1, 27
            xi = c + 0.5d0*gp(:,kk)
            CALL CHI_Q1_MAP(nds, xi, pnt, jc, dj)
            total = total + 0.125d0*gw(kk)*ABS(dj)*pnt(1)*pnt(2)*pnt(3)
          END DO
        END DO
      END DO
    END DO
  END FUNCTION subdivided_integral

END PROGRAM test_chi_geometry
