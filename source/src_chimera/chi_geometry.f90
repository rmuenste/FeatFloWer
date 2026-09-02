!=========================================================================
! CHI_GEOMETRY - dependency-free geometric kernel of the Chimera
! component (design: chimera-integration-design.md v3, Layer M).
!
! Contents: the Q1 trilinear forward map of a FeatFloWer hexahedron
! (functional twin of EL_Q1_MAP in source/src_el/el_quadrature.f90,
! pinned against it by test_chi_geometry - el_quadrature itself is NOT
! reused here because it drags in EL_HALO/EL_FIELDS/MPI state), a damped
! Newton inverse map (idiom of GetPointFromElement,
! source/src_particles/part_step.f90), a 3x3x3 Gauss cubature rule and a
! 3x3 inverse.
!
! Vertex ordering convention (FeatFloWer/FEAT kvert):
!   1:(-,-,-) 2:(+,-,-) 3:(+,+,-) 4:(-,+,-)
!   5:(-,-,+) 6:(+,-,+) 7:(+,+,+) 8:(-,+,+)
! in reference coordinates xi in [-1,1]^3.
!
! No COMMON blocks, no module state, no USE of solver modules: every
! routine is pure in the reentrancy sense (all data through arguments).
!=========================================================================
MODULE CHI_GEOMETRY

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CHI_Q1_MAP
  PUBLIC :: CHI_INVERSE_MAP
  PUBLIC :: CHI_GAUSS3
  PUBLIC :: CHI_M33INV

  ! Reference-coordinate signs of the 8 hex vertices (FEAT ordering).
  REAL*8, PARAMETER :: chi_vertex_signs(3,8) = RESHAPE((/ &
    -1d0,-1d0,-1d0,  1d0,-1d0,-1d0,  1d0, 1d0,-1d0, -1d0, 1d0,-1d0, &
    -1d0,-1d0, 1d0,  1d0,-1d0, 1d0,  1d0, 1d0, 1d0, -1d0, 1d0, 1d0 /), (/3,8/))

  INTEGER, PARAMETER :: chi_newton_maxiter = 50
  REAL*8,  PARAMETER :: chi_newton_steptol = 1d-13
  REAL*8,  PARAMETER :: chi_newton_divergebound = 1.5d0

CONTAINS

  !-----------------------------------------------------------------------
  ! Trilinear (Q1) forward map: reference point xi -> physical point,
  ! Jacobian d(point)/d(xi) and its determinant.  nodes(1:3,1:8) are the
  ! vertex coordinates in FEAT ordering.
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_Q1_MAP(nodes, xi, point, jac, detj)
    REAL*8, INTENT(IN)  :: nodes(3,8), xi(3)
    REAL*8, INTENT(OUT) :: point(3), jac(3,3), detj

    REAL*8 :: shp, der(3)
    INTEGER :: i

    point = 0d0
    jac = 0d0
    DO i = 1, 8
      shp = 0.125d0*(1d0+chi_vertex_signs(1,i)*xi(1)) &
                   *(1d0+chi_vertex_signs(2,i)*xi(2)) &
                   *(1d0+chi_vertex_signs(3,i)*xi(3))
      der(1) = 0.125d0*chi_vertex_signs(1,i) &
                   *(1d0+chi_vertex_signs(2,i)*xi(2)) &
                   *(1d0+chi_vertex_signs(3,i)*xi(3))
      der(2) = 0.125d0*chi_vertex_signs(2,i) &
                   *(1d0+chi_vertex_signs(1,i)*xi(1)) &
                   *(1d0+chi_vertex_signs(3,i)*xi(3))
      der(3) = 0.125d0*chi_vertex_signs(3,i) &
                   *(1d0+chi_vertex_signs(1,i)*xi(1)) &
                   *(1d0+chi_vertex_signs(2,i)*xi(2))
      point = point + shp*nodes(:,i)
      jac(:,1) = jac(:,1) + der(1)*nodes(:,i)
      jac(:,2) = jac(:,2) + der(2)*nodes(:,i)
      jac(:,3) = jac(:,3) + der(3)*nodes(:,i)
    END DO
    detj = jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3)) &
         - jac(2,1)*(jac(1,2)*jac(3,3)-jac(3,2)*jac(1,3)) &
         + jac(3,1)*(jac(1,2)*jac(2,3)-jac(2,2)*jac(1,3))
  END SUBROUTINE CHI_Q1_MAP

  !-----------------------------------------------------------------------
  ! Newton inverse of the trilinear map: physical point p -> reference
  ! coordinates xi.  converged = .TRUE. iff the iteration converged AND
  ! the point lies inside the reference cube (|xi_i| <= 1 + tol).
  ! Divergence guard: |xi| beyond 1.5 after a few iterations aborts, as
  ! in GetPointFromElement.
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_INVERSE_MAP(nodes, p, xi, converged, inside_tol)
    REAL*8, INTENT(IN)  :: nodes(3,8), p(3)
    REAL*8, INTENT(OUT) :: xi(3)
    LOGICAL, INTENT(OUT) :: converged
    REAL*8, INTENT(IN), OPTIONAL :: inside_tol

    REAL*8 :: point(3), jac(3,3), jacinv(3,3), detj, res(3), step(3)
    REAL*8 :: tol_inside, stepnorm2
    LOGICAL :: ok
    INTEGER :: iter

    tol_inside = 1d-2
    IF (PRESENT(inside_tol)) tol_inside = inside_tol

    xi = 0d0
    converged = .FALSE.

    DO iter = 1, chi_newton_maxiter
      CALL CHI_Q1_MAP(nodes, xi, point, jac, detj)
      CALL CHI_M33INV(jac, jacinv, ok)
      IF (.NOT. ok) RETURN
      res = p - point
      step = MATMUL(jacinv, res)
      xi = xi + step
      stepnorm2 = step(1)*step(1) + step(2)*step(2) + step(3)*step(3)
      IF (iter .GT. 4 .AND. &
          (ABS(xi(1)) .GT. chi_newton_divergebound .OR. &
           ABS(xi(2)) .GT. chi_newton_divergebound .OR. &
           ABS(xi(3)) .GT. chi_newton_divergebound)) RETURN
      IF (stepnorm2 .LT. chi_newton_steptol*chi_newton_steptol) EXIT
    END DO

    IF (ABS(xi(1)) .GT. 1d0+tol_inside .OR. &
        ABS(xi(2)) .GT. 1d0+tol_inside .OR. &
        ABS(xi(3)) .GT. 1d0+tol_inside) RETURN

    converged = .TRUE.
  END SUBROUTINE CHI_INVERSE_MAP

  !-----------------------------------------------------------------------
  ! 3x3x3 tensor-product Gauss rule on [-1,1]^3 (exact through degree 5
  ! per direction).  pts(1:3,1:27), wts(1:27); sum of weights = 8.
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_GAUSS3(pts, wts)
    REAL*8, INTENT(OUT) :: pts(3,27), wts(27)

    REAL*8, PARAMETER :: g = 0.774596669241483377035853079956d0  ! sqrt(3/5)
    REAL*8 :: p1(3), w1(3)
    INTEGER :: i, j, k, n

    p1 = (/ -g, 0d0, g /)
    w1 = (/ 5d0/9d0, 8d0/9d0, 5d0/9d0 /)

    n = 0
    DO k = 1, 3
      DO j = 1, 3
        DO i = 1, 3
          n = n + 1
          pts(1,n) = p1(i)
          pts(2,n) = p1(j)
          pts(3,n) = p1(k)
          wts(n) = w1(i)*w1(j)*w1(k)
        END DO
      END DO
    END DO
  END SUBROUTINE CHI_GAUSS3

  !-----------------------------------------------------------------------
  ! Inverse of a 3x3 matrix; ok = .FALSE. on (near-)singular input.
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_M33INV(a, ainv, ok)
    REAL*8, INTENT(IN)  :: a(3,3)
    REAL*8, INTENT(OUT) :: ainv(3,3)
    LOGICAL, INTENT(OUT) :: ok

    REAL*8 :: det, cof(3,3)

    cof(1,1) =  (a(2,2)*a(3,3) - a(2,3)*a(3,2))
    cof(1,2) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))
    cof(1,3) =  (a(2,1)*a(3,2) - a(2,2)*a(3,1))
    cof(2,1) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))
    cof(2,2) =  (a(1,1)*a(3,3) - a(1,3)*a(3,1))
    cof(2,3) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))
    cof(3,1) =  (a(1,2)*a(2,3) - a(1,3)*a(2,2))
    cof(3,2) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))
    cof(3,3) =  (a(1,1)*a(2,2) - a(1,2)*a(2,1))

    det = a(1,1)*cof(1,1) + a(1,2)*cof(1,2) + a(1,3)*cof(1,3)

    IF (ABS(det) .LT. TINY(1d0)*64d0) THEN
      ok = .FALSE.
      ainv = 0d0
      RETURN
    END IF

    ainv = TRANSPOSE(cof)/det
    ok = .TRUE.
  END SUBROUTINE CHI_M33INV

END MODULE CHI_GEOMETRY
