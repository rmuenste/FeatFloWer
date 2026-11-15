!=========================================================================
! EL_MM_OPERATOR - solid-fraction operators of momentum model A (operator form)
!
! Model A weights the pressure gradient and the viscous stress of the fluid momentum equation
! with the porosity: -eps grad p + div(eps tau). With phi = 1 - eps the weak forms split into
! the solver's own unweighted operators and a solid-fraction part,
!     int p div(eps_h N)        = int p div N        - int p div(phi_h N)
!     int eps_h c grad u:grad N = int c grad u:grad N - int phi_h c grad u:grad N ,
! and this module evaluates the two solid-fraction parts element by element:
!     fp(c,i) = int_e p_h d/dx_c (phi_h N_i) dV        fv(c,i) = int_e phi_h grad u_c . grad N_i dV
! (N_i the Q2 velocity basis, p_h the discontinuous P1 pressure p_1 + p_{2:4}.(x - x_c) with
! x_c the centre dof, phi_h = 1 - eps_h with the continuous transport porosity eps_h of the
! conservative continuity form). The fluid adapter (Transport_Q2P1, QuadSc_el_mm.f90) subtracts
! them from what the unweighted operators contribute; no shared operator is modified.
!
! Properties (unit tests): fp is the transpose of the solid-fraction part of the conservative
! divergence, sum_i u_i . fp_i(p) = sum_m p_m [F_m(u; 1) - F_m(u; eps_h)], so pressure gradient
! and continuity form one saddle point; a constant pressure gives zero force at every interior
! node for any continuous phi_h (gauge invariance); fv annihilates constant velocities.
!
! Geometry and quadrature: isoparametric Q2 map and the 3 x 3 x 3 Gauss rule of EL_VA_FLUX.
!=========================================================================
MODULE EL_MM_OPERATOR

  USE TYPES, ONLY: tMultiMesh
  USE EL_VA_FLUX, ONLY: tElVaFluxTables, EL_VA_FLUX_TABLES, EL_VA_ELEMENT_DOFS, NQ2, NGPV

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: EL_MM_ELEMENT_PHI_OPERATORS, EL_MM_PHI_OPERATORS

  TYPE(tElVaFluxTables), SAVE :: tables

CONTAINS

  ! One element: x(3,27) dof coordinates, phi(27) nodal solid fraction, p(4) pressure dofs,
  ! u(3,27) nodal velocity. fv is computed only when with_viscous is set (zero otherwise).
  PURE SUBROUTINE EL_MM_ELEMENT_PHI_OPERATORS(t, x, phi, p, u, with_viscous, fp, fv)

    TYPE(tElVaFluxTables), INTENT(IN) :: t
    REAL*8, INTENT(IN) :: x(3,NQ2), phi(NQ2), p(4), u(3,NQ2)
    LOGICAL, INTENT(IN) :: with_viscous
    REAL*8, INTENT(OUT) :: fp(3,NQ2), fv(3,NQ2)

    REAL*8 :: jac(3,3), inv(3,3), detj, xq(3), xc(3), dn(3,NQ2), phiq, gphi(3), gu(3,3)
    REAL*8 :: pq, w
    INTEGER :: q, l, k, c

    fp = 0.0d0
    fv = 0.0d0
    xc = x(:,NQ2)
    DO q = 1, NGPV
      jac = 0.0d0
      xq = 0.0d0
      DO l = 1, NQ2
        DO k = 1, 3
          jac(:,k) = jac(:,k) + x(:,l)*t%vder(k,l,q)
        END DO
        xq = xq + x(:,l)*t%vval(l,q)
      END DO
      CALL invert3(jac, inv, detj)
      w = t%vw(q)*ABS(detj)
      phiq = 0.0d0
      gphi = 0.0d0
      gu = 0.0d0
      DO l = 1, NQ2
        ! d N_l / d x_j = sum_k dref(k,l) inv(k,j)
        DO k = 1, 3
          dn(k,l) = t%vder(1,l,q)*inv(1,k) + t%vder(2,l,q)*inv(2,k) + t%vder(3,l,q)*inv(3,k)
        END DO
        phiq = phiq + phi(l)*t%vval(l,q)
        gphi = gphi + phi(l)*dn(:,l)
        IF (with_viscous) THEN
          DO k = 1, 3
            gu(:,k) = gu(:,k) + u(:,l)*dn(k,l)
          END DO
        END IF
      END DO
      pq = p(1) + DOT_PRODUCT(p(2:4), xq - xc)
      DO l = 1, NQ2
        fp(:,l) = fp(:,l) + pq*(phiq*dn(:,l) + t%vval(l,q)*gphi)*w
        IF (with_viscous) THEN
          DO c = 1, 3
            fv(c,l) = fv(c,l) + phiq*DOT_PRODUCT(gu(c,:), dn(:,l))*w
          END DO
        END IF
      END DO
    END DO

  CONTAINS

    PURE SUBROUTINE invert3(a, b, det)
      REAL*8, INTENT(IN) :: a(3,3)
      REAL*8, INTENT(OUT) :: b(3,3), det
      det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
          + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      b(1,1) =  (a(2,2)*a(3,3)-a(2,3)*a(3,2))/det
      b(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))/det
      b(1,3) =  (a(1,2)*a(2,3)-a(1,3)*a(2,2))/det
      b(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))/det
      b(2,2) =  (a(1,1)*a(3,3)-a(1,3)*a(3,1))/det
      b(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))/det
      b(3,1) =  (a(2,1)*a(3,2)-a(2,2)*a(3,1))/det
      b(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))/det
      b(3,3) =  (a(1,1)*a(2,2)-a(1,2)*a(2,1))/det
    END SUBROUTINE invert3

  END SUBROUTINE EL_MM_ELEMENT_PHI_OPERATORS

  ! Assembled (rank-local, like every element assembly: the caller sums across ranks) nodal
  ! vectors of the two operators over the elements of this rank.
  !   coor(3,ndof) Q2 dof coordinates, eps_q2(ndof) transport porosity, pressure(4*nel) the P1
  !   pressure dofs, val_u/v/w nodal velocity; fp(3,ndof), fv(3,ndof).
  SUBROUTINE EL_MM_PHI_OPERATORS(mesh, ilev, coor, eps_q2, pressure, val_u, val_v, val_w, &
                                 with_viscous, fp, fv)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: coor(:,:), eps_q2(:), pressure(:), val_u(:), val_v(:), val_w(:)
    LOGICAL, INTENT(IN) :: with_viscous
    REAL*8, INTENT(OUT) :: fp(:,:), fv(:,:)

    INTEGER :: iel, l, dofs(NQ2)
    REAL*8 :: x(3,NQ2), phi(NQ2), u(3,NQ2), efp(3,NQ2), efv(3,NQ2)

    IF (.NOT.tables%ready) CALL EL_VA_FLUX_TABLES(tables)
    fp = 0.0d0
    fv = 0.0d0
    DO iel = 1, mesh%level(ilev)%nel
      CALL EL_VA_ELEMENT_DOFS(mesh, ilev, iel, dofs)
      DO l = 1, NQ2
        x(:,l) = coor(:,dofs(l))
        phi(l) = 1.0d0 - eps_q2(dofs(l))
        u(1,l) = val_u(dofs(l))
        u(2,l) = val_v(dofs(l))
        u(3,l) = val_w(dofs(l))
      END DO
      CALL EL_MM_ELEMENT_PHI_OPERATORS(tables, x, phi, pressure(4*(iel-1)+1:4*iel), u, &
        with_viscous, efp, efv)
      DO l = 1, NQ2
        fp(:,dofs(l)) = fp(:,dofs(l)) + efp(:,l)
        fv(:,dofs(l)) = fv(:,dofs(l)) + efv(:,l)
      END DO
    END DO

  END SUBROUTINE EL_MM_PHI_OPERATORS

END MODULE EL_MM_OPERATOR
