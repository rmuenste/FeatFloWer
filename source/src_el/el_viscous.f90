!=========================================================================
! EL_VISCOUS - kernel average of the viscous force density at a particle
!
! Momentum models Bfull and A give the particle the viscous force F_v = V div tau of the
! undisturbed (volume-averaged) fluid stress, CFDEM's viscForce. The fluid operator of this
! solver on the Newtonian path is the component Laplacian, tau = c grad u with the coefficient
! c the diffusion matrix carries (EL_FORCES: EL_FLUID_OPERATOR_VISCOSITY), so the force is
!     F_v = V c <lap u>_W ,
! with the kernel average taken in the weak form (no second derivatives of the Q2 field, the
! gradient jumps between elements are included automatically):
!     <lap u_c>_W = ( - int_Omega grad W . grad u_c dV + int_dOmega W du_c/dn dS ) / int_Omega W dV.
! The boundary integral runs over the element faces on a physical Dirichlet boundary inside the
! kernel support (the kernel is cut by the wall there; without the term the average would be
! biased for wall-adjacent particles). Subdomain interfaces are not boundaries: every rank
! integrates its own elements and the owner sums the contributions (EL_HALO), so the volume
! term is additive and interior faces cancel in the continuous sense.
!
! Geometry: isoparametric Q2 map and the Gauss rules of EL_VA_FLUX (3 x 3 x 3 volume, 3 x 3 per
! face), the same element description as the conservative continuity fluxes.
!
! Discrete consistency. With the kernel as wide as a few elements the Gauss rule does not
! integrate grad W exactly, so the plain weak form would return a spurious force for a uniform
! velocity gradient. Since - int grad W dV + int_dOmega W n dS = 0 exactly, the same quantity
! evaluated numerically (the defect B) is subtracted with the kernel-averaged gradient G:
!     <lap u_c>_W = ( A_c - B . G_c ) / N ,
!     A_c = - int grad W . grad u_c + int_wall W grad u_c . n      B = - int grad W + int_wall W n
!     G_c = int W grad u_c / N                                     N = int W
! which is exactly zero for a linear velocity field on any mesh and unchanged in the continuum.
!
! EL_INTEGRATE_PARTICLE_VISCOUS returns the additive pieces (A, B, int W grad u, N; EL_VIS_SIZE
! numbers) for the elements of one rank; they are summed at the owner and combined by
! EL_VISCOUS_LAPLACIAN.
!=========================================================================
MODULE EL_VISCOUS

  USE TYPES, ONLY: tMultiMesh
  USE EL_KERNEL_FUNCTIONS, ONLY: EL_KERNEL_VALUE, EL_KERNEL_DERIVATIVE
  USE EL_VA_FLUX, ONLY: tElVaFluxTables, EL_VA_FLUX_TABLES, EL_VA_ELEMENT_DOFS, &
                        EL_VA_FACE_GEOMETRY, NQ2, NFACE, NGPF, NGPV
  USE EL_QUADRATURE, ONLY: EL_ELEMENT_INTERSECTS_SUPPORT

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: EL_INTEGRATE_PARTICLE_VISCOUS, EL_VISCOUS_ELEMENT, EL_VISCOUS_LAPLACIAN, EL_VIS_SIZE

  ! layout of the additive record: A(1:3), B(4:6), int W grad u (7:15, row c at 7+3(c-1)), N (16)
  INTEGER, PARAMETER :: EL_VIS_SIZE = 16

  TYPE(tElVaFluxTables), SAVE :: tables
  INTEGER, SAVE :: face_dof(NFACE) = 0     ! local Q2 dof at the centre of each reference face

CONTAINS

  SUBROUTINE EL_VISCOUS_SETUP()
    INTEGER :: f
    IF (tables%ready) RETURN
    CALL EL_VA_FLUX_TABLES(tables)
    ! the face-centre dof is the one whose basis function is 1 at the face centre, i.e. at
    ! the middle face point (point 5 of the 3 x 3 rule sits at the centre)
    DO f = 1, NFACE
      face_dof(f) = MAXLOC(tables%fval(:,5,f), DIM=1)
    END DO
  END SUBROUTINE EL_VISCOUS_SETUP

  ! Contribution of one element with dof coordinates x and nodal velocities u to the numerator
  ! of <lap u>_W for a kernel centred at `position` with support radius `width`; wall(f) marks
  ! the faces on a physical Dirichlet boundary.
  PURE SUBROUTINE EL_VISCOUS_ELEMENT(t, x, u, wall, position, width, vis)

    TYPE(tElVaFluxTables), INTENT(IN) :: t
    REAL*8, INTENT(IN) :: x(3,NQ2), u(3,NQ2), position(3), width
    LOGICAL, INTENT(IN) :: wall(NFACE)
    REAL*8, INTENT(OUT) :: vis(EL_VIS_SIZE)

    REAL*8 :: jac(3,3), inv(3,3), detj, xq(3), gu(3,3), rvec(3), rlen, gw(3), w, nds(3), dv
    INTEGER :: p, l, k, f

    vis = 0.0d0
    ! volume terms
    DO p = 1, NGPV
      jac = 0.0d0
      xq = 0.0d0
      DO l = 1, NQ2
        DO k = 1, 3
          jac(:,k) = jac(:,k) + x(:,l)*t%vder(k,l,p)
        END DO
        xq = xq + x(:,l)*t%vval(l,p)
      END DO
      rvec = xq - position
      rlen = SQRT(SUM(rvec*rvec))
      IF (rlen.GE.width) CYCLE
      CALL invert3(jac, inv, detj)
      dv = t%vw(p)*ABS(detj)
      w = EL_KERNEL_VALUE(rlen, width)
      gw = 0.0d0
      IF (rlen.GT.0.0d0) gw = EL_KERNEL_DERIVATIVE(rlen, width)*rvec/rlen
      CALL velocity_gradient(u, t%vder(:,:,p), inv, gu)
      DO k = 1, 3
        vis(k) = vis(k) - DOT_PRODUCT(gu(k,:), gw)*dv
        vis(7+3*(k-1):9+3*(k-1)) = vis(7+3*(k-1):9+3*(k-1)) + w*gu(k,:)*dv
      END DO
      vis(4:6) = vis(4:6) - gw*dv
      vis(16) = vis(16) + w*dv
    END DO
    ! boundary terms on the wall faces
    DO f = 1, NFACE
      IF (.NOT.wall(f)) CYCLE
      DO p = 1, NGPF
        xq = 0.0d0
        DO l = 1, NQ2
          xq = xq + x(:,l)*t%fval(l,p,f)
        END DO
        rlen = SQRT(SUM((xq - position)**2))
        w = EL_KERNEL_VALUE(rlen, width)
        IF (w.LE.0.0d0) CYCLE
        CALL EL_VA_FACE_GEOMETRY(t, f, p, x, jac, detj, nds)
        CALL invert3(jac, inv, detj)
        CALL velocity_gradient(u, t%fder(:,:,p,f), inv, gu)
        DO k = 1, 3
          vis(k) = vis(k) + w*DOT_PRODUCT(gu(k,:), nds)
        END DO
        vis(4:6) = vis(4:6) + w*nds
      END DO
    END DO

  CONTAINS

    ! gu(c,j) = d u_c / d x_j from the reference derivatives dref(k,l) and inv = (dx/dxi)^-1
    PURE SUBROUTINE velocity_gradient(u, dref, inv, gu)
      REAL*8, INTENT(IN) :: u(3,NQ2), dref(3,NQ2), inv(3,3)
      REAL*8, INTENT(OUT) :: gu(3,3)
      REAL*8 :: dphys(3)
      INTEGER :: l, j
      gu = 0.0d0
      DO l = 1, NQ2
        ! d N_l / d x_j = sum_k dref(k,l) * dxi_k/dx_j = sum_k dref(k,l) inv(k,j)
        DO j = 1, 3
          dphys(j) = dref(1,l)*inv(1,j) + dref(2,l)*inv(2,j) + dref(3,l)*inv(3,j)
        END DO
        DO j = 1, 3
          gu(:,j) = gu(:,j) + u(:,l)*dphys(j)
        END DO
      END DO
    END SUBROUTINE velocity_gradient

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

  END SUBROUTINE EL_VISCOUS_ELEMENT

  ! Additive pieces of <lap u>_W over the elements of this rank. coor = Q2 dof coordinates of the
  ! level (mesh%level(ilev)%dcorvg), wall_dof = Dirichlet-boundary mask of the velocity dofs.
  SUBROUTINE EL_INTEGRATE_PARTICLE_VISCOUS(mesh, ilev, coor, val_u, val_v, val_w, wall_dof, &
                                           position, width, vis)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: coor(:,:), val_u(:), val_v(:), val_w(:), position(3), width
    LOGICAL, INTENT(IN) :: wall_dof(:)
    REAL*8, INTENT(OUT) :: vis(EL_VIS_SIZE)

    INTEGER :: iel, l, f, dofs(NQ2)
    REAL*8 :: x(3,NQ2), u(3,NQ2), part(EL_VIS_SIZE)
    LOGICAL :: wall(NFACE)

    vis = 0.0d0
    IF (width.LE.0.0d0) RETURN
    CALL EL_VISCOUS_SETUP()
    DO iel = 1, mesh%level(ilev)%nel
      IF (.NOT.EL_ELEMENT_INTERSECTS_SUPPORT(mesh, ilev, iel, position, width)) CYCLE
      CALL EL_VA_ELEMENT_DOFS(mesh, ilev, iel, dofs)
      DO l = 1, NQ2
        x(:,l) = coor(:,dofs(l))
        u(1,l) = val_u(dofs(l))
        u(2,l) = val_v(dofs(l))
        u(3,l) = val_w(dofs(l))
      END DO
      DO f = 1, NFACE
        wall(f) = wall_dof(dofs(face_dof(f)))
      END DO
      CALL EL_VISCOUS_ELEMENT(tables, x, u, wall, position, width, part)
      vis = vis + part
    END DO

  END SUBROUTINE EL_INTEGRATE_PARTICLE_VISCOUS

  ! <lap u_c>_W from the summed record (see the module header); zero without kernel support.
  PURE FUNCTION EL_VISCOUS_LAPLACIAN(vis) RESULT(lap)
    REAL*8, INTENT(IN) :: vis(EL_VIS_SIZE)
    REAL*8 :: lap(3)
    INTEGER :: c
    lap = 0.0d0
    IF (vis(16).LE.0.0d0) RETURN
    DO c = 1, 3
      lap(c) = (vis(c) - DOT_PRODUCT(vis(4:6), vis(7+3*(c-1):9+3*(c-1)))/vis(16))/vis(16)
    END DO
  END FUNCTION EL_VISCOUS_LAPLACIAN

END MODULE EL_VISCOUS
