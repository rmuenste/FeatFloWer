!=========================================================================
! CHI_FORCES - surface-stress force/torque integration on a submesh
! boundary, Layer M (design: chimera-integration-design.md v3, section 5).
!
! Normative convention (paper eqs. (4a)/(4b) with the v1-review
! correction): with n_out the OUTWARD normal of the atmosphere domain
! (which on the inner surface points INTO the particle, n_out = -n_B),
!
!   F = - surfint sigma . n_out dS
!   T = - surfint (x - X) x (sigma . n_out) dS
!
! which equals the classical + surfint sigma . n_B over the particle
! surface.  sigma = -p I + mu (grad u + grad u^T); the pressure is the
! discontinuous P1 of the adjacent element (centroid-linear), gradients
! are Q2.
!=========================================================================
MODULE CHI_FORCES

  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP, CHI_M33INV
  USE CHI_FEM_EVAL, ONLY: CHI_Q2_BASIS, CHI_Q2_DOFMAP, CHI_EVAL_P1
  USE CHI_KERNELS, ONLY: CHI_FACE_RULE, CHI_FACE_GEOM

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CHI_COMPUTE_FORCES

CONTAINS

  !-----------------------------------------------------------------------
  ! Integrate force and torque over a boundary-face list (typically the
  ! inner surface).  xref: torque reference point (particle center).
  ! u/v/w: Q2 velocity dofs; p: 4*nel discontinuous P1 pressure dofs.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_COMPUTE_FORCES(faces, nfaces, nvt, net, nat, nel, &
                                kvert, kedge, karea, dcorvg, &
                                u, v, w, p, mu, xref, F, T)
    INTEGER, INTENT(IN) :: faces(2,*), nfaces, nvt, net, nat, nel
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*)
    REAL*8,  INTENT(IN) :: dcorvg(3,*)
    REAL*8,  INTENT(IN) :: u(*), v(*), w(*), p(*)
    REAL*8,  INTENT(IN) :: mu, xref(3)
    REAL*8,  INTENT(OUT) :: F(3), T(3)

    REAL*8 :: xiq(3,9), wq(9), phi(27), dphi(3,27)
    REAL*8 :: nodes(3,8), jac(3,3), jacinv(3,3), detj
    REAL*8 :: x(3), nrm(3), sj, dSw, xc(3)
    REAL*8 :: gphi(3,27), gradU(3,3), sigma(3,3), tr(3), pq, rel(3)
    INTEGER :: idx(27), ifc, e, fl, q, i, j, a, b
    LOGICAL :: ok

    F = 0d0
    T = 0d0

    DO ifc = 1, nfaces
      e = faces(1,ifc)
      fl = faces(2,ifc)
      DO i = 1, 8
        nodes(:,i) = dcorvg(:,kvert(i,e))
      END DO
      CALL CHI_Q2_DOFMAP(e, kvert, kedge, karea, nvt, net, nat, idx)
      xc = 0.125d0*SUM(nodes, DIM=2)
      CALL CHI_FACE_RULE(fl, xiq, wq)

      DO q = 1, 9
        CALL CHI_Q2_BASIS(xiq(:,q), phi, dphi)
        CALL CHI_FACE_GEOM(nodes, fl, xiq(:,q), x, nrm, sj)
        dSw = wq(q)*sj
        CALL CHI_Q1_MAP(nodes, xiq(:,q), x, jac, detj)
        CALL CHI_M33INV(jac, jacinv, ok)
        IF (.NOT. ok) THEN
          WRITE(*,*) 'CHI_COMPUTE_FORCES: singular Jacobian, el', e
          STOP 1
        END IF
        DO j = 1, 27
          gphi(:,j) = MATMUL(TRANSPOSE(jacinv), dphi(:,j))
        END DO
        ! gradU(a,:) = grad of velocity component a
        gradU = 0d0
        DO j = 1, 27
          gradU(1,:) = gradU(1,:) + u(idx(j))*gphi(:,j)
          gradU(2,:) = gradU(2,:) + v(idx(j))*gphi(:,j)
          gradU(3,:) = gradU(3,:) + w(idx(j))*gphi(:,j)
        END DO
        pq = CHI_EVAL_P1(p(4*(e-1)+1:4*(e-1)+4), x, xc)
        DO a = 1, 3
          DO b = 1, 3
            sigma(a,b) = mu*(gradU(a,b) + gradU(b,a))
          END DO
          sigma(a,a) = sigma(a,a) - pq
        END DO
        tr = MATMUL(sigma, nrm)
        F = F - dSw*tr
        rel = x - xref
        T(1) = T(1) - dSw*(rel(2)*tr(3) - rel(3)*tr(2))
        T(2) = T(2) - dSw*(rel(3)*tr(1) - rel(1)*tr(3))
        T(3) = T(3) - dSw*(rel(1)*tr(2) - rel(2)*tr(1))
      END DO
    END DO
  END SUBROUTINE CHI_COMPUTE_FORCES

END MODULE CHI_FORCES
