!=========================================================================
! CHI_FEM_EVAL - mesh-agnostic Q2/P1 evaluation primitives of the
! Chimera component (design: chimera-integration-design.md v3, Layer M).
!
! The 27-node triquadratic basis reproduces EXACTLY the FeatFloWer Q2
! local ordering (verified against the hard-coded basis in RETURN_Velo,
! source/src_particles/part_step.f90:1711, and E013):
!   1..8   vertices        (FEAT kvert ordering, see CHI_GEOMETRY)
!   9..20  edge midpoints  (FEAT kedge ordering: 4 bottom, 4 vertical,
!                           4 top - edges (1,2)(2,3)(3,4)(4,1),
!                           (1,5)(2,6)(3,7)(4,8), (5,6)(6,7)(7,8)(8,5))
!   21..26 face centers    (FEAT karea ordering: bottom (1234),
!                           front (1265), right (2376), back (3487),
!                           left (4158), top (5678))
!   27     element center
! Global Q2 DOF numbering (SetUp_myQ2Coor convention):
!   1..nvt vertices, nvt+1..nvt+net edge DOFs, +nat face DOFs, +iel
!   element-center DOFs.
!
! P1 pressure is discontinuous linear in PHYSICAL coordinates about the
! element centroid: p(x) = p1 + p2*(x-xc) + p3*(y-yc) + p4*(z-zc)
! (convention of IntP1toQ2, source/src_quadLS/QuadSc_def.f90).
!
! No COMMON blocks, no module state.
!=========================================================================
MODULE CHI_FEM_EVAL

  USE CHI_GEOMETRY, ONLY: CHI_M33INV, CHI_Q1_MAP

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CHI_Q2_BASIS
  PUBLIC :: CHI_Q2_DOFMAP
  PUBLIC :: CHI_Q2_REFNODES
  PUBLIC :: CHI_EVAL_Q2_SCALAR
  PUBLIC :: CHI_EVAL_Q2_GRADIENT
  PUBLIC :: CHI_EVAL_P1
  PUBLIC :: CHI_EVAL_FIELD_AT

  ! 1D node index (-1, 0, +1) of each of the 27 local DOFs per direction,
  ! matching the FeatFloWer local ordering documented above.
  INTEGER, PARAMETER :: chi_q2_nodes(3,27) = RESHAPE((/ &
    -1,-1,-1,   1,-1,-1,   1, 1,-1,  -1, 1,-1, &   ! vertices 1-4 (bottom)
    -1,-1, 1,   1,-1, 1,   1, 1, 1,  -1, 1, 1, &   ! vertices 5-8 (top)
     0,-1,-1,   1, 0,-1,   0, 1,-1,  -1, 0,-1, &   ! edges 9-12 (bottom)
    -1,-1, 0,   1,-1, 0,   1, 1, 0,  -1, 1, 0, &   ! edges 13-16 (vertical)
     0,-1, 1,   1, 0, 1,   0, 1, 1,  -1, 0, 1, &   ! edges 17-20 (top)
     0, 0,-1,   0,-1, 0,   1, 0, 0,   0, 1, 0, &   ! faces 21-24
    -1, 0, 0,   0, 0, 1,                        &   ! faces 25-26
     0, 0, 0 /), (/3,27/))                          ! center 27

CONTAINS

  !-----------------------------------------------------------------------
  ! Q2 basis values and reference derivatives at an arbitrary reference
  ! point xi in [-1,1]^3.  Tensor product of the 1D quadratic Lagrange
  ! functions on nodes {-1, 0, +1}:
  !   L_{-1}(t) = t(t-1)/2,  L_0(t) = 1-t^2,  L_{+1}(t) = t(t+1)/2
  ! which reproduces the RETURN_Velo/E013 basis exactly.
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_Q2_BASIS(xi, phi, dphi)
    REAL*8, INTENT(IN)  :: xi(3)
    REAL*8, INTENT(OUT) :: phi(27), dphi(3,27)

    REAL*8 :: l(3,-1:1), dl(3,-1:1)
    INTEGER :: d, i, a, b, c

    DO d = 1, 3
      l(d,-1) = 0.5d0*xi(d)*(xi(d)-1d0)
      l(d, 0) = 1d0 - xi(d)*xi(d)
      l(d, 1) = 0.5d0*xi(d)*(xi(d)+1d0)
      dl(d,-1) = xi(d) - 0.5d0
      dl(d, 0) = -2d0*xi(d)
      dl(d, 1) = xi(d) + 0.5d0
    END DO

    DO i = 1, 27
      a = chi_q2_nodes(1,i)
      b = chi_q2_nodes(2,i)
      c = chi_q2_nodes(3,i)
      phi(i)    =  l(1,a)* l(2,b)* l(3,c)
      dphi(1,i) = dl(1,a)* l(2,b)* l(3,c)
      dphi(2,i) =  l(1,a)*dl(2,b)* l(3,c)
      dphi(3,i) =  l(1,a)* l(2,b)*dl(3,c)
    END DO
  END SUBROUTINE CHI_Q2_BASIS

  !-----------------------------------------------------------------------
  ! Reference coordinates of the 27 local DOFs (for tests and node-
  ! coordinate generation).
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_Q2_REFNODES(refxi)
    REAL*8, INTENT(OUT) :: refxi(3,27)
    INTEGER :: i
    DO i = 1, 27
      refxi(:,i) = DBLE(chi_q2_nodes(:,i))
    END DO
  END SUBROUTINE CHI_Q2_REFNODES

  !-----------------------------------------------------------------------
  ! Local (1..27) -> global Q2 DOF indices of element iel.
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_Q2_DOFMAP(iel, kvert, kedge, karea, nvt, net, nat, idx)
    INTEGER, INTENT(IN)  :: iel, nvt, net, nat
    INTEGER, INTENT(IN)  :: kvert(8,*), kedge(12,*), karea(6,*)
    INTEGER, INTENT(OUT) :: idx(27)

    INTEGER :: i

    DO i = 1, 8
      idx(i) = kvert(i,iel)
    END DO
    DO i = 1, 12
      idx(8+i) = nvt + kedge(i,iel)
    END DO
    DO i = 1, 6
      idx(20+i) = nvt + net + karea(i,iel)
    END DO
    idx(27) = nvt + net + nat + iel
  END SUBROUTINE CHI_Q2_DOFMAP

  !-----------------------------------------------------------------------
  ! Value of a Q2 field from its 27 element-local nodal values.
  !-----------------------------------------------------------------------
  PURE FUNCTION CHI_EVAL_Q2_SCALAR(vals, phi) RESULT(val)
    REAL*8, INTENT(IN) :: vals(27), phi(27)
    REAL*8 :: val
    val = SUM(vals*phi)
  END FUNCTION CHI_EVAL_Q2_SCALAR

  !-----------------------------------------------------------------------
  ! Physical gradient of a Q2 field: grad_x = J^{-T} grad_xi, with the
  ! Q1 geometry Jacobian jac from CHI_Q1_MAP (FeatFloWer hexes carry Q1
  ! geometry; cf. the EL_Q1_MAP commentary in el_quadrature.f90).
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_EVAL_Q2_GRADIENT(vals, dphi, jac, grad, ok)
    REAL*8, INTENT(IN)  :: vals(27), dphi(3,27), jac(3,3)
    REAL*8, INTENT(OUT) :: grad(3)
    LOGICAL, INTENT(OUT) :: ok

    REAL*8 :: jacinv(3,3), gref(3)
    INTEGER :: d

    DO d = 1, 3
      gref(d) = SUM(vals*dphi(d,:))
    END DO
    CALL CHI_M33INV(jac, jacinv, ok)
    IF (.NOT. ok) THEN
      grad = 0d0
      RETURN
    END IF
    ! grad_x = J^{-T} * gref
    grad = MATMUL(TRANSPOSE(jacinv), gref)
  END SUBROUTINE CHI_EVAL_Q2_GRADIENT

  !-----------------------------------------------------------------------
  ! Discontinuous P1 pressure at physical point x, from the 4 element
  ! DOFs and the element centroid xc.
  !-----------------------------------------------------------------------
  PURE FUNCTION CHI_EVAL_P1(pdofs, x, xc) RESULT(p)
    REAL*8, INTENT(IN) :: pdofs(4), x(3), xc(3)
    REAL*8 :: p
    p = pdofs(1) + pdofs(2)*(x(1)-xc(1)) + pdofs(3)*(x(2)-xc(2)) &
                 + pdofs(4)*(x(3)-xc(3))
  END FUNCTION CHI_EVAL_P1

  !-----------------------------------------------------------------------
  ! Velocity, velocity gradient and P1 pressure of a Q2/P1 field at the
  ! reference point xi of element iel (Phase 3: background evaluation for
  ! the Robin data and submesh evaluation for fringe values).
  ! gradu(a,b) = d u_a / d x_b.  p holds 4 dofs per element (FeatFloWer
  ! P1 layout: value at the centroid + 3 slopes).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_EVAL_FIELD_AT(iel, xi, kvert, kedge, karea, nvt, net, nat, &
                               dcorvg, u, v, w, p, uval, gradu, pval, ok)
    INTEGER, INTENT(IN) :: iel, kvert(8,*), kedge(12,*), karea(6,*)
    INTEGER, INTENT(IN) :: nvt, net, nat
    REAL*8,  INTENT(IN) :: xi(3), dcorvg(3,*), u(*), v(*), w(*), p(*)
    REAL*8,  INTENT(OUT) :: uval(3), gradu(3,3), pval
    LOGICAL, INTENT(OUT) :: ok

    REAL*8 :: nodes(3,8), phi(27), dphi(3,27), jac(3,3), detj
    REAL*8 :: x(3), xc(3), vals(27), g(3)
    INTEGER :: idx(27), i, a
    LOGICAL :: gok

    ok = .FALSE.
    uval = 0d0
    gradu = 0d0
    pval = 0d0
    DO i = 1, 8
      nodes(:,i) = dcorvg(:,kvert(i,iel))
    END DO
    CALL CHI_Q2_DOFMAP(iel, kvert, kedge, karea, nvt, net, nat, idx)
    CALL CHI_Q2_BASIS(xi, phi, dphi)
    CALL CHI_Q1_MAP(nodes, xi, x, jac, detj)
    DO a = 1, 3
      SELECT CASE (a)
      CASE (1)
        vals = u(idx)
      CASE (2)
        vals = v(idx)
      CASE DEFAULT
        vals = w(idx)
      END SELECT
      uval(a) = CHI_EVAL_Q2_SCALAR(vals, phi)
      CALL CHI_EVAL_Q2_GRADIENT(vals, dphi, jac, g, gok)
      IF (.NOT. gok) RETURN
      gradu(a,:) = g
    END DO
    CALL CHI_Q1_MAP(nodes, (/0d0,0d0,0d0/), xc, jac, detj)
    pval = CHI_EVAL_P1(p(4*(iel-1)+1:4*iel), x, xc)
    ok = .TRUE.
  END SUBROUTINE CHI_EVAL_FIELD_AT

END MODULE CHI_FEM_EVAL
