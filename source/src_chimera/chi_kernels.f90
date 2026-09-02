!=========================================================================
! CHI_KERNELS - reentrant Chimera assembly kernels, Layer M
! (design: chimera-integration-design.md v3, section 3; replaces any use
! of the legacy F77 assembly kernels for submesh work - those read
! /ELEM/ /CUB/ /COAUX1/ /TRIAD/ and var_QuadScalar and are therefore not
! reentrant; v2-review blocker 1).
!
! Everything is passed explicitly: connectivity, coordinates, physical
! parameters, Picard velocity, and the target CSR. No COMMON blocks, no
! module state, no USE of solver modules.
!
! Discretization: Q2 velocity / discontinuous P1 pressure (physical
! centroid-linear, IntP1toQ2 convention), Q1 (trilinear) geometry,
! 3x3x3 Gauss volume rule, 3x3 Gauss face rule.
!
! Weak form (test function v = phi_I e_a, trial u = phi_J e_b):
!   momentum block (a,b):
!     dtinv*rho*(phi_J, phi_I) delta_ab            [time term, 0 = steady]
!     + rho*(((u_k - u_mesh) . grad phi_J), phi_I) delta_ab   [Picard]
!     + mu*( delta_ab grad phi_J . grad phi_I + d_a phi_J d_b phi_I )
!                                    [deformation form, 2 mu D(u):D(v)]
!   momentum-pressure:  -(psi_kp, d_a phi_I)       [-(p, div v)]
!   continuity:         +(psi_ip, d_b phi_J)       [(q, div u)]
!
! Robin outer boundary (design section 5, paper eq. (5d), n = n_Gamma =
! outward normal of the atmosphere):
!   LHS += - alpha (u_k . n) (phi_J phi_I) delta_ab   [Picard-linearized]
!   RHS += + (h, phi_I),  h = sigma(u_bg,p_bg) n - alpha (u_bg.n) u_bg
!
! Global unknown layout (block order): u = 1..ndof, v = ndof+1..2*ndof,
! w = 2*ndof+1..3*ndof, pressure = 3*ndof + 4*(iel-1) + (1..4).
!=========================================================================
MODULE CHI_KERNELS

  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP, CHI_GAUSS3, CHI_M33INV
  USE CHI_FEM_EVAL, ONLY: CHI_Q2_BASIS, CHI_Q2_DOFMAP

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CHI_BUILD_SADDLE_CSR
  PUBLIC :: CHI_ASM_SADDLE
  PUBLIC :: CHI_ASM_ROBIN
  PUBLIC :: CHI_APPLY_DIRICHLET_ROW
  PUBLIC :: CHI_FACE_RULE
  PUBLIC :: CHI_FACE_GEOM
  PUBLIC :: chi_bc_velocity, chi_robin_data

  ! Reference-face parametrization, consistent with the FEAT karea
  ! ordering (bottom, front, right, back, left, top).  For face f the
  ! coordinate chi_face_fixed_axis(f) is fixed at chi_face_fixed_val(f);
  ! the free axes are (axisA, axisB) and the OUTWARD normal of the
  ! element is chi_face_sign(f) * (dx/dA x dx/dB).
  INTEGER, PARAMETER :: chi_face_fixed_axis(6) = (/ 3, 2, 1, 2, 1, 3 /)
  REAL*8,  PARAMETER :: chi_face_fixed_val(6)  = &
    (/ -1d0, -1d0, 1d0, 1d0, -1d0, 1d0 /)
  INTEGER, PARAMETER :: chi_face_axisA(6) = (/ 1, 1, 2, 1, 2, 1 /)
  INTEGER, PARAMETER :: chi_face_axisB(6) = (/ 2, 3, 3, 3, 3, 2 /)
  REAL*8,  PARAMETER :: chi_face_sign(6)  = &
    (/ -1d0, 1d0, 1d0, -1d0, -1d0, 1d0 /)

  ABSTRACT INTERFACE
    ! Dirichlet velocity value at a (projected) boundary Q2 node.
    SUBROUTINE chi_bc_velocity(x, u)
      REAL*8, INTENT(IN)  :: x(3)
      REAL*8, INTENT(OUT) :: u(3)
    END SUBROUTINE chi_bc_velocity
    ! Robin data h at an outer-boundary quadrature point with outward
    ! normal n (design section 5).
    SUBROUTINE chi_robin_data(x, n, h)
      REAL*8, INTENT(IN)  :: x(3), n(3)
      REAL*8, INTENT(OUT) :: h(3)
    END SUBROUTINE chi_robin_data
  END INTERFACE

CONTAINS

  !-----------------------------------------------------------------------
  ! Build the sparsity pattern of the monolithic saddle-point matrix
  ! (n = 3*ndof + 4*nel).  Columns are emitted in ascending order per
  ! row, which the assembly scatter (binary search) relies on.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_BUILD_SADDLE_CSR(nel, nvt, net, nat, kvert, kedge, karea, &
                                  n, LdA, ColA)
    INTEGER, INTENT(IN) :: nel, nvt, net, nat
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*)
    INTEGER, INTENT(OUT) :: n
    INTEGER, ALLOCATABLE, INTENT(OUT) :: LdA(:), ColA(:)

    INTEGER :: ndof, e, i, j, k, b, r, pos, nnz, deg
    INTEGER :: idx(27), row27(27)
    INTEGER, ALLOCATABLE :: eldofs(:,:)
    INTEGER, ALLOCATABLE :: d2ePtr(:), d2e(:), cnt(:)
    INTEGER, ALLOCATABLE :: adjPtr(:), adj(:), stamp(:)
    INTEGER :: adjnnz, je, jj, col

    ndof = nvt + net + nat + nel
    n = 3*ndof + 4*nel

    ALLOCATE(eldofs(27,nel))
    DO e = 1, nel
      CALL CHI_Q2_DOFMAP(e, kvert, kedge, karea, nvt, net, nat, idx)
      eldofs(:,e) = idx
    END DO

    ! --- dof -> element incidence (CSR) --------------------------------
    ALLOCATE(cnt(ndof))
    cnt = 0
    DO e = 1, nel
      DO k = 1, 27
        cnt(eldofs(k,e)) = cnt(eldofs(k,e)) + 1
      END DO
    END DO
    ALLOCATE(d2ePtr(ndof+1))
    d2ePtr(1) = 1
    DO i = 1, ndof
      d2ePtr(i+1) = d2ePtr(i) + cnt(i)
    END DO
    ALLOCATE(d2e(d2ePtr(ndof+1)-1))
    cnt = 0
    DO e = 1, nel          ! ascending e => element lists sorted
      DO k = 1, 27
        i = eldofs(k,e)
        d2e(d2ePtr(i) + cnt(i)) = e
        cnt(i) = cnt(i) + 1
      END DO
    END DO
    ! dedupe repeated element entries cannot occur (27 distinct dofs/el)

    ! --- scalar Q2 adjacency (two-pass, marker dedupe, sorted rows) ----
    ALLOCATE(stamp(ndof))
    stamp = 0
    ALLOCATE(adjPtr(ndof+1))
    adjPtr(1) = 1
    DO i = 1, ndof
      deg = 0
      DO je = d2ePtr(i), d2ePtr(i+1)-1
        e = d2e(je)
        DO k = 1, 27
          j = eldofs(k,e)
          IF (stamp(j) .NE. i) THEN
            stamp(j) = i
            deg = deg + 1
          END IF
        END DO
      END DO
      adjPtr(i+1) = adjPtr(i) + deg
    END DO
    ALLOCATE(adj(adjPtr(ndof+1)-1))
    stamp = 0
    DO i = 1, ndof
      deg = 0
      DO je = d2ePtr(i), d2ePtr(i+1)-1
        e = d2e(je)
        DO k = 1, 27
          j = eldofs(k,e)
          IF (stamp(j) .NE. i) THEN
            stamp(j) = i
            adj(adjPtr(i) + deg) = j
            deg = deg + 1
          END IF
        END DO
      END DO
      CALL sort_ascending(adj(adjPtr(i):adjPtr(i+1)-1), deg)
    END DO

    ! --- global CSR ----------------------------------------------------
    ALLOCATE(LdA(n+1))
    LdA(1) = 1
    DO b = 1, 3
      DO i = 1, ndof
        r = (b-1)*ndof + i
        LdA(r+1) = LdA(r) + 3*(adjPtr(i+1)-adjPtr(i)) &
                          + 4*(d2ePtr(i+1)-d2ePtr(i))
      END DO
    END DO
    DO e = 1, nel
      DO k = 1, 4
        r = 3*ndof + 4*(e-1) + k
        LdA(r+1) = LdA(r) + 3*27 + 4
      END DO
    END DO
    nnz = LdA(n+1) - 1
    ALLOCATE(ColA(nnz))

    DO b = 1, 3
      DO i = 1, ndof
        r = (b-1)*ndof + i
        pos = LdA(r)
        DO k = 1, 3
          DO jj = adjPtr(i), adjPtr(i+1)-1
            ColA(pos) = (k-1)*ndof + adj(jj)
            pos = pos + 1
          END DO
        END DO
        DO je = d2ePtr(i), d2ePtr(i+1)-1
          e = d2e(je)
          DO k = 1, 4
            ColA(pos) = 3*ndof + 4*(e-1) + k
            pos = pos + 1
          END DO
        END DO
      END DO
    END DO
    DO e = 1, nel
      row27 = eldofs(:,e)
      CALL sort_ascending(row27, 27)
      DO k = 1, 4
        r = 3*ndof + 4*(e-1) + k
        pos = LdA(r)
        DO b = 1, 3
          DO jj = 1, 27
            ColA(pos) = (b-1)*ndof + row27(jj)
            pos = pos + 1
          END DO
        END DO
        DO col = 1, 4
          ColA(pos) = 3*ndof + 4*(e-1) + col
          pos = pos + 1
        END DO
      END DO
    END DO

    DEALLOCATE(eldofs, d2ePtr, d2e, cnt, adjPtr, adj, stamp)
  END SUBROUTINE CHI_BUILD_SADDLE_CSR

  !-----------------------------------------------------------------------
  ! Assemble the volume terms of the saddle-point system into Avals
  ! (zeroed here).  uk/vk/wk: Picard velocity at the Q2 dofs; umesh:
  ! constant submesh (ALE) velocity; dtinv = 1/dt (0 = steady).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_ASM_SADDLE(nel, nvt, net, nat, kvert, kedge, karea, &
                            dcorvg, n, LdA, ColA, Avals, &
                            rho, mu, dtinv, uk, vk, wk, umesh)
    INTEGER, INTENT(IN) :: nel, nvt, net, nat
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*)
    REAL*8,  INTENT(IN) :: dcorvg(3,*)
    INTEGER, INTENT(IN) :: n, LdA(*), ColA(*)
    REAL*8,  INTENT(OUT) :: Avals(*)
    REAL*8,  INTENT(IN) :: rho, mu, dtinv
    REAL*8,  INTENT(IN) :: uk(*), vk(*), wk(*), umesh(3)

    REAL*8 :: gp(3,27), gw(27)
    REAL*8 :: phi(27), dphi(3,27), phig(27,27), dphig(3,27,27)
    REAL*8 :: nodes(3,8), jac(3,3), jacinv(3,3), detj, xq(3), xc(3)
    REAL*8 :: gphi(3,27), ukq(3), conv(27), psi(4), w
    REAL*8 :: Aloc(27,3,27,3), Gloc(27,3,4), Dloc(4,27,3)
    INTEGER :: idx(27), ndof, e, q, i, j, a, b, k, r, pos
    LOGICAL :: ok
    REAL*8 :: diag_ab, lap

    ndof = nvt + net + nat + nel
    Avals(1:LdA(n+1)-1) = 0d0

    CALL CHI_GAUSS3(gp, gw)
    DO q = 1, 27
      CALL CHI_Q2_BASIS(gp(:,q), phi, dphi)
      phig(:,q) = phi
      dphig(:,:,q) = dphi
    END DO

    DO e = 1, nel
      DO i = 1, 8
        nodes(:,i) = dcorvg(:,kvert(i,e))
      END DO
      CALL CHI_Q2_DOFMAP(e, kvert, kedge, karea, nvt, net, nat, idx)
      xc = 0.125d0*SUM(nodes, DIM=2)

      Aloc = 0d0
      Gloc = 0d0
      Dloc = 0d0

      DO q = 1, 27
        CALL CHI_Q1_MAP(nodes, gp(:,q), xq, jac, detj)
        CALL CHI_M33INV(jac, jacinv, ok)
        IF (.NOT. ok) THEN
          WRITE(*,*) 'CHI_ASM_SADDLE: singular geometry Jacobian, el', e
          STOP 1
        END IF
        w = gw(q)*ABS(detj)
        ! physical gradients: grad_x = J^{-T} grad_xi
        DO j = 1, 27
          gphi(:,j) = MATMUL(TRANSPOSE(jacinv), dphig(:,j,q))
        END DO
        ! Picard transport velocity at the point
        ukq = 0d0
        DO j = 1, 27
          ukq(1) = ukq(1) + phig(j,q)*uk(idx(j))
          ukq(2) = ukq(2) + phig(j,q)*vk(idx(j))
          ukq(3) = ukq(3) + phig(j,q)*wk(idx(j))
        END DO
        ukq = ukq - umesh
        DO j = 1, 27
          conv(j) = rho*(ukq(1)*gphi(1,j) + ukq(2)*gphi(2,j) + &
                         ukq(3)*gphi(3,j))
        END DO
        psi(1) = 1d0
        psi(2) = xq(1) - xc(1)
        psi(3) = xq(2) - xc(2)
        psi(4) = xq(3) - xc(3)

        DO j = 1, 27
          DO i = 1, 27
            lap = gphi(1,i)*gphi(1,j) + gphi(2,i)*gphi(2,j) + &
                  gphi(3,i)*gphi(3,j)
            diag_ab = dtinv*rho*phig(i,q)*phig(j,q) &
                    + phig(i,q)*conv(j) + mu*lap
            DO a = 1, 3
              Aloc(i,a,j,a) = Aloc(i,a,j,a) + w*diag_ab
              DO b = 1, 3
                Aloc(i,a,j,b) = Aloc(i,a,j,b) + w*mu*gphi(a,j)*gphi(b,i)
              END DO
            END DO
          END DO
        END DO
        DO k = 1, 4
          DO i = 1, 27
            DO a = 1, 3
              Gloc(i,a,k) = Gloc(i,a,k) - w*psi(k)*gphi(a,i)
              Dloc(k,i,a) = Dloc(k,i,a) + w*psi(k)*gphi(a,i)
            END DO
          END DO
        END DO
      END DO

      ! --- scatter -----------------------------------------------------
      DO a = 1, 3
        DO i = 1, 27
          r = (a-1)*ndof + idx(i)
          DO b = 1, 3
            DO j = 1, 27
              pos = csr_pos(LdA, ColA, r, (b-1)*ndof + idx(j))
              Avals(pos) = Avals(pos) + Aloc(i,a,j,b)
            END DO
          END DO
          DO k = 1, 4
            pos = csr_pos(LdA, ColA, r, 3*ndof + 4*(e-1) + k)
            Avals(pos) = Avals(pos) + Gloc(i,a,k)
          END DO
        END DO
      END DO
      DO k = 1, 4
        r = 3*ndof + 4*(e-1) + k
        DO a = 1, 3
          DO j = 1, 27
            pos = csr_pos(LdA, ColA, r, (a-1)*ndof + idx(j))
            Avals(pos) = Avals(pos) + Dloc(k,j,a)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE CHI_ASM_SADDLE

  !-----------------------------------------------------------------------
  ! 3x3 Gauss rule on a reference face: returns the 9 volumetric
  ! reference coordinates and 2D weights for face lface.
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_FACE_RULE(lface, xiq, wq)
    INTEGER, INTENT(IN)  :: lface
    REAL*8,  INTENT(OUT) :: xiq(3,9), wq(9)

    REAL*8, PARAMETER :: g = 0.774596669241483377035853079956d0
    REAL*8 :: p1(3), w1(3)
    INTEGER :: i, j, q

    p1 = (/ -g, 0d0, g /)
    w1 = (/ 5d0/9d0, 8d0/9d0, 5d0/9d0 /)
    q = 0
    DO j = 1, 3
      DO i = 1, 3
        q = q + 1
        xiq(chi_face_fixed_axis(lface), q) = chi_face_fixed_val(lface)
        xiq(chi_face_axisA(lface), q) = p1(i)
        xiq(chi_face_axisB(lface), q) = p1(j)
        wq(q) = w1(i)*w1(j)
      END DO
    END DO
  END SUBROUTINE CHI_FACE_RULE

  !-----------------------------------------------------------------------
  ! Face geometry at a reference point: physical point, OUTWARD unit
  ! normal of the element/domain, and the surface Jacobian |t_A x t_B|
  ! (multiply by the 2D weight for the surface measure).
  !
  ! Orientation: the reference sign table is only valid for positively
  ! oriented elements, but the FEAT 1:8 refinement produces children of
  ! MIXED orientation (negative detJ; harmless in the production code,
  ! which uses |detJ| and never derives normals from reference tables).
  ! The outward direction is therefore fixed geometrically: the normal
  ! must point from the element centroid toward the face point.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_FACE_GEOM(nodes, lface, xi, x, nrm, surfjac)
    REAL*8,  INTENT(IN)  :: nodes(3,8), xi(3)
    INTEGER, INTENT(IN)  :: lface
    REAL*8,  INTENT(OUT) :: x(3), nrm(3), surfjac

    REAL*8 :: jac(3,3), detj, tA(3), tB(3), cr(3), outv(3)

    CALL CHI_Q1_MAP(nodes, xi, x, jac, detj)
    tA = jac(:,chi_face_axisA(lface))
    tB = jac(:,chi_face_axisB(lface))
    cr(1) = tA(2)*tB(3) - tA(3)*tB(2)
    cr(2) = tA(3)*tB(1) - tA(1)*tB(3)
    cr(3) = tA(1)*tB(2) - tA(2)*tB(1)
    surfjac = SQRT(cr(1)*cr(1) + cr(2)*cr(2) + cr(3)*cr(3))
    IF (surfjac .LE. 0d0) THEN
      WRITE(*,*) 'CHI_FACE_GEOM: degenerate face'
      STOP 1
    END IF
    nrm = chi_face_sign(lface)*cr/surfjac
    ! geometric outward orientation (element centroid -> face point)
    outv = x - 0.125d0*SUM(nodes, DIM=2)
    IF (nrm(1)*outv(1) + nrm(2)*outv(2) + nrm(3)*outv(3) .LT. 0d0) THEN
      nrm = -nrm
    END IF
  END SUBROUTINE CHI_FACE_GEOM

  !-----------------------------------------------------------------------
  ! Robin boundary contribution on a list of outer faces (design
  ! section 5): matrix += -alpha (u_k.n) phi_J phi_I (block diagonal),
  ! rhs += (h, phi_I) with h from the callback.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_ASM_ROBIN(faces, nfaces, nvt, net, nat, nel, &
                           kvert, kedge, karea, dcorvg, &
                           n, LdA, ColA, Avals, rhs, &
                           alpha, uk, vk, wk, robin_h)
    INTEGER, INTENT(IN) :: faces(2,*), nfaces, nvt, net, nat, nel
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*)
    REAL*8,  INTENT(IN) :: dcorvg(3,*)
    INTEGER, INTENT(IN) :: n, LdA(*), ColA(*)
    REAL*8,  INTENT(INOUT) :: Avals(*), rhs(*)
    REAL*8,  INTENT(IN) :: alpha, uk(*), vk(*), wk(*)
    PROCEDURE(chi_robin_data) :: robin_h

    REAL*8 :: xiq(3,9), wq(9), phi(27), dphi(3,27)
    REAL*8 :: nodes(3,8), x(3), nrm(3), sj, dSw
    REAL*8 :: ukq(3), un, h(3)
    INTEGER :: idx(27), ndof, ifc, e, f, q, i, j, a, r, pos

    ndof = nvt + net + nat + nel

    DO ifc = 1, nfaces
      e = faces(1,ifc)
      f = faces(2,ifc)
      DO i = 1, 8
        nodes(:,i) = dcorvg(:,kvert(i,e))
      END DO
      CALL CHI_Q2_DOFMAP(e, kvert, kedge, karea, nvt, net, nat, idx)
      CALL CHI_FACE_RULE(f, xiq, wq)
      DO q = 1, 9
        CALL CHI_Q2_BASIS(xiq(:,q), phi, dphi)
        CALL CHI_FACE_GEOM(nodes, f, xiq(:,q), x, nrm, sj)
        dSw = wq(q)*sj
        ukq = 0d0
        DO j = 1, 27
          ukq(1) = ukq(1) + phi(j)*uk(idx(j))
          ukq(2) = ukq(2) + phi(j)*vk(idx(j))
          ukq(3) = ukq(3) + phi(j)*wk(idx(j))
        END DO
        un = ukq(1)*nrm(1) + ukq(2)*nrm(2) + ukq(3)*nrm(3)
        CALL robin_h(x, nrm, h)
        DO i = 1, 27
          DO a = 1, 3
            r = (a-1)*ndof + idx(i)
            rhs(r) = rhs(r) + dSw*h(a)*phi(i)
            DO j = 1, 27
              pos = csr_pos(LdA, ColA, r, (a-1)*ndof + idx(j))
              Avals(pos) = Avals(pos) - dSw*alpha*un*phi(j)*phi(i)
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE CHI_ASM_ROBIN

  !-----------------------------------------------------------------------
  ! Replace row r by the identity row with right-hand side value.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_APPLY_DIRICHLET_ROW(r, value, LdA, ColA, Avals, rhs)
    INTEGER, INTENT(IN) :: r, LdA(*), ColA(*)
    REAL*8,  INTENT(INOUT) :: Avals(*), rhs(*)
    REAL*8,  INTENT(IN) :: value

    INTEGER :: pos
    LOGICAL :: founddiag

    founddiag = .FALSE.
    DO pos = LdA(r), LdA(r+1)-1
      IF (ColA(pos) .EQ. r) THEN
        Avals(pos) = 1d0
        founddiag = .TRUE.
      ELSE
        Avals(pos) = 0d0
      END IF
    END DO
    IF (.NOT. founddiag) THEN
      WRITE(*,*) 'CHI_APPLY_DIRICHLET_ROW: no diagonal in row', r
      STOP 1
    END IF
    rhs(r) = value
  END SUBROUTINE CHI_APPLY_DIRICHLET_ROW

  !-----------------------------------------------------------------------
  ! Position of (row, col) in the sorted CSR; aborts on a structural
  ! miss (assembly bug guard).
  !-----------------------------------------------------------------------
  INTEGER FUNCTION csr_pos(LdA, ColA, r, c)
    INTEGER, INTENT(IN) :: LdA(*), ColA(*), r, c
    INTEGER :: lo, hi, mid
    lo = LdA(r)
    hi = LdA(r+1) - 1
    DO WHILE (lo .LE. hi)
      mid = (lo + hi)/2
      IF (ColA(mid) .EQ. c) THEN
        csr_pos = mid
        RETURN
      ELSE IF (ColA(mid) .LT. c) THEN
        lo = mid + 1
      ELSE
        hi = mid - 1
      END IF
    END DO
    WRITE(*,*) 'CHI_KERNELS csr_pos: entry not in pattern, row/col', r, c
    STOP 1
    csr_pos = -1
  END FUNCTION csr_pos

  PURE SUBROUTINE sort_ascending(a, n)
    INTEGER, INTENT(INOUT) :: a(*)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, k, key
    DO i = 2, n
      key = a(i)
      k = i - 1
      DO WHILE (k .GE. 1)
        IF (a(k) .LE. key) EXIT
        a(k+1) = a(k)
        k = k - 1
      END DO
      a(k+1) = key
    END DO
  END SUBROUTINE sort_ascending

END MODULE CHI_KERNELS
