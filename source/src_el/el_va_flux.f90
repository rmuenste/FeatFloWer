MODULE EL_VA_FLUX

  ! Conservative flux moments of the Euler-Lagrange volume-averaged continuity constraint
  ! (SimPar@ELVAContinuityForm = conservative; design note
  ! docs/md_docs/el_va_conservative_form_2026-11-05.md).
  !
  ! The constraint div(eps u) = -d eps/dt is enforced per element in the weak form of the four
  ! P1 pressure test functions q_1 = 1, q_m = x_m - xc_m (xc = the element-centre Q2 node, the
  ! convention of the B matrices). For a continuous transport porosity eps_h and the Q2
  ! velocity u_h this module evaluates, per element,
  !
  !   flux_moments(1) = sum_faces  int_f eps_h u_h.n dS                       (outward n)
  !   flux_moments(m) = sum_faces  int_f (x_m - xc_m) eps_h u_h.n dS
  !                     - int_e eps_h u_h,m dV                                  (m = 2..4)
  !
  ! plus the volume integrals int_e eps_h u_h dV, int_e u_h dV and the geometric moments
  ! int_e 1 dV, int_e (x_m - xc_m) dV. Both elements adjacent to an interior face evaluate the
  ! same nodal data at the same face points with opposite normals, so the sum of
  ! flux_moments(1) over all elements telescopes to the boundary flux as an algebraic identity
  ! (no neighbour lookup, no face exchange; partition faces included), and
  ! sum_e [xc_m flux_moments(1) + flux_moments(1+m)] = -int eps_h u_h,m dV + boundary term.
  !
  ! Transport porosity: eps_h is the trilinear (Q1) interpolant of vertex values, each vertex
  ! value the volume-weighted mean of the P0 element porosity of the elements sharing the
  ! vertex (EL_VA_RECONSTRUCT_EPS_LOCAL + a cross-rank nodal sum by the caller +
  ! EL_VA_EMBED_Q1_IN_Q2). It is bounded by the element values, continuous across ranks, and
  ! a reconstruction only: the P0 porosity epsilon_f stays the conserved storage variable.
  !
  ! Quadrature: 3-point Gauss per direction (faces 3x3, volume 3x3x3), the Q2 isoparametric
  ! geometry map of the B matrices (ILINT = 2, all 27 nodes). On trilinear cells the face
  ! integrands (Q1 eps x Q2 u x degree-2 cofactor) are of degree <= 5 and integrated exactly.
  ! The outward normal is n dS = sign_f sgn(det J) (dx/dxi_a x dx/dxi_b) dxi_a dxi_b with the
  ! cyclic tangent pair (a,b) of the fixed reference coordinate of face f (tables below) and the
  ! sign of the Jacobian determinant, so that left-handed (mirrored) elements, which FEAT
  ! accepts since it only uses |det J|, get the outward normal too: explicit, no centre test.
  !
  ! The module is MPI-free and holds no mesh state: the reference tables depend on nothing but
  ! the element type and are built once (EL_VA_FLUX_TABLES); the per-element kernel
  ! EL_VA_ELEMENT_FLUX_MOMENTS takes plain arrays (nodal coordinates, porosity, velocity) and
  ! is unit-tested on the synthetic mesh of el_test_helpers (test_el_va_flux). Errors are
  ! returned to the caller (nbad, ierr), which decides on a rank-uniform stop.
  !
  ! Geometry restriction: the EL meshes are straight-edged (trilinear) hexes, as the kernel
  ! deposition (EL_Q1_MAP) already assumes; EL_VA_CHECK_DOF_LAYOUT relies on it (every edge,
  ! face and centre node at the mean of its vertices) and would report a curved Q2 mesh as a
  ! layout mismatch.
  !
  ! Local Q2 dof numbering (FEAT E013): 1-8 vertices, 9-20 edges, 21-26 faces, 27 centre, in
  ! the order of kvert / kedge / karea / element; global index of local l:
  ! kvert(l,iel), nvt+kedge(l-8,iel), nvt+net+karea(l-20,iel), nvt+net+nat+iel
  ! (EL_VA_ELEMENT_DOFS; the same layout NDFGL produces for element type 13).

  USE TYPES, ONLY: tMultiMesh

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: EL_VA_RECONSTRUCT_EPS_LOCAL, EL_VA_EMBED_Q1_IN_Q2
  PUBLIC :: EL_VA_FLUX_MOMENTS, EL_VA_ELEMENT_FLUX_MOMENTS, EL_VA_ELEMENT_FACE_FLUX
  PUBLIC :: EL_VA_ELEMENT_DOFS, EL_VA_Q2_BASIS, EL_VA_FLUX_TABLES, tElVaFluxTables
  PUBLIC :: EL_VA_CHECK_DOF_LAYOUT, EL_VA_CONSERVATIVE_ROWS
  PUBLIC :: NQ2, NFACE, NGPF, NGPV

  INTEGER, PARAMETER :: NQ2 = 27      ! Q2 dofs per element
  INTEGER, PARAMETER :: NFACE = 6
  INTEGER, PARAMETER :: NGP1 = 3      ! Gauss points per direction
  INTEGER, PARAMETER :: NGPF = 9      ! per face
  INTEGER, PARAMETER :: NGPV = 27     ! per volume

  ! Reference node of each local Q2 dof as a triple in {-1,0,+1} (FEAT E013 numbering).
  INTEGER, PARAMETER :: node_ijk(3,NQ2) = RESHAPE((/ &
    -1,-1,-1,   1,-1,-1,   1, 1,-1,  -1, 1,-1, &     ! 1-4   bottom vertices
    -1,-1, 1,   1,-1, 1,   1, 1, 1,  -1, 1, 1, &     ! 5-8   top vertices
     0,-1,-1,   1, 0,-1,   0, 1,-1,  -1, 0,-1, &     ! 9-12  bottom edges
    -1,-1, 0,   1,-1, 0,   1, 1, 0,  -1, 1, 0, &     ! 13-16 vertical edges
     0,-1, 1,   1, 0, 1,   0, 1, 1,  -1, 0, 1, &     ! 17-20 top edges
     0, 0,-1,   0,-1, 0,   1, 0, 0,   0, 1, 0, &     ! 21-24 faces 1-4
    -1, 0, 0,   0, 0, 1,   0, 0, 0 /), (/3,NQ2/))    ! 25-26 faces 5-6, 27 centre

  ! Reference faces (FEAT): 1 zeta=-1, 2 eta=-1, 3 xi=+1, 4 eta=+1, 5 xi=-1, 6 zeta=+1.
  INTEGER, PARAMETER :: face_axis(NFACE) = (/ 3, 2, 1, 2, 1, 3 /)
  REAL*8, PARAMETER :: face_coord(NFACE) = (/ -1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, 1.0d0 /)
  ! cyclic tangent pair (a,b) with e_a x e_b = e_axis, and the outward sign
  INTEGER, PARAMETER :: face_ta(NFACE) = (/ 1, 3, 2, 3, 2, 1 /)
  INTEGER, PARAMETER :: face_tb(NFACE) = (/ 2, 1, 3, 1, 3, 2 /)
  REAL*8, PARAMETER :: face_sign(NFACE) = (/ -1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, 1.0d0 /)

  ! Vertices (local 1-8) defining the edge, face and centre dofs: the trilinear interpolant
  ! takes the mean of these vertex values at the node (edge midpoint, face centre, centre).
  INTEGER, PARAMETER :: edge_vertices(2,12) = RESHAPE((/ &
    1,2, 2,3, 4,3, 1,4, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 8,7, 5,8 /), (/2,12/))
  INTEGER, PARAMETER :: face_vertices(4,NFACE) = RESHAPE((/ &
    1,2,3,4,  1,2,6,5,  2,3,7,6,  4,3,7,8,  1,4,8,5,  5,6,7,8 /), (/4,NFACE/))

  TYPE tElVaFluxTables
    LOGICAL :: ready = .FALSE.
    REAL*8 :: fxi(3,NGPF,NFACE)          ! reference coordinates of the face points
    REAL*8 :: fw(NGPF)                   ! face weights (sum 4 = reference face area)
    REAL*8 :: fval(NQ2,NGPF,NFACE)       ! basis values at the face points
    REAL*8 :: fder(3,NQ2,NGPF,NFACE)     ! reference derivatives at the face points
    REAL*8 :: vxi(3,NGPV)                ! volume points
    REAL*8 :: vw(NGPV)                   ! volume weights (sum 8)
    REAL*8 :: vval(NQ2,NGPV)
    REAL*8 :: vder(3,NQ2,NGPV)
  END TYPE tElVaFluxTables

  TYPE(tElVaFluxTables), SAVE :: tables

CONTAINS

  ! Self-check of the dof layout assumed by EL_VA_ELEMENT_DOFS against the coordinates:
  ! every edge, face and centre dof must sit at the mean of its defining vertices
  ! (straight-edged hexes). Returns the number of offending dofs of the first offending
  ! element (0 = consistent) and prints the real pairing of that element; the caller
  ! (EL_TRANSFER, once per mesh, rank-uniform via a reduction) stops the run, because the
  ! flux kernel would otherwise gather the wrong nodes silently.
  INTEGER FUNCTION EL_VA_CHECK_DOF_LAYOUT(mesh, ilev, coor) RESULT(nbad)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: coor(:,:)

    INTEGER :: iel, l, dofs(NQ2), a, b, c
    REAL*8 :: v(3,8), mid(3), h, tol

    nbad = 0
    DO iel = 1, mesh%level(ilev)%nel
      CALL EL_VA_ELEMENT_DOFS(mesh, ilev, iel, dofs)
      DO l = 1, 8
        v(:,l) = coor(:,dofs(l))
      END DO
      h = MAXVAL(MAXVAL(v,DIM=2) - MINVAL(v,DIM=2))
      tol = 1.0d-6*h
      DO l = 1, 12
        mid = 0.5d0*(v(:,edge_vertices(1,l)) + v(:,edge_vertices(2,l)))
        IF (MAXVAL(ABS(coor(:,dofs(8+l)) - mid)).GT.tol) nbad = nbad + 1
      END DO
      DO l = 1, NFACE
        mid = 0.25d0*(v(:,face_vertices(1,l)) + v(:,face_vertices(2,l)) + &
                      v(:,face_vertices(3,l)) + v(:,face_vertices(4,l)))
        IF (MAXVAL(ABS(coor(:,dofs(20+l)) - mid)).GT.tol) nbad = nbad + 1
      END DO
      mid = 0.125d0*SUM(v, DIM=2)
      IF (MAXVAL(ABS(coor(:,dofs(NQ2)) - mid)).GT.tol) nbad = nbad + 1
      IF (nbad.GT.0) THEN
        WRITE(*,'(A,I0,A,I0)') 'EL_VA_FLUX dof layout mismatch in element ', iel, &
          ', bad dofs: ', nbad
        DO l = 1, 12
          CALL find_pair(coor(:,dofs(8+l)), v, tol, a, b)
          WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '  edge ', l, ': assumed ', &
            edge_vertices(1,l), '-', edge_vertices(2,l), ' found ', a, '-', b
        END DO
        DO l = 1, NFACE
          CALL find_face(coor(:,dofs(20+l)), v, tol, c)
          WRITE(*,'(A,I0,A,I0,A,I0)') '  face ', l, ': assumed ', l, ' found ', c
        END DO
        RETURN
      END IF
    END DO

  CONTAINS

    SUBROUTINE find_pair(x, v, tol, a, b)
      REAL*8, INTENT(IN) :: x(3), v(3,8), tol
      INTEGER, INTENT(OUT) :: a, b
      INTEGER :: i, j
      a = 0
      b = 0
      DO i = 1, 8
        DO j = i+1, 8
          IF (MAXVAL(ABS(x - 0.5d0*(v(:,i)+v(:,j)))).LE.tol) THEN
            a = i
            b = j
          END IF
        END DO
      END DO
    END SUBROUTINE find_pair

    SUBROUTINE find_face(x, v, tol, c)
      REAL*8, INTENT(IN) :: x(3), v(3,8), tol
      INTEGER, INTENT(OUT) :: c
      INTEGER :: f
      c = 0
      DO f = 1, NFACE
        IF (MAXVAL(ABS(x - 0.25d0*(v(:,face_vertices(1,f)) + v(:,face_vertices(2,f)) + &
            v(:,face_vertices(3,f)) + v(:,face_vertices(4,f))))).LE.tol) c = f
      END DO
    END SUBROUTINE find_face

  END FUNCTION EL_VA_CHECK_DOF_LAYOUT

  ! Row algebra of the conservative constraint for one element (design note): with the flux
  ! moments F_m of the current velocity, the storage rate d_cont, the storage moments Q_m, the
  ! assembled divergence moments D_m and the preconditioning coefficient eps_e,
  !   r_m  = F_m + d_cont Q_m          unsinked residual
  !   gV_1 = D_1 - r_1/eps_e           volume integral of the row-1 source (sink bookkeeping)
  ! and the projection right-hand side is -r_m/(eps_e dt) (row 1 minus the sink), applied by
  ! the caller. At the fixed point r_m = 0 (compatible storage), i.e. gV_1 = D_1.
  PURE SUBROUTINE EL_VA_CONSERVATIVE_ROWS(flux, d_cont, q, d, eps, r, gv1)

    REAL*8, INTENT(IN) :: flux(4), d_cont, q(4), d(4), eps
    REAL*8, INTENT(OUT) :: r(4), gv1

    r = flux + d_cont*q
    gv1 = d(1) - r(1)/eps

  END SUBROUTINE EL_VA_CONSERVATIVE_ROWS

  ! ---- Q2 Lagrange basis on [-1,1]^3 in the FEAT E013 local numbering ---------------------
  ! val(l) and dref(:,l) = d val(l) / d xi. Tensor products of the 1D quadratic Lagrange
  ! polynomials on the nodes -1, 0, +1; reproduces E013/E013A (checked by the unit test).
  PURE SUBROUTINE EL_VA_Q2_BASIS(xi, val, dref)

    REAL*8, INTENT(IN) :: xi(3)
    REAL*8, INTENT(OUT) :: val(NQ2), dref(3,NQ2)

    REAL*8 :: l1(3,-1:1), d1(3,-1:1)
    INTEGER :: k, l, i, j, m

    DO k = 1, 3
      l1(k,-1) = 0.5d0*xi(k)*(xi(k)-1.0d0)
      l1(k, 0) = 1.0d0 - xi(k)*xi(k)
      l1(k, 1) = 0.5d0*xi(k)*(xi(k)+1.0d0)
      d1(k,-1) = xi(k) - 0.5d0
      d1(k, 0) = -2.0d0*xi(k)
      d1(k, 1) = xi(k) + 0.5d0
    END DO
    DO l = 1, NQ2
      i = node_ijk(1,l)
      j = node_ijk(2,l)
      m = node_ijk(3,l)
      val(l) = l1(1,i)*l1(2,j)*l1(3,m)
      dref(1,l) = d1(1,i)*l1(2,j)*l1(3,m)
      dref(2,l) = l1(1,i)*d1(2,j)*l1(3,m)
      dref(3,l) = l1(1,i)*l1(2,j)*d1(3,m)
    END DO

  END SUBROUTINE EL_VA_Q2_BASIS

  ! Reference tables: Gauss points, weights and basis values/derivatives at the face and
  ! volume points. Independent of the mesh, built once.
  SUBROUTINE EL_VA_FLUX_TABLES(t)

    TYPE(tElVaFluxTables), INTENT(OUT) :: t

    REAL*8, PARAMETER :: g = 0.774596669241483377d0
    REAL*8 :: gp(NGP1), gw(NGP1), xi(3)
    INTEGER :: f, i, j, k, p, ax, a, b

    gp = (/ -g, 0.0d0, g /)
    gw = (/ 5.0d0/9.0d0, 8.0d0/9.0d0, 5.0d0/9.0d0 /)

    DO f = 1, NFACE
      ax = face_axis(f)
      ! the two free axes in ascending order
      a = MOD(ax,3) + 1
      b = MOD(ax+1,3) + 1
      IF (a.GT.b) THEN
        k = a
        a = b
        b = k
      END IF
      p = 0
      DO j = 1, NGP1
        DO i = 1, NGP1
          p = p + 1
          xi(ax) = face_coord(f)
          xi(a) = gp(i)
          xi(b) = gp(j)
          t%fxi(:,p,f) = xi
          t%fw(p) = gw(i)*gw(j)
          CALL EL_VA_Q2_BASIS(xi, t%fval(:,p,f), t%fder(:,:,p,f))
        END DO
      END DO
    END DO

    p = 0
    DO k = 1, NGP1
      DO j = 1, NGP1
        DO i = 1, NGP1
          p = p + 1
          xi = (/ gp(i), gp(j), gp(k) /)
          t%vxi(:,p) = xi
          t%vw(p) = gw(i)*gw(j)*gw(k)
          CALL EL_VA_Q2_BASIS(xi, t%vval(:,p), t%vder(:,:,p))
        END DO
      END DO
    END DO
    t%ready = .TRUE.

  END SUBROUTINE EL_VA_FLUX_TABLES

  ! Global Q2 dof of every local dof of element iel (layout of NDFGL, element type 13).
  PURE SUBROUTINE EL_VA_ELEMENT_DOFS(mesh, ilev, iel, dofs)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, iel
    INTEGER, INTENT(OUT) :: dofs(NQ2)

    INTEGER :: l, nvt, net, nat

    nvt = mesh%level(ilev)%nvt
    net = mesh%level(ilev)%net
    nat = mesh%level(ilev)%nat
    DO l = 1, 8
      dofs(l) = mesh%level(ilev)%kvert(l,iel)
    END DO
    DO l = 1, 12
      dofs(8+l) = nvt + mesh%level(ilev)%kedge(l,iel)
    END DO
    DO l = 1, NFACE
      dofs(20+l) = nvt + net + mesh%level(ilev)%karea(l,iel)
    END DO
    dofs(NQ2) = nvt + net + nat + iel

  END SUBROUTINE EL_VA_ELEMENT_DOFS

  ! Flux moments of one face f of an element: moments(1) = int_f eps u.n dS (outward),
  ! moments(1+m) = int_f (x_m - xc_m) eps u.n dS, and the unweighted flux int_f u.n dS.
  ! Exposed for the pairwise-cancellation test.
  PURE SUBROUTINE EL_VA_ELEMENT_FACE_FLUX(t, f, x, eps, u, xc, moments, flux_unweighted)

    TYPE(tElVaFluxTables), INTENT(IN) :: t
    INTEGER, INTENT(IN) :: f
    REAL*8, INTENT(IN) :: x(3,NQ2), eps(NQ2), u(3,NQ2), xc(3)
    REAL*8, INTENT(OUT) :: moments(4), flux_unweighted

    REAL*8 :: jac(3,3), nds(3), xq(3), uq(3), epsq, un, un0, detj
    INTEGER :: p, l, k, a, b

    moments = 0.0d0
    flux_unweighted = 0.0d0
    a = face_ta(f)
    b = face_tb(f)
    DO p = 1, NGPF
      jac = 0.0d0
      xq = 0.0d0
      uq = 0.0d0
      epsq = 0.0d0
      DO l = 1, NQ2
        DO k = 1, 3
          jac(:,k) = jac(:,k) + x(:,l)*t%fder(k,l,p,f)
        END DO
        xq = xq + x(:,l)*t%fval(l,p,f)
        uq = uq + u(:,l)*t%fval(l,p,f)
        epsq = epsq + eps(l)*t%fval(l,p,f)
      END DO
      ! outward n dS = sign sgn(det J) (dx/dxi_a x dx/dxi_b) w
      detj = jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3)) &
           - jac(2,1)*(jac(1,2)*jac(3,3)-jac(3,2)*jac(1,3)) &
           + jac(3,1)*(jac(1,2)*jac(2,3)-jac(2,2)*jac(1,3))
      nds(1) = jac(2,a)*jac(3,b) - jac(3,a)*jac(2,b)
      nds(2) = jac(3,a)*jac(1,b) - jac(1,a)*jac(3,b)
      nds(3) = jac(1,a)*jac(2,b) - jac(2,a)*jac(1,b)
      nds = face_sign(f)*SIGN(1.0d0, detj)*t%fw(p)*nds
      un0 = DOT_PRODUCT(uq, nds)
      un = epsq*un0
      flux_unweighted = flux_unweighted + un0
      moments(1) = moments(1) + un
      moments(2:4) = moments(2:4) + (xq - xc)*un
    END DO

  END SUBROUTINE EL_VA_ELEMENT_FACE_FLUX

  ! Flux moments of one element (see the module header).
  !   x(3,27)   nodal coordinates (Q2 geometry), eps(27) transport porosity, u(3,27) velocity
  !   flux(4)   conservative weak-form moments
  !   eps_u(3)  int_e eps_h u_h dV,  vol_u(3) int_e u_h dV
  !   geom(0:3) int_e 1 dV, int_e (x_m - xc_m) dV
  !   div_u     sum_faces int_f u_h.n dS (= int_e div u_h, the consistency check against B)
  PURE SUBROUTINE EL_VA_ELEMENT_FLUX_MOMENTS(t, x, eps, u, flux, eps_u, vol_u, geom, div_u)

    TYPE(tElVaFluxTables), INTENT(IN) :: t
    REAL*8, INTENT(IN) :: x(3,NQ2), eps(NQ2), u(3,NQ2)
    REAL*8, INTENT(OUT) :: flux(4), eps_u(3), vol_u(3), geom(0:3), div_u

    REAL*8 :: xc(3), moments(4), jac(3,3), detj, w, xq(3), uq(3), epsq, un0
    INTEGER :: f, p, l, k

    xc = x(:,NQ2)
    flux = 0.0d0
    div_u = 0.0d0
    DO f = 1, NFACE
      CALL EL_VA_ELEMENT_FACE_FLUX(t, f, x, eps, u, xc, moments, un0)
      flux = flux + moments
      div_u = div_u + un0
    END DO

    eps_u = 0.0d0
    vol_u = 0.0d0
    geom = 0.0d0
    DO p = 1, NGPV
      jac = 0.0d0
      xq = 0.0d0
      uq = 0.0d0
      epsq = 0.0d0
      DO l = 1, NQ2
        DO k = 1, 3
          jac(:,k) = jac(:,k) + x(:,l)*t%vder(k,l,p)
        END DO
        xq = xq + x(:,l)*t%vval(l,p)
        uq = uq + u(:,l)*t%vval(l,p)
        epsq = epsq + eps(l)*t%vval(l,p)
      END DO
      detj = jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3)) &
           - jac(2,1)*(jac(1,2)*jac(3,3)-jac(3,2)*jac(1,3)) &
           + jac(3,1)*(jac(1,2)*jac(2,3)-jac(2,2)*jac(1,3))
      w = t%vw(p)*ABS(detj)
      geom(0) = geom(0) + w
      geom(1:3) = geom(1:3) + (xq - xc)*w
      vol_u = vol_u + uq*w
      eps_u = eps_u + epsq*uq*w
    END DO
    flux(2:4) = flux(2:4) - eps_u

  END SUBROUTINE EL_VA_ELEMENT_FLUX_MOMENTS

  ! Driver over the local elements. coor(3,:) must hold the coordinates of every Q2 dof
  ! (dcorvg of the EL coupling level, or the BUILD_DOF_COORDINATES array of the tests).
  ! ierr = 1 when an input array is shorter than the Q2 dof count (nothing evaluated).
  SUBROUTINE EL_VA_FLUX_MOMENTS(mesh, ilev, coor, eps_q2, val_u, val_v, val_w, &
                                flux_moments, eps_velocity, velocity, geom_moments, div_u, ierr)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: coor(:,:), eps_q2(:), val_u(:), val_v(:), val_w(:)
    REAL*8, INTENT(OUT) :: flux_moments(:,:)    ! (4,nel)
    REAL*8, INTENT(OUT) :: eps_velocity(:,:)    ! (3,nel)
    REAL*8, INTENT(OUT) :: velocity(:,:)        ! (3,nel)
    REAL*8, INTENT(OUT) :: geom_moments(0:,:)   ! (0:3,nel)
    REAL*8, INTENT(OUT) :: div_u(:)             ! (nel)
    INTEGER, INTENT(OUT) :: ierr

    INTEGER :: iel, l, dofs(NQ2), ndof
    REAL*8 :: x(3,NQ2), eps(NQ2), u(3,NQ2)

    IF (.NOT.tables%ready) CALL EL_VA_FLUX_TABLES(tables)

    ierr = 0
    ndof = mesh%level(ilev)%nvt + mesh%level(ilev)%net + mesh%level(ilev)%nat + &
           mesh%level(ilev)%nel
    IF (SIZE(coor,2).LT.ndof .OR. SIZE(eps_q2).LT.ndof .OR. SIZE(val_u).LT.ndof .OR. &
        SIZE(val_v).LT.ndof .OR. SIZE(val_w).LT.ndof) THEN
      ierr = 1
      RETURN
    END IF

    DO iel = 1, mesh%level(ilev)%nel
      CALL EL_VA_ELEMENT_DOFS(mesh, ilev, iel, dofs)
      DO l = 1, NQ2
        x(:,l) = coor(:,dofs(l))
        eps(l) = eps_q2(dofs(l))
        u(1,l) = val_u(dofs(l))
        u(2,l) = val_v(dofs(l))
        u(3,l) = val_w(dofs(l))
      END DO
      CALL EL_VA_ELEMENT_FLUX_MOMENTS(tables, x, eps, u, flux_moments(:,iel), &
        eps_velocity(:,iel), velocity(:,iel), geom_moments(:,iel), div_u(iel))
    END DO

  END SUBROUTINE EL_VA_FLUX_MOMENTS

  ! Local part of the vertex reconstruction: num(v) = sum_e eps_e |e|, den(v) = sum_e |e|
  ! over the local elements sharing vertex v (|e| = dvol, the storage volume of the
  ! deposition). Vertex entries of Q2-sized arrays; the other entries stay zero. The caller
  ! sums both arrays across ranks (E013Sum2) and divides.
  PURE SUBROUTINE EL_VA_RECONSTRUCT_EPS_LOCAL(mesh, ilev, epsilon_f, num, den)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: epsilon_f(:)
    REAL*8, INTENT(OUT) :: num(:), den(:)

    INTEGER :: iel, l, iv
    REAL*8 :: vol

    num = 0.0d0
    den = 0.0d0
    DO iel = 1, mesh%level(ilev)%nel
      vol = mesh%level(ilev)%dvol(iel)
      DO l = 1, 8
        iv = mesh%level(ilev)%kvert(l,iel)
        num(iv) = num(iv) + epsilon_f(iel)*vol
        den(iv) = den(iv) + vol
      END DO
    END DO

  END SUBROUTINE EL_VA_RECONSTRUCT_EPS_LOCAL

  ! Trilinear embedding: with the vertex entries of eps_q2 set, fill the edge, face and
  ! centre entries with the means of their defining vertices (exact Q1 interpolant in Q2).
  PURE SUBROUTINE EL_VA_EMBED_Q1_IN_Q2(mesh, ilev, eps_q2)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(INOUT) :: eps_q2(:)

    INTEGER :: iel, l, dofs(NQ2)
    REAL*8 :: v(8)

    DO iel = 1, mesh%level(ilev)%nel
      CALL EL_VA_ELEMENT_DOFS(mesh, ilev, iel, dofs)
      DO l = 1, 8
        v(l) = eps_q2(dofs(l))
      END DO
      DO l = 1, 12
        eps_q2(dofs(8+l)) = 0.5d0*(v(edge_vertices(1,l)) + v(edge_vertices(2,l)))
      END DO
      DO l = 1, NFACE
        eps_q2(dofs(20+l)) = 0.25d0*(v(face_vertices(1,l)) + v(face_vertices(2,l)) + &
                                     v(face_vertices(3,l)) + v(face_vertices(4,l)))
      END DO
      eps_q2(dofs(NQ2)) = 0.125d0*SUM(v)
    END DO

  END SUBROUTINE EL_VA_EMBED_Q1_IN_Q2

END MODULE EL_VA_FLUX
