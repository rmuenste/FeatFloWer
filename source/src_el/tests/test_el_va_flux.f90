PROGRAM TEST_EL_VA_FLUX

  ! Unit tests of EL_VA_FLUX (conservative continuity form) on the synthetic structured mesh
  ! of el_test_helpers. MPI-free kernel; MPI_Init only for the shared test harness.
  !   1 basis numbering: EL_VA_Q2_BASIS reproduces the FEAT E013A values and derivatives
  !   2 one-cell geometry: face areas, outward orientation, sum of fluxes of u = x is 3|e|
  !   3 shared face: the face fluxes of the two adjacent elements cancel pairwise
  !   4 telescoping and the three coordinate-moment identities (brick and warped mesh)
  !   5 constant eps and linear velocity: F_1 = eps |e|_geom div u exactly
  !   6 reconstruction: constant and linear fields reproduced, bounds across a jump
  !   7 dof-layout check: consistent mesh passes, a swapped edge pair is reported
  !   8 row algebra: compatible flux gives zero residual, incompatible storage is exposed

  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE EL_VA_FLUX, ONLY: EL_VA_Q2_BASIS, EL_VA_FLUX_TABLES, tElVaFluxTables, &
                        EL_VA_ELEMENT_DOFS, EL_VA_ELEMENT_FACE_FLUX, &
                        EL_VA_ELEMENT_FLUX_MOMENTS, EL_VA_FLUX_MOMENTS, &
                        EL_VA_RECONSTRUCT_EPS_LOCAL, EL_VA_EMBED_Q1_IN_Q2, NQ2, &
                        EL_VA_CHECK_DOF_LAYOUT, EL_VA_CONSERVATIVE_ROWS
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, EL_TEST_CHECK, EL_TEST_REPORT, &
                             EL_TEST_RELEASE_MESH, BUILD_STRUCTURED_MESH, &
                             BUILD_DOF_COORDINATES, nel, nvt, net, nat

  IMPLICIT NONE

  INTEGER :: rank, nproc, ierr, failures
  TYPE(tElVaFluxTables) :: tab

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  CALL EL_TEST_INIT()
  CALL EL_VA_FLUX_TABLES(tab)

  failures = 0
  CALL test_basis_numbering(failures)
  CALL test_one_cell_geometry(failures, .FALSE.)
  CALL test_one_cell_geometry(failures, .TRUE.)
  CALL test_mirrored_cell(failures)
  CALL test_shared_face_cancellation(failures)
  CALL test_telescoping_and_identities(failures, .FALSE.)
  CALL test_telescoping_and_identities(failures, .TRUE.)
  CALL test_constant_eps_linear_velocity(failures)
  CALL test_reconstruction(failures)
  CALL test_layout_check(failures)
  CALL test_conservative_rows(failures)

  CALL EL_TEST_REPORT(failures, rank, nproc, 'EL VA flux')
  CALL MPI_Finalize(ierr)

  IF (failures.GT.0) ERROR STOP 1

CONTAINS

  ! Mesh of n^3 cells on the unit cube; warped moves every interior vertex by a smooth,
  ! deterministic displacement of up to 0.2 h. Q2 dof coordinates from the test helper.
  SUBROUTINE setup_mesh(n, warped, mesh, coor, ndof)

    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: warped
    TYPE(tMultiMesh), INTENT(INOUT) :: mesh
    REAL*8, ALLOCATABLE, INTENT(OUT) :: coor(:,:)
    INTEGER, INTENT(OUT) :: ndof

    INTEGER :: iv
    REAL*8 :: x(3), h

    CALL BUILD_STRUCTURED_MESH(mesh, n)
    ndof = nvt + net + nat + nel
    h = 1.0d0/DBLE(n)
    IF (warped) THEN
      DO iv = 1, nvt
        x = mesh%level(1)%dcorvg(:,iv)
        IF (MINVAL(x).GT.1.0d-12 .AND. MAXVAL(x).LT.1.0d0-1.0d-12) THEN
          mesh%level(1)%dcorvg(1,iv) = x(1) + 0.2d0*h*SIN(3.1d0*x(2)+1.3d0*x(3))
          mesh%level(1)%dcorvg(2,iv) = x(2) + 0.15d0*h*COS(2.2d0*x(1)+0.7d0*x(3))
          mesh%level(1)%dcorvg(3,iv) = x(3) + 0.18d0*h*SIN(1.7d0*x(1)-2.1d0*x(2))
        END IF
      END DO
    END IF
    ALLOCATE(coor(3,ndof))
    CALL BUILD_DOF_COORDINATES(mesh, coor)

  END SUBROUTINE setup_mesh

  ! Nodal Q2 velocity with zero normal component on the walls of the unit cube.
  SUBROUTINE wall_tangent_velocity(coor, u, v, w)

    REAL*8, INTENT(IN) :: coor(:,:)
    REAL*8, INTENT(OUT) :: u(:), v(:), w(:)
    INTEGER :: i
    REAL*8 :: x, y, z

    DO i = 1, SIZE(coor,2)
      x = coor(1,i)
      y = coor(2,i)
      z = coor(3,i)
      u(i) = x*(1.0d0-x)*(0.3d0 + y*z + 0.5d0*SIN(2.0d0*z))
      v(i) = y*(1.0d0-y)*(0.7d0*x - z*z + 0.2d0*COS(3.0d0*x))
      w(i) = z*(1.0d0-z)*(1.0d0 + x*y - 0.4d0*x*x)
    END DO

  END SUBROUTINE wall_tangent_velocity

  ! Deterministic pseudo-random nodal porosity in [0.4, 1].
  SUBROUTINE random_eps(n, eps)

    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(OUT) :: eps(:)
    INTEGER :: i

    DO i = 1, n
      eps(i) = 0.7d0 + 0.3d0*SIN(12.9898d0*DBLE(i) + 78.233d0*DBLE(MOD(i,7)))
    END DO

  END SUBROUTINE random_eps

  SUBROUTINE gather_element(mesh, iel, coor, eps, u, v, w, x, e, uu)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: iel
    REAL*8, INTENT(IN) :: coor(:,:), eps(:), u(:), v(:), w(:)
    REAL*8, INTENT(OUT) :: x(3,NQ2), e(NQ2), uu(3,NQ2)
    INTEGER :: dofs(NQ2), l

    CALL EL_VA_ELEMENT_DOFS(mesh, 1, iel, dofs)
    DO l = 1, NQ2
      x(:,l) = coor(:,dofs(l))
      e(l) = eps(dofs(l))
      uu(1,l) = u(dofs(l))
      uu(2,l) = v(dofs(l))
      uu(3,l) = w(dofs(l))
    END DO

  END SUBROUTINE gather_element

  ! 1. The tensor-product basis of the module against the FEAT E013A (values and reference
  !    derivatives), i.e. the local numbering assumed by EL_VA_ELEMENT_DOFS.
  SUBROUTINE test_basis_numbering(nf)

    INTEGER, INTENT(INOUT) :: nf

    INTEGER, PARAMETER :: NNBAS = 27, NNDER = 10, NNVE = 8
    REAL*8 :: dx(NNVE), dy(NNVE), dz(NNVE), djac(3,3), detj, dbas(3,NNBAS,NNDER)
    LOGICAL :: bder(NNDER)
    INTEGER :: kve(NNVE), iel_c, ndim
    COMMON /ELEM/ dx, dy, dz, djac, detj, dbas, bder, kve, iel_c, ndim
    REAL*8 :: dhelp(NNBAS,4,1), xi(3), val(NQ2), dref(3,NQ2), err
    INTEGER :: ip, k
    REAL*8, PARAMETER :: pts(3,4) = RESHAPE((/ 0.31d0,-0.52d0,0.77d0, -0.9d0,0.1d0,-0.2d0, &
                                               0.0d0,0.0d0,0.0d0, 0.66d0,0.66d0,-0.66d0 /), (/3,4/))
    EXTERNAL E013A

    bder = .FALSE.
    bder(1:4) = .TRUE.
    err = 0.0d0
    DO ip = 1, 4
      xi = pts(:,ip)
      CALL E013A(xi(1), xi(2), xi(3), dhelp, 1)
      CALL EL_VA_Q2_BASIS(xi, val, dref)
      err = MAX(err, MAXVAL(ABS(val - dhelp(:,1,1))))
      DO k = 1, 3
        err = MAX(err, MAXVAL(ABS(dref(k,:) - dhelp(:,1+k,1))))
      END DO
    END DO
    CALL EL_TEST_CHECK(err.LT.1.0d-13, 'Q2 basis and numbering match E013A', nf, rank)

  END SUBROUTINE test_basis_numbering

  ! 2. One cell: each face area and outward orientation from a uniform velocity, and the
  !    divergence theorem for u = x (sum of face fluxes = 3 |e|_geom, exact on trilinear cells).
  SUBROUTINE test_one_cell_geometry(nf, warped)

    INTEGER, INTENT(INOUT) :: nf
    LOGICAL, INTENT(IN) :: warped

    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:), eps(:), u(:), v(:), w(:)
    INTEGER :: ndof, f, k, l
    REAL*8 :: x(3,NQ2), e(NQ2), uu(3,NQ2), mom(4), flux(4), eps_u(3), vol_u(3), geom(0:3)
    REAL*8 :: err_area, err_sign, expected, un0, div_u
    ! expected outward unit-normal direction per face (brick): -z, -y, +x, +y, -x, +z
    INTEGER, PARAMETER :: naxis(6) = (/ 3, 2, 1, 2, 1, 3 /)
    REAL*8, PARAMETER :: nsign(6) = (/ -1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, 1.0d0 /)
    CHARACTER(LEN=8) :: tag

    tag = 'brick'
    IF (warped) tag = 'warped'
    ! a single cell has no interior vertex; warp the cell by hand instead
    CALL setup_mesh(1, .FALSE., mesh, coor, ndof)
    IF (warped) THEN
      ! vertex numbering of the helper: vid(i,j,k) = 1 + i + 2j + 4k; move two corners
      mesh%level(1)%dcorvg(:,8) = (/ 1.1d0, 0.9d0, 1.15d0 /)    ! (1,1,1)
      mesh%level(1)%dcorvg(:,4) = (/ 1.05d0, 0.95d0, 0.1d0 /)   ! (1,1,0)
      CALL BUILD_DOF_COORDINATES(mesh, coor)
    END IF
    ALLOCATE(eps(ndof), u(ndof), v(ndof), w(ndof))
    eps = 1.0d0
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    CALL gather_element(mesh, 1, coor, eps, u, v, w, x, e, uu)

    IF (.NOT.warped) THEN
      err_area = 0.0d0
      err_sign = 0.0d0
      DO f = 1, 6
        DO k = 1, 3
          uu = 0.0d0
          uu(k,:) = 1.0d0
          CALL EL_VA_ELEMENT_FACE_FLUX(tab, f, x, e, uu, x(:,NQ2), mom, un0)
          expected = 0.0d0
          IF (k.EQ.naxis(f)) expected = nsign(f)
          err_sign = MAX(err_sign, ABS(mom(1) - expected))
        END DO
      END DO
      CALL EL_TEST_CHECK(err_sign.LT.1.0d-13, 'unit cell: face areas and outward normals', nf, rank)
    END IF

    DO l = 1, NQ2
      uu(:,l) = x(:,l)
    END DO
    CALL EL_VA_ELEMENT_FLUX_MOMENTS(tab, x, e, uu, flux, eps_u, vol_u, geom, div_u)
    CALL EL_TEST_CHECK(ABS(flux(1) - 3.0d0*geom(0)).LT.1.0d-13*MAX(1.0d0,geom(0)) .AND. &
      ABS(div_u - flux(1)).LT.1.0d-13, &
      'unit cell '//TRIM(tag)//': divergence theorem for u = x', nf, rank)
    CALL EL_TEST_CHECK(ABS(vol_u(1) - geom(1) - x(1,NQ2)*geom(0)).LT.1.0d-13, &
      'unit cell '//TRIM(tag)//': volume moments consistent', nf, rank)

    CALL EL_TEST_RELEASE_MESH(mesh)
    DEALLOCATE(coor, eps, u, v, w)

  END SUBROUTINE test_one_cell_geometry

  ! 2b. A left-handed (mirrored, det J < 0) cell: FEAT accepts such elements, so the outward
  !     normal must not depend on the orientation. Mirror x -> 1 - x of the unit cell.
  SUBROUTINE test_mirrored_cell(nf)

    INTEGER, INTENT(INOUT) :: nf

    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:), eps(:), u(:), v(:), w(:)
    INTEGER :: ndof, f, k, l
    REAL*8 :: x(3,NQ2), e(NQ2), uu(3,NQ2), mom(4), flux(4), eps_u(3), vol_u(3), geom(0:3)
    REAL*8 :: err_sign, expected, un0, div_u
    ! outward unit normals of the mirrored cell: faces 3 and 5 swap their direction
    INTEGER, PARAMETER :: naxis(6) = (/ 3, 2, 1, 2, 1, 3 /)
    REAL*8, PARAMETER :: nsign(6) = (/ -1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0 /)

    CALL setup_mesh(1, .FALSE., mesh, coor, ndof)
    mesh%level(1)%dcorvg(1,:) = 1.0d0 - mesh%level(1)%dcorvg(1,:)
    CALL BUILD_DOF_COORDINATES(mesh, coor)
    ALLOCATE(eps(ndof), u(ndof), v(ndof), w(ndof))
    eps = 1.0d0
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    CALL gather_element(mesh, 1, coor, eps, u, v, w, x, e, uu)
    err_sign = 0.0d0
    DO f = 1, 6
      DO k = 1, 3
        uu = 0.0d0
        uu(k,:) = 1.0d0
        CALL EL_VA_ELEMENT_FACE_FLUX(tab, f, x, e, uu, x(:,NQ2), mom, un0)
        expected = 0.0d0
        IF (k.EQ.naxis(f)) expected = nsign(f)
        err_sign = MAX(err_sign, ABS(mom(1) - expected))
      END DO
    END DO
    CALL EL_TEST_CHECK(err_sign.LT.1.0d-13, 'mirrored cell: outward normals', nf, rank)
    DO l = 1, NQ2
      uu(:,l) = x(:,l)
    END DO
    CALL EL_VA_ELEMENT_FLUX_MOMENTS(tab, x, e, uu, flux, eps_u, vol_u, geom, div_u)
    CALL EL_TEST_CHECK(ABS(flux(1) - 3.0d0*geom(0)).LT.1.0d-13 .AND. ABS(geom(0)-1.0d0).LT.1.0d-13, &
      'mirrored cell: divergence theorem for u = x', nf, rank)
    CALL EL_TEST_RELEASE_MESH(mesh)
    DEALLOCATE(coor, eps, u, v, w)

  END SUBROUTINE test_mirrored_cell

  ! 3. The flux through a shared face evaluated from both adjacent elements cancels exactly
  !    (same nodal data, same points, opposite normals), brick and warped.
  SUBROUTINE test_shared_face_cancellation(nf)

    INTEGER, INTENT(INOUT) :: nf

    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:), eps(:), u(:), v(:), w(:)
    INTEGER :: ndof, iw
    REAL*8 :: x1(3,NQ2), e1(NQ2), u1(3,NQ2), x2(3,NQ2), e2(NQ2), u2(3,NQ2), m1(4), m2(4)
    REAL*8 :: err, scale, un1, un2
    LOGICAL :: warped

    DO iw = 0, 1
      warped = (iw.EQ.1)
      CALL setup_mesh(3, warped, mesh, coor, ndof)
      ALLOCATE(eps(ndof), u(ndof), v(ndof), w(ndof))
      CALL random_eps(ndof, eps)
      CALL wall_tangent_velocity(coor, u, v, w)
      ! element 14 is the centre cell of the 3x3x3 mesh, element 15 its +x neighbour:
      ! face 3 of 14 is face 5 of 15
      CALL gather_element(mesh, 14, coor, eps, u, v, w, x1, e1, u1)
      CALL gather_element(mesh, 15, coor, eps, u, v, w, x2, e2, u2)
      CALL EL_VA_ELEMENT_FACE_FLUX(tab, 3, x1, e1, u1, x1(:,NQ2), m1, un1)
      CALL EL_VA_ELEMENT_FACE_FLUX(tab, 5, x2, e2, u2, x2(:,NQ2), m2, un2)
      scale = MAX(ABS(m1(1)), ABS(m2(1)), 1.0d-30)
      err = ABS(m1(1) + m2(1))/scale
      CALL EL_TEST_CHECK(ABS(m1(1)).GT.1.0d-6 .AND. err.LT.1.0d-12, &
        MERGE('warped', 'brick ', warped)//' shared face: fluxes cancel pairwise', nf, rank)
      DEALLOCATE(eps, u, v, w)
      CALL EL_TEST_RELEASE_MESH(mesh)
      DEALLOCATE(coor)
    END DO

  END SUBROUTINE test_shared_face_cancellation

  ! 4. Closed box: sum_e F_1 = 0 (telescoping) and sum_e [xc_m F_1 + F_1+m] + int eps_h u_m = 0
  !    for the three coordinates, with random nodal eps and a wall-tangent velocity.
  SUBROUTINE test_telescoping_and_identities(nf, warped)

    INTEGER, INTENT(INOUT) :: nf
    LOGICAL, INTENT(IN) :: warped

    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:), eps(:), u(:), v(:), w(:)
    REAL*8, ALLOCATABLE :: flux(:,:), eps_u(:,:), vol_u(:,:), geom(:,:), div_u(:)
    INTEGER :: ndof, iel, m, ioff, ierr
    REAL*8 :: phi, scale, ident(3), iscale(3), xc(3)
    CHARACTER(LEN=8) :: tag

    tag = 'brick'
    IF (warped) tag = 'warped'
    CALL setup_mesh(3, warped, mesh, coor, ndof)
    ALLOCATE(eps(ndof), u(ndof), v(ndof), w(ndof))
    ALLOCATE(flux(4,nel), eps_u(3,nel), vol_u(3,nel), geom(0:3,nel), div_u(nel))
    CALL random_eps(ndof, eps)
    CALL wall_tangent_velocity(coor, u, v, w)
    CALL EL_VA_FLUX_MOMENTS(mesh, 1, coor, eps, u, v, w, flux, eps_u, vol_u, geom, div_u, ierr)

    ioff = nvt + net + nat
    phi = 0.0d0
    scale = 0.0d0
    ident = 0.0d0
    iscale = 0.0d0
    DO iel = 1, nel
      xc = coor(:,ioff+iel)
      phi = phi + flux(1,iel)
      scale = scale + ABS(flux(1,iel))
      DO m = 1, 3
        ident(m) = ident(m) + xc(m)*flux(1,iel) + flux(1+m,iel) + eps_u(m,iel)
        iscale(m) = iscale(m) + ABS(eps_u(m,iel)) + ABS(xc(m)*flux(1,iel))
      END DO
    END DO
    CALL EL_TEST_CHECK(scale.GT.1.0d-6 .AND. ABS(phi).LT.1.0d-12*scale, &
      TRIM(tag)//' box: sum of element fluxes telescopes to zero', nf, rank)
    CALL EL_TEST_CHECK(ALL(ABS(ident).LT.1.0d-12*iscale), &
      TRIM(tag)//' box: coordinate-moment identities', nf, rank)
    CALL EL_TEST_CHECK(ABS(SUM(geom(0,:)) - 1.0d0).LT.1.0d-12, &
      TRIM(tag)//' box: quadrature volumes sum to the box', nf, rank)

    DEALLOCATE(eps, u, v, w, flux, eps_u, vol_u, geom, div_u)
    CALL EL_TEST_RELEASE_MESH(mesh)
    DEALLOCATE(coor)

  END SUBROUTINE test_telescoping_and_identities

  ! 5. Constant eps and a linear velocity u = A x + b: F_1 = eps |e|_geom tr(A) per element.
  SUBROUTINE test_constant_eps_linear_velocity(nf)

    INTEGER, INTENT(INOUT) :: nf

    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:), eps(:), u(:), v(:), w(:)
    REAL*8, ALLOCATABLE :: flux(:,:), eps_u(:,:), vol_u(:,:), geom(:,:), div_u(:)
    INTEGER :: ndof, i, iel, ierr
    REAL*8 :: err, tra
    REAL*8, PARAMETER :: a(3,3) = RESHAPE((/ 0.3d0, -0.2d0, 0.5d0, 0.1d0, 0.8d0, -0.4d0, &
                                             -0.6d0, 0.25d0, -0.35d0 /), (/3,3/))
    REAL*8, PARAMETER :: b(3) = (/ 0.2d0, -0.1d0, 0.4d0 /)

    CALL setup_mesh(2, .TRUE., mesh, coor, ndof)
    ALLOCATE(eps(ndof), u(ndof), v(ndof), w(ndof))
    ALLOCATE(flux(4,nel), eps_u(3,nel), vol_u(3,nel), geom(0:3,nel), div_u(nel))
    eps = 0.7d0
    DO i = 1, ndof
      u(i) = DOT_PRODUCT(a(1,:), coor(:,i)) + b(1)
      v(i) = DOT_PRODUCT(a(2,:), coor(:,i)) + b(2)
      w(i) = DOT_PRODUCT(a(3,:), coor(:,i)) + b(3)
    END DO
    CALL EL_VA_FLUX_MOMENTS(mesh, 1, coor, eps, u, v, w, flux, eps_u, vol_u, geom, div_u, ierr)
    tra = a(1,1) + a(2,2) + a(3,3)
    err = 0.0d0
    DO iel = 1, nel
      err = MAX(err, ABS(flux(1,iel) - 0.7d0*geom(0,iel)*tra)/geom(0,iel))
    END DO
    CALL EL_TEST_CHECK(err.LT.1.0d-13, 'constant eps, linear velocity: F_1 = eps |e| div u', nf, rank)

    DEALLOCATE(eps, u, v, w, flux, eps_u, vol_u, geom, div_u)
    CALL EL_TEST_RELEASE_MESH(mesh)
    DEALLOCATE(coor)

  END SUBROUTINE test_constant_eps_linear_velocity

  ! 6. Reconstruction: a linear element field is reproduced at every dof whose defining
  !    vertices are interior (the wall vertices average one-sidedly); across a jump every
  !    nodal value stays within the element range (boundedness of the Q1 embedding).
  SUBROUTINE test_reconstruction(nf)

    INTEGER, INTENT(INOUT) :: nf

    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:), eps_e(:), num(:), den(:), eps_q2(:)
    INTEGER :: ndof, iel, iv, i
    REAL*8 :: xc(3), err, x(3), h
    LOGICAL :: ok

    CALL setup_mesh(4, .FALSE., mesh, coor, ndof)
    h = 0.25d0
    ALLOCATE(eps_e(nel), num(ndof), den(ndof), eps_q2(ndof))

    ! linear field eps = 0.5 + 0.3 x - 0.1 y + 0.2 z evaluated at the element centres
    DO iel = 1, nel
      xc = coor(:, nvt+net+nat+iel)
      eps_e(iel) = 0.5d0 + 0.3d0*xc(1) - 0.1d0*xc(2) + 0.2d0*xc(3)
    END DO
    CALL EL_VA_RECONSTRUCT_EPS_LOCAL(mesh, 1, eps_e, num, den)
    eps_q2 = 0.0d0
    DO iv = 1, nvt
      IF (den(iv).GT.0.0d0) eps_q2(iv) = num(iv)/den(iv)
    END DO
    CALL EL_VA_EMBED_Q1_IN_Q2(mesh, 1, eps_q2)
    err = 0.0d0
    DO i = 1, ndof
      x = coor(:,i)
      IF (MINVAL(x).GT.h-1.0d-12 .AND. MAXVAL(x).LT.1.0d0-h+1.0d-12) &
        err = MAX(err, ABS(eps_q2(i) - (0.5d0 + 0.3d0*x(1) - 0.1d0*x(2) + 0.2d0*x(3))))
    END DO
    CALL EL_TEST_CHECK(err.LT.1.0d-13, 'reconstruction reproduces a linear field away from the walls', nf, rank)

    ! jump: 0.4 below x = 0.5, 1 above
    DO iel = 1, nel
      xc = coor(:, nvt+net+nat+iel)
      eps_e(iel) = MERGE(0.4d0, 1.0d0, xc(1).LT.0.5d0)
    END DO
    CALL EL_VA_RECONSTRUCT_EPS_LOCAL(mesh, 1, eps_e, num, den)
    eps_q2 = 0.0d0
    DO iv = 1, nvt
      IF (den(iv).GT.0.0d0) eps_q2(iv) = num(iv)/den(iv)
    END DO
    CALL EL_VA_EMBED_Q1_IN_Q2(mesh, 1, eps_q2)
    ok = MINVAL(eps_q2).GE.0.4d0-1.0d-14 .AND. MAXVAL(eps_q2).LE.1.0d0+1.0d-14
    ! the interface vertices carry the mean 0.7, the others the element value
    DO iv = 1, nvt
      x = coor(:,iv)
      IF (ABS(x(1)-0.5d0).LT.1.0d-12) THEN
        ok = ok .AND. ABS(eps_q2(iv)-0.7d0).LT.1.0d-14
      ELSE IF (x(1).LT.0.5d0) THEN
        ok = ok .AND. ABS(eps_q2(iv)-0.4d0).LT.1.0d-14
      ELSE
        ok = ok .AND. ABS(eps_q2(iv)-1.0d0).LT.1.0d-14
      END IF
    END DO
    CALL EL_TEST_CHECK(ok, 'reconstruction bounded across a jump', nf, rank)

    DEALLOCATE(eps_e, num, den, eps_q2)
    CALL EL_TEST_RELEASE_MESH(mesh)
    DEALLOCATE(coor)

  END SUBROUTINE test_reconstruction

  ! 7. The layout self-check: the consistent helper mesh passes; swapping the two edges of
  !    one element (the midpoint coordinates no longer match) is reported.
  SUBROUTINE test_layout_check(nf)

    INTEGER, INTENT(INOUT) :: nf

    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:)
    INTEGER :: ndof, nbad, k

    CALL setup_mesh(2, .TRUE., mesh, coor, ndof)
    nbad = EL_VA_CHECK_DOF_LAYOUT(mesh, 1, coor)
    CALL EL_TEST_CHECK(nbad.EQ.0, 'layout check: consistent mesh passes', nf, rank)
    k = mesh%level(1)%kedge(1,3)
    mesh%level(1)%kedge(1,3) = mesh%level(1)%kedge(2,3)
    mesh%level(1)%kedge(2,3) = k
    WRITE(*,'(A)') '  (the next layout report is expected)'
    nbad = EL_VA_CHECK_DOF_LAYOUT(mesh, 1, coor)
    CALL EL_TEST_CHECK(nbad.EQ.2, 'layout check: swapped edge pair reported', nf, rank)
    CALL EL_TEST_RELEASE_MESH(mesh)
    DEALLOCATE(coor)

  END SUBROUTINE test_layout_check

  ! 8. Row algebra: a flux that balances the storage (F = -d_cont Q) gives r = 0 and
  !    gV_1 = D_1; an incompatible storage rate shows up in r_1 and is not absorbed.
  SUBROUTINE test_conservative_rows(nf)

    INTEGER, INTENT(INOUT) :: nf

    REAL*8 :: flux(4), q(4), d(4), r(4), gv1, d_cont, eps
    REAL*8, PARAMETER :: tol = 1.0d-15

    q = (/ 2.0d-3, 1.0d-5, -3.0d-5, 2.0d-5 /)
    d = (/ -7.0d-7, 1.0d-8, 2.0d-8, -4.0d-8 /)
    d_cont = 3.5d-4
    eps = 0.93d0
    flux = -d_cont*q
    CALL EL_VA_CONSERVATIVE_ROWS(flux, d_cont, q, d, eps, r, gv1)
    CALL EL_TEST_CHECK(MAXVAL(ABS(r)).LT.tol .AND. ABS(gv1 - d(1)).LT.tol, &
      'rows: balanced flux gives zero residual and gV_1 = D_1', nf, rank)
    flux(1) = flux(1) + 1.0d-6
    CALL EL_VA_CONSERVATIVE_ROWS(flux, d_cont, q, d, eps, r, gv1)
    CALL EL_TEST_CHECK(ABS(r(1) - 1.0d-6).LT.tol .AND. ABS(gv1 - (d(1) - 1.0d-6/eps)).LT.tol, &
      'rows: incompatible storage exposed in r_1', nf, rank)

  END SUBROUTINE test_conservative_rows

END PROGRAM TEST_EL_VA_FLUX
