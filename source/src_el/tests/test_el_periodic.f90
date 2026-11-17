PROGRAM TEST_EL_PERIODIC

  ! Unit tests of the periodic image treatment of the particle-mesh coupling
  ! (SimPar@ELPeriodicKernel = wrap; EL_PERIODIC, EL_IMAGE_TRANSFER) on the synthetic structured
  ! mesh of el_test_helpers: n^3 cells on the unit cube, taken as one cell of a lattice with
  ! period 1 on the periodic axes. One rank holds the whole cell, i.e. it touches both periodic
  ! planes of every periodic axis (the self-image case). MPI_Init only for the test harness.
  !
  !   1 image enumeration: interior 1, face 2, edge 4, corner 8 images; a particle exactly on a
  !     plane or outside the cell; a box that touches one plane only; a kernel wider than half
  !     a period; the radius limit
  !   2 halo predicates: rank boxes and particles coupled across the periodic planes
  !   3 default path: with the lattice off the record routines are the single-kernel routines,
  !     bit for bit
  !   4 translation invariance: a particle on / near a periodic plane, edge or corner gives the
  !     fields of the same configuration shifted into the interior (sample, solid volume and
  !     its gradient, reaction force, drag coefficient), whatever the position convention of
  !     the wrapped particle
  !   5 normalisation and conservation: kernel integral as in the interior, deposited volume,
  !     reaction force and drag coefficient conserved; without the lattice the same particle
  !     loses kernel support (the test discriminates)
  !   6 continuity: a particle carried across a periodic plane in small steps changes the
  !     solid-volume field per step exactly as the same motion in the interior does
  !   7 viscous force: kernel-averaged Laplacian of a periodic field, seam against interior;
  !     the boundary-mask check of the transfer pass
  !   8 wide kernel: a kernel wider than half a period (a particle overlapping its own image)
  !     keeps translation invariance and conservation

  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_kernel, el_kernel_width_factor
  USE EL_FIELDS, ONLY: el_field_data, EL_ENSURE_FIELDS, EL_BEGIN_FIELD_UPDATE, EL_FINALIZE
  USE EL_HALO, ONLY: tElParticleRecord
  USE EL_PERIODIC, ONLY: EL_SET_PERIODIC_LATTICE, EL_CLEAR_PERIODIC_LATTICE, &
                         EL_PERIODIC_ACTIVE, EL_PERIODIC_KERNEL_FITS, EL_IMAGE_SHIFTS, &
                         EL_MAX_IMAGES, EL_PERIODIC_SPHERE_REACHES_BOX, &
                         EL_PERIODIC_BOXES_WITHIN_DELTA, EL_PERIODIC_MARKED_ON_PLANES
  USE EL_IMAGE_TRANSFER, ONLY: EL_SAMPLE_RECORD, EL_SAMPLE_RECORD_VISCOUS, EL_DEPOSIT_RECORD
  USE EL_QUADRATURE, ONLY: EL_SAMPLE_SIZE, EL_SAMPLE_NORMALIZATION, EL_INTEGRATE_PARTICLE, &
                           EL_DEPOSIT_PARTICLE
  USE EL_VISCOUS, ONLY: EL_VISCOUS_LAPLACIAN, EL_VIS_SIZE
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, EL_TEST_CHECK, EL_TEST_REPORT, &
                             EL_TEST_RELEASE_MESH, BUILD_STRUCTURED_MESH, &
                             BUILD_DOF_COORDINATES, nel, nvt, net, nat

  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 6                       ! cells per axis, h = 1/6
  REAL*8, PARAMETER :: unit_box(6) = (/ 0.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0 /)
  REAL*8, PARAMETER :: one3(3) = (/ 1.0d0, 1.0d0, 1.0d0 /)
  REAL*8, PARAMETER :: radius_std = 0.04d0          ! kernel radius 0.2 = 1.2 h

  ! the fields of one coupled configuration, gathered on lattice indices for comparison
  TYPE tCoupled
    REAL*8 :: sample(EL_SAMPLE_SIZE)
    REAL*8 :: alpha(0:n-1,0:n-1,0:n-1)
    REAL*8 :: grad_alpha(3,0:n-1,0:n-1,0:n-1)
    ! Q2 dof values on the half-step lattice; the dofs of a periodic plane pair are added
    ! (index 2n folded onto 0 on a periodic axis), as the parallel sum of the solver does
    REAL*8, ALLOCATABLE :: force(:,:,:,:)             ! (3,0:2n,0:2n,0:2n)
    REAL*8, ALLOCATABLE :: drag_b(:,:,:)              ! (0:2n,0:2n,0:2n)
    REAL*8 :: volume
  END TYPE tCoupled

  TYPE(tMultiMesh) :: mesh
  REAL*8, ALLOCATABLE :: coor(:,:), val_u(:), val_v(:), val_w(:), val_p(:), eps(:)
  INTEGER :: ndof, rank, nproc, ierr, failures
  REAL*8 :: pi, h

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  CALL EL_TEST_INIT()
  pi = ACOS(-1.0d0)
  h = 1.0d0/DBLE(n)
  el_kernel = 'deen_poly'
  el_kernel_width_factor = 2.5d0

  CALL BUILD_STRUCTURED_MESH(mesh, n)
  ndof = nvt + net + nat + nel
  ALLOCATE(coor(3,ndof), val_u(ndof), val_v(ndof), val_w(ndof), val_p(4*nel), eps(nel))
  CALL BUILD_DOF_COORDINATES(mesh, coor)
  CALL EL_ENSURE_FIELDS(nel, ndof)

  failures = 0
  CALL test_image_enumeration(failures)
  CALL test_halo_predicates(failures)
  CALL test_default_path(failures)
  CALL test_translation_invariance(failures)
  CALL test_normalisation_and_conservation(failures)
  CALL test_continuity_across_plane(failures)
  CALL test_viscous_average(failures)
  CALL test_wide_kernel(failures)

  CALL EL_CLEAR_PERIODIC_LATTICE()
  CALL EL_TEST_REPORT(failures, rank, nproc, 'EL periodic kernel')
  CALL EL_TEST_RELEASE_MESH(mesh)
  CALL EL_FINALIZE()
  DEALLOCATE(coor, val_u, val_v, val_w, val_p, eps)
  CALL MPI_Finalize(ierr)

  IF (failures.GT.0) ERROR STOP 1

CONTAINS

  !----------------------------------------------------------------------
  ! 1 image enumeration
  !----------------------------------------------------------------------
  SUBROUTINE test_image_enumeration(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: s(3,EL_MAX_IMAGES), r, narrow(6)
    INTEGER :: m
    LOGICAL :: ok

    r = 0.2d0
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3)
    CALL EL_TEST_CHECK(EL_PERIODIC_ACTIVE(), 'images: lattice xy is active', nf, rank)

    CALL EL_IMAGE_SHIFTS((/ 0.5d0, 0.5d0, 0.5d0 /), r, unit_box, s, m)
    CALL EL_TEST_CHECK(m.EQ.1 .AND. ALL(s(:,1).EQ.0.0d0), &
      'images: interior particle is its only image, unshifted', nf, rank)

    CALL EL_IMAGE_SHIFTS((/ 0.95d0, 0.5d0, 0.5d0 /), r, unit_box, s, m)
    ok = (m.EQ.2) .AND. has_shift(s, m, (/ 0.0d0, 0.0d0, 0.0d0 /)) .AND. &
         has_shift(s, m, (/ -1.0d0, 0.0d0, 0.0d0 /))
    CALL EL_TEST_CHECK(ok, 'images: near a periodic face 2 images (0 and -L)', nf, rank)

    CALL EL_IMAGE_SHIFTS((/ 0.95d0, 0.05d0, 0.5d0 /), r, unit_box, s, m)
    ok = (m.EQ.4) .AND. has_shift(s, m, (/ -1.0d0, 1.0d0, 0.0d0 /)) .AND. &
         has_shift(s, m, (/ 0.0d0, 1.0d0, 0.0d0 /)) .AND. has_shift(s, m, (/ -1.0d0, 0.0d0, 0.0d0 /))
    CALL EL_TEST_CHECK(ok, 'images: near a periodic edge 4 images', nf, rank)

    ! the wall axis is never shifted, however close the particle is to the wall
    CALL EL_IMAGE_SHIFTS((/ 0.5d0, 0.5d0, 0.02d0 /), r, unit_box, s, m)
    CALL EL_TEST_CHECK(m.EQ.1 .AND. ALL(s(:,1).EQ.0.0d0), &
      'images: no image across a non-periodic boundary', nf, rank)

    ! exactly on a plane: both sides, from either representation of the plane
    CALL EL_IMAGE_SHIFTS((/ 0.0d0, 0.5d0, 0.5d0 /), r, unit_box, s, m)
    ok = (m.EQ.2) .AND. has_shift(s, m, (/ 0.0d0, 0.0d0, 0.0d0 /)) .AND. &
         has_shift(s, m, (/ 1.0d0, 0.0d0, 0.0d0 /))
    CALL EL_IMAGE_SHIFTS((/ 1.0d0, 0.5d0, 0.5d0 /), r, unit_box, s, m)
    ok = ok .AND. (m.EQ.2) .AND. has_shift(s, m, (/ 0.0d0, 0.0d0, 0.0d0 /)) .AND. &
         has_shift(s, m, (/ -1.0d0, 0.0d0, 0.0d0 /))
    CALL EL_TEST_CHECK(ok, 'images: particle exactly on a plane (x = 0 and x = L)', nf, rank)

    ! a position outside the cell (not yet wrapped by the particle engine)
    CALL EL_IMAGE_SHIFTS((/ -0.03d0, 0.5d0, 0.5d0 /), r, unit_box, s, m)
    ok = (m.EQ.2) .AND. has_shift(s, m, (/ 0.0d0, 0.0d0, 0.0d0 /)) .AND. &
         has_shift(s, m, (/ 1.0d0, 0.0d0, 0.0d0 /))
    CALL EL_IMAGE_SHIFTS((/ 1.7d0, 0.5d0, 0.5d0 /), r, unit_box, s, m)
    ok = ok .AND. (m.EQ.1) .AND. has_shift(s, m, (/ -1.0d0, 0.0d0, 0.0d0 /))
    CALL EL_IMAGE_SHIFTS((/ -2.4d0, 0.5d0, 0.5d0 /), r, unit_box, s, m)
    ok = ok .AND. (m.EQ.1) .AND. has_shift(s, m, (/ 3.0d0, 0.0d0, 0.0d0 /))
    CALL EL_TEST_CHECK(ok, 'images: unwrapped positions outside the cell', nf, rank)

    ! a rank box that touches the plane x = 0 only: the far particle arrives as one image
    narrow = (/ 0.0d0, 0.3d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0 /)
    CALL EL_IMAGE_SHIFTS((/ 0.95d0, 0.5d0, 0.5d0 /), r, narrow, s, m)
    ok = (m.EQ.1) .AND. has_shift(s, m, (/ -1.0d0, 0.0d0, 0.0d0 /))
    CALL EL_IMAGE_SHIFTS((/ 0.6d0, 0.5d0, 0.5d0 /), r, narrow, s, m)
    ok = ok .AND. (m.EQ.0)
    CALL EL_TEST_CHECK(ok, 'images: box on one plane only (1 image, or none out of reach)', &
      nf, rank)

    ! corner of a fully periodic cell
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .TRUE. /), one3)
    CALL EL_IMAGE_SHIFTS((/ 0.95d0, 0.05d0, 0.97d0 /), r, unit_box, s, m)
    ok = (m.EQ.8) .AND. has_shift(s, m, (/ -1.0d0, 1.0d0, -1.0d0 /)) .AND. &
         has_shift(s, m, (/ 0.0d0, 0.0d0, 0.0d0 /))
    CALL EL_TEST_CHECK(ok, 'images: near a periodic corner 8 images', nf, rank)
    ! candidates beyond the kernel radius are not images. 0.12 from three planes: the face
    ! images (gap 0.12) and the edge images (0.12 sqrt 2 = 0.170) reach the cell, the corner
    ! image (0.12 sqrt 3 = 0.208 > 0.2) does not
    CALL EL_IMAGE_SHIFTS((/ 0.88d0, 0.12d0, 0.88d0 /), r, unit_box, s, m)
    CALL EL_TEST_CHECK(m.EQ.7 .AND. .NOT.has_shift(s, m, (/ -1.0d0, 1.0d0, -1.0d0 /)), &
      'images: the corner image out of kernel reach is dropped (7)', nf, rank)
    ! 0.15 from three planes: the edge images (0.212) are out of reach as well
    CALL EL_IMAGE_SHIFTS((/ 0.85d0, 0.15d0, 0.85d0 /), r, unit_box, s, m)
    CALL EL_TEST_CHECK(m.EQ.4 .AND. .NOT.has_shift(s, m, (/ -1.0d0, 1.0d0, 0.0d0 /)), &
      'images: edge and corner images out of kernel reach are dropped (4)', nf, rank)

    ! kernel wider than half a period: the particle and both neighbours on every axis
    CALL EL_IMAGE_SHIFTS((/ 0.5d0, 0.5d0, 0.5d0 /), 0.9d0, unit_box, s, m)
    CALL EL_TEST_CHECK(m.EQ.27, 'images: kernel radius 0.9 L reaches the cell 27 times', &
      nf, rank)
    ok = EL_PERIODIC_KERNEL_FITS(0.999d0) .AND. .NOT.EL_PERIODIC_KERNEL_FITS(1.0d0)
    CALL EL_SET_PERIODIC_LATTICE((/ .FALSE., .TRUE., .FALSE. /), (/ 9.0d9, 0.5d0, 9.0d9 /))
    ok = ok .AND. EL_PERIODIC_KERNEL_FITS(0.49d0) .AND. .NOT.EL_PERIODIC_KERNEL_FITS(0.5d0)
    CALL EL_TEST_CHECK(ok, 'images: kernel radius must be below every period', nf, rank)

    CALL EL_CLEAR_PERIODIC_LATTICE()
    CALL EL_IMAGE_SHIFTS((/ 0.95d0, 0.05d0, 0.5d0 /), r, unit_box, s, m)
    CALL EL_TEST_CHECK(m.EQ.1 .AND. ALL(s(:,1).EQ.0.0d0) .AND. .NOT.EL_PERIODIC_ACTIVE(), &
      'images: no lattice, no images', nf, rank)
  END SUBROUTINE test_image_enumeration

  LOGICAL FUNCTION has_shift(s, m, target)
    REAL*8, INTENT(IN) :: s(3,EL_MAX_IMAGES), target(3)
    INTEGER, INTENT(IN) :: m
    INTEGER :: k
    has_shift = .FALSE.
    DO k = 1, m
      IF (ALL(s(:,k).EQ.target)) has_shift = .TRUE.
    END DO
  END FUNCTION has_shift

  !----------------------------------------------------------------------
  ! 2 halo predicates
  !----------------------------------------------------------------------
  SUBROUTINE test_halo_predicates(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: a(6), b(6), c(6)
    LOGICAL :: ok

    a = (/ 0.0d0, 0.3d0, 0.0d0, 0.5d0, 0.0d0, 1.0d0 /)
    b = (/ 0.7d0, 1.0d0, 0.0d0, 0.5d0, 0.0d0, 1.0d0 /)      ! plain gap 0.4 in x
    c = (/ 0.7d0, 1.0d0, 0.5d0, 1.0d0, 0.0d0, 1.0d0 /)      ! diagonal neighbour across x
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .FALSE., .FALSE. /), one3)
    ok = EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, 0.2d0) .AND. EL_PERIODIC_BOXES_WITHIN_DELTA(b, a, 0.2d0)
    ok = ok .AND. EL_PERIODIC_BOXES_WITHIN_DELTA(a, c, 0.0d0)
    ok = ok .AND. EL_PERIODIC_SPHERE_REACHES_BOX((/ 0.95d0, 0.25d0, 0.5d0 /), 0.1d0, a)
    ok = ok .AND. .NOT.EL_PERIODIC_SPHERE_REACHES_BOX((/ 0.85d0, 0.25d0, 0.5d0 /), 0.1d0, a)
    ! the same pair with the lattice on another axis: not coupled
    CALL EL_SET_PERIODIC_LATTICE((/ .FALSE., .TRUE., .FALSE. /), one3)
    ok = ok .AND. .NOT.EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, 0.2d0)
    ok = ok .AND. EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, 0.4d0)
    ok = ok .AND. .NOT.EL_PERIODIC_SPHERE_REACHES_BOX((/ 0.95d0, 0.25d0, 0.5d0 /), 0.1d0, a)
    ! xy lattice: boxes at opposite corners of the cell are neighbours through both planes
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3)
    a = (/ 0.0d0, 0.3d0, 0.0d0, 0.3d0, 0.0d0, 1.0d0 /)
    c = (/ 0.7d0, 1.0d0, 0.7d0, 1.0d0, 0.0d0, 1.0d0 /)
    ok = ok .AND. EL_PERIODIC_BOXES_WITHIN_DELTA(a, c, 0.05d0)
    CALL EL_TEST_CHECK(ok, 'halo: boxes and particles coupled across periodic planes only', &
      nf, rank)
    CALL EL_CLEAR_PERIODIC_LATTICE()
  END SUBROUTINE test_halo_predicates

  !----------------------------------------------------------------------
  ! 3 default path
  !----------------------------------------------------------------------
  SUBROUTINE test_default_path(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: s0(EL_SAMPLE_SIZE), s1(EL_SAMPLE_SIZE), force(3)
    REAL*8, ALLOCATABLE :: a0(:), f0(:,:), b0(:), g0(:,:)
    LOGICAL :: ok

    CALL EL_CLEAR_PERIODIC_LATTICE()
    CALL set_fields((/ 0, 0, 0 /))
    p = particle((/ 0.95d0, 0.03d0, 0.5d0 /), radius_std)
    force = (/ 0.3d0, -1.1d0, 0.7d0 /)

    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, eps, p, s0)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, s0(1), force, el_field_data, drag_b=2.5d0)
    a0 = el_field_data%alpha_p
    g0 = el_field_data%grad_alpha_p
    f0 = el_field_data%force_rhs
    b0 = el_field_data%drag_B_ml

    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_SAMPLE_RECORD(mesh, 1, val_u, val_v, val_w, val_p, eps, p, unit_box, s1)
    CALL EL_DEPOSIT_RECORD(mesh, 1, p, unit_box, s1(1), force, el_field_data, drag_b=2.5d0)
    ok = ALL(s0.EQ.s1) .AND. ALL(a0.EQ.el_field_data%alpha_p) .AND. &
         ALL(g0.EQ.el_field_data%grad_alpha_p) .AND. ALL(f0.EQ.el_field_data%force_rhs) .AND. &
         ALL(b0.EQ.el_field_data%drag_B_ml)
    CALL EL_TEST_CHECK(ok, 'default: record routines equal the single-kernel routines bitwise', &
      nf, rank)
  END SUBROUTINE test_default_path

  !----------------------------------------------------------------------
  ! 4 translation invariance
  !----------------------------------------------------------------------
  SUBROUTINE test_translation_invariance(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: base(3)

    ! reference position in the interior, off every mesh symmetry
    base = (/ 0.47d0, 0.53d0, 0.51d0 /)

    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3)
    ! 3 cells in x: x = 0.97, the kernel crosses the plane x = L
    CALL check_shift(base, (/ 3, 0, 0 /), radius_std, .FALSE., &
      'invariance xy: face x (particle at 0.97)', nf)
    ! 3 cells in x and y, position left unwrapped (1.03 > L): edge, both conventions
    CALL check_shift(base + (/ 0.06d0, 0.0d0, 0.0d0 /), (/ 3, 3, 0 /), radius_std, .FALSE., &
      'invariance xy: edge, position not wrapped (1.03, 1.03)', nf)
    CALL check_shift(base + (/ 0.06d0, 0.0d0, 0.0d0 /), (/ 3, 3, 0 /), radius_std, .TRUE., &
      'invariance xy: edge, position wrapped into [0,L)', nf)
    CALL check_shift(base, (/ -3, 3, 0 /), radius_std, .TRUE., &
      'invariance xy: edge (-x, +y)', nf)
    ! exactly on the plane x = L and x = 0: interior twin at x = 0.5
    CALL check_shift((/ 0.5d0, 0.53d0, 0.51d0 /), (/ 3, 0, 0 /), radius_std, .FALSE., &
      'invariance xy: particle exactly on the plane x = L', nf)
    CALL check_shift((/ 0.5d0, 0.53d0, 0.51d0 /), (/ -3, 0, 0 /), radius_std, .FALSE., &
      'invariance xy: particle exactly on the plane x = 0', nf)
    ! exactly on the edge x = y = 0
    CALL check_shift((/ 0.5d0, 0.5d0, 0.51d0 /), (/ -3, -3, 0 /), radius_std, .FALSE., &
      'invariance xy: particle exactly on the edge x = y = 0', nf)
    ! next to the wall z = 0 and on the periodic edge: the wall cut is the same in both
    CALL check_shift((/ 0.47d0, 0.53d0, 0.07d0 /), (/ 3, 3, 0 /), radius_std, .TRUE., &
      'invariance xy: periodic edge next to a wall (kernel cut by z = 0)', nf)

    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .TRUE. /), one3)
    CALL check_shift(base, (/ 3, 3, 3 /), radius_std, .TRUE., &
      'invariance xyz: corner (8 images)', nf)
    CALL check_shift(base, (/ 3, -3, 3 /), radius_std, .FALSE., &
      'invariance xyz: corner, position not wrapped', nf)
    CALL check_shift((/ 0.5d0, 0.5d0, 0.5d0 /), (/ 3, 3, 3 /), radius_std, .FALSE., &
      'invariance xyz: particle exactly on the corner', nf)
    CALL EL_CLEAR_PERIODIC_LATTICE()
  END SUBROUTINE test_translation_invariance

  ! Configuration A: particle at pos_ref, fields F(x). Configuration B: particle at
  ! pos_ref + cells*h (wrapped into [0,1) if asked), fields F(x - cells*h). B is A translated
  ! by a whole number of cells along periodic axes, so every field of B is the field of A at
  ! the translated lattice index, to rounding.
  SUBROUTINE check_shift(pos_ref, cells, radius, wrap, label, nf)
    REAL*8, INTENT(IN) :: pos_ref(3), radius
    INTEGER, INTENT(IN) :: cells(3)
    LOGICAL, INTENT(IN) :: wrap
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tCoupled) :: ca, cb
    REAL*8 :: pos(3), force(3), drag_b, err
    INTEGER :: a

    force = (/ 0.3d0, -1.1d0, 0.7d0 /)
    drag_b = 2.5d0
    CALL set_fields((/ 0, 0, 0 /))
    CALL couple(particle(pos_ref, radius), force, drag_b, ca)
    pos = pos_ref + DBLE(cells)*h
    IF (wrap) THEN
      DO a = 1, 3
        IF (cells(a).NE.0) pos(a) = MODULO(pos(a), 1.0d0)
      END DO
    END IF
    CALL set_fields(cells)
    CALL couple(particle(pos, radius), force, drag_b, cb)
    err = coupled_difference(ca, cb, cells)
    IF (rank.EQ.0) WRITE(*,'(A,ES10.2,A,A)') '    max difference ', err, '  ', label
    CALL EL_TEST_CHECK(err.LE.1.0d-11, label, nf, rank)
  END SUBROUTINE check_shift

  !----------------------------------------------------------------------
  ! 5 normalisation and conservation
  !----------------------------------------------------------------------
  SUBROUTINE test_normalisation_and_conservation(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tCoupled) :: c_in, c_seam, c_off
    REAL*8 :: force(3), drag_b, vp, pos_in(3), pos_seam(3), err
    INTEGER :: c

    force = (/ 0.3d0, -1.1d0, 0.7d0 /)
    drag_b = 2.5d0
    vp = 4.0d0*pi*radius_std**3/3.0d0
    pos_in = (/ 0.47d0, 0.53d0, 0.51d0 /)
    pos_seam = (/ 0.97d0, 0.03d0, 0.51d0 /)             ! 3 cells in x and y: periodic edge
    CALL set_fields((/ 0, 0, 0 /))

    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3)
    CALL couple(particle(pos_in, radius_std), force, drag_b, c_in)
    CALL couple(particle(pos_seam, radius_std), force, drag_b, c_seam)
    CALL EL_CLEAR_PERIODIC_LATTICE()
    CALL couple(particle(pos_seam, radius_std), force, drag_b, c_off)

    err = ABS(c_seam%sample(EL_SAMPLE_NORMALIZATION)/c_in%sample(EL_SAMPLE_NORMALIZATION) - 1.0d0)
    CALL EL_TEST_CHECK(err.LE.1.0d-12, &
      'weights: kernel integral on a periodic edge equals the interior one', nf, rank)
    err = c_off%sample(EL_SAMPLE_NORMALIZATION)/c_in%sample(EL_SAMPLE_NORMALIZATION)
    IF (rank.EQ.0) WRITE(*,'(A,F8.5)') &
      '    kernel integral on the edge without the lattice / interior = ', err
    CALL EL_TEST_CHECK(err.LT.0.5d0, &
      'weights: without the lattice the same particle loses kernel support (discriminates)', &
      nf, rank)

    CALL EL_TEST_CHECK(ABS(c_seam%volume/vp - 1.0d0).LE.1.0d-12, &
      'conservation: deposited solid volume = particle volume on the edge', nf, rank)
    err = 0.0d0
    DO c = 1, 3
      err = MAX(err, ABS(SUM(c_seam%force(c,:,:,:)) + force(c)))
    END DO
    CALL EL_TEST_CHECK(err.LE.1.0d-12, &
      'conservation: reaction force sums to -F on the edge', nf, rank)
    CALL EL_TEST_CHECK(ABS(SUM(c_seam%drag_b)/drag_b - 1.0d0).LE.1.0d-12, &
      'conservation: drag coefficient sums to B on the edge', nf, rank)
    ! the volume left of the plane x = L and the one right of x = 0 are both there
    err = SUM(c_seam%alpha(n-1,:,:))/SUM(c_seam%alpha(0,:,:))
    IF (rank.EQ.0) WRITE(*,'(A,F8.5)') &
      '    solid volume in the cell layer below x = L / above x = 0 = ', err
    CALL EL_TEST_CHECK(SUM(c_seam%alpha(0,:,:)).GT.0.0d0 .AND. SUM(c_seam%alpha(n-1,:,:)).GT.0.0d0 &
      .AND. SUM(c_off%alpha(0,:,:)).EQ.0.0d0, &
      'deposit: both sides of the plane receive solid volume (one side without the lattice)', &
      nf, rank)
  END SUBROUTINE test_normalisation_and_conservation

  !----------------------------------------------------------------------
  ! 6 continuity across a periodic plane
  !----------------------------------------------------------------------
  SUBROUTINE test_continuity_across_plane(nf)
    INTEGER, INTENT(INOUT) :: nf
    INTEGER, PARAMETER :: nstep = 12
    TYPE(tCoupled) :: c_prev, c_now
    REAL*8 :: d_seam(nstep), d_in(nstep), d_off(nstep), force(3), x
    INTEGER :: k, pass

    force = 0.0d0
    CALL set_fields((/ 0, 0, 0 /))
    ! the particle moves from x = 0.94 to 1.06 in steps of 0.01 and is wrapped at x = L
    ! (pass 1), the same motion 3 cells further in the interior (pass 2), and pass 1 again
    ! without the lattice (pass 3)
    DO pass = 1, 3
      IF (pass.EQ.3) THEN
        CALL EL_CLEAR_PERIODIC_LATTICE()
      ELSE
        CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3)
      END IF
      DO k = 0, nstep
        x = 0.94d0 + 0.01d0*DBLE(k)
        IF (pass.EQ.2) x = x - 3.0d0*h
        x = MODULO(x, 1.0d0)
        CALL couple(particle((/ x, 0.53d0, 0.51d0 /), radius_std), force, 0.0d0, c_now)
        IF (k.GT.0) THEN
          SELECT CASE (pass)
          CASE (1)
            d_seam(k) = MAXVAL(ABS(c_now%alpha - c_prev%alpha))
          CASE (2)
            d_in(k) = MAXVAL(ABS(c_now%alpha - c_prev%alpha))
          CASE (3)
            d_off(k) = MAXVAL(ABS(c_now%alpha - c_prev%alpha))
          END SELECT
        END IF
        c_prev = c_now
      END DO
    END DO
    CALL EL_CLEAR_PERIODIC_LATTICE()
    IF (rank.EQ.0) THEN
      WRITE(*,'(A,ES11.3,A,ES11.3)') '    max change of alpha_p per step: across the plane ', &
        MAXVAL(d_seam), ', interior ', MAXVAL(d_in)
      WRITE(*,'(A,ES11.3,A,I0)') '    without the lattice               : ', MAXVAL(d_off), &
        ' at step ', MAXLOC(d_off, DIM=1)
    END IF
    CALL EL_TEST_CHECK(MAXVAL(ABS(d_seam - d_in)).LE.1.0d-11*MAXVAL(d_in), &
      'continuity: per-step change of alpha_p across the plane equals the interior one', &
      nf, rank)
    CALL EL_TEST_CHECK(MAXVAL(d_off).GT.10.0d0*MAXVAL(d_in), &
      'continuity: without the lattice the footprint jumps at the crossing (discriminates)', &
      nf, rank)
  END SUBROUTINE test_continuity_across_plane

  !----------------------------------------------------------------------
  ! 7 viscous force
  !----------------------------------------------------------------------
  SUBROUTINE test_viscous_average(nf)
    INTEGER, INTENT(INOUT) :: nf
    LOGICAL, ALLOCATABLE :: wall(:), surface(:), wall_faces(:)
    REAL*8 :: vis(EL_VIS_SIZE), lap_in(3), lap_seam(3), lap_off(3), err
    TYPE(tElParticleRecord) :: p
    INTEGER :: i

    ! xy lattice: the faces on z = 0 and z = 1 are walls, the periodic planes are not
    ALLOCATE(wall(ndof))
    DO i = 1, ndof
      wall(i) = (coor(3,i).LT.1.0d-12 .OR. coor(3,i).GT.1.0d0-1.0d-12)
    END DO
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3)
    ! the check of the transfer pass on a boundary mask: the wall mask holds no face on a
    ! periodic plane, a mask of the whole box surface (what a rank holding both planes of an
    ! axis gets from the fluid) does. Face-centre dofs only, as in the mask of the adapter.
    ALLOCATE(surface(ndof), wall_faces(ndof))
    surface = .FALSE.
    wall_faces = .FALSE.
    DO i = nvt+net+1, nvt+net+nat
      surface(i) = (MINVAL(coor(:,i)).LT.1.0d-12 .OR. MAXVAL(coor(:,i)).GT.1.0d0-1.0d-12)
      wall_faces(i) = wall(i)
    END DO
    CALL EL_TEST_CHECK(EL_PERIODIC_MARKED_ON_PLANES(coor, wall_faces, unit_box).EQ.0 .AND. &
      EL_PERIODIC_MARKED_ON_PLANES(coor, surface, unit_box).EQ.4*n*n, &
      'viscous: boundary-mask check finds the faces on periodic planes (4 n^2), none in the wall mask', &
      nf, rank)
    DEALLOCATE(surface, wall_faces)
    ! kernel radius 2.5 h = 0.417 (the resolution of the viscous unit test); next to z = 0
    p = particle((/ 0.47d0, 0.53d0, 0.6d0*h /), 2.5d0*h/5.0d0)
    CALL set_fields((/ 0, 0, 0 /))
    CALL EL_SAMPLE_RECORD_VISCOUS(mesh, 1, coor, val_u, val_v, val_w, wall, p, unit_box, vis)
    lap_in = EL_VISCOUS_LAPLACIAN(vis)
    p%position = (/ 0.97d0, 0.03d0, 0.6d0*h /)
    CALL set_fields((/ 3, 3, 0 /))
    CALL EL_SAMPLE_RECORD_VISCOUS(mesh, 1, coor, val_u, val_v, val_w, wall, p, unit_box, vis)
    lap_seam = EL_VISCOUS_LAPLACIAN(vis)
    CALL EL_CLEAR_PERIODIC_LATTICE()
    CALL EL_SAMPLE_RECORD_VISCOUS(mesh, 1, coor, val_u, val_v, val_w, wall, p, unit_box, vis)
    lap_off = EL_VISCOUS_LAPLACIAN(vis)

    err = MAXVAL(ABS(lap_seam - lap_in))/MAXVAL(ABS(lap_in))
    IF (rank.EQ.0) THEN
      WRITE(*,'(A,3ES13.5)') '    viscous interior lap      = ', lap_in
      WRITE(*,'(A,3ES13.5,A,ES10.2)') '    viscous periodic edge lap = ', lap_seam, &
        '  rel. difference ', err
      WRITE(*,'(A,3ES13.5)') '    without the lattice       = ', lap_off
    END IF
    CALL EL_TEST_CHECK(err.LE.1.0d-10, &
      'viscous: Laplacian average on a periodic edge equals the interior one', nf, rank)
    CALL EL_TEST_CHECK(MAXVAL(ABS(lap_off - lap_in))/MAXVAL(ABS(lap_in)).GT.1.0d-2, &
      'viscous: without the lattice the edge value differs (discriminates)', nf, rank)
    DEALLOCATE(wall)
  END SUBROUTINE test_viscous_average

  !----------------------------------------------------------------------
  ! 8 kernel wider than half a period
  !----------------------------------------------------------------------
  SUBROUTINE test_wide_kernel(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tCoupled) :: c
    REAL*8 :: radius, vp, force(3), err
    INTEGER :: k

    radius = 0.13d0                                     ! kernel radius 0.65 L
    force = (/ 0.3d0, -1.1d0, 0.7d0 /)
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .TRUE. /), one3)
    CALL check_shift((/ 0.47d0, 0.53d0, 0.51d0 /), (/ 3, 2, -1 /), radius, .TRUE., &
      'wide kernel xyz (radius 0.65 L): translation invariance', nf)
    CALL check_shift((/ 0.47d0, 0.53d0, 0.51d0 /), (/ 1, 0, 0 /), radius, .FALSE., &
      'wide kernel xyz: translation by one cell', nf)
    CALL set_fields((/ 0, 0, 0 /))
    CALL couple(particle((/ 0.97d0, 0.03d0, 0.51d0 /), radius), force, 2.5d0, c)
    vp = 4.0d0*pi*radius**3/3.0d0
    err = ABS(c%volume/vp - 1.0d0)
    DO k = 1, 3
      err = MAX(err, ABS(SUM(c%force(k,:,:,:)) + force(k)))
    END DO
    CALL EL_TEST_CHECK(err.LE.1.0d-12, &
      'wide kernel xyz: volume and reaction force conserved', nf, rank)
    CALL EL_CLEAR_PERIODIC_LATTICE()
  END SUBROUTINE test_wide_kernel

  !----------------------------------------------------------------------
  ! Helpers
  !----------------------------------------------------------------------
  FUNCTION particle(position, radius) RESULT(p)
    REAL*8, INTENT(IN) :: position(3), radius
    TYPE(tElParticleRecord) :: p
    p%position = position
    p%radius = radius
    p%density = 1.0d0
    p%velocity = (/ 0.1d0, -0.2d0, 0.3d0 /)
  END FUNCTION particle

  ! Smooth fields of period 1 on every axis, evaluated at x - cells*h: velocity on the Q2
  ! dofs, pressure (value and gradient at the cell centre, the P1 layout of the sample
  ! routine) and porosity per cell.
  SUBROUTINE set_fields(cells)
    INTEGER, INTENT(IN) :: cells(3)
    REAL*8 :: x(3), tp
    INTEGER :: i, iel, jj

    tp = 2.0d0*pi
    DO i = 1, ndof
      x = coor(:,i) - DBLE(cells)*h
      val_u(i) = SIN(tp*x(1))*COS(tp*x(2)) + 0.5d0*SIN(tp*x(3))
      val_v(i) = COS(tp*x(1)) + SIN(tp*x(2))*SIN(tp*x(3))
      val_w(i) = 1.0d0 + 0.3d0*COS(tp*(x(1)+x(2))) + 0.2d0*SIN(tp*x(3))
    END DO
    DO iel = 1, nel
      x = coor(:, nvt+net+nat+iel) - DBLE(cells)*h
      jj = 4*(iel-1) + 1
      val_p(jj) = COS(tp*x(1))*SIN(tp*x(2)) + 0.1d0*COS(tp*x(3))
      val_p(jj+1) = -tp*SIN(tp*x(1))*SIN(tp*x(2))
      val_p(jj+2) = tp*COS(tp*x(1))*COS(tp*x(2))
      val_p(jj+3) = -0.1d0*tp*SIN(tp*x(3))
      eps(iel) = 0.8d0 + 0.1d0*SIN(tp*x(1))*COS(tp*x(3)) + 0.05d0*COS(tp*x(2))
    END DO
  END SUBROUTINE set_fields

  ! One transfer pass for one particle on the single rank of the test: sample, then deposit
  ! with the sampled normalisation (what the owner sends back), gathered on lattice indices.
  SUBROUTINE couple(p, force, drag_b, c)
    TYPE(tElParticleRecord), INTENT(IN) :: p
    REAL*8, INTENT(IN) :: force(3), drag_b
    TYPE(tCoupled), INTENT(OUT) :: c
    LOGICAL :: axis(3)
    INTEGER :: i, iel, ix(3), a

    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_SAMPLE_RECORD(mesh, 1, val_u, val_v, val_w, val_p, eps, p, unit_box, c%sample)
    CALL EL_DEPOSIT_RECORD(mesh, 1, p, unit_box, c%sample(EL_SAMPLE_NORMALIZATION), force, &
                           el_field_data, drag_b=drag_b)
    c%volume = SUM(el_field_data%alpha_p*mesh%level(1)%dvol)
    DO iel = 1, nel
      ix(1) = MOD(iel-1, n)
      ix(2) = MOD((iel-1)/n, n)
      ix(3) = (iel-1)/(n*n)
      c%alpha(ix(1),ix(2),ix(3)) = el_field_data%alpha_p(iel)
      c%grad_alpha(:,ix(1),ix(2),ix(3)) = el_field_data%grad_alpha_p(:,iel)
    END DO
    axis = lattice_axes()
    ALLOCATE(c%force(3,0:2*n,0:2*n,0:2*n), c%drag_b(0:2*n,0:2*n,0:2*n))
    c%force = 0.0d0
    c%drag_b = 0.0d0
    DO i = 1, ndof
      ix = NINT(2.0d0*DBLE(n)*coor(:,i))
      DO a = 1, 3
        IF (axis(a) .AND. ix(a).EQ.2*n) ix(a) = 0
      END DO
      c%force(:,ix(1),ix(2),ix(3)) = c%force(:,ix(1),ix(2),ix(3)) + el_field_data%force_rhs(:,i)
      c%drag_b(ix(1),ix(2),ix(3)) = c%drag_b(ix(1),ix(2),ix(3)) + el_field_data%drag_B_ml(i)
    END DO
  END SUBROUTINE couple

  FUNCTION lattice_axes() RESULT(axis)
    USE EL_PERIODIC, ONLY: EL_PERIODIC_AXES
    LOGICAL :: axis(3)
    axis = EL_PERIODIC_AXES()
  END FUNCTION lattice_axes

  ! Largest difference between configuration a and configuration b read at the lattice index
  ! translated by `cells` (b is a moved by that many cells; periodic axes only).
  REAL*8 FUNCTION coupled_difference(a, b, cells) RESULT(err)
    TYPE(tCoupled), INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: cells(3)
    INTEGER :: i, j, k, ib, jb, kb, m(3)

    err = MAXVAL(ABS(a%sample - b%sample))
    err = MAX(err, ABS(a%volume - b%volume)/a%volume)
    DO k = 0, n-1
      DO j = 0, n-1
        DO i = 0, n-1
          ib = MODULO(i + cells(1), n)
          jb = MODULO(j + cells(2), n)
          kb = MODULO(k + cells(3), n)
          err = MAX(err, ABS(a%alpha(i,j,k) - b%alpha(ib,jb,kb)))
          ! the gradient is of order alpha/h: compared relative to that scale
          err = MAX(err, MAXVAL(ABS(a%grad_alpha(:,i,j,k) - b%grad_alpha(:,ib,jb,kb)))*h)
        END DO
      END DO
    END DO
    ! dof lattice: a periodic axis has the indices 0..2n-1 (2n folded onto 0), a wall axis
    ! 0..2n and no translation
    m = 2*n
    WHERE (cells.EQ.0) m = 2*n + 1
    DO k = 0, m(3)-1
      DO j = 0, m(2)-1
        DO i = 0, m(1)-1
          ib = MODULO(i + 2*cells(1), m(1))
          jb = MODULO(j + 2*cells(2), m(2))
          kb = MODULO(k + 2*cells(3), m(3))
          err = MAX(err, MAXVAL(ABS(a%force(:,i,j,k) - b%force(:,ib,jb,kb))))
          err = MAX(err, ABS(a%drag_b(i,j,k) - b%drag_b(ib,jb,kb)))
        END DO
      END DO
    END DO
  END FUNCTION coupled_difference

END PROGRAM TEST_EL_PERIODIC
