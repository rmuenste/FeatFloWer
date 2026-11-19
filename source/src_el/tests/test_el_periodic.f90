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
  !   9 neighbour test: f(a,b) and f(b,a) agree bit for bit on pairs whose gap equals the
  !     kernel radius (the ties on which two ranks could disagree), named cases and a scan
  !  10 wall lift (saffman_mei_wall): a periodic plane is not a wall; the lift does not change
  !     with the position relative to a periodic plane, and a fully periodic box gives the
  !     unbounded form
  !  11 refusals and robustness, at predicate level: invalid period, lattice against domain,
  !     rank box spanning a period, mask and coordinate sizes, positions far outside the cell
  !     and positions that are not finite
  ! The distributed path (graph build, record exchange, owner reduction on several ranks) is
  ! in test_el_periodic_mpi.f90.

  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_kernel, el_kernel_width_factor
  USE EL_FIELDS, ONLY: el_field_data, EL_ENSURE_FIELDS, EL_BEGIN_FIELD_UPDATE, EL_FINALIZE
  USE EL_HALO, ONLY: tElParticleRecord
  USE EL_PERIODIC, ONLY: EL_SET_PERIODIC_LATTICE, EL_CLEAR_PERIODIC_LATTICE, &
                         EL_PERIODIC_ACTIVE, EL_PERIODIC_KERNEL_FITS, EL_IMAGE_SHIFTS, &
                         EL_MAX_IMAGES, EL_PERIODIC_SPHERE_REACHES_BOX, &
                         EL_PERIODIC_BOXES_WITHIN_DELTA, EL_PERIODIC_MARKED_ON_PLANES, &
                         EL_PERIODIC_LATTICE_FITS_DOMAIN, EL_PERIODIC_BOX_SPANS_AXIS
  USE EL_VISCOUS, ONLY: EL_INTEGRATE_PARTICLE_VISCOUS
  USE EL_CONFIG, ONLY: el_lift_model, EL_KERNEL_SUPPORT_RADIUS
  USE EL_FORCES, ONLY: EL_LIFT_CLOSURE
  USE EL_GEOMETRY, ONLY: EL_SET_DOMAIN_BOX, el_domain_set, el_domain_box_set
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
  LOGICAL :: lattice_ok

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
  CALL test_neighbour_symmetry(failures)
  CALL test_wall_lift(failures)
  CALL test_refusals_and_robustness(failures)

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
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3, lattice_ok)
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
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .TRUE. /), one3, lattice_ok)
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
    CALL EL_SET_PERIODIC_LATTICE((/ .FALSE., .TRUE., .FALSE. /), (/ 9.0d9, 0.5d0, 9.0d9 /), lattice_ok)
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
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .FALSE., .FALSE. /), one3, lattice_ok)
    ok = EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, 0.2d0) .AND. EL_PERIODIC_BOXES_WITHIN_DELTA(b, a, 0.2d0)
    ok = ok .AND. EL_PERIODIC_BOXES_WITHIN_DELTA(a, c, 0.0d0)
    ok = ok .AND. EL_PERIODIC_SPHERE_REACHES_BOX((/ 0.95d0, 0.25d0, 0.5d0 /), 0.1d0, a)
    ok = ok .AND. .NOT.EL_PERIODIC_SPHERE_REACHES_BOX((/ 0.85d0, 0.25d0, 0.5d0 /), 0.1d0, a)
    ! the same pair with the lattice on another axis: not coupled
    CALL EL_SET_PERIODIC_LATTICE((/ .FALSE., .TRUE., .FALSE. /), one3, lattice_ok)
    ok = ok .AND. .NOT.EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, 0.2d0)
    ok = ok .AND. EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, 0.4d0)
    ok = ok .AND. .NOT.EL_PERIODIC_SPHERE_REACHES_BOX((/ 0.95d0, 0.25d0, 0.5d0 /), 0.1d0, a)
    ! xy lattice: boxes at opposite corners of the cell are neighbours through both planes
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3, lattice_ok)
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
    REAL*8 :: v0(EL_VIS_SIZE), v1(EL_VIS_SIZE)
    REAL*8, ALLOCATABLE :: a0(:), f0(:,:), b0(:), g0(:,:)
    LOGICAL, ALLOCATABLE :: wall(:)
    LOGICAL :: ok
    INTEGER :: i

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
    ! bit patterns, not floating-point equality (which would take -0 for +0)
    ok = same_bits(s0, s1) .AND. same_bits(a0, el_field_data%alpha_p) .AND. &
         same_bits(RESHAPE(g0, (/ SIZE(g0) /)), &
                   RESHAPE(el_field_data%grad_alpha_p, (/ SIZE(g0) /))) .AND. &
         same_bits(RESHAPE(f0, (/ SIZE(f0) /)), &
                   RESHAPE(el_field_data%force_rhs, (/ SIZE(f0) /))) .AND. &
         same_bits(b0, el_field_data%drag_B_ml)
    CALL EL_TEST_CHECK(ok, 'default: sample and deposit of a record have the bit patterns of the single-kernel routines', &
      nf, rank)

    ! the viscous record against the single-kernel routine, next to a wall and a box face
    ALLOCATE(wall(ndof))
    DO i = 1, ndof
      wall(i) = (coor(3,i).LT.1.0d-12 .OR. coor(3,i).GT.1.0d0-1.0d-12)
    END DO
    p = particle((/ 0.97d0, 0.03d0, 0.6d0*h /), 2.5d0*h/5.0d0)
    CALL EL_INTEGRATE_PARTICLE_VISCOUS(mesh, 1, coor, val_u, val_v, val_w, wall, p%position, &
      EL_KERNEL_SUPPORT_RADIUS(p%radius), v0)
    CALL EL_SAMPLE_RECORD_VISCOUS(mesh, 1, coor, val_u, val_v, val_w, wall, p, unit_box, v1)
    CALL EL_TEST_CHECK(same_bits(v0, v1) .AND. MAXVAL(ABS(v0)).GT.0.0d0, &
      'default: viscous record has the bit pattern of the single-kernel routine', nf, rank)
    DEALLOCATE(wall)
  END SUBROUTINE test_default_path

  !----------------------------------------------------------------------
  ! 4 translation invariance
  !----------------------------------------------------------------------
  SUBROUTINE test_translation_invariance(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: base(3)

    ! reference position in the interior, off every mesh symmetry
    base = (/ 0.47d0, 0.53d0, 0.51d0 /)

    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3, lattice_ok)
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

    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .TRUE. /), one3, lattice_ok)
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

    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3, lattice_ok)
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
        CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3, lattice_ok)
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
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3, lattice_ok)
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
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .TRUE. /), one3, lattice_ok)
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
  ! 9 neighbour test: bitwise symmetry on ties
  !----------------------------------------------------------------------
  SUBROUTINE test_neighbour_symmetry(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: a(6), b(6), slab(6,8), dh, delta, len(3), lo, w, d(4)
    INTEGER :: i, j, k, it, c, seed, n_asym, n_old_asym, n_checked
    LOGICAL :: ok

    ! review case 1: L = 1, a = [0.1,0.2], b = [0.8,0.9] in x, delta = 0.2 (the gap across
    ! the plane is 0.2 up to rounding; the previous expression gave different verdicts)
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .FALSE., .FALSE. /), one3, lattice_ok)
    a = (/ 0.1d0, 0.2d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0 /)
    b = (/ 0.8d0, 0.9d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0 /)
    ok = (EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, 0.2d0) .EQV. EL_PERIODIC_BOXES_WITHIN_DELTA(b, a, 0.2d0))
    CALL EL_TEST_CHECK(ok .AND. (old_within(a, b, 0.2d0, one3) .NEQV. old_within(b, a, 0.2d0, one3)), &
      'neighbours: L = 1, [0.1,0.2] and [0.8,0.9], delta 0.2: same verdict both ways (the previous expression disagreed)', &
      nf, rank)

    ! review case 2: L = 0.02, eight slabs, every pair, delta = 0.005 written two ways
    len = (/ 0.02d0, 1.0d0, 1.0d0 /)
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .FALSE., .FALSE. /), len, lattice_ok)
    dh = 0.02d0/8.0d0
    DO k = 1, 8
      slab(:,k) = (/ DBLE(k-1)*dh, DBLE(k)*dh, 0.0d0, 1.0d0, 0.0d0, 1.0d0 /)
    END DO
    n_asym = 0
    n_old_asym = 0
    DO it = 1, 2
      delta = 0.005d0
      IF (it.EQ.2) delta = 2.0d0*2.5d0*0.001d0
      DO i = 1, 8
        DO j = 1, 8
          IF (EL_PERIODIC_BOXES_WITHIN_DELTA(slab(:,i), slab(:,j), delta) .NEQV. &
              EL_PERIODIC_BOXES_WITHIN_DELTA(slab(:,j), slab(:,i), delta)) n_asym = n_asym + 1
          IF (old_within(slab(:,i), slab(:,j), delta, len) .NEQV. &
              old_within(slab(:,j), slab(:,i), delta, len)) n_old_asym = n_old_asym + 1
        END DO
      END DO
    END DO
    IF (rank.EQ.0) WRITE(*,'(A,I0,A,I0)') '    eight slabs of L = 0.02, delta = 0.005: asymmetric pairs ', &
      n_asym, ', with the previous expression ', n_old_asym
    CALL EL_TEST_CHECK(n_asym.EQ.0, &
      'neighbours: eight slabs of L = 0.02, tie distance 0.005: every pair symmetric', nf, rank)
    ! the pair named in the review, [0.005,0.0075] and [0.0175,0.02]: coupled across the plane
    CALL EL_TEST_CHECK(EL_PERIODIC_BOXES_WITHIN_DELTA(slab(:,3), slab(:,8), 0.0051d0) .AND. &
      .NOT.EL_PERIODIC_BOXES_WITHIN_DELTA(slab(:,3), slab(:,8), 0.0049d0), &
      'neighbours: slabs 3 and 8 of L = 0.02 are 0.005 apart across the plane', nf, rank)

    ! scan: pseudo-random box pairs on three lattices; delta is put on the tie (the gap as
    ! the previous expression computed it in either order, and its two neighbours in the
    ! floating-point line)
    n_asym = 0
    n_old_asym = 0
    n_checked = 0
    seed = 12345
    DO it = 1, 3
      len = (/ 1.0d0, 0.02d0, 0.3d0 /)
      len = CSHIFT(len, it-1)
      CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., it.EQ.3 /), len, lattice_ok)
      DO k = 1, 4000
        DO c = 1, 3
          lo = lcg(seed)*0.9d0*len(c)
          w = (0.02d0 + 0.08d0*lcg(seed))*len(c)
          a(2*c-1) = lo
          a(2*c) = MIN(lo + w, len(c))
          lo = lcg(seed)*0.9d0*len(c)
          w = (0.02d0 + 0.08d0*lcg(seed))*len(c)
          b(2*c-1) = lo
          b(2*c) = MIN(lo + w, len(c))
        END DO
        ! two of three pairs differ along x only, where a tie is a single rounding
        IF (MOD(k,3).NE.0) THEN
          b(3:6) = a(3:6)
        END IF
        d(1) = old_gap(a, b, len)
        d(2) = old_gap(b, a, len)
        d(3) = NEAREST(d(1), 1.0d0)
        d(4) = NEAREST(d(2), -1.0d0)
        DO c = 1, 4
          n_checked = n_checked + 1
          IF (EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, d(c)) .NEQV. &
              EL_PERIODIC_BOXES_WITHIN_DELTA(b, a, d(c))) n_asym = n_asym + 1
          IF (old_within(a, b, d(c), len) .NEQV. old_within(b, a, d(c), len)) &
            n_old_asym = n_old_asym + 1
        END DO
      END DO
    END DO
    IF (rank.EQ.0) WRITE(*,'(A,I0,A,I0,A,I0)') '    scan: ', n_checked, &
      ' ties, asymmetric verdicts ', n_asym, ', with the previous expression ', n_old_asym
    CALL EL_TEST_CHECK(n_asym.EQ.0, 'neighbours: scan of tie distances, every pair symmetric', &
      nf, rank)
    CALL EL_TEST_CHECK(n_old_asym.GT.0, &
      'neighbours: the scan contains ties on which the previous expression disagreed (discriminates)', &
      nf, rank)
    CALL EL_CLEAR_PERIODIC_LATTICE()
  END SUBROUTINE test_neighbour_symmetry

  ! the neighbour test as it was before the review (box b shifted, then subtracted)
  REAL*8 FUNCTION old_gap(a, b, len) RESULT(dist)
    USE EL_PERIODIC, ONLY: EL_PERIODIC_AXES
    REAL*8, INTENT(IN) :: a(6), b(6), len(3)
    REAL*8 :: gap(3), sh
    LOGICAL :: axis(3)
    INTEGER :: k, m
    axis = EL_PERIODIC_AXES()
    DO k = 1, 3
      gap(k) = MAX(0.0d0, MAX(a(2*k-1)-b(2*k), b(2*k-1)-a(2*k)))
      IF (.NOT.axis(k)) CYCLE
      DO m = -1, 1, 2
        sh = DBLE(m)*len(k)
        gap(k) = MIN(gap(k), MAX(0.0d0, MAX(a(2*k-1)-(b(2*k)+sh), (b(2*k-1)+sh)-a(2*k))))
      END DO
    END DO
    dist = SQRT(SUM(gap*gap))
  END FUNCTION old_gap

  LOGICAL FUNCTION old_within(a, b, delta, len)
    REAL*8, INTENT(IN) :: a(6), b(6), delta, len(3)
    old_within = old_gap(a, b, len).LE.delta
  END FUNCTION old_within

  ! deterministic pseudo-random number in [0,1)
  REAL*8 FUNCTION lcg(seed)
    INTEGER, INTENT(INOUT) :: seed
    seed = MOD(1103515245_8*INT(seed,8) + 12345_8, 2147483648_8)
    lcg = DBLE(seed)/2147483648.0d0
  END FUNCTION lcg

  !----------------------------------------------------------------------
  ! 10 wall lift: a periodic plane is not a wall
  !----------------------------------------------------------------------
  SUBROUTINE test_wall_lift(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: state(8), u(3), grad_u(3,3), lift_a(3), lift_b(3), lift_c(3), lift_off(3)
    REAL*8 :: lift_free(3), lift_xyz(3)
    CHARACTER(LEN=32) :: saved_model
    REAL*8, PARAMETER :: rho = 1.0d0, nu = 0.01d0

    saved_model = el_lift_model
    CALL EL_SET_DOMAIN_BOX(0.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0)
    ! slip 1, diameter 0.04: Re = 4, the wall-bounded branch of the closure
    state = 0.0d0
    state(7) = 0.02d0
    state(8) = 2.0d0
    u = (/ 1.0d0, 0.0d0, 0.0d0 /)
    grad_u = 0.0d0
    grad_u(1,3) = 3.0d0
    grad_u(2,1) = -1.0d0

    el_lift_model = 'saffman_mei_wall'
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3, lattice_ok)
    state(1:3) = (/ 0.02d0, 0.5d0, 0.3d0 /)          ! 0.02 from the periodic plane x = 0
    CALL EL_LIFT_CLOSURE(state, u, grad_u, rho, nu, 1.0d0, lift_a)
    state(1:3) = (/ 0.52d0, 0.5d0, 0.3d0 /)          ! the same, moved by half a period
    CALL EL_LIFT_CLOSURE(state, u, grad_u, rho, nu, 1.0d0, lift_b)
    state(1:3) = (/ 0.98d0, 0.97d0, 0.3d0 /)         ! on the other side of the plane, near y = L
    CALL EL_LIFT_CLOSURE(state, u, grad_u, rho, nu, 1.0d0, lift_c)
    CALL EL_TEST_CHECK(same_bits(lift_a, lift_b) .AND. same_bits(lift_a, lift_c) .AND. &
      lift_a(3).NE.0.0d0 .AND. lift_a(1).EQ.0.0d0 .AND. lift_a(2).EQ.0.0d0, &
      'wall lift xy: the same on either side of a periodic plane, directed along the wall normal z', &
      nf, rank)

    CALL EL_CLEAR_PERIODIC_LATTICE()
    state(1:3) = (/ 0.02d0, 0.5d0, 0.3d0 /)
    CALL EL_LIFT_CLOSURE(state, u, grad_u, rho, nu, 1.0d0, lift_off)
    CALL EL_TEST_CHECK(lift_off(1).NE.0.0d0 .AND. lift_off(3).EQ.0.0d0, &
      'wall lift: without the lattice the plane x = 0 acts as the wall (discriminates)', nf, rank)

    ! every face periodic: no wall, the unbounded form
    el_lift_model = 'saffman_mei'
    CALL EL_LIFT_CLOSURE(state, u, grad_u, rho, nu, 1.0d0, lift_free)
    el_lift_model = 'saffman_mei_wall'
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .TRUE. /), one3, lattice_ok)
    CALL EL_LIFT_CLOSURE(state, u, grad_u, rho, nu, 1.0d0, lift_xyz)
    CALL EL_TEST_CHECK(same_bits(lift_xyz, lift_free) .AND. MAXVAL(ABS(lift_free)).GT.0.0d0, &
      'wall lift xyz: no wall left, the unbounded saffman_mei lift', nf, rank)

    CALL EL_CLEAR_PERIODIC_LATTICE()
    el_lift_model = saved_model
    el_domain_set = .FALSE.
    el_domain_box_set = .FALSE.
  END SUBROUTINE test_wall_lift

  !----------------------------------------------------------------------
  ! 11 refusals (predicate level) and robustness of the image enumeration
  !----------------------------------------------------------------------
  SUBROUTINE test_refusals_and_robustness(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: s(3,EL_MAX_IMAGES), box(6), big, nan, inf, zero
    LOGICAL :: ok, ok2
    LOGICAL, ALLOCATABLE :: mask(:)
    INTEGER :: m, bad

    ! a period that is not a positive finite number sets no lattice
    zero = 0.0d0
    nan = zero/zero
    inf = HUGE(1.0d0)
    inf = inf*2.0d0
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .FALSE., .FALSE. /), (/ 0.0d0, 1.0d0, 1.0d0 /), ok)
    ok2 = .NOT.ok .AND. .NOT.EL_PERIODIC_ACTIVE()
    CALL EL_SET_PERIODIC_LATTICE((/ .FALSE., .TRUE., .FALSE. /), (/ 1.0d0, -2.0d0, 1.0d0 /), ok)
    ok2 = ok2 .AND. .NOT.ok .AND. .NOT.EL_PERIODIC_ACTIVE()
    CALL EL_SET_PERIODIC_LATTICE((/ .FALSE., .FALSE., .TRUE. /), (/ 1.0d0, 1.0d0, nan /), ok)
    ok2 = ok2 .AND. .NOT.ok .AND. .NOT.EL_PERIODIC_ACTIVE()
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .FALSE., .FALSE. /), (/ inf, 1.0d0, 1.0d0 /), ok)
    ok2 = ok2 .AND. .NOT.ok .AND. .NOT.EL_PERIODIC_ACTIVE()
    ! the length of an axis that is not periodic is not read
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .FALSE., .FALSE. /), (/ 1.0d0, -1.0d0, nan /), ok)
    ok2 = ok2 .AND. ok .AND. EL_PERIODIC_ACTIVE()
    CALL EL_TEST_CHECK(ok2, 'refusal: a period that is zero, negative or not finite sets no lattice', &
      nf, rank)

    ! lattice against the fluid domain (start-up check 913)
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), one3, lattice_ok)
    ok = EL_PERIODIC_LATTICE_FITS_DOMAIN(unit_box, bad) .AND. bad.EQ.0
    ok = ok .AND. EL_PERIODIC_LATTICE_FITS_DOMAIN((/ 2.0d0, 3.0d0, -0.5d0, 0.5d0, 0.0d0, 7.0d0 /), bad)
    ok = ok .AND. EL_PERIODIC_LATTICE_FITS_DOMAIN((/ 0.0d0, 1.0d0+5.0d-9, 0.0d0, 1.0d0, 0.0d0, 1.0d0 /), bad)
    ok = ok .AND. .NOT.EL_PERIODIC_LATTICE_FITS_DOMAIN((/ 0.0d0, 1.0d0+5.0d-8, 0.0d0, 1.0d0, 0.0d0, 1.0d0 /), bad)
    ok = ok .AND. bad.EQ.1
    ! the batch-settling box with the default period 1 (the deck without PeriodicLength)
    ok = ok .AND. .NOT.EL_PERIODIC_LATTICE_FITS_DOMAIN((/ 0.0d0, 0.02d0, 0.0d0, 0.02d0, 0.0d0, 0.2d0 /), bad)
    ok = ok .AND. .NOT.EL_PERIODIC_LATTICE_FITS_DOMAIN((/ 0.0d0, 1.0d0, 0.0d0, 0.5d0, 0.0d0, 1.0d0 /), bad)
    ok = ok .AND. bad.EQ.2
    ok = ok .AND. .NOT.EL_PERIODIC_LATTICE_FITS_DOMAIN((/ 0.0d0, nan, 0.0d0, 1.0d0, 0.0d0, 1.0d0 /), bad)
    CALL EL_TEST_CHECK(ok, 'refusal: the lattice must be the fluid domain (1e-8 periods)', nf, rank)

    ! a rank box that spans a period (halo check 914)
    ok = EL_PERIODIC_BOX_SPANS_AXIS(unit_box)
    ok = ok .AND. .NOT.EL_PERIODIC_BOX_SPANS_AXIS((/ 0.0d0, 0.5d0, 0.5d0, 1.0d0, 0.0d0, 1.0d0 /))
    ok = ok .AND. EL_PERIODIC_BOX_SPANS_AXIS((/ 0.0d0, 0.5d0, 0.0d0, 1.0d0-1.0d-9, 0.0d0, 0.2d0 /))
    ok = ok .AND. .NOT.EL_PERIODIC_BOX_SPANS_AXIS((/ 0.0d0, 0.5d0, 0.0d0, 1.0d0-1.0d-7, 0.0d0, 1.0d0 /))
    ok = ok .AND. EL_PERIODIC_BOX_SPANS_AXIS((/ 0.0d0, nan, 0.0d0, 0.5d0, 0.0d0, 1.0d0 /))
    CALL EL_TEST_CHECK(ok, 'refusal: a rank box that spans a period is recognised (z, not periodic, may)', &
      nf, rank)

    ! boundary mask longer than the coordinates (check 912)
    ALLOCATE(mask(ndof+1))
    mask = .FALSE.
    ok = EL_PERIODIC_MARKED_ON_PLANES(coor, mask, unit_box).EQ.-1
    ok = ok .AND. EL_PERIODIC_MARKED_ON_PLANES(coor, mask(1:ndof), unit_box).EQ.0
    ok = ok .AND. EL_PERIODIC_MARKED_ON_PLANES(coor, mask(1:nvt), unit_box).EQ.0
    DEALLOCATE(mask)
    CALL EL_TEST_CHECK(ok, 'refusal: a boundary mask that the coordinates do not cover is reported', &
      nf, rank)

    ! image enumeration far outside the cell: no integer conversion, never more than 27
    box = (/ 0.0d0, 0.3d0, 0.0d0, 0.4d0, 0.0d0, 1.0d0 /)
    CALL EL_IMAGE_SHIFTS((/ 3.0d9 + 0.25d0, 0.2d0, 0.5d0 /), 0.2d0, box, s, m)
    ok = (m.EQ.1)
    IF (ok) ok = (s(1,1).EQ.-3.0d9 .AND. s(2,1).EQ.0.0d0)
    CALL EL_IMAGE_SHIFTS((/ -7.0d12 + 0.95d0, 0.2d0, 0.5d0 /), 0.2d0, box, s, m)
    ok = ok .AND. (m.EQ.1)
    IF (ok) ok = (s(1,1).EQ.7.0d12 - 1.0d0)
    big = 1.0d300
    CALL EL_TEST_CHECK(ok, 'images: positions 3e9 and 7e12 periods outside the cell', nf, rank)
    ! beyond 2**52 periods the lattice index is not exact: no image, not three copies of one
    CALL EL_IMAGE_SHIFTS((/ big, -big, 0.5d0 /), 0.2d0, box, s, m)
    ok = (m.EQ.0)
    CALL EL_IMAGE_SHIFTS((/ 0.1d0, HUGE(1.0d0), 0.5d0 /), 0.2d0, box, s, m)
    ok = ok .AND. (m.EQ.0)
    CALL EL_TEST_CHECK(ok, 'images: a position beyond 2**52 periods has no image', nf, rank)
    CALL EL_IMAGE_SHIFTS((/ nan, 0.2d0, 0.5d0 /), 0.2d0, box, s, m)
    ok = (m.EQ.0)
    CALL EL_IMAGE_SHIFTS((/ 0.1d0, inf, 0.5d0 /), 0.2d0, box, s, m)
    ok = ok .AND. (m.EQ.0)
    CALL EL_IMAGE_SHIFTS((/ 0.1d0, 0.2d0, -inf /), 0.2d0, box, s, m)
    ok = ok .AND. (m.EQ.0)
    ok = ok .AND. .NOT.EL_PERIODIC_SPHERE_REACHES_BOX((/ nan, 0.2d0, 0.5d0 /), 0.2d0, box)
    CALL EL_TEST_CHECK(ok, 'images: a position that is not finite has no image', nf, rank)
    ! the bound of 27 does not depend on the preconditions: a radius of ten periods
    CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .TRUE. /), one3, lattice_ok)
    CALL EL_IMAGE_SHIFTS((/ 0.5d0, 0.5d0, 0.5d0 /), 10.0d0, unit_box, s, m)
    CALL EL_TEST_CHECK(m.EQ.EL_MAX_IMAGES, &
      'images: at most 27 whatever the radius (a radius beyond a period is refused by the halo)', &
      nf, rank)
    CALL EL_CLEAR_PERIODIC_LATTICE()
  END SUBROUTINE test_refusals_and_robustness

  !----------------------------------------------------------------------
  ! Helpers
  !----------------------------------------------------------------------
  ! the two arrays have the same bit patterns (and the same size)
  LOGICAL FUNCTION same_bits(a, b)
    REAL*8, INTENT(IN) :: a(:), b(:)
    same_bits = .FALSE.
    IF (SIZE(a).NE.SIZE(b)) RETURN
    same_bits = ALL(TRANSFER(a, (/ 0_8 /)).EQ.TRANSFER(b, (/ 0_8 /)))
  END FUNCTION same_bits

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
