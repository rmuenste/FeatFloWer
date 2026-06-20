PROGRAM TEST_EL_CONVERGENCE

  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_kernel, el_kernel_width_factor
  USE EL_FIELDS, ONLY: el_field_data, EL_ENSURE_FIELDS, &
                       EL_BEGIN_FIELD_UPDATE, EL_FINALIZE
  USE EL_HALO, ONLY: tElParticleRecord
  USE EL_QUADRATURE, ONLY: EL_SAMPLE_SIZE, EL_SAMPLE_NORMALIZATION, &
                           EL_SAMPLE_UF_BEGIN, EL_SAMPLE_UF_END, &
                           EL_SAMPLE_GRAD_U_BEGIN, EL_SAMPLE_GRAD_U_END, &
                           EL_SAMPLE_PRESSURE, EL_SAMPLE_GRAD_P_BEGIN, &
                           EL_SAMPLE_GRAD_P_END, EL_INTEGRATE_PARTICLE, &
                           EL_DEPOSIT_PARTICLE
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, EL_TEST_CHECK, EL_TEST_REPORT, &
                             EL_TEST_RELEASE_MESH, BUILD_STRUCTURED_MESH, &
                             BUILD_DOF_COORDINATES, nel, nvt, net, nat

  IMPLICIT NONE

  INTEGER :: rank, nproc, ierr, failures

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  CALL EL_TEST_INIT()

  failures = 0
  el_kernel = 'deen_poly'
  el_kernel_width_factor = 2.5d0

  CALL test_exact_polynomial_gather(failures)
  CALL test_smooth_field_convergence(failures)
  CALL test_kernel_width_study(failures)
  CALL test_spread_moments(failures)
  CALL test_gather_scatter_adjoint(failures)

  CALL EL_TEST_REPORT(failures, rank, nproc, 'EL convergence')
  CALL EL_FINALIZE()
  CALL MPI_Finalize(ierr)

  IF (failures.GT.0) ERROR STOP 1

CONTAINS

  SUBROUTINE test_exact_polynomial_gather(nf)

    INTEGER, INTENT(INOUT) :: nf
    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coord(:,:), u(:), v(:), w(:), p(:)
    TYPE(tElParticleRecord) :: particle
    REAL*8 :: sample(EL_SAMPLE_SIZE), carrier(3), grad_u(3,3), grad_p(3)
    REAL*8 :: expected_carrier(3), expected_grad_u(3,3), expected_grad_p(3)
    REAL*8 :: expected_pressure, pressure_value
    INTEGER :: ndof, i
    LOGICAL :: ok

    CALL setup_case(8, mesh, coord, u, v, w, p, ndof)
    particle = make_particle(0.5d0, 0.5d0, 0.5d0, 0.35d0/8.0d0)

    DO i = 1, ndof
      u(i) = 0.3d0 + 0.5d0*coord(1,i) - 0.25d0*coord(2,i)
      v(i) = -0.7d0 + 0.2d0*coord(2,i) + 0.125d0*coord(3,i)
      w(i) = 0.1d0 - 0.4d0*coord(1,i) + 0.3d0*coord(3,i)
    END DO
    CALL fill_linear_pressure(mesh, p, (/0.7d0, 0.2d0, -0.3d0, 0.4d0/))

    CALL EL_INTEGRATE_PARTICLE(mesh, 1, u, v, w, p, el_field_data%epsilon_f, &
                               particle, sample)
    carrier = sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) / &
              sample(EL_SAMPLE_NORMALIZATION)
    grad_u(1,:) = sample(EL_SAMPLE_GRAD_U_BEGIN:EL_SAMPLE_GRAD_U_BEGIN+2) / &
                  sample(EL_SAMPLE_NORMALIZATION)
    grad_u(2,:) = sample(EL_SAMPLE_GRAD_U_BEGIN+3:EL_SAMPLE_GRAD_U_BEGIN+5) / &
                  sample(EL_SAMPLE_NORMALIZATION)
    grad_u(3,:) = sample(EL_SAMPLE_GRAD_U_BEGIN+6:EL_SAMPLE_GRAD_U_END) / &
                  sample(EL_SAMPLE_NORMALIZATION)
    pressure_value = sample(EL_SAMPLE_PRESSURE) / sample(EL_SAMPLE_NORMALIZATION)
    grad_p = sample(EL_SAMPLE_GRAD_P_BEGIN:EL_SAMPLE_GRAD_P_END) / &
             sample(EL_SAMPLE_NORMALIZATION)

    expected_carrier = (/ &
      0.3d0 + 0.5d0*particle%position(1) - 0.25d0*particle%position(2), &
     -0.7d0 + 0.2d0*particle%position(2) + 0.125d0*particle%position(3), &
      0.1d0 - 0.4d0*particle%position(1) + 0.3d0*particle%position(3) /)
    expected_grad_u = 0.0d0
    expected_grad_u(1,:) = (/0.5d0, -0.25d0, 0.0d0/)
    expected_grad_u(2,:) = (/0.0d0, 0.2d0, 0.125d0/)
    expected_grad_u(3,:) = (/-0.4d0, 0.0d0, 0.3d0/)
    expected_pressure = 0.7d0 + 0.2d0*particle%position(1) - &
      0.3d0*particle%position(2) + 0.4d0*particle%position(3)
    expected_grad_p = (/0.2d0, -0.3d0, 0.4d0/)

    ok = sample(EL_SAMPLE_NORMALIZATION).GT.0.0d0
    ok = ok .AND. MAXVAL(ABS(carrier-expected_carrier)).LT.1.0d-10
    ok = ok .AND. MAXVAL(ABS(grad_u-expected_grad_u)).LT.1.0d-10
    ok = ok .AND. ABS(pressure_value-expected_pressure).LT.1.0d-12
    ok = ok .AND. MAXVAL(ABS(grad_p-expected_grad_p)).LT.1.0d-14
    CALL EL_TEST_CHECK(ok, 'exact constant/linear velocity and pressure gather', &
                       nf, rank)

    CALL cleanup_case(mesh, coord, u, v, w, p)

  END SUBROUTINE test_exact_polynomial_gather

  SUBROUTINE test_smooth_field_convergence(nf)

    INTEGER, INTENT(INOUT) :: nf
    INTEGER, PARAMETER :: NC(4) = (/4, 8, 16, 32/)
    REAL*8 :: value_error(4), grad_error(4), h
    REAL*8 :: value_order_1, value_order_2, grad_order_1, grad_order_2
    INTEGER :: i
    LOGICAL :: ok

    DO i = 1, 4
      h = 1.0d0 / DBLE(NC(i))
      CALL smooth_error_at_resolution(NC(i), 0.35d0*h, value_error(i), &
                                      grad_error(i))
    END DO

    value_order_1 = LOG(value_error(1)/value_error(2)) / LOG(2.0d0)
    value_order_2 = LOG(value_error(2)/value_error(3)) / LOG(2.0d0)
    grad_order_1 = LOG(grad_error(1)/grad_error(2)) / LOG(2.0d0)
    grad_order_2 = LOG(grad_error(2)/grad_error(3)) / LOG(2.0d0)

    ok = value_error(4).LT.value_error(3) .AND. value_error(3).LT.value_error(2)
    ok = ok .AND. grad_error(4).LT.grad_error(3) .AND. &
      grad_error(3).LT.grad_error(2)
    ok = ok .AND. MIN(value_order_1,value_order_2).GT.1.2d0
    ok = ok .AND. MIN(grad_order_1,grad_order_2).GT.1.2d0
    CALL EL_TEST_CHECK(ok, 'smooth gather errors decrease with measured EOC', &
                       nf, rank)

  END SUBROUTINE test_smooth_field_convergence

  SUBROUTINE test_kernel_width_study(nf)

    INTEGER, INTENT(INOUT) :: nf
    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coord(:,:), u(:), v(:), w(:), p(:)
    TYPE(tElParticleRecord) :: particle
    REAL*8 :: factors(3), sampled(3), norms(3), center_value, old_factor
    REAL*8 :: sample(EL_SAMPLE_SIZE)
    INTEGER :: ndof, i, g
    LOGICAL :: ok

    old_factor = el_kernel_width_factor
    factors = (/1.25d0, 2.5d0, 3.75d0/)
    CALL setup_case(16, mesh, coord, u, v, w, p, ndof)
    particle = make_particle(0.5d0, 0.5d0, 0.5d0, 0.02d0)
    center_value = SUM(particle%position**2)

    DO g = 1, ndof
      u(g) = SUM(coord(:,g)**2)
    END DO
    v = 0.0d0
    w = 0.0d0
    p = 0.0d0

    DO i = 1, 3
      el_kernel_width_factor = factors(i)
      CALL EL_INTEGRATE_PARTICLE(mesh, 1, u, v, w, p, el_field_data%epsilon_f, &
                                 particle, sample)
      norms(i) = sample(EL_SAMPLE_NORMALIZATION)
      sampled(i) = sample(EL_SAMPLE_UF_BEGIN) / norms(i)
    END DO
    el_kernel_width_factor = old_factor

    ok = MINVAL(norms).GT.0.0d0
    ok = ok .AND. MINVAL(sampled-center_value).GT.0.0d0
    ok = ok .AND. sampled(3).GT.sampled(2) .AND. sampled(2).GT.sampled(1)
    CALL EL_TEST_CHECK(ok, 'kernel width preserves normalization and smoothing trend', &
                       nf, rank)
    CALL cleanup_case(mesh, coord, u, v, w, p)

  END SUBROUTINE test_kernel_width_study

  SUBROUTINE test_spread_moments(nf)

    INTEGER, INTENT(INOUT) :: nf
    INTEGER, PARAMETER :: NC(4) = (/4, 8, 16, 32/)
    REAL*8 :: moment_error(4), zeroth_error(4)
    INTEGER :: i
    LOGICAL :: ok

    DO i = 1, 4
      CALL spread_moment_error(NC(i), moment_error(i), zeroth_error(i))
    END DO

    ok = MAXVAL(zeroth_error).LT.1.0d-10
    CALL EL_TEST_CHECK(ok, 'P2G zeroth moment exact across refinements', &
                       nf, rank)
    ok = MAXVAL(moment_error).LT.1.0d-10
    CALL EL_TEST_CHECK(ok, 'P2G first moment is centered to roundoff', &
                       nf, rank)

  END SUBROUTINE test_spread_moments

  SUBROUTINE test_gather_scatter_adjoint(nf)

    INTEGER, INTENT(INOUT) :: nf
    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coord(:,:), u(:), v(:), w(:), p(:)
    TYPE(tElParticleRecord) :: particle
    REAL*8 :: sample(EL_SAMPLE_SIZE), force(3), carrier(3)
    REAL*8 :: gather_dot, scatter_dot
    INTEGER :: ndof, g
    LOGICAL :: ok

    CALL setup_case(12, mesh, coord, u, v, w, p, ndof)
    particle = make_particle(0.5d0, 0.5d0, 0.5d0, 0.35d0/12.0d0)
    force = (/0.4d0, -0.8d0, 1.1d0/)
    DO g = 1, ndof
      u(g) = 0.2d0 + 0.1d0*coord(1,g) + 0.3d0*coord(2,g)**2
      v(g) = -0.4d0 + 0.2d0*coord(2,g) - 0.1d0*coord(3,g)**2
      w(g) = 0.5d0 - 0.3d0*coord(1,g)*coord(3,g)
    END DO
    p = 0.0d0

    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, u, v, w, p, el_field_data%epsilon_f, &
                               particle, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, particle, sample(EL_SAMPLE_NORMALIZATION), &
                             force, el_field_data)
    carrier = sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) / &
              sample(EL_SAMPLE_NORMALIZATION)
    gather_dot = DOT_PRODUCT(force, carrier)
    scatter_dot = 0.0d0
    DO g = 1, ndof
      scatter_dot = scatter_dot + el_field_data%force_rhs(1,g)*u(g) + &
        el_field_data%force_rhs(2,g)*v(g) + el_field_data%force_rhs(3,g)*w(g)
    END DO
    ok = ABS(gather_dot + scatter_dot).LT.1.0d-11* &
      MAX(1.0d0,ABS(gather_dot))
    CALL EL_TEST_CHECK(ok, 'gather/scatter weighted transpose identity', nf, rank)
    CALL cleanup_case(mesh, coord, u, v, w, p)

  END SUBROUTINE test_gather_scatter_adjoint

  SUBROUTINE smooth_error_at_resolution(ncell, radius, value_error, grad_error)

    INTEGER, INTENT(IN) :: ncell
    REAL*8, INTENT(IN) :: radius
    REAL*8, INTENT(OUT) :: value_error, grad_error
    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coord(:,:), u(:), v(:), w(:), p(:)
    TYPE(tElParticleRecord) :: particle
    REAL*8 :: sample(EL_SAMPLE_SIZE), carrier(3), grad_u(3,3)
    REAL*8 :: exact_value(3), exact_grad(3,3), x(3), pi
    INTEGER :: ndof, g

    CALL setup_case(ncell, mesh, coord, u, v, w, p, ndof)
    particle = make_particle(0.5d0, 0.5d0, 0.5d0, radius)
    pi = ACOS(-1.0d0)
    DO g = 1, ndof
      x = coord(:,g)
      u(g) = SIN(2.0d0*pi*x(1)) + 0.25d0*COS(2.0d0*pi*x(2))
      v(g) = COS(2.0d0*pi*x(2)) + 0.10d0*SIN(2.0d0*pi*x(3))
      w(g) = SIN(2.0d0*pi*x(3)) + 0.15d0*COS(2.0d0*pi*x(1))
    END DO
    p = 0.0d0

    CALL EL_INTEGRATE_PARTICLE(mesh, 1, u, v, w, p, el_field_data%epsilon_f, &
                               particle, sample)
    carrier = sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) / &
              sample(EL_SAMPLE_NORMALIZATION)
    grad_u(1,:) = sample(EL_SAMPLE_GRAD_U_BEGIN:EL_SAMPLE_GRAD_U_BEGIN+2) / &
                  sample(EL_SAMPLE_NORMALIZATION)
    grad_u(2,:) = sample(EL_SAMPLE_GRAD_U_BEGIN+3:EL_SAMPLE_GRAD_U_BEGIN+5) / &
                  sample(EL_SAMPLE_NORMALIZATION)
    grad_u(3,:) = sample(EL_SAMPLE_GRAD_U_BEGIN+6:EL_SAMPLE_GRAD_U_END) / &
                  sample(EL_SAMPLE_NORMALIZATION)

    x = particle%position
    exact_value(1) = SIN(2.0d0*pi*x(1)) + 0.25d0*COS(2.0d0*pi*x(2))
    exact_value(2) = COS(2.0d0*pi*x(2)) + 0.10d0*SIN(2.0d0*pi*x(3))
    exact_value(3) = SIN(2.0d0*pi*x(3)) + 0.15d0*COS(2.0d0*pi*x(1))
    exact_grad = 0.0d0
    exact_grad(1,1) = 2.0d0*pi*COS(2.0d0*pi*x(1))
    exact_grad(1,2) = -0.5d0*pi*SIN(2.0d0*pi*x(2))
    exact_grad(2,2) = -2.0d0*pi*SIN(2.0d0*pi*x(2))
    exact_grad(2,3) = 0.2d0*pi*COS(2.0d0*pi*x(3))
    exact_grad(3,3) = 2.0d0*pi*COS(2.0d0*pi*x(3))
    exact_grad(3,1) = -0.3d0*pi*SIN(2.0d0*pi*x(1))

    value_error = SQRT(SUM((carrier-exact_value)**2))
    grad_error = SQRT(SUM((grad_u-exact_grad)**2))
    CALL cleanup_case(mesh, coord, u, v, w, p)

  END SUBROUTINE smooth_error_at_resolution

  SUBROUTINE spread_moment_error(ncell, moment_error, zeroth_error)

    INTEGER, INTENT(IN) :: ncell
    REAL*8, INTENT(OUT) :: moment_error, zeroth_error
    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coord(:,:), u(:), v(:), w(:), p(:)
    TYPE(tElParticleRecord) :: particle
    REAL*8 :: sample(EL_SAMPLE_SIZE), force(3), zeroth(3), moment(3,3), target(3,3)
    INTEGER :: ndof, i, j, g

    CALL setup_case(ncell, mesh, coord, u, v, w, p, ndof)
    particle = make_particle(0.5d0, 0.5d0, 0.5d0, 0.35d0/DBLE(ncell))
    force = (/0.7d0, -1.3d0, 2.1d0/)
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    p = 0.0d0

    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, u, v, w, p, el_field_data%epsilon_f, &
                               particle, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, particle, sample(EL_SAMPLE_NORMALIZATION), &
                             force, el_field_data)
    zeroth = 0.0d0
    moment = 0.0d0
    DO g = 1, ndof
      zeroth = zeroth + el_field_data%force_rhs(:,g)
      DO i = 1, 3
        DO j = 1, 3
          moment(i,j) = moment(i,j) + el_field_data%force_rhs(i,g)*coord(j,g)
        END DO
      END DO
    END DO
    DO i = 1, 3
      target(i,:) = -force(i)*particle%position
    END DO
    zeroth_error = MAXVAL(ABS(zeroth + force))
    moment_error = SQRT(SUM((moment-target)**2))
    CALL cleanup_case(mesh, coord, u, v, w, p)

  END SUBROUTINE spread_moment_error

  SUBROUTINE setup_case(ncell, mesh, coord, u, v, w, p, ndof)

    INTEGER, INTENT(IN) :: ncell
    TYPE(tMultiMesh), INTENT(INOUT) :: mesh
    REAL*8, ALLOCATABLE, INTENT(OUT) :: coord(:,:), u(:), v(:), w(:), p(:)
    INTEGER, INTENT(OUT) :: ndof

    CALL BUILD_STRUCTURED_MESH(mesh, ncell)
    ndof = nvt + net + nat + nel
    ALLOCATE(coord(3,ndof), u(ndof), v(ndof), w(ndof), p(4*nel))
    CALL BUILD_DOF_COORDINATES(mesh, coord)
    CALL EL_ENSURE_FIELDS(nel, ndof)
    el_field_data%epsilon_f = 1.0d0

  END SUBROUTINE setup_case

  SUBROUTINE cleanup_case(mesh, coord, u, v, w, p)

    TYPE(tMultiMesh), INTENT(INOUT) :: mesh
    REAL*8, ALLOCATABLE, INTENT(INOUT) :: coord(:,:), u(:), v(:), w(:), p(:)

    IF (ALLOCATED(coord)) DEALLOCATE(coord)
    IF (ALLOCATED(u)) DEALLOCATE(u)
    IF (ALLOCATED(v)) DEALLOCATE(v)
    IF (ALLOCATED(w)) DEALLOCATE(w)
    IF (ALLOCATED(p)) DEALLOCATE(p)
    CALL EL_TEST_RELEASE_MESH(mesh)

  END SUBROUTINE cleanup_case

  FUNCTION make_particle(x, y, z, radius) RESULT(particle)

    REAL*8, INTENT(IN) :: x, y, z, radius
    TYPE(tElParticleRecord) :: particle

    particle%position = (/x, y, z/)
    particle%velocity = 0.0d0
    particle%radius = radius
    particle%density = 1.0d0

  END FUNCTION make_particle

  SUBROUTINE fill_linear_pressure(mesh, pressure, coefficients)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    REAL*8, INTENT(OUT) :: pressure(:)
    REAL*8, INTENT(IN) :: coefficients(4)
    REAL*8 :: center(3)
    INTEGER :: iel, i, j

    DO iel = 1, mesh%level(1)%nel
      center = 0.0d0
      DO i = 1, 8
        center = center + 0.125d0* &
          mesh%level(1)%dcorvg(:, mesh%level(1)%kvert(i,iel))
      END DO
      j = 4*(iel-1) + 1
      pressure(j) = coefficients(1) + DOT_PRODUCT(coefficients(2:4), center)
      pressure(j+1:j+3) = coefficients(2:4)
    END DO

  END SUBROUTINE fill_linear_pressure

END PROGRAM TEST_EL_CONVERGENCE
