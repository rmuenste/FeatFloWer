PROGRAM TEST_EL_MOMENTUM_MODEL

  ! Unit tests of the Euler-Lagrange momentum models (EL_MOMENTUM_MODEL, EL_VISCOUS,
  ! EL_MM_OPERATOR, the model parts of EL_CONFIG and EL_FORCES) on the synthetic structured
  ! mesh of el_test_helpers. Not covered here (they live in Transport_Q2P1 and need the
  ! solver): the coefficient update with its restriction and mass rebuild, the two model A
  ! corrections as applied, the refusals (STOP paths; exercised by negative runs), the
  ! multi-rank sum of the viscous record.
  !   1 policy table: drag scale, residual divisor, viscous force / feedback per model
  !   2 key resolution: implied values, registry of explicit keys
  !   3 force split per model: particle and feedback sums against a hand table
  !   4 viscous kernel average: zero for a linear velocity (interior and at a wall, brick and
  !     warped mesh), the Laplacian of a quadratic field in the interior and at a wall
  !   5 start-up balance: the drag scale of each model in the 1-D slip equation gives the
  !     time-constant ratios 1.25 (unweighted inertia) and 0.62 (model B). This pins the policy
  !     and the algebra of the prediction only; that the solver realises the inertia weights is
  !     shown by the runs of validation case s06_startup_inertia, not by a unit test.
  !   6 model A operators: fp is the transpose of the solid-fraction part of the conservative
  !     divergence (element identity against EL_VA_ELEMENT_FLUX_MOMENTS), a constant pressure
  !     gives no force at interior nodes for a continuous solid fraction, fv is the
  !     phi-weighted stiffness (constants annihilated, energy of a linear field)

  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_momentum_model, el_pressure_force, el_pressure_feedback, &
                       el_volume_averaged, el_va_continuity_form, el_drag_coupling, &
                       el_apply_fluid_feedback, el_fluid_gravity, el_va_source_correctors, &
                       el_drag_model, el_kernel, EL_VALIDATE_CONFIG, EL_NOTE_KEY_SET, &
                       el_n_keys_set, el_n_implied, EL_MM_ID, el_fluid_viscosity, &
                       el_closure_nu, EL_CLOSURE_VISCOSITY
  USE EL_MOMENTUM_MODEL, ONLY: EL_MM_ACTIVE, EL_MM_DRAG_SCALE, EL_MM_RESIDUAL_DIVISOR, &
                               EL_MM_VISCOUS_FORCE, EL_MM_VISCOUS_FEEDBACK, &
                               EL_MM_WEIGHTED_INERTIA, EL_MM_OPERATOR_FORM, EL_MM_NONE, EL_MM_A
  USE EL_FORCES, ONLY: tElForceResult, EL_COMPUTE_PARTICLE_FORCES, EL_FLUID_OPERATOR_VISCOSITY
  USE EL_QUADRATURE, ONLY: EL_SAMPLE_SIZE, EL_SAMPLE_NORMALIZATION, EL_SAMPLE_UF_BEGIN, &
                           EL_SAMPLE_GRAD_P_BEGIN, EL_SAMPLE_EPSILON_F
  USE EL_VISCOUS, ONLY: EL_INTEGRATE_PARTICLE_VISCOUS, EL_VISCOUS_LAPLACIAN, EL_VIS_SIZE
  USE EL_MM_OPERATOR, ONLY: EL_MM_ELEMENT_PHI_OPERATORS, EL_MM_PHI_OPERATORS
  USE EL_VA_FLUX, ONLY: tElVaFluxTables, EL_VA_FLUX_TABLES, EL_VA_ELEMENT_DOFS, &
                        EL_VA_ELEMENT_FLUX_MOMENTS, NQ2
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
  CALL test_policy_table(failures)
  CALL test_key_resolution(failures)
  CALL test_viscosity_convention(failures)
  CALL test_force_split(failures)
  CALL test_viscous_average(failures, .FALSE.)
  CALL test_viscous_average(failures, .TRUE.)
  CALL test_startup_balance(failures)
  CALL test_operator_form(failures, .FALSE.)
  CALL test_operator_form(failures, .TRUE.)

  CALL EL_TEST_REPORT(failures, rank, nproc, 'EL momentum model')
  CALL MPI_Finalize(ierr)

  IF (failures.GT.0) ERROR STOP 1

CONTAINS

  ! Select a model the way a deck would: defaults, no explicit keys, then the resolution.
  SUBROUTINE select_model(name)
    CHARACTER(LEN=*), INTENT(IN) :: name
    el_momentum_model = name
    el_pressure_force = .TRUE.
    el_pressure_feedback = .FALSE.
    el_volume_averaged = 'none'
    el_va_continuity_form = 'advective'
    el_drag_coupling = 'explicit'
    el_apply_fluid_feedback = .FALSE.
    el_fluid_gravity = .TRUE.
    el_va_source_correctors = 0
    el_drag_model = 'difelice_std'
    el_n_keys_set = 0
    el_n_implied = 0
    CALL EL_VALIDATE_CONFIG()
  END SUBROUTINE select_model

  ! SimPar@ELFluidViscosity: what the closures and the viscous force receive in the two
  ! conventions (the conversion of the property itself is done by the application).
  SUBROUTINE test_viscosity_convention(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8, PARAMETER :: rho = 10.0d0, nu = 1.5d-4

    el_fluid_viscosity = 'kinematic'
    CALL select_model('none')
    CALL EL_TEST_CHECK(EL_CLOSURE_VISCOSITY(nu).EQ.nu .AND. &
      EL_FLUID_OPERATOR_VISCOSITY(rho, nu).EQ.nu, &
      'viscosity: default = campaign convention (closures nu, operator nu)', nf, rank)
    el_fluid_viscosity = 'Dynamic'
    CALL select_model('none')
    el_closure_nu = nu
    ! the property now holds mu = rho nu; the closures get the recorded nu back bit for bit
    CALL EL_TEST_CHECK(EL_CLOSURE_VISCOSITY(rho*nu).EQ.nu .AND. &
      EL_FLUID_OPERATOR_VISCOSITY(rho, EL_CLOSURE_VISCOSITY(rho*nu)).EQ.rho*nu, &
      'viscosity: dynamic (closures nu, operator and viscous force rho nu)', nf, rank)
    el_fluid_viscosity = 'kinematic'
    el_closure_nu = -1.0d0
    CALL select_model('none')
    CALL EL_TEST_CHECK(EL_FLUID_OPERATOR_VISCOSITY(rho, nu).EQ.nu, &
      'viscosity: back to the default after re-validation', nf, rank)
  END SUBROUTINE test_viscosity_convention

  SUBROUTINE test_policy_table(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8, PARAMETER :: eps = 0.6d0

    CALL select_model('none')
    CALL EL_TEST_CHECK(.NOT.EL_MM_ACTIVE() .AND. .NOT.EL_MM_WEIGHTED_INERTIA() .AND. &
      EL_MM_DRAG_SCALE(eps).EQ.1.0d0 .AND. EL_MM_RESIDUAL_DIVISOR(eps).EQ.eps .AND. &
      .NOT.EL_MM_VISCOUS_FORCE(), 'policy: none = campaign set (divisor eps, no scale)', nf, rank)
    CALL select_model('Bfull')
    CALL EL_TEST_CHECK(EL_MM_WEIGHTED_INERTIA() .AND. EL_MM_DRAG_SCALE(eps).EQ.1.0d0 .AND. &
      EL_MM_RESIDUAL_DIVISOR(eps).EQ.1.0d0 .AND. EL_MM_VISCOUS_FORCE() .AND. &
      EL_MM_VISCOUS_FEEDBACK(), 'policy: Bfull (divisor 1, viscous force fed back)', nf, rank)
    CALL select_model('B')
    CALL EL_TEST_CHECK(EL_MM_WEIGHTED_INERTIA() .AND. &
      ABS(EL_MM_DRAG_SCALE(eps)-1.0d0/eps).LT.1.0d-15 .AND. &
      EL_MM_RESIDUAL_DIVISOR(eps).EQ.1.0d0 .AND. .NOT.EL_MM_VISCOUS_FORCE(), &
      'policy: B (drag / eps, divisor 1, no viscous force)', nf, rank)
    CALL select_model('A')
    CALL EL_TEST_CHECK(EL_MM_WEIGHTED_INERTIA() .AND. EL_MM_OPERATOR_FORM() .AND. &
      EL_MM_DRAG_SCALE(eps).EQ.1.0d0 .AND. EL_MM_RESIDUAL_DIVISOR(eps).EQ.eps .AND. &
      EL_MM_VISCOUS_FORCE() .AND. .NOT.EL_MM_VISCOUS_FEEDBACK(), &
      'policy: A (operator form, divisor eps, viscous force without feedback)', nf, rank)
    CALL select_model('bfull')
    CALL EL_TEST_CHECK(.NOT.EL_MM_OPERATOR_FORM(), 'policy: Bfull is not the operator form', &
      nf, rank)
    CALL select_model('none')
    CALL EL_TEST_CHECK(.NOT.EL_MM_OPERATOR_FORM() .AND. EL_MM_ID().EQ.EL_MM_NONE, &
      'policy: none has no operator form, id reset by the validation', nf, rank)
  END SUBROUTINE test_policy_table

  SUBROUTINE test_key_resolution(nf)
    INTEGER, INTENT(INOUT) :: nf

    CALL select_model('b')
    CALL EL_TEST_CHECK(TRIM(el_volume_averaged).EQ.'continuity' .AND. &
      TRIM(el_va_continuity_form).EQ.'conservative' .AND. &
      TRIM(el_drag_coupling).EQ.'semi_implicit' .AND. el_apply_fluid_feedback .AND. &
      .NOT.el_pressure_force .AND. .NOT.el_pressure_feedback .AND. .NOT.el_fluid_gravity &
      .AND. el_va_source_correctors.EQ.1, 'keys: model B implies its set', nf, rank)
    CALL EL_TEST_CHECK(el_n_implied.EQ.8, 'keys: eight implied lines for the echo', nf, rank)
    CALL select_model('bfull')
    CALL EL_TEST_CHECK(el_pressure_force .AND. el_pressure_feedback, &
      'keys: Bfull implies pressure force and feedback', nf, rank)
    CALL select_model('a')
    CALL EL_TEST_CHECK(el_pressure_force .AND. .NOT.el_pressure_feedback .AND. &
      EL_MM_ID().EQ.EL_MM_A .AND. TRIM(el_va_continuity_form).EQ.'conservative', &
      'keys: A implies the pressure force without feedback', nf, rank)
    ! an explicit key with the implied value is accepted and a larger corrector count kept
    el_momentum_model = 'bfull'
    el_n_keys_set = 0
    CALL EL_NOTE_KEY_SET('ELPressureForce')
    CALL EL_NOTE_KEY_SET('ELVASourceCorrectors')
    el_va_source_correctors = 3
    CALL EL_VALIDATE_CONFIG()
    CALL EL_TEST_CHECK(el_va_source_correctors.EQ.3 .AND. EL_MM_ID().EQ.1, &
      'keys: explicit consistent values are kept', nf, rank)
    CALL select_model('none')
    CALL EL_TEST_CHECK(el_pressure_force .AND. .NOT.el_pressure_feedback .AND. &
      TRIM(el_volume_averaged).EQ.'none', 'keys: none leaves the campaign keys alone', nf, rank)
  END SUBROUTINE test_key_resolution

  SUBROUTINE test_force_split(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: state(8), sample(EL_SAMPLE_SIZE), gravity(3), lap(3)
    REAL*8, PARAMETER :: rho = 1000.0d0, nu = 1.0d-6, eps = 0.6d0
    TYPE(tElForceResult) :: r0, rb, rf, ra
    REAL*8 :: vol, fv(3)

    state = 0.0d0
    state(7) = 5.0d-4
    state(8) = 1100.0d0
    sample = 0.0d0
    sample(EL_SAMPLE_NORMALIZATION) = 1.0d0
    sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_BEGIN+2) = (/ 0.0d0, 0.0d0, 5.0d-3 /)
    sample(EL_SAMPLE_GRAD_P_BEGIN:EL_SAMPLE_GRAD_P_BEGIN+2) = (/ 0.0d0, 0.0d0, -40.0d0 /)
    sample(EL_SAMPLE_EPSILON_F) = eps
    gravity = (/ 0.0d0, 0.0d0, -9.81d0 /)
    lap = (/ 0.0d0, 0.0d0, 3.0d0 /)
    vol = 4.0d0*ACOS(-1.0d0)*state(7)**3/3.0d0

    CALL select_model('none')
    el_pressure_feedback = .TRUE.
    CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, rho, nu, gravity, r0, lap)
    CALL EL_TEST_CHECK(ALL(r0%viscous.EQ.0.0d0) .AND. &
      ALL(r0%particle_total.EQ.r0%drag+r0%pressure+r0%lift+r0%grav_buoy) .AND. &
      ALL(r0%particle_nondrag.EQ.r0%pressure+r0%lift+r0%grav_buoy) .AND. &
      ALL(r0%feedback_force.EQ.r0%drag+r0%lift+r0%pressure), &
      'forces: none ignores the viscous sample, sums as before', nf, rank)

    CALL select_model('b')
    CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, rho, nu, gravity, rb, lap)
    CALL EL_TEST_CHECK(MAXVAL(ABS(rb%drag-r0%drag/eps)).LE.1.0d-14*MAXVAL(ABS(r0%drag)) .AND. &
      ABS(rb%drag_B-r0%drag_B/eps).LE.1.0d-14*r0%drag_B .AND. ALL(rb%pressure.EQ.0.0d0) .AND. &
      ALL(rb%viscous.EQ.0.0d0) .AND. ALL(rb%feedback_force.EQ.rb%drag+rb%lift), &
      'forces: B = drag / eps, no pressure or viscous force', nf, rank)

    CALL select_model('bfull')
    CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, rho, nu, gravity, rf, lap)
    fv = vol*EL_FLUID_OPERATOR_VISCOSITY(rho, nu)*lap
    CALL EL_TEST_CHECK(MAXVAL(ABS(rf%viscous-fv)).LE.1.0d-30 .AND. ALL(rf%drag.EQ.r0%drag) .AND. &
      ALL(rf%pressure.EQ.r0%pressure) .AND. &
      MAXVAL(ABS(rf%particle_total-(r0%particle_total+fv))).LE.1.0d-20 .AND. &
      MAXVAL(ABS(rf%feedback_force-(r0%feedback_force+fv))).LE.1.0d-20 .AND. &
      MAXVAL(ABS(rf%feedback_nondrag-(rf%pressure+rf%lift+fv))).LE.1.0d-20 .AND. &
      MAXVAL(ABS(rf%particle_nondrag-(rf%pressure+rf%lift+rf%grav_buoy+fv))).LE.1.0d-20, &
      'forces: Bfull adds V c lap u to particle and feedback', nf, rank)
    CALL select_model('a')
    CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, rho, nu, gravity, ra, lap)
    CALL EL_TEST_CHECK(MAXVAL(ABS(ra%viscous-fv)).LE.1.0d-30 .AND. &
      ALL(ra%pressure.EQ.r0%pressure) .AND. ALL(ra%drag.EQ.r0%drag) .AND. &
      MAXVAL(ABS(ra%particle_total-rf%particle_total)).LE.1.0d-20 .AND. &
      ALL(ra%feedback_force.EQ.ra%drag+ra%lift) .AND. ALL(ra%feedback_nondrag.EQ.ra%lift), &
      'forces: A gives the particle pressure and viscous force, feeds back the drag only', &
      nf, rank)
    CALL select_model('none')
  END SUBROUTINE test_force_split

  SUBROUTINE test_viscous_average(nf, warped)
    INTEGER, INTENT(INOUT) :: nf
    LOGICAL, INTENT(IN) :: warped
    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:), u(:), v(:), w(:)
    LOGICAL, ALLOCATABLE :: wall(:)
    INTEGER, PARAMETER :: n = 6
    INTEGER :: ndof, i
    REAL*8 :: vis(EL_VIS_SIZE), lap(3), x(3), h, pos(3), width, err
    CHARACTER(LEN=8) :: tag

    tag = 'brick'
    IF (warped) tag = 'warped'
    el_kernel = 'deen_poly'
    CALL BUILD_STRUCTURED_MESH(mesh, n)
    ndof = nvt + net + nat + nel
    h = 1.0d0/DBLE(n)
    IF (warped) THEN
      DO i = 1, nvt
        x = mesh%level(1)%dcorvg(:,i)
        IF (MINVAL(x).GT.1.0d-12 .AND. MAXVAL(x).LT.1.0d0-1.0d-12) THEN
          mesh%level(1)%dcorvg(1,i) = x(1) + 0.12d0*h*SIN(3.1d0*x(2)+1.3d0*x(3))
          mesh%level(1)%dcorvg(2,i) = x(2) + 0.10d0*h*COS(2.2d0*x(1)+0.7d0*x(3))
          mesh%level(1)%dcorvg(3,i) = x(3) + 0.11d0*h*SIN(1.7d0*x(1)-2.1d0*x(2))
        END IF
      END DO
    END IF
    ALLOCATE(coor(3,ndof), u(ndof), v(ndof), w(ndof), wall(ndof))
    CALL BUILD_DOF_COORDINATES(mesh, coor)
    DO i = 1, ndof
      wall(i) = (MINVAL(coor(:,i)).LT.1.0d-12 .OR. MAXVAL(coor(:,i)).GT.1.0d0-1.0d-12)
    END DO
    width = 2.5d0*h

    ! linear field: zero to rounding, interior and with the kernel cut by a wall
    DO i = 1, ndof
      u(i) = 1.0d0 + 2.0d0*coor(1,i) - 3.0d0*coor(2,i) + 0.5d0*coor(3,i)
      v(i) = -0.7d0*coor(1,i) + 0.2d0*coor(3,i)
      w(i) = 4.0d0*coor(2,i) - coor(3,i)
    END DO
    pos = (/ 0.52d0, 0.47d0, 0.51d0 /)
    CALL EL_INTEGRATE_PARTICLE_VISCOUS(mesh, 1, coor, u, v, w, wall, pos, width, vis)
    lap = EL_VISCOUS_LAPLACIAN(vis)
    CALL EL_TEST_CHECK(MAXVAL(ABS(lap)).LE.1.0d-9, 'viscous '//TRIM(tag)// &
      ': linear field gives zero in the interior', nf, rank)
    pos = (/ 0.52d0, 0.47d0, 0.6d0*h /)
    CALL EL_INTEGRATE_PARTICLE_VISCOUS(mesh, 1, coor, u, v, w, wall, pos, width, vis)
    lap = EL_VISCOUS_LAPLACIAN(vis)
    CALL EL_TEST_CHECK(MAXVAL(ABS(lap)).LE.1.0d-9, 'viscous '//TRIM(tag)// &
      ': linear field gives zero at a wall', nf, rank)

    ! quadratic field (exact in Q2 on the brick): lap u = (2, -4, 6)
    DO i = 1, ndof
      u(i) = coor(3,i)**2 + coor(1,i)*coor(2,i)
      v(i) = -2.0d0*coor(1,i)**2 + coor(3,i)
      w(i) = coor(1,i)**2 + coor(2,i)**2 + coor(3,i)**2 + coor(2,i)*coor(3,i)
    END DO
    pos = (/ 0.52d0, 0.47d0, 0.51d0 /)
    CALL EL_INTEGRATE_PARTICLE_VISCOUS(mesh, 1, coor, u, v, w, wall, pos, width, vis)
    lap = EL_VISCOUS_LAPLACIAN(vis)
    err = MAXVAL(ABS(lap - (/ 2.0d0, -4.0d0, 6.0d0 /)))/6.0d0
    IF (rank.EQ.0) WRITE(*,'(A,A,A,3ES13.5,A,ES10.2)') '    viscous ', TRIM(tag), &
      ' interior lap = ', lap, '  rel. error ', err
    CALL EL_TEST_CHECK(err.LE.0.003d0, 'viscous '//TRIM(tag)// &
      ': Laplacian of a quadratic field, interior (0.3 %)', nf, rank)
    pos = (/ 0.52d0, 0.47d0, 0.6d0*h /)
    CALL EL_INTEGRATE_PARTICLE_VISCOUS(mesh, 1, coor, u, v, w, wall, pos, width, vis)
    lap = EL_VISCOUS_LAPLACIAN(vis)
    err = MAXVAL(ABS(lap - (/ 2.0d0, -4.0d0, 6.0d0 /)))/6.0d0
    IF (rank.EQ.0) WRITE(*,'(A,A,A,3ES13.5,A,ES10.2)') '    viscous ', TRIM(tag), &
      ' wall     lap = ', lap, '  rel. error ', err
    CALL EL_TEST_CHECK(err.LE.0.003d0, 'viscous '//TRIM(tag)// &
      ': Laplacian of a quadratic field, kernel cut by a wall (0.3 %)', nf, rank)
    ! without the boundary term the wall value is badly wrong: the term matters
    wall = .FALSE.
    CALL EL_INTEGRATE_PARTICLE_VISCOUS(mesh, 1, coor, u, v, w, wall, pos, width, vis)
    lap = EL_VISCOUS_LAPLACIAN(vis)
    IF (rank.EQ.0) WRITE(*,'(A,A,A,3ES13.5,A,3ES12.4)') '    viscous ', TRIM(tag), &
      ' no-wall  lap = ', lap, '  B = ', vis(4:6)/vis(16)
    CALL EL_TEST_CHECK(MAXVAL(ABS(lap - (/ 2.0d0, -4.0d0, 6.0d0 /)))/6.0d0.GT.0.1d0, &
      'viscous '//TRIM(tag)//': dropping the wall term is detected', nf, rank)

    DEALLOCATE(coor, u, v, w, wall)
    CALL EL_TEST_RELEASE_MESH(mesh)
  END SUBROUTINE test_viscous_average

  ! 1-D start-up of a homogeneous suspension, zero net flux: with the fluid inertia weight
  ! m_f (eps for weighted, 1 for unweighted inertia), the drag scale and the pressure split of
  ! each model, the slip obeys  M ds/dt = R - K s  with a linear drag K. The time constants
  ! follow from the same policy functions the solver uses.
  SUBROUTINE test_startup_balance(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8, PARAMETER :: phi = 0.4d0, eps = 0.6d0, rf = 1000.0d0, rp = 1100.0d0, beta = 2.0d4
    REAL*8 :: tau_a, tau_unw, tau_b, scale

    ! fluid: m_f rho_f du/dt = -w_p grad p - S beta s ; particle: phi rho_p dv/dt =
    ! -phi_p grad p + S beta s - phi drho g, with u = phi s, v = -eps s.
    ! pressure force on the particle and its feedback (none+pf+fb, Bfull): w_p = eps, phi_p = phi
    CALL select_model('bfull')
    scale = EL_MM_DRAG_SCALE(eps)
    tau_a = startup_tau(eps, .TRUE., scale, phi, eps, rf, rp, beta)
    CALL select_model('none')
    scale = EL_MM_DRAG_SCALE(eps)
    tau_unw = startup_tau(1.0d0, .TRUE., scale, phi, eps, rf, rp, beta)
    CALL select_model('b')
    scale = EL_MM_DRAG_SCALE(eps)
    tau_b = startup_tau(eps, .FALSE., scale, phi, eps, rf, rp, beta)
    CALL select_model('none')
    CALL EL_TEST_CHECK(ABS(tau_a - eps*phi*(rf*phi+rp*eps)/beta).LE.1.0d-12*tau_a, &
      'start-up: weighted inertia gives tau_A', nf, rank)
    CALL EL_TEST_CHECK(ABS(tau_unw/tau_a - 1.2516d0).LE.5.0d-4, &
      'start-up: unweighted inertia is 1.25 x', nf, rank)
    CALL EL_TEST_CHECK(ABS(tau_b/tau_a - 0.6226d0).LE.5.0d-4, &
      'start-up: model B is 0.62 x', nf, rank)

  END SUBROUTINE test_startup_balance

  SUBROUTINE test_operator_form(nf, warped)
    INTEGER, INTENT(INOUT) :: nf
    LOGICAL, INTENT(IN) :: warped
    TYPE(tMultiMesh) :: mesh
    TYPE(tElVaFluxTables) :: tab
    REAL*8, ALLOCATABLE :: coor(:,:), eps(:), u(:), v(:), w(:), pr(:), fp(:,:), fv(:,:)
    INTEGER, PARAMETER :: n = 4
    INTEGER :: ndof, i, iel, l, dofs(NQ2)
    REAL*8 :: x(3), h, xe(3,NQ2), ee(NQ2), ue(3,NQ2), one(NQ2), pe(4), efp(3,NQ2), efv(3,NQ2)
    REAL*8 :: f1(4), fe(4), eu(3), vu(3), geom(0:3), divu, lhs, rhs, worst, scale, tol, s
    CHARACTER(LEN=8) :: tag

    tag = 'brick'
    IF (warped) tag = 'warped'
    CALL EL_VA_FLUX_TABLES(tab)
    CALL BUILD_STRUCTURED_MESH(mesh, n)
    ndof = nvt + net + nat + nel
    h = 1.0d0/DBLE(n)
    IF (warped) THEN
      DO i = 1, nvt
        x = mesh%level(1)%dcorvg(:,i)
        IF (MINVAL(x).GT.1.0d-12 .AND. MAXVAL(x).LT.1.0d0-1.0d-12) THEN
          mesh%level(1)%dcorvg(1,i) = x(1) + 0.12d0*h*SIN(3.1d0*x(2)+1.3d0*x(3))
          mesh%level(1)%dcorvg(2,i) = x(2) + 0.10d0*h*COS(2.2d0*x(1)+0.7d0*x(3))
          mesh%level(1)%dcorvg(3,i) = x(3) + 0.11d0*h*SIN(1.7d0*x(1)-2.1d0*x(2))
        END IF
      END DO
    END IF
    ALLOCATE(coor(3,ndof), eps(ndof), u(ndof), v(ndof), w(ndof), pr(4*nel), fp(3,ndof), fv(3,ndof))
    CALL BUILD_DOF_COORDINATES(mesh, coor)
    ! trilinear porosity (the transport porosity is Q1 embedded in Q2), smooth Q2 velocity
    DO i = 1, ndof
      x = coor(:,i)
      eps(i) = 0.55d0 + 0.2d0*x(1) - 0.15d0*x(2) + 0.1d0*x(3) + 0.12d0*x(1)*x(2)*x(3)
      u(i) = SIN(1.3d0*x(1)+0.4d0) + x(2)*x(3)
      v(i) = x(1)**2 - 0.5d0*x(2) + COS(0.9d0*x(3))
      w(i) = 0.3d0 + x(1)*x(2) - x(3)**2
    END DO
    IF (warped) THEN
      ! on the warped mesh the nodal interpolant of the trilinear function is not trilinear in
      ! the reference element; the identity then holds to quadrature accuracy only (7e-9 here)
      tol = 1.0d-6
    ELSE
      tol = 1.0d-12
    END IF

    ! element identity  sum_i u_i . fp_i(p) = sum_m p_m [F_m(u;1) - F_m(u;eps)]
    one = 1.0d0
    worst = 0.0d0
    DO iel = 1, nel
      CALL EL_VA_ELEMENT_DOFS(mesh, 1, iel, dofs)
      DO l = 1, NQ2
        xe(:,l) = coor(:,dofs(l))
        ee(l) = eps(dofs(l))
        ue(:,l) = (/ u(dofs(l)), v(dofs(l)), w(dofs(l)) /)
      END DO
      pe = (/ 1.7d0, -0.6d0, 0.9d0, 2.3d0 /) + 0.01d0*DBLE(iel)
      CALL EL_MM_ELEMENT_PHI_OPERATORS(tab, xe, 1.0d0-ee, pe, ue, .FALSE., efp, efv)
      CALL EL_VA_ELEMENT_FLUX_MOMENTS(tab, xe, one, ue, f1, eu, vu, geom, divu)
      CALL EL_VA_ELEMENT_FLUX_MOMENTS(tab, xe, ee, ue, fe, eu, vu, geom, divu)
      lhs = SUM(ue*efp)
      rhs = DOT_PRODUCT(pe, f1 - fe)
      scale = MAX(ABS(lhs), ABS(rhs), 1.0d-30)
      worst = MAX(worst, ABS(lhs-rhs)/scale)
    END DO
    IF (rank.EQ.0) WRITE(*,'(A,A,A,ES10.2)') '    operator ', TRIM(tag), &
      ' transpose identity, worst relative difference ', worst
    CALL EL_TEST_CHECK(worst.LE.tol, 'operator '//TRIM(tag)// &
      ': fp is the transpose of the solid-fraction divergence', nf, rank)

    ! constant pressure: no force at interior nodes
    pr = 0.0d0
    DO iel = 1, nel
      pr(4*(iel-1)+1) = 7.0d0
    END DO
    CALL EL_MM_PHI_OPERATORS(mesh, 1, coor, eps, pr, u, v, w, .TRUE., fp, fv)
    worst = 0.0d0
    scale = 0.0d0
    DO i = 1, ndof
      scale = MAX(scale, MAXVAL(ABS(fp(:,i))))
      IF (MINVAL(coor(:,i)).GT.1.0d-12 .AND. MAXVAL(coor(:,i)).LT.1.0d0-1.0d-12) &
        worst = MAX(worst, MAXVAL(ABS(fp(:,i))))
    END DO
    IF (rank.EQ.0) WRITE(*,'(A,A,A,2ES10.2)') '    operator ', TRIM(tag), &
      ' constant pressure: interior max / boundary scale ', worst, scale
    CALL EL_TEST_CHECK(worst.LE.MAX(tol,1.0d-12)*scale .AND. scale.GT.0.0d0, 'operator '// &
      TRIM(tag)//': constant pressure gives no interior force', nf, rank)

    ! fv: constants annihilated; energy of u = (x,0,0) with constant phi is phi * volume
    u = 3.0d0
    v = -1.0d0
    w = 0.5d0
    CALL EL_MM_PHI_OPERATORS(mesh, 1, coor, eps, pr, u, v, w, .TRUE., fp, fv)
    CALL EL_TEST_CHECK(MAXVAL(ABS(fv)).LE.1.0d-12, 'operator '//TRIM(tag)// &
      ': fv annihilates a constant velocity', nf, rank)
    eps = 0.6d0
    u = coor(1,:)
    v = 0.0d0
    w = 0.0d0
    CALL EL_MM_PHI_OPERATORS(mesh, 1, coor, eps, pr, u, v, w, .TRUE., fp, fv)
    s = SUM(u*fv(1,:))
    CALL EL_TEST_CHECK(ABS(s - 0.4d0).LE.1.0d-10 .AND. MAXVAL(ABS(fv(2:3,:))).LE.1.0d-12, &
      'operator '//TRIM(tag)//': fv is the phi-weighted stiffness', nf, rank)

    DEALLOCATE(coor, eps, u, v, w, pr, fp, fv)
    CALL EL_TEST_RELEASE_MESH(mesh)
  END SUBROUTINE test_operator_form

  ! time constant mass / K of the start-up balance (see test_startup_balance)
  REAL*8 FUNCTION startup_tau(mf, pressure_split, drag_scale, phi, eps, rf, rp, beta)
    REAL*8, INTENT(IN) :: mf, drag_scale, phi, eps, rf, rp, beta
    LOGICAL, INTENT(IN) :: pressure_split
    REAL*8 :: mass, k
    IF (pressure_split) THEN
      ! eliminate grad p between fluid (weight eps) and particle (weight phi)
      mass = mf*rf*phi*phi/eps + rp*eps*phi
      k = drag_scale*beta*(phi/eps + 1.0d0)
    ELSE
      ! particle equation alone
      mass = rp*eps*phi
      k = drag_scale*beta
    END IF
    startup_tau = mass/k
  END FUNCTION startup_tau

END PROGRAM TEST_EL_MOMENTUM_MODEL
