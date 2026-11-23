PROGRAM TEST_EL_TRANSFER
!========================================================================
! Scaffolded MPI test for the Euler-Lagrange particle-mesh transfer.
!
! Phase 1 of the E-L plan makes the conservation/interpolation test suite
! the acceptance gate for the transfer module. This program covers the
! parts that can be exercised standalone against a synthetic structured
! Q2/P1 hex mesh, on 1/2/8 ranks:
!
!   * positive compact kernel support and cut-off,
!   * exact constant-field kernel interpolation,
!   * second-order-clean linear-field interpolation away from walls,
!   * central-particle volume conservation  (sum_K |K| alpha_p = V_i),
!   * wall/truncation-adjacent volume conservation (renormalization),
!   * componentwise conservative Q2 force spreading (sum_a b_a = -F_i),
!   * SimPar@ELDragImplicitWeights = hrz: nodal HRZ coefficients against an independent
!     implementation (own Q2 basis, own trilinear map, 27-point Gauss) node by node on the
!     regular and on a distorted mesh; the explicit partner B_i u_i^n node by node; a steady
!     non-uniform field gives the fluid exactly minus the particle drag impulse (review
!     finding 2; the former partner -B u_f with the HRZ weights failed this, job 164663),
!   * zero-particle / zero-force behaviour,
!   * rank-independence of the transfer math (1/2/8 ranks),
!   * halo bounding-box / sphere overlap predicates.
!
! NOTE (scaffold boundary): the genuinely *distributed* transfer path -
! EL_REFRESH_COUPLING_HALO -> EL_INTEGRATE_PARTICLE (partial integrals) ->
! EL_REDUCE_TO_OWNERS -> EL_DEPOSIT_PARTICLE -> E013Sum3 - with a particle
! whose kernel straddles a partition boundary cannot be reproduced here:
! it needs MPI_COMM_SUBS, the partitioned production mesh, the shared-DOF
! reduction tables, and PE ownership, all of which only exist after the
! full application init. Those cases are marked TODO below and must be
! driven through a small q2p1_el_pipeflow regression harness on the staged
! unit_cube_27_case (face/edge/corner straddling particles on 2 and 8
! ranks). This file deliberately keeps the assertable transfer-math
! coverage so CI is green on 1/2/8 ranks today.
!========================================================================

  USE MPI
  USE TYPES,             ONLY: tMultiMesh
  USE EL_CONFIG,         ONLY: el_kernel, el_kernel_width_factor
  USE EL_KERNEL_FUNCTIONS, ONLY: EL_KERNEL_VALUE
  USE EL_FIELDS,         ONLY: el_field_data, EL_ENSURE_FIELDS, &
                               EL_BEGIN_FIELD_UPDATE, &
                               EL_CAPTURE_FLUID_FEEDBACK_SOURCE, &
                               EL_APPLY_FLUID_FEEDBACK_SOURCE, &
                               EL_APPLY_DRAG_SPLIT_PARTNER, EL_FINALIZE
  USE EL_HALO,           ONLY: tElParticleRecord, EL_BOXES_WITHIN_DELTA, &
                               EL_SPHERE_INTERSECTS_BOX
  USE EL_QUADRATURE,     ONLY: EL_SAMPLE_SIZE, EL_INTEGRATE_PARTICLE, &
                               EL_DEPOSIT_PARTICLE
  USE EL_TEST_HELPERS,   ONLY: EL_TEST_INIT, EL_TEST_RELEASE_MESH, &
                               BUILD_STRUCTURED_MESH, BUILD_DOF_COORDINATES, &
                               nel, nvt, net, nat

  IMPLICIT NONE

  TYPE(tMultiMesh) :: mesh
  REAL*8, ALLOCATABLE :: val_u(:), val_v(:), val_w(:), val_p(:)
  REAL*8, ALLOCATABLE :: coord(:,:)
  INTEGER :: ndof, ncell, rank, nproc, ierr, failures
  REAL*8  :: pi, h

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  CALL EL_TEST_INIT()
  pi     = ACOS(-1.0d0)
  failures = 0

  ! 4x4x4 hexes on the unit cube; h = 0.25.
  ncell = 4
  h = 1.0d0 / DBLE(ncell)
  el_kernel = 'deen_poly'
  el_kernel_width_factor = 2.5d0

  CALL build_structured_mesh(mesh, ncell)
  ndof = nvt + net + nat + nel
  ALLOCATE(val_u(ndof), val_v(ndof), val_w(ndof), val_p(4*nel))
  ALLOCATE(coord(3, ndof))
  CALL build_dof_coordinates(mesh, coord)
  CALL EL_ENSURE_FIELDS(nel, ndof)

  CALL test_kernel_positive_compact(failures)
  CALL test_constant_interpolation(failures)
  CALL test_linear_interpolation(failures)
  CALL test_central_volume_conservation(failures)
  CALL test_wall_volume_conservation(failures)
  CALL test_force_spread_conservation(failures)
  CALL test_drag_B_spread_conservation(failures)
  CALL test_drag_B_hrz_weights(failures)
  CALL test_hrz_node_by_node(failures)
  CALL test_feedback_source_capture(failures)
  CALL test_feedback_source_apply(failures)
  CALL test_zero_force(failures)
  CALL test_halo_predicates(failures)
  CALL test_rank_independence(failures)
  CALL test_partitioned_feedback_no_double_count(failures)

  CALL report(failures, rank, nproc)

  CALL EL_TEST_RELEASE_MESH(mesh)
  CALL EL_FINALIZE()
  DEALLOCATE(val_u, val_v, val_w, val_p, coord)
  CALL MPI_Finalize(ierr)

  IF (failures.GT.0) ERROR STOP 1

CONTAINS

  !----------------------------------------------------------------------
  ! Test cases
  !----------------------------------------------------------------------
  SUBROUTINE test_kernel_positive_compact(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: width, d
    INTEGER :: i
    LOGICAL :: ok
    width = 0.2d0
    ok = .TRUE.
    DO i = 0, 100
      d = width*DBLE(i)/100.0d0
      IF (EL_KERNEL_VALUE(d, width).LT.0.0d0) ok = .FALSE.
    END DO
    ! Compact support: zero at and beyond the cut-off radius.
    IF (EL_KERNEL_VALUE(width, width).NE.0.0d0) ok = .FALSE.
    IF (EL_KERNEL_VALUE(1.5d0*width, width).NE.0.0d0) ok = .FALSE.
    CALL check(ok, 'kernel positive + compact support', nf)
  END SUBROUTINE test_kernel_positive_compact

  SUBROUTINE test_constant_interpolation(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: sample(EL_SAMPLE_SIZE), carrier(3)
    p = central_particle()
    val_u = 1.0d0; val_v = 2.0d0; val_w = 3.0d0; val_p = 0.0d0
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL check(sample(1).GT.0.0d0, 'constant: particle has fluid support', nf)
    carrier = sample(2:4) / sample(1)
    CALL check(MAXVAL(ABS(carrier - (/1.0d0, 2.0d0, 3.0d0/))).LT.1.0d-10, &
               'constant-field interpolation is exact', nf)
  END SUBROUTINE test_constant_interpolation

  SUBROUTINE test_linear_interpolation(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: sample(EL_SAMPLE_SIZE), carrier(3), expected(3)
    INTEGER :: g
    p = central_particle()
    ! u = 0.3 + 0.5 x,  v = -0.1 + 0.25 y,  w = 0.2 + 0.75 z
    DO g = 1, ndof
      val_u(g) = 0.3d0  + 0.5d0 *coord(1, g)
      val_v(g) = -0.1d0 + 0.25d0*coord(2, g)
      val_w(g) = 0.2d0  + 0.75d0*coord(3, g)
    END DO
    val_p = 0.0d0
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    carrier = sample(2:4) / sample(1)
    expected(1) = 0.3d0  + 0.5d0 *p%position(1)
    expected(2) = -0.1d0 + 0.25d0*p%position(2)
    expected(3) = 0.2d0  + 0.75d0*p%position(3)
    ! Symmetric interior support => kernel-weighted mean of a linear field
    ! reproduces it at X_i to within the cubature/discretisation error.
    ! Tolerance is intentionally loose; tighten once two-width O(delta^2)
    ! convergence is added.
    CALL check(MAXVAL(ABS(carrier - expected)).LT.1.0d-6, &
               'linear-field interpolation away from walls', nf)
  END SUBROUTINE test_linear_interpolation

  SUBROUTINE test_central_volume_conservation(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: deposited, expected
    deposited = deposit_volume(central_particle(), expected)
    CALL check(ABS(deposited - expected).LE.1.0d-10*expected, &
               'central particle: sum |K| alpha_p = V_i', nf)
  END SUBROUTINE test_central_volume_conservation

  SUBROUTINE test_wall_volume_conservation(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: deposited, expected
    ! Particle near the x=0 wall: support is truncated by the domain
    ! boundary and renormalised over the fluid elements only.
    deposited = deposit_volume(wall_particle(), expected)
    CALL check(ABS(deposited - expected).LE.1.0d-10*expected, &
               'wall-adjacent particle: conservation after renorm', nf)
  END SUBROUTINE test_wall_volume_conservation

  SUBROUTINE test_force_spread_conservation(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: sample(EL_SAMPLE_SIZE), force(3), sumf(3)
    INTEGER :: c
    LOGICAL :: ok
    p = central_particle()
    force = (/0.7d0, -1.3d0, 2.1d0/)
    val_u = 0.0d0; val_v = 0.0d0; val_w = 0.0d0; val_p = 0.0d0
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), force, el_field_data)
    ok = .TRUE.
    DO c = 1, 3
      sumf(c) = SUM(el_field_data%force_rhs(c, :))
      IF (ABS(sumf(c) + force(c)).GT.1.0d-10*(ABS(force(c)) + 1.0d0)) ok = .FALSE.
    END DO
    CALL check(ok, 'force spreading: sum_a b_a = -F_i (componentwise)', nf)
  END SUBROUTINE test_force_spread_conservation

  SUBROUTINE test_drag_B_spread_conservation(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: sample(EL_SAMPLE_SIZE), force(3), drag_b, deposited
    p = central_particle()
    force = 0.0d0
    drag_b = 3.75d0
    val_u = 0.0d0; val_v = 0.0d0; val_w = 0.0d0; val_p = 0.0d0
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), force, el_field_data, &
                             drag_b=drag_b)
    deposited = SUM(el_field_data%drag_B_ml)
    CALL check(ABS(deposited-drag_b).LE.1.0d-10*(ABS(drag_b)+1.0d0), &
               'drag_B spreading: sum_a B_a = B_i', nf)
  END SUBROUTINE test_drag_B_spread_conservation

  ! SimPar@ELDragImplicitWeights = hrz: a small kernel inside one element, centred where some Q2
  ! basis functions are negative. Row-sum weights give a negative nodal drag coefficient there
  ! (the anti-damping of case 11); the HRZ weights are non-negative and keep the total
  ! (EL_DEPOSIT_PARTICLE, argument hrz), and the explicit force keeps the row-sum weights.
  SUBROUTINE test_drag_B_hrz_weights(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: sample(EL_SAMPLE_SIZE), force(3), drag_b, rowsum_min, rowsum_total
    REAL*8, ALLOCATABLE :: rowsum_f(:,:)
    REAL*8, ALLOCATABLE :: rowsum_b(:)
    p%position = (/0.32d0, 0.33d0, 0.34d0/)
    p%radius = 0.02d0                 ! kernel radius 0.1, mostly inside the element [0.25, 0.5]^3
    p%density = 1.0d0
    force = (/0.1d0, -0.2d0, 0.7d0/)
    drag_b = 2.0d0
    val_u = 0.0d0; val_v = 0.0d0; val_w = 0.0d0; val_p = 0.0d0
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), force, el_field_data, drag_b=drag_b)
    ALLOCATE(rowsum_b(SIZE(el_field_data%drag_B_ml)))
    ALLOCATE(rowsum_f(3,SIZE(el_field_data%drag_B_ml)))
    rowsum_b = el_field_data%drag_B_ml
    rowsum_f = el_field_data%force_rhs
    rowsum_min = MINVAL(rowsum_b)
    rowsum_total = SUM(rowsum_b)
    CALL check(rowsum_min.LT.0.0d0, 'drag_B row-sum weights: negative nodal coefficient (the defect)', nf)
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), force, el_field_data, drag_b=drag_b, &
                             hrz=.TRUE.)
    CALL check(MINVAL(el_field_data%drag_B_ml).GE.0.0d0, 'drag_B hrz weights: non-negative', nf)
    CALL check(ABS(SUM(el_field_data%drag_B_ml) - rowsum_total).LE.1.0d-12*ABS(rowsum_total), &
               'drag_B hrz weights: same total as row-sum', nf)
    CALL check(ABS(SUM(el_field_data%drag_B_ml) - drag_b).LE.1.0d-10*drag_b, &
               'drag_B hrz weights: sum_a B_a = B_i', nf)
    CALL check(MAXVAL(ABS(el_field_data%force_rhs - rowsum_f)).EQ.0.0d0, &
               'hrz: the explicit force keeps the row-sum weights (identical deposit)', nf)
    WRITE(*,'(A,ES12.4,A,ES12.4)') '  hrz test: row-sum min B_a ', rowsum_min, &
      ', hrz min B_a ', MINVAL(el_field_data%drag_B_ml)
    ! diagnostic control SimPar@ELDragHrzPartner = sampled (check revived from 77ffd070): the
    ! split partner is deposited with the hrz weights, force_rhs = -split * w_hrz: total -split,
    ! every entry >= 0 for split(3) < 0, the coefficient as without split
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), (/0.0d0, 0.0d0, 0.0d0/), el_field_data, &
                             drag_b=drag_b, hrz=.TRUE., split=(/0.0d0, 0.0d0, -0.7d0/))
    CALL check(ABS(SUM(el_field_data%force_rhs(3,:)) - 0.7d0).LE.1.0d-12 .AND. &
               MINVAL(el_field_data%force_rhs(3,:)).GE.0.0d0, &
               'split partner (ELDragHrzPartner = sampled): deposited with the hrz weights (total -split, one sign)', nf)
    CALL check(ABS(SUM(el_field_data%drag_B_ml) - drag_b).LE.1.0d-10*drag_b .AND. &
               MINVAL(el_field_data%drag_B_ml).GE.0.0d0, &
               'split partner (ELDragHrzPartner = sampled): hrz coefficient unchanged', nf)
    DEALLOCATE(rowsum_b, rowsum_f)
  END SUBROUTINE test_drag_B_hrz_weights

  ! SimPar@ELDragImplicitWeights = hrz (review finding 2): on the regular mesh with a small and a
  ! wide kernel, and on a distorted mesh with a wide kernel,
  !   (a) the nodal coefficient B_i of EL_DEPOSIT_PARTICLE equals an independent HRZ
  !       computation node by node (the row-sum coefficient as well, which checks the
  !       reference), and the partner of EL_APPLY_DRAG_SPLIT_PARTNER is dt B_i u_i node by node,
  !   (b) in a steady non-uniform field the drag terms of one theta step give the fluid exactly
  !       minus the particle drag impulse (job 164663: the former partner -B u_f with the HRZ
  !       weights left B (u_f - sum_i h_i u_i), -1.9 % of F_z here).
  SUBROUTINE test_hrz_node_by_node(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: dcor_saved(3, nvt)
    INTEGER :: ivt, k
    REAL*8 :: x(3)

    p%density = 1.0d0
    p%position = (/0.32d0, 0.33d0, 0.34d0/)
    p%radius = 0.02d0
    CALL hrz_case(p, 'regular mesh, small kernel', nf)
    p%position = (/0.41d0, 0.46d0, 0.53d0/)
    p%radius = 0.03d0
    CALL hrz_case(p, 'regular mesh, wide kernel', nf)

    ! distorted mesh: smooth displacement of the interior vertices (|d| <= 0.04 = h/6.25)
    dcor_saved = mesh%level(1)%dcorvg(:, 1:nvt)
    DO ivt = 1, nvt
      x = mesh%level(1)%dcorvg(:, ivt)
      IF (ANY(x.LT.1.0d-12) .OR. ANY(x.GT.1.0d0-1.0d-12)) CYCLE
      DO k = 1, 3
        mesh%level(1)%dcorvg(k, ivt) = x(k) + 0.04d0*SIN(pi*DBLE(k)*x(MOD(k,3)+1) + 0.3d0)* &
          SIN(2.0d0*pi*x(MOD(k+1,3)+1))
      END DO
    END DO
    CALL build_dof_coordinates(mesh, coord)
    CALL hrz_case(p, 'distorted mesh, wide kernel', nf)
    mesh%level(1)%dcorvg(:, 1:nvt) = dcor_saved
    CALL build_dof_coordinates(mesh, coord)
  END SUBROUTINE test_hrz_node_by_node

  SUBROUTINE hrz_case(p, label, nf)
    TYPE(tElParticleRecord), INTENT(IN) :: p
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nf
    REAL*8, PARAMETER :: dt = 0.125d0, theta = 0.5d0
    REAL*8 :: sample(EL_SAMPLE_SIZE), uf(3), vp(3), force(3), drag_b, fsum(3), err, nref
    REAL*8, ALLOCATABLE :: bh(:), br(:), du(:), dv(:), dw(:), pu(:), pv(:), pw(:)
    INTEGER :: g

    ALLOCATE(bh(ndof), br(ndof), du(ndof), dv(ndof), dw(ndof), pu(ndof), pv(ndof), pw(ndof))
    DO g = 1, ndof
      val_u(g) = 0.1d0 + 0.2d0*coord(1,g)*coord(2,g)
      val_v(g) = -0.3d0*coord(3,g)**2
      val_w(g) = 0.5d0 + coord(1,g) - 0.7d0*coord(2,g)*coord(3,g)
    END DO
    val_p = 0.0d0
    drag_b = 2.0d0
    CALL hrz_reference(p, drag_b, bh, br, nref)

    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL check(ABS(sample(1)/nref - 1.0d0).LE.1.0d-12, 'hrz '//label// &
      ': reference kernel integral = sampler normalisation', nf)
    uf = sample(2:4)/sample(1)
    vp = (/0.0d0, 0.0d0, 0.1d0/)
    force = drag_b*(uf - vp)
    ! row-sum coefficient against the reference (checks the reference itself)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), force, el_field_data, drag_b=drag_b)
    err = MAXVAL(ABS(el_field_data%drag_B_ml - br))/drag_b
    CALL check(err.LE.1.0d-12, 'hrz '//label//': row-sum B_i = independent reference, every node', nf)
    ! (a) HRZ coefficient node by node
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), force, el_field_data, drag_b=drag_b, &
                             hrz=.TRUE.)
    err = MAXVAL(ABS(el_field_data%drag_B_ml - bh))/drag_b
    WRITE(*,'(A,A,A,ES10.2,A,ES11.3,A,ES11.3)') '  hrz ', label, ': max |B_i - ref| / B ', err, &
      ', row-sum min B_i / B ', MINVAL(br)/drag_b, ', hrz min ', MINVAL(bh)/drag_b
    CALL check(err.LE.1.0d-12 .AND. MINVAL(el_field_data%drag_B_ml).GE.0.0d0, &
      'hrz '//label//': B_i = independent HRZ reference, every node, all >= 0', nf)
    CALL EL_CAPTURE_FLUID_FEEDBACK_SOURCE()
    pu = 0.0d0; pv = 0.0d0; pw = 0.0d0
    CALL EL_APPLY_DRAG_SPLIT_PARTNER(pu, pv, pw, val_u, val_v, val_w, dt)
    err = MAX(MAXVAL(ABS(pu - dt*el_field_data%drag_B_source*val_u)), &
              MAXVAL(ABS(pv - dt*el_field_data%drag_B_source*val_v)), &
              MAXVAL(ABS(pw - dt*el_field_data%drag_B_source*val_w)))
    CALL check(err.LE.4.0d0*EPSILON(1.0d0)*dt*drag_b*MAX(MAXVAL(ABS(val_u)), MAXVAL(ABS(val_v)), &
      MAXVAL(ABS(val_w))), 'hrz '//label//': partner = dt B_i u_i^n, every node', nf)
    ! (b) one theta step in a steady field (u^{n+1} = u^n = val): right-hand side dt (source +
    ! partner), drag terms of QuadSc_def theta dt B_i u_i^{n+1} (matrix) and (1-theta) dt B_i u_i^n
    ! (old-velocity defect); their difference per dof is the fluid drag impulse of the dof
    du = 0.0d0; dv = 0.0d0; dw = 0.0d0
    CALL EL_APPLY_FLUID_FEEDBACK_SOURCE(du, dv, dw, dt)
    CALL EL_APPLY_DRAG_SPLIT_PARTNER(du, dv, dw, val_u, val_v, val_w, dt)
    du = du - dt*(theta*el_field_data%drag_B_source*val_u + (1.0d0-theta)*el_field_data%drag_B_source*val_u)
    dv = dv - dt*(theta*el_field_data%drag_B_source*val_v + (1.0d0-theta)*el_field_data%drag_B_source*val_v)
    dw = dw - dt*(theta*el_field_data%drag_B_source*val_w + (1.0d0-theta)*el_field_data%drag_B_source*val_w)
    ! the split correction alone (partner minus the two drag parts of the theta scheme)
    pu = pu - dt*(theta*el_field_data%drag_B_source*val_u + (1.0d0-theta)*el_field_data%drag_B_source*val_u)
    pv = pv - dt*(theta*el_field_data%drag_B_source*val_v + (1.0d0-theta)*el_field_data%drag_B_source*val_v)
    pw = pw - dt*(theta*el_field_data%drag_B_source*val_w + (1.0d0-theta)*el_field_data%drag_B_source*val_w)
    err = MAX(MAXVAL(ABS(pu)), MAXVAL(ABS(pv)), MAXVAL(ABS(pw)))
    CALL check(err.LE.4.0d0*EPSILON(1.0d0)*dt*drag_b*MAX(MAXVAL(ABS(val_u)), MAXVAL(ABS(val_v)), &
      MAXVAL(ABS(val_w))), 'hrz '//label//': steady field, split correction vanishes on every node', nf)
    err = MAX(MAXVAL(ABS(du - dt*el_field_data%fluid_feedback_source(1,:))), &
              MAXVAL(ABS(dv - dt*el_field_data%fluid_feedback_source(2,:))), &
              MAXVAL(ABS(dw - dt*el_field_data%fluid_feedback_source(3,:))))
    CALL check(err.LE.1.0d-14*dt*(MAXVAL(ABS(force)) + drag_b*MAXVAL(ABS(val_w))), &
      'hrz '//label//': steady field, fluid drag force per node = deposited reaction', nf)
    fsum = (/ SUM(du), SUM(dv), SUM(dw) /)/dt
    WRITE(*,'(A,A,A,3ES12.4,A,3ES12.4)') '  hrz ', label, ': fluid drag force + F ', fsum + force, &
      '  F ', force
    CALL check(MAXVAL(ABS(fsum + force)).LE.1.0d-12*MAXVAL(ABS(force)), &
      'hrz '//label//': steady non-uniform field, fluid drag impulse = -particle drag impulse', nf)
    DEALLOCATE(bh, br, du, dv, dw, pu, pv, pw)
  END SUBROUTINE hrz_case

  ! Independent nodal drag coefficients of one particle (drag_b): row-sum br_i = B sum_e int K
  ! phi_i / norm and HRZ bh_i = B sum_e int K phi_i^2 (sum_j int K phi_j)/(sum_j int K phi_j^2) /
  ! norm, norm = int K. Own implementation: Q2 basis as products of 1-D quadratic Lagrange
  ! polynomials at the reference positions of the 27 dofs (vertices in the kvert order, edge
  ! midpoints and face centres as in BUILD_DOF_COORDINATES), own trilinear map, 3 x 3 x 3 Gauss.
  SUBROUTINE hrz_reference(p, drag_b, bh, br, norm)
    TYPE(tElParticleRecord), INTENT(IN) :: p
    REAL*8, INTENT(IN) :: drag_b
    REAL*8, INTENT(OUT) :: bh(:), br(:), norm
    INTEGER, PARAMETER :: ed(2,12) = RESHAPE((/ &
      1,2, 2,3, 4,3, 1,4, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 8,7, 5,8 /), (/2,12/))
    INTEGER, PARAMETER :: fc(4,6) = RESHAPE((/ &
      1,2,3,4,  1,2,6,5,  2,3,7,6,  4,3,7,8,  1,4,8,5,  5,6,7,8 /), (/4,6/))
    REAL*8, PARAMETER :: sg(3,8) = RESHAPE((/ &
      -1.0d0,-1.0d0,-1.0d0,  1.0d0,-1.0d0,-1.0d0,  1.0d0, 1.0d0,-1.0d0, -1.0d0, 1.0d0,-1.0d0, &
      -1.0d0,-1.0d0, 1.0d0,  1.0d0,-1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, -1.0d0, 1.0d0, 1.0d0 /), &
      (/3,8/))
    REAL*8 :: ref(3,27), xv(3,8), gp(3), gw(3), xi(3), xq(3), jac(3,3), detj, w, kern
    REAL*8 :: phi(27), rs(27), sq(27), nv, dn(3), width
    INTEGER :: gl(27), iel, i, a, b, c, k

    gp = (/ -SQRT(0.6d0), 0.0d0, SQRT(0.6d0) /)
    gw = (/ 5.0d0/9.0d0, 8.0d0/9.0d0, 5.0d0/9.0d0 /)
    width = 2.0d0*el_kernel_width_factor*p%radius
    bh = 0.0d0
    br = 0.0d0
    norm = 0.0d0
    DO iel = 1, nel
      DO i = 1, 8
        ref(:,i) = sg(:,i)
        gl(i) = mesh%level(1)%kvert(i,iel)
        xv(:,i) = mesh%level(1)%dcorvg(:, gl(i))
      END DO
      DO i = 1, 12
        ref(:,8+i) = 0.5d0*(sg(:,ed(1,i)) + sg(:,ed(2,i)))
        gl(8+i) = nvt + mesh%level(1)%kedge(i,iel)
      END DO
      DO i = 1, 6
        ref(:,20+i) = 0.25d0*(sg(:,fc(1,i)) + sg(:,fc(2,i)) + sg(:,fc(3,i)) + sg(:,fc(4,i)))
        gl(20+i) = nvt + net + mesh%level(1)%karea(i,iel)
      END DO
      ref(:,27) = 0.0d0
      gl(27) = nvt + net + nat + iel
      rs = 0.0d0
      sq = 0.0d0
      DO c = 1, 3
        DO b = 1, 3
          DO a = 1, 3
            xi = (/ gp(a), gp(b), gp(c) /)
            xq = 0.0d0
            jac = 0.0d0
            DO i = 1, 8
              nv = 0.125d0*(1.0d0 + sg(1,i)*xi(1))*(1.0d0 + sg(2,i)*xi(2))*(1.0d0 + sg(3,i)*xi(3))
              dn(1) = 0.125d0*sg(1,i)*(1.0d0 + sg(2,i)*xi(2))*(1.0d0 + sg(3,i)*xi(3))
              dn(2) = 0.125d0*sg(2,i)*(1.0d0 + sg(1,i)*xi(1))*(1.0d0 + sg(3,i)*xi(3))
              dn(3) = 0.125d0*sg(3,i)*(1.0d0 + sg(1,i)*xi(1))*(1.0d0 + sg(2,i)*xi(2))
              xq = xq + nv*xv(:,i)
              DO k = 1, 3
                jac(:,k) = jac(:,k) + dn(k)*xv(:,i)
              END DO
            END DO
            detj = jac(1,1)*(jac(2,2)*jac(3,3) - jac(3,2)*jac(2,3)) - &
                   jac(2,1)*(jac(1,2)*jac(3,3) - jac(3,2)*jac(1,3)) + &
                   jac(3,1)*(jac(1,2)*jac(2,3) - jac(2,2)*jac(1,3))
            w = gw(a)*gw(b)*gw(c)*ABS(detj)
            kern = EL_KERNEL_VALUE(SQRT(SUM((xq - p%position)**2)), width)
            IF (kern.LE.0.0d0) CYCLE
            DO i = 1, 27
              phi(i) = lag(ref(1,i), xi(1))*lag(ref(2,i), xi(2))*lag(ref(3,i), xi(3))
            END DO
            rs = rs + phi*kern*w
            sq = sq + phi*phi*kern*w
            norm = norm + kern*w
          END DO
        END DO
      END DO
      IF (SUM(sq).LE.0.0d0) CYCLE
      DO i = 1, 27
        br(gl(i)) = br(gl(i)) + drag_b*rs(i)
        bh(gl(i)) = bh(gl(i)) + drag_b*sq(i)*SUM(rs)/SUM(sq)
      END DO
    END DO
    br = br/norm
    bh = bh/norm
  END SUBROUTINE hrz_reference

  ! 1-D quadratic Lagrange polynomial of the node at r (-1, 0, 1), evaluated at t
  PURE REAL*8 FUNCTION lag(r, t)
    REAL*8, INTENT(IN) :: r, t
    IF (r.LT.-0.5d0) THEN
      lag = 0.5d0*t*(t - 1.0d0)
    ELSE IF (r.GT.0.5d0) THEN
      lag = 0.5d0*t*(t + 1.0d0)
    ELSE
      lag = 1.0d0 - t*t
    END IF
  END FUNCTION lag

  SUBROUTINE test_feedback_source_capture(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: sample(EL_SAMPLE_SIZE), feedback(3), expected(3), saved(3)
    REAL*8 :: drag_b, saved_drag_b
    INTEGER :: c
    LOGICAL :: ok
    p = central_particle()
    feedback = (/0.4d0, 0.0d0, -0.9d0/)
    drag_b = 2.5d0
    val_u = 0.0d0; val_v = 0.0d0; val_w = 0.0d0; val_p = 0.0d0
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), feedback, el_field_data, &
                             drag_b=drag_b)
    CALL EL_CAPTURE_FLUID_FEEDBACK_SOURCE()
    ok = .TRUE.
    DO c = 1, 3
      expected(c) = -feedback(c)
      saved(c) = SUM(el_field_data%fluid_feedback_source(c,:))
      IF (ABS(saved(c)-expected(c)).GT.1.0d-10*(ABS(feedback(c))+1.0d0)) &
        ok = .FALSE.
    END DO
    saved_drag_b = SUM(el_field_data%drag_B_source)
    IF (ABS(saved_drag_b-drag_b).GT.1.0d-10*(ABS(drag_b)+1.0d0)) ok = .FALSE.
    CALL check(ok, 'fluid source captures deposited reaction sign', nf)

    ! The post-advance diagnostic pass clears/reuses force_rhs but must not
    ! erase the pre-advance source snapshot used by the fluid solve.
    CALL EL_BEGIN_FIELD_UPDATE(.FALSE.)
    ok = MAXVAL(ABS(el_field_data%force_rhs)).EQ.0.0d0
    DO c = 1, 3
      saved(c) = SUM(el_field_data%fluid_feedback_source(c,:))
      IF (ABS(saved(c)-expected(c)).GT.1.0d-10*(ABS(feedback(c))+1.0d0)) &
        ok = .FALSE.
    END DO
    saved_drag_b = SUM(el_field_data%drag_B_source)
    IF (ABS(saved_drag_b-drag_b).GT.1.0d-10*(ABS(drag_b)+1.0d0)) ok = .FALSE.
    CALL check(ok, 'fluid source survives diagnostic field refresh', nf)
  END SUBROUTINE test_feedback_source_capture

  SUBROUTINE test_feedback_source_apply(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8, ALLOCATABLE :: du(:), dv(:), dw(:)
    REAL*8 :: dt
    LOGICAL :: ok
    INTEGER :: i
    dt = 0.125d0
    ALLOCATE(du(ndof), dv(ndof), dw(ndof))
    du = 1.0d0
    dv = -2.0d0
    dw = 0.5d0
    el_field_data%fluid_feedback_source = 0.0d0
    DO i = 1, ndof
      el_field_data%fluid_feedback_source(1,i) = 0.25d0
      el_field_data%fluid_feedback_source(2,i) = -0.5d0
      el_field_data%fluid_feedback_source(3,i) = 1.0d0
    END DO
    CALL EL_APPLY_FLUID_FEEDBACK_SOURCE(du, dv, dw, dt)
    ok = MAXVAL(ABS(du - (1.0d0 + dt*0.25d0))).LT.1.0d-14
    ok = ok .AND. MAXVAL(ABS(dv - (-2.0d0 - dt*0.5d0))).LT.1.0d-14
    ok = ok .AND. MAXVAL(ABS(dw - (0.5d0 + dt))).LT.1.0d-14
    CALL check(ok, 'fluid source insertion applies exact tstep scaling', nf)
    DEALLOCATE(du, dv, dw)
  END SUBROUTINE test_feedback_source_apply

  SUBROUTINE test_zero_force(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: sample(EL_SAMPLE_SIZE)
    p = central_particle()
    val_u = 0.0d0; val_v = 0.0d0; val_w = 0.0d0; val_p = 0.0d0
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), (/0.0d0,0.0d0,0.0d0/), &
                             el_field_data)
    CALL check(MAXVAL(ABS(el_field_data%force_rhs)).EQ.0.0d0, &
               'zero force produces zero reaction field', nf)
  END SUBROUTINE test_zero_force

  SUBROUTINE test_halo_predicates(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: a(6), b(6)
    LOGICAL :: ok
    a = (/0.0d0,1.0d0, 0.0d0,1.0d0, 0.0d0,1.0d0/)
    b = (/1.2d0,2.0d0, 0.0d0,1.0d0, 0.0d0,1.0d0/)
    ok = .TRUE.
    IF (.NOT.EL_BOXES_WITHIN_DELTA(a, b, 0.3d0)) ok = .FALSE.   ! gap 0.2 <= 0.3
    IF (EL_BOXES_WITHIN_DELTA(a, b, 0.1d0))      ok = .FALSE.   ! gap 0.2 >  0.1
    IF (.NOT.EL_SPHERE_INTERSECTS_BOX((/1.1d0,0.5d0,0.5d0/), 0.2d0, a)) ok = .FALSE.
    IF (EL_SPHERE_INTERSECTS_BOX((/1.4d0,0.5d0,0.5d0/), 0.2d0, a))      ok = .FALSE.
    CALL check(ok, 'halo overlap predicates (box/sphere)', nf)
  END SUBROUTINE test_halo_predicates

  SUBROUTINE test_rank_independence(nf)
    ! Each rank builds an identical mesh + particle and runs the same
    ! transfer math; the conserved quantities must agree across ranks to
    ! floating-point reduction tolerance. This is the assertable part of
    ! the 1/2/8-rank requirement. TODO: extend to a true partition-
    ! straddling particle once the distributed halo path is driven from a
    ! q2p1_el_pipeflow regression harness (see file header).
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: deposited, expected, lo, hi, spread
    deposited = deposit_volume(central_particle(), expected)
    CALL MPI_Allreduce(deposited, lo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                       MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(deposited, hi, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)
    spread = hi - lo
    CALL check(spread.LE.1.0d-12*MAX(1.0d0, expected), &
               'rank-independent deposited volume', nf)
  END SUBROUTINE test_rank_independence

  SUBROUTINE test_partitioned_feedback_no_double_count(nf)
    ! Partitioned-mesh feedback regression (the genuinely distributed case).
    !
    ! The fluid feedback source must be captured from the DISTRIBUTED force_rhs
    ! (before E013Sum3), because the momentum RHS it is added to is distributed
    ! and is summed exactly once by the fluid solver. If it were captured AFTER
    ! E013Sum3 (consistent form), partition-boundary DOFs would be summed twice.
    !
    ! A full FEAT partition cannot be exercised here (E013Sum3 needs the parallel
    ! coupling tables built by the real app), so this models a particle whose
    ! kernel straddles a partition with the REAL EL_CAPTURE_FLUID_FEEDBACK_SOURCE
    ! and field storage:
    !   DOF 1 -> interior to rank 0, DOF 3 -> interior to the last rank,
    !   DOF 2 -> shared by all ranks (each holds a partial S/nproc).
    ! The global reaction summed over unique DOFs is A + S + B == -feedback.
    !
    ! LIMITATIONS (what this test does and does NOT cover):
    !   * It guards the CONVENTION/invariant and the EL_CAPTURE_FLUID_FEEDBACK_
    !     SOURCE semantics: a distributed source sums (once) to -feedback, while a
    !     consistent (pre-summed) source double-counts shared DOFs.
    !   * It does NOT regression-test the production line ordering inside
    !     EL_PARTICLE_MESH_PASS (capture BEFORE E013Sum3). This test builds
    !     force_rhs itself and calls EL_CAPTURE directly, so it would still pass
    !     if the production capture were (wrongly) moved after E013Sum3. The
    !     standalone harness cannot call E013Sum3 because that needs the parallel
    !     coupling tables built by the real app.
    !   * The "partition" is a 3-DOF model, not a real FEAT mesh decomposition:
    !     no actual deposit kernel, element geometry, or shared-DOF topology is
    !     exercised; the shared-DOF partials are prescribed, not produced by
    !     EL_DEPOSIT_PARTICLE straddling a true rank boundary.
    !   * Closing the gap requires an end-to-end partitioned q2p1_el_pipeflow
    !     regression (the standing Phase-1 distributed-straddling TODO), which
    !     drives EL_PARTICLE_MESH_PASS -> E013Sum3 on a real partitioned mesh.
    INTEGER, INTENT(INOUT) :: nf
    INTEGER, PARAMETER :: NL = 3
    REAL*8 :: A(3), B(3), S(3), feedback(3)
    REAL*8 :: local_sum(3), global_sum(3), shared(3), shared_glob(3)
    INTEGER :: c
    LOGICAL :: ok

    A = (/-0.20d0,  0.10d0,  0.05d0/)
    B = (/-0.30d0,  0.40d0, -0.15d0/)
    S = (/-0.50d0, -0.20d0,  0.60d0/)
    feedback = -(A + S + B)

    CALL EL_ENSURE_FIELDS(1, NL)
    el_field_data%force_rhs = 0.0d0
    IF (rank.EQ.0)         el_field_data%force_rhs(:,1) = A
    el_field_data%force_rhs(:,2) = S / DBLE(nproc)
    IF (rank.EQ.nproc-1)   el_field_data%force_rhs(:,3) = B

    ! FIX path: capture distributed force_rhs; the solver single sum (modeled
    ! by the cross-rank Allreduce of local DOF sums) reconstructs -feedback.
    CALL EL_CAPTURE_FLUID_FEEDBACK_SOURCE()
    DO c = 1, 3
      local_sum(c) = SUM(el_field_data%fluid_feedback_source(c,:))
    END DO
    CALL MPI_Allreduce(local_sum, global_sum, 3, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    ok = MAXVAL(ABS(global_sum + feedback)).LT.1.0d-12
    CALL check(ok, 'partitioned feedback: distributed capture sums to -F', nf)

    ! BUG guard (only meaningful when DOFs are actually shared): make the shared
    ! DOF consistent first (an E013Sum), THEN capture. Summed once more by the
    ! solver, the shared DOF is over-counted by (nproc-1)*S, so the global
    ! reaction no longer equals -feedback.
    IF (nproc.GE.2) THEN
      shared = el_field_data%force_rhs(:,2)
      CALL MPI_Allreduce(shared, shared_glob, 3, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      el_field_data%force_rhs(:,2) = shared_glob
      CALL EL_CAPTURE_FLUID_FEEDBACK_SOURCE()
      DO c = 1, 3
        local_sum(c) = SUM(el_field_data%fluid_feedback_source(c,:))
      END DO
      CALL MPI_Allreduce(local_sum, global_sum, 3, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      ok = MAXVAL(ABS(global_sum + feedback)).GT.1.0d-6
      CALL check(ok, &
        'partitioned feedback: consistent capture double-counts (guard)', nf)
    END IF

    ! restore the working-mesh field sizes
    CALL EL_ENSURE_FIELDS(nel, ndof)
  END SUBROUTINE test_partitioned_feedback_no_double_count

  !----------------------------------------------------------------------
  ! Helpers
  !----------------------------------------------------------------------
  FUNCTION central_particle() RESULT(p)
    TYPE(tElParticleRecord) :: p
    p%position = (/0.5d0, 0.5d0, 0.5d0/)
    p%radius   = 0.04d0    ! width = 2*2.5*r = 0.2, support interior
    p%density  = 1.0d0
  END FUNCTION central_particle

  FUNCTION wall_particle() RESULT(p)
    TYPE(tElParticleRecord) :: p
    p%position = (/0.05d0, 0.5d0, 0.5d0/)   ! support reaches past x=0
    p%radius   = 0.04d0
    p%density  = 1.0d0
  END FUNCTION wall_particle

  FUNCTION deposit_volume(p, expected) RESULT(deposited)
    TYPE(tElParticleRecord), INTENT(IN) :: p
    REAL*8, INTENT(OUT) :: expected
    REAL*8 :: deposited
    REAL*8 :: sample(EL_SAMPLE_SIZE)
    val_u = 0.0d0; val_v = 0.0d0; val_w = 0.0d0; val_p = 0.0d0
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), (/0.0d0,0.0d0,0.0d0/), &
                             el_field_data)
    deposited = SUM(el_field_data%alpha_p * mesh%level(1)%dvol)
    expected  = 4.0d0*pi*p%radius**3/3.0d0
  END FUNCTION deposit_volume

  SUBROUTINE check(ok, label, nf)
    LOGICAL, INTENT(IN) :: ok
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nf
    IF (.NOT.ok) THEN
      nf = nf + 1
      IF (rank.EQ.0) WRITE(*,'(A,A)') '  [FAIL] ', label
    ELSE
      IF (rank.EQ.0) WRITE(*,'(A,A)') '  [ ok ] ', label
    END IF
  END SUBROUTINE check

  SUBROUTINE report(nf, r, np)
    INTEGER, INTENT(IN) :: nf, r, np
    INTEGER :: total, e
    CALL MPI_Allreduce(nf, total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, e)
    IF (r.EQ.0) THEN
      IF (total.EQ.0) THEN
        WRITE(*,'(A,I0,A)') 'EL transfer tests passed on ', np, ' rank(s).'
      ELSE
        WRITE(*,'(A,I0,A,I0,A)') 'EL transfer tests FAILED: ', total, &
          ' failure(s) across ', np, ' rank(s).'
      END IF
    END IF
  END SUBROUTINE report

END PROGRAM TEST_EL_TRANSFER
