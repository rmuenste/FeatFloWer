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
                               EL_APPLY_FLUID_FEEDBACK_SOURCE, EL_FINALIZE
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

  SUBROUTINE test_feedback_source_capture(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tElParticleRecord) :: p
    REAL*8 :: sample(EL_SAMPLE_SIZE), feedback(3), expected(3), saved(3)
    INTEGER :: c
    LOGICAL :: ok
    p = central_particle()
    feedback = (/0.4d0, 0.0d0, -0.9d0/)
    val_u = 0.0d0; val_v = 0.0d0; val_w = 0.0d0; val_p = 0.0d0
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    CALL EL_INTEGRATE_PARTICLE(mesh, 1, val_u, val_v, val_w, val_p, &
                               el_field_data%epsilon_f, p, sample)
    CALL EL_DEPOSIT_PARTICLE(mesh, 1, p, sample(1), feedback, el_field_data)
    CALL EL_CAPTURE_FLUID_FEEDBACK_SOURCE()
    ok = .TRUE.
    DO c = 1, 3
      expected(c) = -feedback(c)
      saved(c) = SUM(el_field_data%fluid_feedback_source(c,:))
      IF (ABS(saved(c)-expected(c)).GT.1.0d-10*(ABS(feedback(c))+1.0d0)) &
        ok = .FALSE.
    END DO
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
