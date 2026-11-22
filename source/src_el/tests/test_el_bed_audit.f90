!=========================================================================
! test_el_bed_audit - helpers of the bed audit line (EL_BED_AUDIT; SimPar@ELBedAuditFreq).
! Case 11; docs/md_docs/el_min_fluidization_prep_2026-11-19.md.
!
! Checks (exit code = the failing check):
!   50-53  sphere-slab volume (full, half, cap, outside, sum of a split)
!   60-62  histogram quantile
!   70-73  plane mean of the P1 pressure on the structured unit-cube mesh (affine field: exact;
!          element-wise discontinuous field: the face-centroid sum)
!   74-75  window gradient -<dp/dz>: affine field over a window that cuts elements; element-
!          wise gradients weighted with the cut volumes
!   80-81  window margin MAX(2 d, kernel support radius)
!   82-85  floor_bal_z: forces armed (balance), force-off contact-free fixture (stationary
!          particles, no pe force, the CFD side computed the buoyant weight but did not arm it:
!          unavailable, not an upward floor force), balance of another step, dt = 0
!   86-89  kernel-sampled window means: recorded samples, bin-centre selection at both window
!          ends, sums; unavailable value is NaN
!   92     ELDragImplicitWeights = hrz needs the drag diagonal in the momentum defect (review r2
!          finding 4): refused for FlowType = non-Newtonian with MatrixRenewal S = 0 only
!=========================================================================
PROGRAM TEST_EL_BED_AUDIT

  USE TYPES, ONLY: tMultiMesh
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
  USE EL_BED_AUDIT, ONLY: EL_SPHERE_SLAB_VOLUME, EL_HISTOGRAM_QUANTILE, EL_PLANE_PRESSURE_LOCAL, &
                          EL_WINDOW_GRADIENT_LOCAL, EL_BED_WINDOW_MARGIN, EL_BED_FLOOR_BALANCE, &
                          EL_BED_WINDOW_SAMPLE_SUMS, EL_BED_UNAVAILABLE, &
                          EL_BED_AUDIT_BEGIN_SAMPLES, EL_BED_AUDIT_RECORD_SAMPLE
  USE EL_CONFIG, ONLY: el_kernel_width_factor, EL_DRAG_HRZ_DEFECT_SUPPORTED
  USE EL_GEOMETRY, ONLY: EL_SET_DOMAIN_BOX
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, BUILD_STRUCTURED_MESH, EL_TEST_RELEASE_MESH

  IMPLICIT NONE

  TYPE(tMultiMesh) :: mesh
  INTEGER :: hist(10), iel, n, i, j
  REAL*8 :: pi, r, v, area, pint, a0, b(3), c(3), expected, sums(3), wz, dtt
  REAL*8, ALLOCATABLE :: p(:)

  pi = ACOS(-1.0d0)

  ! ---- 50-53: sphere-slab volume ------------------------------------------------------------
  r = 0.435d-3
  v = 4.0d0*pi*r**3/3.0d0
  IF (ABS(EL_SPHERE_SLAB_VOLUME(0.0d0, r, -1.0d0, 1.0d0)/v - 1.0d0).GT.1.0d-14) STOP 50
  IF (ABS(EL_SPHERE_SLAB_VOLUME(0.0d0, r, 0.0d0, 1.0d0)/(0.5d0*v) - 1.0d0).GT.1.0d-14) STOP 51
  ! cap of height h = r/2: pi h^2 (3r - h)/3
  IF (ABS(EL_SPHERE_SLAB_VOLUME(0.0d0, r, 0.5d0*r, 1.0d0)/ &
          (pi*(0.5d0*r)**2*(3.0d0*r - 0.5d0*r)/3.0d0) - 1.0d0).GT.1.0d-13) STOP 52
  IF (EL_SPHERE_SLAB_VOLUME(0.0d0, r, 2.0d0*r, 1.0d0).NE.0.0d0) STOP 53
  IF (ABS(EL_SPHERE_SLAB_VOLUME(0.1d0, r, 0.0d0, 0.1d0 + 0.3d0*r) + &
          EL_SPHERE_SLAB_VOLUME(0.1d0, r, 0.1d0 + 0.3d0*r, 1.0d0) - v)/v.GT.1.0d-13) STOP 53

  ! ---- 60-62: histogram quantile --------------------------------------------------------------
  hist = 0
  hist(3) = 50
  hist(5) = 50          ! bins of 0.1 over [0, 1]
  IF (ABS(EL_HISTOGRAM_QUANTILE(hist, 0.0d0, 1.0d0, 0.5d0) - 0.3d0).GT.1.0d-14) STOP 60
  IF (ABS(EL_HISTOGRAM_QUANTILE(hist, 0.0d0, 1.0d0, 0.99d0) - 0.498d0).GT.1.0d-14) STOP 61
  IF (ABS(EL_HISTOGRAM_QUANTILE(hist, 0.0d0, 1.0d0, 0.25d0) - 0.25d0).GT.1.0d-14) STOP 62

  ! ---- 70-73: plane pressure on the unit cube --------------------------------------------------
  CALL EL_TEST_INIT()
  n = 4
  CALL BUILD_STRUCTURED_MESH(mesh, n)
  ALLOCATE(p(4*mesh%level(1)%nel))
  a0 = 3.0d0
  b = (/0.7d0, -1.3d0, 2.9d0/)
  DO iel = 1, mesh%level(1)%nel
    c = 0.0d0
    DO i = 1, 8
      c = c + 0.125d0*mesh%level(1)%dcorvg(:, mesh%level(1)%kvert(i, iel))
    END DO
    p(4*(iel-1)+1) = a0 + DOT_PRODUCT(b, c)
    p(4*(iel-1)+2:4*(iel-1)+4) = b
  END DO
  area = 0.0d0
  pint = 0.0d0
  CALL EL_PLANE_PRESSURE_LOCAL(mesh, 1, p, 0.0d0, 1.0d-9, area, pint)
  IF (ABS(area - 1.0d0).GT.1.0d-13) STOP 70
  IF (ABS(pint/area - (a0 + 0.5d0*b(1) + 0.5d0*b(2))).GT.1.0d-13) STOP 71
  area = 0.0d0
  pint = 0.0d0
  CALL EL_PLANE_PRESSURE_LOCAL(mesh, 1, p, 1.0d0, 1.0d-9, area, pint)
  IF (ABS(area - 1.0d0).GT.1.0d-13 .OR. &
      ABS(pint/area - (a0 + 0.5d0*b(1) + 0.5d0*b(2) + b(3))).GT.1.0d-13) STOP 72
  ! discontinuous: value i + 10 j in the bottom element (i, j), gradient (1, 2, 3): the face mean
  ! of an element is p0 + grad . (face centroid - centre) = p0 - 3 h/2
  expected = 0.0d0
  DO iel = 1, mesh%level(1)%nel
    i = MOD(iel - 1, n)
    j = MOD((iel - 1)/n, n)
    p(4*(iel-1)+1) = DBLE(i) + 10.0d0*DBLE(j)
    p(4*(iel-1)+2:4*(iel-1)+4) = (/1.0d0, 2.0d0, 3.0d0/)
    IF (iel.LE.n*n) expected = expected + (DBLE(i) + 10.0d0*DBLE(j) - 1.5d0/DBLE(n))/DBLE(n*n)
  END DO
  area = 0.0d0
  pint = 0.0d0
  CALL EL_PLANE_PRESSURE_LOCAL(mesh, 1, p, 0.0d0, 1.0d-9, area, pint)
  IF (ABS(pint/area - expected).GT.1.0d-12) STOP 73
  ! window [0.3, 0.7] on the 4 x 4 x 4 mesh (layers of 0.25): layers k = 1 [0.25, 0.5] and
  ! k = 2 [0.5, 0.75] each contribute a volume 0.2 (cut elements); with the z gradient 3 + 10 k
  ! in layer k the volume mean is (0.2 * 13 + 0.2 * 23)/0.4 = 18, so -<dp/dz> = -18
  DO iel = 1, mesh%level(1)%nel
    p(4*(iel-1)+4) = 3.0d0 + 10.0d0*DBLE((iel - 1)/(n*n))
  END DO
  area = 0.0d0
  pint = 0.0d0
  CALL EL_WINDOW_GRADIENT_LOCAL(mesh, 1, p, 0.3d0, 0.7d0, area, pint)
  IF (ABS(area - 0.4d0).GT.1.0d-14) STOP 74
  IF (ABS(pint/area + 18.0d0).GT.1.0d-12) STOP 75
  CALL EL_TEST_RELEASE_MESH(mesh)

  ! ---- 80-81: window margin ------------------------------------------------------------------
  r = 0.435d-3
  el_kernel_width_factor = 2.5d0          ! support radius 5 r = 2.5 d > 2 d (case 11)
  IF (ABS(EL_BED_WINDOW_MARGIN(r) - 5.0d0*r).GT.1.0d-18) STOP 80
  el_kernel_width_factor = 0.5d0          ! support radius r < 2 d
  IF (ABS(EL_BED_WINDOW_MARGIN(r) - 4.0d0*r).GT.1.0d-18) STOP 81
  el_kernel_width_factor = 2.5d0

  ! ---- 82-85: floor balance ------------------------------------------------------------------
  dtt = 2.5d-4
  wz = -3.0d-3                            ! buoyant bed weight (z), N
  ! forces armed, bed at rest on the floor: dp = 0, expected = dt W -> floor force -W (upward)
  IF (ABS(EL_BED_FLOOR_BALANCE(7, 7, 0.0d0, dtt*wz, dtt, .TRUE.) + wz).GT.1.0d-15) STOP 82
  ! force-off, contact-free fixture: nothing armed (ELApplyParticleForces = No), no pe force, no
  ! contact: dp = 0 while expected still holds the computed weight impulse. The former
  ! (dp - expected)/dt reported an upward floor force -W here; now unavailable.
  IF (.NOT.IEEE_IS_NAN(EL_BED_FLOOR_BALANCE(7, 7, 0.0d0, dtt*wz, dtt, .FALSE.))) STOP 83
  ! the recorded balance belongs to another step, or no step length
  IF (.NOT.IEEE_IS_NAN(EL_BED_FLOOR_BALANCE(6, 7, 0.0d0, dtt*wz, dtt, .TRUE.))) STOP 84
  IF (.NOT.IEEE_IS_NAN(EL_BED_FLOOR_BALANCE(7, 7, 0.0d0, dtt*wz, 0.0d0, .TRUE.))) STOP 85

  ! ---- 86-89: kernel-sampled window means ----------------------------------------------------
  CALL EL_SET_DOMAIN_BOX(0.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0)   ! 4000 bins of 2.5e-4
  CALL EL_BED_AUDIT_BEGIN_SAMPLES()
  CALL EL_BED_AUDIT_RECORD_SAMPLE(0.10d0, 0.50d0, 0.20d0)     ! below the window
  CALL EL_BED_AUDIT_RECORD_SAMPLE(0.50d0, 0.40d0, 0.10d0)
  CALL EL_BED_AUDIT_RECORD_SAMPLE(0.50d0, 0.38d0, 0.12d0)     ! same bin
  CALL EL_BED_AUDIT_RECORD_SAMPLE(0.19995d0, 0.90d0, 9.0d0)   ! bin centre 0.199875 < z_lo
  CALL EL_BED_AUDIT_RECORD_SAMPLE(0.20010d0, 0.36d0, 0.30d0)  ! bin centre 0.200125 >= z_lo
  CALL EL_BED_AUDIT_RECORD_SAMPLE(0.80010d0, 0.80d0, 7.0d0)   ! bin centre 0.800125 > z_hi
  CALL EL_BED_AUDIT_RECORD_SAMPLE(1.50d0, 0.99d0, 5.0d0)      ! above the box: last bin
  CALL EL_BED_WINDOW_SAMPLE_SUMS(0.0d0, 1.0d0, 0.2d0, 0.8d0, sums)
  IF (ABS(sums(3) - 3.0d0).GT.0.0d0) STOP 86
  IF (ABS(sums(1) - (0.40d0 + 0.38d0 + 0.36d0)).GT.1.0d-15) STOP 87
  IF (ABS(sums(2) - (0.10d0 + 0.12d0 + 0.30d0)).GT.1.0d-15) STOP 88
  IF (.NOT.IEEE_IS_NAN(EL_BED_UNAVAILABLE())) STOP 89
  CALL EL_BED_AUDIT_BEGIN_SAMPLES()
  CALL EL_BED_WINDOW_SAMPLE_SUMS(0.0d0, 1.0d0, 0.2d0, 0.8d0, sums)
  IF (ANY(sums.NE.0.0d0)) STOP 89

  ! ---- 92: hrz and the momentum defect --------------------------------------------------------
  IF (.NOT.EL_DRAG_HRZ_DEFECT_SUPPORTED(.FALSE., 0)) STOP 92
  IF (.NOT.EL_DRAG_HRZ_DEFECT_SUPPORTED(.FALSE., 1)) STOP 92
  IF (.NOT.EL_DRAG_HRZ_DEFECT_SUPPORTED(.TRUE., 1)) STOP 92
  IF (EL_DRAG_HRZ_DEFECT_SUPPORTED(.TRUE., 0)) STOP 92

  WRITE(*,'(A)') 'test_el_bed_audit: all checks passed'

END PROGRAM TEST_EL_BED_AUDIT
