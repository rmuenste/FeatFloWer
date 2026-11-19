!=========================================================================
! test_el_bed_audit - helpers of the bed audit line (EL_BED_AUDIT; SimPar@ELBedAuditFreq).
! Case 11; docs/md_docs/el_min_fluidization_prep_2026-11-19.md.
!
! Checks (exit code = the failing check):
!   50-53  sphere-slab volume (full, half, cap, outside, sum of a split)
!   60-62  histogram quantile
!   70-73  plane mean of the P1 pressure on the structured unit-cube mesh (affine field: exact;
!          element-wise discontinuous field: the face-centroid sum)
!=========================================================================
PROGRAM TEST_EL_BED_AUDIT

  USE TYPES, ONLY: tMultiMesh
  USE EL_BED_AUDIT, ONLY: EL_SPHERE_SLAB_VOLUME, EL_HISTOGRAM_QUANTILE, EL_PLANE_PRESSURE_LOCAL
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, BUILD_STRUCTURED_MESH, EL_TEST_RELEASE_MESH

  IMPLICIT NONE

  TYPE(tMultiMesh) :: mesh
  INTEGER :: hist(10), iel, n, i, j
  REAL*8 :: pi, r, v, area, pint, a0, b(3), c(3), expected
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
  CALL EL_TEST_RELEASE_MESH(mesh)

  WRITE(*,'(A)') 'test_el_bed_audit: all checks passed'

END PROGRAM TEST_EL_BED_AUDIT
