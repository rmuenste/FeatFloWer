!=========================================================================
! CHI_LOCATOR - instance-based point-to-element locator of the Chimera
! component (design: chimera-integration-design.md v3, Layer M).
!
! A uniform bucket grid over element bounding boxes: each grid cell
! stores (CSR-style) the elements whose padded bbox overlaps it.  A
! query walks the candidate list of the cell containing the point,
! prefilters by bbox and confirms by the Newton inverse map.
!
! Deliberately instance-based (TYPE with BUILD/LOCATE/RELEASE): several
! locators coexist (background partition + one per submesh).  The
! module-singleton OctTreeSearch is left untouched.
!=========================================================================
MODULE CHI_LOCATOR

  USE CHI_GEOMETRY, ONLY: CHI_INVERSE_MAP

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tChimeraLocator
  PUBLIC :: CHI_LOCATOR_BUILD
  PUBLIC :: CHI_LOCATE
  PUBLIC :: CHI_LOCATOR_RELEASE

  TYPE tChimeraLocator
    LOGICAL :: initialized = .FALSE.
    INTEGER :: nel = 0
    INTEGER :: nx = 0, ny = 0, nz = 0
    REAL*8  :: xmin(3) = 0d0, dcell(3) = 1d0
    ! CSR cell -> element lists, cell index c = ix + nx*(iy + ny*iz)
    INTEGER, ALLOCATABLE :: cellStart(:)   ! size nx*ny*nz + 1
    INTEGER, ALLOCATABLE :: cellElems(:)
    ! per-element padded bounding boxes for the prefilter
    REAL*8, ALLOCATABLE :: boxlo(:,:), boxhi(:,:)   ! (3,nel)
  END TYPE tChimeraLocator

  REAL*8, PARAMETER :: chi_box_pad_rel = 1d-8

CONTAINS

  !-----------------------------------------------------------------------
  ! Build the bucket grid for a hexahedral mesh given by dcorvg(3,nvt)
  ! and kvert(8,nel).  Grid resolution targets ~1 element per cell per
  ! direction, proportional to the mesh extent.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_LOCATOR_BUILD(loc, dcorvg, kvert, nel, nvt)
    TYPE(tChimeraLocator), INTENT(INOUT) :: loc
    REAL*8,  INTENT(IN) :: dcorvg(3,*)
    INTEGER, INTENT(IN) :: kvert(8,*), nel, nvt

    REAL*8 :: gmin(3), gmax(3), ext(3), pad, diam
    INTEGER :: iel, iv, d, ncells, ic
    INTEGER :: lo(3), hi(3), ix, iy, iz, c
    INTEGER, ALLOCATABLE :: counts(:), fill(:)
    REAL*8 :: p(3)

    CALL CHI_LOCATOR_RELEASE(loc)
    IF (nel .LE. 0 .OR. nvt .LE. 0) RETURN

    loc%nel = nel
    ALLOCATE(loc%boxlo(3,nel), loc%boxhi(3,nel))

    ! element bboxes and global bounds
    gmin =  HUGE(1d0)
    gmax = -HUGE(1d0)
    DO iel = 1, nel
      loc%boxlo(:,iel) =  HUGE(1d0)
      loc%boxhi(:,iel) = -HUGE(1d0)
      DO iv = 1, 8
        p = dcorvg(:,kvert(iv,iel))
        DO d = 1, 3
          loc%boxlo(d,iel) = MIN(loc%boxlo(d,iel), p(d))
          loc%boxhi(d,iel) = MAX(loc%boxhi(d,iel), p(d))
        END DO
      END DO
      DO d = 1, 3
        gmin(d) = MIN(gmin(d), loc%boxlo(d,iel))
        gmax(d) = MAX(gmax(d), loc%boxhi(d,iel))
      END DO
    END DO

    ext = gmax - gmin
    diam = MAX(ext(1), ext(2), ext(3))
    pad = chi_box_pad_rel*MAX(diam, 1d0)
    DO iel = 1, nel
      loc%boxlo(:,iel) = loc%boxlo(:,iel) - pad
      loc%boxhi(:,iel) = loc%boxhi(:,iel) + pad
    END DO
    loc%xmin = gmin - pad

    ! grid resolution: cbrt(nel) cells scaled by relative extent
    DO d = 1, 3
      IF (ext(d) .LE. 0d0) THEN
        lo(d) = 1
      ELSE
        lo(d) = MAX(1, INT(DBLE(nel)**(1d0/3d0)*ext(d)/diam + 0.5d0))
      END IF
    END DO
    loc%nx = lo(1)
    loc%ny = lo(2)
    loc%nz = lo(3)
    DO d = 1, 3
      IF (ext(d) .LE. 0d0) THEN
        loc%dcell(d) = 1d0
      ELSE
        loc%dcell(d) = (ext(d) + 2d0*pad)/DBLE(lo(d))
      END IF
    END DO

    ncells = loc%nx*loc%ny*loc%nz
    ALLOCATE(counts(ncells))
    counts = 0

    ! pass 1: count element-cell overlaps
    DO iel = 1, nel
      CALL cell_range(loc, loc%boxlo(:,iel), loc%boxhi(:,iel), lo, hi)
      DO iz = lo(3), hi(3)
        DO iy = lo(2), hi(2)
          DO ix = lo(1), hi(1)
            c = 1 + ix + loc%nx*(iy + loc%ny*iz)
            counts(c) = counts(c) + 1
          END DO
        END DO
      END DO
    END DO

    ALLOCATE(loc%cellStart(ncells+1))
    loc%cellStart(1) = 1
    DO ic = 1, ncells
      loc%cellStart(ic+1) = loc%cellStart(ic) + counts(ic)
    END DO
    ALLOCATE(loc%cellElems(loc%cellStart(ncells+1)-1))
    ALLOCATE(fill(ncells))
    fill = 0

    ! pass 2: fill
    DO iel = 1, nel
      CALL cell_range(loc, loc%boxlo(:,iel), loc%boxhi(:,iel), lo, hi)
      DO iz = lo(3), hi(3)
        DO iy = lo(2), hi(2)
          DO ix = lo(1), hi(1)
            c = 1 + ix + loc%nx*(iy + loc%ny*iz)
            loc%cellElems(loc%cellStart(c) + fill(c)) = iel
            fill(c) = fill(c) + 1
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE(counts, fill)
    loc%initialized = .TRUE.
  END SUBROUTINE CHI_LOCATOR_BUILD

  !-----------------------------------------------------------------------
  ! Locate physical point p: returns the containing element iel and its
  ! reference coordinates xi, or found = .FALSE. (point outside the mesh
  ! or outside the grid).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_LOCATE(loc, dcorvg, kvert, p, iel, xi, found)
    TYPE(tChimeraLocator), INTENT(IN) :: loc
    REAL*8,  INTENT(IN)  :: dcorvg(3,*), p(3)
    INTEGER, INTENT(IN)  :: kvert(8,*)
    INTEGER, INTENT(OUT) :: iel
    REAL*8,  INTENT(OUT) :: xi(3)
    LOGICAL, INTENT(OUT) :: found

    INTEGER :: ix, iy, iz, c, k, je, iv
    REAL*8 :: nodes(3,8)
    LOGICAL :: conv

    iel = 0
    xi = 0d0
    found = .FALSE.
    IF (.NOT. loc%initialized) RETURN

    ix = INT((p(1) - loc%xmin(1))/loc%dcell(1))
    iy = INT((p(2) - loc%xmin(2))/loc%dcell(2))
    iz = INT((p(3) - loc%xmin(3))/loc%dcell(3))
    IF (ix .LT. 0 .OR. ix .GE. loc%nx) RETURN
    IF (iy .LT. 0 .OR. iy .GE. loc%ny) RETURN
    IF (iz .LT. 0 .OR. iz .GE. loc%nz) RETURN

    c = 1 + ix + loc%nx*(iy + loc%ny*iz)
    DO k = loc%cellStart(c), loc%cellStart(c+1)-1
      je = loc%cellElems(k)
      IF (p(1) .LT. loc%boxlo(1,je) .OR. p(1) .GT. loc%boxhi(1,je)) CYCLE
      IF (p(2) .LT. loc%boxlo(2,je) .OR. p(2) .GT. loc%boxhi(2,je)) CYCLE
      IF (p(3) .LT. loc%boxlo(3,je) .OR. p(3) .GT. loc%boxhi(3,je)) CYCLE
      DO iv = 1, 8
        nodes(:,iv) = dcorvg(:,kvert(iv,je))
      END DO
      CALL CHI_INVERSE_MAP(nodes, p, xi, conv)
      IF (conv) THEN
        iel = je
        found = .TRUE.
        RETURN
      END IF
    END DO
  END SUBROUTINE CHI_LOCATE

  SUBROUTINE CHI_LOCATOR_RELEASE(loc)
    TYPE(tChimeraLocator), INTENT(INOUT) :: loc
    IF (ALLOCATED(loc%cellStart)) DEALLOCATE(loc%cellStart)
    IF (ALLOCATED(loc%cellElems)) DEALLOCATE(loc%cellElems)
    IF (ALLOCATED(loc%boxlo)) DEALLOCATE(loc%boxlo)
    IF (ALLOCATED(loc%boxhi)) DEALLOCATE(loc%boxhi)
    loc%initialized = .FALSE.
    loc%nel = 0
    loc%nx = 0
    loc%ny = 0
    loc%nz = 0
  END SUBROUTINE CHI_LOCATOR_RELEASE

  !-----------------------------------------------------------------------
  ! Grid-cell index range overlapped by a bbox (clamped to the grid).
  !-----------------------------------------------------------------------
  SUBROUTINE cell_range(loc, blo, bhi, lo, hi)
    TYPE(tChimeraLocator), INTENT(IN) :: loc
    REAL*8,  INTENT(IN)  :: blo(3), bhi(3)
    INTEGER, INTENT(OUT) :: lo(3), hi(3)

    INTEGER :: d, n(3)

    n = (/ loc%nx, loc%ny, loc%nz /)
    DO d = 1, 3
      lo(d) = MAX(0, MIN(n(d)-1, INT((blo(d) - loc%xmin(d))/loc%dcell(d))))
      hi(d) = MAX(0, MIN(n(d)-1, INT((bhi(d) - loc%xmin(d))/loc%dcell(d))))
    END DO
  END SUBROUTINE cell_range

END MODULE CHI_LOCATOR
