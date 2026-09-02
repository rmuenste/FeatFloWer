!=========================================================================
! CHI_LEGACY_MESH_ADAPTER - Layer H bridge to the legacy mesh loader
! (design: chimera-integration-design.md v3, section 3; v2-review
! blocker 1 disposition).
!
! readTriCoarse / refineMeshLevel / genMeshStructures (in ff_mesh) USE
! var_QuadScalar and PP3D_MPI, so every call into them lives HERE, in a
! Layer-H compilation unit, and never in the COMMON-free ff_chimera
! library.  This is a temporary bridge by contract: mesh I/O only, never
! assembly.
!
! The loading pipeline mirrors refineMesh (mesh_refine.f90) but inserts
! the Chimera boundary classification + radial projection BETWEEN
! refinement levels (the FeatFloWer treatment of curved boundaries:
! vertices are placed on the true surface on every level, so the next
! refinement bisects projected geometry, not chords):
!
!   readTriCoarse -> structures(1) -> classify level 1 (generator-exact)
!   for l = 1 .. nlmax-1:
!     refineMeshLevel(l -> l+1)
!     propagate the surface bitmask (hierarchical numbering)
!     project new inner/outer-surface vertices radially
!     structures(l+1)
!   re-point coarse dcorvg at the finest array (refineMesh convention)
!   finalize (dof masks, projected Q2 coordinates, boundary face lists)
!=========================================================================
MODULE CHI_LEGACY_MESH_ADAPTER

  USE CHI_SUBMESH, ONLY: tChimeraSubmesh, CHI_CLASSIFY_GEOMETRIC, &
    CHI_PROPAGATE_CLASSIFICATION, CHI_PROJECT_VERTICES, &
    CHI_SUBMESH_FINALIZE_LEVEL
  USE mesh_structures, ONLY: readTriCoarse, refineMeshLevel, &
    genMeshStructures, getNumberOfEdgesOnVerts

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CHI_LOAD_SUBMESH
  PUBLIC :: CHI_RELEASE_SUBMESH_MESH

CONTAINS

  !-----------------------------------------------------------------------
  ! Load the coarse shell .tri, refine to sub%nlmax with per-level
  ! classification/projection, finalize the finest level.  The caller
  ! sets sub%shape/center/radii(/zlo/zhi) and sub%nlmax beforehand.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_LOAD_SUBMESH(sub, trifile, ok)
    TYPE(tChimeraSubmesh), INTENT(INOUT) :: sub
    CHARACTER(*), INTENT(IN) :: trifile
    LOGICAL, INTENT(OUT) :: ok

    INTEGER :: noe, l, nvtc
    INTEGER, ALLOCATABLE :: mask(:), maskc(:)

    ok = .FALSE.
    IF (sub%nlmax .LT. 1) RETURN

    IF (.NOT. ALLOCATED(sub%mesh%level)) THEN
      ALLOCATE(sub%mesh%level(sub%nlmax))
    END IF
    sub%mesh%nlmin = 1
    sub%mesh%nlmax = sub%nlmax
    sub%mesh%maxlevel = sub%nlmax
    sub%mesh%level(1)%nel = 0

    CALL readTriCoarse(trifile, sub%mesh)
    IF (sub%mesh%level(1)%nel .LE. 0) THEN
      WRITE(*,*) 'CHI_LOAD_SUBMESH: could not read ', TRIM(trifile)
      RETURN
    END IF

    CALL getNumberOfEdgesOnVerts(sub%mesh%level(1), noe)
    CALL genMeshStructures(sub%mesh, .FALSE., 1, noe)

    ALLOCATE(mask(sub%mesh%level(1)%nvt))
    CALL CHI_CLASSIFY_GEOMETRIC(sub, sub%mesh%level(1)%dcorvg, &
      sub%mesh%level(1)%nvt, mask)
    CALL CHI_PROJECT_VERTICES(sub, sub%mesh%level(1)%dcorvg, &
      sub%mesh%level(1)%nvt, mask)

    DO l = 1, sub%nlmax - 1
      CALL refineMeshLevel(sub%mesh%level(l), sub%mesh%level(l+1))
      nvtc = sub%mesh%level(l+1)%nvt
      ALLOCATE(maskc(nvtc))
      CALL CHI_PROPAGATE_CLASSIFICATION(mask, sub%mesh%level(l)%nvt, &
        sub%mesh%level(l)%net, sub%mesh%level(l)%nat, &
        sub%mesh%level(l)%nel, sub%mesh%level(l)%kved, &
        sub%mesh%level(l)%kvar, maskc, nvtc)
      CALL CHI_PROJECT_VERTICES(sub, sub%mesh%level(l+1)%dcorvg, &
        nvtc, maskc)
      CALL genMeshStructures(sub%mesh, .FALSE., l+1, noe)
      CALL MOVE_ALLOC(maskc, mask)
    END DO

    ! refineMesh convention: coarse levels share the finest coordinates.
    DO l = 1, sub%nlmax - 1
      DEALLOCATE(sub%mesh%level(l)%dcorvg)
      sub%mesh%level(l)%dcorvg => sub%mesh%level(sub%nlmax)%dcorvg
    END DO

    IF (ALLOCATED(sub%vertmask)) DEALLOCATE(sub%vertmask)
    CALL MOVE_ALLOC(mask, sub%vertmask)

    CALL CHI_SUBMESH_FINALIZE_LEVEL(sub, &
      sub%mesh%level(sub%nlmax)%dcorvg, sub%mesh%level(sub%nlmax)%kvert, &
      sub%mesh%level(sub%nlmax)%kved, sub%mesh%level(sub%nlmax)%kvar, &
      sub%mesh%level(sub%nlmax)%kadj, &
      sub%mesh%level(sub%nlmax)%nvt, sub%mesh%level(sub%nlmax)%net, &
      sub%mesh%level(sub%nlmax)%nat, sub%mesh%level(sub%nlmax)%nel)

    ok = .TRUE.
  END SUBROUTINE CHI_LOAD_SUBMESH

  !-----------------------------------------------------------------------
  ! Release the mesh hierarchy that CHI_LOAD_SUBMESH allocated (the
  ! finest dcorvg once - coarse levels only point at it).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_RELEASE_SUBMESH_MESH(sub)
    TYPE(tChimeraSubmesh), INTENT(INOUT) :: sub
    INTEGER :: l
    IF (.NOT. ALLOCATED(sub%mesh%level)) RETURN
    DO l = 1, SIZE(sub%mesh%level)
      IF (l .LT. SIZE(sub%mesh%level)) THEN
        NULLIFY(sub%mesh%level(l)%dcorvg)
      END IF
    END DO
    IF (ASSOCIATED(sub%mesh%level(SIZE(sub%mesh%level))%dcorvg)) THEN
      DEALLOCATE(sub%mesh%level(SIZE(sub%mesh%level))%dcorvg)
    END IF
    DEALLOCATE(sub%mesh%level)
  END SUBROUTINE CHI_RELEASE_SUBMESH_MESH

END MODULE CHI_LEGACY_MESH_ADAPTER
