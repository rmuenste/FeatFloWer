!=========================================================================
! CHI_SUBMESH - the Chimera atmosphere submesh subsystem, Layer M
! (design: chimera-integration-design.md v3, sections 3 and 7).
!
! A tChimeraSubmesh wraps its OWN tMultiMesh instance (never mg_mesh)
! plus placement, boundary classification and Q2 node coordinates.
! Everything here is COMMON-free and instance-based; the calls into the
! legacy mesh loader/refiner (readTriCoarse/refineMeshLevel/
! genMeshStructures, which USE var_QuadScalar) live in the Layer-H
! adapter chi_legacy_mesh_adapter.f90, which drives the classify/project
! routines of this module between refinement levels.
!
! Boundary classification uses a per-vertex surface BITMASK
! (CHI_SURF_INNER/OUTER/ZLO/ZHI) so that rim nodes lying on two surfaces
! (e.g. inner mantle meets the z=0 face) propagate correctly:
!  - level 1: geometric classification (the meshgen generator places
!    boundary vertices EXACTLY on the analytic surfaces),
!  - level l -> l+1 (hierarchical refinement numbering of
!    refineMeshLevel: new ids nvt+edge, nvt+net+face, nvt+net+nat+elem):
!    edge midpoint mask = IAND of the two kved endpoints, face-center
!    mask = IAND of the four kvar corners, cell centers = 0.
!    Assumption (guaranteed by the structured shell generators): no
!    interior "chord" edge/face has all its vertices on one surface.
!  - after propagation, nodes with the inner/outer bit are radially
!    projected onto the true surface (Q1 vertex projection per level -
!    the FeatFloWer treatment of curved boundaries).
!
! Normal convention (design section 5): boundary-face quadrature returns
! n_out = outward normal of the SUBMESH DOMAIN (the atmosphere). On the
! inner surface this points INTO the particle (n_out = -n_B).
!=========================================================================
MODULE CHI_SUBMESH

  USE types, ONLY: tMultiMesh

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tChimeraSubmesh
  PUBLIC :: CHI_SHAPE_CYLINDER_Z, CHI_SHAPE_SPHERE
  PUBLIC :: CHI_SURF_INNER, CHI_SURF_OUTER, CHI_SURF_ZLO, CHI_SURF_ZHI
  PUBLIC :: CHI_CLASSIFY_GEOMETRIC
  PUBLIC :: CHI_PROPAGATE_CLASSIFICATION
  PUBLIC :: CHI_PROJECT_VERTICES
  PUBLIC :: CHI_SUBMESH_FINALIZE_LEVEL
  PUBLIC :: CHI_SUBMESH_RELEASE

  INTEGER, PARAMETER :: CHI_SHAPE_CYLINDER_Z = 1
  INTEGER, PARAMETER :: CHI_SHAPE_SPHERE     = 2

  ! surface bitmask values
  INTEGER, PARAMETER :: CHI_SURF_INNER = 1
  INTEGER, PARAMETER :: CHI_SURF_OUTER = 2
  INTEGER, PARAMETER :: CHI_SURF_ZLO   = 4
  INTEGER, PARAMETER :: CHI_SURF_ZHI   = 8

  ! FEAT local face -> 4 local corner vertices (karea convention:
  ! bottom(1234), front(1265), right(2376), back(3487), left(4158),
  ! top(5678)).
  INTEGER, PARAMETER :: chi_face_verts(4,6) = RESHAPE((/ &
    1,2,3,4,  1,2,6,5,  2,3,7,6,  3,4,8,7,  4,1,5,8,  5,6,7,8 /), (/4,6/))

  TYPE tChimeraSubmesh
    LOGICAL :: initialized = .FALSE.
    INTEGER :: shape = CHI_SHAPE_CYLINDER_Z
    REAL*8  :: center(3) = 0d0
    REAL*8  :: radius_inner = 0d0, radius_outer = 0d0
    REAL*8  :: zlo = 0d0, zhi = 0d0        ! cylinder shape only
    INTEGER :: nlmax = 0
    TYPE(tMultiMesh) :: mesh               ! own hierarchy instance
    ! finest-level data (filled by CHI_SUBMESH_FINALIZE_LEVEL):
    INTEGER :: ndof = 0                    ! Q2 dofs of the finest level
    INTEGER, ALLOCATABLE :: vertmask(:)    ! per-vertex surface bitmask
    INTEGER, ALLOCATABLE :: dofmask(:)     ! per-Q2-dof surface bitmask
    REAL*8,  ALLOCATABLE :: q2coor(:,:)    ! (3,ndof), boundary dofs projected
    INTEGER, ALLOCATABLE :: innerFaces(:,:) ! (2,n): (iel, local face)
    INTEGER, ALLOCATABLE :: outerFaces(:,:)
  END TYPE tChimeraSubmesh

CONTAINS

  !-----------------------------------------------------------------------
  ! Geometric classification of level-1 vertices (generator-exact
  ! coordinates). tol is absolute.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_CLASSIFY_GEOMETRIC(sub, dcorvg, nvt, mask)
    TYPE(tChimeraSubmesh), INTENT(IN) :: sub
    REAL*8,  INTENT(IN)  :: dcorvg(3,*)
    INTEGER, INTENT(IN)  :: nvt
    INTEGER, INTENT(OUT) :: mask(nvt)

    REAL*8 :: tol, r, d(3)
    INTEGER :: i

    tol = 1d-8*MAX(sub%radius_outer, 1d0)
    mask = 0
    DO i = 1, nvt
      d = dcorvg(:,i) - sub%center
      IF (sub%shape .EQ. CHI_SHAPE_SPHERE) THEN
        r = SQRT(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
      ELSE
        r = SQRT(d(1)*d(1) + d(2)*d(2))
        IF (ABS(dcorvg(3,i) - sub%zlo) .LT. tol) mask(i) = mask(i) + CHI_SURF_ZLO
        IF (ABS(dcorvg(3,i) - sub%zhi) .LT. tol) mask(i) = mask(i) + CHI_SURF_ZHI
      END IF
      IF (ABS(r - sub%radius_inner) .LT. tol) mask(i) = mask(i) + CHI_SURF_INNER
      IF (ABS(r - sub%radius_outer) .LT. tol) mask(i) = mask(i) + CHI_SURF_OUTER
    END DO
  END SUBROUTINE CHI_CLASSIFY_GEOMETRIC

  !-----------------------------------------------------------------------
  ! Propagate the vertex mask from a parent level to its refinement
  ! (hierarchical numbering of refineMeshLevel).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_PROPAGATE_CLASSIFICATION(maskp, nvtp, netp, natp, nelp, &
                                          kved, kvar, maskc, nvtc)
    INTEGER, INTENT(IN)  :: nvtp, netp, natp, nelp, nvtc
    INTEGER, INTENT(IN)  :: maskp(nvtp), kved(2,*), kvar(4,*)
    INTEGER, INTENT(OUT) :: maskc(nvtc)

    INTEGER :: e, a, m

    maskc = 0
    maskc(1:nvtp) = maskp
    DO e = 1, netp
      maskc(nvtp + e) = IAND(maskp(kved(1,e)), maskp(kved(2,e)))
    END DO
    DO a = 1, natp
      m = IAND(maskp(kvar(1,a)), maskp(kvar(2,a)))
      m = IAND(m, maskp(kvar(3,a)))
      m = IAND(m, maskp(kvar(4,a)))
      maskc(nvtp + netp + a) = m
    END DO
    ! cell centers (nvtp+netp+natp+1 .. +nelp) stay 0
    IF (nvtp + netp + natp + nelp .NE. nvtc) THEN
      WRITE(*,*) 'CHI_PROPAGATE_CLASSIFICATION: numbering mismatch', &
        nvtp, netp, natp, nelp, nvtc
      STOP 1
    END IF
  END SUBROUTINE CHI_PROPAGATE_CLASSIFICATION

  !-----------------------------------------------------------------------
  ! Radially project vertices carrying the inner/outer bit onto the true
  ! surface (cylinder: scale x,y about the axis; sphere: scale the full
  ! offset vector). z-face nodes need no projection (planar).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_PROJECT_VERTICES(sub, dcorvg, nvt, mask)
    TYPE(tChimeraSubmesh), INTENT(IN) :: sub
    REAL*8,  INTENT(INOUT) :: dcorvg(3,*)
    INTEGER, INTENT(IN)    :: nvt, mask(nvt)

    INTEGER :: i
    REAL*8 :: rtarget

    DO i = 1, nvt
      IF (IAND(mask(i), CHI_SURF_INNER) .NE. 0) THEN
        rtarget = sub%radius_inner
      ELSE IF (IAND(mask(i), CHI_SURF_OUTER) .NE. 0) THEN
        rtarget = sub%radius_outer
      ELSE
        CYCLE
      END IF
      CALL project_point(sub, dcorvg(:,i), rtarget)
    END DO
  END SUBROUTINE CHI_PROJECT_VERTICES

  SUBROUTINE project_point(sub, p, rtarget)
    TYPE(tChimeraSubmesh), INTENT(IN) :: sub
    REAL*8, INTENT(INOUT) :: p(3)
    REAL*8, INTENT(IN) :: rtarget
    REAL*8 :: d(3), r
    d = p - sub%center
    IF (sub%shape .EQ. CHI_SHAPE_SPHERE) THEN
      r = SQRT(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
      IF (r .GT. 0d0) p = sub%center + d*(rtarget/r)
    ELSE
      r = SQRT(d(1)*d(1) + d(2)*d(2))
      IF (r .GT. 0d0) THEN
        p(1) = sub%center(1) + d(1)*(rtarget/r)
        p(2) = sub%center(2) + d(2)*(rtarget/r)
      END IF
    END IF
  END SUBROUTINE project_point

  !-----------------------------------------------------------------------
  ! Finalize the finest level: per-Q2-dof surface mask, projected Q2
  ! node coordinates, and the inner/outer boundary face lists.
  ! Needs the finest level's kvert/kedge/karea/kved/kvar/kadj and the
  ! finest vertex mask (stored into sub%vertmask by the adapter).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_SUBMESH_FINALIZE_LEVEL(sub, dcorvg, kvert, &
                                        kved, kvar, kadj, nvt, net, nat, nel)
    TYPE(tChimeraSubmesh), INTENT(INOUT) :: sub
    REAL*8,  INTENT(IN) :: dcorvg(3,*)
    INTEGER, INTENT(IN) :: kvert(8,*)
    INTEGER, INTENT(IN) :: kved(2,*), kvar(4,*), kadj(6,*)
    INTEGER, INTENT(IN) :: nvt, net, nat, nel

    INTEGER :: i, e, a, f, m, ninner, nouter, cnt(2)
    REAL*8 :: rtarget

    sub%ndof = nvt + net + nat + nel

    IF (ALLOCATED(sub%dofmask)) DEALLOCATE(sub%dofmask)
    IF (ALLOCATED(sub%q2coor)) DEALLOCATE(sub%q2coor)
    ALLOCATE(sub%dofmask(sub%ndof), sub%q2coor(3,sub%ndof))

    ! --- dof masks: vertices / edge dofs / face dofs / cell dofs -------
    sub%dofmask = 0
    sub%dofmask(1:nvt) = sub%vertmask(1:nvt)
    DO e = 1, net
      sub%dofmask(nvt + e) = IAND(sub%vertmask(kved(1,e)), &
                                  sub%vertmask(kved(2,e)))
    END DO
    DO a = 1, nat
      m = IAND(sub%vertmask(kvar(1,a)), sub%vertmask(kvar(2,a)))
      m = IAND(m, sub%vertmask(kvar(3,a)))
      m = IAND(m, sub%vertmask(kvar(4,a)))
      sub%dofmask(nvt + net + a) = m
    END DO

    ! --- Q2 node coordinates (SetUp_myQ2Coor layout), then project the
    !     boundary dofs onto the true surfaces ---------------------------
    DO i = 1, nvt
      sub%q2coor(:,i) = dcorvg(:,i)
    END DO
    DO e = 1, net
      sub%q2coor(:,nvt+e) = 0.5d0*(dcorvg(:,kved(1,e)) + dcorvg(:,kved(2,e)))
    END DO
    DO a = 1, nat
      sub%q2coor(:,nvt+net+a) = 0.25d0*(dcorvg(:,kvar(1,a)) + &
        dcorvg(:,kvar(2,a)) + dcorvg(:,kvar(3,a)) + dcorvg(:,kvar(4,a)))
    END DO
    DO e = 1, nel
      sub%q2coor(:,nvt+net+nat+e) = 0.125d0*( &
        dcorvg(:,kvert(1,e)) + dcorvg(:,kvert(2,e)) + &
        dcorvg(:,kvert(3,e)) + dcorvg(:,kvert(4,e)) + &
        dcorvg(:,kvert(5,e)) + dcorvg(:,kvert(6,e)) + &
        dcorvg(:,kvert(7,e)) + dcorvg(:,kvert(8,e)))
    END DO
    DO i = 1, sub%ndof
      IF (IAND(sub%dofmask(i), CHI_SURF_INNER) .NE. 0) THEN
        rtarget = sub%radius_inner
      ELSE IF (IAND(sub%dofmask(i), CHI_SURF_OUTER) .NE. 0) THEN
        rtarget = sub%radius_outer
      ELSE
        CYCLE
      END IF
      CALL project_point(sub, sub%q2coor(:,i), rtarget)
    END DO

    ! --- boundary face lists (true boundary faces only: kadj == 0) -----
    cnt = 0
    DO e = 1, nel
      DO f = 1, 6
        IF (kadj(f,e) .NE. 0) CYCLE
        m = face_mask(sub, kvert, e, f)
        IF (IAND(m, CHI_SURF_INNER) .NE. 0) cnt(1) = cnt(1) + 1
        IF (IAND(m, CHI_SURF_OUTER) .NE. 0) cnt(2) = cnt(2) + 1
      END DO
    END DO
    IF (ALLOCATED(sub%innerFaces)) DEALLOCATE(sub%innerFaces)
    IF (ALLOCATED(sub%outerFaces)) DEALLOCATE(sub%outerFaces)
    ALLOCATE(sub%innerFaces(2,MAX(cnt(1),1)), sub%outerFaces(2,MAX(cnt(2),1)))
    ninner = 0
    nouter = 0
    DO e = 1, nel
      DO f = 1, 6
        IF (kadj(f,e) .NE. 0) CYCLE
        m = face_mask(sub, kvert, e, f)
        IF (IAND(m, CHI_SURF_INNER) .NE. 0) THEN
          ninner = ninner + 1
          sub%innerFaces(:,ninner) = (/ e, f /)
        END IF
        IF (IAND(m, CHI_SURF_OUTER) .NE. 0) THEN
          nouter = nouter + 1
          sub%outerFaces(:,nouter) = (/ e, f /)
        END IF
      END DO
    END DO
    ! shrink-to-fit via count bookkeeping: store counts in column 0 style
    ! is not possible; callers use SIZE(...,2) so reallocate exactly.
    CALL shrink_faces(sub%innerFaces, ninner)
    CALL shrink_faces(sub%outerFaces, nouter)

    sub%initialized = .TRUE.
  END SUBROUTINE CHI_SUBMESH_FINALIZE_LEVEL

  INTEGER FUNCTION face_mask(sub, kvert, e, f)
    TYPE(tChimeraSubmesh), INTENT(IN) :: sub
    INTEGER, INTENT(IN) :: kvert(8,*), e, f
    INTEGER :: m, j
    m = sub%vertmask(kvert(chi_face_verts(1,f), e))
    DO j = 2, 4
      m = IAND(m, sub%vertmask(kvert(chi_face_verts(j,f), e)))
    END DO
    face_mask = m
  END FUNCTION face_mask

  SUBROUTINE shrink_faces(list, n)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: list(:,:)
    INTEGER, INTENT(IN) :: n
    INTEGER, ALLOCATABLE :: tmp(:,:)
    ALLOCATE(tmp(2,MAX(n,0)))
    IF (n .GT. 0) tmp(:,1:n) = list(:,1:n)
    CALL MOVE_ALLOC(tmp, list)
  END SUBROUTINE shrink_faces

  !-----------------------------------------------------------------------
  ! Release the finest-level arrays (the tMultiMesh content is released
  ! by the Layer-H adapter, which allocated it).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_SUBMESH_RELEASE(sub)
    TYPE(tChimeraSubmesh), INTENT(INOUT) :: sub
    IF (ALLOCATED(sub%vertmask)) DEALLOCATE(sub%vertmask)
    IF (ALLOCATED(sub%dofmask)) DEALLOCATE(sub%dofmask)
    IF (ALLOCATED(sub%q2coor)) DEALLOCATE(sub%q2coor)
    IF (ALLOCATED(sub%innerFaces)) DEALLOCATE(sub%innerFaces)
    IF (ALLOCATED(sub%outerFaces)) DEALLOCATE(sub%outerFaces)
    sub%initialized = .FALSE.
    sub%ndof = 0
  END SUBROUTINE CHI_SUBMESH_RELEASE

END MODULE CHI_SUBMESH
