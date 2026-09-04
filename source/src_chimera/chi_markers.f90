!=========================================================================
! CHI_MARKERS - Chimera-S hole/fringe classification kernel, Layer M
! (design: chimera-integration-design.md v3, section 3 "chi_coupling",
! marker representation; paper Section 6).
!
! Two-array marker scheme (never a signed encoding):
!   kind(i) = 0 free, 1 fringe, 2 hole      pid(i) = body id where kind>0
! Normative definition (paper Section 6): a background node x_i is a
! HOLE node if x_i lies in the body B_k; for a background cell crossed
! by the body surface (some of its 27 Q2 nodes inside, some outside)
! the outside nodes are FRINGE nodes.  Hole > fringe > free precedence
! is a numeric MAX on kind - which is exactly what the parallel
! E013Max-style synchronisation computes on partition interfaces; the
! body id follows the winning kind (CHI_MERGE_MARKER).
!
! Pure geometry: bodies are analytic (z-aligned cylinder, treated as
! infinite along z in M1; sphere).  COMMON-free, var_QuadScalar-free.
!=========================================================================
MODULE CHI_MARKERS

  USE CHI_FEM_EVAL, ONLY: CHI_Q2_DOFMAP
  USE CHI_PERIODIC, ONLY: tChiPeriodic, CHI_PER_ACTIVE, CHI_PER_DELTA

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tChiBody
  PUBLIC :: CHI_BODY_CYLINDER_Z, CHI_BODY_SPHERE
  PUBLIC :: CHI_MARK_FREE, CHI_MARK_FRINGE, CHI_MARK_HOLE
  PUBLIC :: CHI_POINT_IN_BODY
  PUBLIC :: CHI_BODY_NEAR_BOX
  PUBLIC :: CHI_CLASSIFY_MARKERS
  PUBLIC :: CHI_MERGE_MARKER

  INTEGER, PARAMETER :: CHI_BODY_CYLINDER_Z = 1
  INTEGER, PARAMETER :: CHI_BODY_SPHERE     = 2

  INTEGER, PARAMETER :: CHI_MARK_FREE   = 0
  INTEGER, PARAMETER :: CHI_MARK_FRINGE = 1
  INTEGER, PARAMETER :: CHI_MARK_HOLE   = 2

  TYPE tChiBody
    INTEGER :: shape = CHI_BODY_CYLINDER_Z
    REAL*8  :: center(3) = 0d0
    REAL*8  :: radius = 0d0
  END TYPE tChiBody

CONTAINS

  !-----------------------------------------------------------------------
  ! Point-in-body test (closed body: surface points count as inside).
  !-----------------------------------------------------------------------
  PURE LOGICAL FUNCTION CHI_POINT_IN_BODY(body, x, pb)
    TYPE(tChiBody), INTENT(IN) :: body
    REAL*8, INTENT(IN) :: x(3)
    TYPE(tChiPeriodic), INTENT(IN), OPTIONAL :: pb   ! Phase 5: minimum image
    REAL*8 :: d(3), r2
    IF (PRESENT(pb)) THEN
      d = CHI_PER_DELTA(pb, x, body%center)
    ELSE
      d = x - body%center
    END IF
    IF (body%shape .EQ. CHI_BODY_SPHERE) THEN
      r2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
    ELSE
      r2 = d(1)*d(1) + d(2)*d(2)
    END IF
    CHI_POINT_IN_BODY = (r2 .LE. body%radius*body%radius)
  END FUNCTION CHI_POINT_IN_BODY

  !-----------------------------------------------------------------------
  ! Conservative box test: can a body's neighbourhood of radius `reach`
  ! (radius, or radius + atmosphere width) touch the axis-aligned box
  ! [lo, hi]?  Without an active periodic box this is the plain
  ! bounding-box overlap test (unchanged Phase-3/4 arithmetic); with one
  ! it uses the minimum-image displacement of the box centre, so periodic
  ! images of the body are seen.  Cylinders ignore the z-axis.
  !-----------------------------------------------------------------------
  PURE LOGICAL FUNCTION CHI_BODY_NEAR_BOX(body, lo, hi, reach, pb)
    TYPE(tChiBody), INTENT(IN) :: body
    REAL*8, INTENT(IN) :: lo(3), hi(3), reach
    TYPE(tChiPeriodic), INTENT(IN), OPTIONAL :: pb
    REAL*8 :: mid(3), half(3), d(3), tol
    LOGICAL :: periodic
    periodic = .FALSE.
    IF (PRESENT(pb)) periodic = CHI_PER_ACTIVE(pb)
    IF (periodic) THEN
      mid = 0.5d0*(lo + hi)
      half = 0.5d0*(hi - lo)
      d = CHI_PER_DELTA(pb, mid, body%center)
      tol = reach*(1d0 + 1d-9)
      IF (body%shape .EQ. CHI_BODY_CYLINDER_Z) THEN
        CHI_BODY_NEAR_BOX = (ABS(d(1)) .LE. tol + half(1)) .AND. &
                            (ABS(d(2)) .LE. tol + half(2))
      ELSE
        CHI_BODY_NEAR_BOX = ALL(ABS(d) .LE. tol + half)
      END IF
    ELSE
      IF (body%shape .EQ. CHI_BODY_CYLINDER_Z) THEN
        CHI_BODY_NEAR_BOX = (lo(1) .LE. body%center(1)+reach .AND. &
                             hi(1) .GE. body%center(1)-reach .AND. &
                             lo(2) .LE. body%center(2)+reach .AND. &
                             hi(2) .GE. body%center(2)-reach)
      ELSE
        CHI_BODY_NEAR_BOX = ALL(lo .LE. body%center+reach) .AND. &
                            ALL(hi .GE. body%center-reach)
      END IF
    END IF
  END FUNCTION CHI_BODY_NEAR_BOX

  !-----------------------------------------------------------------------
  ! Classify the Q2 dofs of a (partition of a) background mesh.  q2coor
  ! holds the Q2 node coordinates (3, nvt+net+nat+nel).  kind/pid are
  ! reset to free first.  Precedence between bodies: hole > fringe; the
  ! body id is the first body that raised the node to its final kind
  ! (unambiguous for non-overlapping atmospheres, the M1 assumption).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_CLASSIFY_MARKERS(nel, nvt, net, nat, kvert, kedge, karea, &
                                  q2coor, nbody, bodies, kind, pid, pb)
    INTEGER, INTENT(IN) :: nel, nvt, net, nat
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*)
    REAL*8,  INTENT(IN) :: q2coor(3,*)
    INTEGER, INTENT(IN) :: nbody
    TYPE(tChiBody), INTENT(IN) :: bodies(*)
    INTEGER, INTENT(OUT) :: kind(*), pid(*)
    TYPE(tChiPeriodic), INTENT(IN), OPTIONAL :: pb   ! Phase 5: periodic box

    INTEGER :: ndof, e, k, i, idx(27), nin
    LOGICAL :: inside(27)
    REAL*8 :: lo(3), hi(3)

    ndof = nvt + net + nat + nel
    kind(1:ndof) = CHI_MARK_FREE
    pid(1:ndof) = 0

    DO e = 1, nel
      CALL CHI_Q2_DOFMAP(e, kvert, kedge, karea, nvt, net, nat, idx)
      lo = q2coor(:,idx(1))
      hi = lo
      DO i = 2, 8
        lo = MIN(lo, q2coor(:,idx(i)))
        hi = MAX(hi, q2coor(:,idx(i)))
      END DO
      DO k = 1, nbody
        ! bounding-box reject (arrays: O(nel) per body instead of 27x)
        IF (.NOT. CHI_BODY_NEAR_BOX(bodies(k), lo, hi, bodies(k)%radius, pb)) CYCLE
        nin = 0
        DO i = 1, 27
          inside(i) = CHI_POINT_IN_BODY(bodies(k), q2coor(:,idx(i)), pb)
          IF (inside(i)) nin = nin + 1
        END DO
        IF (nin .EQ. 0) CYCLE
        DO i = 1, 27
          IF (inside(i)) THEN
            CALL raise_marker(kind(idx(i)), pid(idx(i)), CHI_MARK_HOLE, k)
          ELSE IF (nin .LT. 27) THEN
            CALL raise_marker(kind(idx(i)), pid(idx(i)), CHI_MARK_FRINGE, k)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE CHI_CLASSIFY_MARKERS

  SUBROUTINE raise_marker(kind, pid, newkind, newpid)
    INTEGER, INTENT(INOUT) :: kind, pid
    INTEGER, INTENT(IN) :: newkind, newpid
    IF (newkind .GT. kind) THEN
      kind = newkind
      pid = newpid
    END IF
  END SUBROUTINE raise_marker

  !-----------------------------------------------------------------------
  ! Merge rule for a node shared between partitions (the semantics the
  ! two-pass parallel synchronisation implements): kind = MAX, pid
  ! follows the winning kind; on a tie the ids agree under the
  ! non-overlap assumption, otherwise the smaller id wins
  ! deterministically.
  !-----------------------------------------------------------------------
  PURE SUBROUTINE CHI_MERGE_MARKER(kind, pid, kind_remote, pid_remote)
    INTEGER, INTENT(INOUT) :: kind, pid
    INTEGER, INTENT(IN) :: kind_remote, pid_remote
    IF (kind_remote .GT. kind) THEN
      kind = kind_remote
      pid = pid_remote
    ELSE IF (kind_remote .EQ. kind .AND. kind .GT. 0) THEN
      IF (pid_remote .GT. 0 .AND. (pid .EQ. 0 .OR. pid_remote .LT. pid)) pid = pid_remote
    END IF
  END SUBROUTINE CHI_MERGE_MARKER

END MODULE CHI_MARKERS
