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

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tChiBody
  PUBLIC :: CHI_BODY_CYLINDER_Z, CHI_BODY_SPHERE
  PUBLIC :: CHI_MARK_FREE, CHI_MARK_FRINGE, CHI_MARK_HOLE
  PUBLIC :: CHI_POINT_IN_BODY
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
  PURE LOGICAL FUNCTION CHI_POINT_IN_BODY(body, x)
    TYPE(tChiBody), INTENT(IN) :: body
    REAL*8, INTENT(IN) :: x(3)
    REAL*8 :: d(3), r2
    d = x - body%center
    IF (body%shape .EQ. CHI_BODY_SPHERE) THEN
      r2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
    ELSE
      r2 = d(1)*d(1) + d(2)*d(2)
    END IF
    CHI_POINT_IN_BODY = (r2 .LE. body%radius*body%radius)
  END FUNCTION CHI_POINT_IN_BODY

  !-----------------------------------------------------------------------
  ! Classify the Q2 dofs of a (partition of a) background mesh.  q2coor
  ! holds the Q2 node coordinates (3, nvt+net+nat+nel).  kind/pid are
  ! reset to free first.  Precedence between bodies: hole > fringe; the
  ! body id is the first body that raised the node to its final kind
  ! (unambiguous for non-overlapping atmospheres, the M1 assumption).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_CLASSIFY_MARKERS(nel, nvt, net, nat, kvert, kedge, karea, &
                                  q2coor, nbody, bodies, kind, pid)
    INTEGER, INTENT(IN) :: nel, nvt, net, nat
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*)
    REAL*8,  INTENT(IN) :: q2coor(3,*)
    INTEGER, INTENT(IN) :: nbody
    TYPE(tChiBody), INTENT(IN) :: bodies(*)
    INTEGER, INTENT(OUT) :: kind(*), pid(*)

    INTEGER :: ndof, e, k, i, idx(27), nin
    LOGICAL :: inside(27)

    ndof = nvt + net + nat + nel
    kind(1:ndof) = CHI_MARK_FREE
    pid(1:ndof) = 0

    DO e = 1, nel
      CALL CHI_Q2_DOFMAP(e, kvert, kedge, karea, nvt, net, nat, idx)
      DO k = 1, nbody
        nin = 0
        DO i = 1, 27
          inside(i) = CHI_POINT_IN_BODY(bodies(k), q2coor(:,idx(i)))
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
