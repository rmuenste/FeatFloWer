!=========================================================================
! CHI_PERIODIC - periodic-box geometry for the Chimera component
! (design: chimera-integration-design.md v3, Phase 5 "Arrays +
! periodicity"; Layer M, COMMON-free, instance-based).
!
! FeatFloWer realises periodicity by pairing the opposite-face dofs of a
! box that spans exactly one period per periodic axis (dPeriodicity in
! PP3D_MPI).  For the Chimera geometry every "x - X_k" (point-in-body,
! penalty ramp, donor search in an atmosphere) therefore becomes the
! MINIMUM IMAGE displacement, and every submesh sample point that leaves
! the box (an atmosphere straddling a periodic face) is WRAPPED back into
! the box before the background is evaluated.  Both operations are exact
! identities when no axis is periodic, so the non-periodic paths are
! arithmetically unchanged.
!
! The box is described by its lower corner `lo` and the period `len`
! per axis; `per(d)` switches the axis on.  The type has value semantics
! and is passed explicitly (OPTIONAL in the Layer-M kernels).
!=========================================================================
MODULE CHI_PERIODIC

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tChiPeriodic
  PUBLIC :: CHI_PER_ACTIVE
  PUBLIC :: CHI_PER_DELTA
  PUBLIC :: CHI_PER_WRAP
  PUBLIC :: CHI_PER_DIST

  TYPE tChiPeriodic
    LOGICAL :: per(3) = .FALSE.
    REAL*8  :: len(3) = 0d0
    REAL*8  :: lo(3)  = 0d0
  END TYPE tChiPeriodic

CONTAINS

  PURE LOGICAL FUNCTION CHI_PER_ACTIVE(pb)
    TYPE(tChiPeriodic), INTENT(IN) :: pb
    CHI_PER_ACTIVE = ANY(pb%per)
  END FUNCTION CHI_PER_ACTIVE

  !-----------------------------------------------------------------------
  ! Minimum-image displacement x - c.  Non-periodic axes: plain difference.
  !-----------------------------------------------------------------------
  PURE FUNCTION CHI_PER_DELTA(pb, x, c) RESULT(d)
    TYPE(tChiPeriodic), INTENT(IN) :: pb
    REAL*8, INTENT(IN) :: x(3), c(3)
    REAL*8 :: d(3)
    INTEGER :: i
    d = x - c
    DO i = 1, 3
      IF (pb%per(i) .AND. pb%len(i) .GT. 0d0) THEN
        d(i) = d(i) - pb%len(i)*ANINT(d(i)/pb%len(i))
      END IF
    END DO
  END FUNCTION CHI_PER_DELTA

  !-----------------------------------------------------------------------
  ! Minimum-image distance |x - c|.
  !-----------------------------------------------------------------------
  PURE REAL*8 FUNCTION CHI_PER_DIST(pb, x, c)
    TYPE(tChiPeriodic), INTENT(IN) :: pb
    REAL*8, INTENT(IN) :: x(3), c(3)
    REAL*8 :: d(3)
    d = CHI_PER_DELTA(pb, x, c)
    CHI_PER_DIST = SQRT(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
  END FUNCTION CHI_PER_DIST

  !-----------------------------------------------------------------------
  ! Wrap x into [lo, lo+len) per periodic axis (x = lo+len maps to lo).
  !-----------------------------------------------------------------------
  PURE FUNCTION CHI_PER_WRAP(pb, x) RESULT(xw)
    TYPE(tChiPeriodic), INTENT(IN) :: pb
    REAL*8, INTENT(IN) :: x(3)
    REAL*8 :: xw(3)
    INTEGER :: i
    xw = x
    DO i = 1, 3
      IF (pb%per(i) .AND. pb%len(i) .GT. 0d0) THEN
        xw(i) = x(i) - pb%len(i)*FLOOR((x(i) - pb%lo(i))/pb%len(i))
        ! round-off: a point just below lo wraps to lo+len exactly
        IF (xw(i) .GE. pb%lo(i) + pb%len(i)) xw(i) = pb%lo(i)
      END IF
    END DO
  END FUNCTION CHI_PER_WRAP

END MODULE CHI_PERIODIC
