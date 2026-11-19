!=========================================================================
! pressure_coarse_key.f90
!
! Values of the deck key Pres@PeriodicCoarseOperator (shared pressure parser,
! GetPresParameters). The key chooses the operator of the coarse pressure solve of the master
! when validated periodicity is active (decision and refusals: QuadSc_pco_operator.f90):
!
!   legacy    the coarse Schur matrix the master assembles on its unpartitioned mesh, without
!             the couplings across periodic planes (the arithmetic of before; replays the
!             certified periodic logs bit for bit)
!   gathered  the operator the workers apply (owned rows of C plus the ghost columns of the
!             parallel part CP), collected into the master numbering
!   wrapped   an assembly on the master with the velocity dofs of each periodic class
!             identified (masses summed over the class, then transformed once)
!
! PCO_KEY_UNSET means the key is absent: gathered when validated periodicity is active, with a
! logged fallback to legacy when the configuration is not supported. A value given in the
! deck is binding: an unsupported configuration is refused instead of falling back.
! On a run without periodicity every value leaves the run as it was.
!
! Extension point: a further operator is a new constant here, a name in PCO_KEY_NAME /
! PCO_KEY_PARSE and a branch in the coordinator (PCO_REFRESH_RAW_OPERATOR).
!=========================================================================
MODULE pressure_coarse_key

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: PCO_KEY_INVALID  = -1
  INTEGER, PARAMETER, PUBLIC :: PCO_KEY_UNSET    = 0
  INTEGER, PARAMETER, PUBLIC :: PCO_KEY_LEGACY   = 1
  INTEGER, PARAMETER, PUBLIC :: PCO_KEY_GATHERED = 2
  INTEGER, PARAMETER, PUBLIC :: PCO_KEY_WRAPPED  = 3

  PUBLIC :: PCO_KEY_PARSE, PCO_KEY_NAME

CONTAINS

  ! Deck text -> key value (case-insensitive); PCO_KEY_INVALID for anything else.
  FUNCTION PCO_KEY_PARSE(text) RESULT(key)
    CHARACTER(len=*), INTENT(IN) :: text
    INTEGER :: key
    CHARACTER(len=32) :: t
    INTEGER :: i, c

    t = ADJUSTL(text)
    DO i = 1, LEN_TRIM(t)
      c = IACHAR(t(i:i))
      IF (c.GE.IACHAR('A') .AND. c.LE.IACHAR('Z')) t(i:i) = ACHAR(c + 32)
    END DO
    SELECT CASE (TRIM(t))
    CASE ('legacy')
      key = PCO_KEY_LEGACY
    CASE ('gathered')
      key = PCO_KEY_GATHERED
    CASE ('wrapped')
      key = PCO_KEY_WRAPPED
    CASE DEFAULT
      key = PCO_KEY_INVALID
    END SELECT
  END FUNCTION PCO_KEY_PARSE

  FUNCTION PCO_KEY_NAME(key) RESULT(name)
    INTEGER, INTENT(IN) :: key
    CHARACTER(len=8) :: name
    SELECT CASE (key)
    CASE (PCO_KEY_UNSET)
      name = 'unset'
    CASE (PCO_KEY_LEGACY)
      name = 'legacy'
    CASE (PCO_KEY_GATHERED)
      name = 'gathered'
    CASE (PCO_KEY_WRAPPED)
      name = 'wrapped'
    CASE DEFAULT
      name = 'invalid'
    END SELECT
  END FUNCTION PCO_KEY_NAME

END MODULE pressure_coarse_key
