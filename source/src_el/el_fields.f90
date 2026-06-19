MODULE EL_FIELDS

  IMPLICIT NONE

  TYPE tElFieldStorage
    INTEGER :: nel = 0
    INTEGER :: ndof = 0
    REAL*8, ALLOCATABLE :: alpha_p(:)
    REAL*8, ALLOCATABLE :: epsilon_f(:)
    REAL*8, ALLOCATABLE :: epsilon_f_old(:)
    REAL*8, ALLOCATABLE :: deps_f_dt(:)
    REAL*8, ALLOCATABLE :: force_rhs(:,:)
    REAL*8, ALLOCATABLE :: fluid_feedback_source(:,:)
  END TYPE tElFieldStorage

  TYPE(tElFieldStorage), SAVE :: el_field_data

CONTAINS

  SUBROUTINE EL_INITIALIZE(nel, ndof)

    INTEGER, INTENT(IN) :: nel, ndof

    CALL EL_ENSURE_FIELDS(nel, ndof)

  END SUBROUTINE EL_INITIALIZE

  SUBROUTINE EL_ENSURE_FIELDS(nel, ndof)

    INTEGER, INTENT(IN) :: nel, ndof

    IF (el_field_data%nel.EQ.nel .AND. el_field_data%ndof.EQ.ndof) RETURN

    CALL EL_RELEASE_FIELDS()

    el_field_data%nel = nel
    el_field_data%ndof = ndof
    ALLOCATE(el_field_data%alpha_p(nel))
    ALLOCATE(el_field_data%epsilon_f(nel))
    ALLOCATE(el_field_data%epsilon_f_old(nel))
    ALLOCATE(el_field_data%deps_f_dt(nel))
    ALLOCATE(el_field_data%force_rhs(3,ndof))
    ALLOCATE(el_field_data%fluid_feedback_source(3,ndof))

    el_field_data%alpha_p = 0.0d0
    el_field_data%epsilon_f = 1.0d0
    el_field_data%epsilon_f_old = 1.0d0
    el_field_data%deps_f_dt = 0.0d0
    el_field_data%force_rhs = 0.0d0
    el_field_data%fluid_feedback_source = 0.0d0

  END SUBROUTINE EL_ENSURE_FIELDS

  SUBROUTINE EL_BEGIN_FIELD_UPDATE(advance_history)

    LOGICAL, INTENT(IN), OPTIONAL :: advance_history
    LOGICAL :: update_history

    IF (.NOT.ALLOCATED(el_field_data%alpha_p)) RETURN
    update_history = .TRUE.
    IF (PRESENT(advance_history)) update_history = advance_history
    IF (update_history) el_field_data%epsilon_f_old = el_field_data%epsilon_f
    el_field_data%alpha_p = 0.0d0
    el_field_data%force_rhs = 0.0d0

  END SUBROUTINE EL_BEGIN_FIELD_UPDATE

  SUBROUTINE EL_CAPTURE_FLUID_FEEDBACK_SOURCE()

    IF (.NOT.ALLOCATED(el_field_data%force_rhs)) RETURN
    IF (.NOT.ALLOCATED(el_field_data%fluid_feedback_source)) RETURN
    el_field_data%fluid_feedback_source = el_field_data%force_rhs

  END SUBROUTINE EL_CAPTURE_FLUID_FEEDBACK_SOURCE

  SUBROUTINE EL_APPLY_FLUID_FEEDBACK_SOURCE(def_u, def_v, def_w, dt)

    REAL*8, INTENT(INOUT) :: def_u(:), def_v(:), def_w(:)
    REAL*8, INTENT(IN) :: dt
    INTEGER :: ndof

    IF (.NOT.ALLOCATED(el_field_data%fluid_feedback_source)) RETURN
    ndof = SIZE(el_field_data%fluid_feedback_source,2)
    IF (SIZE(def_u).LT.ndof .OR. SIZE(def_v).LT.ndof .OR. SIZE(def_w).LT.ndof) THEN
      ERROR STOP 'E-L fluid feedback source/velocity RHS size mismatch.'
    END IF

    def_u(1:ndof) = def_u(1:ndof) + dt*el_field_data%fluid_feedback_source(1,:)
    def_v(1:ndof) = def_v(1:ndof) + dt*el_field_data%fluid_feedback_source(2,:)
    def_w(1:ndof) = def_w(1:ndof) + dt*el_field_data%fluid_feedback_source(3,:)

  END SUBROUTINE EL_APPLY_FLUID_FEEDBACK_SOURCE

  SUBROUTINE EL_RELEASE_FIELDS()

    IF (ALLOCATED(el_field_data%alpha_p)) DEALLOCATE(el_field_data%alpha_p)
    IF (ALLOCATED(el_field_data%epsilon_f)) DEALLOCATE(el_field_data%epsilon_f)
    IF (ALLOCATED(el_field_data%epsilon_f_old)) DEALLOCATE(el_field_data%epsilon_f_old)
    IF (ALLOCATED(el_field_data%deps_f_dt)) DEALLOCATE(el_field_data%deps_f_dt)
    IF (ALLOCATED(el_field_data%force_rhs)) DEALLOCATE(el_field_data%force_rhs)
    IF (ALLOCATED(el_field_data%fluid_feedback_source)) &
      DEALLOCATE(el_field_data%fluid_feedback_source)
    el_field_data%nel = 0
    el_field_data%ndof = 0

  END SUBROUTINE EL_RELEASE_FIELDS

  SUBROUTINE EL_WRITE_RESTART(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, PARAMETER :: restart_version = 1
    INTEGER :: unit

    IF (.NOT.ALLOCATED(el_field_data%alpha_p)) RETURN
    OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACCESS='STREAM', &
         FORM='UNFORMATTED', ACTION='WRITE')
    WRITE(unit) restart_version, el_field_data%nel, el_field_data%ndof
    WRITE(unit) el_field_data%alpha_p
    WRITE(unit) el_field_data%epsilon_f
    WRITE(unit) el_field_data%epsilon_f_old
    WRITE(unit) el_field_data%deps_f_dt
    WRITE(unit) el_field_data%force_rhs
    CLOSE(unit)

  END SUBROUTINE EL_WRITE_RESTART

  SUBROUTINE EL_READ_RESTART(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, PARAMETER :: restart_version = 1
    INTEGER :: unit, version, nel, ndof

    OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='OLD', ACCESS='STREAM', &
         FORM='UNFORMATTED', ACTION='READ')
    READ(unit) version, nel, ndof
    IF (version.NE.restart_version) THEN
      CLOSE(unit)
      ERROR STOP 'Unsupported E-L restart version.'
    END IF
    CALL EL_ENSURE_FIELDS(nel, ndof)
    READ(unit) el_field_data%alpha_p
    READ(unit) el_field_data%epsilon_f
    READ(unit) el_field_data%epsilon_f_old
    READ(unit) el_field_data%deps_f_dt
    READ(unit) el_field_data%force_rhs
    CLOSE(unit)

  END SUBROUTINE EL_READ_RESTART

  SUBROUTINE EL_FINALIZE()

    CALL EL_RELEASE_FIELDS()

  END SUBROUTINE EL_FINALIZE

END MODULE EL_FIELDS
