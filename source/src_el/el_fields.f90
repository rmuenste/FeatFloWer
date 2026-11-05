MODULE EL_FIELDS

  IMPLICIT NONE

  TYPE tElFieldStorage
    INTEGER :: nel = 0
    INTEGER :: ndof = 0
    ! .FALSE. until the first void-fraction finalisation: the initial
    ! epsilon_f_old = 1 would otherwise produce a spurious d eps/dt = -alpha/dt
    ! kick at step 1 (consumed by the volume-averaged continuity source).
    LOGICAL :: history_valid = .FALSE.
    ! BDF2 continuity source (ELVASourceTimeScheme = BDF2): deps_f_dt of the previous
    ! interval, saved when the history advances; rate_valid marks a genuine backward
    ! difference in deps_f_dt (formed over the interval rate_dt), prev_rate_valid and
    ! rate_dt_prev the same for deps_f_dt_prev. Flags are cleared whenever deps_f_dt is
    ! not a genuine difference (first finalisation, dt <= 0).
    LOGICAL :: rate_valid = .FALSE.
    LOGICAL :: prev_rate_valid = .FALSE.
    REAL*8 :: rate_dt = 0.0d0
    REAL*8 :: rate_dt_prev = 0.0d0
    REAL*8, ALLOCATABLE :: deps_f_dt_prev(:)
    REAL*8, ALLOCATABLE :: alpha_p(:)
    REAL*8, ALLOCATABLE :: grad_alpha_p(:,:)   ! cell-averaged kernel gradient of alpha_p (3,nel)
    REAL*8, ALLOCATABLE :: epsilon_f(:)
    REAL*8, ALLOCATABLE :: epsilon_f_old(:)
    REAL*8, ALLOCATABLE :: deps_f_dt(:)
    REAL*8, ALLOCATABLE :: force_rhs(:,:)
    REAL*8, ALLOCATABLE :: fluid_feedback_source(:,:)
    REAL*8, ALLOCATABLE :: drag_B_ml(:)
    REAL*8, ALLOCATABLE :: drag_B_source(:)
    ! Transport porosity of the conservative continuity form (ELVAContinuityForm =
    ! conservative): the Q1 vertex reconstruction of epsilon_f embedded in the Q2 dofs
    ! (EL_VA_FLUX). A reconstruction, not a conserved variable: rebuilt in every transfer
    ! pass after the void-fraction finalisation (cross-rank summed once there), never
    ! restarted. Allocated by EL_ENSURE_TRANSPORT_EPS only under that form.
    REAL*8, ALLOCATABLE :: epsilon_transport_q2(:)
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
    el_field_data%history_valid = .FALSE.
    ALLOCATE(el_field_data%alpha_p(nel))
    ALLOCATE(el_field_data%grad_alpha_p(3,nel))
    ALLOCATE(el_field_data%epsilon_f(nel))
    ALLOCATE(el_field_data%epsilon_f_old(nel))
    ALLOCATE(el_field_data%deps_f_dt(nel))
    ALLOCATE(el_field_data%deps_f_dt_prev(nel))
    el_field_data%deps_f_dt_prev = 0.0d0
    el_field_data%rate_valid = .FALSE.
    el_field_data%prev_rate_valid = .FALSE.
    el_field_data%rate_dt = 0.0d0
    el_field_data%rate_dt_prev = 0.0d0
    ALLOCATE(el_field_data%force_rhs(3,ndof))
    ALLOCATE(el_field_data%fluid_feedback_source(3,ndof))
    ALLOCATE(el_field_data%drag_B_ml(ndof))
    ALLOCATE(el_field_data%drag_B_source(ndof))

    el_field_data%alpha_p = 0.0d0
    el_field_data%grad_alpha_p = 0.0d0
    el_field_data%epsilon_f = 1.0d0
    el_field_data%epsilon_f_old = 1.0d0
    el_field_data%deps_f_dt = 0.0d0
    el_field_data%force_rhs = 0.0d0
    el_field_data%fluid_feedback_source = 0.0d0
    el_field_data%drag_B_ml = 0.0d0
    el_field_data%drag_B_source = 0.0d0

  END SUBROUTINE EL_ENSURE_FIELDS

  ! Allocation of the transport porosity, separate from EL_ENSURE_FIELDS (whose size-based
  ! early return must not skip it when the form is switched on at an existing size).
  SUBROUTINE EL_ENSURE_TRANSPORT_EPS(ndof)

    INTEGER, INTENT(IN) :: ndof

    IF (ALLOCATED(el_field_data%epsilon_transport_q2)) THEN
      IF (SIZE(el_field_data%epsilon_transport_q2).EQ.ndof) RETURN
      DEALLOCATE(el_field_data%epsilon_transport_q2)
    END IF
    ALLOCATE(el_field_data%epsilon_transport_q2(ndof))
    el_field_data%epsilon_transport_q2 = 1.0d0

  END SUBROUTINE EL_ENSURE_TRANSPORT_EPS

  SUBROUTINE EL_BEGIN_FIELD_UPDATE(advance_history)

    LOGICAL, INTENT(IN), OPTIONAL :: advance_history
    LOGICAL :: update_history

    IF (.NOT.ALLOCATED(el_field_data%alpha_p)) RETURN
    update_history = .TRUE.
    IF (PRESENT(advance_history)) update_history = advance_history
    IF (update_history) THEN
      el_field_data%epsilon_f_old = el_field_data%epsilon_f
      ! once per time step (pre-advance pass): keep the finished interval's rate for BDF2
      IF (ALLOCATED(el_field_data%deps_f_dt_prev)) THEN
        el_field_data%deps_f_dt_prev = el_field_data%deps_f_dt
        el_field_data%prev_rate_valid = el_field_data%rate_valid
        el_field_data%rate_dt_prev = el_field_data%rate_dt
      END IF
    END IF
    el_field_data%alpha_p = 0.0d0
    IF (ALLOCATED(el_field_data%grad_alpha_p)) el_field_data%grad_alpha_p = 0.0d0
    el_field_data%force_rhs = 0.0d0
    el_field_data%drag_B_ml = 0.0d0

  END SUBROUTINE EL_BEGIN_FIELD_UPDATE

  SUBROUTINE EL_CAPTURE_FLUID_FEEDBACK_SOURCE()

    IF (.NOT.ALLOCATED(el_field_data%force_rhs)) RETURN
    IF (.NOT.ALLOCATED(el_field_data%fluid_feedback_source)) RETURN
    IF (.NOT.ALLOCATED(el_field_data%drag_B_ml)) RETURN
    IF (.NOT.ALLOCATED(el_field_data%drag_B_source)) RETURN
    el_field_data%fluid_feedback_source = el_field_data%force_rhs
    el_field_data%drag_B_source = el_field_data%drag_B_ml

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
    IF (ALLOCATED(el_field_data%grad_alpha_p)) DEALLOCATE(el_field_data%grad_alpha_p)
    IF (ALLOCATED(el_field_data%epsilon_f)) DEALLOCATE(el_field_data%epsilon_f)
    IF (ALLOCATED(el_field_data%epsilon_f_old)) DEALLOCATE(el_field_data%epsilon_f_old)
    IF (ALLOCATED(el_field_data%deps_f_dt)) DEALLOCATE(el_field_data%deps_f_dt)
    IF (ALLOCATED(el_field_data%deps_f_dt_prev)) DEALLOCATE(el_field_data%deps_f_dt_prev)
    el_field_data%rate_valid = .FALSE.
    el_field_data%prev_rate_valid = .FALSE.
    el_field_data%rate_dt = 0.0d0
    el_field_data%rate_dt_prev = 0.0d0
    IF (ALLOCATED(el_field_data%force_rhs)) DEALLOCATE(el_field_data%force_rhs)
    IF (ALLOCATED(el_field_data%fluid_feedback_source)) &
      DEALLOCATE(el_field_data%fluid_feedback_source)
    IF (ALLOCATED(el_field_data%drag_B_ml)) DEALLOCATE(el_field_data%drag_B_ml)
    IF (ALLOCATED(el_field_data%drag_B_source)) DEALLOCATE(el_field_data%drag_B_source)
    IF (ALLOCATED(el_field_data%epsilon_transport_q2)) &
      DEALLOCATE(el_field_data%epsilon_transport_q2)
    el_field_data%nel = 0
    el_field_data%ndof = 0
    el_field_data%history_valid = .FALSE.

  END SUBROUTINE EL_RELEASE_FIELDS

  SUBROUTINE EL_WRITE_RESTART(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    ! version 3 (2026-10-28): + history/rate validity flags, rate intervals and the
    ! previous-interval rate of the BDF2 source. Older binaries reject version-3 files.
    INTEGER, PARAMETER :: restart_version = 3
    INTEGER :: unit, iflag(3)

    IF (.NOT.ALLOCATED(el_field_data%alpha_p)) RETURN
    OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACCESS='STREAM', &
         FORM='UNFORMATTED', ACTION='WRITE')
    WRITE(unit) restart_version, el_field_data%nel, el_field_data%ndof
    WRITE(unit) el_field_data%alpha_p
    WRITE(unit) el_field_data%epsilon_f
    WRITE(unit) el_field_data%epsilon_f_old
    WRITE(unit) el_field_data%deps_f_dt
    WRITE(unit) el_field_data%force_rhs
    WRITE(unit) el_field_data%drag_B_ml
    WRITE(unit) el_field_data%drag_B_source
    iflag(1) = MERGE(1, 0, el_field_data%history_valid)
    iflag(2) = MERGE(1, 0, el_field_data%rate_valid)
    iflag(3) = MERGE(1, 0, el_field_data%prev_rate_valid)
    WRITE(unit) iflag
    WRITE(unit) el_field_data%rate_dt, el_field_data%rate_dt_prev
    WRITE(unit) el_field_data%deps_f_dt_prev
    CLOSE(unit)

  END SUBROUTINE EL_WRITE_RESTART

  SUBROUTINE EL_READ_RESTART(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, PARAMETER :: restart_version = 3
    INTEGER :: unit, version, nel, ndof, iflag(3)

    OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='OLD', ACCESS='STREAM', &
         FORM='UNFORMATTED', ACTION='READ')
    READ(unit) version, nel, ndof
    IF (version.LT.1 .OR. version.GT.restart_version) THEN
      CLOSE(unit)
      ERROR STOP 'Unsupported E-L restart version.'
    END IF
    CALL EL_ENSURE_FIELDS(nel, ndof)
    READ(unit) el_field_data%alpha_p
    READ(unit) el_field_data%epsilon_f
    READ(unit) el_field_data%epsilon_f_old
    READ(unit) el_field_data%deps_f_dt
    READ(unit) el_field_data%force_rhs
    IF (version.GE.2) THEN
      READ(unit) el_field_data%drag_B_ml
      READ(unit) el_field_data%drag_B_source
    ELSE
      el_field_data%drag_B_ml = 0.0d0
      el_field_data%drag_B_source = 0.0d0
    END IF
    IF (version.GE.3) THEN
      READ(unit) iflag
      el_field_data%history_valid = iflag(1).EQ.1
      el_field_data%rate_valid = iflag(2).EQ.1
      el_field_data%prev_rate_valid = iflag(3).EQ.1
      READ(unit) el_field_data%rate_dt, el_field_data%rate_dt_prev
      READ(unit) el_field_data%deps_f_dt_prev
    ELSE
      ! versions 1-2 carry no validity information: the next finalisation reseeds the
      ! history (as before this version) and the BDF2 source starts with a BE step
      el_field_data%history_valid = .FALSE.
      el_field_data%rate_valid = .FALSE.
      el_field_data%prev_rate_valid = .FALSE.
      el_field_data%rate_dt = 0.0d0
      el_field_data%rate_dt_prev = 0.0d0
      el_field_data%deps_f_dt_prev = 0.0d0
    END IF
    CLOSE(unit)

  END SUBROUTINE EL_READ_RESTART

  SUBROUTINE EL_FINALIZE()

    CALL EL_RELEASE_FIELDS()

  END SUBROUTINE EL_FINALIZE

END MODULE EL_FIELDS
