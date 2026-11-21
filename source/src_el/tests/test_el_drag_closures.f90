!=========================================================================
! test_el_drag_closures - the porosity closures of the drag table (EL_CONFIG el_drag_names):
! gidaspow, wenyu, ergun (case 11, minimum fluidisation; docs/md_docs/
! el_min_fluidization_prep_2026-11-19.md).
!
! Checks (exit code = the failing check):
!   10-19  per-particle coefficient drag_B = V_p beta/(1-eps) against numbers computed by hand
!          (python, double precision) from the formulas of case spec 11 at eps 0.4 and 0.9,
!          slip 0.5 m/s, the case-11 particle (d 0.87 mm, rho_f 1.0, mu 1.8e-5)
!   20-24  zero-slip limit: finite, equal to the analytic limit, continuous in the slip
!   25     Wen-Yu above Re 1000 (C_d = 0.44)
!   30-33  the switch at eps 0.8: Ergun at 0.8 itself, Wen-Yu above; the jump is printed
!   40-41  fixed-bed identity: n drag_B U/eps / eps = the Ergun pressure gradient at
!          superficial velocity U (model-A drag + pressure-gradient force = model-B drag)
!   50-52  U_mf from weight gradient = drag, gidaspow at eps 0.40 (spec 11 table, column
!          Ergun: 0.210 m/s) and wenyu (column Wen-Yu: 0.170 m/s)
!   60-61  closure table: lookup and validation of the new names, id resolved at validation
!   62     every registered closure has a dispatch case (finite positive coefficient)
!=========================================================================
PROGRAM TEST_EL_DRAG_CLOSURES

  USE EL_CONFIG, ONLY: el_drag_model, el_drag_coupling, el_apply_fluid_feedback, &
                       EL_DRAG_LOOKUP, el_drag_names, EL_DRAG_GIDASPOW, EL_DRAG_WENYU, &
                       EL_DRAG_ERGUN, EL_VALIDATE_CONFIG, EL_RESOLVE_DRAG_CLOSURE, &
                       el_drag_closure_id
  USE EL_FORCES, ONLY: EL_DRAG_CLOSURE

  IMPLICIT NONE

  REAL*8, PARAMETER :: D = 0.87d-3, RHO = 1.0d0, MU = 1.8d-5, RHOP = 860.0d0, G = 9.81d0
  REAL*8, PARAMETER :: NU = MU/RHO
  REAL*8 :: pi, vp, b, b0, b1, ratio0, ratio5, u_mf
  INTEGER :: i

  pi = ACOS(-1.0d0)
  vp = pi*D**3/6.0d0

  ! ---- 10-19: values at eps 0.4 and 0.9, slip 0.5 m/s ----------------------------------
  CALL CHECK('ergun', 0.4d0, 0.5d0, 2.191673209914664d-06, 10)
  CALL CHECK('wenyu', 0.4d0, 0.5d0, 2.8662507350071774d-06, 11)
  CALL CHECK('gidaspow', 0.4d0, 0.5d0, 2.191673209914664d-06, 12)
  CALL CHECK('ergun', 0.9d0, 0.5d0, 4.834322045252142d-07, 13)
  CALL CHECK('wenyu', 0.9d0, 0.5d0, 4.3792679251462653d-07, 14)
  CALL CHECK('gidaspow', 0.9d0, 0.5d0, 4.3792679251462653d-07, 15)

  ! ---- 20-24: zero slip ----------------------------------------------------------------
  CALL CHECK('ergun', 0.4d0, 0.0d0, 1.844900285820606d-06, 20)
  CALL CHECK('wenyu', 0.4d0, 0.0d0, 1.6734159654587396d-06, 21)
  CALL CHECK('wenyu', 0.9d0, 0.0d0, 1.9512828841817462d-07, 22)
  ! analytic limits: Ergun 150 (1-eps) mu V_p/(eps d^2) = 25 pi mu d (1-eps)/eps,
  ! Wen-Yu 3 pi mu d eps^-2.65
  b = DRAG_B('ergun', 0.4d0, 0.0d0)
  IF (ABS(b/(25.0d0*pi*MU*D*0.6d0/0.4d0) - 1.0d0).GT.1.0d-13) STOP 23
  b = DRAG_B('wenyu', 0.9d0, 0.0d0)
  b1 = DRAG_B('wenyu', 0.9d0, 1.0d-12)
  IF (ABS(b/(3.0d0*pi*MU*D*0.9d0**(-2.65d0)) - 1.0d0).GT.1.0d-13) STOP 24
  ! continuity: at slip 1e-12 the coefficient exceeds the limit by exactly the Schiller-Naumann
  ! factor 0.15 Re^0.687 (1e-8 here: the approach to the limit is continuous, not smooth)
  IF (ABS(b1/b - (1.0d0 + 0.15d0*(RHO*0.9d0*1.0d-12*D/MU)**0.687d0)).GT.1.0d-14) STOP 24
  IF (ABS(b1/b - 1.0d0).GT.1.0d-7) STOP 24

  ! ---- 25: Wen-Yu above Re 1000 (slip 30 m/s at eps 0.9: Re 1305) ----------------------
  CALL CHECK('wenyu', 0.9d0, 30.0d0, 4.668444300404826d-06, 25)

  ! ---- 30-33: switch at eps 0.8 ----------------------------------------------------------
  b0 = DRAG_B('gidaspow', 0.8d0, 0.5d0)
  IF (b0.NE.DRAG_B('ergun', 0.8d0, 0.5d0)) STOP 30
  b1 = DRAG_B('gidaspow', 0.8d0 + 1.0d-12, 0.5d0)
  IF (ABS(b1/DRAG_B('wenyu', 0.8d0 + 1.0d-12, 0.5d0) - 1.0d0).GT.0.0d0) STOP 31
  ratio5 = DRAG_B('ergun', 0.8d0, 0.5d0)/DRAG_B('wenyu', 0.8d0, 0.5d0)
  ratio0 = DRAG_B('ergun', 0.8d0, 0.0d0)/DRAG_B('wenyu', 0.8d0, 0.0d0)
  WRITE(*,'(A,F8.4,A,F8.4)') 'gidaspow jump at eps 0.8, Ergun/Wen-Yu: slip 0 ', ratio0, &
    ', slip 0.5 m/s ', ratio5
  ! 6.25 pi mu d against 3 pi mu d 0.8^-2.65 at Re 0
  IF (ABS(ratio0 - 6.25d0/(3.0d0*0.8d0**(-2.65d0))).GT.1.0d-12) STOP 32
  IF (ABS(ratio5 - 1.1426773062166042d0).GT.1.0d-12) STOP 33

  ! ---- 40-41: fixed-bed identity against the Ergun equation -----------------------------
  CALL FIXED_BED('ergun', 0.40d0, 0.10d0, 40)
  CALL FIXED_BED('gidaspow', 0.36d0, 0.25d0, 41)

  ! ---- 50-52: minimum fluidisation velocity ------------------------------------------------
  u_mf = UMF('gidaspow', 0.40d0)
  WRITE(*,'(A,F10.6,A)') 'U_mf gidaspow at eps 0.40: ', u_mf, ' m/s (spec 11 table, Ergun: 0.210)'
  IF (ABS(u_mf - 0.210d0).GT.5.0d-4) STOP 50
  ! the same bisection in python (umf_reference.py logic): 0.21038271148569176
  IF (ABS(u_mf/0.21038271148569176d0 - 1.0d0).GT.1.0d-9) STOP 51
  u_mf = UMF('wenyu', 0.40d0)
  WRITE(*,'(A,F10.6,A)') 'U_mf wenyu    at eps 0.40: ', u_mf, ' m/s (spec 11 table, Wen-Yu: 0.170)'
  IF (ABS(u_mf - 0.170d0).GT.5.0d-4 .OR. ABS(u_mf/0.16972276906909284d0 - 1.0d0).GT.1.0d-9) STOP 52

  ! ---- 60-61: table ------------------------------------------------------------------------
  IF (EL_DRAG_LOOKUP('gidaspow').NE.EL_DRAG_GIDASPOW) STOP 60
  IF (EL_DRAG_LOOKUP('wenyu').NE.EL_DRAG_WENYU) STOP 60
  IF (EL_DRAG_LOOKUP('ergun').NE.EL_DRAG_ERGUN) STOP 60
  IF (EL_DRAG_LOOKUP('wen_yu').NE.-1 .OR. EL_DRAG_LOOKUP(' ').NE.-1) STOP 60
  DO i=LBOUND(el_drag_names,1),UBOUND(el_drag_names,1)
    IF (EL_DRAG_LOOKUP(el_drag_names(i)).NE.i) STOP 60
  END DO
  el_drag_coupling = 'explicit'
  el_apply_fluid_feedback = .FALSE.
  DO i=LBOUND(el_drag_names,1),UBOUND(el_drag_names,1)
    el_drag_model = el_drag_names(i)
    CALL EL_VALIDATE_CONFIG()          ! stops with 1 on an invalid name
    ! the validation resolves the id once (review finding 14)
    IF (el_drag_closure_id.NE.i) STOP 61
  END DO
  el_drag_model = 'GIDASPOW'
  CALL EL_VALIDATE_CONFIG()
  IF (TRIM(el_drag_model).NE.'gidaspow' .OR. el_drag_closure_id.NE.EL_DRAG_GIDASPOW) STOP 61

  ! ---- 62: every registered closure has a dispatch case --------------------------------------
  ! EL_DRAG_CLOSURE stops the run (MPI_Abort) on an id without a CASE, so a closure added to the
  ! table but not to the dispatcher fails here instead of falling back to another law: every
  ! table entry gives a finite positive coefficient at eps 0.5, slip 0.1 m/s.
  DO i=LBOUND(el_drag_names,1),UBOUND(el_drag_names,1)
    b = DRAG_B(TRIM(el_drag_names(i)), 0.5d0, 0.1d0)
    IF (.NOT.(b.GT.0.0d0 .AND. b.LT.HUGE(1.0d0))) STOP 62
  END DO

  WRITE(*,'(A)') 'test_el_drag_closures: all checks passed'

CONTAINS

  ! drag_B of closure `name` for a particle at rest, fluid velocity (0,0,slip), porosity eps
  REAL*8 FUNCTION DRAG_B(name, eps, slip)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL*8, INTENT(IN) :: eps, slip
    REAL*8 :: state(8), force(3)
    state = 0.0d0
    state(7) = 0.5d0*D
    state(8) = RHOP
    el_drag_model = name
    CALL EL_RESOLVE_DRAG_CLOSURE()
    CALL EL_DRAG_CLOSURE(state, (/0.0d0, 0.0d0, slip/), RHO, NU, eps, force, DRAG_B)
    IF (ABS(force(3) - DRAG_B*slip).GT.1.0d-15*ABS(DRAG_B*slip)) STOP 9
  END FUNCTION DRAG_B

  SUBROUTINE CHECK(name, eps, slip, expected, code)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL*8, INTENT(IN) :: eps, slip, expected
    INTEGER, INTENT(IN) :: code
    REAL*8 :: value
    value = DRAG_B(name, eps, slip)
    IF (ABS(value/expected - 1.0d0).GT.1.0d-12) THEN
      WRITE(*,'(A,A,2F8.3,2ES24.16)') 'mismatch ', name, eps, slip, value, expected
      STOP code
    END IF
  END SUBROUTINE CHECK

  ! Fixed bed: particles at rest, superficial velocity u_sup, interstitial slip u_sup/eps.
  ! Model-A drag per unit volume n F_d plus the pressure-gradient force (1-eps)|dp/dz| gives
  ! |dp/dz| = n F_d/eps (fluid balance eps dp/dz = n F_d); compare with the Ergun equation.
  SUBROUTINE FIXED_BED(name, eps, u_sup, code)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL*8, INTENT(IN) :: eps, u_sup
    INTEGER, INTENT(IN) :: code
    REAL*8 :: n, drag_per_volume, gradient, ergun_gradient
    n = (1.0d0 - eps)/vp
    drag_per_volume = n*DRAG_B(name, eps, u_sup/eps)*(u_sup/eps)
    gradient = drag_per_volume/eps
    ergun_gradient = 150.0d0*MU*(1.0d0-eps)**2*u_sup/(eps**3*D**2) + &
                     1.75d0*RHO*(1.0d0-eps)*u_sup**2/(eps**3*D)
    WRITE(*,'(A,A,A,F5.2,A,F6.3,A,ES16.8,A,ES16.8)') 'fixed bed ', name, ' eps ', eps, &
      ' U ', u_sup, ': n F_d/eps ', gradient, '  Ergun ', ergun_gradient
    IF (ABS(gradient/ergun_gradient - 1.0d0).GT.1.0d-12) STOP code
  END SUBROUTINE FIXED_BED

  ! U_mf: (1-eps)(rho_p - rho_f) g = n F_d(U/eps)/eps, bisection as umf_reference.py
  REAL*8 FUNCTION UMF(name, eps)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL*8, INTENT(IN) :: eps
    REAL*8 :: a, c, m, target, gradient
    INTEGER :: k
    target = (1.0d0 - eps)*(RHOP - RHO)*G
    a = 1.0d-6
    c = 10.0d0
    DO k=1,200
      m = 0.5d0*(a + c)
      gradient = (1.0d0 - eps)/vp*DRAG_B(name, eps, m/eps)*(m/eps)/eps
      IF (gradient.LT.target) THEN
        a = m
      ELSE
        c = m
      END IF
    END DO
    UMF = 0.5d0*(a + c)
  END FUNCTION UMF

END PROGRAM TEST_EL_DRAG_CLOSURES
