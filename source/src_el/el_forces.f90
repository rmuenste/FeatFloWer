MODULE EL_FORCES

  USE EL_CONFIG, ONLY: el_drag_model

  IMPLICIT NONE

CONTAINS

  SUBROUTINE EL_DRAG_FORCE(state, carrier_velocity, density, viscosity, force)

    REAL*8, INTENT(IN) :: state(8), carrier_velocity(3)
    REAL*8, INTENT(IN) :: density, viscosity
    REAL*8, INTENT(OUT) :: force(3)

    REAL*8 :: diameter, slip(3), slip_norm, reynolds, correction
    REAL*8 :: dynamic_viscosity

    force = 0.0d0
    diameter = 2.0d0*state(7)
    dynamic_viscosity = density*viscosity
    IF (diameter.LE.0.0d0 .OR. dynamic_viscosity.LE.0.0d0) RETURN

    slip = carrier_velocity - state(4:6)
    slip_norm = SQRT(SUM(slip*slip))
    reynolds = density*diameter*slip_norm/dynamic_viscosity

    correction = 1.0d0
    IF (TRIM(el_drag_model).EQ.'schiller_naumann') THEN
      IF (reynolds.GT.0.0d0 .AND. reynolds.LT.1000.0d0) THEN
        correction = 1.0d0 + 0.15d0*reynolds**0.687d0
      ELSE IF (reynolds.GE.1000.0d0) THEN
        correction = 0.44d0*reynolds/24.0d0
      END IF
    END IF

    force = 3.0d0*ACOS(-1.0d0)*dynamic_viscosity*diameter*correction*slip

  END SUBROUTINE EL_DRAG_FORCE

END MODULE EL_FORCES
