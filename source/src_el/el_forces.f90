MODULE EL_FORCES

  USE EL_CONFIG, ONLY: el_drag_model, el_eps_f_min, el_pressure_force, &
                       el_lift_model
  USE EL_QUADRATURE, ONLY: EL_SAMPLE_SIZE, EL_SAMPLE_NORMALIZATION, &
                           EL_SAMPLE_UF_BEGIN, EL_SAMPLE_UF_END, &
                           EL_SAMPLE_GRAD_U_BEGIN, EL_SAMPLE_GRAD_U_END, &
                           EL_SAMPLE_PRESSURE, EL_SAMPLE_GRAD_P_BEGIN, &
                           EL_SAMPLE_GRAD_P_END, EL_SAMPLE_EPSILON_F

  IMPLICIT NONE

  TYPE tElForceResult
    REAL*8 :: drag(3) = 0.0d0
    REAL*8 :: pressure(3) = 0.0d0
    REAL*8 :: lift(3) = 0.0d0
    REAL*8 :: grav_buoy(3) = 0.0d0
    REAL*8 :: particle_total(3) = 0.0d0
    REAL*8 :: feedback_force(3) = 0.0d0
    REAL*8 :: drag_B = 0.0d0
    REAL*8 :: u_f(3) = 0.0d0
  END TYPE tElForceResult

CONTAINS

  SUBROUTINE EL_COMPUTE_PARTICLE_FORCES(state, sample, density, viscosity, &
      gravity, result)

    REAL*8, INTENT(IN) :: state(8), sample(EL_SAMPLE_SIZE)
    REAL*8, INTENT(IN) :: density, viscosity, gravity(3)
    TYPE(tElForceResult), INTENT(OUT) :: result

    REAL*8 :: norm, grad_u(3,3), grad_p(3), eps_f, volume

    result = tElForceResult()
    norm = sample(EL_SAMPLE_NORMALIZATION)
    IF (norm.LE.0.0d0) RETURN

    result%u_f = sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END)/norm
    grad_u(1,:) = sample(EL_SAMPLE_GRAD_U_BEGIN:EL_SAMPLE_GRAD_U_BEGIN+2)/norm
    grad_u(2,:) = sample(EL_SAMPLE_GRAD_U_BEGIN+3:EL_SAMPLE_GRAD_U_BEGIN+5)/norm
    grad_u(3,:) = sample(EL_SAMPLE_GRAD_U_BEGIN+6:EL_SAMPLE_GRAD_U_END)/norm
    grad_p = sample(EL_SAMPLE_GRAD_P_BEGIN:EL_SAMPLE_GRAD_P_END)/norm
    eps_f = MIN(1.0d0,MAX(el_eps_f_min,sample(EL_SAMPLE_EPSILON_F)/norm))
    volume = EL_PARTICLE_VOLUME(state(7))

    CALL EL_DRAG_CLOSURE(state, result%u_f, density, viscosity, eps_f, &
                         result%drag, result%drag_B)
    IF (el_pressure_force) result%pressure = -volume*grad_p
    CALL EL_LIFT_CLOSURE(state, result%u_f, grad_u, density, viscosity, &
                         eps_f, result%lift)
    result%grav_buoy = (state(8)-density)*volume*gravity

    result%particle_total = result%drag + result%pressure + result%lift + &
                            result%grav_buoy
    result%feedback_force = result%drag + result%lift

  END SUBROUTINE EL_COMPUTE_PARTICLE_FORCES

  SUBROUTINE EL_DRAG_FORCE(state, carrier_velocity, density, viscosity, force)

    REAL*8, INTENT(IN) :: state(8), carrier_velocity(3)
    REAL*8, INTENT(IN) :: density, viscosity
    REAL*8, INTENT(OUT) :: force(3)
    REAL*8 :: drag_B

    CALL EL_DRAG_CLOSURE(state, carrier_velocity, density, viscosity, 1.0d0, &
                         force, drag_B)

  END SUBROUTINE EL_DRAG_FORCE

  SUBROUTINE EL_DRAG_CLOSURE(state, carrier_velocity, density, viscosity, &
                             epsilon_f, force, drag_B)

    REAL*8, INTENT(IN) :: state(8), carrier_velocity(3)
    REAL*8, INTENT(IN) :: density, viscosity, epsilon_f
    REAL*8, INTENT(OUT) :: force(3), drag_B

    REAL*8 :: diameter, slip(3), slip_norm, reynolds, correction
    REAL*8 :: dynamic_viscosity, eps, cd, area, chi

    force = 0.0d0
    drag_B = 0.0d0
    diameter = 2.0d0*state(7)
    dynamic_viscosity = density*viscosity
    IF (diameter.LE.0.0d0 .OR. dynamic_viscosity.LE.0.0d0) RETURN

    slip = carrier_velocity - state(4:6)
    slip_norm = SQRT(SUM(slip*slip))
    reynolds = MAX(1.0d-12,density*diameter*slip_norm/dynamic_viscosity)
    eps = MIN(1.0d0,MAX(1.0d-12,epsilon_f))

    SELECT CASE (TRIM(el_drag_model))
    CASE ('difelice')
      cd = EL_SPHERE_CD(reynolds)
      chi = 3.7d0 - 0.65d0*EXP(-0.5d0*(1.5d0-LOG10(reynolds))**2)
      area = ACOS(-1.0d0)*diameter*diameter/4.0d0
      IF (slip_norm.GT.0.0d0) THEN
        drag_B = 0.5d0*density*cd*area*slip_norm*eps**(2.0d0-chi)
      ELSE
        drag_B = 3.0d0*ACOS(-1.0d0)*dynamic_viscosity*diameter* &
                 eps**(2.0d0-chi)
      END IF
    CASE DEFAULT
      correction = 1.0d0
      IF (TRIM(el_drag_model).EQ.'schiller_naumann') THEN
        IF (reynolds.LT.1000.0d0) THEN
          correction = 1.0d0 + 0.15d0*reynolds**0.687d0
        ELSE
          correction = 0.44d0*reynolds/24.0d0
        END IF
      END IF
      drag_B = 3.0d0*ACOS(-1.0d0)*dynamic_viscosity*diameter*correction
    END SELECT

    force = drag_B*slip

  END SUBROUTINE EL_DRAG_CLOSURE

  SUBROUTINE EL_LIFT_CLOSURE(state, carrier_velocity, grad_u, density, &
                             viscosity, epsilon_f, lift)

    REAL*8, INTENT(IN) :: state(8), carrier_velocity(3), grad_u(3,3)
    REAL*8, INTENT(IN) :: density, viscosity, epsilon_f
    REAL*8, INTENT(OUT) :: lift(3)
    REAL*8 :: slip(3), omega(3), cross_term(3), omega_norm, coeff
    REAL*8 :: dynamic_viscosity, wall_factor

    lift = 0.0d0
    IF (TRIM(el_lift_model).EQ.'none') RETURN

    dynamic_viscosity = density*viscosity
    IF (dynamic_viscosity.LE.0.0d0 .OR. state(7).LE.0.0d0) RETURN

    slip = carrier_velocity - state(4:6)
    omega(1) = grad_u(3,2) - grad_u(2,3)
    omega(2) = grad_u(1,3) - grad_u(3,1)
    omega(3) = grad_u(2,1) - grad_u(1,2)
    omega_norm = SQRT(SUM(omega*omega))
    IF (omega_norm.LE.0.0d0) RETURN

    cross_term = EL_CROSS(slip,omega)
    coeff = 1.615d0*(2.0d0*state(7))**2* &
            SQRT(density*dynamic_viscosity*omega_norm)
    wall_factor = 1.0d0
    IF (TRIM(el_lift_model).EQ.'saffman_mei_wall') THEN
      wall_factor = MIN(1.0d0,MAX(0.0d0,epsilon_f))
    END IF
    lift = wall_factor*coeff*cross_term/MAX(omega_norm,TINY(1.0d0))

  END SUBROUTINE EL_LIFT_CLOSURE

  REAL*8 FUNCTION EL_PARTICLE_VOLUME(radius)

    REAL*8, INTENT(IN) :: radius

    IF (radius.GT.0.0d0) THEN
      EL_PARTICLE_VOLUME = 4.0d0*ACOS(-1.0d0)*radius**3/3.0d0
    ELSE
      EL_PARTICLE_VOLUME = 0.0d0
    END IF

  END FUNCTION EL_PARTICLE_VOLUME

  REAL*8 FUNCTION EL_SPHERE_CD(reynolds)

    REAL*8, INTENT(IN) :: reynolds

    IF (reynolds.LT.1000.0d0) THEN
      EL_SPHERE_CD = 24.0d0/reynolds*(1.0d0 + 0.15d0*reynolds**0.687d0)
    ELSE
      EL_SPHERE_CD = 0.44d0
    END IF

  END FUNCTION EL_SPHERE_CD

  PURE FUNCTION EL_CROSS(a,b) RESULT(c)

    REAL*8, INTENT(IN) :: a(3), b(3)
    REAL*8 :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  END FUNCTION EL_CROSS

END MODULE EL_FORCES
