MODULE EL_FORCES

  USE EL_CONFIG, ONLY: el_drag_model, el_eps_f_min, el_pressure_force, &
                       el_lift_model
  USE EL_GEOMETRY, ONLY: el_domain_box_set, EL_WALL_DISTANCE_AND_NORMAL
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

    REAL*8 :: diameter, slip(3), slip_norm, reynolds, reynolds_raw, correction
    REAL*8 :: dynamic_viscosity, eps, cd, area, chi

    force = 0.0d0
    drag_B = 0.0d0
    diameter = 2.0d0*state(7)
    dynamic_viscosity = density*viscosity
    IF (diameter.LE.0.0d0 .OR. dynamic_viscosity.LE.0.0d0) RETURN

    slip = carrier_velocity - state(4:6)
    slip_norm = SQRT(SUM(slip*slip))
    eps = MIN(1.0d0,MAX(1.0d-12,epsilon_f))
    reynolds_raw = density*diameter*slip_norm/dynamic_viscosity
    reynolds = MAX(1.0d-12,reynolds_raw)

    SELECT CASE (TRIM(el_drag_model))
    CASE ('difelice')
      ! ----------------------------------------------------------------------
      ! Di Felice drag is treated as one consistent closure package:
      ! interstitial slip from the EL sampler, Re_p = rho_f eps_f d |u_rel|/mu,
      ! Dallavalle-type Cd = (0.63 + 4.8/sqrt(Re_p))**2, and voidage factor
      ! eps_f**(-chi). The older note that suggested dropping to 1-chi applies
      ! to the Zhou set-A / Schiller-Naumann convention, not to the original
      ! Di Felice form used here. Richardson-Zaki validation remains the
      ! empirical arbiter for dense suspensions.
      !
      ! Since Cd(Re->0) = 23.04/Re, this closure's dilute low-Re drag is
      ! 2.88*pi*mu*d rather than the Stokes 3*pi*mu*d, about 4% lower. The
      ! zero-slip drag_B fallback is therefore the slip->0 limit
      ! 2.88*pi*mu*d*eps_f**(-chi-1).
      ! ----------------------------------------------------------------------
      reynolds_raw = density*eps*diameter*slip_norm/dynamic_viscosity
      reynolds = MAX(1.0d-12,reynolds_raw)
      cd = (0.63d0 + 4.8d0/SQRT(reynolds))**2
      chi = 3.7d0 - 0.65d0*EXP(-0.5d0*(1.5d0-LOG10(reynolds))**2)
      area = ACOS(-1.0d0)*diameter*diameter/4.0d0
      IF (reynolds_raw.GT.1.0d-12) THEN
        drag_B = 0.5d0*density*cd*area*slip_norm*eps**(-chi)
      ELSE
        drag_B = 2.88d0*ACOS(-1.0d0)*dynamic_viscosity*diameter* &
                 eps**(-chi-1.0d0)
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
    REAL*8 :: dynamic_viscosity, diameter, slip_norm, reynolds, re_shear
    REAL*8 :: beta, mei_factor, saffman_lift(3), zeng_lift(3)
    REAL*8 :: wall_distance, wall_normal(3), l_over_d, cl_zeng, area, blend
    REAL*8 :: zeng_re
    LOGICAL, SAVE :: warned_unset_box = .FALSE.
    LOGICAL, SAVE :: warned_zeng_clamp = .FALSE.

    lift = 0.0d0
    IF (TRIM(el_lift_model).EQ.'none') RETURN

    dynamic_viscosity = density*viscosity
    IF (dynamic_viscosity.LE.0.0d0 .OR. state(7).LE.0.0d0) RETURN

    slip = carrier_velocity - state(4:6)
    slip_norm = SQRT(SUM(slip*slip))
    IF (slip_norm.LE.0.0d0) RETURN

    diameter = 2.0d0*state(7)
    omega(1) = grad_u(3,2) - grad_u(2,3)
    omega(2) = grad_u(1,3) - grad_u(3,1)
    omega(3) = grad_u(2,1) - grad_u(1,2)
    omega_norm = SQRT(SUM(omega*omega))
    IF (omega_norm.LE.0.0d0) RETURN

    cross_term = EL_CROSS(slip,omega)
    coeff = 1.615d0*diameter**2*SQRT(density*dynamic_viscosity*omega_norm)
    saffman_lift = coeff*cross_term/MAX(omega_norm,TINY(1.0d0))

    reynolds = density*diameter*slip_norm/dynamic_viscosity
    re_shear = density*diameter**2*omega_norm/dynamic_viscosity
    beta = 0.5d0*re_shear/MAX(reynolds,TINY(1.0d0))
    beta = MIN(0.4d0,MAX(0.005d0,beta))
    mei_factor = EL_MEI_LIFT_FACTOR(reynolds,beta)

    SELECT CASE (TRIM(el_lift_model))
    CASE ('saffman')
      lift = saffman_lift
    CASE ('saffman_mei')
      lift = mei_factor*saffman_lift
    CASE ('saffman_mei_wall')
      IF (.NOT.el_domain_box_set) THEN
        IF (.NOT.warned_unset_box) THEN
          WRITE(*,'(A)') 'WARNING: saffman_mei_wall requested without EL_SET_DOMAIN_BOX; using saffman_mei.'
          warned_unset_box = .TRUE.
        END IF
        lift = mei_factor*saffman_lift
        RETURN
      END IF

      IF (reynolds.LE.0.5d0) THEN
        lift = mei_factor*saffman_lift
        RETURN
      END IF

      CALL EL_WALL_DISTANCE_AND_NORMAL(state(1:3), wall_distance, wall_normal)
      l_over_d = MAX(wall_distance/diameter,0.5d0)
      zeng_re = reynolds
      IF (zeng_re.GT.200.0d0) THEN
        IF (.NOT.warned_zeng_clamp) THEN
          WRITE(*,'(A)') 'WARNING: Zeng wall-lift Re exceeds 200; clamping correlation input.'
          warned_zeng_clamp = .TRUE.
        END IF
        zeng_re = 200.0d0
      END IF

      cl_zeng = EL_ZENG_SHEAR_LIFT_COEFFICIENT(MAX(zeng_re,1.0d-12), l_over_d)
      area = ACOS(-1.0d0)*diameter*diameter/4.0d0
      zeng_lift = cl_zeng*0.5d0*density*slip_norm*slip_norm*area*wall_normal

      IF (reynolds.GE.2.0d0) THEN
        lift = zeng_lift
      ELSE
        blend = LOG(reynolds/0.5d0)/LOG(4.0d0)
        blend = MIN(1.0d0,MAX(0.0d0,blend))
        lift = (1.0d0-blend)*mei_factor*saffman_lift + blend*zeng_lift
      END IF
    END SELECT

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

  REAL*8 FUNCTION EL_MEI_LIFT_FACTOR(reynolds, beta)

    REAL*8, INTENT(IN) :: reynolds, beta
    REAL*8 :: re_use, beta_use

    re_use = MAX(0.0d0,reynolds)
    beta_use = MAX(0.0d0,beta)
    IF (re_use.LE.40.0d0) THEN
      EL_MEI_LIFT_FACTOR = (1.0d0 - 0.3314d0*SQRT(beta_use))* &
                           EXP(-re_use/10.0d0) + 0.3314d0*SQRT(beta_use)
    ELSE
      EL_MEI_LIFT_FACTOR = 0.0524d0*SQRT(beta_use*re_use)
    END IF

  END FUNCTION EL_MEI_LIFT_FACTOR

  REAL*8 FUNCTION EL_ZENG_CONTACT_LIFT_COEFFICIENT(reynolds)

    REAL*8, INTENT(IN) :: reynolds
    REAL*8 :: re_use

    re_use = MAX(0.0d0,reynolds)
    EL_ZENG_CONTACT_LIFT_COEFFICIENT = &
      3.663d0/(re_use*re_use + 0.1173d0)**0.22d0

  END FUNCTION EL_ZENG_CONTACT_LIFT_COEFFICIENT

  REAL*8 FUNCTION EL_ZENG_SHEAR_LIFT_COEFFICIENT(reynolds, l_over_d)

    REAL*8, INTENT(IN) :: reynolds, l_over_d
    REAL*8 :: re_use, delta, alpha_sl, beta_sl, lambda_sl

    ! Zeng, Najjar, Balachandar & Fischer, Phys. Fluids 21, 033302 (2009):
    ! Eq. (19) gives the contact limit CLs,w; Eqs. (28)-(29) give the
    ! stationary-particle wall-bounded shear composite. Here L is center-wall
    ! distance normalized by particle diameter, so contact is L=0.5 and
    ! delta = L - 0.5. The coefficient is normalized by
    ! F_L = CLs * 0.5*rho*|u_slip|**2*pi*d**2/4 and is positive away from
    ! the nearest wall.
    re_use = MAX(0.0d0,reynolds)
    delta = MAX(0.0d0,l_over_d - 0.5d0)
    alpha_sl = -EXP(-0.3d0 + 0.025d0*re_use)
    beta_sl = 0.8d0 + 0.01d0*re_use
    lambda_sl = (1.0d0 - EXP(-delta))*(re_use/250.0d0)**2.5d0
    EL_ZENG_SHEAR_LIFT_COEFFICIENT = &
      EL_ZENG_CONTACT_LIFT_COEFFICIENT(re_use)* &
      EXP(-0.5d0*delta*(re_use/250.0d0)**(4.0d0/3.0d0))* &
      (EXP(alpha_sl*delta**beta_sl) - lambda_sl)

  END FUNCTION EL_ZENG_SHEAR_LIFT_COEFFICIENT

  PURE FUNCTION EL_CROSS(a,b) RESULT(c)

    REAL*8, INTENT(IN) :: a(3), b(3)
    REAL*8 :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  END FUNCTION EL_CROSS

  !----------------------------------------------------------------------------
  ! EL_SEMI_IMPLICIT_DRAG_IMPULSE
  !
  ! Reproduce the drag impulse PE actually applies to a hydrodynamic particle
  ! over one CFD/PE coupling step, so the fluid can be fed exactly what the
  ! particle's drag removed (issue-D fix). PE advances the body with n_sub
  ! semi-implicit sub-steps of dt/n_sub on the SAME frozen hydro state
  ! (u_f, drag_B, f_other); this mirrors pe::elSemiImplicitVelocity
  ! (libs/pe/pe/core/collisionsystem/ELSemiImplicitDrag.h) exactly.
  !
  !   dt_sub = dt / n_sub
  !   U      = U_old; repeat n_sub: U = (U + dt_sub/m*(B*u_f + f_other))/(1+dt_sub*B/m)
  !   I_drag = m*(U - U_old) - f_other*dt        (exact: f_other contributes f_other*dt)
  !   F_drag_effective = I_drag / dt
  !
  ! f_other (= pressure + lift + grav_buoy) is subtracted out, so the result is
  ! the pure DRAG impulse; lift is re-added explicitly by the caller. Returns 0
  ! for degenerate inputs (no drag / zero step / massless).
  !----------------------------------------------------------------------------
  SUBROUTINE EL_SEMI_IMPLICIT_DRAG_IMPULSE(u_old, m_p, drag_B, u_f, f_other, dt, &
                                           n_sub, drag_force_eff)

    REAL*8, INTENT(IN) :: u_old(3), m_p, drag_B, u_f(3), f_other(3), dt
    INTEGER, INTENT(IN) :: n_sub
    REAL*8, INTENT(OUT) :: drag_force_eff(3)

    REAL*8 :: u(3), dt_sub, denom
    INTEGER :: k, nsub_use

    drag_force_eff = 0.0d0
    IF (dt.LE.0.0d0 .OR. m_p.LE.0.0d0) RETURN

    nsub_use = MAX(1, n_sub)
    dt_sub = dt / DBLE(nsub_use)
    denom = 1.0d0 + dt_sub*drag_B/m_p

    u = u_old
    DO k = 1, nsub_use
      u = (u + (dt_sub/m_p)*(drag_B*u_f + f_other)) / denom
    END DO

    drag_force_eff = (m_p*(u - u_old) - f_other*dt) / dt

  END SUBROUTINE EL_SEMI_IMPLICIT_DRAG_IMPULSE

END MODULE EL_FORCES
