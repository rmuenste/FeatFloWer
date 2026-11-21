MODULE EL_FORCES

  USE MPI, ONLY: MPI_COMM_WORLD
  USE EL_CONFIG, ONLY: el_fluid_visc_dynamic, el_drag_model, el_eps_f_min, el_pressure_force, &
                       el_drag_closure_id, EL_DRAG_DIFELICE, EL_DRAG_DIFELICE_STD, &
                       EL_DRAG_STOKES, EL_DRAG_SCHILLER_NAUMANN, EL_DRAG_GIDASPOW, &
                       EL_DRAG_WENYU, EL_DRAG_ERGUN, &
                       el_pressure_feedback, &
                       el_lift_model, el_inertial_lift, &
                       el_inertial_lift_umax, el_prescribed_umax, &
                       el_cylinder_center, el_cylinder_radius, &
                       el_cylinder_axis
  USE EL_GEOMETRY, ONLY: el_domain_set, EL_WALL_DISTANCE_AND_NORMAL
  USE EL_MOMENTUM_MODEL, ONLY: EL_MM_DRAG_SCALE, EL_MM_VISCOUS_FORCE, EL_MM_VISCOUS_FEEDBACK
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
    ! viscous force V c <lap u> of the momentum models Bfull / A (zero otherwise)
    REAL*8 :: viscous(3) = 0.0d0
    ! everything the particle receives besides the drag (pressure + lift + buoyant weight
    ! + viscous) and everything besides the drag whose reaction the fluid receives
    REAL*8 :: particle_nondrag(3) = 0.0d0
    REAL*8 :: feedback_nondrag(3) = 0.0d0
    REAL*8 :: grav_buoy(3) = 0.0d0
    REAL*8 :: particle_total(3) = 0.0d0
    REAL*8 :: feedback_force(3) = 0.0d0
    REAL*8 :: drag_B = 0.0d0
    REAL*8 :: u_f(3) = 0.0d0
  END TYPE tElForceResult

CONTAINS

  ! visc_laplacian (optional): kernel average <lap u> at the particle (EL_VISCOUS), supplied by
  ! the transfer pass when the momentum model gives the particle a viscous force.
  SUBROUTINE EL_COMPUTE_PARTICLE_FORCES(state, sample, density, viscosity, &
      gravity, result, visc_laplacian)

    REAL*8, INTENT(IN) :: state(8), sample(EL_SAMPLE_SIZE)
    REAL*8, INTENT(IN) :: density, viscosity, gravity(3)
    TYPE(tElForceResult), INTENT(OUT) :: result
    REAL*8, INTENT(IN), OPTIONAL :: visc_laplacian(3)

    REAL*8 :: norm, grad_u(3,3), grad_p(3), eps_f, volume
    REAL*8 :: inertial_lift(3), drag_scale

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
    ! Momentum model B carries the whole interaction in the drag: force and coefficient are
    ! divided by eps once, here, before every particle and feedback consumer (factor 1 and
    ! therefore bitwise unchanged for every other model).
    drag_scale = EL_MM_DRAG_SCALE(eps_f)
    result%drag = drag_scale*result%drag
    result%drag_B = drag_scale*result%drag_B
    IF (el_pressure_force) result%pressure = -volume*grad_p
    CALL EL_LIFT_CLOSURE(state, result%u_f, grad_u, density, viscosity, &
                         eps_f, result%lift)
    IF (TRIM(el_inertial_lift).NE.'none') THEN
      CALL EL_INERTIAL_LIFT_FORCE(state, density, inertial_lift)
      result%lift = result%lift + inertial_lift
    END IF
    result%grav_buoy = (state(8)-density)*volume*gravity

    result%particle_total = result%drag + result%pressure + result%lift + &
                            result%grav_buoy
    result%feedback_force = result%drag + result%lift
    ! Model-A reaction: the fluid loses what the particle gains from grad p
    ! (cell-averaged -eps grad p). Keyed so certified logs stay bitwise.
    IF (el_pressure_feedback) result%feedback_force = result%feedback_force + &
                                                      result%pressure
    ! Non-drag parts in the summation order the transfer pass has always used (bitwise for
    ! the campaign set): what the particle receives, and what is fed back to the fluid.
    result%particle_nondrag = result%pressure + result%lift + result%grav_buoy
    result%feedback_nondrag = result%lift
    IF (el_pressure_feedback) result%feedback_nondrag = result%feedback_nondrag + &
                                                        result%pressure
    ! Viscous force of the momentum models Bfull / A, with the coefficient of the fluid's own
    ! diffusion operator; fed back under Bfull.
    IF (PRESENT(visc_laplacian) .AND. EL_MM_VISCOUS_FORCE()) THEN
      result%viscous = volume*EL_FLUID_OPERATOR_VISCOSITY(density, viscosity)*visc_laplacian
      result%particle_total = result%particle_total + result%viscous
      result%particle_nondrag = result%particle_nondrag + result%viscous
      IF (EL_MM_VISCOUS_FEEDBACK()) THEN
        result%feedback_force = result%feedback_force + result%viscous
        result%feedback_nondrag = result%feedback_nondrag + result%viscous
      END IF
    END IF

  END SUBROUTINE EL_COMPUTE_PARTICLE_FORCES

  ! Coefficient c of the viscous term c lap u in the fluid momentum equation as the solver
  ! assembles it on the Newtonian path (DIFFQ2_NEWT: Properties%Viscosity(1), no density
  ! factor, against a mass matrix that carries the density). The particle viscous force and
  ! model A use the same coefficient so that particle and fluid see one stress. The argument
  ! is the kinematic viscosity of the closures: under ELFluidViscosity = dynamic (default) the
  ! operator carries mu = rho nu, under 'kinematic' (campaign behaviour) that value itself,
  ! i.e. mu/rho (EL_CONFIG). This function is the single place the viscous force follows it.
  PURE REAL*8 FUNCTION EL_FLUID_OPERATOR_VISCOSITY(density, viscosity)
    REAL*8, INTENT(IN) :: density, viscosity
    IF (el_fluid_visc_dynamic) THEN
      EL_FLUID_OPERATOR_VISCOSITY = density*viscosity
    ELSE
      EL_FLUID_OPERATOR_VISCOSITY = viscosity
    END IF
  END FUNCTION EL_FLUID_OPERATOR_VISCOSITY

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

    REAL*8 :: diameter, slip(3), slip_norm, reynolds, reynolds_raw
    REAL*8 :: dynamic_viscosity, eps
    INTEGER :: closure, ierr

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

    ! Closure table: EL_CONFIG el_drag_names. Every closure returns the per-particle
    ! coefficient drag_B of the model-A form (F_d = drag_B (u_f - v), see the header of
    ! EL_DRAG_B_GIDASPOW); model B divides it by eps afterwards (EL_COMPUTE_PARTICLE_FORCES).
    ! The id is resolved once from the deck name (EL_CONFIG:EL_RESOLVE_DRAG_CLOSURE); every
    ! registered closure has its CASE, an id without one stops the run (a closure added to
    ! the table but not here can never fall back to another law).
    closure = el_drag_closure_id
    SELECT CASE (closure)
    CASE (EL_DRAG_DIFELICE, EL_DRAG_DIFELICE_STD)
      drag_B = EL_DRAG_B_DIFELICE(density, dynamic_viscosity, diameter, slip_norm, eps, &
                                  closure.EQ.EL_DRAG_DIFELICE_STD)
    CASE (EL_DRAG_GIDASPOW)
      drag_B = EL_DRAG_B_GIDASPOW(density, dynamic_viscosity, diameter, slip_norm, eps)
    CASE (EL_DRAG_WENYU)
      drag_B = EL_DRAG_B_WENYU(density, dynamic_viscosity, diameter, slip_norm, eps)
    CASE (EL_DRAG_ERGUN)
      drag_B = EL_DRAG_B_ERGUN(density, dynamic_viscosity, diameter, slip_norm, eps)
    CASE (EL_DRAG_STOKES, EL_DRAG_SCHILLER_NAUMANN)
      drag_B = EL_DRAG_B_SINGLE_SPHERE(dynamic_viscosity, diameter, reynolds, &
                                       closure.EQ.EL_DRAG_SCHILLER_NAUMANN)
    CASE DEFAULT
      WRITE(*,'(A,I0,A,A,A)') 'Fatal: EL_DRAG_CLOSURE: closure id ', closure, ' (ELDragModel = ', &
        TRIM(el_drag_model), ') has no dispatch case (unresolved or not implemented).'
      CALL MPI_Abort(MPI_COMM_WORLD, 915, ierr)
    END SELECT

    force = drag_B*slip

  END SUBROUTINE EL_DRAG_CLOSURE

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
  ! (Arithmetic moved verbatim from the former SELECT CASE of EL_DRAG_CLOSURE.)
  ! ----------------------------------------------------------------------
  PURE REAL*8 FUNCTION EL_DRAG_B_DIFELICE(density, dynamic_viscosity, diameter, slip_norm, &
                                          eps, standard) RESULT(drag_B)
    REAL*8, INTENT(IN) :: density, dynamic_viscosity, diameter, slip_norm, eps
    LOGICAL, INTENT(IN) :: standard
    REAL*8 :: reynolds_raw, reynolds, cd, chi, area

    reynolds_raw = density*eps*diameter*slip_norm/dynamic_viscosity
    reynolds = MAX(1.0d-12,reynolds_raw)
    cd = (0.63d0 + 4.8d0/SQRT(reynolds))**2
    chi = 3.7d0 - 0.65d0*EXP(-0.5d0*(1.5d0-LOG10(reynolds))**2)
    area = ACOS(-1.0d0)*diameter*diameter/4.0d0
    IF (reynolds_raw.GT.1.0d-12) THEN
      drag_B = 0.5d0*density*cd*area*slip_norm*eps**(-chi)
      ! 'difelice_std': the Kafui (2002) / Zhu et al. (2007) / CFDEM form
      ! carries eps^2 on the squared slip, i.e. eps^(2-chi) overall;
      ! 'difelice' (campaign form) is eps^(-chi). The two differ by eps^2.
      IF (standard) drag_B = drag_B*eps*eps
    ELSE
      drag_B = 2.88d0*ACOS(-1.0d0)*dynamic_viscosity*diameter* &
               eps**(-chi-1.0d0)
      IF (standard) drag_B = drag_B*eps*eps
    END IF
  END FUNCTION EL_DRAG_B_DIFELICE

  ! Single sphere without porosity: Stokes, or Schiller-Naumann with the single-particle
  ! Reynolds number (arithmetic moved verbatim from the former CASE DEFAULT).
  PURE REAL*8 FUNCTION EL_DRAG_B_SINGLE_SPHERE(dynamic_viscosity, diameter, reynolds, &
                                               schiller_naumann) RESULT(drag_B)
    REAL*8, INTENT(IN) :: dynamic_viscosity, diameter, reynolds
    LOGICAL, INTENT(IN) :: schiller_naumann
    REAL*8 :: correction

    correction = 1.0d0
    IF (schiller_naumann) THEN
      IF (reynolds.LT.1000.0d0) THEN
        correction = 1.0d0 + 0.15d0*reynolds**0.687d0
      ELSE
        correction = 0.44d0*reynolds/24.0d0
      END IF
    END IF
    drag_B = 3.0d0*ACOS(-1.0d0)*dynamic_viscosity*diameter*correction
  END FUNCTION EL_DRAG_B_SINGLE_SPHERE

  ! ----------------------------------------------------------------------
  ! Gidaspow (1994) drag: Ergun for eps <= 0.8, Wen-Yu above.
  !
  ! Convention. beta is the interphase momentum-exchange coefficient per unit mixture volume
  ! of the two-fluid model A (Gidaspow): the particle phase receives beta (u - v) per unit
  ! volume besides its share -(1 - eps) grad p of the pressure gradient, u - v the
  ! INTERSTITIAL slip. With n = (1 - eps)/V_p particles per unit volume the per-particle
  ! coefficient is
  !     drag_B = V_p beta / (1 - eps),     F_d = drag_B (u_f - v),
  ! the model-A form of every closure in this module (difelice_std is the same form). The
  ! momentum models consume it as follows (EL_CONFIG model table, EL_MM_DRAG_SCALE):
  !     A, Bfull : F_d as above; the particle also receives -V_p grad p (pressure force)
  !     B        : F_d / eps, no pressure force (the whole interaction is in the drag)
  !     none     : F_d as above, pressure force per the campaign keys
  ! In a fixed bed at superficial velocity U the model-A drag per unit volume is
  ! n F_d = beta U/eps and the pressure gradient carries the rest: n (F_d + V_p |grad p|) =
  ! n F_d / eps = beta U/eps^2 = the Ergun gradient (unit test test_el_drag_closures).
  !
  ! Verified against the sources (2026-11-19, liggghts-dem tree, read-only):
  !   CFDEM GidaspowDrag.C: betaP = beta/phiP (per particle volume), dragCoefficient = Vs betaP,
  !     divided by voidfraction only for modelType "B"; dense branch for voidfraction <=
  !     switchingVoidfraction (0.8), dilute (Wen-Yu, Schiller-Naumann C_d without a Re cap)
  !     above.
  !   MFiX drag_gs.f DRAG_GIDASPOW: lDgA = beta/(1-eps), Ergun for EPg <= 0.8, Wen-Yu with
  !     C_d = 0.44 above Re 1000; des/drag_gp_des.f: F_gp = DgA PVOL (model A), / EP_g (model B).
  !     MFiX returns zero drag at exactly zero slip (RE == 0); here the slip -> 0 limit is used.
  ! Both codes therefore use exactly this per-particle form for model A and the 1/eps form for
  ! model B. The switch: Ergun at eps = 0.8 itself, as in both codes. The two branches do not
  ! join at 0.8 (Ergun/Wen-Yu = 1.153 at Re 0, 1.143 at slip 0.5 m/s for the case-11 particles).
  ! ----------------------------------------------------------------------
  PURE REAL*8 FUNCTION EL_DRAG_B_GIDASPOW(density, dynamic_viscosity, diameter, slip_norm, &
                                          eps) RESULT(drag_B)
    REAL*8, INTENT(IN) :: density, dynamic_viscosity, diameter, slip_norm, eps
    REAL*8, PARAMETER :: EPS_SWITCH = 0.8d0
    IF (eps.LE.EPS_SWITCH) THEN
      drag_B = EL_DRAG_B_ERGUN(density, dynamic_viscosity, diameter, slip_norm, eps)
    ELSE
      drag_B = EL_DRAG_B_WENYU(density, dynamic_viscosity, diameter, slip_norm, eps)
    END IF
  END FUNCTION EL_DRAG_B_GIDASPOW

  ! Ergun: beta = 150 (1-eps)^2 mu / (eps d^2) + 1.75 (1-eps) rho |u-v| / d, so
  ! drag_B = V_p [150 (1-eps) mu / (eps d^2) + 1.75 rho |u-v| / d]. Finite at zero slip.
  PURE REAL*8 FUNCTION EL_DRAG_B_ERGUN(density, dynamic_viscosity, diameter, slip_norm, &
                                       eps) RESULT(drag_B)
    REAL*8, INTENT(IN) :: density, dynamic_viscosity, diameter, slip_norm, eps
    REAL*8 :: volume
    volume = ACOS(-1.0d0)*diameter**3/6.0d0
    drag_B = volume*(150.0d0*(1.0d0-eps)*dynamic_viscosity/(eps*diameter*diameter) + &
                     1.75d0*density*slip_norm/diameter)
  END FUNCTION EL_DRAG_B_ERGUN

  ! Wen-Yu: beta = 3/4 C_d eps (1-eps) rho |u-v| / d eps^-2.65, C_d = 24/Re (1 + 0.15 Re^0.687)
  ! for Re < 1000 and 0.44 above, Re = rho eps |u-v| d / mu; so
  ! drag_B = V_p 3/4 (C_d |u-v|) rho eps eps^-2.65 / d. C_d |u-v| is formed without the
  ! division by Re (24 mu/(rho eps d) (1 + 0.15 Re^0.687)), as in CFDEM: the zero-slip limit
  ! 3 pi mu d eps^-2.65 is reached continuously, no special branch.
  PURE REAL*8 FUNCTION EL_DRAG_B_WENYU(density, dynamic_viscosity, diameter, slip_norm, &
                                       eps) RESULT(drag_B)
    REAL*8, INTENT(IN) :: density, dynamic_viscosity, diameter, slip_norm, eps
    REAL*8 :: volume, reynolds, cd_slip
    volume = ACOS(-1.0d0)*diameter**3/6.0d0
    reynolds = density*eps*slip_norm*diameter/dynamic_viscosity
    IF (reynolds.LT.1000.0d0) THEN
      cd_slip = 24.0d0*dynamic_viscosity/(density*eps*diameter)* &
                (1.0d0 + 0.15d0*reynolds**0.687d0)
    ELSE
      cd_slip = 0.44d0*slip_norm
    END IF
    drag_B = volume*0.75d0*cd_slip*density*eps*eps**(-2.65d0)/diameter
  END FUNCTION EL_DRAG_B_WENYU

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
      IF (.NOT.el_domain_set) THEN
        IF (.NOT.warned_unset_box) THEN
          WRITE(*,'(A)') 'WARNING: saffman_mei_wall requested without EL domain; using saffman_mei.'
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
      ! no wall (every face of the box is a periodic plane of the coupling): unbounded form
      IF (wall_distance.GE.HUGE(1.0d0)) THEN
        lift = mei_factor*saffman_lift
        RETURN
      END IF
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

  SUBROUTINE EL_INERTIAL_LIFT_FORCE(state, density, force)

    ! Neutrally-buoyant inertial (Segre-Silberberg) lift for tube Poiseuille
    ! flow. Matched-asymptotics lift profile digitized from Matas, Morris &
    ! Guazzelli (2004) J. Fluid Mech. 515, 171-195, figure 14, Rc = 30 curve
    ! (Rc = Umax*D/nu, Asmolov's channel Reynolds number; Schonberg-Hinch/
    ! Asmolov theory, valid for Re_p = Rc*(a/D)^2 << 1). With the paper's
    ! normalization Fl*Rc^{1/2}*eps^{-3}/(eta*Umax*a), eps = sqrt(Re_p/2),
    ! the dimensional force collapses to
    !   F_r(s) = fhat(s) * rho_f * Umax^2 * a^4 / (8*sqrt(2)*R^2),  s = r/R,
    ! applied along the radial unit vector (positive fhat pushes outward).
    ! Channel profile applied radially in the tube (the paper's own channel
    ! -> pipe mapping, 2z/l <-> r/R); zero crossing at s_eq ~ 0.68.
    REAL*8, INTENT(IN) :: state(8)
    REAL*8, INTENT(IN) :: density
    REAL*8, INTENT(OUT) :: force(3)

    INTEGER, PARAMETER :: NTAB = 19
    ! s = 0.00 .. 0.90 in steps of 0.05
    REAL*8, PARAMETER :: FHAT(NTAB) = (/ 0.00d0, 1.10d0, 1.70d0, 2.30d0, &
      2.90d0, 3.50d0, 3.85d0, 4.15d0, 4.20d0, 3.90d0, 3.50d0, 2.85d0, &
      2.00d0, 0.80d0, -0.80d0, -2.60d0, -5.20d0, -8.20d0, -12.10d0 /)
    REAL*8, PARAMETER :: DS = 0.05d0

    REAL*8 :: axial(3), radial(3), dx(3), r, s, umax, fhat_s, coeff, frac
    INTEGER :: ibin

    force = 0.0d0
    IF (el_cylinder_radius.LE.0.0d0) RETURN

    umax = el_inertial_lift_umax
    IF (umax.LE.0.0d0) umax = el_prescribed_umax
    IF (umax.LE.0.0d0) RETURN

    axial = 0.0d0
    SELECT CASE (TRIM(el_cylinder_axis))
    CASE ('x')
      axial(1) = 1.0d0
    CASE ('y')
      axial(2) = 1.0d0
    CASE DEFAULT
      axial(3) = 1.0d0
    END SELECT

    dx = state(1:3) - el_cylinder_center
    radial = dx - SUM(dx*axial)*axial
    r = SQRT(SUM(radial*radial))
    IF (r.LE.1.0d-14*el_cylinder_radius) RETURN
    radial = radial/r
    s = r/el_cylinder_radius

    ! Linear interpolation on the table; linear extrapolation of the last
    ! segment beyond s = 0.90 (steep near-wall repulsion; the particle
    ! center is geometrically capped at s = 1 - a/R).
    IF (s.GE.DBLE(NTAB-1)*DS) THEN
      fhat_s = FHAT(NTAB) + (FHAT(NTAB)-FHAT(NTAB-1))/DS * &
               (s - DBLE(NTAB-1)*DS)
    ELSE
      ibin = INT(s/DS) + 1
      frac = s/DS - DBLE(ibin-1)
      fhat_s = (1.0d0-frac)*FHAT(ibin) + frac*FHAT(ibin+1)
    END IF

    coeff = density*umax*umax*state(7)**4 / &
            (8.0d0*SQRT(2.0d0)*el_cylinder_radius**2)
    force = fhat_s*coeff*radial

  END SUBROUTINE EL_INERTIAL_LIFT_FORCE

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
