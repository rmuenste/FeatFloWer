PROGRAM TEST_EL_KERNEL_FORCES

  USE EL_CONFIG, ONLY: el_kernel, el_drag_model, el_pressure_force, &
                       el_lift_model, el_drag_coupling, &
                       el_apply_fluid_feedback, el_drag_semi_implicit, &
                       el_inertial_lift, el_inertial_lift_umax, &
                       el_cylinder_center, el_cylinder_radius, &
                       el_cylinder_axis, &
                       EL_VALIDATE_CONFIG
  USE EL_KERNEL_FUNCTIONS, ONLY: EL_KERNEL_VALUE
  USE EL_FORCES, ONLY: tElForceResult, EL_DRAG_FORCE, &
                       EL_COMPUTE_PARTICLE_FORCES, EL_DRAG_CLOSURE, &
                       EL_MEI_LIFT_FACTOR, &
                       EL_ZENG_CONTACT_LIFT_COEFFICIENT, &
                       EL_ZENG_SHEAR_LIFT_COEFFICIENT, &
                       EL_INERTIAL_LIFT_FORCE
  USE EL_GEOMETRY, ONLY: EL_SET_DOMAIN_BOX
  USE EL_QUADRATURE, ONLY: EL_SAMPLE_SIZE, EL_SAMPLE_NORMALIZATION, &
                           EL_SAMPLE_UF_BEGIN, EL_SAMPLE_UF_END, &
                           EL_SAMPLE_GRAD_U_BEGIN, EL_SAMPLE_GRAD_U_END, &
                           EL_SAMPLE_PRESSURE, EL_SAMPLE_GRAD_P_BEGIN, &
                           EL_SAMPLE_GRAD_P_END, EL_SAMPLE_EPSILON_F
  USE EL_FIELDS, ONLY: el_field_data, EL_INITIALIZE, EL_WRITE_RESTART, &
                       EL_READ_RESTART, EL_FINALIZE

  IMPLICIT NONE

  INTEGER :: i, j, restart_unit
  REAL*8 :: distance, width
  REAL*8 :: state(8), carrier_velocity(3), force(3), expected(3)
  REAL*8 :: sample(EL_SAMPLE_SIZE), gravity(3), slip(3), force_eps1(3)
  TYPE(tElForceResult) :: force_result
  REAL*8 :: relative_error, diameter, dynamic_viscosity, reynolds_values(5)
  REAL*8 :: reynolds, slip_norm, correction, cd, area, chi, eps, drag_B
  REAL*8 :: expected_ratio, stokes_norm, difelice_norm
  REAL*8 :: grad_u(3,3), omega_norm, lift_coeff, mei_factor
  REAL*8 :: lift_low(3), lift_high(3), zeng_value, dns_value

  el_drag_coupling = 'explicit'
  el_apply_fluid_feedback = .FALSE.
  CALL EL_VALIDATE_CONFIG()
  IF (el_drag_semi_implicit) STOP 70
  el_drag_coupling = 'semi_implicit'
  el_apply_fluid_feedback = .TRUE.
  CALL EL_VALIDATE_CONFIG()
  IF (.NOT.el_drag_semi_implicit) STOP 71
  el_drag_coupling = 'explicit'
  el_apply_fluid_feedback = .FALSE.
  CALL EL_VALIDATE_CONFIG()

  width = 0.25d0
  DO i=0,100
    distance = width*DBLE(i)/100.0d0
    IF (EL_KERNEL_VALUE(distance,width).LT.0.0d0) STOP 1
  END DO
  IF (EL_KERNEL_VALUE(width,width).NE.0.0d0) STOP 2

  el_kernel = 'gaussian'
  DO i=0,100
    distance = width*DBLE(i)/100.0d0
    IF (EL_KERNEL_VALUE(distance,width).LT.0.0d0) STOP 3
  END DO

  state = 0.0d0
  state(7) = 0.005d0
  carrier_velocity = (/1.0d-5, -2.0d-5, 3.0d-5/)
  el_drag_model = 'stokes'
  CALL EL_DRAG_FORCE(state, carrier_velocity, 1.0d0, 1.0d-3, force)

  expected = 3.0d0*ACOS(-1.0d0)*1.0d-3*0.01d0*carrier_velocity
  relative_error = SQRT(SUM((force-expected)**2))/SQRT(SUM(expected**2))
  IF (relative_error.GT.1.0d-12) STOP 4

  el_drag_model = 'difelice'
  carrier_velocity = (/1.0d-9, -2.0d-9, 3.0d-9/)
  CALL EL_DRAG_FORCE(state, carrier_velocity, 1.0d0, 1.0d-3, force)
  diameter = 2.0d0*state(7)
  dynamic_viscosity = 1.0d0*1.0d-3
  area = ACOS(-1.0d0)*diameter*diameter/4.0d0
  slip_norm = SQRT(SUM(carrier_velocity*carrier_velocity))
  reynolds = MAX(1.0d-12,diameter*slip_norm/1.0d-3)
  cd = (0.63d0 + 4.8d0/SQRT(reynolds))**2
  expected = 0.5d0*1.0d0*cd*area*slip_norm*carrier_velocity
  relative_error = SQRT(SUM((force-expected)**2))/SQRT(SUM(expected**2))
  IF (relative_error.GT.1.0d-12) STOP 40

  reynolds = 1.0d-6
  slip_norm = reynolds*dynamic_viscosity/(1.0d0*diameter)
  carrier_velocity = (/slip_norm, 0.0d0, 0.0d0/)
  CALL EL_DRAG_FORCE(state, carrier_velocity, 1.0d0, 1.0d-3, force)
  difelice_norm = SQRT(SUM(force**2))
  el_drag_model = 'stokes'
  CALL EL_DRAG_FORCE(state, carrier_velocity, 1.0d0, 1.0d-3, force)
  stokes_norm = SQRT(SUM(force**2))
  expected_ratio = (0.63d0*SQRT(reynolds) + 4.8d0)**2 / 24.0d0
  IF (ABS(difelice_norm/stokes_norm - expected_ratio).GT.1.0d-3) STOP 54

  diameter = 2.0d0*state(7)
  dynamic_viscosity = 1.0d0*1.0d-3
  area = ACOS(-1.0d0)*diameter*diameter/4.0d0
  reynolds_values = (/1.0d-3, 1.0d-1, 1.0d0, 1.0d1, 5.0d2/)
  DO j = 1, SIZE(reynolds_values)
    reynolds = reynolds_values(j)
    slip_norm = reynolds*dynamic_viscosity/(1.0d0*diameter)
    carrier_velocity = (/slip_norm, 0.0d0, 0.0d0/)

    el_drag_model = 'stokes'
    CALL EL_DRAG_FORCE(state, carrier_velocity, 1.0d0, 1.0d-3, force)
    expected = 3.0d0*ACOS(-1.0d0)*dynamic_viscosity*diameter* &
               carrier_velocity
    relative_error = SQRT(SUM((force-expected)**2))/MAX(SQRT(SUM(expected**2)), &
                     TINY(1.0d0))
    IF (relative_error.GT.1.0d-12) STOP 50

    el_drag_model = 'schiller_naumann'
    correction = 1.0d0 + 0.15d0*reynolds**0.687d0
    CALL EL_DRAG_FORCE(state, carrier_velocity, 1.0d0, 1.0d-3, force)
    expected = 3.0d0*ACOS(-1.0d0)*dynamic_viscosity*diameter* &
               correction*carrier_velocity
    relative_error = SQRT(SUM((force-expected)**2))/MAX(SQRT(SUM(expected**2)), &
                     TINY(1.0d0))
    IF (relative_error.GT.1.0d-12) STOP 51

    el_drag_model = 'difelice'
    cd = (0.63d0 + 4.8d0/SQRT(reynolds))**2
    CALL EL_DRAG_FORCE(state, carrier_velocity, 1.0d0, 1.0d-3, force)
    expected = 0.5d0*1.0d0*cd*area*slip_norm*carrier_velocity
    relative_error = SQRT(SUM((force-expected)**2))/MAX(SQRT(SUM(expected**2)), &
                     TINY(1.0d0))
    IF (relative_error.GT.1.0d-12) STOP 52
  END DO

  sample = 0.0d0
  sample(EL_SAMPLE_NORMALIZATION) = 2.0d0
  sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) = 2.0d0*(/0.2d0,0.1d0,-0.3d0/)
  sample(EL_SAMPLE_GRAD_P_BEGIN:EL_SAMPLE_GRAD_P_END) = 2.0d0*(/4.0d0,-5.0d0,6.0d0/)
  sample(EL_SAMPLE_EPSILON_F) = 2.0d0
  gravity = (/0.0d0,0.0d0,-9.81d0/)
  state = 0.0d0
  state(7) = 0.005d0
  state(8) = 2.5d0
  el_pressure_force = .TRUE.
  el_lift_model = 'none'
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  expected = -(4.0d0*ACOS(-1.0d0)*state(7)**3/3.0d0)*(/4.0d0,-5.0d0,6.0d0/)
  IF (MAXVAL(ABS(force_result%pressure-expected)).GT.1.0d-14) STOP 41
  expected = (state(8)-1.0d0)*(4.0d0*ACOS(-1.0d0)*state(7)**3/3.0d0)*gravity
  IF (MAXVAL(ABS(force_result%grav_buoy-expected)).GT.1.0d-14) STOP 42
  IF (MAXVAL(ABS(force_result%lift)).NE.0.0d0) STOP 43
  IF (MAXVAL(ABS(force_result%particle_total - &
      (force_result%drag+force_result%pressure+force_result%grav_buoy))).GT.1.0d-14) STOP 44
  IF (MAXVAL(ABS(force_result%feedback_force-force_result%drag)).GT.1.0d-14) STOP 45
  slip = force_result%u_f - state(4:6)
  IF (MAXVAL(ABS(force_result%drag_B*slip-force_result%drag)).GT.1.0d-14) STOP 46

  force_eps1 = force_result%drag
  sample(EL_SAMPLE_EPSILON_F) = 2.0d0*0.5d0
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  IF (SQRT(SUM(force_result%drag**2)).LE.SQRT(SUM(force_eps1**2))) STOP 47
  slip = force_result%u_f - state(4:6)
  slip_norm = SQRT(SUM(slip*slip))
  eps = 0.5d0
  reynolds = MAX(1.0d-12, 1.0d0*eps*diameter*slip_norm/dynamic_viscosity)
  cd = (0.63d0 + 4.8d0/SQRT(reynolds))**2
  chi = 3.7d0 - 0.65d0*EXP(-0.5d0*(1.5d0-LOG10(reynolds))**2)
  expected = 0.5d0*1.0d0*cd*area*slip_norm*eps**(-chi)*slip
  relative_error = SQRT(SUM((force_result%drag-expected)**2))/ &
                   MAX(SQRT(SUM(expected**2)),TINY(1.0d0))
  IF (relative_error.GT.1.0d-12) STOP 53

  state(4:6) = 0.0d0
  carrier_velocity = 0.0d0
  CALL EL_DRAG_CLOSURE(state, carrier_velocity, 1.0d0, 1.0d-3, 0.5d0, &
                       force, drag_B)
  chi = 3.7d0 - 0.65d0*EXP(-0.5d0*(1.5d0-LOG10(1.0d-12))**2)
  expected_ratio = 2.88d0*ACOS(-1.0d0)*dynamic_viscosity*diameter* &
                   0.5d0**(-chi-1.0d0)
  IF (ABS(drag_B-expected_ratio)/expected_ratio.GT.1.0d-12) STOP 55

  carrier_velocity = (/1.0d-14, 0.0d0, 0.0d0/)
  CALL EL_DRAG_CLOSURE(state, carrier_velocity, 1.0d0, 1.0d-3, 0.5d0, &
                       force, drag_B)
  expected = expected_ratio*carrier_velocity
  relative_error = SQRT(SUM((force-expected)**2))/ &
                   MAX(SQRT(SUM(expected**2)),TINY(1.0d0))
  IF (relative_error.GT.1.0d-12) STOP 56

  sample(EL_SAMPLE_EPSILON_F) = 2.0d0
  state(4:6) = 0.0d0
  sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) = 2.0d0*(/1.0d0, 0.0d0, 0.0d0/)
  grad_u = 0.0d0
  grad_u(1,2) = 2.0d0
  sample(EL_SAMPLE_GRAD_U_BEGIN:EL_SAMPLE_GRAD_U_END) = &
    2.0d0*(/grad_u(1,:), grad_u(2,:), grad_u(3,:)/)
  el_lift_model = 'saffman'
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  omega_norm = 2.0d0
  lift_coeff = 1.615d0*diameter*diameter* &
               SQRT(1.0d0*dynamic_viscosity*omega_norm)
  expected = (/0.0d0, lift_coeff, 0.0d0/)
  relative_error = SQRT(SUM((force_result%lift-expected)**2))/ &
                   MAX(SQRT(SUM(expected**2)),TINY(1.0d0))
  IF (relative_error.GT.1.0d-12) STOP 48

  IF (ABS(EL_MEI_LIFT_FACTOR(0.0d0,0.1d0)-1.0d0).GT.1.0d-12) STOP 57
  IF (ABS(EL_MEI_LIFT_FACTOR(100.0d0,0.1d0)- &
      0.1657033493928231d0).GT.1.0d-12) STOP 58
  IF (ABS(EL_MEI_LIFT_FACTOR(40.0d0,1.0d0)- &
      ((1.0d0-0.3314d0)*EXP(-4.0d0)+0.3314d0)).GT.1.0d-12) STOP 59

  el_lift_model = 'saffman_mei'
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  mei_factor = EL_MEI_LIFT_FACTOR(10.0d0,0.01d0)
  expected = mei_factor*(/0.0d0, lift_coeff, 0.0d0/)
  relative_error = SQRT(SUM((force_result%lift-expected)**2))/ &
                   MAX(SQRT(SUM(expected**2)),TINY(1.0d0))
  IF (relative_error.GT.1.0d-12) STOP 60

  sample(EL_SAMPLE_EPSILON_F) = 2.0d0*0.5d0
  el_lift_model = 'saffman_mei_wall'
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  expected = mei_factor*(/0.0d0, lift_coeff, 0.0d0/)
  relative_error = SQRT(SUM((force_result%lift-expected)**2))/ &
                   MAX(SQRT(SUM(expected**2)),TINY(1.0d0))
  IF (relative_error.GT.1.0d-12) STOP 49

  IF (ABS(EL_ZENG_CONTACT_LIFT_COEFFICIENT(0.0d0)-5.86936821706617d0) &
      .GT.1.0d-12) STOP 61
  dns_value = 2.653d0
  zeng_value = EL_ZENG_CONTACT_LIFT_COEFFICIENT(2.0d0)
  IF (ABS(zeng_value-dns_value)/dns_value.GT.0.06d0) STOP 62
  dns_value = 1.305d0
  zeng_value = EL_ZENG_CONTACT_LIFT_COEFFICIENT(10.0d0)
  IF (ABS(zeng_value-dns_value)/dns_value.GT.0.06d0) STOP 63
  dns_value = 0.3384d0
  zeng_value = EL_ZENG_CONTACT_LIFT_COEFFICIENT(200.0d0)
  IF (ABS(zeng_value-dns_value)/dns_value.GT.0.06d0) STOP 64
  IF (ABS(EL_ZENG_SHEAR_LIFT_COEFFICIENT(20.0d0,0.5d0)- &
      EL_ZENG_CONTACT_LIFT_COEFFICIENT(20.0d0)).GT.1.0d-12) STOP 65
  IF (EL_ZENG_SHEAR_LIFT_COEFFICIENT(20.0d0,4.0d0).LE.0.0d0) STOP 66
  IF (EL_ZENG_SHEAR_LIFT_COEFFICIENT(30.0d0,4.0d0).GE.0.0d0) STOP 67

  CALL EL_SET_DOMAIN_BOX(0.0d0,1.0d0,0.0d0,1.0d0,0.0d0,1.0d0)
  sample(EL_SAMPLE_EPSILON_F) = 2.0d0
  el_lift_model = 'saffman_mei_wall'
  state(1:3) = (/0.5d0,0.5d0,0.75d0/)

  state(4:6) = 0.0d0
  sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) = &
    2.0d0*(/0.05d0, 0.0d0, 0.0d0/)
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  lift_low = force_result%lift
  sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) = &
    2.0d0*(/0.0500001d0, 0.0d0, 0.0d0/)
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  IF (SQRT(SUM((force_result%lift-lift_low)**2)).GT.1.0d-10) STOP 68

  sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) = &
    2.0d0*(/0.2d0, 0.0d0, 0.0d0/)
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  lift_high = force_result%lift
  expected = 0.0d0
  expected(3) = -EL_ZENG_SHEAR_LIFT_COEFFICIENT(2.0d0,25.0d0)* &
                0.5d0*1.0d0*0.2d0*0.2d0*area
  relative_error = SQRT(SUM((lift_high-expected)**2))/ &
                   MAX(SQRT(SUM(expected**2)),TINY(1.0d0))
  IF (relative_error.GT.1.0d-12) STOP 69

  ! Matas-Asmolov inertial lift: table anchors, interpolation, zero crossing,
  ! and radial direction (F = fhat(r/R)*rho*Umax^2*a^4/(8*sqrt(2)*R^2)*e_r).
  el_inertial_lift = 'matas_asmolov'
  el_inertial_lift_umax = 0.6d0
  el_cylinder_center = 0.0d0
  el_cylinder_radius = 0.5d0
  el_cylinder_axis = 'z'
  state = 0.0d0
  state(7) = 0.025d0

  ! Grid point s = 0.40 (fhat = 4.20), particle on +x: outward push.
  state(1:3) = (/0.2d0, 0.0d0, 1.0d0/)
  CALL EL_INERTIAL_LIFT_FORCE(state, 1.0d0, force)
  expected = 0.0d0
  expected(1) = 4.20d0*1.0d0*0.36d0*0.025d0**4/ &
                (8.0d0*SQRT(2.0d0)*0.25d0)
  relative_error = SQRT(SUM((force-expected)**2))/SQRT(SUM(expected**2))
  IF (relative_error.GT.1.0d-12) STOP 74

  ! s = 0.80 (fhat = -5.20): inward push, placed on -y so e_r = (0,-1,0).
  state(1:3) = (/0.0d0, -0.4d0, 1.0d0/)
  CALL EL_INERTIAL_LIFT_FORCE(state, 1.0d0, force)
  IF (force(2).LE.0.0d0) STOP 75
  IF (ABS(force(1)).GT.1.0d-18 .OR. ABS(force(3)).GT.1.0d-18) STOP 75

  ! Zero crossing of the interpolated table: between s = 0.65 (+0.80) and
  ! s = 0.70 (-0.80) the linear interpolant vanishes at s = 0.675.
  state(1:3) = (/0.3375d0, 0.0d0, 1.0d0/)
  CALL EL_INERTIAL_LIFT_FORCE(state, 1.0d0, force)
  IF (SQRT(SUM(force**2)).GT.1.0d-18) STOP 76

  ! Off-grid interpolation at s = 0.42: 0.6*4.20 + 0.4*3.90 = 4.08.
  state(1:3) = (/0.21d0, 0.0d0, 1.0d0/)
  CALL EL_INERTIAL_LIFT_FORCE(state, 1.0d0, force)
  expected(1) = 4.08d0*1.0d0*0.36d0*0.025d0**4/ &
                (8.0d0*SQRT(2.0d0)*0.25d0)
  relative_error = ABS(force(1)-expected(1))/expected(1)
  IF (relative_error.GT.1.0d-12) STOP 77
  el_inertial_lift = 'none'
  el_inertial_lift_umax = 0.0d0
  el_cylinder_radius = -1.0d0

  CALL EL_INITIALIZE(3,5)
  el_field_data%alpha_p = (/0.1d0,0.2d0,0.3d0/)
  el_field_data%epsilon_f = (/0.9d0,0.8d0,0.7d0/)
  el_field_data%epsilon_f_old = (/0.95d0,0.85d0,0.75d0/)
  el_field_data%deps_f_dt = (/-0.5d0,-0.5d0,-0.5d0/)
  el_field_data%force_rhs = RESHAPE((/(DBLE(i),i=1,15)/),(/3,5/))
  el_field_data%drag_B_ml = (/(0.25d0*DBLE(i),i=1,5)/)
  el_field_data%drag_B_source = (/(0.5d0*DBLE(i),i=1,5)/)
  CALL EL_WRITE_RESTART('test_el_fields.restart')
  CALL EL_FINALIZE()
  CALL EL_READ_RESTART('test_el_fields.restart')
  IF (MAXVAL(ABS(el_field_data%epsilon_f_old- &
      (/0.95d0,0.85d0,0.75d0/))).GT.1.0d-14) STOP 5
  IF (MAXVAL(ABS(el_field_data%deps_f_dt+0.5d0)).GT.1.0d-14) STOP 6
  IF (MAXVAL(ABS(el_field_data%drag_B_ml- &
      (/(0.25d0*DBLE(i),i=1,5)/))).GT.1.0d-14) STOP 72
  IF (MAXVAL(ABS(el_field_data%drag_B_source- &
      (/(0.5d0*DBLE(i),i=1,5)/))).GT.1.0d-14) STOP 73
  CALL EL_FINALIZE()
  OPEN(NEWUNIT=restart_unit,FILE='test_el_fields.restart',STATUS='OLD')
  CLOSE(restart_unit,STATUS='DELETE')

  WRITE(*,*) 'EL kernel and force tests passed.'

END PROGRAM TEST_EL_KERNEL_FORCES
