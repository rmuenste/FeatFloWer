PROGRAM TEST_EL_KERNEL_FORCES

  USE EL_CONFIG, ONLY: el_kernel, el_drag_model, el_pressure_force, &
                       el_lift_model
  USE EL_KERNEL_FUNCTIONS, ONLY: EL_KERNEL_VALUE
  USE EL_FORCES, ONLY: tElForceResult, EL_DRAG_FORCE, &
                       EL_COMPUTE_PARTICLE_FORCES, EL_SPHERE_CD
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
  REAL*8 :: reynolds, slip_norm, correction, cd, area, chi, eps_ratio
  REAL*8 :: grad_u(3,3), omega_norm, lift_coeff

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
  expected = 3.0d0*ACOS(-1.0d0)*1.0d-3*0.01d0*carrier_velocity
  relative_error = SQRT(SUM((force-expected)**2))/SQRT(SUM(expected**2))
  IF (relative_error.GT.1.0d-2) STOP 40

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
    cd = EL_SPHERE_CD(reynolds)
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
  reynolds = MAX(1.0d-12, 1.0d0*diameter*slip_norm/dynamic_viscosity)
  chi = 3.7d0 - 0.65d0*EXP(-0.5d0*(1.5d0-LOG10(reynolds))**2)
  eps_ratio = 0.5d0**(2.0d0-chi)
  relative_error = ABS(SQRT(SUM(force_result%drag**2))/ &
    SQRT(SUM(force_eps1**2)) - eps_ratio) / eps_ratio
  IF (relative_error.GT.1.0d-12) STOP 53

  sample(EL_SAMPLE_EPSILON_F) = 2.0d0
  state(4:6) = 0.0d0
  sample(EL_SAMPLE_UF_BEGIN:EL_SAMPLE_UF_END) = 2.0d0*(/1.0d0, 0.0d0, 0.0d0/)
  grad_u = 0.0d0
  grad_u(1,2) = 2.0d0
  sample(EL_SAMPLE_GRAD_U_BEGIN:EL_SAMPLE_GRAD_U_END) = &
    2.0d0*(/grad_u(1,:), grad_u(2,:), grad_u(3,:)/)
  el_lift_model = 'saffman_mei'
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  omega_norm = 2.0d0
  lift_coeff = 1.615d0*diameter*diameter* &
               SQRT(1.0d0*dynamic_viscosity*omega_norm)
  expected = (/0.0d0, lift_coeff, 0.0d0/)
  relative_error = SQRT(SUM((force_result%lift-expected)**2))/ &
                   MAX(SQRT(SUM(expected**2)),TINY(1.0d0))
  IF (relative_error.GT.1.0d-12) STOP 48

  sample(EL_SAMPLE_EPSILON_F) = 2.0d0*0.5d0
  el_lift_model = 'saffman_mei_wall'
  CALL EL_COMPUTE_PARTICLE_FORCES(state, sample, 1.0d0, 1.0d-3, gravity, &
                                  force_result)
  expected = 0.5d0*(/0.0d0, lift_coeff, 0.0d0/)
  relative_error = SQRT(SUM((force_result%lift-expected)**2))/ &
                   MAX(SQRT(SUM(expected**2)),TINY(1.0d0))
  IF (relative_error.GT.1.0d-12) STOP 49

  CALL EL_INITIALIZE(3,5)
  el_field_data%alpha_p = (/0.1d0,0.2d0,0.3d0/)
  el_field_data%epsilon_f = (/0.9d0,0.8d0,0.7d0/)
  el_field_data%epsilon_f_old = (/0.95d0,0.85d0,0.75d0/)
  el_field_data%deps_f_dt = (/-0.5d0,-0.5d0,-0.5d0/)
  el_field_data%force_rhs = RESHAPE((/(DBLE(i),i=1,15)/),(/3,5/))
  CALL EL_WRITE_RESTART('test_el_fields.restart')
  CALL EL_FINALIZE()
  CALL EL_READ_RESTART('test_el_fields.restart')
  IF (MAXVAL(ABS(el_field_data%epsilon_f_old- &
      (/0.95d0,0.85d0,0.75d0/))).GT.1.0d-14) STOP 5
  IF (MAXVAL(ABS(el_field_data%deps_f_dt+0.5d0)).GT.1.0d-14) STOP 6
  CALL EL_FINALIZE()
  OPEN(NEWUNIT=restart_unit,FILE='test_el_fields.restart',STATUS='OLD')
  CLOSE(restart_unit,STATUS='DELETE')

  WRITE(*,*) 'EL kernel and force tests passed.'

END PROGRAM TEST_EL_KERNEL_FORCES
