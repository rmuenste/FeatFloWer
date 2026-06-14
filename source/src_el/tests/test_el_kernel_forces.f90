PROGRAM TEST_EL_KERNEL_FORCES

  USE EL_CONFIG, ONLY: el_kernel, el_drag_model
  USE EL_KERNEL_FUNCTIONS, ONLY: EL_KERNEL_VALUE
  USE EL_FORCES, ONLY: EL_DRAG_FORCE
  USE EL_FIELDS, ONLY: el_field_data, EL_INITIALIZE, EL_WRITE_RESTART, &
                       EL_READ_RESTART, EL_FINALIZE

  IMPLICIT NONE

  INTEGER :: i, restart_unit
  REAL*8 :: distance, width
  REAL*8 :: state(8), carrier_velocity(3), force(3), expected(3)
  REAL*8 :: relative_error

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
