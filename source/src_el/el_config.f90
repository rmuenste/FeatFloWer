MODULE EL_CONFIG

  IMPLICIT NONE

  CHARACTER(LEN=32) :: el_kernel = 'deen_poly'
  CHARACTER(LEN=32) :: el_drag_model = 'difelice'
  CHARACTER(LEN=32) :: el_lift_model = 'none'
  CHARACTER(LEN=32) :: el_prescribed_field = 'none'
  REAL*8 :: el_kernel_width_factor = 2.5d0
  REAL*8 :: el_eps_f_min = 0.4d0
  REAL*8 :: el_eps_f_relax = 0.0d0
  REAL*8 :: el_shear_rate = 0.0d0
  LOGICAL :: el_apply_particle_forces = .TRUE.
  LOGICAL :: el_apply_fluid_feedback = .FALSE.
  LOGICAL :: el_pressure_force = .TRUE.
  LOGICAL :: el_magnus = .FALSE.
  LOGICAL :: el_write_diagnostics = .TRUE.
  INTEGER :: el_momentum_audit_freq = 100

CONTAINS

  SUBROUTINE EL_VALIDATE_CONFIG()

    CALL EL_LOWERCASE(el_kernel)
    CALL EL_LOWERCASE(el_drag_model)
    CALL EL_LOWERCASE(el_lift_model)
    CALL EL_LOWERCASE(el_prescribed_field)

    SELECT CASE (TRIM(el_kernel))
    CASE ('deen_poly', 'gaussian')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELKernel: ', TRIM(el_kernel)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_drag_model))
    CASE ('difelice', 'stokes', 'schiller_naumann')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELDragModel: ', TRIM(el_drag_model)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_lift_model))
    CASE ('none', 'saffman', 'saffman_mei', 'saffman_mei_wall')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELLiftModel: ', TRIM(el_lift_model)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_prescribed_field))
    CASE ('none', 'linear_shear')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELPrescribedField: ', TRIM(el_prescribed_field)
      STOP 1
    END SELECT

    IF (el_magnus) THEN
      WRITE(*,*) 'ELMagnus is not implemented in Phase 2.'
      STOP 1
    END IF

    IF (el_kernel_width_factor.LE.0.0d0) THEN
      WRITE(*,*) 'ELKernelWidthFactor must be positive.'
      STOP 1
    END IF

    IF (el_eps_f_min.LE.0.0d0 .OR. el_eps_f_min.GT.1.0d0) THEN
      WRITE(*,*) 'ELEpsFMin must be in (0,1].'
      STOP 1
    END IF

    IF (el_eps_f_relax.LT.0.0d0 .OR. el_eps_f_relax.GE.1.0d0) THEN
      WRITE(*,*) 'ELEpsFRelax must be in [0,1).'
      STOP 1
    END IF

    IF (el_momentum_audit_freq.LE.0) THEN
      WRITE(*,*) 'ELMomentumAuditFreq must be positive.'
      STOP 1
    END IF

  END SUBROUTINE EL_VALIDATE_CONFIG

  SUBROUTINE EL_PRINT_CONFIG(mfile, mterm)

    INTEGER, INTENT(IN) :: mfile, mterm

    WRITE(mfile,'(A,A)') 'EL kernel                 = ', TRIM(el_kernel)
    WRITE(mfile,'(A,ES14.6)') 'EL kernel width factor    = ', el_kernel_width_factor
    WRITE(mfile,'(A,ES14.6)') 'EL epsilon minimum        = ', el_eps_f_min
    WRITE(mfile,'(A,ES14.6)') 'EL epsilon relaxation     = ', el_eps_f_relax
    WRITE(mfile,'(A,A)') 'EL drag model             = ', TRIM(el_drag_model)
    WRITE(mfile,'(A,L1)') 'EL pressure force         = ', el_pressure_force
    WRITE(mfile,'(A,A)') 'EL lift model             = ', TRIM(el_lift_model)
    WRITE(mfile,'(A,A)') 'EL prescribed field       = ', TRIM(el_prescribed_field)
    WRITE(mfile,'(A,ES14.6)') 'EL shear rate             = ', el_shear_rate
    WRITE(mfile,'(A,L1)') 'EL Magnus                 = ', el_magnus
    WRITE(mfile,'(A,L1)') 'EL apply particle forces  = ', el_apply_particle_forces
    WRITE(mfile,'(A,L1)') 'EL apply fluid feedback   = ', el_apply_fluid_feedback
    WRITE(mfile,'(A,I0)') 'EL momentum audit freq    = ', el_momentum_audit_freq

    IF (mterm.NE.mfile) THEN
      WRITE(mterm,'(A,A)') 'EL kernel                 = ', TRIM(el_kernel)
      WRITE(mterm,'(A,ES14.6)') 'EL kernel width factor    = ', el_kernel_width_factor
      WRITE(mterm,'(A,ES14.6)') 'EL epsilon minimum        = ', el_eps_f_min
      WRITE(mterm,'(A,ES14.6)') 'EL epsilon relaxation     = ', el_eps_f_relax
      WRITE(mterm,'(A,A)') 'EL drag model             = ', TRIM(el_drag_model)
      WRITE(mterm,'(A,L1)') 'EL pressure force         = ', el_pressure_force
      WRITE(mterm,'(A,A)') 'EL lift model             = ', TRIM(el_lift_model)
      WRITE(mterm,'(A,A)') 'EL prescribed field       = ', TRIM(el_prescribed_field)
      WRITE(mterm,'(A,ES14.6)') 'EL shear rate             = ', el_shear_rate
      WRITE(mterm,'(A,L1)') 'EL Magnus                 = ', el_magnus
      WRITE(mterm,'(A,L1)') 'EL apply particle forces  = ', el_apply_particle_forces
      WRITE(mterm,'(A,L1)') 'EL apply fluid feedback   = ', el_apply_fluid_feedback
      WRITE(mterm,'(A,I0)') 'EL momentum audit freq    = ', el_momentum_audit_freq
    END IF

  END SUBROUTINE EL_PRINT_CONFIG

  SUBROUTINE EL_LOWERCASE(text)

    CHARACTER(LEN=*), INTENT(INOUT) :: text
    INTEGER :: i, code

    DO i=1,LEN_TRIM(text)
      code = IACHAR(text(i:i))
      IF (code.GE.IACHAR('A') .AND. code.LE.IACHAR('Z')) THEN
        text(i:i) = ACHAR(code + IACHAR('a') - IACHAR('A'))
      END IF
    END DO

  END SUBROUTINE EL_LOWERCASE

END MODULE EL_CONFIG
