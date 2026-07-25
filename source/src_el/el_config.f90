MODULE EL_CONFIG

  IMPLICIT NONE

  CHARACTER(LEN=32) :: el_kernel = 'deen_poly'
  CHARACTER(LEN=32) :: el_drag_model = 'difelice'
  CHARACTER(LEN=32) :: el_drag_coupling = 'explicit'
  CHARACTER(LEN=32) :: el_lift_model = 'none'
  CHARACTER(LEN=32) :: el_inertial_lift = 'none'
  CHARACTER(LEN=32) :: el_prescribed_field = 'none'
  CHARACTER(LEN=32) :: el_domain_type = 'box'
  CHARACTER(LEN=8) :: el_cylinder_axis = 'x'
  REAL*8 :: el_kernel_width_factor = 2.5d0
  REAL*8 :: el_eps_f_min = 0.4d0
  REAL*8 :: el_eps_f_relax = 0.0d0
  REAL*8 :: el_shear_rate = 0.0d0
  REAL*8 :: el_prescribed_umax = 0.0d0
  REAL*8 :: el_inertial_lift_umax = 0.0d0
  REAL*8 :: el_cylinder_center(3) = (/0.0d0,0.0d0,0.0d0/)
  REAL*8 :: el_cylinder_radius = -1.0d0
  LOGICAL :: el_apply_particle_forces = .TRUE.
  LOGICAL :: el_apply_fluid_feedback = .FALSE.
  ! Apply Prop@Gravity as a fluid body force (AddGravForce). Set to No for
  ! fully periodic suspensions: without boundaries there is no hydrostatic
  ! pressure to absorb rho_f*g, so fluid gravity must be dropped (absorbed
  ! into the modified pressure) while particles keep the analytic
  ! grav_buoy = (rho_p-rho_f)*V*g. Default Yes preserves bounded-domain
  ! behavior (tier2 terminal, ten Cate).
  LOGICAL :: el_fluid_gravity = .TRUE.
  LOGICAL :: el_drag_semi_implicit = .FALSE.
  LOGICAL :: el_pressure_force = .TRUE.
  LOGICAL :: el_magnus = .FALSE.
  LOGICAL :: el_write_diagnostics = .TRUE.
  ! Measured-leak momentum compensator (SimPar@ELMomentumFix): each step
  ! the EL_FLUID_PAIR audit measures the fluid-internal momentum error of
  ! the previous fluid solve (dominated by the Galerkin convective form,
  ! -rho*int(u_i div u)); when enabled, its negative is applied as a
  ! uniform body force through the Grav_QuadSc pathway one step lagged
  ! (deadbeat: cumulative drift stays bounded by one step's leak).
  ! Use with the legacy convective form; the momentum-exact divergence
  ! form (ELConvectionForm=divergence) is not energy-robust in developed
  ! pseudo-turbulence (NaN blowup, see v2_rz_settling RUNBOOK part 4).
  LOGICAL :: el_momentum_fix = .FALSE.
  INTEGER :: el_momentum_audit_freq = 100
  REAL*8 :: el_tavg_window = 0.5d0
  ! Mesoscale lane-mode filter for periodic settling boxes: each step,
  ! subtract from u_z its zero-mean x-binned and y-binned plane averages
  ! (the gravest "lane" modes of the two-way-coupled concentration
  ! instability, which lateral container walls would suppress in the
  ! wall-bounded experiments RZ describes). The subtracted field is
  ! z-independent (divergence-free on the axis-aligned box) and has its
  ! global mean removed (momentum-neutral; any lumping residue lands in
  ! the EL_FLUID_PAIR mismatch and is absorbed by ELMomentumFix).
  ! Periodic-box suspension runs only.
  LOGICAL :: el_meso_filter = .FALSE.
  INTEGER :: el_meso_filter_bins = 16

CONTAINS

  SUBROUTINE EL_VALIDATE_CONFIG()

    CALL EL_LOWERCASE(el_kernel)
    CALL EL_LOWERCASE(el_drag_model)
    CALL EL_LOWERCASE(el_drag_coupling)
    CALL EL_LOWERCASE(el_lift_model)
    CALL EL_LOWERCASE(el_inertial_lift)
    CALL EL_LOWERCASE(el_prescribed_field)
    CALL EL_LOWERCASE(el_domain_type)
    CALL EL_LOWERCASE(el_cylinder_axis)

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

    SELECT CASE (TRIM(el_drag_coupling))
    CASE ('explicit')
      el_drag_semi_implicit = .FALSE.
    CASE ('semi_implicit')
      el_drag_semi_implicit = .TRUE.
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELDragCoupling: ', TRIM(el_drag_coupling)
      STOP 1
    END SELECT

    IF (el_drag_semi_implicit .AND. .NOT.el_apply_fluid_feedback) THEN
      WRITE(*,*) 'ELDragCoupling=semi_implicit requires ELApplyFluidFeedback=Yes.'
      STOP 1
    END IF

    SELECT CASE (TRIM(el_lift_model))
    CASE ('none', 'saffman', 'saffman_mei', 'saffman_mei_wall')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELLiftModel: ', TRIM(el_lift_model)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_inertial_lift))
    CASE ('none', 'matas_asmolov')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELInertialLift: ', TRIM(el_inertial_lift)
      STOP 1
    END SELECT

    IF (TRIM(el_inertial_lift).NE.'none') THEN
      IF (TRIM(el_domain_type).NE.'cylinder' .OR. &
          el_cylinder_radius.LE.0.0d0) THEN
        WRITE(*,*) 'ELInertialLift requires ELDomainType=cylinder with ', &
          'ELCylinderRadius > 0.'
        STOP 1
      END IF
      IF (el_inertial_lift_umax.LE.0.0d0 .AND. &
          el_prescribed_umax.LE.0.0d0) THEN
        WRITE(*,*) 'ELInertialLift requires ELInertialLiftUmax > 0 ', &
          '(or ELPrescribedUmax > 0 in frozen runs).'
        STOP 1
      END IF
    END IF

    SELECT CASE (TRIM(el_prescribed_field))
    CASE ('none', 'linear_shear', 'poiseuille')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELPrescribedField: ', TRIM(el_prescribed_field)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_domain_type))
    CASE ('box', 'cylinder')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELDomainType: ', TRIM(el_domain_type)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_cylinder_axis))
    CASE ('x', 'y', 'z')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELCylinderAxis: ', TRIM(el_cylinder_axis)
      STOP 1
    END SELECT

    IF (TRIM(el_domain_type).EQ.'cylinder' .AND. &
        el_cylinder_radius.LE.0.0d0) THEN
      WRITE(*,*) 'ELDomainType=cylinder requires ELCylinderRadius > 0.'
      STOP 1
    END IF

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

    IF (el_tavg_window.LE.0.0d0 .OR. el_tavg_window.GT.1.0d0) THEN
      WRITE(*,*) 'ELTAvgWindow must be in (0,1].'
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
    WRITE(mfile,'(A,A)') 'EL drag coupling          = ', TRIM(el_drag_coupling)
    WRITE(mfile,'(A,L1)') 'EL pressure force         = ', el_pressure_force
    WRITE(mfile,'(A,A)') 'EL lift model             = ', TRIM(el_lift_model)
    WRITE(mfile,'(A,A)') 'EL inertial lift          = ', TRIM(el_inertial_lift)
    WRITE(mfile,'(A,ES14.6)') 'EL inertial lift umax     = ', el_inertial_lift_umax
    WRITE(mfile,'(A,A)') 'EL prescribed field       = ', TRIM(el_prescribed_field)
    WRITE(mfile,'(A,ES14.6)') 'EL shear rate             = ', el_shear_rate
    WRITE(mfile,'(A,ES14.6)') 'EL prescribed umax        = ', el_prescribed_umax
    WRITE(mfile,'(A,A)') 'EL domain type            = ', TRIM(el_domain_type)
    WRITE(mfile,'(A,3ES14.6)') 'EL cylinder center        = ', el_cylinder_center
    WRITE(mfile,'(A,ES14.6)') 'EL cylinder radius        = ', el_cylinder_radius
    WRITE(mfile,'(A,A)') 'EL cylinder axis          = ', TRIM(el_cylinder_axis)
    WRITE(mfile,'(A,L1)') 'EL Magnus                 = ', el_magnus
    WRITE(mfile,'(A,L1)') 'EL apply particle forces  = ', el_apply_particle_forces
    WRITE(mfile,'(A,L1)') 'EL apply fluid feedback   = ', el_apply_fluid_feedback
    WRITE(mfile,'(A,L1)') 'EL fluid gravity          = ', el_fluid_gravity
    WRITE(mfile,'(A,L1)') 'EL momentum fix           = ', el_momentum_fix
    WRITE(mfile,'(A,L1,A,I0,A)') 'EL meso filter            = ', el_meso_filter, &
      ' (bins ', el_meso_filter_bins, ')'
    WRITE(mfile,'(A,I0)') 'EL momentum audit freq    = ', el_momentum_audit_freq
    WRITE(mfile,'(A,ES14.6)') 'EL tavg window            = ', el_tavg_window

    IF (mterm.NE.mfile) THEN
      WRITE(mterm,'(A,A)') 'EL kernel                 = ', TRIM(el_kernel)
      WRITE(mterm,'(A,ES14.6)') 'EL kernel width factor    = ', el_kernel_width_factor
      WRITE(mterm,'(A,ES14.6)') 'EL epsilon minimum        = ', el_eps_f_min
      WRITE(mterm,'(A,ES14.6)') 'EL epsilon relaxation     = ', el_eps_f_relax
      WRITE(mterm,'(A,A)') 'EL drag model             = ', TRIM(el_drag_model)
      WRITE(mterm,'(A,A)') 'EL drag coupling          = ', TRIM(el_drag_coupling)
      WRITE(mterm,'(A,L1)') 'EL pressure force         = ', el_pressure_force
      WRITE(mterm,'(A,A)') 'EL lift model             = ', TRIM(el_lift_model)
      WRITE(mterm,'(A,A)') 'EL inertial lift          = ', TRIM(el_inertial_lift)
      WRITE(mterm,'(A,ES14.6)') 'EL inertial lift umax     = ', el_inertial_lift_umax
      WRITE(mterm,'(A,A)') 'EL prescribed field       = ', TRIM(el_prescribed_field)
      WRITE(mterm,'(A,ES14.6)') 'EL shear rate             = ', el_shear_rate
      WRITE(mterm,'(A,ES14.6)') 'EL prescribed umax        = ', el_prescribed_umax
      WRITE(mterm,'(A,A)') 'EL domain type            = ', TRIM(el_domain_type)
      WRITE(mterm,'(A,3ES14.6)') 'EL cylinder center        = ', el_cylinder_center
      WRITE(mterm,'(A,ES14.6)') 'EL cylinder radius        = ', el_cylinder_radius
      WRITE(mterm,'(A,A)') 'EL cylinder axis          = ', TRIM(el_cylinder_axis)
      WRITE(mterm,'(A,L1)') 'EL Magnus                 = ', el_magnus
      WRITE(mterm,'(A,L1)') 'EL apply particle forces  = ', el_apply_particle_forces
      WRITE(mterm,'(A,L1)') 'EL apply fluid feedback   = ', el_apply_fluid_feedback
      WRITE(mterm,'(A,L1)') 'EL fluid gravity          = ', el_fluid_gravity
      WRITE(mterm,'(A,L1)') 'EL momentum fix           = ', el_momentum_fix
      WRITE(mterm,'(A,L1,A,I0,A)') 'EL meso filter            = ', el_meso_filter, &
        ' (bins ', el_meso_filter_bins, ')'
      WRITE(mterm,'(A,I0)') 'EL momentum audit freq    = ', el_momentum_audit_freq
      WRITE(mterm,'(A,ES14.6)') 'EL tavg window            = ', el_tavg_window
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
