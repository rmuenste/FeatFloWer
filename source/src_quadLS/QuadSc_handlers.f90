! ----------------------------------------------
!
SUBROUTINE Init_QuadScalar(mfile)
INTEGER I,J,ndof,mfile

 ! Init pe handlers
 call Init_PE_Handlers()

 QuadSc%cName = "Velo"
 LinSc%cName = "Pres"

 CALL GetVeloParameters(QuadSc%prm,QuadSc%cName,mfile)
 CALL GetPresParameters(LinSc%prm,LinSc%cName,mfile)

END SUBROUTINE Init_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE Init_Die_Handlers()
! In this function we set the function handlers
! for Die FBM simulations --> as normal but return surface distances as well
implicit none

 fbm_geom_handler_ptr => GetFictKnpr_Die
 fbm_vel_bc_handler_ptr => FictKnpr_velBC
 fbm_vel_bc_handler_ptr => fbm_velBC

END SUBROUTINE Init_Die_Handlers
!
! ----------------------------------------------
!
SUBROUTINE Init_Wangen_Handlers()
! In this function we set the function handlers
! for Die FBM simulations --> as normal but return surface distances as well
implicit none

 fbm_geom_handler_ptr => GetFictKnpr_Wangen
 fbm_vel_bc_handler_ptr => FictKnpr_velBC_Wangen

END SUBROUTINE Init_Wangen_Handlers
!
! ----------------------------------------------
!
SUBROUTINE Init_Laser_Handlers()
implicit none

 fbm_up_handler_ptr => fbm_updateLaser

END SUBROUTINE Init_Laser_Handlers
!
! ----------------------------------------------
!
SUBROUTINE Init_Default_Handlers()
! In this function we set the function handlers
! for FBM, etc. to their default values
implicit none

 fbm_up_handler_ptr => fbm_updateDefault
 fbm_geom_handler_ptr => fbm_getFictKnpr
 fbm_vel_bc_handler_ptr => fbm_velBC

 fbm_force_handler_ptr => fbm_updateForcesFC2

END SUBROUTINE Init_Default_Handlers
!
! ----------------------------------------------
!
SUBROUTINE Init_PE_Handlers()
! In this function we set the function handlers
! for FBM, etc. to their default values
implicit none

 fbm_up_handler_ptr => fbm_updateDefaultFC2
 fbm_geom_handler_ptr => fbm_getFictKnprFC2
 fbm_vel_bc_handler_ptr => fbm_velBCFC2

 fbm_force_handler_ptr => fbm_updateForcesFC2

END SUBROUTINE Init_PE_Handlers

