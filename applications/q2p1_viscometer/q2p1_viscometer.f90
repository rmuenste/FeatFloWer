!===============================================================================
! q2p1_viscometer - D5.1 numerical viscometer
!
! Searle-type rotational viscometer: stationary outer vessel, rotating
! inner bob.  The bob is NOT meshed and NOT an FBM body - it is a blind
! coaxial bore in the mesh whose surfaces carry the rotating Dirichlet
! condition u = Omega x r (boundary type 'Inflow14', see
! source/src_quadLS/QuadSc_user.f90).  The same bore surfaces are tagged
! 'WallF' by an overlay .par so that they, and only they, form the
! boundary-force mask BndrForce used by both torque estimators.
!
! Observable: eta_rel(phi) = T_z(phi)/T_z(0) at fixed Omega.
!
! Two independent torque estimators are printed every time step, raw
! (unnormalised - normalisation is done offline):
!
!   VISC_TORQUE_RES time= <t> Tz= <value>   primary, residual-based
!                                           (EvaluateTorque_residual,
!                                            BenchForce-style reaction)
!   VISC_TORQUE_DNA time= <t> Tz= <value>   secondary, grad(alpha)-band
!                                           volume integral (GetVolumeFormTorque)
!
! Both calls live here, in the application, on purpose: the shared
! per-step code path in QuadSc_main.f90 is left untouched so that no
! other application changes behaviour.
!
! Cloned from applications/q2p1_dkt (serial PE + HardContactAndFluid +
! external particles.xyz packing).
!===============================================================================
PROGRAM Q2P1_VISCOMETER

  include 'defs_include.h'

  use solution_io, only: postprocessing_app

  use app_initialization, only: init_q2p1_app

  ! D5.1 boundary-torque estimators (both are app-local calls)
  use Transport_Q2P1, only: VISC_GetTorqueResidual, VISC_GetTorqueVolume

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile, ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0

  !-------INIT PHASE-------

  call init_q2p1_app(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dout = Real(INT(timens/dtgmv)+1)*dtgmv

  !-------MAIN LOOP-------

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
  CALL Transport_q2p1_UxyzP_fc_ext(ufile,inonln_u,itns)

  inonln_t = 2

  !-------D5.1 BOUNDARY TORQUE (both estimators, raw values)-------
  ! Evaluated after the converged NS step, on the BndrForce ('WallF')
  ! masked bore surface.  Newtonian only (deck: SimPar@FlowType = Newtonian).
  CALL VISC_GetTorqueResidual(ufile)
  CALL VISC_GetTorqueVolume(ufile)

  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1

  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  call sim_finalize(tt0,ufile)

END PROGRAM Q2P1_VISCOMETER
