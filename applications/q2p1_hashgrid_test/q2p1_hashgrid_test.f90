PROGRAM Q2P1_HASHGRID_GATE
  ! Layer-2 FBM classification gate. Real partitioned mesh, real Q2 DOF sweep, pe's
  ! accelerated point query against the exhaustive one on every DOF of every rank.
  ! There is deliberately NO fluid solve in this program: the only production routine
  ! it borrows from the solver loop is one particle step, and that for a reason given
  ! below. The stages themselves are in hashgrid_gate_driver.f90; README.md explains
  ! the construction they rely on.

  include 'defs_include.h'

  USE PP3D_MPI, ONLY : myid
  use app_initialization, only: init_q2p1_app
  use var_QuadScalar, only: Properties, QuadSc, LinSc
  use Transport_Q2P1, only: fbm_up_handler_ptr
  use fbm, only: fbm_updateFBM
  use hashgrid_gate_driver, only: gate_preflight, gate_run_stages
  use post_utils, only: sim_finalize

  integer :: ufile
  real    :: tt0 = 0.0

  !-------INIT PHASE-------
  ! The shared initialisation: mesh hierarchy, Q2 numbering and comms, FictKNPR /
  ! FictKNPR_uint64 (Init_QuadScalar_Stuctures - the reduced variant does not allocate
  ! them), and the pe world through this app's General_init_ext -> commf2c_hashgrid_test
  ! -> setupHashGridTest, which reads pe_user_config.json with packingMethod_ gate_probe.
  call init_q2p1_app(ufile)

  CALL ZTIME(tt0)

  if (myid == 1) write(*,'(A)') '> q2p1_hashgrid_test: Layer-2 FBM classification gate, no fluid solve'

  ! One particle step and nothing else of the production loop. pe's hash grid activates
  ! inside findContacts, i.e. inside step_simulation(), and fbm_updateFBM is also what
  ! sets collision_pipeline_initialized - without which fbm_getFictKnprFC2 never takes
  ! the accelerated branch and every comparison would be vacuous. Gravity is zero and
  ! every body starts at rest, so the step moves nothing.
  dt = tstep
  timens = timens + dt
  call fbm_updateFBM(Properties%Density(1), tstep, timens, &
                     Properties%Gravity, ufile, myid, &
                     QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                     LinSc%valP(NLMAX)%x, &
                     fbm_up_handler_ptr)

  ! Assert the gate is exercising what it claims to (grid active on every worker, span
  ! as computed, accelerated path on, comparison compiled in). Aborts otherwise.
  call gate_preflight(ufile)

  !-------THE STAGES-------
  call gate_run_stages(timens, ufile)

  call sim_finalize(tt0, ufile)

END PROGRAM Q2P1_HASHGRID_GATE
