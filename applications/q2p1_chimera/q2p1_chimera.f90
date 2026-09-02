PROGRAM Q2P1_CHIMERA

  ! Chimera overlapping-mesh application (chimera-integration-design.md
  ! v3, section 8): the standard Q2/P1 flow solve (shared fluid core,
  ! FBM path disabled) with the particle handled by the Chimera
  ! component.  Hooks H5/H6 (Chimera_Initialize/Chimera_Finalize) are
  ! app-local by design; the per-step coupling (H1) and the hole/fringe
  ! constraints (H2-H4) live behind the CHIMERA_API facade inside the
  ! shared solver.  With SimPar@ChimeraEnable = No this program is a
  ! plain flow solver (no FBM, no forces) - the reference application
  ! for the standard mode remains q2p1_fc_ext.

  ! Include definitions
  include 'defs_include.h'

  use solution_io, only: postprocessing_app

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize
  USE var_QuadScalar, ONLY :  myTimer
  USE Transport_Q2P1, ONLY : Transport_q2p1_UxyzP_fluid_core
  USE chimera_api, ONLY : Chimera_Initialize, Chimera_Finalize

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  real               :: dtt10 = 0.0

  character(len=100) :: arg
  character(len=100) :: version_string
  character(len=100) :: git_commit_hash_trim
  logical :: show_version
#include "./version.h"

  !-------INIT PHASE-------

  ! Read command line arguments
  if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      show_version = trim(arg) == "-v"
  else
      show_version = .false.
  endif

  ! Display version information if "-v" argument is given and exit
  if (show_version) then
      version_string = "Version: " // trim(PROJECT_VERSION)
      git_commit_hash_trim = "Git Commit Hash: " // trim(GIT_COMMIT_HASH)
      print *, version_string
      print *, git_commit_hash_trim
      stop
  endif

  call init_q2p1_app(ufile)

  ! Hook H5 (app-local): submeshes, markers, donor caches.
  CALL Chimera_Initialize(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dout = Real(INT(timens/dtgmv)+1)*dtgmv

  !-------MAIN LOOP-------

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  ! Solve Navier-Stokes: shared fluid core, FBM disabled (the body is
  ! represented by the Chimera hole/fringe constraints).
  CALL Transport_q2p1_UxyzP_fluid_core(ufile,inonln_u,itns,.FALSE.)

  IF (bTracer) THEN
    ! Solve transport equation for linear scalar
    CALL Transport_LinScalar(ufile,inonln_t)
  ELSE
    inonln_t = 2
  END IF

  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  ! Hook H6 (app-local): release the component (idempotent).
  CALL Chimera_Finalize()

  CALL Output_MPI_Timings()

  call sim_finalize(tt0,ufile)


END PROGRAM Q2P1_CHIMERA
