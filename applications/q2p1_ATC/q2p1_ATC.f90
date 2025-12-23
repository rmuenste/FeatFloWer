PROGRAM Q2P1_ATC

  ! Include definitions
  include 'defs_include.h'

  use solution_io, only: postprocessing_app

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize
  USE var_QuadScalar, ONLY :  myTimer

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

  !-------SOLVER PHASE-------
  call solve_q2p1_app(ufile)

  !-------FINALIZATION PHASE-------
  call sim_finalize()

END PROGRAM Q2P1_ATC
