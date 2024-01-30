PROGRAM Q2P1_SSE

  include 'defs_include.h'

  use solution_io, only: postprocessing_sse

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize_sse

  use Transport_q2p1, only : Transport_q2p1_UxyzP_sse,DetermineIfGoalsWereReached
  use Sigma_User, only : bKTPRelease,mySetup
  use var_QuadScalar, only : SSE_HAS_ANGLE, extruder_angle,DivergedSolution,myErrorCode
  use var_QuadScalar, only : myErrorCode
  USE PP3D_MPI, ONLY: MPI_COMM_WORLD
  use f90getopt
 
  
  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  character(len=60)  :: arg
  real               :: dout = 0.0
  integer            :: ufile,ilog,i,ierr
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  LOGICAL bGoalsReached

  character(len=*), parameter :: version = '1.0'
  real*8                      :: angle = 0.0
  integer                     :: iangle = 0
  type(option_s)              :: opts(3)

  integer id,mem
 
  opts(1) = option_s('angle', .true.,  'a')
  opts(2) = option_s('version',  .false., 'v')
  opts(3) = option_s('help',  .false., 'h')

  ! check the command line arguments
  do
      select case (getopt('a:vh', opts))
          case (char(0))
              exit
          case ('a')
              read(optarg,*) angle
              iangle = int(angle)
              extruder_angle = angle
              SSE_HAS_ANGLE=.true.
!               write(*,*)'got angle', extruder_angle
          case ('-a')
              read(optarg,*) angle
              iangle = int(angle)
              extruder_angle = angle
              SSE_HAS_ANGLE=.true.
!               write(*,*)'got angle ifc', extruder_angle
          case ('v')
              print '(a, f3.1)', 'version ', version
              call exit(0)
          case ('h')
            call print_help()
            call exit(0)
          case default
      end select
  end do
  !-------INIT PHASE-------

  call init_q2p1_ext(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dout = Real(INT(timens/dtgmv)+1)*dtgmv

  !-------MAIN LOOP-------
#if !defined WIN32
  IF (bKTPRelease) CALL xSEND_START()
#endif

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  if (itns.eq.1) CALL MemoryPrint(1,'w','t0')
  
  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
  call Transport_q2p1_UxyzP_sse(ufile,inl_u, itns)

  if (DivergedSolution .eqv. .true.) EXIT
  
  if (mySetup%bPressureConvergence) THEN
   timens=timens+dtgmv
   call postprocessing_sse(dout, inonln_u, inonln_t,ufile)
   exit
  end if

  IF (bTracer) THEN
    ! Solve transport equation for linear scalar
    CALL Transport_LinScalar(ufile,inonln_t)
  ELSE
    inonln_t = 2
  END IF

  CALL MemoryPrint(1,'w','tnn')
  
  call postprocessing_sse(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT
  
  END DO

#if !defined WIN32
  IF (bKTPRelease) CALL xSEND_FINISH()
#endif

  if (DivergedSolution) call MPI_Abort(MPI_COMM_WORLD, myErrorCode%DIVERGENCE_U, ierr)
  CALL DetermineIfGoalsWereReached(bGoalsReached)
  if (.not.bGoalsReached)  call MPI_Abort(MPI_COMM_WORLD, myErrorCode%RUNOUTOFTIMESTEPS, ierr)
  
  call sim_finalize_sse(tt0,ufile)
 
  contains

    subroutine print_help()
        print '(a, /)', 'command-line options:'
        print '(a)',    '  -v, --version     print version information and exit'
        print '(a)',    '  -h, --help        print usage information and exit'
        print '(a)',    '  -a, --angle       the angle position for the simulation'
    end subroutine print_help  

END PROGRAM Q2P1_SSE
