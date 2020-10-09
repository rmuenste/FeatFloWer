PROGRAM HEAT

  include 'defs_include.h'

  use solution_io, only: postprocessing_app_heat,write_sol_to_file

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize
                         
  use Transport_Q2P1,  only:        updateFBMGeometry       
  use Transport_Q1, ONLY : AddSource_EWIKON,Boundary_LinSc_Val_EWIKON,Transport_LinScalar_EWIKON
  use var_QuadScalar, only : DivergedSolution,ConvergedSolution
  USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI


  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog,iXX
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0

  !-------INIT PHASE-------

  call init_q2p1_ext(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dout = Real(INT(timens/dtgmv)+1)*dtgmv

!   CALL updateFBMGeometry()
  !-------MAIN LOOP-------

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt
  inonln_u = 2
  inonln_t = 2

  ! Solve transport equation for linear scalar
  CALL Transport_LinScalar_EWIKON(Boundary_LinSc_Val_EWIKON,AddSource_EWIKON,ufile,inonln_t)

  call postprocessing_app_heat(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

   if (DivergedSolution) EXIT
   
   IF (ConvergedSolution) EXIT
   
  END DO

  
  if (.not.DivergedSolution) THEN
  
   IF (ConvergedSolution) THEN
    IF (myid.eq.showid) THEN
      WRITE(*,*) "HEAT_APP has successfully finished before reaching timestep limit. "
      WRITE(ufile,*) "HEAT_APP has successfully finished before reaching timestep limit. "
    END IF
   ELSE
    IF (myid.eq.showid) THEN
      WRITE(*,*) "HEAT_APP has successfully finished. "
      WRITE(ufile,*) "HEAT_APP has successfully finished. "
    END IF
   END IF   
   
  ELSE
   IF (myid.eq.showid) THEN
     WRITE(*,*) "HEAT_APP has been stopped due to recognition of a diverged solution ..."
     WRITE(ufile,*) "HEAT_APP has been stopped due to recognition of a diverged solution ..."
   END IF
  END IF

  CALL Barrier_myMPI()
  CALL MPI_Finalize(ierr)
  
  if (DivergedSolution) STOP 1

END PROGRAM HEAT
