PROGRAM Q1_GenScalar

  include 'defs_include.h'

  use solution_io, only: postprocessing_app
  use solution_io, only: postprocessing_sse_q1_scalar

  use Transport_Q1, only : Reinit_GenLinSc_Q1
  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  logical            :: bRestartTime =.true.

  !-------INIT PHASE-------

  call init_q1_scalar(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  IF (bRestartTime) THEN
   istep_ns = 0
   timens = 0d0
  END IF
  
  dout = Real(INT(timens/dtgmv)+1)*dtgmv
!   pause
  
  !-------MAIN LOOP-------

!   CALL Reinit_GenLinSc_Q1()
  
  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
!   CALL Transport_q2p1_UxyzP_fc_ext(ufile,inonln_u,itns)
 
  IF (bTracer) THEN
    ! Solve transport equation for linear scalar
    CALL Transport_GenLinSc_Q1_Multimat(ufile,inonln_t)
  ELSE
    inonln_t = 2
  END IF

  call postprocessing_sse_q1_scalar(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  call sim_finalize(tt0,ufile)

END PROGRAM Q1_GenScalar
