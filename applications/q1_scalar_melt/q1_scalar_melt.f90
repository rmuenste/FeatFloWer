PROGRAM Q1_GenScalar

  include 'defs_include.h'

  use solution_io, only: postprocessing_app
  use solution_io, only: postprocessing_sse_q1_scalar

  use Transport_Q1, only : Reinit_GenLinSc_Q1,Correct_GenLinSc_Q1_ALPHA,EstimateAlphaTimeStepSize
  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize_sse

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  logical            :: bRestartTime =.true.
  real*8             :: Orig_tsep

  !-------INIT PHASE-------

  call init_q1_scalar(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  IF (bRestartTime) THEN
   istep_ns = 0
   timens = 0d0
  END IF
  
  Orig_tsep = tstep
  
  dout = Real(INT(timens/dtgmv)+1)*dtgmv
  
  !-------MAIN LOOP-------

  DO itns=1,nitns

1 CONTINUE  

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  ! Solve transport equation for linear scalar
  CALL Transport_GenLinSc_Q1_Melt(ufile,inonln_t)

  call postprocessing_sse_q1_scalar(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

!   write(*,*) '0',bAlphaConverged
!  CALL Correct_GenLinSc_Q1_ALPHA(ufile)
  timens = timens + dtgmv
  itns = max(itns,2)
  call postprocessing_sse_q1_scalar(dout, inonln_u, inonln_t,ufile)

2 CONTINUE
  
  call sim_finalize_sse(tt0,ufile)

END PROGRAM Q1_GenScalar
