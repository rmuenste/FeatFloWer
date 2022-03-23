PROGRAM Q1_GenScalar

  include 'defs_include.h'

  use solution_io, only: postprocessing_general
  use var_QuadScalar, only: bAlphaConverged,DivergedSolution,myErrorCode,AlphaControl
  USE PP3D_MPI, ONLY: MPI_COMM_WORLD

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
  integer            :: iChange,nChange=8

  !-------INIT PHASE-------

  call init_q1_scalar(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  IF (bRestartTime) THEN
   istep_ns = 0
   timens = 0d0
  END IF
  
  CALL EstimateAlphaTimeStepSize(ufile)
  Orig_tsep = tstep
  
  dout = Real(INT(timens/dtgmv)+1)*dtgmv
  
  !-------MAIN LOOP-------

  bAlphaConverged = .false.
  DivergedSolution = .false.
  
  iChange = 0
  
  DO itns=1,nitns

1 CONTINUE  

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  ! Solve transport equation for linear scalar
  CALL Transport_GenLinSc_Q1_Multimat(ufile,inonln_t)

  if (DivergedSolution .eqv. .true.) THEN
   timens=timens-dt
   tstep = tstep*0.5d0
   if (tstep.lt.Orig_tsep*(2d0/9d0)) then
    if (myid.eq.1) write(MTERM,*) 'Q1_scalar simulation has diverged after multiple timestep reductions!'
    if (myid.eq.1) write(ufile,*) 'Q1_scalar simulation has diverged after multiple timestep reductions!'
    GOTO 2
   ELSE
    DivergedSolution = .false.
    if (myid.eq.1) write(MTERM,*) 'Timestep reduction due to divergence in Q1_scalar solver '
    if (myid.eq.1) write(ufile,*) 'Timestep reduction due to divergence in Q1_scalar solver '
    GOTO 1
   END IF
  END IF
  
  !!!!!!!!!!!!!!!! TimestepControl   !!!!!!!!!!!!!!!
  if (inonln_t.lt.3.and.iChange.gt.nChange) then
   tstep = min((3d0/2d0)*tstep,Orig_tsep*(3d0/1d0))
   iChange = -1
  end if
  if (inonln_t.gt.5.and.iChange.gt.nChange) then
   tstep = max((2d0/3d0)*tstep,Orig_tsep*(1d0/3d0))
   iChange = -1
  end if
  iChange = iChange + 1
  !!!!!!!!!!!!!!!! TimestepControl   !!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!! ConvergenceControl   !!!!!!!!!!!!!!!
  IF (itns.gt.INT(0.5d0*DBLE(nitns)).and.AlphaControl.lt.0.5d0) THEN
   DivergedSolution = .TRUE.
   GOTO 2
  END IF
  !!!!!!!!!!!!!!!! ConvergenceControl   !!!!!!!!!!!!!!!

  call postprocessing_general(dout, inonln_u, inonln_t,ufile,'v,p,q,t')
!   call postprocessing_sse_q1_scalar(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  
  ! Exit if done
  IF (bAlphaConverged) EXIT
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  if (.not.bAlphaConverged)  call MPI_Abort(MPI_COMM_WORLD, myErrorCode%RUNOUTOFTIMESTEPS, ierr)
  
!   write(*,*) '0',bAlphaConverged
  CALL Correct_GenLinSc_Q1_ALPHA(ufile)
  timens = timens + dtgmv
  itns = max(itns,2)
  call postprocessing_general(dout, inonln_u, inonln_t,ufile,'v,p,q,t')
!   call postprocessing_sse_q1_scalar(dout, inonln_u, inonln_t,ufile)

2 CONTINUE

  if (DivergedSolution)  call MPI_Abort(MPI_COMM_WORLD, myErrorCode%DIVERGENCE_ALPHA, ierr)
  
  call sim_finalize_sse(tt0,ufile)

END PROGRAM Q1_GenScalar
