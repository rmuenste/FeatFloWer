PROGRAM Q2P1_xParticles
  USE PP3D_MPI, ONLY : Barrier_myMPI,myid

  include 'defs_include.h'

  IMPLICIT NONE
  real               :: tt0 = 0.0
  real               :: tt1 = 0.0
  integer            :: ufile,ierr
  
  CALL ZTIME(tt0)
  call init_q2p1_xParicles(ufile)
  CALL ZTIME(tt1)
  if (myid.eq.1) write(*,'(A,ES12.4,A)') "initialization took ", dble(tt1-tt0),"s"
  
  
  CALL Barrier_myMPI()
  
  CALL ZTIME(tt0)
  CALL Transport_xParticles_MPI(ufile)
  CALL ZTIME(tt1)
  if (myid.eq.1) write(*,'(A,ES12.4,A)') "computation took ", dble(tt1-tt0),"s"
  
  
  CALL Barrier_myMPI()
  CALL OUTPUT_xParticleVTK(0)  

  CALL MPI_Finalize(ierr)

END PROGRAM Q2P1_xParticles
