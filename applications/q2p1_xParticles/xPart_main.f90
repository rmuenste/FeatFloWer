PROGRAM Q2P1_xParticles
  USE PP3D_MPI, ONLY : Barrier_myMPI,myid

  include 'defs_include.h'

  IMPLICIT NONE
  real               :: tt0 = 0.0
  real               :: tt1 = 0.0
  real*8             :: mpi_tt0 = 0d0
  real*8             :: mpi_tt1 = 0d0
  real*8             :: MPI_WTIME
  integer            :: ufile,ierr
  
  CALL ZTIME(tt0)
  call init_q2p1_xParicles(ufile)
  CALL ZTIME(tt1)
  if (myid.eq.1) write(*,'(A,ES12.4,A)') "initialization took ", dble(tt1-tt0),"s"
  
  
  CALL Barrier_myMPI()
  mpi_tt0 = MPI_WTIME()
  CALL Transport_xParticles_MPI(ufile)
  CALL Barrier_myMPI()
  mpi_tt1 = MPI_WTIME()
  if (myid.eq.1) write(*,'(A,ES12.4,A)') "computation took ", mpi_tt1-mpi_tt0,"s"
  
  
  CALL Barrier_myMPI()
  CALL OUTPUT_xParticleDUMP(0)
  
  CALL OUTPUT_xParticleVTK(0)  
  

  CALL MPI_Finalize(ierr)

END PROGRAM Q2P1_xParticles
