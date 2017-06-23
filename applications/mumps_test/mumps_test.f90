PROGRAM MUMPS_TEST

  include 'defs_include.h'

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile=111, uterm,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0

  !-------INIT PHASE-------

  call init_q2p1_ext(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dout = Real(INT(timens/dtgmv)+1)*dtgmv

  !-------MAIN LOOP-------

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
  CALL Transport_q2p1_UxyzP(ufile,inonln_u,itns)


  inonln_t = 2


  call postprocessing_fc_ext(dout, iogmv, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

  call sim_finalize(tt0,ufile)
  

!   IMPLICIT NONE
!   INCLUDE 'mpif.h'
!   INCLUDE 'dmumps_struc.h'
!   TYPE (DMUMPS_STRUC) mumps_par
!   INTEGER IERR, I
!   CHARACTER cFile*(20)
!   integer :: myFile = 55
!   integer :: istat
! 
!   CALL MPI_INIT(IERR)
!   ! Define a communicator for the package.
!   mumps_par%COMM = MPI_COMM_WORLD
! 
!   !  Initialize an instance of the package
!   !  for L U factorization (sym = 0, with working host)
!   mumps_par%JOB = -1
!   mumps_par%SYM = 0
!   mumps_par%PAR = 0
!   CALL DMUMPS(mumps_par)
! 
!   IF (mumps_par%INFOG(1).LT.0) THEN
!     WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
!     "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
!     "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
!     GOTO 500
!   END IF
! 
!   WRITE(cFile,'(A,I2.2,A)')  'input_',mumps_par%MYID,'.txt'
! 
!   OPEN(FILE=ADJUSTL(TRIM(cFile)),UNIT=myFile, &
!   action="read",iostat=istat)
! 
!   !  Define problem on the host (processor 0)
!   IF ( mumps_par%MYID .eq. 0 ) THEN
!     READ(myFile,*) mumps_par%N
!     READ(myFile,*) mumps_par%NZ
!     ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
!     ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
!     
!     !         ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )
!     !mumps_par%A  = 0d0
! 
!     DO I = 1, mumps_par%NZ
!     READ(myFile,*) mumps_par%IRN(I),mumps_par%JCN(I)
!     END DO
!     ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
!     DO I = 1, mumps_par%N
!     READ(myFile,*) mumps_par%RHS(I)
!     !          WRITE(*,*) mumps_par%rhs(I)
!     END DO
!   ELSE
! 
!     READ(myFile,*) mumps_par%NZ_loc
!     ALLOCATE( mumps_par%IRN_loc ( mumps_par%NZ_loc) )
!     ALLOCATE( mumps_par%JCN_loc ( mumps_par%NZ_loc) )
!     ALLOCATE( mumps_par%A_loc( mumps_par%NZ_loc) )
!     DO I = 1, mumps_par%NZ_loc
!     READ(myFile,*) mumps_par%IRN_loc(I),mumps_par%JCN_loc(I),&
!     mumps_par%A_loc(I)
!     END DO
! 
!   END IF
! 
!   CLOSE(myFile)
! 
!   !  Call package for solution
!   MUMPS_PAR%icntl(6)  = 7
!   !     pivot order (automatic)
!   MUMPS_PAR%icntl(7)  = 7
!   !     scaling (automatic)
!   MUMPS_PAR%icntl(8)  = 7
!   !     no transpose
!   MUMPS_PAR%icntl(9)  = 1
!   !     max steps for iterative refinement
!   MUMPS_PAR%icntl(10) = 0
!   !     statistics info
!   MUMPS_PAR%icntl(11) = 0 
!   !     controls parallelism
!   MUMPS_PAR%icntl(12) = 0
!   !     use ScaLAPACK for root node
!   MUMPS_PAR%icntl(13) = 0
!   !     percentage increase in estimated workspace
!   MUMPS_PAR%icntl(14) = 100
! 
!   mumps_par%ICNTL(5)=0
!   mumps_par%ICNTL(18)=3
!   !     mumps_par%WRITE_PROBLEM='matrix'
! 
!   mumps_par%JOB = 6
!   CALL DMUMPS(mumps_par)
! 
!   !       mumps_par%JOB = 2
!   !       CALL DMUMPS(mumps_par)
!   ! 
!   !       mumps_par%JOB = 3
!   !       CALL DMUMPS(mumps_par)
! 
!   IF (mumps_par%INFOG(1).LT.0) THEN
!     WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
!     "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),& 
!     "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
!     GOTO 500
!   END IF
! 
!   !  Solution has been assembled on the host
!   IF ( mumps_par%MYID .eq. 0 ) THEN
!     WRITE( 6, * ) ' Solution is '
!     DO I=1,mumps_par%N
!     WRITE( 6, * ) mumps_par%RHS(I)
!     END DO
!   END IF
! 
!   !  Deallocate user data
!   IF ( mumps_par%MYID .eq. 0 )THEN
!     !         DEALLOCATE( mumps_par%A )
!     DEALLOCATE( mumps_par%RHS )
!   END IF
! 
!   !  Destroy the instance (deallocate internal data structures)
!   mumps_par%JOB = -2
!   CALL DMUMPS(mumps_par)
!   IF (mumps_par%INFOG(1).LT.0) THEN
!     WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
!     "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),& 
!     "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
!     GOTO 500
!   END IF
! 
!   500  CALL MPI_FINALIZE(IERR)
!   STOP
! 
!   !       write(*,*) mumps_par%MYID,"... is done!"
!   !       pause

END PROGRAM MUMPS_TEST
