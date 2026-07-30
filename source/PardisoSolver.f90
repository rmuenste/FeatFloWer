MODULE PardisoSolver
! MKL PARDISO direct solver for the gathered coarse pressure system
! (Pres@MGCrsSolverType = 9). Mirrors the UMFPackSolver
! factorize/solve/free split: the factorization is built once in
! Setup_UMFPACK_CoarseSolver and reused across coarse solves until
! myPardiso_Free is called. Runs on the master rank only, on the same
! gathered UMF_CMat/UMF_lMat matrix as solver type 2.
! Fatal paths must call MPI_Abort, not STOP: a master-only stop would leave
! worker ranks blocked in MPI.
!
! Compiled only with USE_MKL_PARDISO=ON (links sequential MKL; the master
! rank is bound to one core during the coarse solve).
IMPLICIT NONE

INTEGER*8 :: pardiso_pt(64) = 0
INTEGER :: pardiso_iparm(64) = 0
INTEGER, PARAMETER :: pardiso_mtype = 11   ! real, nonsymmetric
INTEGER, PARAMETER :: pardiso_maxfct = 1, pardiso_mnum = 1
INTEGER :: pardiso_n = 0
LOGICAL :: bPardisoFactorized = .FALSE.

CONTAINS

SUBROUTINE myPardiso_Factorize(Ax,Mat)
USE var_QuadScalar
USE PP3D_MPI, ONLY: MPI_COMM_WORLD
IMPLICIT NONE

 REAL*8 Ax(*)
 TYPE(TMatrix) Mat
 INTEGER i,j1,j2,ierror,idum(1),abortError
 REAL*8 ddum(1)

 ! PARDISO requires ascending column indices within each row
 DO i=1,Mat%nu
  j1 = Mat%LdA(i)
  j2 = Mat%LdA(i+1)-1
  CALL PardisoSortRow(j2-j1+1,Mat%ColA(j1),Ax(j1))
 END DO

 IF (bPardisoFactorized) CALL myPardiso_Free()

 pardiso_pt = 0
 pardiso_iparm = 0
 pardiso_iparm(1)  = 1    ! user-supplied iparm
 pardiso_iparm(2)  = 2    ! METIS nested-dissection ordering
 pardiso_iparm(8)  = 2    ! max iterative refinement steps
 pardiso_iparm(10) = 13   ! pivot perturbation 1e-13
 pardiso_iparm(11) = 1    ! scaling
 pardiso_iparm(13) = 1    ! weighted matching
 pardiso_iparm(35) = 0    ! 1-based (Fortran) indexing
 pardiso_n = Mat%nu

 ! Analysis + numerical factorization (phase 12), 1-based CSR
 CALL pardiso(pardiso_pt,pardiso_maxfct,pardiso_mnum,pardiso_mtype,12,&
      pardiso_n,Ax,Mat%LdA,Mat%ColA,idum,1,pardiso_iparm,0,ddum,ddum,ierror)

 IF (ierror.NE.0) THEN
  WRITE(*,*) 'PARDISO factorization failed with error ',ierror
  FLUSH(6)
  CALL MPI_Abort(MPI_COMM_WORLD,myErrorCode%SOLVER_RUNTIME_FAILURE,abortError)
 END IF

 bPardisoFactorized = .TRUE.

END SUBROUTINE myPardiso_Factorize
!
! -----------------------------------------------------
!
SUBROUTINE myPardiso_Solve(x,b,Ax,Mat)
USE var_QuadScalar
USE PP3D_MPI, ONLY: MPI_COMM_WORLD
IMPLICIT NONE

 REAL*8 x(*),b(*),Ax(*)
 TYPE(TMatrix) Mat
 INTEGER ierror,idum(1),abortError

 IF (.NOT.bPardisoFactorized) THEN
  WRITE(*,*) 'PARDISO solve called without factorization'
  FLUSH(6)
  CALL MPI_Abort(MPI_COMM_WORLD,myErrorCode%SOLVER_RUNTIME_FAILURE,abortError)
 END IF

 ! Solve + iterative refinement (phase 33); b is unchanged, x receives
 ! the solution
 CALL pardiso(pardiso_pt,pardiso_maxfct,pardiso_mnum,pardiso_mtype,33,&
      pardiso_n,Ax,Mat%LdA,Mat%ColA,idum,1,pardiso_iparm,0,b,x,ierror)

 IF (ierror.NE.0) THEN
  WRITE(*,*) 'PARDISO solve failed with error ',ierror
  FLUSH(6)
  CALL MPI_Abort(MPI_COMM_WORLD,myErrorCode%SOLVER_RUNTIME_FAILURE,abortError)
 END IF

END SUBROUTINE myPardiso_Solve
!
! -----------------------------------------------------
!
SUBROUTINE myPardiso_Free()
IMPLICIT NONE

 INTEGER ierror,idum(1)
 REAL*8 ddum(1)

 IF (.NOT.bPardisoFactorized) RETURN

 ! Release all internal memory (phase -1)
 CALL pardiso(pardiso_pt,pardiso_maxfct,pardiso_mnum,pardiso_mtype,-1,&
      pardiso_n,ddum,idum,idum,idum,1,pardiso_iparm,0,ddum,ddum,ierror)

 bPardisoFactorized = .FALSE.

END SUBROUTINE myPardiso_Free
!
! -----------------------------------------------------
!
SUBROUTINE PardisoSortRow(n,K,A)
IMPLICIT NONE
REAL*8 A(*)
INTEGER n,K(*)
INTEGER i,j,iWA
REAL*8 dWA

DO I=2,N
 DO J=N,I,-1
  IF (K(J).LT.K(J-1)) THEN
   iWA     = K(J)
   dWA     = A(J)
   K(J)    = K(J-1)
   A(J)    = A(J-1)
   K(J-1)  = iWA
   A(J-1)  = dWA
  END IF
 END DO
END DO

END SUBROUTINE PardisoSortRow

END MODULE PardisoSolver
