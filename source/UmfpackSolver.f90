MODULE UMFPackSolver
!-------------------------------------------------------------------------------------------------
! Two interfaces to UMFPACK (umf4* wrapper, extern/libraries/umfpack4/modsrc):
!
! 1. myUmfPack_Factorize / _Solve / _Free (legacy, one global factorisation). The arrays are
!    the row-compressed matrix of the solver, handed to UMFPACK as column-compressed: UMFPACK
!    therefore factorises and solves with the TRANSPOSE of the stored matrix (equal for the
!    symmetric matrices it is used for; kept as it is, the certified runs depend on it). No
!    status is checked. Each handle is tracked on its own (legacy_symbolic_live,
!    legacy_numeric_live: the wrapper returns -1 for a failed step); numeric factorisation is
!    skipped after a failed symbolic analysis, and _Free releases only handles that exist.
!
! 2. tUmfCSC (explicit orientation, used by the periodic coarse operator, QuadSc_pco_solver):
!    the CSR matrix A of the caller is converted once into the compressed-COLUMN storage of A
!    (the own convention of UMFPACK), so umfpack solves A x = b (sys = UMFPACK_A = 0) and not the
!    transpose. The structure and the symbolic analysis are kept while the pattern is fixed;
!    refreshing the values redoes only the numeric factorisation. Every UMFPACK call returns
!    its status (Info(UMFPACK_STATUS+1)): 0 ok, 1 singular (warning), < 0 error; the caller
!    decides how to stop. Each tUmfCSC owns its own symbolic/numeric handles.
!-------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER numeric, symbolic,n,nz,sys
REAL*8 control (20), info (90)
LOGICAL :: legacy_symbolic_live = .FALSE., legacy_numeric_live = .FALSE.

! UMFPACK constants (umfpack.h)
INTEGER, PARAMETER :: UMF_SYS_A = 0          ! UMFPACK_A: solve A x = b
INTEGER, PARAMETER :: UMF_STATUS_OK = 0
INTEGER, PARAMETER :: UMF_WARNING_SINGULAR = 1

TYPE tUmfCSC
 INTEGER :: n = 0, nz = 0
 INTEGER, ALLOCATABLE :: ap(:), ai(:)     ! 0-based column pointers / row indices (ascending)
 REAL*8 , ALLOCATABLE :: ax(:)
 INTEGER, ALLOCATABLE :: from_csr(:)      ! CSC slot k holds CSR entry from_csr(k)
 INTEGER :: symbolic = -1, numeric = -1   ! wrapper handles, -1 = none
 REAL*8  :: control(20), info(90)
END TYPE tUmfCSC

PUBLIC :: tUmfCSC, UmfCSC_SetStructure, UmfCSC_SetValues, UmfCSC_Analyse, UmfCSC_Factorize, &
          UmfCSC_Solve, UmfCSC_Release

CONTAINS

SUBROUTINE myUmfPack_Solve(x,b,Ax,Mat,ItRef)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8 x(*),b(*),Ax(*)
 TYPE(TMatrix) Mat
 INTEGER ItRef

 n  = Mat%nu
 nz = Mat%na

 ! solving by a certain solution technique
 IF     (ItRef.EQ.1) THEN
  sys = 0
  call umf4solr (sys, Mat%LdA, Mat%ColA, Ax, x, b, numeric, control, info)
 ELSEIF (ItRef.EQ.2) THEN
  sys = 0
  call umf4sol (sys, x, b, numeric, control, info)
 ELSE
  WRITE(*,*) "no available solver"
 END IF

END SUBROUTINE myUmfPack_Solve
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_Factorize(Ax,Mat)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8 Ax(*)
 TYPE(TMatrix) Mat
! LOCAL
 INTEGER i,j,j1,j2

 n  = Mat%nu
 nz = Mat%na

 ! Reorganize the matrix into ascending format
 DO i=1,n
  j1 = Mat%LdA(i)
  j2 = Mat%LdA(i+1)-1
  CALL SortIt(j2-j1+1,Mat%ColA(j1),Ax(j1))
 END DO

 ! Reorganize the matrix to a 0-based one
 DO i = 1, n+1
   Mat%LdA(i) = Mat%LdA (i) - 1
 END DO
 DO i = 1, nz
   Mat%ColA (i) = Mat%ColA(i) - 1
 END DO

 ! Set up the default solver parameters
 call umf4def (control)

 ! Symbolic factorization
 call umf4sym (n, n, Mat%LdA, Mat%ColA, Ax, symbolic, control, info)
 legacy_symbolic_live = (symbolic.GE.0)

 ! Numeric factorization (not on a failed symbolic analysis: its handle is -1)
 IF (.NOT.legacy_symbolic_live) THEN
  numeric = -1
  legacy_numeric_live = .FALSE.
  RETURN
 END IF
 call umf4num (Mat%LdA, Mat%ColA, Ax, symbolic, numeric, control, info)
 legacy_numeric_live = (numeric.GE.0)

END SUBROUTINE myUmfPack_Factorize
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_Free()

 ! free the symbolic analysis
 IF (legacy_symbolic_live) call umf4fsym (symbolic)
 legacy_symbolic_live = .FALSE.

! free the numeric factorization
 IF (legacy_numeric_live) call umf4fnum (numeric)
 legacy_numeric_live = .FALSE.

END SUBROUTINE myUmfPack_Free
!
! -----------------------------------------------------
! Explicitly oriented interface (tUmfCSC), see the module header.
! -----------------------------------------------------
!
! Structure of A (n x n, CSR with 1-based ldA/colA, duplicates not allowed) as CSC of A.
! Releases a previous factorisation of F. Within a column the row indices ascend because
! the rows are visited in ascending order.
SUBROUTINE UmfCSC_SetStructure(F,nrow,ldA,colA)
IMPLICIT NONE
 TYPE(tUmfCSC), INTENT(INOUT) :: F
 INTEGER, INTENT(IN) :: nrow, ldA(:), colA(:)
 INTEGER :: i, k, j, nzz
 INTEGER, ALLOCATABLE :: fill(:)

 CALL UmfCSC_Release(F)
 nzz = ldA(nrow+1) - 1
 F%n  = nrow
 F%nz = nzz
 IF (ALLOCATED(F%ap)) DEALLOCATE(F%ap)
 IF (ALLOCATED(F%ai)) DEALLOCATE(F%ai)
 IF (ALLOCATED(F%ax)) DEALLOCATE(F%ax)
 IF (ALLOCATED(F%from_csr)) DEALLOCATE(F%from_csr)
 ALLOCATE(F%ap(nrow+1), F%ai(MAX(nzz,1)), F%ax(MAX(nzz,1)), F%from_csr(MAX(nzz,1)))
 ALLOCATE(fill(nrow+1))
 fill = 0
 DO k = 1, nzz
  j = colA(k)
  fill(j+1) = fill(j+1) + 1
 END DO
 F%ap(1) = 0
 DO j = 1, nrow
  F%ap(j+1) = F%ap(j) + fill(j+1)
 END DO
 fill(1:nrow) = F%ap(1:nrow)          ! next free 0-based slot per column
 DO i = 1, nrow
  DO k = ldA(i), ldA(i+1)-1
   j = colA(k)
   F%ai(fill(j)+1) = i - 1
   F%from_csr(fill(j)+1) = k
   fill(j) = fill(j) + 1
  END DO
 END DO
 DEALLOCATE(fill)
 F%ax = 0d0

END SUBROUTINE UmfCSC_SetStructure
!
! Values of A in the CSR order of UmfCSC_SetStructure.
SUBROUTINE UmfCSC_SetValues(F,a)
IMPLICIT NONE
 TYPE(tUmfCSC), INTENT(INOUT) :: F
 REAL*8, INTENT(IN) :: a(:)
 INTEGER :: k
 DO k = 1, F%nz
  F%ax(k) = a(F%from_csr(k))
 END DO
END SUBROUTINE UmfCSC_SetValues
!
! Symbolic analysis of the current structure (and values). status as returned by UMFPACK.
SUBROUTINE UmfCSC_Analyse(F,status)
IMPLICIT NONE
 TYPE(tUmfCSC), INTENT(INOUT) :: F
 INTEGER, INTENT(OUT) :: status
 IF (F%numeric.GE.0) CALL umf4fnum(F%numeric)
 F%numeric = -1
 IF (F%symbolic.GE.0) CALL umf4fsym(F%symbolic)
 F%symbolic = -1
 CALL umf4def(F%control)
 CALL umf4sym(F%n, F%n, F%ap, F%ai, F%ax, F%symbolic, F%control, F%info)
 status = NINT(F%info(1))
 IF (F%symbolic.LT.0 .AND. status.GE.0) status = -999
END SUBROUTINE UmfCSC_Analyse
!
! Numeric factorisation of the current values (symbolic analysis must exist).
SUBROUTINE UmfCSC_Factorize(F,status)
IMPLICIT NONE
 TYPE(tUmfCSC), INTENT(INOUT) :: F
 INTEGER, INTENT(OUT) :: status
 IF (F%symbolic.LT.0) THEN
  status = -998
  RETURN
 END IF
 IF (F%numeric.GE.0) CALL umf4fnum(F%numeric)
 F%numeric = -1
 CALL umf4num(F%ap, F%ai, F%ax, F%symbolic, F%numeric, F%control, F%info)
 status = NINT(F%info(1))
 IF (F%numeric.LT.0 .AND. status.GE.0) status = -997
END SUBROUTINE UmfCSC_Factorize
!
! Solve A x = b with iterative refinement (the stored A, not its transpose).
SUBROUTINE UmfCSC_Solve(F,x,b,status)
IMPLICIT NONE
 TYPE(tUmfCSC), INTENT(INOUT) :: F
 REAL*8, INTENT(INOUT) :: x(*)
 REAL*8, INTENT(IN) :: b(*)
 INTEGER, INTENT(OUT) :: status
 INTEGER :: isys
 IF (F%numeric.LT.0) THEN
  status = -996
  RETURN
 END IF
 isys = UMF_SYS_A
 CALL umf4solr(isys, F%ap, F%ai, F%ax, x, b, F%numeric, F%control, F%info)
 status = NINT(F%info(1))
END SUBROUTINE UmfCSC_Solve
!
SUBROUTINE UmfCSC_Release(F)
IMPLICIT NONE
 TYPE(tUmfCSC), INTENT(INOUT) :: F
 IF (F%numeric.GE.0) CALL umf4fnum(F%numeric)
 IF (F%symbolic.GE.0) CALL umf4fsym(F%symbolic)
 F%numeric = -1
 F%symbolic = -1
END SUBROUTINE UmfCSC_Release
!
! -----------------------------------------------------
!
SUBROUTINE SortIt(n,K,A)
IMPLICIT NONE
REAL*8 A(*)
INTEGER n,K(*)
INTEGER i,j,iWA
REAL*8 dWA

! WRITE(*,*) n
! WRITE(*,*) (K(i),i=1,n)
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
! WRITE(*,*) (K(i),i=1,n)

END SUBROUTINE SortIt

END MODULE UMFPackSolver