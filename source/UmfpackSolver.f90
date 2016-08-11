MODULE UMFPackSolver
IMPLICIT NONE

INTEGER numeric, symbolic,n,nz,sys
REAL*8 control (20), info (90)

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

 ! Numeric factorization
 call umf4num (Mat%LdA, Mat%ColA, Ax, symbolic, numeric, control, info)

END SUBROUTINE myUmfPack_Factorize
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_Free()

! free the symbolic analysis
 call umf4fsym (symbolic)

! free the numeric factorization
 call umf4fnum (numeric)

END SUBROUTINE myUmfPack_Free
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