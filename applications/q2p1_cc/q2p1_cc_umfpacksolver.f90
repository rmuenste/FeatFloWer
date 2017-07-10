MODULE UMFPackSolver_CC
IMPLICIT NONE

INTEGER numeric, symbolic,n,nz,sys
REAL*8 control (20), info (90)

CONTAINS

SUBROUTINE myUmfPack_crsCCReshape(Ax,Mat)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8 Ax(*)
 TYPE(TMatrix) Mat
 INTEGER sym,num
! LOCAL
 INTEGER i,j,j1,j2

 n  = Mat%nu
 nz = Mat%na

 ! Reorganize the matrix into ascending format
 DO i=1,n
  j1 = Mat%LdA(i)
  j2 = Mat%LdA(i+1)-1
  CALL SortIt_CC(j2-j1+1,Mat%ColA(j1),Ax(j1))
 END DO

 ! Reorganize the matrix to a 0-based one
 DO i = 1, n+1
   Mat%LdA(i) = Mat%LdA (i) - 1
 END DO
 DO i = 1, nz
   Mat%ColA (i) = Mat%ColA(i) - 1
 END DO

END SUBROUTINE myUmfPack_crsCCReshape
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_crsCCSymFactorize(Ax,Mat,sym)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8 Ax(*)
 TYPE(TMatrix) Mat
 INTEGER sym

 n  = Mat%nu
 nz = Mat%na

 ! Set up the default solver parameters
 call umf4def (control)

 ! Symbolic factorization
 call umf4sym (n, n, Mat%LdA, Mat%ColA, Ax, sym, control, info)

END SUBROUTINE myUmfPack_crsCCSymFactorize
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_crsCCNumFactorize(Ax,Mat,sym,num)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8 Ax(*)
 TYPE(TMatrix) Mat
 INTEGER sym,num
! LOCAL
 INTEGER i,j,j1,j2

 n  = Mat%nu
 nz = Mat%na

 ! Numeric factorization
 call umf4num (Mat%LdA, Mat%ColA, Ax, sym, num, control, info)

END SUBROUTINE myUmfPack_crsCCNumFactorize
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_CCSolveMaster(x,b,Ax,MatLdA,MatColA,sym,num,nEq)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8 Ax(*),x(*),b(*)
 INTEGER MatLdA(*),MatColA(*)
 INTEGER sym,num,nEq

  sys = 2
  call umf4solr (sys, MatLdA, MatColA, Ax, x, b, num, control, info)

END SUBROUTINE myUmfPack_CCSolveMaster
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_CCSolve(x,b,Ax,MatLdA,MatColA,sym,num,nEq)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8 Ax(*),x(*),b(*)
 INTEGER MatLdA(*),MatColA(*)
 INTEGER sym,num,nEq

  sys = 0
  call umf4solr (sys, MatLdA, MatColA, Ax, x, b, num, control, info)

END SUBROUTINE myUmfPack_CCSolve
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_CCFactorizeLocalMat(Ax,Mat,sym,num)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8  Ax(*)
 TYPE(TMatrix) Mat
 INTEGER sym,num
! LOCAL
 INTEGER i,j,j1,j2,n,nz

 n  = Mat%nu
 nz = Mat%na

 ! Reorganize the matrix into ascending format
 DO i=1,n
  j1 = Mat%LdA(i)
  j2 = Mat%LdA(i+1)-1
  CALL SortIt_CC(j2-j1+1,Mat%ColA(j1),Ax(j1))
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
 call umf4sym (n, n, Mat%LdA, Mat%ColA, Ax, sym, control, info)

 ! Numeric factorization
 call umf4num (Mat%LdA, Mat%ColA, Ax, sym, num, control, info)

END SUBROUTINE myUmfPack_CCFactorizeLocalMat
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_CCSolveLocalMat(x,b,Ax,Mat,sym,num,nEq)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8  Ax(*),x(*),b(*)
 TYPE(TMatrix) Mat
 INTEGER sym,num,nEq

  sys = 2
  call umf4solr (sys, Mat%LdA, Mat%ColA, Ax, x, b, num, control, info)

END SUBROUTINE myUmfPack_CCSolveLocalMat
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_CCFactorize(Ax,MatLdA,MatColA,sym,num,nEq)
USE var_QuadScalar
IMPLICIT NONE

 ! INPUT
 REAL*8  Ax(*)
 INTEGER MatLdA(*),MatColA(*)
 INTEGER sym,num,nEq
! LOCAL
 INTEGER i,j,j1,j2

 ! Set up the default solver parameters
 call umf4def (control)

 ! Symbolic factorization
 call umf4sym (nEq, nEq, MatLdA, MatColA, Ax, sym, control, info)

 ! Numeric factorization
 call umf4num (MatLdA, MatColA, Ax, sym, num, control, info)

END SUBROUTINE myUmfPack_CCFactorize
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_CCcrsFree(num)
 INTEGER num

! free the numeric factorization
 call umf4fnum (num)

END SUBROUTINE myUmfPack_CCcrsFree
!
! -----------------------------------------------------
!
SUBROUTINE myUmfPack_CCFree(sym, num)
 INTEGER sym, num

! free the symbolic analysis
 call umf4fsym (sym)

! free the numeric factorization
 call umf4fnum (num)

END SUBROUTINE myUmfPack_CCFree
!
! -----------------------------------------------------
!
SUBROUTINE SortIt_CC(n,K,A)
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

END SUBROUTINE SortIt_CC

END MODULE UMFPackSolver_CC
