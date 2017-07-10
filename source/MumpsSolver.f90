MODULE MumpsSolver
USE PP3D_MPI
USE var_QuadScalar

IMPLICIT NONE
INCLUDE 'dmumps_struc.h'
TYPE (DMUMPS_STRUC) mumps_par

 CONTAINS
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_Init()
IMPLICIT NONE

mumps_par%COMM = MPI_COMM_WORLD
mumps_par%MYID = myid

mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 0

MUMPS_PAR%icntl(1)  = 0
MUMPS_PAR%icntl(2)  = 0
MUMPS_PAR%icntl(3)  = 0

CALL DMUMPS(mumps_par)

IF (mumps_par%INFOG(1).LT.0) THEN
 WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
 STOP
END IF

END SUBROUTINE MUMPS_Init
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_SetUpP1_MASTER(L,B)
IMPLICIT NONE
REAL*8 B(*)
TYPE(tMatrix) L
INTEGER I,J,K

 mumps_par%N  = L%nu
 mumps_par%NZ = L%na
 ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
 ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
 K = 0
 DO I = 1, L%nu
  DO J = L%LdA(i),L%LdA(i+1)-1
   K = K + 1
   mumps_par%IRN(k) = I
   mumps_par%JCN(k) = L%ColA(J) 
  END DO
 END DO
 ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
 DO I = 1, mumps_par%N
  mumps_par%RHS(I) = B(I)
 END DO

END SUBROUTINE MUMPS_SetUpP1_MASTER
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_SetUpQ2_1_MASTER(L,B)
IMPLICIT NONE
INTEGER NJ
REAL*8 B(*)
TYPE(tMatrix) L
INTEGER I,J,K

 NJ = L%nu
 mumps_par%N  = 3*L%nu
 mumps_par%NZ = 3*L%na
 ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
 ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
 K = 0
 DO I = 1, NJ
  DO J = L%LdA(i),L%LdA(i+1)-1
   K = K + 1
   mumps_par%IRN(k) = 0*NJ + I
   mumps_par%JCN(k) = 0*NJ + L%ColA(J) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   K = K + 1
   mumps_par%IRN(k) = 1*NJ + I
   mumps_par%JCN(k) = 1*NJ + L%ColA(J) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   K = K + 1
   mumps_par%IRN(k) = 2*NJ + I
   mumps_par%JCN(k) = 2*NJ + L%ColA(J) 
  END DO
 END DO
 
 ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
 DO I = 1, mumps_par%N
  mumps_par%RHS(I) = B(I)
 END DO
END SUBROUTINE MUMPS_SetUpQ2_1_MASTER
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_SetUpQ2_3_MASTER(L,B)
IMPLICIT NONE
REAL*8 B(*)
TYPE(tMatrix) L
INTEGER I,J,K
INTEGER NJ

 NJ = L%nu
 mumps_par%N  = 3*L%nu
 mumps_par%NZ = 9*L%na
 ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
 ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
 K = 0
 DO I = 1, NJ
  DO J = L%LdA(i),L%LdA(i+1)-1
   K = K + 1
   mumps_par%IRN(k) = 0*NJ + I
   mumps_par%JCN(k) = 0*NJ + L%ColA(J) 

   K = K + 1
   mumps_par%IRN(k) = 0*NJ + I
   mumps_par%JCN(k) = 1*NJ + L%ColA(J) 

   K = K + 1
   mumps_par%IRN(k) = 0*NJ + I
   mumps_par%JCN(k) = 2*NJ + L%ColA(J) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   K = K + 1
   mumps_par%IRN(k) = 1*NJ + I
   mumps_par%JCN(k) = 0*NJ + L%ColA(J) 

   K = K + 1
   mumps_par%IRN(k) = 1*NJ + I
   mumps_par%JCN(k) = 1*NJ + L%ColA(J) 

   K = K + 1
   mumps_par%IRN(k) = 1*NJ + I
   mumps_par%JCN(k) = 2*NJ + L%ColA(J) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   K = K + 1
   mumps_par%IRN(k) = 2*NJ + I
   mumps_par%JCN(k) = 0*NJ + L%ColA(J) 

   K = K + 1
   mumps_par%IRN(k) = 2*NJ + I
   mumps_par%JCN(k) = 1*NJ + L%ColA(J) 

   K = K + 1
   mumps_par%IRN(k) = 2*NJ + I
   mumps_par%JCN(k) = 2*NJ + L%ColA(J) 
  END DO
 END DO
 
 ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
 DO I = 1, mumps_par%N
  mumps_par%RHS(I) = B(I)
 END DO

END SUBROUTINE MUMPS_SetUpQ2_3_MASTER
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_SetUpP1_SLAVE(A,AP,L,LP,X,XP,N)
IMPLICIT NONE
INTEGER N
REAL*8 A(*),AP(*),X(*),XP(*)
TYPE(tMatrix) L,LP
INTEGER I,J,K,iG,jG

 DO i=1,N
  X(i) = DBLE(GlobalNumberingP1(i))
 END DO
 
 CALL GetParPressure(X,XP)

 mumps_par%NZ_loc = L%na + LP%na
 ALLOCATE( mumps_par%IRN_loc ( mumps_par%NZ_loc) )
 ALLOCATE( mumps_par%JCN_loc ( mumps_par%NZ_loc) )
 ALLOCATE( mumps_par%A_loc   ( mumps_par%NZ_loc) )
 K = 0
 DO I = 1, L%nu
  iG = GlobalNumberingP1(I) 

  DO J = L%LdA(i),L%LdA(i+1)-1
   jG = GlobalNumberingP1(L%ColA(J)) 
   K = K + 1
   mumps_par%IRN_loc(k) = iG
   mumps_par%JCN_loc(k) = jG
   mumps_par%A_loc(k)   = A(J)
  END DO

  DO J = LP%LdA(i),LP%LdA(i+1)-1
   jG = INT(XP(LP%ColA(J)))
   K = K + 1
   mumps_par%IRN_loc(k) = iG
   mumps_par%JCN_loc(k) = jG
   mumps_par%A_loc(k)   = AP(J)
  END DO
 END DO
     
END SUBROUTINE MUMPS_SetUpP1_SLAVE
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_SetUpQ2_1_SLAVE(A11,a22,a33,L,X,N)
IMPLICIT NONE
INTEGER N
REAL*8 A11(*),A22(*),A33(*),X(*)
TYPE(tMatrix) L
INTEGER I,J,K,iG,jG

 mumps_par%NZ_loc = 3*L%na
 ALLOCATE( mumps_par%IRN_loc ( mumps_par%NZ_loc) )
 ALLOCATE( mumps_par%JCN_loc ( mumps_par%NZ_loc) )
 ALLOCATE( mumps_par%A_loc   ( mumps_par%NZ_loc) )

 K = 0

 DO I = 1, L%nu
  iG = 0*myGlobal_ndof + GlobalNumberingQ2(I) 
  DO J = L%LdA(i),L%LdA(i+1)-1
   jG = 0*myGlobal_ndof + GlobalNumberingQ2(L%ColA(J)) 
   K = K + 1
   mumps_par%IRN_loc(k) = iG
   mumps_par%JCN_loc(k) = jG
   mumps_par%A_loc(k)   = A11(J)
  END DO
 END DO
     
 DO I = 1, L%nu
  iG = 1*myGlobal_ndof + GlobalNumberingQ2(I) 
  DO J = L%LdA(i),L%LdA(i+1)-1
   jG = 1*myGlobal_ndof + GlobalNumberingQ2(L%ColA(J)) 
   K = K + 1
   mumps_par%IRN_loc(k) = iG
   mumps_par%JCN_loc(k) = jG
   mumps_par%A_loc(k)   = A22(J)
  END DO
 END DO

 DO I = 1, L%nu
  iG = 2*myGlobal_ndof + GlobalNumberingQ2(I) 
  DO J = L%LdA(i),L%LdA(i+1)-1
   jG = 2*myGlobal_ndof + GlobalNumberingQ2(L%ColA(J)) 
   K = K + 1
   mumps_par%IRN_loc(k) = iG
   mumps_par%JCN_loc(k) = jG
   mumps_par%A_loc(k)   = A33(J)
  END DO
 END DO
 
END SUBROUTINE MUMPS_SetUpQ2_1_SLAVE
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_SetUpQ2_3_SLAVE(A11,A12,A13,A21,A22,A23,A31,A32,A33,L,X,N)
IMPLICIT NONE
INTEGER N
REAL*8 A11(*),A22(*),A33(*),X(*)
REAL*8 A12(*),A13(*),A21(*),A23(*),A31(*),A32(*)
TYPE(tMatrix) L
INTEGER I,J,K,iG,jG
INTEGER NJ,NI

 NJ = L%nu
 NI = myGlobal_ndof

 mumps_par%NZ_loc = 9*L%na
 ALLOCATE( mumps_par%IRN_loc ( mumps_par%NZ_loc) )
 ALLOCATE( mumps_par%JCN_loc ( mumps_par%NZ_loc) )
 ALLOCATE( mumps_par%A_loc   ( mumps_par%NZ_loc) )

 K = 0

 DO I = 1, NJ
  iG = GlobalNumberingQ2(I) 
  DO J = L%LdA(i),L%LdA(i+1)-1
   jG = GlobalNumberingQ2(L%ColA(J)) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   K = K + 1
   mumps_par%IRN_loc(k) = 0*NI + iG
   mumps_par%JCN_loc(k) = 0*NI + jG
   mumps_par%A_loc(k)   = A11(J)

   K = K + 1
   mumps_par%IRN_loc(k) = 0*NI + iG
   mumps_par%JCN_loc(k) = 1*NI + jG
   mumps_par%A_loc(k)   = A12(J)
   
   K = K + 1
   mumps_par%IRN_loc(k) = 0*NI + iG
   mumps_par%JCN_loc(k) = 2*NI + jG
   mumps_par%A_loc(k)   = A13(J)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   K = K + 1
   mumps_par%IRN_loc(k) = 1*NI + iG
   mumps_par%JCN_loc(k) = 0*NI + jG
   mumps_par%A_loc(k)   = A21(J)

   K = K + 1
   mumps_par%IRN_loc(k) = 1*NI + iG
   mumps_par%JCN_loc(k) = 1*NI + jG
   mumps_par%A_loc(k)   = A22(J)
   
   K = K + 1
   mumps_par%IRN_loc(k) = 1*NI + iG
   mumps_par%JCN_loc(k) = 2*NI + jG
   mumps_par%A_loc(k)   = A23(J)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   K = K + 1
   mumps_par%IRN_loc(k) = 2*NI + iG
   mumps_par%JCN_loc(k) = 0*NI + jG
   mumps_par%A_loc(k)   = A31(J)

   K = K + 1
   mumps_par%IRN_loc(k) = 2*NI + iG
   mumps_par%JCN_loc(k) = 1*NI + jG
   mumps_par%A_loc(k)   = A32(J)
   
   K = K + 1
   mumps_par%IRN_loc(k) = 2*NI + iG
   mumps_par%JCN_loc(k) = 2*NI + jG
   mumps_par%A_loc(k)   = A33(J)
  END DO
 END DO
    
END SUBROUTINE MUMPS_SetUpQ2_3_SLAVE
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_Solve(X)
IMPLICIT NONE
REAL*8 X(*)
INTEGER I

MUMPS_PAR%icntl(1)  = 0 
MUMPS_PAR%icntl(2)  = 0
MUMPS_PAR%icntl(3)  = 0
!  Call package for solution
MUMPS_PAR%icntl(6)  = 7
!     pivot order (automatic)
MUMPS_PAR%icntl(7)  = 7
!     scaling (automatic)
MUMPS_PAR%icntl(8)  = 77
!     no transpose
MUMPS_PAR%icntl(9)  = 1
!     max steps for iterative refinement
MUMPS_PAR%icntl(10) = 0
!     statistics info
MUMPS_PAR%icntl(11) = 0
!     controls parallelism
MUMPS_PAR%icntl(12) = 0
!     use ScaLAPACK for root node
MUMPS_PAR%icntl(13) = 0
!     percentage increase in estimated workspace
MUMPS_PAR%icntl(14) = 100

! MUMPS_PAR%icntl(4)  = 0
! mumps_par%ICNTL(5)  = 0
mumps_par%ICNTL(18) = 3

mumps_par%JOB = 6
CALL DMUMPS(mumps_par)

!       mumps_par%JOB = 2
!       CALL DMUMPS(mumps_par)
! 
!       mumps_par%JOB = 3
!       CALL DMUMPS(mumps_par)

IF (mumps_par%INFOG(1).LT.0) THEN
 WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
 STOP
END IF


!  Solution has been assembled on the host
IF ( MYID .eq. 0 ) THEN
!   WRITE( 6, * ) ' Solution is '
 DO I=1,mumps_par%N
!   WRITE( 6, * ) mumps_par%RHS(I)
  X(i) = mumps_par%RHS(I)
 END DO
END IF

! pause
END SUBROUTINE MUMPS_Solve
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_CleanUp
IMPLICIT NONE

!  Deallocate user data
IF ( mumps_par%MYID .eq. 0 )THEN
  DEALLOCATE( mumps_par%RHS )
END IF

mumps_par%JOB = -2
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1).LT.0) THEN
 WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
 STOP
END IF

END SUBROUTINE MUMPS_CleanUp
!
! ------------------------------------------------------------------------
!
SUBROUTINE MUMPS_solver_Distributed()
USE var_QuadScalar, ONLY:myCrsMat
USE PP3D_MPI
IMPLICIT NONE
! INCLUDE 'mpif.h'
INCLUDE 'dmumps_struc.h'
TYPE (DMUMPS_STRUC) mumps_par
INTEGER I
CHARACTER cFile*(20)

!CALL MPI_INIT(IERR)

! Define a communicator for the package.
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%MYID = myid
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)

mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 0

MUMPS_PAR%icntl(1)  = 0
MUMPS_PAR%icntl(2)  = 0
MUMPS_PAR%icntl(3)  = 0

CALL DMUMPS(mumps_par)

IF (mumps_par%INFOG(1).LT.0) THEN
 WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
 GOTO 500
END IF

! WRITE(cFile,'(A,I2.2,A)')  'input_',mumps_par%MYID,'.txt'
! OPEN(FILE=ADJUSTL(TRIM(cFile)),UNIT=55)

!  Define problem on the host (processor 0)
IF ( mumps_par%MYID .eq. 0 ) THEN

  mumps_par%N  = myCrsMat%nu
  mumps_par%NZ = myCrsMat%na
  ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
  ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
!         ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )
!   mumps_par%A  = 0d0
  DO I = 1, mumps_par%NZ
    mumps_par%IRN(I) = myCrsMat%Row(i)
    mumps_par%JCN(I) = myCrsMat%Col(i)
  END DO
  ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
  DO I = 1, mumps_par%N
   mumps_par%RHS(I) = myCrsMat%D(I)
  END DO
ELSE
  
  mumps_par%NZ_loc = myCrsMat%na
  ALLOCATE( mumps_par%IRN_loc ( mumps_par%NZ_loc) )
  ALLOCATE( mumps_par%JCN_loc ( mumps_par%NZ_loc) )
  ALLOCATE( mumps_par%A_loc( mumps_par%NZ_loc) )
  DO I = 1, mumps_par%NZ_loc
    mumps_par%IRN_loc(I) = myCrsMat%Row(i)
    mumps_par%JCN_loc(I) = myCrsMat%Col(i)
    mumps_par%A_loc(I)   = myCrsMat%A(i)
  END DO
END IF

! CLOSE(55)  
     
MUMPS_PAR%icntl(1)  = 0 
MUMPS_PAR%icntl(2)  = 0
MUMPS_PAR%icntl(3)  = 0
!  Call package for solution
MUMPS_PAR%icntl(6)  = 7
!     pivot order (automatic)
MUMPS_PAR%icntl(7)  = 7
!     scaling (automatic)
MUMPS_PAR%icntl(8)  = 77
!     no transpose
MUMPS_PAR%icntl(9)  = 1
!     max steps for iterative refinement
MUMPS_PAR%icntl(10) = 0
!     statistics info
MUMPS_PAR%icntl(11) = 0
!     controls parallelism
MUMPS_PAR%icntl(12) = 0
!     use ScaLAPACK for root node
MUMPS_PAR%icntl(13) = 0
!     percentage increase in estimated workspace
MUMPS_PAR%icntl(14) = 100

! MUMPS_PAR%icntl(4)  = 0
! mumps_par%ICNTL(5)  = 0
mumps_par%ICNTL(18) = 3

mumps_par%JOB = 6
CALL DMUMPS(mumps_par)

!       mumps_par%JOB = 2
!       CALL DMUMPS(mumps_par)
! 
!       mumps_par%JOB = 3
!       CALL DMUMPS(mumps_par)

IF (mumps_par%INFOG(1).LT.0) THEN
 WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
 GOTO 500
END IF
!  Solution has been assembled on the hostCALL MPI_FINALIZE(IERR)
IF ( mumps_par%MYID .eq. 0 ) THEN
!   WRITE( 6, * ) ' Solution is '
 DO I=1,mumps_par%N
!   WRITE( 6, * ) mumps_par%RHS(I)
  myCrsMat%D(I) = mumps_par%RHS(I)
 END DO
END IF

! write(*,*) myid,' is done!' 
! pause

!  Deallocate user data
IF ( mumps_par%MYID .eq. 0 )THEN
!         DEALLOCATE( mumps_par%A )
  DEALLOCATE( mumps_par%RHS )
END IF
!  Destroy the instance (deallocate internal data structures)

mumps_par%JOB = -2
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1).LT.0) THEN
 WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
 GOTO 500
END IF

500  CONTINUE
!500  CALL MPI_FINALIZE(IERR)

END SUBROUTINE MUMPS_solver_Distributed
END MODULE MumpsSolver
