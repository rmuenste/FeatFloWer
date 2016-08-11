C
      SUBROUTINE JacobiSolver(VA,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA)
C
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      DO 1 ITE=1,NIT
C
      CALL LCL1(DD,NEQ)
      DO 2 IEQ=1,NEQ
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     DD(IEQ)=DD(IEQ)-DBLE(VA(ICOL))*DX(KCOL(ICOL))
2     CONTINUE
C
      CALL E011Sum(DD)
C
      DO 4 IEQ=1,NEQ
4     DD(IEQ)=DD(IEQ)+DB(IEQ)
      DO 5 IEQ=1,NEQ
5     DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *         OMEGA*DD(IEQ)/DBLE(E011_UE(IEQ))
C
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE SSORSolver(VA,KCOL,KLD,DX,DB,DD,KNPR,NEQ,
     *          NIT,OMEGA)
C
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*),KNPR(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      DO 1 ITE=1,NIT
C
!     The real SOR part is here
!     ----------------------------------------------
      DO 2 IEQ=1,NEQ
       IF (KNPR(IEQ).NE.0) GOTO 2
       AUX=DB(IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        AUX=AUX-DBLE(VA(ICOL))*DX(KCOL(ICOL))
       END DO
       AUX=OMEGA*(AUX/DBLE(VA(KLD(IEQ)))-DX(IEQ))+DX(IEQ)
       DX(IEQ)=AUX
2     CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 3 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
       DD(IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        DD(IEQ)=DD(IEQ)-DBLE(VA(ICOL))*DX(KCOL(ICOL))
       END DO
3     CONTINUE
C
      CALL E011Sum(DD)
C
      DO 4 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *          OMEGA*(DD(IEQ)+DB(IEQ))/DBLE(E011_UE(IEQ))
4     CONTINUE
!     ----------------------------------------------
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO 22 IEQ=NEQ-1,1,-1
       IF (KNPR(IEQ).NE.0) GOTO 22
       AUX=DB(IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        AUX=AUX-DBLE(VA(ICOL))*DX(KCOL(ICOL))
       END DO
       AUX=OMEGA*(AUX/DBLE(VA(KLD(IEQ)))-DX(IEQ))+DX(IEQ)
       DX(IEQ)=AUX
22    CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 33 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 33
       DD(IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        DD(IEQ)=DD(IEQ)-DBLE(VA(ICOL))*DX(KCOL(ICOL))
       END DO
33    CONTINUE
C
      CALL E011Sum(DD)
C
      DO 44 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 44
       DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *          OMEGA*(DD(IEQ)+DB(IEQ))/DBLE(E011_UE(IEQ))
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE GetDefNorm(VA,KCOL,KLD,DX,DB,DD,NEQ,DEF)
C
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (myid.ne.0) THEN
       DO IEQ=1,NEQ
        DD(IEQ) = 0d0
        DO ICOL=KLD(IEQ),KLD(IEQ+1)-1
         DD(IEQ)=DD(IEQ)-DBLE(VA(ICOL))*DX(KCOL(ICOL))
        END DO
       END DO
C
       CALL E011Sum(DD)
C
       DO IEQ=1,NEQ
        DD(IEQ)=DD(IEQ)+DB(IEQ)
       END DO
C
       CALL LL21(DD,NEQ,DEF)
      END IF

      CALL COMM_Maximum(DEF)
C
      END
