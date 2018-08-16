      SUBROUTINE E011_SSORSolver(DA11,DA22,DA33,KCOL,KLD,DX,
     *           DB,DD,KNPR,NEQ,NIT,OMEGA,DEF)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE011,myid,E011SUM
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCOL(*),KLD(*),DX(*),DB(*),DD(*),KNPR(*)
      DIMENSION DA11(*),DA22(*),DA33(*)
      COMMON /ERRCTL/ IER,ICHECK
      INTEGER MEQ1,MEQ2,MEQ3
      SAVE /ERRCTL/
C
      MEQ1 = 0
      MEQ2 = NEQ
      MEQ3 = 2*NEQ
C
      DO 1 ITE=1,NIT
C
!     The real SOR part is here
!     ----------------------------------------------
      DO 2 IEQ=1,NEQ
       IF (KNPR(IEQ).NE.0) GOTO 2
       AUX1=DB(MEQ1+IEQ)
       AUX2=DB(MEQ2+IEQ)
       AUX3=DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J=KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
2     CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 3 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
       DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
       DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
       DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
       END DO
3     CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
C
      DO 4 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *          OMEGA*DD(MEQ1+IEQ)/MGE011(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *          OMEGA*DD(MEQ2+IEQ)/MGE011(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *          OMEGA*DD(MEQ3+IEQ)/MGE011(ILEV)%UE33(IEQ)
4     CONTINUE
!     ----------------------------------------------
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO 22 IEQ=NEQ-1,1,-1
       IF (KNPR(IEQ).NE.0) GOTO 22
       AUX1=DB(MEQ1+IEQ)
       AUX2=DB(MEQ2+IEQ)
       AUX3=DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
22    CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 33 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 33
       DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
       DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
       DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
       END DO
33    CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
C
      DO 44 IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *          OMEGA*DD(MEQ1+IEQ)/MGE011(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *          OMEGA*DD(MEQ2+IEQ)/MGE011(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *          OMEGA*DD(MEQ3+IEQ)/MGE011(ILEV)%UE33(IEQ)
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      DO 122 IEQ=1,NEQ
      DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
      DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
      DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
      DO 133 ICOL=KLD(IEQ),KLD(IEQ+1)-1
      J=KCOL(ICOL)
      DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
      DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
      DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
133    CONTINUE
122    CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
      CALL LL21(DD,3*NEQ,DEF)
C
!       write(*,*) myid,DEF
!       IF (NIT.NE.0) THEN
! !       IF (myid.eq.1) write(*,'(2I5,6D12.4)') 
! !      * 1,1,DX(108),Db(108),Dd(108),DX(271),Db(271),Dd(271)
! !       IF (myid.eq.2) write(*,'(2I5,6D12.4)') 
! !      * 2,1,DX(94),Db(94),Dd(94),DX(261),Db(261),Dd(261)
!        pause
!       END IF
      END
C
C
C
      SUBROUTINE E011_SSORSmoother(DA11,DA22,DA33,KCOL,KLD,DX,DB,
     *           DD,KNPR,NEQ,NIT,OMEGA)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE011,myid,E011SUM
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCOL(*),KLD(*),DX(*),DB(*),DD(*),KNPR(*)
      DIMENSION DA11(*),DA22(*),DA33(*)
      COMMON /ERRCTL/ IER,ICHECK
      INTEGER MEQ1,MEQ2,MEQ3
      SAVE /ERRCTL/
C
      MEQ1 = 0
      MEQ2 = NEQ
      MEQ3 = 2*NEQ
C
      DO 1 ITE=1,NIT
C
!     The real SOR part is here
!     ----------------------------------------------
      DO 2 IEQ=1,NEQ
       IF (KNPR(IEQ).NE.0) GOTO 2
       AUX1=DB(MEQ1+IEQ)
       AUX2=DB(MEQ2+IEQ)
       AUX3=DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J=KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
2     CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 3 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
       DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
       DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
       DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
       END DO
3     CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
C
      DO 4 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *          OMEGA*DD(MEQ1+IEQ)/MGE011(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *          OMEGA*DD(MEQ2+IEQ)/MGE011(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *          OMEGA*DD(MEQ3+IEQ)/MGE011(ILEV)%UE33(IEQ)
4     CONTINUE
!     ----------------------------------------------
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO 22 IEQ=NEQ-1,1,-1
       IF (KNPR(IEQ).NE.0) GOTO 22
       AUX1=DB(MEQ1+IEQ)
       AUX2=DB(MEQ2+IEQ)
       AUX3=DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
22    CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 33 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 33
       DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
       DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
       DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
       END DO
33    CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
C
      DO 44 IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *          OMEGA*DD(MEQ1+IEQ)/MGE011(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *          OMEGA*DD(MEQ2+IEQ)/MGE011(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *          OMEGA*DD(MEQ3+IEQ)/MGE011(ILEV)%UE33(IEQ)
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      END
C
C
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
      SUBROUTINE SSORSolverX(DA,KCOL,KLD,DX,DB,DD,KNPR,NEQ,
     *          NIT,OMEGA,iC)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE011,E011_UE,myid,E011SUM
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*),KNPR(*)
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
        AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
       END DO
       AUX=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
       DX(IEQ)=AUX
2     CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 3 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
       DD(IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
       END DO
3     CONTINUE
C
      CALL E011Sum(DD)
C
      DO 4 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       IF (iC.eq.1) DIAG=MGE011(ILEV)%UE11(IEQ)
       IF (iC.eq.2) DIAG=MGE011(ILEV)%UE22(IEQ)
       IF (iC.eq.3) DIAG=MGE011(ILEV)%UE33(IEQ)
!        DIAG = E011_UE(IEQ)
       DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *          OMEGA*(DD(IEQ)+DB(IEQ))/DIAG
4     CONTINUE
!     ----------------------------------------------
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO 22 IEQ=NEQ-1,1,-1
       IF (KNPR(IEQ).NE.0) GOTO 22
       AUX=DB(IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
       END DO
       AUX=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
       DX(IEQ)=AUX
22    CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 33 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 33
       DD(IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
       END DO
33    CONTINUE
C
      CALL E011Sum(DD)
C
      DO 44 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 44
       IF (iC.eq.1) DIAG=MGE011(ILEV)%UE11(IEQ)
       IF (iC.eq.2) DIAG=MGE011(ILEV)%UE22(IEQ)
       IF (iC.eq.3) DIAG=MGE011(ILEV)%UE33(IEQ)
!        DIAG = E011_UE(IEQ)
       DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *          OMEGA*(DD(IEQ)+DB(IEQ))/DIAG
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE SSORSolver(DA,KCOL,KLD,DX,DB,DD,KNPR,NEQ,
     *          NIT,OMEGA)
C
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*4 DA(*)
      DIMENSION KCOL(*),KLD(*),DX(*),DB(*),DD(*),KNPR(*)
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
        AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
       END DO
       AUX=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
       DX(IEQ)=AUX
2     CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 3 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
       DD(IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
       END DO
3     CONTINUE
C
      CALL E011Sum(DD)
C
      DO 4 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *          OMEGA*(DD(IEQ)+DB(IEQ))/E011_UE(IEQ)
4     CONTINUE
!     ----------------------------------------------
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO 22 IEQ=NEQ-1,1,-1
       IF (KNPR(IEQ).NE.0) GOTO 22
       AUX=DB(IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
       END DO
       AUX=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
       DX(IEQ)=AUX
22    CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 33 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 33
       DD(IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
       END DO
33    CONTINUE
C
      CALL E011Sum(DD)
C
      DO 44 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 44
       DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *          OMEGA*(DD(IEQ)+DB(IEQ))/E011_UE(IEQ)
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      END
C
c
C
      SUBROUTINE SSORSolverXXX(DA11,DA22,DA33,KCOL,KLD,
     *           DX,DB,DD,KNPR,NEQ,NIT,OMEGA)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE011,E011_UE,myid,E011SUM
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCOL(*),KLD(*),DX(*),DB(*),DD(*),KNPR(*)
      DIMENSION DA11(*),DA22(*),DA33(*)
      INTEGER MEQ1,MEQ2,MEQ3
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      MEQ1 = 0
      MEQ2 = NEQ
      MEQ3 = 2*NEQ
C
      DO 1 ITE=1,NIT
C
!     The real SOR part is here
!     ----------------------------------------------
      DO 2 IEQ=1,NEQ
       IF (KNPR(IEQ).NE.0) GOTO 2
       AUX1=DB(MEQ1+IEQ)
       AUX2=DB(MEQ2+IEQ)
       AUX3=DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J=KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
2     CONTINUE
!     ----------------------------------------------
C
!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 3 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
       DD(MEQ1+IEQ) = 0d0
       DD(MEQ2+IEQ) = 0d0
       DD(MEQ3+IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
       END DO
3     CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
! C
      DO 4 IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     * OMEGA*(DD(MEQ1+IEQ)+DB(MEQ1+IEQ))/MGE011(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     * OMEGA*(DD(MEQ2+IEQ)+DB(MEQ2+IEQ))/MGE011(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     * OMEGA*(DD(MEQ3+IEQ)+DB(MEQ3+IEQ))/MGE011(ILEV)%UE33(IEQ)
4     CONTINUE
!     ----------------------------------------------
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO 22 IEQ=NEQ,1,-1
       IF (KNPR(IEQ).NE.0) GOTO 22
       AUX1=DB(MEQ1+IEQ)
       AUX2=DB(MEQ2+IEQ)
       AUX3=DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J=KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
22    CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 33 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 33
       DD(MEQ1+IEQ) = 0d0
       DD(MEQ2+IEQ) = 0d0
       DD(MEQ3+IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
       END DO
33    CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
! C
C
      DO 44 IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     * OMEGA*(DD(MEQ1+IEQ)+DB(MEQ1+IEQ))/MGE011(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     * OMEGA*(DD(MEQ2+IEQ)+DB(MEQ2+IEQ))/MGE011(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     * OMEGA*(DD(MEQ3+IEQ)+DB(MEQ3+IEQ))/MGE011(ILEV)%UE33(IEQ)
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE SSORSolverYYY(DA11,DA22,DA33,KCOL,KLD,
     *           DX,DB,DD,KNPR,NEQ,NIT,OMEGA)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE011,E011_UE,myid,E011SUM
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCOL(*),KLD(*),DX(*),DB(*),DD(*),KNPR(*)
      DIMENSION DA11(*),DA22(*),DA33(*)
      INTEGER MEQ1,MEQ2,MEQ3
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      MEQ1 = 0
      MEQ2 = NEQ
      MEQ3 = 2*NEQ
C
      DO 1 ITE=1,NIT
C
!     The real SOR part is here
!     ----------------------------------------------
      DO 2 IEQ=1,NEQ
       IF (KNPR(IEQ).NE.0) GOTO 2
       AUX1=DB(MEQ1+IEQ)
       AUX2=DB(MEQ2+IEQ)
       AUX3=DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J=KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
2     CONTINUE
!     ----------------------------------------------
C
!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 3 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
       DD(MEQ1+IEQ) = 0d0
       DD(MEQ2+IEQ) = 0d0
       DD(MEQ3+IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
       END DO
3     CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
C
      DO 4 IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     * OMEGA*(DD(MEQ1+IEQ)+DB(MEQ1+IEQ))/MGE011(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     * OMEGA*(DD(MEQ2+IEQ)+DB(MEQ2+IEQ))/MGE011(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     * OMEGA*(DD(MEQ3+IEQ)+DB(MEQ3+IEQ))/MGE011(ILEV)%UE33(IEQ)
4     CONTINUE
!     ----------------------------------------------
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO 22 IEQ=NEQ,1,-1
       IF (KNPR(IEQ).NE.0) GOTO 22
       AUX1=DB(MEQ1+IEQ)
       AUX2=DB(MEQ2+IEQ)
       AUX3=DB(MEQ3+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J=KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
22    CONTINUE
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO 33 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 33
       DD(MEQ1+IEQ) = 0d0
       DD(MEQ2+IEQ) = 0d0
       DD(MEQ3+IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
       END DO
33    CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
C
      DO 44 IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     * OMEGA*(DD(MEQ1+IEQ)+DB(MEQ1+IEQ))/MGE011(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     * OMEGA*(DD(MEQ2+IEQ)+DB(MEQ2+IEQ))/MGE011(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     * OMEGA*(DD(MEQ3+IEQ)+DB(MEQ3+IEQ))/MGE011(ILEV)%UE33(IEQ)
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      DO 122 IEQ=1,NEQ
      DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
      DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
      DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
      DO 133 ICOL=KLD(IEQ),KLD(IEQ+1)-1
      J=KCOL(ICOL)
      DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
      DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
      DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
133   CONTINUE
122   CONTINUE
C
      CALL E011Sum(DD(MEQ1+1))
      CALL E011Sum(DD(MEQ2+1))
      CALL E011Sum(DD(MEQ3+1))
      CALL LL21(DD,3*NEQ,DEF)

      END
C
C
C
      SUBROUTINE GetDefNorm(VA,KCOL,KLD,DX,DB,DD,NEQ,DEF)
C
      USE PP3D_MPI, only : e011sum
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
       DO IEQ=1,NEQ
        DD(IEQ) = 0d0
        DO ICOL=KLD(IEQ),KLD(IEQ+1)-1
         DD(IEQ)=DD(IEQ)-REAL(VA(ICOL))*DX(KCOL(ICOL))
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

      END
