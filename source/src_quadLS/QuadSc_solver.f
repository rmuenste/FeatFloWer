C
      SUBROUTINE E013_SSORSmoother9(DA11,DA22,DA33,DA12,DA13,DA23,
     *           DA21,DA31,DA32,KCOL,KLD,DX,DB,DD,KNPR,NEQ,NIT,OMEGA)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE013,myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA11(*),DA22(*),DA33(*)
      DIMENSION DA12(*),DA13(*),DA23(*),DA21(*),DA31(*),DA32(*)
      DIMENSION KCOL(*),KLD(*),DX(*),DB(*),DD(*),KNPR(*)
      COMMON /ERRCTL/ IER,ICHECK
      INTEGER MEQ1,MEQ2,MEQ3
      SAVE /ERRCTL/
C
      MEQ1 = 0
      MEQ2 = NEQ
      MEQ3 = 2*NEQ
C
      DO ITE=1,NIT
C
!     The real SOR part is here
!     ----------------------------------------------
      DO IEQ=1,NEQ
      IF (KNPR(IEQ).EQ.0) THEN
      ICOL=KLD(IEQ)
      J=KCOL(ICOL)
      AUX1=DB(MEQ1+IEQ)-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      AUX2=DB(MEQ2+IEQ)-DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
      AUX3=DB(MEQ3+IEQ)-DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
      DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J=KCOL(ICOL)
      AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)-
     *          DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      AUX2=AUX2-DA21(ICOL)*DX(MEQ1+J)-
     *          DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
      AUX3=AUX3-DA31(ICOL)*DX(MEQ1+J)-
     *          DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
      END DO
      AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
      AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
      AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
      DX(MEQ1+IEQ)=AUX1
      DX(MEQ2+IEQ)=AUX2
      DX(MEQ3+IEQ)=AUX3
      END IF
      END DO
!     ----------------------------------------------

!       CALL ll21(DX(MEQ1+1),NEQ,ddd1)
!       CALL ll21(DX(MEQ2+1),NEQ,ddd2)
!       CALL ll21(DX(MEQ3+1),NEQ,ddd3)
!       WRITE(*,*) myid, ddd1,ddd2,ddd3,ilev
!       pause

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
      ICOL=KLD(IEQ)
      J=KCOL(ICOL)
      DD(MEQ1+IEQ) = DB(MEQ1+IEQ)-
     *               DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      DD(MEQ2+IEQ) = DB(MEQ2+IEQ)-
     *               DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
      DD(MEQ3+IEQ) = DB(MEQ3+IEQ)-
     *               DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
      DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J = KCOL(ICOL)
      DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-
     *DA11(ICOL)*DX(MEQ1+J)-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-
     *DA21(ICOL)*DX(MEQ1+J)-DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
      DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-
     *DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
      END DO
      END DO
C
!       CALL E013Sum(DD(MEQ1+1))
!       CALL E013Sum(DD(MEQ2+1))
!       CALL E013Sum(DD(MEQ3+1))
      CALL E013UVWSum(DD)
C
      DO IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *          OMEGA*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *          OMEGA*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *          OMEGA*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
      END DO
!     ----------------------------------------------

!       CALL ll21(DX,3*NEQ,ddd1)
!       CALL ll21(DX,3*NEQ,ddd2)
!       CALL ll21(DX,3*NEQ,ddd3)
!       WRITE(*,*) myid, ddd1,ddd2,ddd3
!       pause
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO IEQ=NEQ-1,1,-1
      IF (KNPR(IEQ).EQ.0) THEN
      ICOL=KLD(IEQ)
      J=KCOL(ICOL)
      AUX1=DB(MEQ1+IEQ)-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      AUX2=DB(MEQ2+IEQ)-DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
      AUX3=DB(MEQ3+IEQ)-DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
      DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J=KCOL(ICOL)
      AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)-
     *          DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      AUX2=AUX2-DA21(ICOL)*DX(MEQ1+J)-
     *          DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
      AUX3=AUX3-DA31(ICOL)*DX(MEQ1+J)-
     *          DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
      END DO
      AUX1=OMEGA*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
      AUX2=OMEGA*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
      AUX3=OMEGA*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
      DX(MEQ1+IEQ)=AUX1
      DX(MEQ2+IEQ)=AUX2
      DX(MEQ3+IEQ)=AUX3
      END IF
      END DO
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
      ICOL=KLD(IEQ)
      J=KCOL(ICOL)
      DD(MEQ1+IEQ) = DB(MEQ1+IEQ)-
     *               DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      DD(MEQ2+IEQ) = DB(MEQ2+IEQ)-
     *               DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
      DD(MEQ3+IEQ) = DB(MEQ3+IEQ)-
     *               DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
      DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J = KCOL(ICOL)
      DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-
     *DA11(ICOL)*DX(MEQ1+J)-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-
     *DA21(ICOL)*DX(MEQ1+J)-DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
      DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-
     *DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
      END DO
      END DO
C
!       CALL E013Sum(DD(MEQ1+1))
!       CALL E013Sum(DD(MEQ2+1))
!       CALL E013Sum(DD(MEQ3+1))
      CALL E013UVWSum(DD)
C
      DO IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *          OMEGA*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *          OMEGA*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *          OMEGA*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
      END DO
!     ----------------------------------------------
C
      END DO
C
      END
C
C
C
      SUBROUTINE E013_SSORSmoother(DA11,DA22,DA33,KCOL,KLD,DX,DB,
     *           DD,KNPR,NEQ,NIT,OMEGA)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE013,myid
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
!       CALL E013Sum(DD(MEQ1+1))
!       CALL E013Sum(DD(MEQ2+1))
!       CALL E013Sum(DD(MEQ3+1))
      CALL E013UVWSum(DD)
C
      DO 4 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *          OMEGA*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *          OMEGA*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *          OMEGA*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
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
!       CALL E013Sum(DD(MEQ1+1))
!       CALL E013Sum(DD(MEQ2+1))
!       CALL E013Sum(DD(MEQ3+1))
      CALL E013UVWSum(DD)
C
      DO 44 IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *          OMEGA*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *          OMEGA*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *          OMEGA*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      END
C
C
C
       SUBROUTINE E013_SORSmoother(DA,KCOL,KLD,DX,DB,DD,KNPR,NEQ,
     *          NIT,OMEGA)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE013,myid
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
       IF (KNPR(IEQ).EQ.0) GOTO 3
       DD(IEQ) = DB(IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
       END DO
3     CONTINUE
C
      CALL E013Sum(DD)
C
      DO 4 IEQ=1,NEQ
       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *          OMEGA*DD(IEQ)/MGE013(ILEV)%UE(IEQ)
4     CONTINUE
!     ----------------------------------------------
C
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE E013_JacobiSmoother(DA,KCOL,KLD,DX,DB,DD,NEQ,
     *          NIT,OMEGA)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE013,myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      DO 1 ITE=1,NIT
C
      DO 2 IEQ=1,NEQ
      DD(IEQ) = DB(IEQ)
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
2     CONTINUE
C
      CALL E013Sum(DD)
C
!       DO 4 IEQ=1,NEQ
! 4     DD(IEQ)=DD(IEQ)+DB(IEQ)
      DO 5 IEQ=1,NEQ
5     DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *         OMEGA*DD(IEQ)/MGE013(ILEV)%UE(IEQ)
C
1     CONTINUE
C
!       IF (myid.eq.1) write(*,'(2I5,6D12.4)') 
!      * 1,ilev,DX(108),Db(108),Dd(108),DX(271),Db(271),Dd(271)
!       IF (myid.eq.2) write(*,'(2I5,6D12.4)') 
!      * 2,ilev,DX(94),Db(94),Dd(94),DX(261),Db(261),Dd(261)
C
      END
C
C
C
      SUBROUTINE E013_JacobiSolver9(DA11,DA22,DA33,DA12,DA13,DA23,
     *            DA21,DA31,DA32,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA,DEF)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE013,myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA11(*),DA22(*),DA33(*)
      DIMENSION DA12(*),DA13(*),DA23(*)
      DIMENSION DA21(*),DA31(*),DA32(*)
      DIMENSION KCOL(*),KLD(*),DX(*),DB(*),DD(*)
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
      DO 2 IEQ=1,NEQ
       ICOL=KLD(IEQ)
       J=KCOL(ICOL)
       DD(MEQ1+IEQ) = DB(MEQ1+IEQ)-
     *                DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
       DD(MEQ2+IEQ) = DB(MEQ2+IEQ)-
     *                DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
       DD(MEQ3+IEQ) = DB(MEQ3+IEQ)-
     *                DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
       DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-
     * DA11(ICOL)*DX(MEQ1+J)-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-
     * DA21(ICOL)*DX(MEQ1+J)-DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-
     * DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
3     CONTINUE
2     CONTINUE

      CALL E013UVWSum(DD)
C
      DO 5 IEQ=1,NEQ
      DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *         OMEGA*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
      DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *         OMEGA*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
      DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *         OMEGA*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
5     CONTINUE
C
1     CONTINUE
C
      DO 22 IEQ=1,NEQ
      DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
      DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
      DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
      DO 33 ICOL=KLD(IEQ),KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-
     * DA11(ICOL)*DX(MEQ1+J)-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
        DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-
     * DA21(ICOL)*DX(MEQ1+J)-DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
        DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-
     * DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
33    CONTINUE
22    CONTINUE
C
      CALL E013UVWSum(DD)
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
      SUBROUTINE E013_JacobiSolver(DA11,DA22,DA33,KCOL,KLD,DX,
     *           DB,DD,NEQ,NIT,OMEGA,DEF)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE013,myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCOL(*),KLD(*),DX(*),DB(*),DD(*)
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
      DO 2 IEQ=1,NEQ
      DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
      DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
      DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J=KCOL(ICOL)
      DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
      DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
      DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
3     CONTINUE
2     CONTINUE
C
      CALL E013UVWSum(DD)
C
      DO 5 IEQ=1,NEQ
      DX(MEQ1+IEQ)=(1D0-OMEGA)*DX(MEQ1+IEQ)+
     *         OMEGA*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
      DX(MEQ2+IEQ)=(1D0-OMEGA)*DX(MEQ2+IEQ)+
     *         OMEGA*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
      DX(MEQ3+IEQ)=(1D0-OMEGA)*DX(MEQ3+IEQ)+
     *         OMEGA*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
5     CONTINUE
C
1     CONTINUE
C
      DO 22 IEQ=1,NEQ
      DD(MEQ1+IEQ) = DB(MEQ1+IEQ)
      DD(MEQ2+IEQ) = DB(MEQ2+IEQ)
      DD(MEQ3+IEQ) = DB(MEQ3+IEQ)
      DO 33 ICOL=KLD(IEQ),KLD(IEQ+1)-1
      J=KCOL(ICOL)
      DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-DA11(ICOL)*DX(MEQ1+J)
      DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-DA22(ICOL)*DX(MEQ2+J)
      DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-DA33(ICOL)*DX(MEQ3+J)
33    CONTINUE
22    CONTINUE
C
      CALL E013UVWSum(DD)
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
      SUBROUTINE E012_JacobiSolver(DA,KCOL,KLD,DAP,KCOLP,KLDP,
     *           DXP,DX,DB,DD,NEQ,def0,def,NIT,OMEGA)
C
      USE PP3D_MPI, ONLY :myid,COMM_Maximum,COMM_NLComplete
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
      DIMENSION DAP(*),KCOLP(*),KLDP(*),DXP(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (myid.ne.0) THEN

      CALL LL21 (DB,NEQ,def0)
C
      CALL COMM_Maximum(def0)

      DO 1 ITE=1,NIT
C
      INLComplete = 0
      CALL GetParPressure(DX,DXP)
C
      CALL LCL1(DD,NEQ)
      DO 2 IEQ=1,NEQ
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
2     CONTINUE
C
      DO 20 IEQ=1,NEQ
      DO 30 ICOL=KLDP(IEQ),KLDP(IEQ+1)-1
30     DD(IEQ)=DD(IEQ)-DAP(ICOL)*DXP(KCOLP(ICOL))
20     CONTINUE
C
      DO 4 IEQ=1,NEQ
4     DD(IEQ)=DD(IEQ)+DB(IEQ)
      DO 5 IEQ=1,NEQ
      DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+
     *         OMEGA*DD(IEQ)/DA(KLD(IEQ))!
5     CONTINUE
C
      IF (MOD(ITE,100).EQ.0) THEN
C     ! check the convergence
      DO 6 IEQ=1,NEQ
       DD(IEQ)=0d0
       DO 7 ICOL=KLD(IEQ),KLD(IEQ+1)-1
7       DD(IEQ)=DD(IEQ)+DA(ICOL)*DX(KCOL(ICOL))
6     CONTINUE
C
      DO 60 IEQ=1,NEQ
      DO 70 ICOL=KLDP(IEQ),KLDP(IEQ+1)-1
70     DD(IEQ)=DD(IEQ)+DAP(ICOL)*DXP(KCOLP(ICOL))
60    CONTINUE
C
      DO 8 IEQ=1,NEQ
       DD(IEQ)=DB(IEQ)-DD(IEQ)
8     CONTINUE

      CALL LL21 (DD,NEQ,def)
      ddd = def
      IF (def0.lt.1d-32) THEN
       def = 0d0
      ELSE
       def = def/def0
      END IF

      CALL COMM_Maximum(def)
      CALL COMM_Maximum(ddd)
      IF (def.lt.1d-6) INLComplete = 1
      IF (ddd.lt.1d-9) INLComplete = 1
!       write(*,'(2I5,3D12.4)') myid,ITE,def,ddd,def0
      CALL COMM_NLComplete(INLComplete)
      IF (INLComplete.eq.1) GOTO 88

      END IF

1     CONTINUE

      ELSE

      CALL COMM_Maximum(def0)
      DO ITE=1,NIT
       IF (MOD(ITE,100).EQ.0) THEN
        CALL COMM_Maximum(def)
        CALL COMM_Maximum(ddd)
        CALL COMM_NLComplete(INLComplete)
        IF (INLComplete.eq.1) GOTO 88
       END IF
      END DO

      END IF
C
88    CONTINUE
C
      NIT = ITE
C
      END
C
C
C 
      SUBROUTINE E012_BiCGStabSmoother(DX,DXP,DB,NEQ,NIT,DAX0,
     *DCG0C,BNOCON,DR,DR0,DP,DPA,DSA)

      USE PP3D_MPI, ONLY :E011Sum,E011Mean,myid,COMM_SUMM,COMM_Maximum,
     *                    COMM_NLComplete
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DXP(*),DB(*),DR(*),DR0(*)
      DIMENSION DP(*),DPA(*),DSA(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA DEF/0D0/

      IF (myid.ne.0) THEN
C
C *** Initialization
      RHO0  =1D0
      DALPHA=0D0
      OMEGA0=1D0
C
      CALL LCP1(DB,DR,NEQ)
      CALL DAX0(DX,DXP,DR,NEQ,-1D0,1D0)
      IF (BNOCON) CALL DCG0C(DR,DXP,NEQ)
C
      CALL LCP1(DR,DR0,NEQ)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      INLComplete = 0
      CALL LSP1(DR0,DR,NEQ,RHO1)
      CALL COMM_SUMM(RHO1)
      DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
      RHO0 =RHO1
C
      CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
      CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)
C
      CALL DAX0(DP,DXP,DPA,NEQ,1D0,0D0)
      IF (BNOCON) CALL DCG0C(DPA,DXP,NEQ)
C
      CALL LSP1(DR0,DPA,NEQ,DALPHA)
      CALL COMM_SUMM(DALPHA)
      DALPHA=RHO1/DALPHA
C
      CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)
C
      CALL DAX0(DR,DXP,DSA,NEQ,1D0,0D0)
      IF (BNOCON) CALL DCG0C(DSA,DXP,NEQ)
C
      CALL LSP1(DSA,DR ,NEQ,OMEGA1)
      CALL COMM_SUMM(OMEGA1)
      CALL LSP1(DSA,DSA,NEQ,OMEGA2)
      CALL COMM_SUMM(OMEGA2)
      OMEGA0=OMEGA1/OMEGA2
C
      CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
      CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)
C
      CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)
C
100   CONTINUE
C
      ELSE ! Myid.eq.0

      DO 200 ITE=1,NIT
       CALL COMM_SUMM(RHO1)
       CALL COMM_SUMM(DALPHA)
       CALL COMM_SUMM(OMEGA1)
       CALL COMM_SUMM(OMEGA2)

200   CONTINUE
C
      END IF
C
88    CONTINUE
C
      NIT = ITE-1
C
99999 END
C
C
C
      SUBROUTINE E012_BiCGStabSolver(DX,DXP,DB,NEQ,NIT,DAX0,
     *DCG0C,BNOCON,DR,DR0,DP,DPA,DSA,DefDropCrit)

      USE PP3D_MPI, ONLY :E011Sum,E011Mean,myid,COMM_SUMM,COMM_Maximum,
     *                    COMM_NLComplete
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DXP(*),DB(*),DR(*),DR0(*)
      DIMENSION DP(*),DPA(*),DSA(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA DEF/0D0/

      IF (myid.ne.0) THEN
C
!       CALL LCP1(DB,DR,NEQ)
!       IF (BNOCON) CALL DCG0C(DR,NEQ)
C
C *** Initialization
      RHO0  =1D0
      DALPHA=0D0
      OMEGA0=1D0
C
      CALL LCP1(DB,DR,NEQ)
      CALL DAX0(DX,DXP,DR,NEQ,-1D0,1D0)
      IF (BNOCON) CALL DCG0C(DR,DXP,NEQ)
      CALL LL21(DR,NEQ,DEF0)
      CALL COMM_Maximum(DEF0)
C
      CALL LCP1(DR,DR0,NEQ)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      INLComplete = 0
      CALL LSP1(DR0,DR,NEQ,RHO1)
      CALL COMM_SUMM(RHO1)
      DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
      RHO0 =RHO1
C
      CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
      CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)
C
      CALL DAX0(DP,DXP,DPA,NEQ,1D0,0D0)
      IF (BNOCON) CALL DCG0C(DPA,DXP,NEQ)
C
      CALL LSP1(DR0,DPA,NEQ,DALPHA)
      CALL COMM_SUMM(DALPHA)
      DALPHA=RHO1/DALPHA
C
      CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)
C
      CALL DAX0(DR,DXP,DSA,NEQ,1D0,0D0)
      IF (BNOCON) CALL DCG0C(DSA,DXP,NEQ)
C
      CALL LSP1(DSA,DR ,NEQ,OMEGA1)
      CALL COMM_SUMM(OMEGA1)
      CALL LSP1(DSA,DSA,NEQ,OMEGA2)
      CALL COMM_SUMM(OMEGA2)
      OMEGA0=OMEGA1/OMEGA2
C
      CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
      CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)
C
      CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)
C
      CALL LL21(DR,NEQ,DEF)
C
      IF (DEF0.GT.1d-32) THEN
       DefDrop = DEF/DEF0
      ELSE
       DefDrop = 0d0
      END IF
C
      CALL COMM_Maximum(DefDrop)
      CALL COMM_Maximum(Def)
C
      IF (DefDrop.lt.DefDropCrit) INLComplete = 1
!      IF (Def.lt.1d-11) INLComplete = 1
      CALL COMM_NLComplete(INLComplete)
      IF (INLComplete.eq.1) GOTO 88
C
100   CONTINUE
C
      ELSE ! Myid.eq.0

      CALL COMM_Maximum(DEF0)

      DO 200 ITE=1,NIT
       CALL COMM_SUMM(RHO1)
       CALL COMM_SUMM(DALPHA)
       CALL COMM_SUMM(OMEGA1)
       CALL COMM_SUMM(OMEGA2)

       CALL COMM_Maximum(DefDrop)
       CALL COMM_Maximum(Def)

       CALL COMM_NLComplete(INLComplete)
       IF (INLComplete.eq.1) GOTO 88
200   CONTINUE
C
      END IF
C
88    CONTINUE
C
      NIT = ITE-1
C
99999 END
C
C
C
      SUBROUTINE E012_BiCGStabSolverMaster(DX,DB,NEQ,NIT,DAX0,
     *DCG0C,BNOCON,DR,DR0,DP,DPA,DSA,DefDropCrit)

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DB(*),DR(*),DR0(*)
      DIMENSION DP(*),DPA(*),DSA(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA DEF/0D0/

C *** Initialization
      RHO0  =1D0
      DALPHA=0D0
      OMEGA0=1D0
C
      CALL LCP1(DB,DR,NEQ)
      CALL DAX0(DX,DR,NEQ,-1D0,1D0)
      IF (BNOCON) CALL DCG0C(DR,NEQ)
      CALL LL21(DR,NEQ,DEF0)
C
      CALL LCP1(DR,DR0,NEQ)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      INLComplete = 0
      CALL LSP1(DR0,DR,NEQ,RHO1)
      DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
      RHO0 =RHO1
C
      CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
      CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)
C
      CALL DAX0(DP,DPA,NEQ,1D0,0D0)
      IF (BNOCON) CALL DCG0C(DPA,NEQ)
C
      CALL LSP1(DR0,DPA,NEQ,DALPHA)
      DALPHA=RHO1/DALPHA
C
      CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)
C
      CALL DAX0(DR,DSA,NEQ,1D0,0D0)
      IF (BNOCON) CALL DCG0C(DSA,NEQ)
C
      CALL LSP1(DSA,DR ,NEQ,OMEGA1)
      CALL LSP1(DSA,DSA,NEQ,OMEGA2)
      OMEGA0=OMEGA1/OMEGA2
C
      CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
      CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)
C
      CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)
C
      CALL LL21(DR,NEQ,DEF)
C
      IF (DEF0.GT.1d-32) THEN
       DefDrop = DEF/DEF0
      ELSE
       DefDrop = 0d0
      END IF
C
      IF (DefDrop.lt.DefDropCrit) GOTO 88
C
100   CONTINUE
C
88    CONTINUE
C
      NIT = ITE-1
C
!       WRITE(*,*) NIT
99999 END
C
C
C
      SUBROUTINE PARID117(DA,DAP,KCOL,KPCOL,KLD,KPLD,DX,DXP,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*)
      DIMENSION DAP(*),KPCOL(*),KPLD(*),DXP(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID117 ','01/02/89')
C
      DO 1 IEQ=1,NEQ
      AUX=0D0
      ILD=KLD(IEQ)
      DO 4 ICOL=KPLD(IEQ),KPLD(IEQ+1)-1
4     AUX=AUX+DAP(ICOL)*DXP(KPCOL(ICOL))
      DO 5 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).GE.IEQ) GOTO 1
5     AUX=AUX+DA(ICOL)*DX(KCOL(ICOL))
1     DX(IEQ)=(DX(IEQ)-AUX*OMEGA)/DA(ILD)
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      ILD=KLD(IEQ)
      DO 12 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).LE.IEQ) GOTO 12
      AUX=AUX+DA(ICOL)*DX(KCOL(ICOL))
12    CONTINUE
10    DX(IEQ)=DX(IEQ)-AUX*OMEGA/DA(ILD)
C
      END
C
C
C
      SUBROUTINE E012_SOR_Smoother(DD,DA,DAP,KCOL,KPCOL,
     *           KLD,KPLD,DX,DXP,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DD(*)
      DIMENSION DAP(*),KPCOL(*),KPLD(*),DXP(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID117 ','01/02/89')
C
      DO IEQ=1,NEQ
       AUX=DD(IEQ)
       DO ICOL=KPLD(IEQ),KPLD(IEQ+1)-1
        AUX=AUX-DAP(ICOL)*DXP(KPCOL(ICOL))
       END DO
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
       END DO
       DX(IEQ)=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
      END DO
C
      END
C
C
C
      SUBROUTINE crs_E012_SOR_Smoother(DD,DA,KCOL,
     *           KLD,DX,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID117 ','01/02/89')
C
      DO IEQ=1,NEQ
       AUX=DD(IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
       END DO
       DDAUX =  DX(IEQ)
       DX(IEQ)=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
!        WRITE(*,*) DDAUX ,DX(IEQ)
      END DO
C
      END
C
C
C
      SUBROUTINE E012_ROS_Smoother(DD,DA,DAP,KCOL,KPCOL,
     *           KLD,KPLD,DX,DXP,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DD(*)
      DIMENSION DAP(*),KPCOL(*),KPLD(*),DXP(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID117 ','01/02/89')
C
      DO IEQ=NEQ-1,1,-1
       AUX=DD(IEQ)
       DO ICOL=KPLD(IEQ),KPLD(IEQ+1)-1
        AUX=AUX-DAP(ICOL)*DXP(KPCOL(ICOL))
       END DO
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
       END DO
       DX(IEQ)=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
      END DO
C
      END
C
C
C
      SUBROUTINE PARID117_MASTER(DA,KCOL,KLD,DX,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID117 ','01/02/89')
C
      DO 1 IEQ=1,NEQ
      AUX=0D0
      ILD=KLD(IEQ)
      DO 5 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).GE.IEQ) GOTO 1
5     AUX=AUX+DA(ICOL)*DX(KCOL(ICOL))
1     DX(IEQ)=(DX(IEQ)-AUX*OMEGA)/DA(ILD)
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      ILD=KLD(IEQ)
      DO 12 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).LE.IEQ) GOTO 12
      AUX=AUX+DA(ICOL)*DX(KCOL(ICOL))
12    CONTINUE
10    DX(IEQ)=DX(IEQ)-AUX*OMEGA/DA(ILD)
C
      END
C
C
C
      SUBROUTINE E013_BiCGStabSolver(DX,DB,NEQ,NIT,DAX0,
     *DCG0C,BNOCON,DR,DR0,DP,DPA,DSA,DefDropCrit)

      USE PP3D_MPI, ONLY :myid,COMM_SUMM,COMM_Maximum,
     *                    COMM_NLComplete
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DB(*),DR(*),DR0(*)
      DIMENSION DP(*),DPA(*),DSA(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA DEF/0D0/

      IF (myid.ne.0) THEN
C
C *** Initialization
      RHO0  =1D0
      DALPHA=0D0
      OMEGA0=1D0
C
      CALL LCP1(DB,DR,NEQ)
      CALL DAX0(DX,DR,-1D0,1D0)
      CALL E013UVWSum(DR)
      IF (BNOCON) CALL DCG0C(DR,NEQ)
      CALL LL21(DR,NEQ,DEF0)
      CALL COMM_Maximum(DEF0)
C
      CALL LCP1(DR,DR0,NEQ)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      INLComplete = 0
      CALL LSP1(DR0,DR,NEQ,RHO1)
      CALL COMM_SUMM(RHO1)
      DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
      RHO0 =RHO1
C
      CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
      CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)
C
      CALL DAX0(DP,DPA,1D0,0D0)
      CALL E013UVWSum(DPA)
      IF (BNOCON) CALL DCG0C(DPA,NEQ)
C
      CALL LSP1(DR0,DPA,NEQ,DALPHA)
      CALL COMM_SUMM(DALPHA)
      DALPHA=RHO1/DALPHA
C
      CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)
C
      CALL DAX0(DR,DSA,1D0,0D0)
      CALL E013UVWSum(DSA)
      IF (BNOCON) CALL DCG0C(DSA,NEQ)
C
      CALL LSP1(DSA,DR ,NEQ,OMEGA1)
      CALL COMM_SUMM(OMEGA1)
      CALL LSP1(DSA,DSA,NEQ,OMEGA2)
      CALL COMM_SUMM(OMEGA2)
      OMEGA0=OMEGA1/OMEGA2
C
      CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
      CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)
C
      CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)
C
      CALL LL21(DR,NEQ,DEF)
C
      IF (DEF0.GT.1d-32) THEN
       DefDrop = DEF/DEF0
      ELSE
       DefDrop = 0d0
      END IF
C
      CALL COMM_Maximum(DefDrop)
      CALL COMM_Maximum(Def)
      IF (myid.eq.1) then
!        WRITE(*,'(I,10ES12.4)')ite, DefDrop
      end if
C
      IF (DefDrop.lt.DefDropCrit) INLComplete = 1
!      IF (Def.lt.1d-11) INLComplete = 1
      CALL COMM_NLComplete(INLComplete)
      IF (ite.gt.4.and.INLComplete.eq.1) GOTO 88
C
100   CONTINUE
C
      ELSE ! Myid.eq.0

      CALL COMM_Maximum(DEF0)

      DO 200 ITE=1,NIT
       CALL COMM_SUMM(RHO1)
       CALL COMM_SUMM(DALPHA)
       CALL COMM_SUMM(OMEGA1)
       CALL COMM_SUMM(OMEGA2)

       CALL COMM_Maximum(DefDrop)
       CALL COMM_Maximum(Def)

       CALL COMM_NLComplete(INLComplete)
       IF (ite.gt.4.and.INLComplete.eq.1) GOTO 88
200   CONTINUE
C
      END IF
C
88    CONTINUE
C
!       if (myid.eq.1) write(*,*) 'bicgstab: ',ITE,DEF0,DEF
!       NIT = ITE
      
      NIT = ITE-1
C
99999 END
C
C
C
      SUBROUTINE PARID117_E013(DA11,DA22,DA33,KCOL,KLD,DX,DD,
     *           KNPR,NEQ,OMEGA1,OMEGA2)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE013,myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCOL(*),KLD(*),DX(*),DD(*),KNPR(*)
      DIMENSION DA11(*),DA22(*),DA33(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID117 ','01/02/89')
C
!      return
C      
      MEQ1 = 0
      MEQ2 = NEQ
      MEQ3 = 2*NEQ
C
      DO 1 ITE=1,1
C
!     The real SOR part is here
!     ----------------------------------------------
      DO 2 IEQ=1,NEQ
       IF (KNPR(IEQ).NE.0) GOTO 2
       AUX1=0d0
       AUX2=0d0
       AUX3=0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J=KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA2*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA2*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA2*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
       DX(MEQ1+IEQ)=AUX1
       DX(MEQ2+IEQ)=AUX2
       DX(MEQ3+IEQ)=AUX3
2     CONTINUE
!     ----------------------------------------------

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
      CALL E013UVWSum(DD)
C
      DO 4 IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(MEQ1+IEQ)=(1D0-OMEGA1)*DX(MEQ1+IEQ)+
     *          OMEGA1*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA1)*DX(MEQ2+IEQ)+
     *          OMEGA1*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA1)*DX(MEQ3+IEQ)+
     *          OMEGA1*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
4     CONTINUE
!     ----------------------------------------------
C
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO 22 IEQ=NEQ-1,1,-1
       IF (KNPR(IEQ).NE.0) GOTO 22
       AUX1=0d0
       AUX2=0d0
       AUX3=0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)
        AUX2=AUX2-DA22(ICOL)*DX(MEQ2+J)
        AUX3=AUX3-DA33(ICOL)*DX(MEQ3+J)
       END DO
       AUX1=OMEGA2*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
       AUX2=OMEGA2*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
       AUX3=OMEGA2*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
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
      CALL E013UVWSum(DD)
C
      DO 44 IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA1)*DX(MEQ1+IEQ)+
     *          OMEGA1*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA1)*DX(MEQ2+IEQ)+
     *          OMEGA1*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA1)*DX(MEQ3+IEQ)+
     *          OMEGA1*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
44    CONTINUE
!     ----------------------------------------------
C
1     CONTINUE

      END
C
C
C
      SUBROUTINE PARID117_E0139(DA11,DA22,DA33,DA12,DA13,DA23,
     *DA21,DA31,DA32,KCOL,KLD,DX,DD,KNPR,NEQ,OMEGA1,OMEGA2)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE013,myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA11(*),DA22(*),DA33(*)
      DIMENSION DA12(*),DA13(*),DA23(*),DA21(*),DA31(*),DA32(*)
      DIMENSION KCOL(*),KLD(*),DX(*),DD(*),KNPR(*)
      COMMON /ERRCTL/ IER,ICHECK
      INTEGER MEQ1,MEQ2,MEQ3
      SAVE /ERRCTL/
C
!      return
C      
      MEQ1 = 0
      MEQ2 = NEQ
      MEQ3 = 2*NEQ
C
      DO ITE=1,1!NIT
C
!     The real SOR part is here
!     ----------------------------------------------
      DO IEQ=1,NEQ
      IF (KNPR(IEQ).EQ.0) THEN
      ICOL=KLD(IEQ)
      J=KCOL(ICOL)
      AUX1=-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      AUX2=-DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
      AUX3=-DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
      DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J=KCOL(ICOL)
      AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)-
     *          DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      AUX2=AUX2-DA21(ICOL)*DX(MEQ1+J)-
     *          DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
      AUX3=AUX3-DA31(ICOL)*DX(MEQ1+J)-
     *          DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
      END DO
      AUX1=OMEGA2*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
      AUX2=OMEGA2*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
      AUX3=OMEGA2*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
      DX(MEQ1+IEQ)=AUX1
      DX(MEQ2+IEQ)=AUX2
      DX(MEQ3+IEQ)=AUX3
      END IF
      END DO
!
!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
      ICOL=KLD(IEQ)
      J=KCOL(ICOL)
      DD(MEQ1+IEQ) = 0d0-
     *               DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      DD(MEQ2+IEQ) = 0d0-
     *               DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
      DD(MEQ3+IEQ) = 0d0-
     *               DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
      DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J = KCOL(ICOL)
      DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-
     *DA11(ICOL)*DX(MEQ1+J)-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-
     *DA21(ICOL)*DX(MEQ1+J)-DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
      DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-
     *DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
      END DO
      END DO
C
      CALL E013UVWSum(DD)
C
      DO IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 4
       DX(MEQ1+IEQ)=(1D0-OMEGA1)*DX(MEQ1+IEQ)+
     *          OMEGA1*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA1)*DX(MEQ2+IEQ)+
     *          OMEGA1*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA1)*DX(MEQ3+IEQ)+
     *          OMEGA1*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
      END DO
!     ----------------------------------------------
!
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO IEQ=NEQ-1,1,-1
      IF (KNPR(IEQ).EQ.0) THEN
      ICOL=KLD(IEQ)
      J=KCOL(ICOL)
      AUX1=0d0-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      AUX2=0d0-DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
      AUX3=0d0-DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
      DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J=KCOL(ICOL)
      AUX1=AUX1-DA11(ICOL)*DX(MEQ1+J)-
     *          DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      AUX2=AUX2-DA21(ICOL)*DX(MEQ1+J)-
     *          DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
      AUX3=AUX3-DA31(ICOL)*DX(MEQ1+J)-
     *          DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
      END DO
      AUX1=OMEGA2*(AUX1/DA11(KLD(IEQ))-DX(MEQ1+IEQ))+DX(MEQ1+IEQ)
      AUX2=OMEGA2*(AUX2/DA22(KLD(IEQ))-DX(MEQ2+IEQ))+DX(MEQ2+IEQ)
      AUX3=OMEGA2*(AUX3/DA33(KLD(IEQ))-DX(MEQ3+IEQ))+DX(MEQ3+IEQ)
      DX(MEQ1+IEQ)=AUX1
      DX(MEQ2+IEQ)=AUX2
      DX(MEQ3+IEQ)=AUX3
      END IF
      END DO
!     ----------------------------------------------

!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO IEQ=1,NEQ
!       IF (KNPR(IEQ).EQ.0) GOTO 3
      ICOL=KLD(IEQ)
      J=KCOL(ICOL)
      DD(MEQ1+IEQ) = 0d0-
     *               DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      DD(MEQ2+IEQ) = 0d0-
     *               DA21(ICOL)*DX(MEQ1+J)-DA23(ICOL)*DX(MEQ3+J)
      DD(MEQ3+IEQ) = 0d0-
     *               DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)
      DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      J = KCOL(ICOL)
      DD(MEQ1+IEQ)=DD(MEQ1+IEQ)-
     *DA11(ICOL)*DX(MEQ1+J)-DA12(ICOL)*DX(MEQ2+J)-DA13(ICOL)*DX(MEQ3+J)
      DD(MEQ2+IEQ)=DD(MEQ2+IEQ)-
     *DA21(ICOL)*DX(MEQ1+J)-DA22(ICOL)*DX(MEQ2+J)-DA23(ICOL)*DX(MEQ3+J)
      DD(MEQ3+IEQ)=DD(MEQ3+IEQ)-
     *DA31(ICOL)*DX(MEQ1+J)-DA32(ICOL)*DX(MEQ2+J)-DA33(ICOL)*DX(MEQ3+J)
      END DO
      END DO
C
      CALL E013UVWSum(DD)
C
      DO IEQ=1,NEQ
       DX(MEQ1+IEQ)=(1D0-OMEGA1)*DX(MEQ1+IEQ)+
     *          OMEGA1*DD(MEQ1+IEQ)/MGE013(ILEV)%UE11(IEQ)
       DX(MEQ2+IEQ)=(1D0-OMEGA1)*DX(MEQ2+IEQ)+
     *          OMEGA1*DD(MEQ2+IEQ)/MGE013(ILEV)%UE22(IEQ)
       DX(MEQ3+IEQ)=(1D0-OMEGA1)*DX(MEQ3+IEQ)+
     *          OMEGA1*DD(MEQ3+IEQ)/MGE013(ILEV)%UE33(IEQ)
      END DO
!     ----------------------------------------------
C
      END DO
C
      END
C
C
C
