      SUBROUTINE Get_CMat(M,C,CLD,CCOL,BX,BY,BZ,BLD,BCOL,
     *           BTX,BTY,BTZ,BTLD,BTCOL,KNPRU,KNPRV,KNPRW,NEL,NDOF)
      USE PP3D_MPI, ONLY:myid
      IMPLICIT NONE
      REAL*8  C(*),BX(*),BY(*),BZ(*),BTX(*),BTY(*),BTZ(*),M(*)
      INTEGER CLD(*),CCOL(*),BLD(*),BCOL(*),BTLD(*),BTCOL(*)
      INTEGER KNPRU(*),KNPRV(*),KNPRW(*),NDOF,NEL
      INTEGER IEL,JEL,KEL,IDOF,IC,IB,JB,II,JJ
      REAL*8  DBTX,DBTY,DBTZ,DBX,DBY,DBZ,DMU,DMV,DMW,DC,DAUX
      LOGICAL BNEU

      JJ = 0
      DO IEL = 1, NEL
       DO IC = CLD(IEL),CLD(IEL+1)-1
        JEL = CCOL(IC)
!       Filling in the entry (IEL,JEL) in C matrix 
        DC = 0d0
        II=0
        DO IB = BTLD(IEL),BTLD(IEL+1)-1
         IDOF = BTCOL(IB)
          DO JB=BLD(IDOF),BLD(IDOF+1)-1
           KEL =  BCOL(JB)
           IF (KEL.EQ.JEL) GOTO 1
          END DO
          GOTO 2
1         CONTINUE
          II=II+1
          DBTX  = BTX(IB)
          DBTY  = BTY(IB)
          DBTZ  = BTZ(IB)
          IF (KNPRU(IDOF).EQ.1) THEN
           DMU    = 0d0
          ELSE
           DMU    = 1d0/M(IDOF)
          END IF
          IF (KNPRV(IDOF).EQ.1) THEN
           DMV    = 0d0
          ELSE
           DMV    = 1d0/M(IDOF)
          END IF
          IF (KNPRW(IDOF).EQ.1) THEN
           DMW    = 0d0
          ELSE
           DMW    = 1d0/M(IDOF)
          END IF
          DBX  =  BX(JB)
          DBY  =  BY(JB)
          DBZ  =  BZ(JB)
          DC = DC + DBTX*DBX*DMU + DBTY*DBY*DMV + DBTZ*DBZ*DMW
2         CONTINUE
        END DO
        JJ = JJ + II
        C(IC) = DC
       END DO
      END DO

      END
C
C
C
      SUBROUTINE Get_CMatStruct(CLD,CCOL,BLD,BCOL,KVERT,NEL,IN,iP)
      IMPLICIT NONE
      INTEGER iP,IN
      INTEGER CLD(*),CCOL(IN),BLD(*),BCOL(*),KVERT(8,*),NEL
      INTEGER IEL,JEL,KEL,IVT,JVT,IVERT,IB,IC,JC,IEQ,NNEL
      INTEGER ELEMS(100)
      LOGICAL bNEW

      IF (IP.EQ.1) THEN
      CLD(1) = 1
      DO IEL=1,NEL
       NNEL = 0
       ELEMS = 0
       DO IVT=1,8
        IVERT = KVERT(IVT,IEL)
        DO IB=BLD(IVERT),BLD(IVERT+1)-1,4
         JEL = (BCOL(IB)-1)/4+1
         bNEW = .TRUE.
         DO KEL=1,NNEL
          IF (ELEMS(KEL).EQ.JEL) bNEW=.FALSE.
         END DO
         IF (bNEW.AND.JEL.NE.IEL) THEN
          NNEL = NNEL + 1
          ELEMS(NNEL) = JEL
         END IF
        END DO
       END DO
       CALL SORT1D(ELEMS(1:NNEL),NNEL)
       IEQ = 4*(IEL-1)+1

       DO IVT = 1,4
        CLD(IEQ+IVT) = CLD(IEQ+IVT-1) + 4*(NNEL+1)
       END DO

       CCOL(CLD(IEQ+0)+0) = IEQ + 0
       CCOL(CLD(IEQ+0)+1) = IEQ + 1
       CCOL(CLD(IEQ+0)+2) = IEQ + 2
       CCOL(CLD(IEQ+0)+3) = IEQ + 3

       CCOL(CLD(IEQ+1)+0) = IEQ + 1
       CCOL(CLD(IEQ+1)+1) = IEQ + 0
       CCOL(CLD(IEQ+1)+2) = IEQ + 2
       CCOL(CLD(IEQ+1)+3) = IEQ + 3

       CCOL(CLD(IEQ+2)+0) = IEQ + 2
       CCOL(CLD(IEQ+2)+1) = IEQ + 0
       CCOL(CLD(IEQ+2)+2) = IEQ + 1
       CCOL(CLD(IEQ+2)+3) = IEQ + 3

       CCOL(CLD(IEQ+3)+0) = IEQ + 3
       CCOL(CLD(IEQ+3)+1) = IEQ + 0
       CCOL(CLD(IEQ+3)+2) = IEQ + 1
       CCOL(CLD(IEQ+3)+3) = IEQ + 2

       DO IVT=0,3
        DO JEL=1,NNEL
         IC = CLD(IEQ+IVT)+4*JEL
         JC = 4*(ELEMS(JEL)-1)+1
         IF(IC.GT.IN) write(*,*) "1-problem",IN
         DO JVT=0,3
          CCOL(IC+JVT) = JC + JVT
         END DO
        END DO
       END DO
      END DO

      ELSE

      CLD(1) = 1
      DO IEL=1,NEL
       NNEL = 0
       ELEMS = 0
       DO IVT=1,8
        IVERT = KVERT(IVT,IEL)
        DO IB=BLD(IVERT),BLD(IVERT+1)-1,4
         JEL = (BCOL(IB)-1)/4+1
         bNEW = .TRUE.
         DO KEL=1,NNEL
          IF (ELEMS(KEL).EQ.JEL) bNEW=.FALSE.
         END DO
         IF (bNEW) THEN
          NNEL = NNEL + 1
          ELEMS(NNEL) = JEL
         END IF
        END DO
       END DO
       IEQ = 4*(IEL-1)+1
       IF (NNEL.GT.0) THEN
        CALL SORT1D(ELEMS(1:NNEL),NNEL)
        DO IVT = 1,4
         CLD(IEQ+IVT) = CLD(IEQ+IVT-1) + 4*NNEL
        END DO

        DO IVT=0,3
         DO JEL=1,NNEL
          IC = CLD(IEQ+IVT)+4*(JEL-1)
          JC = 4*(ELEMS(JEL)-1)+1
          IF(IC.GT.IN) write(*,*) "2-problem",IN,CLD(ieq+ivt+1)
          DO JVT=0,3
           CCOL(IC+JVT) = JC + JVT
          END DO
         END DO
        END DO
       ELSE
        DO IVT = 1,4
         CLD(IEQ+IVT) = CLD(IEQ+IVT-1)
        END DO
       END IF

      END DO

      END IF

      CONTAINS

      SUBROUTINE SORT1D(LW,N)
      IMPLICIT NONE
      INTEGER LW(N),LWA
      INTEGER I,J,N

      DO I=2,N
      DO J=N,I,-1
      IF (LW(J).LT.LW(J-1)) THEN
       LWA     = LW(J)
       LW(J)   = LW(J-1)
       LW(J-1) = LWA
      END IF
      END DO
      END DO

      END SUBROUTINE SORT1D

      END
C
C
C
      SUBROUTINE Get_CMatLen(BLD,BCOL,KVERT,NEL,NNNEL,iP)
      IMPLICIT NONE
      INTEGER iP
      INTEGER BLD(*),BCOL(*),KVERT(8,*),NEL,NNNEL
      INTEGER IEL,JEL,KEL,IVT,IVERT,IB,NNEL
      INTEGER ELEMS(100)
      LOGICAL bNEW

      IF (iP.EQ.1) THEN
      NNNEL = 0
      DO IEL=1,NEL
       NNEL = 0
       ELEMS = 0
       DO IVT=1,8
        IVERT = KVERT(IVT,IEL)
        DO IB=BLD(IVERT),BLD(IVERT+1)-1,4
         JEL = (BCOL(IB)-1)/4+1
         bNEW = .TRUE.
         DO KEL=1,NNEL
          IF (ELEMS(KEL).EQ.JEL) bNEW=.FALSE.
         END DO
         IF (bNEW.AND.JEL.NE.IEL) THEN
          NNEL = NNEL + 1
          ELEMS(NNEL) = JEL
         END IF
        END DO
       END DO
       NNNEL = NNNEL + NNEL
      END DO

      NNNEL = 4*4*(NNNEL + NEL)

      ELSE

      NNNEL = 0
      DO IEL=1,NEL
       NNEL = 0
       ELEMS = 0
       DO IVT=1,8
        IVERT = KVERT(IVT,IEL)
        DO IB=BLD(IVERT),BLD(IVERT+1)-1,4
         JEL = (BCOL(IB)-1)/4+1
         bNEW = .TRUE.
         DO KEL=1,NNEL
          IF (ELEMS(KEL).EQ.JEL) bNEW=.FALSE.
         END DO
         IF (bNEW) THEN
          NNEL = NNEL + 1
          ELEMS(NNEL) = JEL
         END IF
        END DO
       END DO
       NNNEL = NNNEL + NNEL
      END DO

      NNNEL = 4*4*NNNEL

      END IF
!       WRITE(*,*) NNNEL

      END
C
C
C
************************************************************************
      SUBROUTINE  C_Mul_q (C,LDC,COLC,CP,LDCP,COLCP,Q,QP,D,NU)
      IMPLICIT NONE
      REAL*8 C(*),CP(*),Q(*),QP(*),D(*)
      INTEGER LDC(*),COLC(*),LDCP(*),COLCP(*),NU
      INTEGER I,J,K

      DO I=1,NU
       DO J=LDC(I),LDC(I+1)-1
        K = COLC(J)
        D(I) = D(I) - C(J)*Q(K)
       END DO
       DO J=LDCP(I),LDCP(I+1)-1
        K = COLCP(J)
        D(I) = D(I) - CP(J)*QP(K)
       END DO
      END DO

      END
************************************************************************
      SUBROUTINE   Assemble_CMat(DC,KCOLC,KLDC,DM,
     *             B1,B2,B3,KCOLBT,KLDBT,KCOLB,KLDB,KBT,NVT)
************************************************************************
*     Purpose:  calculates the matrix entries for projection of C
*-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER IVT,JVT,NVT,IDOF,JDOF
      REAL*8  DC(*),B1(*),B2(*),B3(*),DM(*)
      INTEGER KCOLC(*),KLDC(*),KCOLBT(*),KLDBT(*),KBT(*)
      INTEGER ILOCC,ILOCB,JLOCB,ILOCKB,JLOCKB,KCOLB(*),KLDB(*)
      REAL*8  DENTRY,DAUX,DIFF
C
      DO IVT=1,NVT
       DO ILOCC=KLDC(IVT),KLDC(IVT+1)-1
        JVT=KCOLC(ILOCC)

        ! Diagonal entry
        DENTRY = 0D0
        IF (JVT.EQ.IVT) THEN
         DO ILOCB=KLDBT(IVT),KLDBT(IVT+1)-1
          IDOF = KCOLBT(ILOCB)
          ILOCKB = KBT(ILOCB)
          DAUX = B1(ILOCKB)**2+B2(ILOCKB)**2+B3(ILOCKB)**2
          DENTRY = DENTRY + DAUX/DM(IDOF)
!           write(*,*) IVT,KCOLB(ILOCKB)
         END DO

        ! Offdiagonal entry
        ELSE
         DO ILOCB=KLDBT(IVT),KLDBT(IVT+1)-1
          IDOF = KCOLBT(ILOCB)
          ILOCKB = KBT(ILOCB)
          DO JLOCB=KLDBT(JVT),KLDBT(JVT+1)-1
           JDOF = KCOLBT(JLOCB)
           IF (IDOF.EQ.JDOF) THEN
            JLOCKB = KBT(JLOCB)
            DAUX = B1(ILOCKB)*B1(JLOCKB)+B2(ILOCKB)*B2(JLOCKB)+
     *             B3(ILOCKB)*B3(JLOCKB)
            DENTRY = DENTRY + DAUX/DM(IDOF)
!             write(*,*) IVT,JVT,KCOLB(ILOCKB),KCOLB(JLOCKB)
!             GOTO 1
           END IF
          END DO
! 1         CONTINUE
         END DO
        END IF
        DC(ILOCC) = DENTRY
       END DO
      END DO
C
C
C
      RETURN
      DIFF = 0d0
      DO IVT=1,NVT
       IDOF=KLDC(IVT+1)-KLDC(IVT)
       IF (IDOF.EQ.27) THEN
        DAUX = 0d0
        DO ILOCC=KLDC(IVT)+1,KLDC(IVT+1)-1
         DAUX = DAUX + DC(ILOCC)
        END DO
!         DC(ILOCC) = DAUX
!         WRITE(*,*) DAUX
        DIFF = MAX(ABS(DAUX),DIFF)
       END IF
      END DO
      WRITE(*,*) DIFF
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   Assemble_CMat_mod(DC,KCOLC,KLDC,DM,
     *             B1,B2,B3,KCOLB,KLDB,
     *             BT1,BT2,BT3,KCOLBT,KLDBT,NVT)
************************************************************************
*     Purpose:  calculates the matrix entries for projection of C
*-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER IVT,JVT,NVT,IDOF,JDOF
      REAL*8  DC(*),B1(*),B2(*),B3(*),BT1(*),BT2(*),BT3(*),DM(*)
      INTEGER KCOLC(*),KLDC(*),KCOLBT(*),KLDBT(*)
      INTEGER ILOCC,ILOCB,JLOCB,ILOCKB,JLOCKB,KCOLB(*),KLDB(*)
      REAL*8  DENTRY,DAUX,DIFF
C
      DO IVT=1,NVT
       DO ILOCC=KLDC(IVT),KLDC(IVT+1)-1
        JVT=KCOLC(ILOCC)

       END DO
      END DO
C
      END
************************************************************************
      SUBROUTINE   Create_BTMat(KCOLB,KLDB,KCOLBT,KLDBT,KBT,NVT)
************************************************************************

      IMPLICIT NONE

      INTEGER NVT
      INTEGER IVT,JVT,IDOF,JDOF,ILOC,JLOC
      INTEGER KCOLB(*),KLDB(*),KCOLBT(*),KLDBT(*),KBT(*)
C
      DO IVT=1,NVT
       DO ILOC=KLDBT(IVT),KLDBT(IVT+1)-1
        IDOF = KCOLBT(ILOC)
        DO JLOC=KLDB(IDOF),KLDB(IDOF+1)-1
         JVT = KCOLB(JLOC)
         IF (JVT.EQ.IVT) THEN
          KBT(ILOC) = JLOC
          GOTO 1
         END IF
        END DO
        write(*,*) "Entry has not been found ...",IVT,IDOF
        STOP
1       CONTINUE
       END DO
      END DO
      END
C
************************************************************************C
!       SUBROUTINE   B_Mul_q(DM,KCOLB,KLDB,KCOLBT,KLDBT,KBT,B1,B2,B3,Q,
!      *             DR1,DR2,DR3,DD,KNPR,NVT,NDOF)
      SUBROUTINE   B_Mul_q(DM,KCOLB,KLDB,B1,B2,B3,Q,
     *             DR1,DR2,DR3,DD,NVT,NDOF)
      IMPLICIT NONE

      REAL*8 B1(*),B2(*),B3(*),DM(*),DD(*),DR1(*),DR2(*),DR3(*),Q(*)
      INTEGER KCOLB(*),KLDB(*)!,KCOLBT(*),KLDBT(*),KBT(*),KNPR(*)
      INTEGER NDOF,NVT
      INTEGER IVT,IDOF,ILOC,ILOCBT
      REAL*8 Def1,Def2,Def3
      REAL*8 DDD(NDOF)

       CALL LAX19(B1,KCOLB,KLDB,NDOF,Q,DR1,1d0,0d0)
       CALL LAX19(B2,KCOLB,KLDB,NDOF,Q,DR2,1d0,0d0)
       CALL LAX19(B3,KCOLB,KLDB,NDOF,Q,DR3,1d0,0d0)

       DO IDOF=1,NDOF
        DDD(IDOF) = DM(IDOF)
       END DO

       CALL E013Sum(DR1)
       CALL E013Sum(DR2)
       CALL E013Sum(DR3)
       CALL E013Sum(DDD)

       DO IDOF=1,NDOF
        DR1(IDOF) = DR1(IDOF)/DDD(IDOF)
        DR2(IDOF) = DR2(IDOF)/DDD(IDOF)
        DR3(IDOF) = DR3(IDOF)/DDD(IDOF)
       END DO

       CALL LTX19(B1,KCOLB,KLDB,NDOF,NVT,DR1,DD,+1d0,0D0)
       CALL LTX19(B2,KCOLB,KLDB,NDOF,NVT,DR2,DD,+1d0,1D0)
       CALL LTX19(B3,KCOLB,KLDB,NDOF,NVT,DR3,DD,+1d0,1D0)

!       DO IDOF=1,NDOF
! !        IF (KNPR(IDOF).NE.1) THEN
!         DR1(IDOF) = 0d0
!         DR2(IDOF) = 0d0
!         DR3(IDOF) = 0d0
!         DO ILOC=KLDB(IDOF),KLDB(IDOF+1)-1
!          IVT = KCOLB(ILOC)
!          DR1(IDOF) = DR1(IDOF) + B1(ILOC)*Q(IVT)
!          DR2(IDOF) = DR2(IDOF) + B2(ILOC)*Q(IVT)
!          DR3(IDOF) = DR3(IDOF) + B3(ILOC)*Q(IVT)
!         END DO
!         DR1(IDOF) = DR1(IDOF)/DM(IDOF)
!         DR2(IDOF) = DR2(IDOF)/DM(IDOF)
!         DR3(IDOF) = DR3(IDOF)/DM(IDOF)
! !        ENDIF
!       END DO
! 
!       DO IVT=1,NVT
!        DD(IVT) = 0d0
!         DO ILOC=KLDBT(IVT),KLDBT(IVT+1)-1
!          ILOCBT = KBT(ILOC)
!          IDOF   = KCOLBT(ILOC)
! !          IF (KNPR(IDOF).NE.1) THEN
!           DD(IVT) = DD(IVT) - B1(ILOCBT)*DR1(IDOF) -
!      *              B2(ILOCBT)*DR2(IDOF) - B3(ILOCBT)*DR3(IDOF)
! !          END IF
!         END DO
!       END DO

      END

C
************************************************************************C
      SUBROUTINE   B_Mul_q_mod(DM,DMC,VA,KCOLB,KLDB,B1,B2,B3,KCOLBT,
     *             KLDBT,BT1,BT2,BT3,Q,DR1,DR2,DR3,AUX1,AUX2,AUX3,
     *             KCOLA,KLDA,DD,NVT,NDOF)
      IMPLICIT NONE

      REAL*4 VA(*)
      REAL*8 B1(*),B2(*),B3(*),DM(*),DD(*),DR1(*),DR2(*),DR3(*),Q(*)
      REAL*8 BT1(*),BT2(*),BT3(*),AUX1(*),AUX2(*),AUX3(*),DMC(*)
      INTEGER KCOLBT(*),KLDBT(*)
      INTEGER KCOLB(*),KLDB(*),KCOLA(*),KLDA(*)
      INTEGER NDOF,NVT
      INTEGER IVT,IDOF,ILOC,ILOCBT
      REAL*8 Def1,Def2,Def3
      REAL*8 DDD(NDOF)

       CALL LAX19(B1,KCOLB,KLDB,NDOF,Q,AUX1,1d0,0d0)
       CALL LAX19(B2,KCOLB,KLDB,NDOF,Q,AUX2,1d0,0d0)
       CALL LAX19(B3,KCOLB,KLDB,NDOF,Q,AUX3,1d0,0d0)

       CALL INVERSE_MASS(DM,DMC,AUX1,AUX2,AUX3,DR1,DR2,DR3,
     *      KCOLA,KLDA,NDOF)

       CALL LAX19(BT1,KCOLBT,KLDBT,NVT,DR1,DD,-1d0,0D0)
       CALL LAX19(BT2,KCOLBT,KLDBT,NVT,DR2,DD,-1d0,1D0)
       CALL LAX19(BT3,KCOLBT,KLDBT,NVT,DR3,DD,-1d0,1D0)

      END

C
************************************************************************C
C
      SUBROUTINE   INVERSE_MASS(DM,DMC,AUX1,AUX2,AUX3,U1,U2,U3,
     *      KCOLA,KLDA,NDOF)
      USE QuadScalar, ONLY : Boundary_QuadScalar_Def

      IMPLICIT NONE

      REAL*8  U1(*),U2(*),U3(*),DM(*)
      REAL*8  AUX1(*),AUX2(*),AUX3(*),DMC(*)
      INTEGER KCOLA(*),KLDA(*)
      INTEGER NDOF
      INTEGER IDOF,I,ILOC,ITE
      REAL*8  DDD(NDOF)
      REAL*8  DD1(NDOF),DD2(NDOF),DD3(NDOF)
      REAL*8  DAUX1,DAUX2,DAUX3

      ! Parallel Lumped mass matrix
      DO IDOF=1,NDOF
       DDD(IDOF) = DM(IDOF) !DMC(KLDA(IDOF)) !
      END DO
      CALL E013Sum(DDD)

      ! Initial guess
      CALL E013Sum(AUX1)
      CALL E013Sum(AUX2)
      CALL E013Sum(AUX3)
      DO IDOF=1,NDOF
       U1(IDOF) = AUX1(IDOF)/DDD(IDOF)
       U2(IDOF) = AUX2(IDOF)/DDD(IDOF)
       U3(IDOF) = AUX3(IDOF)/DDD(IDOF)
      END DO

      CALL Boundary_QuadScalar_Def()
      RETURN

      DO ITE=1,10

       DO I=1,NDOF
        DAUX1=0d0
        DAUX2=0d0
        DAUX3=0d0
        DO ILOC=KLDA(I),KLDA(I+1)-1
         DAUX1 = DAUX1 + DMC(ILOC)*U1(KCOLA(ILOC))
         DAUX2 = DAUX2 + DMC(ILOC)*U2(KCOLA(ILOC))
         DAUX3 = DAUX3 + DMC(ILOC)*U3(KCOLA(ILOC))
        END DO
        DD1(IDOF) = DAUX1 - DM(I)*U1(KCOLA(I))
        DD2(IDOF) = DAUX2 - DM(I)*U2(KCOLA(I))
        DD3(IDOF) = DAUX3 - DM(I)*U3(KCOLA(I))
       END DO

       CALL E013Sum(DD1)
       CALL E013Sum(DD2)
       CALL E013Sum(DD3)

       DO I=1,NDOF
        U1(IDOF) = (AUX1(IDOF)+DD1(IDOF))/DDD(IDOF)
        U2(IDOF) = (AUX2(IDOF)+DD2(IDOF))/DDD(IDOF)
        U3(IDOF) = (AUX3(IDOF)+DD3(IDOF))/DDD(IDOF)
       END DO

       CALL Boundary_QuadScalar_Def()

      END DO


      END
C
************************************************************************C
      SUBROUTINE   BT_Mul_U(KCOLB,KLDB,B1,B2,B3,
     *             DU1,DU2,DU3,DD,NDOF,NVT,DT)
      IMPLICIT NONE

      REAL*8 B1(*),B2(*),B3(*),DD(*),DU1(*),DU2(*),DU3(*)
      INTEGER KCOLB(*),KLDB(*)
      INTEGER NVT,NDOF
      INTEGER IVT
      REAL*8 DT

      CALL LTX19(B1,KCOLB,KLDB,NDOF,NVT,DU1,DD,1d0/DT,0D0)
      CALL LTX19(B2,KCOLB,KLDB,NDOF,NVT,DU2,DD,1d0/DT,1D0)
      CALL LTX19(B3,KCOLB,KLDB,NDOF,NVT,DU3,DD,1d0/DT,1D0)

!       DO IVT=1,NVT
!         DD(IVT) = 0d0
!       END DO

!       DO IVT=1,NVT
!         if (abs(dd(ivt)).gt.1d-5) write(*,*) ivt,DD(IVT)
!       END DO

      END
C
      SUBROUTINE   BT_Mul_U_mod(KCOLB,KLDB,B1,B2,B3,
     *             DU1,DU2,DU3,DD,NDOF,NVT,DT)
      IMPLICIT NONE

      REAL*8 B1(*),B2(*),B3(*),DD(*),DU1(*),DU2(*),DU3(*)
      INTEGER KCOLB(*),KLDB(*)
      INTEGER NVT,NDOF
      INTEGER IVT
      REAL*8 DT

      CALL LAX19(B1,KCOLB,KLDB,NVT,DU1,DD,-1d0/DT,0D0)
      CALL LAX19(B2,KCOLB,KLDB,NVT,DU2,DD,-1d0/DT,1D0)
      CALL LAX19(B3,KCOLB,KLDB,NVT,DU3,DD,-1d0/DT,1D0)

!       DO IVT=1,NVT
!         DD(IVT) = 0d0
!       END DO

!       DO IVT=1,NVT
!         if (abs(dd(ivt)).gt.1d-5) write(*,*) ivt,DD(IVT)
!       END DO

      END
C
************************************************************************C
************************************************************************C
      SUBROUTINE   B_Mul_U(KCOLB,KLDB,B1,B2,B3,Q,
     *             DR1,DR2,DR3,NDOF,A1,A2)
      IMPLICIT NONE

      REAL*8 B1(*),B2(*),B3(*),DR1(*),DR2(*),DR3(*),Q(*)
      INTEGER KCOLB(*),KLDB(*)
      INTEGER NDOF
      INTEGER I,J
      REAL*8 A1,A2,P

      IF (A2.EQ.1d0) THEN
       DO I=1,NDOF
        DO J=KLDB(I),KLDB(I+1)-1
         P = Q(KCOLB(J))
         DR1(I) = DR1(I) + A1*B1(J)*P
         DR2(I) = DR2(I) + A1*B2(J)*P
         DR3(I) = DR3(I) + A1*B3(J)*P
        END DO
       END DO
      END IF

      IF (A2.EQ.0d0) THEN
       DO I=1,NDOF
        DR1(I) = 0d0
        DR2(I) = 0d0
        DR3(I) = 0d0
        DO J=KLDB(I),KLDB(I+1)-1
         P = Q(KCOLB(J))
         DR1(I) = DR1(I) + A1*B1(J)*P
         DR2(I) = DR2(I) + A1*B2(J)*P
         DR3(I) = DR3(I) + A1*B3(J)*P
        END DO
       END DO
      END IF

!       CALL LAX19(B1,KCOLB,KLDB,NDOF,Q,DR1,A1,A2)
!       CALL LAX19(B2,KCOLB,KLDB,NDOF,Q,DR2,A1,A2)
!       CALL LAX19(B3,KCOLB,KLDB,NDOF,Q,DR3,A1,A2)

      END
C
**********************************************************************
      SUBROUTINE SetBCOnBMat(B1,B2,B3,KLDB,NDOF,KNPR)
      IMPLICIT NONE

      REAL*8 B1(*),B2(*),B3(*)
      INTEGER KLDB(*),KNPR(*)
      INTEGER NDOF
      INTEGER IDOF,ILOC

      DO IDOF=1,NDOF
       IF (KNPR(IDOF).EQ.1) THEN
        DO ILOC=KLDB(IDOF),KLDB(IDOF+1)-1
         B1(ILOC) = 0d0
         B2(ILOC) = 0d0
         B3(ILOC) = 0d0
        END DO
       END IF
      END DO

      END
