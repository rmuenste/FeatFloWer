************************************************************************
*     Purpose: Find the midpoint pairs for the PERIODIC BC
*-----------------------------------------------------------------------
      SUBROUTINE PERFIND(KFBD,KABD,KNPR,DCORVG,KVERT,KEBD,
     *           KERIODk,KERIODP,KERIODV,KCROSS,KADJ,IPER,KPER)
************************************************************************
C  
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION KABD(*),KEBD(*),KFBD(*),KERIODP(NEL,2),KADJ(NNAE,*)
      DIMENSION KCROSS(*),KERIODk(*),KERIODV(2,*)
      DIMENSION KNPR(*),DCORVG(3,*),KVERT(8,*)
      REAL*8 P_I(4,3),P_J(4,3), PI(3),PJ(3),PII(3),PJJ(3)
      INTEGER IVT(4),JVT(4)
      LOGICAL BFound
      SAVE
C
      IPER = 0
      JPER = 0
      KPER = 0
C
      DO 100 IABD=1,NABD
C
      IEL=KEBD(IABD)
      IAT=KFBD(IABD)
      INPR=KNPR(NVT+KABD(IABD))
      IF (INPR.NE.3) GOTO 100
C
      KERIODk(IABD) = 0
      JPER = JPER + 1
      IF (IAT.EQ.1) THEN
       IVT(1)=KVERT(1,IEL)
       IVT(2)=KVERT(2,IEL)
       IVT(3)=KVERT(3,IEL)
       IVT(4)=KVERT(4,IEL)
       GOTO 300
      ENDIF
C
      IF (IAT.EQ.2) THEN
       IVT(1)=KVERT(1,IEL)
       IVT(2)=KVERT(2,IEL)
       IVT(3)=KVERT(6,IEL)
       IVT(4)=KVERT(5,IEL)
       GOTO 300
      ENDIF
C
      IF (IAT.EQ.3) THEN
       IVT(1)=KVERT(2,IEL)
       IVT(2)=KVERT(3,IEL)
       IVT(3)=KVERT(7,IEL)
       IVT(4)=KVERT(6,IEL)
       GOTO 300
      ENDIF
C
      IF (IAT.EQ.4) THEN
       IVT(1)=KVERT(3,IEL)
       IVT(2)=KVERT(4,IEL)
       IVT(3)=KVERT(8,IEL)
       IVT(4)=KVERT(7,IEL)
       GOTO 300
      ENDIF
C
      IF (IAT.EQ.5) THEN
       IVT(1)=KVERT(4,IEL)
       IVT(2)=KVERT(1,IEL)
       IVT(3)=KVERT(5,IEL)
       IVT(4)=KVERT(8,IEL)
       GOTO 300
      ENDIF
C
      IF (IAT.EQ.6) THEN
       IVT(1)=KVERT(5,IEL)
       IVT(2)=KVERT(6,IEL)
       IVT(3)=KVERT(7,IEL)
       IVT(4)=KVERT(8,IEL)
       GOTO 300
      ENDIF
C
C
300   P_I(1,1)=DCORVG(1,IVT(1))
      P_I(1,2)=DCORVG(2,IVT(1))
      P_I(1,3)=DCORVG(3,IVT(1))
      P_I(2,1)=DCORVG(1,IVT(2))
      P_I(2,2)=DCORVG(2,IVT(2))
      P_I(2,3)=DCORVG(3,IVT(2))
      P_I(3,1)=DCORVG(1,IVT(3))
      P_I(3,2)=DCORVG(2,IVT(3))
      P_I(3,3)=DCORVG(3,IVT(3))
      P_I(4,1)=DCORVG(1,IVT(4))
      P_I(4,2)=DCORVG(2,IVT(4))
      P_I(4,3)=DCORVG(3,IVT(4))
C
      PI(1)=(P_I(1,1)+P_I(2,1)+P_I(3,1)+P_I(4,1))*0.25D0
      PI(2)=(P_I(1,2)+P_I(2,2)+P_I(3,2)+P_I(4,2))*0.25D0
      PI(3)=(P_I(1,3)+P_I(2,3)+P_I(3,3)+P_I(4,3))*0.25D0
C
      DO 101 JABD=1,NABD
C
      IF (IABD.EQ.JABD) GOTO 101
C
      JEL=KEBD(JABD) 
      JAT=KFBD(JABD)
      JNPR=KNPR(NVT+KABD(JABD))
      IF (JNPR.NE.3) GOTO 101
C
      IF (JAT.EQ.1) THEN
       JVT(1)=KVERT(1,JEL)
       JVT(2)=KVERT(2,JEL)
       JVT(3)=KVERT(3,JEL)
       JVT(4)=KVERT(4,JEL)
       GOTO 301
      ENDIF
C
      IF (JAT.EQ.2) THEN
       JVT(1)=KVERT(1,JEL)
       JVT(2)=KVERT(2,JEL)
       JVT(3)=KVERT(6,JEL)
       JVT(4)=KVERT(5,JEL)
       GOTO 301
      ENDIF
C
      IF (JAT.EQ.3) THEN
       JVT(1)=KVERT(2,JEL)
       JVT(2)=KVERT(3,JEL)
       JVT(3)=KVERT(7,JEL)
       JVT(4)=KVERT(6,JEL)
       GOTO 301
      ENDIF
C
      IF (JAT.EQ.4) THEN
       JVT(1)=KVERT(3,JEL)
       JVT(2)=KVERT(4,JEL)
       JVT(3)=KVERT(8,JEL)
       JVT(4)=KVERT(7,JEL)
       GOTO 301
      ENDIF
C
      IF (JAT.EQ.5) THEN
       JVT(1)=KVERT(4,JEL)
       JVT(2)=KVERT(1,JEL)
       JVT(3)=KVERT(5,JEL)
       JVT(4)=KVERT(8,JEL)
       GOTO 301
      ENDIF
C
      IF (JAT.EQ.6) THEN
       JVT(1)=KVERT(5,JEL)
       JVT(2)=KVERT(6,JEL)
       JVT(3)=KVERT(7,JEL)
       JVT(4)=KVERT(8,JEL)
       GOTO 301
      ENDIF
C
C
301   P_J(1,1)=DCORVG(1,JVT(1))
      P_J(1,2)=DCORVG(2,JVT(1))
      P_J(1,3)=DCORVG(3,JVT(1))
      P_J(2,1)=DCORVG(1,JVT(2))
      P_J(2,2)=DCORVG(2,JVT(2))
      P_J(2,3)=DCORVG(3,JVT(2))
      P_J(3,1)=DCORVG(1,JVT(3))
      P_J(3,2)=DCORVG(2,JVT(3))
      P_J(3,3)=DCORVG(3,JVT(3))
      P_J(4,1)=DCORVG(1,JVT(4))
      P_J(4,2)=DCORVG(2,JVT(4))
      P_J(4,3)=DCORVG(3,JVT(4))
C
      PJ(1)=(P_J(1,1)+P_J(2,1)+P_J(3,1)+P_J(4,1))*0.25D0
      PJ(2)=(P_J(1,2)+P_J(2,2)+P_J(3,2)+P_J(4,2))*0.25D0
      PJ(3)=(P_J(1,3)+P_J(2,3)+P_J(3,3)+P_J(4,3))*0.25D0
C
      CALL PERIODa(IFound,PI,PJ)
C
      IF (ABS(IFound).EQ.1) THEN
C      Assignement of the midpoints
       KERIODk(IABD) = (+IFound)*KABD(JABD)
       KERIODk(JABD) = (-IFound)*KABD(IABD)
       KNPR(NVT+KABD(IABD)) = 0
       KNPR(NVT+KABD(JABD)) = 0
       IF (KERIODP(IEL,1).EQ.0) THEN
        KERIODP(IEL,1) = (+IFound)*JEL
        KERIODP(IEL,2) = IAT
       ELSE
        IF (KERIODP(IEL,2).GT.IAT) THEN
         IU = KERIODP(IEL,1)
         KERIODP(IEL,1) = (+IFound)*JEL
         KERIODP(IEL,2) = IU
        ELSE
         KERIODP(IEL,2) = (+IFound)*JEL
        END IF
       END IF
       IF (KERIODP(JEL,1).EQ.0) THEN
        KERIODP(JEL,1) = (-IFound)*IEL
        KERIODP(JEL,2) = JAT
       ELSE
        IF (KERIODP(JEL,2).GT.JAT) THEN
         IU = KERIODP(JEL,1)
         KERIODP(JEL,1) = (-IFound)*IEL
         KERIODP(JEL,2) = IU
        ELSE
         KERIODP(JEL,2) = (-IFound)*IEL
        END IF
       END IF
       IPER = IPER + 1
C      Assignement of the vertices
       DO 400 KA=1,4
        DO 401 KB =1,4
         PII(1)=P_I(KA,1)
         PII(2)=P_I(KA,2)
         PII(3)=P_I(KA,3)
         PJJ(1)=P_J(KB,1)
         PJJ(2)=P_J(KB,2)
         PJJ(3)=P_J(KB,3)
         CALL PERIODa(JFound,PII,PJJ)
         IF (ABS(JFound).EQ.1) THEN
          DO KC=1,KPER
           IF (IVT(KA).EQ.KERIODV(1,KC).OR.
     *         IVT(KA).EQ.KERIODV(2,KC)) GOTO 400
          END DO
          IF (JFound.GT.0) THEN
           KERIODV(1,KPER+1)=IVT(KA)
           KERIODV(2,KPER+1)=JVT(KB)
          ELSE
           KERIODV(1,KPER+1)=JVT(KB)
           KERIODV(2,KPER+1)=IVT(KA)
          END IF
          KPER = KPER + 1
          KNPR(IVT(1)) = 0
          KNPR(IVT(2)) = 0
          KNPR(IVT(3)) = 0
          KNPR(IVT(4)) = 0
          KNPR(JVT(1)) = 0
          KNPR(JVT(2)) = 0
          KNPR(JVT(3)) = 0
          KNPR(JVT(4)) = 0
          GOTO 400
         ELSE 
          GOTO 401
         END IF
         WRITE(*,*) "Problems with vertice assignement ... "
         STOP
401     CONTINUE
400    CONTINUE
C
       GOTO 100
      END IF

101   CONTINUE
C
100   CONTINUE

       IF (JPER.NE.IPER) THEN
        WRITE(*,'(A57,I8,I8)') "Periodic vertice assignement wasn't succ
     *esful, #Mid 1,2: ", IPER,JPER
        STOP
       END IF
C
      DO ICOL=1,NAT
       KCROSS(ICOL) = 0
      END DO
C
      DO IABD=1,NABD
       ICOL = KABD(IABD)
       KCROSS(ICOL) = KERIODk(IABD)
      END DO
C
      END

C
C
C
************************************************************************
C     Rebuilds the matrix structure and adds new enties in the matrix
*-----------------------------------------------------------------------
      SUBROUTINE PERMATA(KCOL_O,KLD_O,KERIO,KCOL_N,KLD_N,
     *                   KCROSS,JCOL_N)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)

      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      DIMENSION KCOL_O(*),KLD_O(*),KERIO(*),KCOL_N(*),KLD_N(*)
      DIMENSION LW(20),KCROSS(*)
      DATA BLOOK/.FALSE./

      JCOL_N = 1
      JJJ = 0
      III = 0
      KKK = 0
      IIK = 0
      JJK = 0
C
      DO ICOL=1,NAT
C
       KLD_N(ICOL) = JCOL_N
C
       IF     (KCROSS(ICOL).EQ.0) THEN

        ILW = 0
        KK=KKK
        KCOL_N(JCOL_N) = KCOL_O(KLD_O(ICOL))
        JCOL_N = JCOL_N + 1
        DO JCOL=KLD_O(ICOL)+1,KLD_O(ICOL+1)-1
         ILW = ILW + 1
         LW(ILW) = KCOL_O(JCOL)
         IF (KCROSS(LW(ILW)).NE.0) THEN
          KKK = KKK + 1
          ILW = ILW + 1
          LW(ILW) = ABS(KCROSS(LW(ILW-1)))
         END IF
        END DO
!        WRITE(*,*) ICOL," - ", ILW, KLD_O(ICOL+1)-1-KLD_O(ICOL)-1+1
        CALL SORTING(LW,ILW)
        DO JCOL=1,ILW
         KCOL_N(JCOL_N) = LW(JCOL)
         JCOL_N = JCOL_N + 1
        END DO

        IF (BLOOK) THEN
         IF (KK.NE.KKK) THEN
          WRITE(*,'(A3,12I7)') 
     *    "N ",(KCOL_N(J),J=KLD_N(ICOL),KLD_N(ICOL)+ILW)
          WRITE(*,'(A3,12I7)') 
     *    "N ",(KCOL_O(J),J=KLD_O(ICOL),KLD_O(ICOL+1)-1)
          WRITE(*,*) "--------------------------------------------"
         END IF
        END IF

       ELSEIF (KCROSS(ICOL).LT.0) THEN

        III = III + 1
        ILW = 0
        KCOL_N(JCOL_N) = KCOL_O(KLD_O(ICOL))
        JCOL_N = JCOL_N + 1
        DO JCOL=KLD_O(ICOL)+1,KLD_O(ICOL+1)-1
         ILW = ILW + 1
         LW(ILW) = KCOL_O(JCOL)
        END DO
        DO INEW=0,5
         ILW = ILW + 1
         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+INEW)
        END DO
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+0)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+1)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+2)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+3)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+4)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+5)
        CALL SORTING(LW,ILW)
        DO JCOL=1,ILW
         KCOL_N(JCOL_N) = LW(JCOL)
         JCOL_N = JCOL_N + 1
        END DO
        IIK = IIK + ILW

        IF (BLOOK) THEN
         WRITE(*,'(A3,12I7)') 
     *   "- ",(KCOL_N(J),J=KLD_N(ICOL),KLD_N(ICOL)+ILW)
         WRITE(*,'(A3,12I7)') 
     *   "- ",(KCOL_O(J),J=KLD_O(ICOL),KLD_O(ICOL+1)-1)
         WRITE(*,*) "--------------------------------------------"
        END IF
      ELSEIF (KCROSS(ICOL).GT.0) THEN

       JJJ = JJJ + 1
       ILW = 0
        KCOL_N(JCOL_N) = KCOL_O(KLD_O(ICOL))
        JCOL_N = JCOL_N + 1
        DO JCOL=KLD_O(ICOL)+1,KLD_O(ICOL+1)-1
         ILW = ILW + 1
         LW(ILW) = KCOL_O(JCOL)
        END DO
        DO INEW=0,5
         ILW = ILW + 1
         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+INEW)
        END DO
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+0)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+1)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+2)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+3)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+4)
!         ILW = ILW + 1
!         LW(ILW) = KCOL_O(KLD_O(ABS(KCROSS(ICOL)))+5)
        CALL SORTING(LW,ILW)
        DO JCOL=1,ILW
         KCOL_N(JCOL_N) = LW(JCOL)
         JCOL_N = JCOL_N + 1
        END DO
        JJK = JJK + ILW

        IF (BLOOK) THEN
         WRITE(*,'(A3,12I7)') 
     *   "+ ",(KCOL_N(J),J=KLD_N(ICOL),KLD_N(ICOL)+ILW)
         WRITE(*,'(A3,12I7)') 
     *   "+ ",(KCOL_O(J),J=KLD_O(ICOL),KLD_O(ICOL+1)-1)
         WRITE(*,*) "--------------------------------------------"
        END IF

       END IF

      END DO

      KLD_N(NAT+1) = JCOL_N

      END
C
C
C
************************************************************************
C     Inserts the entires into the global matrix according to 
C     the periodic BC's
*-----------------------------------------------------------------------
      SUBROUTINE BDRPER(D1,D2,D3,VA,KLD,KCOL,KABD,NABD,KERIOA,IDEF)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*4 VA,A1,A2
      DIMENSION D1(*),D2(*),D3(*)
      DIMENSION KLD(*),KCOL(*),KABD(*),KERIOA(*)
      DIMENSION VA(*),A1(12),A2(12)
      DIMENSION L1(12),L2(12)
C
      DO 1 IAT=1,NABD
C
       IABD=KABD(IAT)
C
       IF (KERIOA(IAT).GT.0) THEN
        IPER = IABD
        JPER = ABS(KERIOA(IAT))
        JC = 0
        DO IC = KLD(IPER),KLD(IPER+1)-1
         JC = JC + 1
         A1(JC) = VA(IC)
         L1(JC) = KCOL(IC)
        END DO
        JC = 0
        DO IC = KLD(JPER),KLD(JPER+1)-1
         JC = JC + 1
         A2(JC) = VA(IC)
         L2(JC) = KCOL(IC)
        END DO
        CALL PERSETA(L1,L2,A1,A2,IDEF)
        JC = 0
        DO IC = KLD(IPER),KLD(IPER+1)-1
         JC = JC + 1
         VA(IC) = A1(JC)
        END DO
        JC = 0
        DO IC = KLD(JPER),KLD(JPER+1)-1
         JC = JC + 1
         VA(IC) = A2(JC)
        END DO
        IF (IDEF.GE.1) THEN
         CALL PERSETD(D1(IPER),D2(IPER),D3(IPER),
     *                D1(JPER),D2(JPER),D3(JPER))
        END IF
       END IF
C
1     CONTINUE
      END
C
C
C
************************************************************************
C     Update of U = U~ - k B P according to periodic BC's
C  Division of the correction with ML and update the velocities
*-----------------------------------------------------------------------
      SUBROUTINE DIVPERM(DM,DRHO,D1,D2,D3,DU1,DU2,DU3,KERIOA,KCROSS,
     *           KABD,NAT,NABD,ILEV)
************************************************************************C
      USE PP3D_MPI
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DM(*),D1(*),D2(*),D3(*),DU1(*),DU2(*),DU3(*)
      DIMENSION KERIOA(*),KCROSS(*),KABD(*)
      REAL*8    DRHO(*),RHODM(NAT)
C
      DO I=1,NAT
       RHODM(I) = DM(I)*DRHO(I)
      END DO
C
      DO IUU=1,mg_mpi(ILEV)%NeighNum                                  ! PARALLEL
       mg_mpi(ILEV)%parST(IUU)%i=1                                    ! PARALLEL
      END DO                                                          ! PARALLEL

      CALL CommBMMul(RHODM,ILEV)                                      ! PARALLEL
C
      DO 2 IAT=1,NAT
       JAT =KCROSS(IAT)
       IF      (JAT.NE.0) THEN
         IF (JAT.LT.0) GOTO 2
         DU1(IAT)=DU1(IAT)-D1(IAT)/(RHODM(IAT)+RHODM(JAT))            ! PERIODIC
         DU2(IAT)=DU2(IAT)-D2(IAT)/(RHODM(IAT)+RHODM(JAT))            ! PERIODIC
         DU3(IAT)=DU3(IAT)-D3(IAT)/(RHODM(IAT)+RHODM(JAT))            ! PERIODIC
       ELSE
        DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
         IUU = mg_mpi(ILEV)%parST(pID)%i                              ! PARALLEL
         IF (mg_mpi(ILEV)%parST(pID)%FaceLink(1,IUU).EQ.IAT) THEN     ! PARALLEL
          DM_IUU = mg_mpi(ILEV)%parST(pID)%RDVect(IUU)                ! PARALLEL
          DU1(IAT)=DU1(IAT)-D1(IAT)/(RHODM(IAT)+DM_IUU)               ! PARALLEL
          DU2(IAT)=DU2(IAT)-D2(IAT)/(RHODM(IAT)+DM_IUU)               ! PARALLEL
          DU3(IAT)=DU3(IAT)-D3(IAT)/(RHODM(IAT)+DM_IUU)               ! PARALLEL
          mg_mpi(ILEV)%parST(pID)%i=IUU+1                             ! PARALLEL
          GOTO 66                                                     ! PARALLEL
         END IF                                                       ! PARALLEL
        END DO
        DU1(IAT)=DU1(IAT)-D1(IAT)/RHODM(IAT)                          ! REGULAR
        DU2(IAT)=DU2(IAT)-D2(IAT)/RHODM(IAT)                          ! REGULAR
        DU3(IAT)=DU3(IAT)-D3(IAT)/RHODM(IAT)                          ! REGULAR
66      CONTINUE
       END IF
2     CONTINUE
C
      DO 3 IABD=1,NABD
       IAT=KABD(IABD)
       JAT=KCROSS(IAT)
       IF (JAT.LT.0) THEN
         DU1(IAT)=DU1(ABS(JAT))
         DU2(IAT)=DU2(ABS(JAT))
         DU3(IAT)=DU3(ABS(JAT))
       END IF
C
3     CONTINUE
C
      DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
       IF (mg_mpi(ILEV)%parST(pID)%i-1.NE.                          ! PARALLEL
     *     mg_mpi(ILEV)%parST(pID)%Num) WRITE(*,*)                  ! PARALLEL
     *     "Problem with",myid,"--",mg_mpi(ILEV)%parST(pID)%Neigh   ! PARALLEL
      END DO

      END
C
C
C
************************************************************************
C     Sorts the vector in increasing order of entries
*-----------------------------------------------------------------------
      SUBROUTINE SORTING(LW,N)
************************************************************************C
      INTEGER LW(20),LWA
      INTEGER I,J,N,JJJ
C
!      WRITE(*,*) (LW(JJJ),JJJ=1,N),N
!      IF (N.GT.13) WRITE(*,*) (LW(JJJ),JJJ=1,N)

      DO I=2,N
       DO J=N,I,-1
        IF (LW(J).LT.LW(J-1)) THEN
         LWA     = LW(J)
         LW(J)   = LW(J-1)
         LW(J-1) = LWA
        END IF
       END DO
!      WRITE(*,*) (LW(JJJ),JJJ=1,N),N
      END DO
C
!      IF (N.GT.13) WRITE(*,*) (LW(JJJ),JJJ=1,N)
!      IF (N.GT.13) WRITE(*,*) ("---",JJJ=1,N)
      END
C
C
C
************************************************************************
C     Sets the matrix entries according to the periodic BC
*-----------------------------------------------------------------------
      SUBROUTINE PERSETA(LCOL1,LCOL2,A1,A2,IW)
************************************************************************C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NSIZE = 12)
      INTEGER LCOL1(12),LCOL2(12)
      REAL*4 A1,A2
      DIMENSION A1(12),A2(12)
      DATA BLOOK/.FALSE./
C
      IF (A1(1).EQ.0D0.OR.A2(1).EQ.0E0) RETURN
C
      IF (LCOL2(1).EQ.21101.AND.BLOOK) THEN
       WRITE(*,3) "A1 ", (A1(I),I=1,NSIZE)
       WRITE(*,4) "L1 ", (LCOL1(I),I=1,NSIZE)
       WRITE(*,3) "A2 ", (A2(I),I=1,NSIZE)
       WRITE(*,4) "L2 ", (LCOL2(I),I=1,NSIZE)
      END IF
C
      DO 1 I=2,NSIZE
       DO 2 J=2,NSIZE
        IF (LCOL1(I).EQ.LCOL2(J)) THEN
         A1(I) = A1(I) + A2(J)
         GOTO 1
        END IF
 2     CONTINUE
 1    CONTINUE
      A1(1) = A1(1) + A2(1)
C
      A2(1) = +1D0
      DO I=2,NSIZE
       IF (LCOL1(1).NE.LCOL2(I)) THEN
        A2(I) =  0D0
       ELSE
        A2(I) = -1D0
       END IF
      END DO
C
      IF (LCOL2(1).EQ.21101.AND.BLOOK) THEN
       WRITE(*,3) "A1 ", (A1(I),I=1,NSIZE)
       WRITE(*,4) "L1 ", (LCOL1(I),I=1,NSIZE)
       WRITE(*,3) "A2 ", (A2(I),I=1,NSIZE)
       WRITE(*,4) "L2 ", (LCOL2(I),I=1,NSIZE)
       WRITE(*,*) "****************", IW
      END IF
3     FORMAT(A3,12(D9.2))
4     FORMAT(A3,12(I9))
      END
C
C
C
************************************************************************
C     Sets the right hand side according to the periodic BC
*-----------------------------------------------------------------------
      SUBROUTINE PERSETD(DI1,DI2,DI3,DJ1,DJ2,DJ3)
************************************************************************C
      REAL*8 DI1,DI2,DI3,DJ1,DJ2,DJ3
C
      DI1 = DI1 + DJ1
      DI2 = DI2 + DJ2
      DI3 = DI3 + DJ3
      DJ1 = 0D0
      DJ2 = 0D0
      DJ3 = 0D0
C
      END
C
C
C
************************************************************************
C     Modifies the solution vector according to periodic BC 
C               Needed just for GMV output!!
*-----------------------------------------------------------------------
      SUBROUTINE PERGMV(VSOL,KERIOV,NKER)
************************************************************************C
      REAL*4 V,VSOL
      INTEGER KERIOV
      DIMENSION VSOL(*),KERIOV(2,*)
C
      DO IKER=1,NKER
       V=0.5E0*(VSOL(KERIOV(1,IKER))+VSOL(KERIOV(2,IKER)))
       VSOL(KERIOV(1,IKER)) = V
       VSOL(KERIOV(2,IKER)) = V
      END DO
C
      END
C
C
C
************************************************************************
