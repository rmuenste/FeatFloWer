C
C
************************************************************************
      SUBROUTINE MATDEF_AFC(U1,U2,U3,D1,D2,D3,DRHO,F1,F2,F3,A,
     *           DK,DS,KCOLA,KLDA,NU,KVERT,KAREA,KEDGE,KINT,DCORVG,
     *           DNUT,IELT,DML,DMASS,ISEP,IAUX,INOD,JNOD,AEDGE,IDEF)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*4 A
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
C
      DIMENSION U1(*),U2(*),U3(*),D1(*),D2(*),D3(*),F1(*),F2(*),F3(*)
      DIMENSION A(*),DK(*),DS(*),DML(*),DMASS(*),DNUT(*),DRHO(*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KCOLA(*),KLDA(*),KINT(*),DCORVG(NNDIM,*)
      DIMENSION ISEP(*),IAUX(*),INOD(*),JNOD(*),AEDGE(*)
C
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
C
      EXTERNAL E030,E031,EM31,EM30
C
      SAVE 
C
C     Add the contribution of the mass matrix to the defect
      IF (IDEF.NE.0) THEN
        SIGML=DBLE(SIGN(1,IDEF))
        DO I=1,NU
          DAUX=SIGML*DRHO(I)*DML(I)
          D1(I)=D1(I)+DAUX*U1(I)
          D2(I)=D2(I)+DAUX*U2(I)
          D3(I)=D3(I)+DAUX*U3(I)
        ENDDO
        IEDGE=0
      ENDIF
C
      NA=KLDA(NU+1)-1
C
      IF (IELT.EQ.0) CALL CONVDG(U1,U2,U3,DRHO,DK,NA,KCOLA,KLDA,
     *    KVERT,KAREA,KEDGE,KINT,DCORVG,E031)
      IF (IELT.EQ.1) CALL CONVDG(U1,U2,U3,DRHO,DK,NA,KCOLA,KLDA,
     *    KVERT,KAREA,KEDGE,KINT,DCORVG,E030)
      IF (IELT.EQ.2) CALL CONVNP(U1,U2,U3,DRHO,DK,NA,KCOLA,KLDA,
     *    KVERT,KAREA,KEDGE,KINT,DCORVG,EM31)
      IF (IELT.EQ.3) CALL CONVNP(U1,U2,U3,DRHO,DK,NA,KCOLA,KLDA,
     *    KVERT,KAREA,KEDGE,KINT,DCORVG,EM30)
C
C     Find the location of the first entry a_ij such that j > i
      DO I=1,NU
        ILOC=KLDA(I)
        IF (IDEF.LE.0) A(ILOC)=REAL(DRHO(I)*DML(I))
 25     ILOC=ILOC+1 
        IF (KCOLA(ILOC).LT.I.AND.ILOC.LT.KLDA(I+1)) GOTO 25
        ISEP(I)=ILOC-1
      END DO
      CALL LCP3(ISEP,IAUX,NU)
C
      DO I=1,NU
C
        II_LOC=KLDA(I)
C
        DO IJ_LOC=KLDA(I)+1,ISEP(I)
C
C         Position of entries in global matrices
          J=KCOLA(IJ_LOC)
          IAUX(J)=IAUX(J)+1
          JI_LOC=IAUX(J)
          JJ_LOC=KLDA(J)
!           IF (I.NE.KCOLA(JI_LOC).OR.J.NE.KCOLA(IJ_LOC)) THEN
!            IAUX(J)=IAUX(J)-1
!            IF (IDEF.LE.0) A(IJ_LOC)=0e0
!            GOTO 30
!           END IF
C
C         Discrete diffusion coefficients
          DK_IJ=DK(IJ_LOC)
          DK_JI=DK(JI_LOC)
          D_IJ=MAX(-DK_IJ,0D0,-DK_JI)
C
C         Elimination of negative entries
          DK(II_LOC)=DK(II_LOC)-D_IJ
          DK(IJ_LOC)=DK(IJ_LOC)+D_IJ
          DK(JJ_LOC)=DK(JJ_LOC)-D_IJ
          DK(JI_LOC)=DK(JI_LOC)+D_IJ
C
C         Perform global matrix assembly
          IF (IDEF.LE.0) THEN
           DAUX=DK(IJ_LOC)+DS(IJ_LOC)
           A(IJ_LOC)=-REAL(THSTEP*DAUX)
           A(II_LOC)=A(II_LOC)-A(IJ_LOC)
           DAUX=DK(JI_LOC)+DS(JI_LOC)
           A(JI_LOC)=-REAL(THSTEP*DAUX)
           A(JJ_LOC)=A(JJ_LOC)-A(JI_LOC)
          ENDIF
C
C         Determine which node is located upwind 
          IF (IDEF.NE.0) THEN
           IEDGE=IEDGE+1
           SIGML=DBLE(SIGN(1,IDEF))
           DAUX=SIGML*DMASS(IJ_LOC)
           F1(IEDGE)=DAUX*(U1(J)-U1(I))
           F2(IEDGE)=DAUX*(U2(J)-U2(I))
           F3(IEDGE)=DAUX*(U3(J)-U3(I))
           IF (DK_JI.GT.DK_IJ) THEN
            INOD(IEDGE)=I
            JNOD(IEDGE)=J
            AEDGE(IEDGE)=MIN(D_IJ,DK_JI+D_IJ)
           ELSE
            INOD(IEDGE)=J
            JNOD(IEDGE)=I
            AEDGE(IEDGE)=MIN(D_IJ,DK_IJ+D_IJ)
           ENDIF
          ENDIF
C
30      CONTINUE
C
        ENDDO
C
      ENDDO
C
C     Augment the defect if required
      IF (IDEF.NE.0) THEN
        CALL DEFTVD(DK,KCOLA,KLDA,NU,U1,D1,F1,DS,
     *              IEDGE,INOD,JNOD,AEDGE)
        CALL DEFTVD(DK,KCOLA,KLDA,NU,U2,D2,F2,DS,
     *              IEDGE,INOD,JNOD,AEDGE)
        CALL DEFTVD(DK,KCOLA,KLDA,NU,U3,D3,F3,DS,
     *              IEDGE,INOD,JNOD,AEDGE)
      ENDIF
C
99999 END
C
C 
************************************************************************
      SUBROUTINE MATDEF(U1,U2,U3,D1,D2,D3,DRHO,A,DK,DS,
     *           KCOLA,KLDA,NU,KVERT,KAREA,KEDGE,KINT,
     *           DCORVG,DNUT,IELT,DML,DMASS,IDEF)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*4 A
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
C
      DIMENSION U1(*),U2(*),U3(*),D1(*),D2(*),D3(*),DRHO(*)
      DIMENSION A(*),DK(*),DS(*),DML(*),DMASS(*),DNUT(*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
C
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
C
      EXTERNAL E030,E031,EM31,EM30
C
      SAVE 
C
      NA=KLDA(NU+1)-1
C
!       IF (IELT.EQ.0) CALL CONVDG(U1,U2,U3,DK,NA,KCOLA,KLDA,
!      *    KVERT,KAREA,KEDGE,KINT,DCORVG,E031)
!       IF (IELT.EQ.1) CALL CONVDG(U1,U2,U3,DK,NA,KCOLA,KLDA,
!      *    KVERT,KAREA,KEDGE,KINT,DCORVG,E030)
!       IF (IELT.EQ.2) CALL CONVNP(U1,U2,U3,DK,NA,KCOLA,KLDA,
!      *    KVERT,KAREA,KEDGE,KINT,DCORVG,EM31)
!       IF (IELT.EQ.3) CALL CONVNP(U1,U2,U3,DK,NA,KCOLA,KLDA,
!      *    KVERT,KAREA,KEDGE,KINT,DCORVG,EM30)
C
C     Assembly of the upwinded convection term
      CALL GUPWD(U1,U2,U3,DRHO,DK,NA,KCOLA,KLDA,KVERT,KAREA,DCORVG)
C
C     Perform global matrix assembly
      IF (IDEF.LE.0) THEN
        DO I=1,NU
C
          II_LOC=KLDA(I)
          DAUX=1.0d0*DK(II_LOC)+DS(II_LOC)
          A(II_LOC)=REAL(DRHO(I)*DML(I))-REAL(THSTEP*DAUX)
C
          DO IJ_LOC=KLDA(I)+1,KLDA(I+1)-1
            DAUX=1.0d0*DK(IJ_LOC)+DS(IJ_LOC)
            A(IJ_LOC)=-REAL(THSTEP*DAUX)
          ENDDO
        ENDDO
      END IF
C
C     Augment the defect if required
      IF (IDEF.NE.0) THEN
        SIGML=DBLE(SIGN(1,IDEF))
        DO I=1,NU
         DAUX=SIGML*DRHO(I)*DML(I)
         DD1=(1.0d0*DK(KLDA(I))+DS(KLDA(I)))*U1(I)
         DD2=(1.0d0*DK(KLDA(I))+DS(KLDA(I)))*U2(I)
         DD3=(1.0d0*DK(KLDA(I))+DS(KLDA(I)))*U3(I)
         DO ILOC=KLDA(I)+1,KLDA(I+1)-1
           J=KCOLA(ILOC)
           DD1=DD1+(1.0d0*DK(ILOC)+DS(ILOC))*U1(J)
           DD2=DD2+(1.0d0*DK(ILOC)+DS(ILOC))*U2(J)
           DD3=DD3+(1.0d0*DK(ILOC)+DS(ILOC))*U3(J)
         ENDDO
         D1(I)=D1(I)+DAUX*U1(I)+THSTEP*DD1
         D2(I)=D2(I)+DAUX*U2(I)+THSTEP*DD2
         D3(I)=D3(I)+DAUX*U3(I)+THSTEP*DD3
        ENDDO
      ENDIF
C
99999 END
C
C 
      SUBROUTINE DEFTVD(DK,KCOLA,KLDA,NU,U,D,FLUX,DS,
     *                  NEDGE,INOD,JNOD,AEDGE)
C
      USE PP3D_MPI, ONLY: CommSum
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DK(*),DS(*),KCOLA(*),KLDA(*),U(*),D(*)
      DIMENSION FLUX(*),INOD(*),JNOD(*),AEDGE(*)
C
      DIMENSION PP(NU),PM(NU),QP(NU),QM(NU)
      PARAMETER (DEPS=1D-15)
C
      PARAMETER (NNLEV=9)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
C
C
      DO I=1,NU
C
       PP(I)=0D0
       PM(I)=0D0
       QP(I)=0D0
       QM(I)=0D0
C
       DD=(DK(KLDA(I))+DS(KLDA(I)))*U(I)
       DO ILOC=KLDA(I)+1,KLDA(I+1)-1
         J=KCOLA(ILOC)
         DD=DD+(DK(ILOC)+DS(ILOC))*U(J)
       ENDDO
       D(I)=D(I)+THSTEP*DD
C
      ENDDO
C
      DO IEDGE=1,NEDGE
C
       I=INOD(IEDGE)
       J=JNOD(IEDGE)
C
       DAUX=AEDGE(IEDGE)*(U(I)-U(J))
C 
       PP(I)=PP(I)+MAX(0D0,DAUX)
       PM(I)=PM(I)+MIN(0D0,DAUX)
C
       QP(I)=QP(I)+MAX(0D0,-DAUX)
       QM(I)=QM(I)+MIN(0D0,-DAUX)
C
       QP(J)=QP(J)+MAX(0D0,DAUX)
       QM(J)=QM(J)+MIN(0D0,DAUX)
C
      ENDDO
C
      CALL CommSum(PP,ILEV)          ! PARALLEL
      CALL CommSum(PM,ILEV)          ! PARALLEL
      CALL CommSum(QP,ILEV)          ! PARALLEL
      CALL CommSum(QM,ILEV)          ! PARALLEL

      DO I=1,NU
        IF (PP(I).GT.QP(I)+DEPS) THEN
          PP(I)=QP(I)/PP(I)
        ELSE
          PP(I)=1D0
        ENDIF
        IF (PM(I).LT.QM(I)-DEPS) THEN
          PM(I)=QM(I)/PM(I)
        ELSE
          PM(I)=1D0
        ENDIF
      ENDDO
C
      DO IEDGE=1,NEDGE
C
       I=INOD(IEDGE)
       J=JNOD(IEDGE)
C
       DAUX=THSTEP*AEDGE(IEDGE)*(U(I)-U(J))
C 
       IF (I.GT.J) THEN
         FLUX(IEDGE)=FLUX(IEDGE)+DAUX
       ELSE
         FLUX(IEDGE)=FLUX(IEDGE)-DAUX
       END IF
C
       IF (DAUX.GT.0) THEN
        DAUX=PP(I)*DAUX
       ELSE
        DAUX=PM(I)*DAUX
       ENDIF
C
       IF (I.GT.J) THEN
         FLUX(IEDGE)=FLUX(IEDGE)-DAUX
       ELSE
         FLUX(IEDGE)=FLUX(IEDGE)+DAUX
       END IF
C
      D(I)=D(I)+DAUX
      D(J)=D(J)-DAUX
C
      ENDDO
C
      END
C

      SUBROUTINE DEFFCT(DU,DML,KNPR,NU,NEDGE,FLUX,FLUX0,INOD,JNOD)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DU(*),DML(*),KNPR(*)
      DIMENSION FLUX(*),FLUX0(*),INOD(*),JNOD(*)
C
      DIMENSION PP(NU),PM(NU),QP(NU),QM(NU),D(NU)
      PARAMETER (DEPS=1D-15)
C
C
      DO I=1,NU
       PP(I)=0D0
       PM(I)=0D0
       QP(I)=0D0
       QM(I)=0D0
        D(I)=0D0
      ENDDO
C
      DO IEDGE=1,NEDGE
C
       I=MAX(INOD(IEDGE),JNOD(IEDGE))
       J=MIN(INOD(IEDGE),JNOD(IEDGE))
C
       DAUX=FLUX(IEDGE)+FLUX0(IEDGE)
C
       FLUX(IEDGE)=DAUX
C 
       PP(I)=PP(I)+MAX(0D0,DAUX)
       PM(I)=PM(I)+MIN(0D0,DAUX)
C
       PP(J)=PP(J)+MAX(0D0,-DAUX)
       PM(J)=PM(J)+MIN(0D0,-DAUX)
C
       DAUX=DU(J)-DU(I)
C
       QP(I)=MAX(QP(I),DAUX)
       QM(I)=MIN(QM(I),DAUX)
C
       QP(J)=MAX(QP(J),-DAUX)
       QM(J)=MIN(QM(J),-DAUX)
C
      ENDDO
C
      DO I=1,NU
        QP(I)=DML(I)*QP(I)
        QM(I)=DML(I)*QM(I)
        IF (PP(I).GT.QP(I)+DEPS) THEN
          PP(I)=QP(I)/PP(I)
        ELSE
          PP(I)=1D0
        ENDIF
        IF (PM(I).LT.QM(I)-DEPS) THEN
          PM(I)=QM(I)/PM(I)
        ELSE
          PM(I)=1D0
        ENDIF
      ENDDO
C
      DO IEDGE=1,NEDGE
C
       I=MAX(INOD(IEDGE),JNOD(IEDGE))
       J=MIN(INOD(IEDGE),JNOD(IEDGE))
C
       DAUX=FLUX(IEDGE)
C
       IF (DAUX.GT.0) THEN
        DAUX=MIN(PP(I),PM(J))*DAUX
       ELSE
        DAUX=MIN(PM(I),PP(J))*DAUX
       ENDIF
C
       IF (KNPR(I).NE.0.OR.KNPR(J).NE.0) DAUX=0D0 
C
       D(I)=D(I)+DAUX
       D(J)=D(J)-DAUX
C
      ENDDO
C
      DO I=1,NU
       DU(I)=DU(I)+D(I)/DML(I)
      ENDDO
C
      END
C
C
************************************************************************
      SUBROUTINE XDFKD(KB1,KB2,KB3,KCOLB,KLDB,KFP,KU1,KU2,KU3,NU,NP,
     *                 TSTEPH,LELBD)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      IF (IPRECB.EQ.0) THEN
       CALL  LTX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KFP),1D0/TSTEPH,0D0)
       CALL  LTX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KFP),1D0/TSTEPH,1D0)
       CALL  LTX39 (VWORK(KB3),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU3),DWORK(KFP),1D0/TSTEPH,1D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.1) THEN
       CALL  LTX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KFP),1D0/TSTEPH,0D0)
       CALL  LTX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KFP),1D0/TSTEPH,1D0)
       CALL  LTX19 (DWORK(KB3),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU3),DWORK(KFP),1D0/TSTEPH,1D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.2) THEN
       CALL BTMUL1(KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *             DWORK(KFP),NEL,NVT,NAT,1D0/TSTEPH)
      ENDIF
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE XDFKG(KM1,DRHO,KB1,KB2,KB3,KCOLB,KLDB,KD1,KD2,KD3,
     *                 KU1,KU2,KU3,KP,NU,NP,TSTEPH,ILEV)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
C
      PARAMETER (NNARR=299,NNWORK=1,NNLEV=9)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
      REAL*8    DRHO(*)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /MGPER/  KKERI(NNLEV),KKERP(NNLEV),KKERV(NNLEV),
     *                NVERPER(NNLEV),KKROS(NNLEV),LKERI,LKERP,LKROS
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      IF (IPRECB.EQ.0) THEN
       CALL  LAX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD1),TSTEPH,0D0)
       CALL  LAX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD2),TSTEPH,0D0)
       CALL  LAX39 (VWORK(KB3),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD3),TSTEPH,0D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.1) THEN
       CALL  LAX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD1),TSTEPH,0D0)
       CALL  LAX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD2),TSTEPH,0D0)
       CALL  LAX19 (DWORK(KB3),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD3),TSTEPH,0D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.2) THEN
       CALL BMUL1 (KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KP),DWORK(KD1),DWORK(KD2),
     *             DWORK(KD3),KWORK(L(LKERI)),KWORK(L(LKERP)),
     *             KWORK(L(LKROS)),KWORK(L(LABD)),TSTEPH,0D0,ILEV)
      ENDIF
C
      CALL BDRY0 (DWORK(KD1),KWORK(L(LABD)),NABD,KWORK(L(LNPR)),NVT)
      CALL BDRY0 (DWORK(KD2),KWORK(L(LABD)),NABD,KWORK(L(LNPR)),NVT)
      CALL BDRY0 (DWORK(KD3),KWORK(L(LABD)),NABD,KWORK(L(LNPR)),NVT)
C
      CALL DIVPERM(DWORK(KM1),DRHO,DWORK(KD1),DWORK(KD2),DWORK(KD3),
     *             DWORK(KU1),DWORK(KU2),DWORK(KU3),KWORK(L(LKERI)),
     *             KWORK(L(LKROS)),KWORK(L(LABD)),NU,NABD,ILEV)
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE RESDFK(D1,D2,D3,F1,F2,F3,NU,RESU1,RESU2,RESU3)
************************************************************************
*    Purpose:  Computes the norms RESU
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z)
      DIMENSION D1(*),D2(*),D3(*),F1(*),F2(*),F3(*)
      SAVE 
C
C-----------------------------------------------------------------------
C     Compute the relative l2-norms  RESU1,RESU2,RESU3
C-----------------------------------------------------------------------
      CALL LL21 (F1,NU,RESF1)
      CALL LL21 (F2,NU,RESF2)
      CALL LL21 (F3,NU,RESF3)
      RESF=MAX(1D-12,RESF1,RESF2,RESF3)
C
      CALL LL21 (D1,NU,RESU1)
      CALL LL21 (D2,NU,RESU2)
      CALL LL21 (D3,NU,RESU3)
      RESU1=RESU1/RESF
      RESU2=RESU2/RESF
      RESU3=RESU3/RESF
C
      END
C
C
C
************************************************************************
      SUBROUTINE XDFKDV(KP,KFP,VMP,NP,A1)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VMP(*)
      DIMENSION VWORK(1),KWORK(1)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IP=1,NP
      DWORK(KP+IP-1)=DWORK(KP+IP-1)+A1*NY*DWORK(KFP+IP-1)/DBLE(VMP(IP))
10    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE VADM1(VX,VY,VZ,NX,AY,AZ)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*),VY(*),VZ(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      DO 10 IX=1,NX
10    VX(IX)=AY*VY(IX)+AZ*VZ(IX)
C
      END
C
C
C
************************************************************************
      SUBROUTINE VADM2(VX,VY,NX,AY)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      DO 10 IX=1,NX
10    VX(IX)=AY*VY(IX)
C
      END
C
C
C
************************************************************************
      SUBROUTINE DADM2(DX,DY,NX,AY)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      DO 10 IX=1,NX
10    DX(IX)=AY*DY(IX)
C
      END
C
C
C
************************************************************************
      SUBROUTINE V1LAX7(VA,KCOL,KLD,NEQ,DX1,DX2,DX3,DAX1,DAX2,DAX3,
     *                  A1,A2)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION DX1(*),DX2(*),DX3(*),DAX1(*),DAX2(*),DAX3(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE
C
      DO 1 IEQ=1,NEQ
      DAX1(IEQ)=A2*DAX1(IEQ)
      DAX2(IEQ)=A2*DAX2(IEQ)
      DAX3(IEQ)=A2*DAX3(IEQ)
C
      DO 2 ICOL=KLD(IEQ),KLD(IEQ+1)-1
      JCOL  =KCOL(ICOL)
      DAICOL=A1*DBLE(VA(ICOL))
      IF (DAICOL.EQ.0D0) GOTO 2
C
      DAX1(IEQ)=DAX1(IEQ)+DAICOL*DX1(JCOL)
      DAX2(IEQ)=DAX2(IEQ)+DAICOL*DX2(JCOL)
      DAX3(IEQ)=DAX3(IEQ)+DAICOL*DX3(JCOL)
C
2     CONTINUE
1     CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE D1LAX7(DA,KCOL,KLD,NEQ,DX1,DX2,DX3,DAX1,DAX2,DAX3,
     *                  A1,A2)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*)
      DIMENSION DX1(*),DX2(*),DX3(*),DAX1(*),DAX2(*),DAX3(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE
C
      DO 1 IEQ=1,NEQ
      DAX1(IEQ)=A2*DAX1(IEQ)
      DAX2(IEQ)=A2*DAX2(IEQ)
      DAX3(IEQ)=A2*DAX3(IEQ)
C
      DO 2 ICOL=KLD(IEQ),KLD(IEQ+1)-1
      JCOL  =KCOL(ICOL)
      DAICOL=A1*DA(ICOL)
      IF (DAICOL.EQ.0D0) GOTO 2
C
      DAX1(IEQ)=DAX1(IEQ)+DAICOL*DX1(JCOL)
      DAX2(IEQ)=DAX2(IEQ)+DAICOL*DX2(JCOL)
      DAX3(IEQ)=DAX3(IEQ)+DAICOL*DX3(JCOL)
C
2     CONTINUE
1     CONTINUE
C
      END

