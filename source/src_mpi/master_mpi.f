************************************************************************
      SUBROUTINE NSDEF0(MFILE,MSHOW,BSTOP,BNLEND,INL,DEF)  
************************************************************************
*   Purpose: - solver for the stationary Burgers-equation
*              equations via fixed point defect correction plus
*              multigrid for linear problems
*            - nonlinear version:
*                   - fixed point defect correction as outer iteration
*                   - mg as solver for the linear auxiliary
*                     problems
*                   - nonlinear parameter optimization for the 
*                     correction from the linear solver step
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      USE PP3D_MPI

      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120,CFILE*60
C
C *** Arrays for multigrid modul M010 
      DIMENSION  KOFFX(NNLEV),KOFFD(NNLEV),KOFFB(NNLEV),KNEQ(NNLEV),
     *           KIT(NNLEV),KIT0(NNLEV)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
      COMMON /EXTRA/ KLKX(NNLEV),KLKY(NNLEV),KLKZ(NNLEV),KLDK(NNLEV),
     *               LRHS,LRHS1,LRHS2,LRHS3,LNUT,KLNUT(NNLEV),
     *               LISEP,LIAUX,LINOD,LJNOD,LAEDGE,
     *               LF0U1,LF0U2,LF0U3,LF1U1,LF1U2,LF1U3
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC,TTILU,
     *                TTMGU,TTSU,TTEU,TTDU,TTPU,TTRU,
     *                TTMGP,TTSP,TTEP,TTDP,TTPP,TTRP
      COMMON /NSCOUN/ NNONL,NMGU,NMGP
      COMMON /NSEXL/  ITEXL,LTML,TIML11,TIML12,TIML31,TIML32
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGPER/  KKERI(NNLEV),KKERP(NNLEV),KKERV(NNLEV),
     *                NVERPER(NNLEV),KKROS(NNLEV),LKERI,LKERP,LKROS
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGIEL/  KLINT(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      COMMON /WALBC/  KKFBD(NNLEV),KAFBD(NNLEV),KNFBD(NNLEV)
      COMMON /LESPAR/ D_SMAG,I_SMAG
      SAVE
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Coefficient of stiffness matrix
      EXTERNAL UE
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31
C *** Multigrid components
      EXTERNAL  YAXU,YPROLU,YRESTU,YSMU,YEXU,YEXAU,YDBCU,YSTEPU

      Real :: dtt0 = 0.0
      Real :: dtt1 = 0.0
C
C=======================================================================
C     First generation of the nonlinear block A on level NLMAX
C=======================================================================

      write(*,*)'NSDEF0 deprecated'
      stop
C
      BSTOP =.FALSE.
      BNLEND=.FALSE.
C
C=======================================================================
C *** Loop of nonlinear iteration
C=======================================================================
C
      INLComplete = 0
C
      DO 222  INL=1,INLMAX
      NNONL=NNONL+1
C
      ISETLV=2
      ILEV=NLMAX
      CALL SETLEV (ISETLV)
C
      CALL COMMDISTRU(DWORK(KU1),NU)
      CALL COMMDISTRU(DWORK(KU2),NU)
      CALL COMMDISTRU(DWORK(KU3),NU)
C
      CALL ZTIME(dtt0)
      IF (IUPW.EQ.0) THEN
        CALL MATDEF_AFC(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *       DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *       DAUX,
     *       DWORK(L(LF1U1)),DWORK(L(LF1U2)),DWORK(L(LF1U3)),
     *       VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),
     *       KWORK(KCOLA),KWORK(KLDA),NU,KWORK(L(LVERT)),
     *       KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(NLEV))),
     *       DWORK(L(LCORVG)),DWORK(L(LNUT)),IELT,
     *       DWORK(KM1),DWORK(KMASS1),KWORK(L(LISEP)),
     *       KWORK(L(LIAUX)),KWORK(L(LINOD)),KWORK(L(LJNOD)),
     *       DWORK(L(LAEDGE)),0)
      ELSE
        CALL MATDEF(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *       DWORK(KF1),DWORK(KF2),DWORK(KF3),DAUX,
     *       VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),KWORK(KCOLA),
     *       KWORK(KLDA),NU,KWORK(L(LVERT)),KWORK(L(LAREA)),
     *       KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),
     *       DWORK(L(LNUT)),IELT,DWORK(KM1),DWORK(KMASS1),0)
      END IF
      CALL ZTIME(dtt1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL ZTIME(dtt0)
      CALL BDRYS(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *           DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *           VWORK(KA1),KWORK(KLDA),KWORK(L(LABD)),
     *           NABD,KWORK(L(LNPR)),DWORK(L(KNFBD(ILEV))),
     *           NVT,THSTEP,ISTOK)
C
      CALL ZTIME(dtt1)
      TTBDR=TTBDR+TTT1-TTT0  
C
C=======================================================================
C     matrix sorting ISORTU > 0
C=======================================================================
      CALL ZTIME(TTT0)
      IF ((ISORTU.GT.0).OR.(ISMU.EQ.4).OR.(ISLU.EQ.4)) THEN
C
       ILEV=NLMIN
C
       IF (ISORTU.GT.0) THEN
        CALL ZNEW(KNA(ILEV)  ,-2,LAH  ,'VAH   ')
        CALL ZNEW(KNA(ILEV)  ,-3,LCOLH,'KCOLH ')
        CALL ZNEW(KNU(ILEV)+1,-3,LLDH ,'KLDH  ')
        IF (IER.NE.0) GOTO 99998
C
        CALL ZCPY(KLA(ILEV)   ,'VA    ',LAH  ,'VAH   ')
        CALL ZCPY(KLCOLA(ILEV),'KCOLA ',LCOLH,'KCOLH ')
        CALL ZCPY(KLLDA(ILEV) ,'KLDA  ',LLDH, 'KLDH  ')
        IF (IER.NE.0) GOTO 99998
C
        CALL MTSRTV(VWORK(L(KLA(ILEV)))   ,VWORK(L(LAH)),
     *              KWORK(L(KLCOLA(ILEV))),KWORK(L(LCOLH)),
     *              KWORK(L(KLLDA(ILEV))) ,KWORK(L(LLDH)),
     *              KWORK(L(KLTRA1(ILEV))),KWORK(L(KLTRA2(ILEV))),
     *              KNU(ILEV))
C
        CALL ZDISP(0,LLDH ,'KLDH  ')
        CALL ZDISP(0,LCOLH,'KCOLH ')
        CALL ZDISP(0,LAH  ,'DAH   ')
        IF (IER.NE.0) GOTO 99998
       ENDIF
C
       IF ((ISMU.EQ.4).OR.(ISLU.EQ.3).OR.(ISLU.EQ.4)) THEN
        CALL ZNEW(KNA(ILEV),-2,KLAILU(ILEV),'VVAILU')
        IF (IER.NE.0) GOTO 99998
        CALL ZCPY(KLA(ILEV),'VA    ',KLAILU(ILEV),'VVAILU')
        IF (IER.NE.0) GOTO 99998
C
        TOLILU=1D-12
        ALPILU=0.0D0
        INDILU=1
        CALL IFD27(VWORK(L(KLAILU(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),KNU(ILEV),INDILU,ALPILU,TOLILU)
        IF (IER.NE.0) GOTO 99998
       ENDIF
C
      ENDIF
C
      CALL ZTIME(TTT1)
      TTILU=TTILU+TTT1-TTT0
C
C
C=======================================================================
C *** Initialization of the offset arrays KOFFX,KOFFB,KOFFD and KNEQ
C=======================================================================
C
      ILEV=NLMIN
      KOFFX(ILEV)=L(KLUP  (ILEV))-1
      KOFFB(ILEV)=L(KLF12P(ILEV))-1
      KOFFD(ILEV)=L(KLAUX (ILEV))-1
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
C
      ICYCLE=ICYCU
      IRELMG=1
      ITMG =ILMINU
      EPSUMG=1D99
      IF (ILMAXU.GT.ILMINU) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
      CALL  LCL1 (DWORK(1+KOFFX(NLMAX)),KNEQ(NLMAX))
      CALL  M0110(DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXU,ITMG,DMPUMG,EPSUMG,DEFUMG,
     *            YAXU,YPROLU,YRESTU,YSMU,YSMU,YEXU,YEXAU,YDBCU,YSTEPU,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMG1,BMGU1)
      NMGU=NMGU+ITMG
C
      TTMGU=TTMGU+TTMG
      TTSU=TTSU+TTS
      TTEU=TTEU+TTE
      TTDU=TTDU+TTD
      TTPU=TTPU+TTP
      TTRU=TTRU+TTR
C
C
      ILEV=NLMIN
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+KNU(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+KNU(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+KNU(ILEV)
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
C
      ICYCLE=ICYCU
      IRELMG=1
      ITMG =ILMINU
      EPSUMG=1D99
      IF (ILMAXU.GT.ILMINU) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
      CALL  LCL1 (DWORK(1+KOFFX(NLMAX)),KNEQ(NLMAX))
      CALL  M0110(DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXU,ITMG,DMPUMG,EPSUMG,DEFUMG,
     *            YAXU,YPROLU,YRESTU,YSMU,YSMU,YEXU,YEXAU,YDBCU,YSTEPU,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMG2,BMGU2)
      NMGU=NMGU+ITMG
C
      TTMGU=TTMGU+TTMG
      TTSU=TTSU+TTS
      TTEU=TTEU+TTE
      TTDU=TTDU+TTD
      TTPU=TTPU+TTP
      TTRU=TTRU+TTR
C
      ILEV=NLMIN
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+KNU(ILEV)+KNU(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+KNU(ILEV)+KNU(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+KNU(ILEV)+KNU(ILEV)
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
C
      ICYCLE=ICYCU
      IRELMG=1
      ITMG =ILMINU
      EPSUMG=1D99
      IF (ILMAXU.GT.ILMINU) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
      CALL  LCL1 (DWORK(1+KOFFX(NLMAX)),KNEQ(NLMAX))
      CALL  M0110(DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXU,ITMG,DMPUMG,EPSUMG,DEFUMG,
     *            YAXU,YPROLU,YRESTU,YSMU,YSMU,YEXU,YEXAU,YDBCU,YSTEPU,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMG3,BMGU3)
      NMGU=NMGU+ITMG
C
      TTMGU=TTMGU+TTMG
      TTSU=TTSU+TTS
      TTEU=TTEU+TTE
      TTDU=TTDU+TTD
      TTPU=TTPU+TTP
      TTRU=TTRU+TTR
C
      RHOLMG=MAX(RHOMG1,RHOMG2,RHOMG3)
C
C=======================================================================
C     matrix resorting ISORTU > 0
C=======================================================================
C
      CALL ZTIME(TTT0)
      IF ((ISORTU.GT.0).OR.(ISMU.EQ.4).OR.(ISLU.EQ.3).OR.
     *    (ISLU.EQ.4)) THEN
C
       ILEV=NLMIN
C
       IF (ISORTU.GT.0) THEN
        CALL ZNEW(KNA(ILEV)  ,-3,LCOLH,'KCOLH ')
        CALL ZNEW(KNU(ILEV)+1,-3,LLDH ,'KLDH  ')
        IF (IER.NE.0) GOTO 99998
C
        CALL ZCPY(KLCOLA(ILEV),'KCOLA ',LCOLH,'KCOLH ')
        CALL ZCPY(KLLDA(ILEV) ,'KLDA  ',LLDH, 'KLDH  ')
        IF (IER.NE.0) GOTO 99998
C
        CALL MTSRTR(KWORK(L(KLCOLA(ILEV))),KWORK(L(LCOLH)),
     *              KWORK(L(KLLDA(ILEV))) ,KWORK(L(LLDH)),
     *              KWORK(L(KLTRA1(ILEV))),KWORK(L(KLTRA2(ILEV))),
     *              KNU(ILEV))
C
         CALL ZDISP(0,LLDH ,'KLDH  ')
         CALL ZDISP(0,LCOLH,'KCOLH ')
         IF (IER.NE.0) GOTO 99998
       ENDIF
C
       IF ((ISMU.EQ.4).OR.(ISLU.EQ.3).OR.(ISLU.EQ.4)) THEN
        CALL ZDISP(0,KLAILU(ILEV) ,'VVAILU  ')
        IF (IER.NE.0) GOTO 99998
       ENDIF
C
      ENDIF
C
      CALL ZTIME(TTT1)
      TTILU=TTILU+TTT1-TTT0
C
C=======================================================================
C    End of multigrid iteration
C=======================================================================
C
      CALL COMM_NLComplete(INLComplete)
      if (INLComplete.eq.1) GOTO 221
C
C
222   CONTINUE
C
C=======================================================================
C *** End of the nonlinear loop
C=======================================================================
C
221   CONTINUE
C
C
C
99998 END
C
C
C
      SUBROUTINE  M0110(DX,DB,DD,KOFFX,KOFFB,KOFFD,KNEQ,NIT,ITE,
     *                  EPS1,EPS2,DEF,DAX,DPROL,DREST,DPRSM,DPOSM,
     *                  DEX,DEXA,DBC,DSTEP,KIT0,KIT,IREL,IDEFMG,RHOLMG,
     *                  BMGEND)
C
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION DX(*),DB(*),DD(*),KOFFX(*),KOFFB(*),KOFFD(*)
      DIMENSION KNEQ(*),KIT0(*),KIT(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      SAVE
C
C
      SUB='M0110 '
      IER=0
C
      BMGEND=.FALSE.
      MT0=MT
      MT=0
C
      ILEV=NLMAX
C
      CALL COMM_MGComplete(IMGComplete)
      IF (IMGComplete.EQ.1) GOTO 1
C
      CALL COMM_MGComplete(IMGComplete)
      IF (IMGComplete.EQ.1) GOTO 1
C
      BREPEAT=.true.
      DO WHILE (BREPEAT)
       CALL CommDex(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),
     *              DD(1+KOFFD(NLMIN)),KNEQ(NLMIN),RHOLMG)
       CALL COMM_MGComplete(IMGComplete)
       IF (IMGComplete.EQ.1) BREPEAT=.false.
      END DO
C
1     CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
      CALL DAX (DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *          -1D0,1D0)
      CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
C
      END


