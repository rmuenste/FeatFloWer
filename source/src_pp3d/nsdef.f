************************************************************************
      SUBROUTINE NSDEF(MFILE,MSHOW,BSTOP,BNLEND,INL,DEF)  
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
C
C=======================================================================
C     First generation of the nonlinear block A on level NLMAX
C=======================================================================
C
      BSTOP =.FALSE.
      BNLEND=.FALSE.
C
      ISETLV=2
      ILEV=NLMAX
      CALL SETLEV (ISETLV)
C
      IF (myid.EQ.MASTER) THEN                               ! PARALLEL
        CALL COMMDISTRU(DWORK(KU1),NU)                       ! PARALLEL
        CALL COMMDISTRU(DWORK(KU2),NU)                       ! PARALLEL
        CALL COMMDISTRU(DWORK(KU3),NU)                       ! PARALLEL
      END IF                                                 ! PARALLEL
C
      CALL ZTIME(TTT0)
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
     *       DWORK(L(LAEDGE)),-1)
      ELSE
        CALL MATDEF(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *       DWORK(KF1),DWORK(KF2),DWORK(KF3),DAUX,
     *       VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),KWORK(KCOLA),
     *       KWORK(KLDA),NU,KWORK(L(LVERT)),KWORK(L(LAREA)),
     *       KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),
     *       DWORK(L(LNUT)),IELT,DWORK(KM1),DWORK(KMASS1),-1)
      END IF
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL CommSum(DWORK(KF1),NLEV)          ! PARALLEL
      CALL CommSum(DWORK(KF2),NLEV)          ! PARALLEL
      CALL CommSum(DWORK(KF3),NLEV)          ! PARALLEL
C
      CALL ZTIME(TTT0)
      CALL BDRYS(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *           DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *           VWORK(KA1),KWORK(KLDA),KWORK(L(LABD)),
     *           NABD,KWORK(L(LNPR)),DWORK(L(KNFBD(ILEV))),
     *           NVT,THSTEP,ISTOK)
C
C     Imposing the periodic BC's
      CALL BDRPER(DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *            VWORK(KA1),KWORK(KLDA),KWORK(KCOLA),
     *            KWORK(L(LABD)),NABD,
     *            KWORK(L(KKERI(ILEV))),ISTOK)
C
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0  
C
      CALL ZTIME(TTT0)
      CALL LCP1(DWORK(KU1),DWORK(L(LU1OLD)),NU)
      CALL LCP1(DWORK(KU2),DWORK(L(LU2OLD)),NU)
      CALL LCP1(DWORK(KU3),DWORK(L(LU3OLD)),NU)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C     Calculation of initial defects
C=======================================================================
C
      CALL ZTIME(TTT0)
      CALL RESDFK(DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *            DWORK(L(LRHS1)),DWORK(L(LRHS2)),DWORK(L(LRHS3)),
     *            NU,RESU1,RESU2,RESU3)
C
      DEF=MAX(RESU1,RESU2,RESU3)
C
      RES   =DEF
      RES0  =DEF
      RESOLD=DEF
      EPSRES=DMPUD*DEF
C
C
C=======================================================================
C     Solution
C=======================================================================
C
      INL=0
C
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      IF (MSHOW.GE.2) WRITE(MTERM,1001)
      IF (MSHOW.GE.1) WRITE(MFILE,1001)
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      IF (MSHOW.GE.2) WRITE(MTERM,1002)  INL,RESU1,RESU2,RESU3
      IF (MSHOW.GE.1) WRITE(MFILE,1002)  INL,RESU1,RESU2,RESU3
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      OMEGA=OMGINI
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Loop of nonlinear iteration
C=======================================================================
C
      INLComplete = 0                                          ! PARALLEL
C
      DO 222  INL=1,INLMAX
      NNONL=NNONL+1
C
C=======================================================================
C *** Generate the A blocks for all coarse levels
C=======================================================================
C
      IF (NLMAX.GT.NLMIN)  THEN
       DO 22  ILEV=NLMAX-1,NLMIN,-1
       CALL ZTIME(TTT0)
       ISETLV=2
       CALL  SETLEV (ISETLV)
C
       I1=ILEV+1
       KU1F=L(KLUP(I1))
       KU2F=KU1F+KNU(I1)
       KU3F=KU2F+KNU(I1)
       KVERTF=L(KLVERT(I1))
       KAREAF=L(KLAREA (I1))
       KADJF =L(KLADJ (I1))
       CALL RESTRU(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *             DWORK(KU1F),DWORK(KU2F),DWORK(KU3F),
     *             KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LADJ)),
     *             NU,NP,NVT,KWORK(KVERTF),KWORK(KAREAF),
     *             KWORK(KADJF),KNU(I1),KNP(I1),KNVT(I1),
     *             VWORK(L(KLVOL(ILEV))),1)
C
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
       CALL CommSumHalf(DWORK(KU1),ILEV)                       ! PARALLEL
       CALL CommSumHalf(DWORK(KU2),ILEV)                       ! PARALLEL
       CALL CommSumHalf(DWORK(KU3),ILEV)                       ! PARALLEL
       IF (ILEV.EQ.NLMIN) THEN                                 ! PARALLEL
        CALL COMMDISTRU(DWORK(KU1),NU)                         ! PARALLEL
        CALL COMMDISTRU(DWORK(KU2),NU)                         ! PARALLEL
        CALL COMMDISTRU(DWORK(KU3),NU)                         ! PARALLEL
       END IF                                                  ! PARALLEL
C
       CALL ZTIME(TTT0)
       IF (IUPW.EQ.0) THEN
         CALL MATDEF_AFC(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *        DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *        DAUX,
     *        DWORK(L(LF1U1)),DWORK(L(LF1U2)),DWORK(L(LF1U3)),
     *        VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),
     *        KWORK(KCOLA),KWORK(KLDA),NU,KWORK(L(LVERT)),
     *        KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),
     *        DWORK(L(LCORVG)),DWORK(L(LNUT)),IELT,
     *        DWORK(KM1),DWORK(KMASS1),KWORK(L(LISEP)),
     *        KWORK(L(LIAUX)),KWORK(L(LINOD)),KWORK(L(LJNOD)),
     *        DWORK(L(LAEDGE)),0)
       ELSE
         CALL MATDEF(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *        DWORK(KF1),DWORK(KF2),DWORK(KF3),DAUX,
     *        VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),KWORK(KCOLA),
     *        KWORK(KLDA),NU,KWORK(L(LVERT)),KWORK(L(LAREA)),
     *        KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),
     *        DWORK(L(LNUT)),IELT,DWORK(KM1),DWORK(KMASS1),0)
       END IF
       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       CALL BDRYA(VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            KWORK(L(KLABD(ILEV))),KNABD(ILEV),KWORK(L(LNPR)),NVT)
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
  22   CONTINUE
      ENDIF
C
      DO ILEV=NLMIN,NLMAX                                       ! PARALLEL
       CALL CommMat(VWORK(L(KLA(ILEV))),KWORK(L(KLLDA(ILEV))),  ! PARALLEL
     *      KNU(ILEV),ILEV)                                     ! PARALLEL
      END DO                                                    ! PARALLEL
!      write(*,*) "Signal ",myid,inl
C
C
C=======================================================================
C     matrix sorting ISORTU > 0
C=======================================================================
      CALL ZTIME(TTT0)
      IF ((ISORTU.GT.0).OR.(ISMU.EQ.4).OR.(ISLU.EQ.4)) THEN
C
       DO 11  ILEV=NLMIN,NLMAX
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
        TOLILU=0d0
        ALPILU=0.0D0
        INDILU=1
        CALL IFD27(VWORK(L(KLAILU(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),KNU(ILEV),INDILU,ALPILU,TOLILU)
!        if (ilev.eq.nlmax) write(*,*) myid,inl,ier
!        IF (IER.NE.0) GOTO 99998
       ENDIF
C
11     CONTINUE
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
      DO 12  ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1
      KOFFB(ILEV)=L(KLF12P(ILEV))-1
      KOFFD(ILEV)=L(KLAUX (ILEV))-1
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
12    CONTINUE
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
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
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
      DO 13  ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+KNU(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+KNU(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+KNU(ILEV)
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
13    CONTINUE
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
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
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
      DO 14  ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+KNU(ILEV)+KNU(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+KNU(ILEV)+KNU(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+KNU(ILEV)+KNU(ILEV)
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
14    CONTINUE
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
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
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
       DO 15  ILEV=NLMIN,NLMAX
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
15     CONTINUE
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
C *** Set finest level NLMAX
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
C=======================================================================
C *** Calculate the optimal correction
C=======================================================================
C
      CALL  ZTIME(TTT0)
      CALL XOPTCN(DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(L(LU1OLD)),
     *     DWORK(L(LU2OLD)),DWORK(L(LU3OLD)),NU,DELU1,DELU2,DELU3,OMEGA)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
      IF ((BMGU1.OR.BMGU2.OR.BMGU3).AND.(ABS(OMEGA).GE.1D-2)) 
     *      BSTOP=.TRUE.
C
C=======================================================================
C *** Calculation of defects
C=======================================================================
C
      IF (INLMAX.GT.1) THEN
       CALL  ZTIME(TTT0)
       CALL LCP1 (DWORK(L(LRHS1)),DWORK(KF1),NU)
       CALL LCP1 (DWORK(L(LRHS2)),DWORK(KF2),NU)
       CALL LCP1 (DWORK(L(LRHS3)),DWORK(KF3),NU)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       IF (IUPW.EQ.0) THEN
         CALL MATDEF_AFC(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *        DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *        DAUX,
     *        DWORK(L(LF1U1)),DWORK(L(LF1U2)),DWORK(L(LF1U3)),
     *        VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),
     *        KWORK(KCOLA),KWORK(KLDA),NU,KWORK(L(LVERT)),
     *        KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(NLEV))),
     *        DWORK(L(LCORVG)),DWORK(L(LNUT)),IELT,
     *        DWORK(KM1),DWORK(KMASS1),KWORK(L(LISEP)),
     *        KWORK(L(LIAUX)),KWORK(L(LINOD)),KWORK(L(LJNOD)),
     *        DWORK(L(LAEDGE)),-1)
       ELSE
         CALL MATDEF(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *        DWORK(KF1),DWORK(KF2),DWORK(KF3),DAUX,
     *        VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),KWORK(KCOLA),
     *        KWORK(KLDA),NU,KWORK(L(LVERT)),KWORK(L(LAREA)),
     *        KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),
     *        DWORK(L(LNUT)),IELT,DWORK(KM1),DWORK(KMASS1),-1)
       END IF
       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       CALL BDRYS(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *            DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *            VWORK(KA1),KWORK(KLDA),KWORK(L(LABD)),
     *            NABD,KWORK(L(LNPR)),DWORK(L(KNFBD(ILEV))),
     *            NVT,THSTEP,ISTOK)
C
C     Imposing the periodic BC's
      CALL BDRPER(DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *            VWORK(KA1),KWORK(KLDA),KWORK(KCOLA),
     *            KWORK(L(LABD)),NABD,
     *            KWORK(L(KKERI(ILEV))),ISTOK)
C
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       CALL LCP1(DWORK(KU1),DWORK(L(LU1OLD)),NU)
       CALL LCP1(DWORK(KU2),DWORK(L(LU2OLD)),NU)
       CALL LCP1(DWORK(KU3),DWORK(L(LU3OLD)),NU)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Calculation of defects and norms
C=======================================================================
C
       CALL CommSum(DWORK(KF1),NLEV)          ! PARALLEL
       CALL CommSum(DWORK(KF2),NLEV)          ! PARALLEL
       CALL CommSum(DWORK(KF3),NLEV)          ! PARALLEL
C
       CALL  ZTIME(TTT0)
       CALL  RESDFK(DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *              DWORK(L(LRHS1)),DWORK(L(LRHS2)),DWORK(L(LRHS3)),
     *              NU,RESU1,RESU2,RESU3)
C
!       WRITE(*,*) "receiving ... ", myid,RESU1,RESU2,RESU3
       RES=MAX(RESU1,RESU2,RESU3)
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
!      WRITE(*,*) "critical:",myid,BSTOP,RES,RESOLD
       IF ( RESOLD.GT.1d-12.AND.RES.GT.1D2*RESOLD) BSTOP=.TRUE.
C
       RESOLD=RES
       RHO   =(RES/RES0)**(1D0/DBLE(INL))
       DELU  =MAX(DELU1,DELU2,DELU3)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ELSE
       RES   =0D0
       RESU1 =0D0
       RESU2 =0D0
       RESU3 =0D0
       RESOLD=0D0
       RHO   =0D0
       DELU  =MAX(DELU1,DELU2,DELU3)
      ENDIF
C
C=======================================================================
C *** Control of terminating the nonlinear iteration
C=======================================================================
C
      IF ((DELU.LE.EPSUR).AND.(RES.LE.EPSUD).AND.!(RES.LE.EPSRES).AND.
     *    (INL.GE.INLMIN))  BNLEND=.TRUE.
      IF ((DELU.LE.EPSUR).AND.(INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.1))
     *                      BNLEND=.TRUE.
C
       IF (MSHOW.GE.2) 
     * WRITE(MTERM,1003)  INL,DELU1,DELU2,DELU3,RESU1,RESU2,RESU3,
     *                    RHO,OMEGA,RHOMG1,RHOMG2,RHOMG3
!        WRITE(MTERM,1003)  myid,DELU1,DELU2,DELU3,RESU1,RESU2,RESU3,
!      *                    RHO,OMEGA,RHOMG1,RHOMG2,RHOMG3
      IF (MSHOW.GE.1) 
     * WRITE(MFILE,1003)  INL,DELU1,DELU2,DELU3,RESU1,RESU2,RESU3,
     *                    RHO,OMEGA,RHOMG1,RHOMG2,RHOMG3
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
!      THIS HAS TO BE A QUESTION OF COMMUNICATION !!!!
!       write(*,*) myid,BSTOP
       IF (BSTOP) RETURN
C
C=======================================================================
C *** Autosave
C=======================================================================
C
      IF (IAUSAV.NE.0) THEN
      IF (MOD(INL,IAUSAV).EQ.0) THEN
       CALL ZTIME(TTT0)
       CFILE='#data/#AUTOSAV   '
       CALL  OF0 (39,CFILE,0)
       CALL  OWA1 (DWORK(KU1),'DU12P ',NUP, 39,0)
       REWIND(39)
       CLOSE(39)
       CALL ZTIME(TTT1)
       TTPOST=TTPOST+TTT1-TTT0
      ENDIF
      ENDIF
C
C=======================================================================
C *** Return if BNLEND=true
C=======================================================================
      IF ((BNLEND).OR.((INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.1))
     *            .OR.(ABS(OMEGA).LT.1D-1)) INLComplete = 1
C
      CALL COMM_NLComplete(INLComplete)                            ! PARALLEL
      if (INLComplete.eq.1) GOTO 221
C
C
222   CONTINUE
C
C=======================================================================
C *** End of the nonlinear loop
C=======================================================================
C
221   IF (MSHOW.GE.6) THEN
      WRITE(MTERM,*)
      WRITE(MTERM,*) ' U-MULTIGRID COMPONENTS [in percent]:',TTMGU
      WRITE(MTERM,*) ' smoothing     :', 1.D2*TTSU/TTMGU
      WRITE(MTERM,*) ' solver        :', 1.D2*TTEU/TTMGU
      WRITE(MTERM,*) ' defect calc.  :', 1.D2*TTDU/TTMGU
      WRITE(MTERM,*) ' prolongation  :', 1.D2*TTPU/TTMGU
      WRITE(MTERM,*) ' restriction   :', 1.D2*TTRU/TTMGU
      WRITE(MTERM,1)
      ENDIF
C
      IF (MSHOW.GE.5) THEN
      WRITE(MFILE,*)
      WRITE(MFILE,*) ' U-MULTIGRID COMPONENTS [in percent]:',TTMGU
      WRITE(MFILE,*) ' smoothing     :', 1.D2*TTSU/TTMGU
      WRITE(MFILE,*) ' solver        :', 1.D2*TTEU/TTMGU
      WRITE(MFILE,*) ' defect calc.  :', 1.D2*TTDU/TTMGU
      WRITE(MFILE,*) ' prolongation  :', 1.D2*TTPU/TTMGU
      WRITE(MFILE,*) ' restriction   :', 1.D2*TTRU/TTMGU
      WRITE(MFILE,1)
      ENDIF
C
      IF (INL.GT.INLMAX)  INL=INLMAX
C
   1  FORMAT(80('-'))
1001  FORMAT(' IT REL-U1',3X,'REL-U2',3X,'REL-U3',3X,
     *       'DEF-U1',3X,'DEF-U2',3X,'DEF-U3',3X,
     *       'RHONL ',3X,'OMEGNL',3X,'RHOMG1',3X,'RHOMG2',3X,'RHOMG3')
1002  FORMAT(I3,27X,3(D9.2))
1003  FORMAT(I3,12(D9.2))
C
C
C
99998 END

      SUBROUTINE SHOWMAT(VA,KLD,NEQ)
      REAL*4 VA(*)
      INTEGER NEQ,KLD(*)

      DO I=1,NEQ
      J1=KLD(I)
      J2=KLD(I+1)-1
      WRITE(*,'(A2,15(1XG12.5))') '*|',(VA(J),J=J1,J2)
      END DO

      END
