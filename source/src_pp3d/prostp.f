************************************************************************
      SUBROUTINE PROSTP(MFILE,MSHOW,IPROJ,BSTOP,BNLEND,INONLN,DEF)
************************************************************************
*
*   Purpose: - performs 1 macro step of projection method
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      USE PP3D_MPI

      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNAB=21, NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NBLOCF=3)
      DIMENSION BCONF(NBLOCF),ARRDF(NBLOCF),BSNGLF(NBLOCF)
      DIMENSION KFN(NBLOCF),KF(NNAB,NBLOCF),LF(NBLOCF)
      CHARACTER ARRDF*6
      DATA BCONF /.FALSE.,.FALSE.,.FALSE./,ARRDF/'DF1   ','DF2   ',
     *           'DF3   '/
      DATA BSNGLF/.FALSE.,.FALSE.,.FALSE./,KFN/1,1,1/
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
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSSAV/  INSAV,INSAVN
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IGID,IGMV,IFINIT
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC,TTILU,
     *                TTMGU,TTSU,TTEU,TTDU,TTPU,TTRU,
     *                TTMGP,TTSP,TTEP,TTDP,TTPP,TTRP,TTLES
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
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
      COMMON /WALBC/ KKFBD(NNLEV),KAFBD(NNLEV),KNFBD(NNLEV)
      COMMON /DELTA/ KLDELT,iLDELT,LDSTRS
      COMMON /LESPAR/D_SMAG,I_SMAG
      SAVE
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Parametrization of the domain
      EXTERNAL PARX,PARY,PARZ
C *** Coefficient of stiffness matrix, right hand side, exact solution
      EXTERNAL COEFFN,RHS,UE,PE,UEX,UEY,UEZ
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30
C *** Multigrid components
      EXTERNAL  YAXP,YPROLP,YRESTP,YSMP,YEXP,YEXAP,YDBCP,YSTEPP
C
C=======================================================================
C *** Begin of nonsteady loop
C=======================================================================
C
      VARIABLE=1                              ! PARALLEL
      ISETLV=2
      ILEV=NLEV
      CALL SETLEV(ISETLV)
C
      THSTEP=TSTEP*(1D0-THETA)
      TMSTEP=TSTEP*THETA
C
      BSTOP=.FALSE.
      RETURN
C--------- Assembly of the constant right-hand side ---------
C     Contribution of the mass matrix / transport operator
      CALL ZTIME(TTT0)
      CALL LCL1(DWORK(KF1),3*NU)
C
! C     Add the gravity force term
!       CALL GravForce(DWORK(KM1),DWORK(KF1),DWORK(KF2),DWORK(KF3))
C
C     Add the surface tension force term
!       CALL XSurfTens()
C
      IF (IUPW.EQ.0) THEN
        CALL MATDEF_AFC(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *       DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *       DAUX,
     *       DWORK(L(LF0U1)),DWORK(L(LF0U2)),DWORK(L(LF0U3)),
     *       VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),
     *       KWORK(KCOLA),KWORK(KLDA),NU,KWORK(L(LVERT)),
     *       KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(NLEV))),
     *       DWORK(L(LCORVG)),DWORK(L(LNUT)),IELT,
     *       DWORK(KM1),DWORK(KMASS1),KWORK(L(LISEP)),
     *       KWORK(L(LIAUX)),KWORK(L(LINOD)),KWORK(L(LJNOD)),
     *       DWORK(L(LAEDGE)),1)
      ELSE
        CALL MATDEF(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *       DWORK(KF1),DWORK(KF2),DWORK(KF3),DAUX,
     *       VWORK(KA1),DWORK(L(KLDK(ILEV))),DWORK(KST1),KWORK(KCOLA),
     *       KWORK(KLDA),NU,KWORK(L(LVERT)),KWORK(L(LAREA)),
     *       KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),
     *       DWORK(L(LNUT)),IELT,DWORK(KM1),DWORK(KMASS1),1)
      END IF
C
C
C     Contribution of the pressure gradient
      IF (IPROJ.EQ.1.AND.myid.NE.MASTER) THEN
      CALL BMUL1(KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LADJ)),
     *            DWORK(L(LCORVG)),DWORK(KP),DWORK(KF1),DWORK(KF2),
     *            DWORK(KF3),KWORK(L(LKERI)),KWORK(L(LKERP)),
     *            KWORK(L(LKROS)),KWORK(L(LABD)),-TSTEP,1D0,NLEV)
      ENDIF
C
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C     Contribution of the surface stress
      CALL ZTIME(TTT0)
      CALL BDRYS(DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(KF1),
     *           DWORK(KF2),DWORK(KF3),VWORK(KA1),KWORK(KLDA),
     *           KWORK(L(LABD)),NABD,KWORK(L(LNPR)),
     *           DWORK(L(KNFBD(ILEV))),NVT,THSTEP,0)
C
C     Imposing the periodic BC's
      CALL BDRPER(DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *            VWORK(KA1),KWORK(KLDA),KWORK(KCOLA),
     *            KWORK(L(LABD)),NABD,
     *            KWORK(L(KKERI(ILEV))),1)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
C     Save the constant right-hand side
      CALL LCP1(DWORK(KF1),DWORK(L(LRHS)),3*NU)
C
      CALL ZTIME(TTT0)
C
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL BDRSET(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *            KWORK(L(KLABD(ILEV))),KNABD(ILEV),
     *            DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LAREA)),
     *            KWORK(L(KELBD(ILEV))),UE,
     *            KWORK(L(LNPR)),DWORK(L(KNFBD(ILEV))),NVT)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL LCP1(DWORK(KF1),DWORK(L(LRHS1)),NU)
      CALL LCP1(DWORK(KF2),DWORK(L(LRHS2)),NU)
      CALL LCP1(DWORK(KF3),DWORK(L(LRHS3)),NU)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      THSTEP=TMSTEP
C
C***********************************************************************
C *** fixed point defect correction for stationary Burgers-equ.
C***********************************************************************
C
      IF (myid.eq.MASTER) THEN
       CALL NSDEF0(MFILE,MSHOW,BSTOP,BNLEND,INONLN,DEF)
      ELSE
       CALL NSDEF(MFILE,MSHOW,BSTOP,BNLEND,INONLN,DEF)
      END IF
!      THIS HAS TO BE A QUESTION OF COMMUNICATION !!!!
!      IF (DEF.LT.1D-12) RETURN                   ! PARALLEL
C
      IF (IMASS.EQ.1.AND.IUPW.EQ.0) THEN
        NEDGE=(NA-NU)/2
        CALL DEFFCT(DWORK(KU1),DWORK(KM1),KWORK(L(LNPR)),
     *       NU,NEDGE,DWORK(L(LF1U1)),DWORK(L(LF0U1)),
     *       KWORK(L(LINOD)),KWORK(L(LJNOD)))
        CALL DEFFCT(DWORK(KU2),DWORK(KM1),KWORK(L(LNPR)),
     *       NU,NEDGE,DWORK(L(LF1U2)),DWORK(L(LF0U2)),
     *       KWORK(L(LINOD)),KWORK(L(LJNOD)))
        CALL DEFFCT(DWORK(KU3),DWORK(KM1),KWORK(L(LNPR)),
     *       NU,NEDGE,DWORK(L(LF1U3)),DWORK(L(LF0U3)),
     *       KWORK(L(LINOD)),KWORK(L(LJNOD)))
      ENDIF 
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
      IF (BSTOP) RETURN
C
      VARIABLE=2                       ! PARALLEL
      ILEV=NLEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
C
C***********************************************************************
C *** Calculation of QR = 1/K B^T U~
C***********************************************************************
C
!!!! TEMPORARY !!!!!!!!!
!      GOTO 99999
!!!! TEMPORARY !!!!!!!!!
      CALL ZTIME(TTT0)
      CALL XDFKD(KB1,KB2,KB3,KCOLB,KLDB,KFP,KU1,KU2,KU3,NU,NP,
     *           TSTEP,KELBD(ILEV))
      CALL LL21(DWORK(KFP),NP,DRES_BEF)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL LL21 (DWORK(KU1),NU,DNRMU1)
      CALL LL21 (DWORK(KU2),NU,DNRMU2)
      CALL LL21 (DWORK(KU3),NU,DNRMU3)
      DNORMU=SQRT(DNRMU1*DNRMU1+DNRMU2*DNRMU2+DNRMU3*DNRMU3)
C
C***********************************************************************
C *** Solution of B^T M^-1 B P = FP
C***********************************************************************
C
      CALL LCP1(DWORK(KP),DWORK(L(LPOLD)),NP)
      CALL LCL1(DWORK(KP),NP)
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      DO 30 ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+3*KNU(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+3*KNU(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+3*KNU(ILEV)
      KNEQ (ILEV)=KNP(ILEV)
      KPRSM(ILEV)=NSMP*NSMPFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMP*NSMPFA**(NLMAX-ILEV)
30    CONTINUE
C
      ICYCLE=ICYCP
      IRELMG=1
      ITMGP =ILMINP
      IF (ILMAXP.GT.ILMINP) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
C
      IF (myid.eq.MASTER) THEN
      CALL  M0110(DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXP,ITMGP,DMPPMG,DNORMU*EPSP/TSTEP,DEFPMG,
     *            YAXP,YPROLP,YRESTP,YSMP,YSMP,YEXP,YEXAP,YDBCP,YSTEPP,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMGP,BSTOP)
      ELSE
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXP,ITMGP,DMPPMG,DNORMU*EPSP/TSTEP,DEFPMG,
     *            YAXP,YPROLP,YRESTP,YSMP,YSMP,YEXP,YEXAP,YDBCP,YSTEPP,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMGP,BSTOP)
      END IF
      NMGP=NMGP+ITMGP
C
      TTMGP=TTMGP+TTMG
      TTSP =TTSP+TTS
      TTEP =TTEP+TTE
      TTDP =TTDP+TTD
      TTPP =TTPP+TTP
      TTRP =TTRP+TTR
C
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
      IF (MSHOW.GE.2) WRITE(MTERM,1001)
      IF (MSHOW.GE.1) WRITE(MFILE,1001)
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
      IF (MSHOW.GE.2) WRITE(MTERM,1003) ITMGP,TSTEP*DEFPMG/DNORMU,
     *                                  RHOMGP,TIMENS
      IF (MSHOW.GE.1) WRITE(MFILE,1003) ITMGP,TSTEP*DEFPMG/DNORMU,
     *                                  RHOMGP,TIMENS
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      IF (MSHOW.GE.6) THEN
      WRITE(MTERM,*)
      WRITE(MTERM,*) ' P-MULTIGRID COMPONENTS [in percent]:',TTMGP
      WRITE(MTERM,*) ' smoothing     :', 1.D2*TTSP/TTMGP
      WRITE(MTERM,*) ' solver        :', 1.D2*TTEP/TTMGP
      WRITE(MTERM,*) ' defect calc.  :', 1.D2*TTDP/TTMGP
      WRITE(MTERM,*) ' prolongation  :', 1.D2*TTPP/TTMGP
      WRITE(MTERM,*) ' restriction   :', 1.D2*TTRP/TTMGP
      WRITE(MTERM,1)
      ENDIF
C
      IF (MSHOW.GE.5) THEN
      WRITE(MFILE,*)
      WRITE(MFILE,*) ' P-MULTIGRID COMPONENTS [in percent]:',TTMGP
      WRITE(MFILE,*) ' smoothing     :', 1.D2*TTSP/TTMGP
      WRITE(MFILE,*) ' solver        :', 1.D2*TTEP/TTMGP
      WRITE(MFILE,*) ' defect calc.  :', 1.D2*TTDP/TTMGP
      WRITE(MFILE,*) ' prolongation  :', 1.D2*TTPP/TTMGP
      WRITE(MFILE,*) ' restriction   :', 1.D2*TTRP/TTMGP
      WRITE(MFILE,1)
      ENDIF
C
      IF (IER.GT.0) IER=0
      IF (IER.NE.0) GOTO 99999
C
C***********************************************************************
C
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
      IF (BSTOP) RETURN
C
C***********************************************************************
C *** Update of U = U~ - k B P 
C***********************************************************************
C
      CALL ZTIME(TTT0)
      IF (myid.NE.MASTER) THEN
      CALL XDFKG(KM1,DAUX,KB1,KB2,KB3,KCOLB,KLDB,
     *     L(LD1),L(LD2),L(LD3),KU1,KU2,KU3,KP,NU,NP,TSTEP,NLEV)
      END IF
C
!       CALL XDFKD(KB1,KB2,KB3,KCOLB,KLDB,KFP,KU1,KU2,KU3,NU,NP,
!      *           TSTEP,KELBD(ILEV))
!       CALL LL21(DWORK(KFP),NP,DRES_AFT)
!       WRITE(*,'(A43,2(2XG8.2))') 
!      *"Incompressibility residuum (before/after): ",DRES_BEF,DRES_AFT
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C ************************************************************************
C *** Update of P(n+1)=P(n) + alphap P for second order
C***********************************************************************
C
C
      CALL ZTIME(TTT0)
C
      IF (IPROJ.EQ.1) THEN
       ALPHAT=1D0
       ALPHAP=1D0/THETA
      ELSE
       ALPHAT=0D0
       ALPHAP=1D0
       IF (IPROJ.EQ.-1) IPROJ=1
       IF (IPROJ.LT.-1) IPROJ=IPROJ+1
      ENDIF
C
      CALL LL21(DWORK(KP),NP,RELP2)
      CALL LLI1(DWORK(KP),NP,RELPM,INDMAX)
      CALL LLC1(DWORK(L(LPOLD)),DWORK(KP),NP,ALPHAT,ALPHAP)
      CALL LL21(DWORK(KP),NP,DSXN)
      RELP2=RELP2/DSXN
      CALL LLI1(DWORK(KP),NP,DSXN,INDMAX)
      RELPM=RELPM/DSXN
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,10002)
     * ITNS,TIMENS,TSTEP,RELP2,RELPM
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,*)
C
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,10002)
     *  ITNS,TIMENS,TSTEP,RELP2,RELPM
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
C
C
C
      GOTO 99999
C
C
C

   1  FORMAT(80('-'))
   3  FORMAT(80('+'))
1000  FORMAT (6E12.5)
1001  FORMAT('  IT DIV-L2',3X,'RHOMGP',3X,' TIME ')
1003  FORMAT(I4,3(D9.2))
10002 FORMAT ('#(',I6,')',1X,'TIME=',D10.3,1X,'TSTEP=',
     *        D10.3,1X,'REL2(P)=',D10.3,1X,'RELM(P)=',D10.3)
C
C
C
99999 END
