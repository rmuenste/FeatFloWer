!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               This module contains FEATFLOW declarations                !   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE def_FEAT

  !=======================================================================
  !     Declarations
  !=======================================================================
  IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
  PARAMETER (NNDIM=3,NNARR=299,NNLEV=9,NNAB=21)
  PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6)
  CHARACTER SUB*6,FMT*15,CPARAM*120
  CHARACTER CDATA*60

  !-----------------------------------------------------------------------
  !     C O M M O N S 
  !-----------------------------------------------------------------------
  
  ! *** Standard COMMON blocks
  COMMON /ERRCTL/ IER,ICHECK
  COMMON /CHAR/   SUB,FMT(3),CPARAM
  COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
  COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,    &
                  NVAR,NEAR,NBCT,NVBD,NEBD,NABD
  COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,   &
                  LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,  &
                  LNPR,LBCT,LVBD,LEBD,LABD,LVOL
  COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,     &
                  DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE), &
                  IEL,NDIM
  COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
  COMMON /COAUX1/ KDFG(NNVE),KDFL(NNVE),IDFL

  ! *** COMMON blocks for time discretization
  COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
  COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
  COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL, &
                  EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
  COMMON /NSSAV/  INSAV,INSAVN
  COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,    &
                  IFUSAV,IFPSAV,IFXSAV,IAVS,IGMV,IFINIT       
  COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC,TTILU, &
                  TTMGU,TTSU,TTEU,TTDU,TTPU,TTRU,             &
                  TTMGP,TTSP,TTEP,TTDP,TTPP,TTRP
  COMMON /NSCOUN/ NNONL,NMGU,NMGP
  COMMON /NSEXL/  ITEXL,LTML,TIML11,TIML12,TIML31,TIML32


  ! *** COMMON blocks for multigrid data management
  COMMON /MGPER/  KKERI(NNLEV),KKERP(NNLEV),KKERV(NNLEV),    &
                  NVERPER(NNLEV),KKROS(NNLEV),LKERI,LKERP,LKROS
  COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,                     &
                  ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
  COMMON /MGTRD/  KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),     &
                  KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),    &
                  KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),    &
                  KNABD(NNLEV)
!  COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),       &
!                  KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),       &
!                  KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),     &
!                  KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),    &
!                  KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),    &
!                  KNABD(NNLEV)
  COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
  COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),    &
                  KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV), &
                  KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),    &
                  KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),    &
                  KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),    &
                  KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),    &
                  KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
  COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
  COMMON /MGIEL/  KLINT(NNLEV)
  COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),             &
                  KLM(NNLEV),KLCOLA(NNLEV),                         &
                  KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV), &
                  KLCOLB(NNLEV),KLLDB(NNLEV),                       &
                  KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),           &
                  KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,                &
                  LPOLD,LD1,LD2,LD3,LDP
  COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
  COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),KNUP(NNLEV)

  COMMON /LEVDIM/ NA,NB,NU,NP,NUP
  COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB, &
                  KLDB,KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,  &
                  KAUX3,KAUXP

  ! *** user COMMON block
  INTEGER  VIPARM 
  DIMENSION VIPARM(100)                     
  EQUIVALENCE (IAUSAV,VIPARM)   
  COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,             &
                 IMASS,IMASSL,IUPW,IPRECA,IPRECB,                &
                 ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,           &
                 INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,        &
                 ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP, &
                 IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

  DOUBLE PRECISION VRPARM,NY
  DIMENSION VRPARM(100)
  EQUIVALENCE (NY,VRPARM)                          
  COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,              &
                  EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU, &
                  AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,  &
                  AMINP,AMAXP

  CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
  COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1, &
                 ISTART,MSTART,CSTART,ISOL,MSOL,CSOL

  COMMON /EXTRA/ KLKX(NNLEV),KLKY(NNLEV),KLKZ(NNLEV),KLDK(NNLEV), &
                 LRHS,LRHS1,LRHS2,LRHS3,LNUT,KLNUT(NNLEV), &
                 LISEP,LIAUX,LINOD,LJNOD,LAEDGE, &
                 LF0U1,LF0U2,LF0U3,LF1U1,LF1U2,LF1U3

  SAVE 

  !-----------------------------------------------------------------------
  !     E X T E R N A L S
  !-----------------------------------------------------------------------

  ! *** Parametrization of the domain
  EXTERNAL PARX,PARY,PARZ

  ! *** Coefficient of exact solution
  EXTERNAL UE,PE,UEX,UEY,UEZ

  ! *** definition of finite elements
  EXTERNAL E030,E031,EM30,EM31,E010

  ! *** Control of refinement - here: regular refinement
  EXTERNAL SEDB,SADB

  EXTERNAL ZVALUE1

END MODULE def_FEAT
