************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: J. Harig, P.Schreiber, S.Turek                              *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
* based on the 2-D routine, modificated by P.Schreiber                 *
* XMABM7                                                               *
*                                                                      *
* Purpose  Calculate matrices corresponding to bilinear forms          *
*          (multigrid version)                                         *
*          Successive call of XABM7                                    *
*                                                                      *
* Subroutines/functions called  XABM7                                  *
*                                                                      *
* Version from  21/02/95                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEv     I*4    Number of levels                                     *
*                 for which matrices are generated                     *
* KLA      I*4    NBLOC*NLEV numbers of arrays                         *
* KLCOL    I*4    NLEV numbers of pointer vectors                      * 
* KLLD     I*4                                                         *
* NBLOC    I*4                                                         *
* ICLEAR   I*4                                                         *
* ELE      SUBR                                                        *
* COEFF    SUBR                                                        *
* BCON     L*4                                                         *
* KAB      I*4    See parameters of XABM7 and ABM7                     *
* KABN     I*4                                                         *
* ICUB     I*4                                                         *
* ISYMM    I*4                                                         *
* CARR     C*6    Names of matrices for all levels                     *
* BSNGL    L*4    NBLOC*NLEV logical values                            *
*                 .TRUE. meand that corresponding array is converted   *
*                 to single precision after completion                 *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KLA      I*4    NBLOC*NLEV numbers of arrays on DA                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMABM7(KLA,KLCOL,KLLD,KNA,KNEQ,NBLCA,ICLR,ELE,
     *                  COEFF,BCON,KAB,KABN,ICUB,ISYMM,ILINT,
     *                  BSNGL,CARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CARR*6
C
      PARAMETER (NNARR=299,NNAB=21,NNLEV=9)
      DIMENSION KLA(NBLCA,*),KLCOL(*),KLLD(*)
      DIMENSION KNA(*),KNEQ(*),KAB(2,NNAB,*),KABN(*),BCON(*)
      DIMENSION CARR(NBLCA,*),BSNGL(NBLCA,*)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL ELE,COEFF
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMABM7'
      IF (ICHECK.GE.997) CALL OTRC('XMABM7','21/02/95')
      IER=0
C
C
      DO 10 ILEV=NLMIN,NLMAX
C
      NEL =KNEL(ILEV)    
      NVT =KNVT(ILEV)
      NET =KNET(ILEV)
      NAT =KNAT(ILEV)
      NVE =KNVE(ILEV)
      NEE =KNEE(ILEV)
      NAE =KNAE(ILEV)
      NVEL =KNVEL(ILEV)
      NEEL =KNEEL(ILEV)
      NVED =KNVED(ILEV)
      NVAR =KNVAR(ILEV)
      NEAR =KNEAR(ILEV)
      NBCT =KNBCT(ILEV)
      NVBD =KNVBD(ILEV)
      NEBD =KNEBD(ILEV)
      NABD =KNABD(ILEV)
C
      LCORVG =KLCVG(ILEV) 
      LCORMG =KLCMG(ILEV) 
      LCORAG =KLCAG(ILEV) 
      LVERT  =KLVERT(ILEV)
      LEDGE  =KLEDGE(ILEV) 
      LAREA  =KLAREA(ILEV) 
      LADJ   =KLADJ(ILEV) 
      LVEL   =KLVEL(ILEV) 
      LEEL   =KLEEL(ILEV)
      LAEL   =KLAEL(ILEV)   
      LVED   =KLVED(ILEV) 
      LAED   =KLAED(ILEV) 
      LVAR   =KLVAR(ILEV) 
      LEAR   =KLEAR(ILEV) 
      LEVE   =KLEVE(ILEV) 
      LAVE   =KLAVE(ILEV) 
      LNPR   =KLNPR(ILEV) 
      LBCT   =KLBCT(ILEV) 
      LVBD   =KLVBD(ILEV) 
      LEBD   =KLEBD(ILEV) 
      LABD   =KLABD(ILEV) 
C
      CALL XABM7(KLA(1,ILEV),KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),
     *           KNEQ(ILEV),NBLCA,ICLR,ELE,COEFF,BCON,
     *           KAB,KABN,ICUB,ISYMM,ILINT,
     *           BSNGL,CARR(1,ILEV))
      IF (IER.NE.0) GOTO 99999
C
10    CONTINUE
C     
99999 END
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: J. Harig, P.Schreiber, S. Turek                             *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XABM7                                                                *
*                                                                      *
* Purpose  Allocation of the blocks in vector DA , KOFF and COECON     *
*          Call of ABM7                                                *
*                                                                      *
* Subroutines/functions called  ABM7, ZNEW, ZDISP, ZLEN, ZTYPE, ZCLEAR *
*                                                                      *
* Version from  21/02/95                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* LCOL     I*4    Numbers of the pointer arrays KCOL and KLD           *
* LLD      I*4    calculated by AP7                                    *
* ICLEAR   I*4    ICLEAR=1  all matrices are cleared                   *
* BSNGL    LOG    =.TRUE.  conversion to single precision              *
* ARR      C*6    Names of blocks (for error messages only)            *
* For the description of the remaining parameters see ABM7             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LA       I*4    Vector of numbers of arrays on DA(VA)                *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NA            *
*                                                                      *
************************************************************************
*                                                                      *
* ABM7                                                                 *
*                                                                      *
* Purpose  Calculation of NBLOC matrices corresponding to              *
*          bilinear forms  a(v,u)                                      *
*          Nonparametric quadrilateral elements                        *
*          Same trial and test functions                               *
*          Storage technique 7/8                                       *
*          Each matrix makes use of the same pointer vectors           *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB3H                       *
*                                                                      *
* Version from  21/02/95                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL     I*4    Pointer vectors for each matrix on DA corresponding  *
* KLD      I*4    to storage technique 7/8 (calculated by AP7)         *
* NA       I*4    Number of elements per matrix                        *
* NEQ      I*4    Number of equations per matrix                       *
* NBLOC    I*4    Number of matrices stored on DA                      *
* KOFF     I*4    Matrix IBLOC starts at position KOFF(IBLOC)+1 on DA  *
* KVERT    I*4    Arrays describing the triangulation                  *
* KEDGE    I*4                                                         *
* KAREA    I*4                                                         *
* DCORVG   R*8                                                         *
* ELE      SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the bilinear form*
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 matrix IBLOC                                         *
* KAB      I*4    Pairs of multiindices occuring in the bilinear forms *
*                 specified separately for each matrix                 *
* KABN     I*4    Number of additive terms in each bilinear form       *
* ICUB     I*4    Number of cubature formula in CB3H                   *
* ISYMM    I*4    0  storage technique 7                               *
*                 1  storage technique 8                               *
*                 2  storage technique 7 but for symmetric matrices    *
* ILINT    I*4    0  full trilinear transformation                     *
*                 1  linear transformation only                        *
*                 2  axiparallel grid                                  *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DA       R*8    Calculated matrices                                  *
* IER      I*4    Error indicator                                      *
*                 -116 Wrong value in array KAB                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XABM7(LA,LCOL,LLD,NA,NEQ,NBLOC,ICLEAR,ELE,
     *                 COEFF,BCON,KAB,KABN,ICUB,ISYMM,ILINT,BSNGL,ARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      PARAMETER (NNARR=299,NNAB=21,NNDER=10)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION LA(*),KAB(2,NNAB,*),KABN(*),BCON(*),ARR(*),BSNGL(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL COEFF,ELE
      SAVE /TRIAA/,/OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XABM7'
      IF (ICHECK.GE.997) CALL OTRC('XABM7 ','21/02/95')
      IER=0
C
      DO 1 IBLOC=1,NBLOC
      IF (LA(IBLOC).EQ.0) THEN
       CALL ZNEW(NA,1,LA(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ELSE
C ***  Check input parameter
       CALL ZTYPE(LA(IBLOC),ITYPE)
       IF (ITYPE.NE.1) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-114,'XAB07 ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LA(IBLOC),ILEN)
       IF (ILEN.LT.NA) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XAB07 ')
        GOTO 99999
       ENDIF
       IF (ICLEAR.EQ.1) CALL ZCLEAR(LA(IBLOC),ARR(IBLOC))
      ENDIF
1     CONTINUE
C
      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NBLOC*NNDER**2,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999
C
      DO 3 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
3     CONTINUE
C
      CALL ABM7(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,
     *          NBLOC,KWORK(L(LOFF)),KWORK(L(LVERT)),KWORK(L(LEDGE)),
     *          KWORK(L(LAREA)),DWORK(L(LCORVG)),ELE,COEFF,BCON,
     *          DWORK(L(LOECON)),KAB,KABN,ICUB,ISYMM,ILINT)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')
C
      DO 4 IBLOC=NBLOC,1,-1
      IF (BSNGL(IBLOC)) THEN
       CALL ZCTYPE(2,LA(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ENDIF
4     CONTINUE
C
99999 END
C
C
C
      SUBROUTINE ABM7(DA,KCOL,KLD,NA,NEQ,NBLOC,KOFF,KVERT,KEDGE,KAREA,
     *                DCORVG,ELE,COEFF,BCON,COECON,KAB,KABN,ICUB,ISYMM,
     *                ILINT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNAB=21,NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
      DIMENSION KCOL(*),KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*)
      DIMENSION KOFF(*),KLD(*),KDFG(NNBAS),KDFL(NNBAS),DA(*),DCORVG(3,*)
      DIMENSION BCON(*),KAB(2,NNAB,*),KABN(*),COECON(NNDER,NNDER,*)
      DIMENSION DB(NNDIM),KENTRY(NNBAS,NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      EXTERNAL ELE,COEFF
      SAVE
C
      SUB='ABM7'
      IF (ICHECK.GE.997) CALL OTRC('ABM7  ','21/02/95')
C
C *** Preparation - evaluation of parameters
      IER=0
      NDIM=1
      BSYMM=ISYMM.GE.1
C *** Which derivatives of basis functions are needed?
      DO 1 I=1,NNDER
1     BDER(I)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 I=1,KABN(IBLOC)
      DO 3 J=1,2
      I1=KAB(J,I,IBLOC)
      IF (I1.LE.0.OR.I1.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-116,'AB07  ')
       GOTO 99999
      ENDIF
3     BDER(I1)=.TRUE.
2     CONTINUE
C *** Dummy call of ELE sets number of element
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,0D0,-1,-1,0,BFIRST)
      BCON0=.TRUE.
C
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 I=1,KABN(IBLOC)
       IA=KAB(1,I,IBLOC)
       IB=KAB(2,I,IBLOC)
5      COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,0D0,IA,IB,IBLOC,BFIRST)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
C
C *** Set zero elements of Jacobian for axiparallel grid
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
C *** Loop over all elements
      JDOFE1=1
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLD(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      IF (BSYMM) JDOFE1=JDOFE
      JCOL0=ILD
      DO 111 IDOFE=JDOFE1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOL(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
120   CONTINUE
C
      IF (ILINT.EQ.2) THEN
       DJ11=(DX(2)+DX(4))*Q2
       DJ12=(DY(2)+DY(4))*Q2
       DJ13=(DZ(1)+DZ(5))*Q2
       DJAC(1,1)=(-DX(1)+DX(2))*Q2
       DJAC(2,1)=(-DY(1)+DY(2))*Q2
       DJAC(1,2)=(-DX(1)+DX(4))*Q2
       DJAC(2,2)=(-DY(1)+DY(4))*Q2
       DJAC(3,3)=(-DZ(1)+DZ(5))*Q2
       DETJ=DJAC(3,3)*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2))
      ELSE IF (ILINT.EQ.1) THEN
       DJ11=(DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=(DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=(DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,1)=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJAC(2,1)=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJAC(3,1)=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJAC(1,2)=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,2)=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,2)=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,3)=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,3)=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,3)=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ELSE
       DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
       DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
      ENDIF
C
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
C *** Cubature points on the reference element
      DO 200 ICUBP=1,NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the trilinear mapping onto the reference element
      IF (ILINT.EQ.0) THEN
       DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
       DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
       DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
       DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
       DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
       DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
       DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
       DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
       DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ENDIF
C
C *** Cubature points on the actual element + weights
      IF (ILINT.EQ.2) THEN
       XX=DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2
       YY=DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2
       ZZ=DJ13+DJAC(3,3)*XI3
      ELSE IF (ILINT.EQ.1) THEN
       XX=DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2+DJAC(1,3)*XI3
       YY=DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2+DJAC(2,3)*XI3
       ZZ=DJ13+DJAC(3,1)*XI1+DJAC(3,2)*XI2+DJAC(3,3)*XI3
      ELSE
       XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
       YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
       ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
      ENDIF
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XX,YY,ZZ,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IALBET=1,KABN(IBLOC)
      IA=KAB(1,IALBET,IBLOC)
      IB=KAB(2,IALBET,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,ZZ,IA,IB,IBLOC,BFIRST)*OM
      ELSE
       AUX=COECON(IA,IB,IBLOC)*OM
      ENDIF
C
      DO 302 JDOFE=1,IDFL
      DO 303 IDIM=1,NDIM
303   DB(IDIM)=DBAS(IDIM,KDFL(JDOFE),IA)
      IF (BSYMM) JDOFE1=JDOFE
      DO 302 IDOFE=JDOFE1,IDFL
      JCOLB=KENTRY(JDOFE,IDOFE)+KOFF(IBLOC)
      DO 304 IDIM=1,NDIM
304   DA(JCOLB)=DA(JCOLB)+DB(IDIM)*DBAS(IDIM,KDFL(IDOFE),IB)*AUX
302   CONTINUE
301   CONTINUE
C
      BFIRST=.FALSE.
300   CONTINUE
C
200   CONTINUE
100   CONTINUE
C
      IF (ISYMM.GT.1) THEN
       DO 400 IBLOC=1,NBLOC
C
       DO 401 IROW=1,NEQ
       DO 402 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
       J1=KCOL(ICOL)
       IF (J1.GE.IROW) GOTO 401
       DO 403 ICOL1=KLD(J1+1)-1,KLD(J1),-1
       IF (KCOL(ICOL1).EQ.IROW) GOTO 404
403    CONTINUE
       GOTO 402
404    ICOLB=ICOL+KOFF(IBLOC)
       ICOL1B=ICOL1+KOFF(IBLOC)
C ***  Elements ICOL and ICOL1 in block IBLOC
       DA(ICOLB)=DA(ICOL1B)
402    CONTINUE
401    CONTINUE
C
400    CONTINUE
      ENDIF
C
99999 END
C
C
C
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Turek                                 *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XVB0                                                                 *
*                                                                      *
* Purpose  Allocation of the blocks in vector DB on DWORK              *
*          Call of VB0                                                 *
*                                                                      *
* Subroutines/functions called  VB0, ZNEW, ZDISP, ZLEN, ZTYPE, ZCLEAR  *
*                                                                      *
* Version from  10/01/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LB       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LB(IBLOC)=0    *
* NEQ      I*4    Length of each vector                                *
* ICLEAR   I*4    ICLEAR=1  all vectors are cleared                    *
* BSNGL    LOG                                                         *
* ARR      C*6    Names of blocks (for error messages only)            *
* For the description of the remaining parameters see VB0              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LB       I*4    Vector of numbers of arrays on DB                    *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NEQ           *
*                                                                      *
************************************************************************
*                                                                      *
* VB0                                                                  *
*                                                                      *
* Purpose  Calculation of NBLOC vectors corresponding to               *
*          linear forms  l(u)                                          *
*          Quadrilateral elements - trilinear transformation           *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB3H                       *
*                                                                      *
* Version from  10/01/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NBLOC    I*4    Number of vectors stored on DB                       *
* KOFF     I*4    Vector IBLOC starts at position KOFF(IBLOC)+1 on DB  *
* KVERT    I*4                                                         *
* KEDGE    I*4                                                         *
* KAREA    I*4                                                         *
* DCORVG   R*8    Arrays describing the triangulation                  *
* ELE      SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the linear forms *
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 vector IBLOC                                         *
* KB       I*4    Multiindices occuring in the linear forms            *
*                 specified separately for each vector                 *
* KBN      I*4    Number of additive terms in each linear form         *
* ICUB     I*4    Number of cubature formula in CB3H                   *
* ILINT    I*4    0  full trilinear transformation                     *
*                 1  linear transformation only                        *
*                 2  axiparallel grid                                  *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DB       R*8    Calculated vectors                                   *
* IER      I*4    Error indicator                                      *
*                 -117 Wrong value in array KB                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE XVBM0(LB,NEQ,NBLOC,ICLEAR,ELE,COEFF,BCON,KB,KBN,
     *                 ICUB,ILINT,BSNGL,ARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      PARAMETER (NNARR=299,NNAB=21,NNDER=10)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION LB(*),KB(NNAB,*),KBN(*),BCON(*),ARR(*),BSNGL(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL COEFF,ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAA/
C
      SUB='XVBM0'
      IF (ICHECK.GE.997) CALL OTRC('XVBM0 ','10/01/95')
      IER=0
C
      DO 1 IBLOC=1,NBLOC
      IF (LB(IBLOC).EQ.0) THEN
       CALL ZNEW(NEQ,1,LB(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ELSE
C ***  Check input parameter
       CALL ZTYPE(LB(IBLOC),ITYPE)
       IF (ITYPE.NE.1) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-114,'XVBM0 ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LB(IBLOC),ILEN)
       IF (ILEN.LT.NEQ) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XVBM0 ')
        GOTO 99999
       ENDIF
       IF (ICLEAR.EQ.1) CALL ZCLEAR(LB(IBLOC),ARR(IBLOC))
      ENDIF
1     CONTINUE
C
      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NBLOC*NNDER,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999
C
      DO 2 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LB(IBLOC))-1
2     CONTINUE
C
      CALL VBM0(DWORK(1),NBLOC,KWORK(L(LOFF)),KWORK(L(LVERT)),
     *          KWORK(L(LEDGE)),KWORK(L(LAREA)),DWORK(L(LCORVG)),
     *          ELE,COEFF,BCON,DWORK(L(LOECON)),KB,KBN,ICUB,ILINT)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')
C
      DO 3 IBLOC=NBLOC,1,-1
      IF (BSNGL(IBLOC)) THEN
       CALL ZCTYPE(2,LB(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ENDIF
3     CONTINUE
C
99999 END
C
C
C
      SUBROUTINE VBM0(DB,NBLOC,KOFF,KVERT,KEDGE,KAREA,DCORVG,
     *                ELE,COEFF,BCON,COECON,KB,KBN,ICUB,ILINT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNAB=21,NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
      DIMENSION KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*),KOFF(*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS),DB(*),DCORVG(3,*)
      DIMENSION BCON(*),KB(NNAB,*),KBN(*),COECON(NNDER,*)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      EXTERNAL ELE,COEFF
      SAVE
C
      SUB='VBM0'
      IF (ICHECK.GE.997) CALL OTRC('VBM0  ','28/02/95')
C
C *** Preparation - evaluation of parameters
      IER=0
C
C *** Which derivatives of basis functions are needed?
      DO 1 IDER=1,NNDER
1     BDER(IDER)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 IBN=1,KBN(IBLOC)
      IB=KB(IBN,IBLOC)
      IF (IB.LE.0.OR.IB.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-117,'VBM0  ')
       GOTO 99999
      ENDIF
3     BDER(IB)=.TRUE.
2     CONTINUE
C
C *** Dummy call of ELE sets number of element
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,0D0,-1,0,BFIRST)
C
      BCON0=.TRUE.
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 IBN=1,KBN(IBLOC)
       IB=KB(IBN,IBLOC)
5      COECON(IB,IBLOC)=COEFF(0D0,0D0,0D0,IB,IBLOC,BFIRST)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
************************************************************************
C *** Calculation of the linear form
************************************************************************
C
C *** Set zero elements of Jacobian for axiparallel grid
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices
      DO 110 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
110   CONTINUE
C
      IF (ILINT.EQ.2) THEN
       DJ11=(DX(2)+DX(4))*Q2
       DJ12=(DY(2)+DY(4))*Q2
       DJ13=(DZ(1)+DZ(5))*Q2
       DJAC(1,1)=(-DX(1)+DX(2))*Q2
       DJAC(2,1)=(-DY(1)+DY(2))*Q2
       DJAC(1,2)=(-DX(1)+DX(4))*Q2
       DJAC(2,2)=(-DY(1)+DY(4))*Q2
       DJAC(3,3)=(-DZ(1)+DZ(5))*Q2
       DETJ=DJAC(3,3)*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2))
      ELSE IF (ILINT.EQ.1) THEN
       DJ11=(DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=(DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=(DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,1)=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJAC(2,1)=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJAC(3,1)=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJAC(1,2)=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,2)=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,2)=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,3)=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,3)=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,3)=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ELSE
       DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
       DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
      ENDIF
C
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
C *** Cubature points on the reference element
      DO 200 ICUBP=1,NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the trilinear mapping onto the reference element
      IF (ILINT.EQ.0) THEN
       DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
       DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
       DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
       DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
       DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
       DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
       DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
       DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
       DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ENDIF
C
C *** Cubature points on the actual element + weights
      IF (ILINT.EQ.2) THEN
       XX=DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2
       YY=DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2
       ZZ=DJ13+DJAC(3,3)*XI3
      ELSE IF (ILINT.EQ.1) THEN
       XX=DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2+DJAC(1,3)*XI3
       YY=DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2+DJAC(2,3)*XI3
       ZZ=DJ13+DJAC(3,1)*XI1+DJAC(3,2)*XI2+DJAC(3,3)*XI3
      ELSE
       XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
       YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
       ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
      ENDIF
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XX,YY,ZZ,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IBN=1,KBN(IBLOC)
      IB=KB(IBN,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,ZZ,IB,IBLOC,BFIRST)*OM
      ELSE
       AUX=COECON(IB,IBLOC)*OM
      ENDIF
C
      DO 302 JDOFE=1,IDFL
      IGLOB=KDFG(JDOFE)+KOFF(IBLOC)
      DO 303 IDIM=1,NDIM
303   DB(IGLOB)=DB(IGLOB)+DBAS(IDIM,KDFL(JDOFE),IB)*AUX
302   CONTINUE
301   CONTINUE
C
      BFIRST=.FALSE.
300   CONTINUE
C
200   CONTINUE
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XMAB09                                                               *
*                                                                      *
* Purpose  Calculate matrices corresponding to bilinear forms          *
*          (multigrid version)                                         *
*          Successive call of XAB09                                    *
*                                                                      *
* Subroutines/functions called  XAB09                                  *
*                                                                      *
* Version from  07/21/92                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Number of levels                                     *
*                 for which matrices are generated                     *
* KLA      I*4    NBLOC*NLEV numbers of arrays                         *
* KLCOL    I*4    NLEV numbers of pointer vectors                      * 
* KLLD     I*4                                                         *
* NBLOC    I*4                                                         *
* ICLEAR   I*4                                                         *
* ELE1     SUBR                                                        *
* ELE2     SUBR                                                        *
* ELE3     SUBR                                                        *
* COEFF    SUBR                                                        *
* BCON     L*4                                                         *
* KAB      I*4    See parameters of XAB07 and AB07                     *
* KABN     I*4                                                         *
* ICUB     I*4                                                         *
* CARR     C*6    Names of matrices for all levels                     *
* BSNGL    L*4    NBLOC*NLEV logical values                            *
*                 .TRUE. meand that corresponding array is converted   *
*                 to single precision after completion                 *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KLA      I*4    NBLOC*NLEV numbers of arrays on DA                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMAB09(KLA,KLCOL,KLLD,KNA,KNEQ,NBLOC,ICLEAR,
     *                  ELE1,ELE2,ELE3,
     *                  COEFF,BCON,KAB,KABN,ICUB,ILINT,BSNGL,CARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CARR*6
C
      PARAMETER (NNARR=299,NNAB=21,NNLEV=9)
      DIMENSION KLA(NBLOC,*),KLCOL(*),KLLD(*),KNEQ(*)
      DIMENSION KNA(*),KAB(3,NNAB,*),KABN(*),BCON(*)
      DIMENSION CARR(NBLOC,*),BSNGL(NBLOC,*)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL ELE1,ELE2,ELE3,COEFF
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMAB09'
      IF (ICHECK.GE.997) CALL OTRC('XMAB09','07/21/92')
      IER=0
C
C
      DO 10 ILEV=NLMIN,NLMAX
C
      NEL =KNEL(ILEV)    
      NVT =KNVT(ILEV)
      NET =KNET(ILEV)
      NAT =KNAT(ILEV)
      NVE =KNVE(ILEV)
      NEE =KNEE(ILEV)
      NAE =KNAE(ILEV)
      NVEL =KNVEL(ILEV)
      NEEL =KNEEL(ILEV)
      NVED =KNVED(ILEV)
      NVAR =KNVAR(ILEV)
      NEAR =KNEAR(ILEV)
      NBCT =KNBCT(ILEV)
      NVBD =KNVBD(ILEV)
      NEBD =KNEBD(ILEV)
      NABD =KNABD(ILEV)
C
      LCORVG =KLCVG(ILEV) 
      LCORMG =KLCMG(ILEV) 
      LCORAG =KLCAG(ILEV) 
      LVERT  =KLVERT(ILEV)
      LEDGE  =KLEDGE(ILEV) 
      LAREA  =KLAREA(ILEV) 
      LADJ   =KLADJ(ILEV) 
      LVEL   =KLVEL(ILEV) 
      LEEL   =KLEEL(ILEV)
      LAEL   =KLAEL(ILEV)   
      LVED   =KLVED(ILEV) 
      LAED   =KLAED(ILEV) 
      LVAR   =KLVAR(ILEV) 
      LEAR   =KLEAR(ILEV) 
      LEVE   =KLEVE(ILEV) 
      LAVE   =KLAVE(ILEV) 
      LNPR   =KLNPR(ILEV) 
      LBCT   =KLBCT(ILEV) 
      LVBD   =KLVBD(ILEV) 
      LEBD   =KLEBD(ILEV) 
      LABD   =KLABD(ILEV)
      CALL XAB09(KLA(1,ILEV),KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),
     *           KNEQ(ILEV),NBLOC,ICLEAR,ELE1,ELE2,ELE3,COEFF,
     *           BCON,
     *           KAB,KABN,ICUB,ILINT,BSNGL(1,ILEV),
     *           CARR(1,ILEV))
      IF (IER.NE.0) GOTO 99999
C
10    CONTINUE
C     
99999 END
c
c
c
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XMAP9                                                                *
*                                                                      *
* Purpose  Calculate pointer vectors                                   *
*          (multigrid version)                                         *
*          Successive call of XAP9                                     *
*                                                                      *
* Subroutines/functions called  XAP9                                   *
*                                                                      *
* Version from  07/21/92                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Number of levels                                     *
* ELE      SUBR                                                        *
* ISYMM    I*4                                                         *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KNA      I*4    Number of entries for each level                     *
* KNEQ     I*4    Number of unknowns for each level                    *
* IER      I*4    Error indicator                                      *
*                 Set by ZNEW                                          *
*                                                                      *
************************************************************************
      SUBROUTINE XMAP9(KLCOL,KLLD,KNA,KNEQ,ELE1,ELE2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION KLCOL(*),KLLD(*),KNA(*),KNEQ(*)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMAP9 '
      IF (ICHECK.GE.997) CALL OTRC('XMAP9 ','07/21/92')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
      NEL =KNEL(ILEV)    
      NVT =KNVT(ILEV)
      NET =KNET(ILEV)
      NAT =KNAT(ILEV)
      NVE =KNVE(ILEV)
      NEE =KNEE(ILEV)
      NAE =KNAE(ILEV)
      NVEL =KNVEL(ILEV)
      NEEL =KNEEL(ILEV)
      NVED =KNVED(ILEV)
      NVAR =KNVAR(ILEV)
      NEAR =KNEAR(ILEV)
      NBCT =KNBCT(ILEV)
      NVBD =KNVBD(ILEV)
      NEBD =KNEBD(ILEV)
      NABD =KNABD(ILEV)
C
      LCORVG =KLCVG(ILEV) 
      LCORMG =KLCMG(ILEV) 
      LCORAG =KLCAG(ILEV) 
      LVERT  =KLVERT(ILEV)
      LEDGE  =KLEDGE(ILEV) 
      LAREA  =KLAREA(ILEV) 
      LADJ   =KLADJ(ILEV) 
      LVEL   =KLVEL(ILEV) 
      LEEL   =KLEEL(ILEV)
      LAEL   =KLAEL(ILEV)   
      LVED   =KLVED(ILEV) 
      LAED   =KLAED(ILEV) 
      LVAR   =KLVAR(ILEV) 
      LEAR   =KLEAR(ILEV) 
      LEVE   =KLEVE(ILEV) 
      LAVE   =KLAVE(ILEV) 
      LNPR   =KLNPR(ILEV) 
      LBCT   =KLBCT(ILEV) 
      LVBD   =KLVBD(ILEV) 
      LEBD   =KLEBD(ILEV) 
      LABD   =KLABD(ILEV)
C
      CALL XAP9(KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),KNEQ(ILEV),ELE1,ELE2,
     *          6)
      IF (IER.NE.0) GOTO 99999
C     
10    CONTINUE
C     
99999 END
c
c
c
************************************************************************
*                                                                      *
* XMSB3                                                                *
*                                                                      *
* Purpose  Generate sequence of meshes for multigrid applications      *
*          by successive calls of XSB0X  (two-level ordering)          *
*                                                                      *
* Subroutines/functions called  XSB0X, WERR, ZCPY                      *
*                                                                      *
* Version from  05/11/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    desired number of l                                  *
* ISCAD    I*4    =1 Determine array KADJ from coarse grid             *
* ISE      I*4    =1 Determine number of midpoints                     *
* ISA      I*4    =0 Release array KADJ on return                      *
*                    after determination of the new subdivisions       *
* ISEEL    I*4    =1 Determine numbers of elements meeting at each     *
*                    edge                                              *
* ISAEL    I*4    =1 Determine numbers of elements meeting at each     *
*                    area                                              *
* ISVEL    I*4    =1 Determine numbers of elements meeting at each     *
*                    vertex                                            *
* IDISP    I*4    =1 Release free space on all arrays after completion *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8   Cartesian coordinates of vertices                     *                                                                      *                                                                      *
************************************************************************
C
      SUBROUTINE XMSB3(ISCAD,ISE,ISA,ISVEL,ISEEL,ISAEL,ISVED,ISAED,
     *                 ISVAR,ISEAR,ISEVE,ISAVE,ISVBD,ISEBD,ISABD,
     *                 IDISP,PARX,PARY,PARZ,SEDB,SADB)
C
      USE Parametrization, ONLY: InitParametrization,ParametrizeBndr
      USE PP3D_MPI, ONLY : myid,myMPI_Barrier
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL SEDB,SADB,PARX,PARY,PARZ
      SAVE
C
      SUB='XMSB3 '
      IF (ICHECK.GE.997) CALL OTRC('XMSB3 ','11/05/91')
C
      NLEV=NLMAX
C
      IF (NLEV.GT.NNLEV) THEN
       CALL WERR(-180,'XMSB2 ')
       GOTO 99999
      ENDIF
C
************************************************************************
C
      DO 10 ILEV=1,NLEV
      write(*,*)'level:',ilev, 'NET:',NET
      IF (ILEV.EQ.1) THEN
      CALL XSB0X(0,1,MAX(1,ISE),ISA,ISVEL,ISEEL,ISAEL,
     *           ISVED,ISAED,ISVAR,ISEAR,ISEVE,ISAVE,
     *           ISVBD,ISEBD,ISABD,IDISP,PARX,PARY,PARZ,
     *           SEDB,SADB)
      ELSE
      CALL XSB0X(1,0,MAX(1,ISE),ISA,ISVEL,ISEEL,ISAEL,
     *           ISVED,ISAED,ISVAR,ISEAR,ISEVE,ISAVE,
     *           ISVBD,ISEBD,ISABD,IDISP,PARX,PARY,PARZ,
     *           SEDB,SADB)
      ENDIF
C
      IF (ILEV.GE.NLMIN) THEN
C
C *** Save dimensions for all levels between NLMIN and NLMAX
C
      KNEL(ILEV)  =NEL
      KNVT(ILEV)  =NVT
      KNET(ILEV)  =NET
      KNAT(ILEV)  =NAT
!      KNVE(ILEV)  =NVE
!      KNEE(ILEV)  =NEE
!      KNAE(ILEV)  =NAE
!      KNVEL(ILEV) =NVEL
!      KNEEL(ILEV) =NEEL
!      KNVED(ILEV) =NVED
!      KNVAR(ILEV) =NVAR
!      KNEAR(ILEV) =NEAR
!      KNBCT(ILEV) =NBCT
!      KNVBD(ILEV) =NVBD
!      KNEBD(ILEV) =NEBD
!      KNABD(ILEV) =NABD
C
C
C *** Save arrays for all levels between NLMIN and NLMAX
C
C
!      IF (ILEV.LT.NLEV) THEN
!C
!C       CALL ZCPY(LCORVG,'DCORVG',KLCVG(ILEV),'DCVG0 ')
!C       IF (IER.NE.0) GOTO 99999
!       IF (LCORMG.NE.0) THEN
!        CALL ZCPY(LCORMG,'DCORMG',KLCMG(ILEV),'DCMG0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LCORAG.NE.0) THEN
!        CALL ZCPY(LCORAG,'DCORAG',KLCAG(ILEV),'DCAG0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       CALL ZCPY(LVERT,'KVERT ',KLVERT(ILEV),'KVERT0')
!       IF (IER.NE.0) GOTO 99999
!       IF (LEDGE.NE.0) THEN
!        CALL ZCPY(LEDGE,'KEDGE  ',KLEDGE(ILEV),'KEDG0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LAREA.NE.0) THEN
!        CALL ZCPY(LAREA,'KAREA  ',KLAREA(ILEV),'KARE0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LADJ.NE.0) THEN
!        CALL ZCPY(LADJ,'KADJ  ',KLADJ(ILEV),'KADJ0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LVEL.NE.0) THEN
!        CALL ZCPY(LVEL,'KVEL  ',KLVEL(ILEV),'KVEL0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LEEL.NE.0) THEN
!        CALL ZCPY(LEEL,'KEEL  ',KLEEL(ILEV),'KEEL0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LAEL.NE.0) THEN
!        CALL ZCPY(LAEL,'KAEL  ',KLAEL(ILEV),'KAEL0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LVED.NE.0) THEN
!        CALL ZCPY(LVED,'KVED  ',KLVED(ILEV),'KVED0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LAED.NE.0) THEN
!        CALL ZCPY(LAED,'KAED  ',KLAED(ILEV),'KAED0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LVAR.NE.0) THEN
!        CALL ZCPY(LVAR,'KVAR  ',KLVAR(ILEV),'KVAR0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LEAR.NE.0) THEN
!        CALL ZCPY(LEAR,'KEAR  ',KLEAR(ILEV),'KEAR0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LEVE.NE.0) THEN
!        CALL ZCPY(LEVE,'KEVE  ',KLEVE(ILEV),'KEVE0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LAVE.NE.0) THEN
!        CALL ZCPY(LAVE,'KAVE  ',KLAVE(ILEV),'KAVE0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       CALL ZCPY(LNPR,'KNPR  ',KLNPR(ILEV),'KNPR0 ')
!       IF (IER.NE.0) GOTO 99999
!       IF (LBCT.NE.0) THEN
!        CALL ZCPY(LBCT,'KLBCT ',KLBCT(ILEV),'KBCT0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LVBD.NE.0) THEN
!        CALL ZCPY(LVBD,'KVBD  ',KLVBD(ILEV),'KVBD0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LEBD.NE.0) THEN
!        CALL ZCPY(LEBD,'KEBD  ',KLEBD(ILEV),'KEBD0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF
!       IF (LABD.NE.0) THEN
!        CALL ZCPY(LABD,'KABD  ',KLABD(ILEV),'KABD0 ')
!        IF (IER.NE.0) GOTO 99999
!       ENDIF   
!C
!       ELSE
!C
!       KLCVG(ILEV) =LCORVG
!       KLCMG(ILEV) =LCORMG
!       KLCAG(ILEV) =LCORAG
!       KLVERT(ILEV)=LVERT
!       KLEDGE(ILEV)=LEDGE
!       KLAREA(ILEV)=LAREA
!       KLADJ(ILEV) =LADJ
!       KLVEL(ILEV) =LVEL
!       KLEEL(ILEV) =LEEL
!       KLAEL(ILEV) =LAEL
!       KLVED(ILEV) =LVED
!       KLAED(ILEV) =LAED
!       KLVAR(ILEV) =LVAR
!       KLEAR(ILEV) =LEAR
!       KLEVE(ILEV) =LEVE
!       KLAVE(ILEV) =LAVE
!       KLNPR(ILEV) =LNPR
!       KLBCT(ILEV) =LBCT
!       KLVBD(ILEV) =LVBD
!       KLEBD(ILEV) =LEBD
!       KLABD(ILEV) =LABD
!C     
!      ENDIF
C
      DO 20 ILEVH=1,NLEV-1
!      KLCVG(ILEVH)=LCORVG
20    CONTINUE
C
      ENDIF


      !CALL ZDISP(0,LCORVG,'DCVORG')
      !CALL ZDISP(0,LCORVG,'DCVORG')

C
      IF (ILEV.EQ.1) THEN
!        CALL InitParametrization()
      ELSE
!        CALL ParametrizeBndr()
!       CALL MIDDEF(DWORK(L(KLCVG(ILEV))),KWORK(L(KLVERT(ILEV))),
!     *      KWORK(L(KLAREA(ILEV))),KNEL(ILEV),DWORK(L(KLCAG(ILEV))))
      END IF
C
10    CONTINUE
C
C
************************************************************************
C     Here comes the deformation stuff
!       CALL DEFORM(DWORK(L(LCORVG)),NVT)
C     MIDPOINTS
C
************************************************************************
99999 END
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
* XMORS3                                                                *
*                                                                      *
* Purpose  Read multiple triangulations (multigrid version)           *
*                                                                      *
* Subroutines/functions called  XORS                                   *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* MFILE    I*4    Unit number used for output                          *
* CCFILE   C*(*)  Array of filenames used for output                   *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 Set by OF0 or OWA                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMORS3(IRMESH)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CTRIF(NNLEV)*15
C
      DATA CTRIF/ '#tries/#tria1  ','#tries/#tria2  ','#tries/#tria3  ',
     *            '#tries/#tria4  ','#tries/#tria5  ','#tries/#tria6  ',
     *            '#tries/#tria7  ','#tries/#tria8  ','#tries/#tria9  '/
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
      SUB='XMORS3'
      IF (ICHECK.GE.997) CALL OTRC('XMORS3 ','06/25/95')
      IER=0
C
      NLEV=NLMAX
C
      IF (IRMESH.GT.1) THEN
       IFMT=1
      ELSE
       IFMT=0
      ENDIF
C
      DO 10 ILEV=1,NLEV
C
C *** Read refined mesh in file CTRIF(ILEV) with format IFMT
      CALL XORS(65,CTRIF(ILEV),IFMT)
      REWIND(65)
      CLOSE(65)
C
      KNEL(ILEV) =NEL
      KNVT(ILEV) =NVT
      KNET(ILEV) =NET
      KNAT(ILEV) =NAT
      KNVE(ILEV) =NVE
      KNEE(ILEV) =NEE
      KNAE(ILEV) =NAE
      KNVEL(ILEV) =NVEL
      KNEEL(ILEV) =NEEL
      KNVED(ILEV) =NVED
      KNVAR(ILEV) =NVAR
      KNEAR(ILEV) =NEAR
      KNBCT(ILEV) =NBCT
      KNVBD(ILEV) =NVBD
      KNEBD(ILEV) =NEBD
      KNABD(ILEV) =NABD
C
C
C *** Save arrays for all levels between NLMIN and NLMAX
      IF (ILEV.LT.NLMAX) THEN
C
       IF (LCORMG.NE.0) THEN
        CALL ZCPY(LCORMG,'DCORMG',KLCMG(ILEV),'DCMG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LCORAG.NE.0) THEN
        CALL ZCPY(LCORAG,'DCORAG',KLCAG(ILEV),'DCAG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LVERT.NE.0) THEN
        CALL ZCPY(LVERT,'KVERT ',KLVERT(ILEV),'KVERT0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LEDGE.NE.0) THEN
        CALL ZCPY(LEDGE,'KEDGE  ',KLEDGE(ILEV),'KEDG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LAREA.NE.0) THEN
        CALL ZCPY(LAREA,'KAREA  ',KLAREA(ILEV),'KARE0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LADJ.NE.0) THEN
        CALL ZCPY(LADJ,'KADJ  ',KLADJ(ILEV),'KADJ0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LVEL.NE.0) THEN
        CALL ZCPY(LVEL,'KVEL  ',KLVEL(ILEV),'KVEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LEEL.NE.0) THEN
        CALL ZCPY(LEEL,'KEEL  ',KLEEL(ILEV),'KEEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LAEL.NE.0) THEN
        CALL ZCPY(LAEL,'KAEL  ',KLAEL(ILEV),'KAEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LVED.NE.0) THEN
        CALL ZCPY(LVED,'KVED  ',KLVED(ILEV),'KVED0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LAED.NE.0) THEN
        CALL ZCPY(LAED,'KAED  ',KLAED(ILEV),'KAED0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LVAR.NE.0) THEN
        CALL ZCPY(LVAR,'KVAR  ',KLVAR(ILEV),'KVAR0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LEAR.NE.0) THEN
        CALL ZCPY(LEAR,'KEAR  ',KLEAR(ILEV),'KEAR0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LEVE.NE.0) THEN
        CALL ZCPY(LEVE,'KEVE  ',KLEVE(ILEV),'KEVE0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LAVE.NE.0) THEN
        CALL ZCPY(LAVE,'KAVE  ',KLAVE(ILEV),'KAVE0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF(LNPR.NE.0) THEN
        CALL ZCPY(LNPR,'KNPR  ',KLNPR(ILEV),'KNPR0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LBCT.NE.0) THEN
       CALL ZCPY(LBCT,'KLBCT ',KLBCT(ILEV),'KBCT0 ')
       IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LVBD.NE.0) THEN
        CALL ZCPY(LVBD,'KVBD  ',KLVBD(ILEV),'KVBD0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LEBD.NE.0) THEN
        CALL ZCPY(LEBD,'KEBD  ',KLEBD(ILEV),'KEBD0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LABD.NE.0) THEN
        CALL ZCPY(LABD,'KABD  ',KLABD(ILEV),'KABD0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF   
C
C
      ELSE
C     
C
      KLCVG(ILEV) =LCORVG
      KLCMG(ILEV) =LCORMG
      KLCAG(ILEV) =LCORAG
      KLVERT(ILEV)=LVERT
      KLEDGE(ILEV)=LEDGE
      KLAREA(ILEV)=LAREA
      KLADJ(ILEV) =LADJ
      KLVEL(ILEV) =LVEL
      KLEEL(ILEV) =LEEL
      KLAEL(ILEV) =LAEL  
      KLVED(ILEV) =LVED
      KLAED(ILEV) =LAED 
      KLVAR(ILEV) =LVAR
      KLEAR(ILEV) =LEAR 
      KLEVE(ILEV) =LEVE
      KLAVE(ILEV) =LAVE
      KLNPR(ILEV) =LNPR
      KLBCT(ILEV) =LBCT
      KLVBD(ILEV) =LVBD
      KLEBD(ILEV) =LEBD
      KLABD(ILEV) =LABD 
      ENDIF
C
10    CONTINUE   
C
      DO 20 ILEV=1,NLMAX-1
      KLCVG(ILEV)=LCORVG
20    CONTINUE
C     
99999 END
C
C
C
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.0)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XMORAn                                                               *
*                                                                      *
* Purpose  Get matrix and pointer vectors corresponding to the storage *
*          technique, previously stored by XMOWAn, back on DWORK       *
*          (multigrid version)                                         *
*          Successive call of XORA                                     *
*          Matrix stored in technique  n  (see Reference Manual)       *
*                                                                      *
* Subroutines/functions called  XORA                                   *
*                                                                      *
* Version from  02/18/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* MFILE    I*4    Unit number used for output                          *
* CCFILE   C*(*)  Array of filenames used for output                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 Set by OF0 or ORA                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMORA3(KLA,KLDIA,KLDIAS,KNDIA,MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
      DIMENSION KLA(*),KLDIA(*),KLDIAS(*),KNDIA(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XMORA3'
      IF (ICHECK.GE.997) CALL OTRC('XMORA3','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
C
      LA   =0
      LDIA =0
      LDIAS=0
C
      CALL XORA(LA   ,'DA    ',MFILE,CCFILE(ILEV),IFMT)
      CALL XORA(LDIA ,'KDIA  ',MFILE,CCFILE(ILEV),IFMT)
      CALL XORA(LDIAS,'KDIAS ',MFILE,CCFILE(ILEV),IFMT)
	IF (IFMT.EQ.1) THEN
	 READ(MFILE,'(I6)') NDIA
	ELSE
	 READ(MFILE) NDIA
	ENDIF
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
      KLA(ILEV)   =LA
      KLDIA(ILEV) =LDIA
      KLDIAS(ILEV)=LDIAS
C     
10    CONTINUE   
C     
99999 END
C
C
C
      SUBROUTINE XMORA7(KLA,KLCOL,KLLD,MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
      DIMENSION KLA(*),KLCOL(*),KLLD(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XMORA7'
      IF (ICHECK.GE.997) CALL OTRC('XMORA7','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
C
      LA  =0
      LCOL=0
      LLD =0
C
      CALL XORA(LA  ,'DA    ',MFILE,CCFILE(ILEV),IFMT)
      CALL XORA(LCOL,'KCOL  ',MFILE,CCFILE(ILEV),IFMT)
      CALL XORA(LLD ,'KLD   ',MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
      KLA(ILEV)  =LA
      KLCOL(ILEV)=LCOL
      KLLD(ILEV) =LLD
C     
10    CONTINUE   
C     
99999 END
C
C
C      
      SUBROUTINE XMORA9(KLA,KLCOL,KLLD,NBLOC,MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
      DIMENSION KLA(NBLOC,*),KLCOL(*),KLLD(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XMORA9'
      IF (ICHECK.GE.997) CALL OTRC('XMORA9','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
      DO 20 IBLOC=1,NBLOC
      CALL XORA(LA  ,'DA    ',MFILE,CCFILE(ILEV),IFMT)
      KLA(IBLOC,ILEV)=LA
20    CONTINUE
C
      CALL XORA(LCOL,'KCOL  ',MFILE,CCFILE(ILEV),IFMT)
      CALL XORA(LLD ,'KLD   ',MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
C     
      CLOSE(MFILE)
      KLCOL(ILEV)=LCOL
      KLLD (ILEV)=LLD
10    CONTINUE 
C     
99999 END
c
c
c
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.0)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XMOWAn                                                               *
*                                                                      *
* Purpose  Store matrix and pointer vectors corresponding to the       *
*          storage technique  (multigrid version)                      *
*          Successive call of XOWA                                     *
*          Matrix stored in technique  n  (see Reference Manual)       *
*                                                                      *
* Subroutines/functions called  XOWA                                   *
*                                                                      *
* Version from  02/18/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* MFILE    I*4    Unit number used for output                          *
* CCFILE   C*(*)  Array of filenames used for output                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 Set by OF0 or OWA                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMOWA3(KLA,KLDIA,KLDIAS,KNDIA,MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
      DIMENSION KLA(*),KLDIA(*),KLDIAS(*),KNDIA(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XMOWA3'
      IF (ICHECK.GE.997) CALL OTRC('XMOWA3','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
C
      LA   =KLA(ILEV)
      LDIA =KLDIA(ILEV)
      LDIAS=KLDIAS(ILEV)
	NDIA =KNDIA(ILEV)
C
      CALL XOWA(LA   ,'DA    ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LDIA ,'KDIA  ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LDIAS,'KDIAS ',MFILE,CCFILE(ILEV),IFMT)
	IF (IFMT.EQ.1) THEN
	 WRITE(MFILE,'(I6)') NDIA
	ELSE
	 WRITE(MFILE) NDIA
	ENDIF
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
C     
10    CONTINUE   
C     
99999 END
C
C
C
      SUBROUTINE XMOWA7(KLA,KLCOL,KLLD,MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
      DIMENSION KLA(*),KLCOL(*),KLLD(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XMOWA7'
      IF (ICHECK.GE.997) CALL OTRC('XMOWA7','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
C
      LA  =KLA(ILEV)
      LCOL=KLCOL(ILEV)
      LLD =KLLD(ILEV)
C
      CALL XOWA(LA  ,'DA    ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LCOL,'KCOL  ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LLD ,'KLD   ',MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
C     
10    CONTINUE   
C     
99999 END
C
C
C
      SUBROUTINE XMOWA9(KLA,KLCOL,KLLD,NBLOC,MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
      DIMENSION KLA(NBLOC,*),KLCOL(*),KLLD(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XMOWA9'
      IF (ICHECK.GE.997) CALL OTRC('XMOWA9','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
      DO 20 IBLOC=1,NBLOC
      LA  =KLA(IBLOC,ILEV)
      CALL XOWA(LA  ,'DA    ',MFILE,CCFILE(ILEV),IFMT)
20    CONTINUE
C
      LCOL=KLCOL(ILEV)
      LLD =KLLD(ILEV)
      CALL XOWA(LCOL,'KCOL  ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LLD ,'KLD   ',MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
C     
      CLOSE(MFILE)
10    CONTINUE 
C     
99999 END
C
C
C
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* ORALCn                                                               *
*                                                                      *
* Purpose  Get arrays from file, previously stored by XOWA or OWAn     *
*                                                                      *
* Subroutines/functions called  ORALC0                                 *
*                                                                      *
* Version from  10/11/94                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----    ----                                                        *
* MFILE   I*4    I/O unit - previouly opened (by OF0)                  *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
* A1      R*8    Scaling factor: DX=DX+A1*DX(readin)                   *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* DX      R*8    Double precision array - version ORA1                 *
* VX      R*4    Single precision array - version ORA2                 *
* ARR     C*6    Name of array read from data file                     *
* IER     I*4    Error indicator                                       *
*                -110  Error while reading from unit  MFILE            *
*                -113  Wrong value of  ITYPE  read from  MFILE         *
*                                                                      *
************************************************************************
C
      SUBROUTINE ORALC1(DX,A1,ARR,MFILE,IFMT)
C
C *** Double precision version
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      DIMENSION DX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='ORALC1'
      IF (ICHECK.GE.997) CALL OTRC('ORALC1','10/11/94')
C
      IER=0
C
      IF (IFMT.NE.0) GOTO 99997
C
      READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
C
      IF (ITYPE.NE.1) GOTO 99997
C
      CALL ORALCD(DX,A1,ILEN8,MFILE,JRECL)
      IF (IER.NE.0) GOTO 99999
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORALC1')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORALC')
C
99999 END
C
C
C
      SUBROUTINE ORALC2(VX,A1,ARR,MFILE,IFMT)
C
C *** Single precision version
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      DIMENSION VX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='ORALC2'
      IF (ICHECK.GE.997) CALL OTRC('ORALC2','10/11/94')
C
      IER=0
C
      IF (IFMT.NE.0) GOTO 99997
C
      READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
C
      IF (ITYPE.NE.2) GOTO 99997
C
      CALL ORALCV(VX,A1,ILEN8,MFILE,JRECL)
      IF (IER.NE.0) GOTO 99999
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORALC2')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORALC')
C
99999 END
C
C
C
      SUBROUTINE ORALCD(DX,A1,N,MFILE,IRECL8)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DXH(512)
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /CHAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('ORALCD','10/11/94')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
      READ (MFILE,ERR=99998,END=99998) (DXH(I),I=1,IRECL8)
      DO 1 I=1,IRECL8
      DX(J1+I)=DX(J1+I)+DXH(I)
1     CONTINUE
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      READ (MFILE,ERR=99998,END=99998) (DXH(I),I=1,MOD(N,IRECL8))
      DO 2 I=1,MOD(N,IRECL8)
      DX(J1+I)=DX(J1+I)+DXH(I)
2     CONTINUE
      GOTO 99999
C
99998 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORALCD')
C
99999 END
C
C
C
      SUBROUTINE ORALCV(VX,A1,N,MFILE,IRECL8)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VX(*),VXH(512)
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /CHAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('ORALCV','10/11/94')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
      READ (MFILE,ERR=99998,END=99998) (VXH(I),I=1,IRECL8)
      DO 1 I=1,IRECL8
      VX(J1+I)=VX(J1+I)+VXH(I)
1     CONTINUE
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      READ (MFILE,ERR=99998,END=99998) (VXH(I),I=1,MOD(N,IRECL8))
      DO 2 I=1,MOD(N,IRECL8)
      VX(J1+I)=VX(J1+I)+VXH(I)
2     CONTINUE
      GOTO 99999
C
99998 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORALCV')
C
99999 END
C
C
C
      SUBROUTINE OWA2V(VX,ARR,NX,MFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6,CFORM*15
      DIMENSION VX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='OWA2V'
      IF (ICHECK.GE.997) CALL OTRC('OWA2V ','10/11/94')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       CFORM=FMT(1)
       WRITE (MFILE,'(2A15,2I15)') ARR,CFORM,1,NX
      ELSE
       WRITE (MFILE) ARR,2,NX,NX,IRECL8
      ENDIF
C
      IF (BFMT) THEN
       WRITE (MFILE,CFORM) (VX(I),I=1,NX)
      ELSE
       CALL OWA0V(VX,NX,MFILE,IRECL8)
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(8,'XOWA2V')
C
99999 END
C
C
C
      SUBROUTINE OWA0V(VX,N,MFILE,IRECL8)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('OWA0V ','10/11/94')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
1     WRITE (MFILE) (VX(J1+I),I=1,IRECL8)
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      WRITE (MFILE) (VX(J1+I),I=1,MOD(N,IRECL8))
C
99999 END
C
C
C
      SUBROUTINE XXSB0(NFINE,ISE,ISA,ISVEL,ISEEL,ISAEL,ISVED,ISAED,
     *                ISVAR,ISEAR,ISEVE,ISAVE,ISVBD,ISEBD,ISABD,IDISP,
     *                PARX,PARY,PARZ,SEDB,SADB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNVE=8,NNAE=6,NNEE=12,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAH/  LCORH
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL PARX,PARY,PARZ,SEDB,SADB
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/,/TRIAH/
C
      SUB='XSB0'
      IF (ICHECK.GE.997) CALL OTRC('XSB0  ','12/01/93')
C
      CALL ZLEN(LVERT,ILEN)
      NNEL=ILEN/NNVE
      CALL ZLEN(LADJ,ILEN)
      NNEL1=ILEN/NNAE
      CALL ZLEN(LEDGE,ILEN)
      NNEL2=ILEN/NNEE
      NNEL=MIN0(NNEL,NNEL1,NNEL2)
      CALL ZLEN(LCORVG,ILEN)
      NNVT=ILEN/3
      NNEL1=NNEL
      NNVT1=NNVT
C
      IF (LCORH.NE.0) THEN
       CALL ZLEN(LCORVG,ILEN1)
       CALL ZLEN(LCORH ,ILEN2)
       IF (ILEN1.NE.ILEN2) GOTO 99999
       CALL ZDISP(0,LCORVG,'DCVORG')
       IF (IER.NE.0) GOTO 99999
       LCORVG=LCORH
       LCORH =0
      ENDIF
C
      NVEL=0
      CALL ZLEN(LVEL,ILEN)
      IF (ILEN.GT.0) CALL ZDISP(0,LVEL,'KVEL ')
      CALL ZNEW(NVEL,-3,LVEL,'KVEL  ')
      IF (IER.NE.0) GOTO 99999
C
      CALL SB0(DWORK(L(LCORVG)),DWORK(L(LCORMG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LEDGE)),KWORK(L(LNPR)),
     *         KWORK(L(LVEL)),NFINE,NNEL,NNVT,PARX,PARY,PARZ,
     *         SEDB,SADB)
C
      IF (IER) 99999,1,2
C
1     IF ((ISVEL.GT.0).OR.(ISE.GT.0)) THEN
       CALL SBVEL(KWORK(L(LVERT)),KWORK(L(LVEL)),0)
       CALL SBVEL(KWORK(L(LVERT)),KWORK(L(LVEL)),1)
       CALL ZDISP(NVEL*NVT,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZDISP(0,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
       LVEL=0
      ENDIF
C
      CALL ZNEW(3*NVT,-1,LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99999
      CALL ZCPY (LCORVG,'DCORVG',LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99999
      CALL SVC(DWORK(L(LCORVG)),KWORK(L(LNPR)),NVT,PARX,PARY,PARZ)
      IF (IER.NE.0) GOTO 99999
C
      GOTO 99995
C
C
2     CALL ZDISP(0,LVEL,'KVEL  ')
C
      IF (NNEL.GT.NNEL1) THEN
       CALL ZNEW(NNVE*NNEL,3,LV1,'KVERT ')
       IF (IER.NE.0) GOTO 99999
       IF (LVERT.NE.0) THEN
        CALL ZCPY (LVERT,'KVOLD ',LV1,'KVERT ')
        CALL ZDISP(0,LVERT,'KVOLD ')
       ENDIF
       LVERT=LV1
C
       CALL ZNEW(NNAE*NNEL,3,LA1,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
       IF (LADJ.NE.0) THEN
        CALL ZCPY (LADJ,'KAOLD ',LA1,'KADJ  ')
        CALL ZDISP(0,LADJ,'KAOLD ')
       ENDIF
       LADJ=LA1
C
       CALL ZNEW(NNEE*NNEL,3,LE1,'KEDGE ')
       IF (IER.NE.0) GOTO 99999
       IF (LEDGE.NE.0) THEN
        CALL ZCPY (LEDGE,'KEOLD ',LE1,'KEDGE ')
        CALL ZDISP(0,LEDGE,'KEOLD ')
       ENDIF
       LEDGE=LE1
      ENDIF
C
      IF (NNVT.GT.NNVT1) THEN
       CALL ZNEW(3*NNVT,1,LCV1,'DCORVG')
       IF (IER.NE.0) GOTO 99999
       IF (LCORVG.NE.0) THEN
        CALL ZCPY (LCORVG,'DCVOLD',LCV1,'DCORVG')
        CALL ZDISP(0,LCORVG,'DCVOLD')
       ENDIF
       LCORVG=LCV1
C
       CALL ZNEW(NNVT,3,LNPR1,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
       IF (LNPR.NE.0) THEN
        CALL ZCPY (LNPR,'KNPOLD',LNPR1,'KNPR  ')
        CALL ZDISP(0,LNPR,'KNPOLD')
       ENDIF
       LNPR=LNPR1
       LNPR1=0
      ENDIF
C
      NVEL=0
      CALL ZLEN(LVEL,ILEN)
      IF (ILEN.GT.0) CALL ZDISP(0,LVEL,'KVEL ')
      CALL ZNEW(NVEL,-3,LVEL,'KVEL  ')
      IF (IER.NE.0) GOTO 99999
C
      CALL SB0(DWORK(L(LCORVG)),DWORK(L(LCORMG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LEDGE)),KWORK(L(LNPR)),
     *         KWORK(L(LVEL)),NFINE,NNEL,NNVT,PARX,PARY,PARZ,
     *         SEDB,SADB)
C
C
      IF ((ISVEL.GT.0).OR.(ISE.GT.0)) THEN
       CALL SBVEL(KWORK(L(LVERT)),KWORK(L(LVEL)),0)
       CALL SBVEL(KWORK(L(LVERT)),KWORK(L(LVEL)),1)
       CALL ZDISP(NVEL*NVT,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZDISP(0,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
       LVEL=0
      ENDIF
C
      CALL ZNEW(3*NVT,-1,LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99999
      CALL ZCPY (LCORVG,'DCORVG',LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99999
      CALL SVC(DWORK(L(LCORVG)),KWORK(L(LNPR)),NVT,PARX,PARY,PARZ)
      IF (IER.NE.0) GOTO 99999
C
C
99995 IF (ISE.GT.0) THEN
       IF (ISE.GE.2) THEN
        LNPR1=0
        CALL ZNEW(NVT+8*NET,3,LNPR1,'KNPR1 ')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY (LNPR,'KNPOLD',LNPR1,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
        CALL ZDISP(0,LNPR,'KNPOLD')
        IF (IER.NE.0) GOTO 99999
        LNPR=LNPR1
        LNPR1=0
       ENDIF
C        
       IF (ISE.GE.3) THEN
        IF (LCORMG.NE.0) THEN
         CALL ZDISP(0,LCORMG,'DCORMG')
         IF (IER.NE.0) GOTO 99999
        ENDIF
        CALL ZNEW(3*8*NET,1,LCORMG,'DCOREG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
       CALL SBE(KWORK(L(LVERT)),KWORK(L(LNPR)),KWORK(L(LVEL)),
     *          DWORK(L(LCORH)),KWORK(L(LEDGE)),DWORK(L(LCORMG)),ISE,
     *          PARX,PARY,PARZ,SEDB)
       IF (IER.NE.0) GOTO 99999
C
       IF (ISE.GE.3) THEN
        CALL ZDISP(3*NET,LCORMG,'DCOREG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
       IF (ISE.GE.2) THEN
        CALL ZDISP(NVT+NET,LNPR,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
      ENDIF
C
C
      IF ((ISVEL.EQ.0).AND.(ISE.GT.0)) THEN
       CALL ZDISP(0,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
       LVEL=0
      ENDIF
C
C
      IF ((ISEEL.GT.0).AND.(ISE.GT.0)) THEN
       NEEL=0
       CALL ZLEN(LEEL,ILEN)
       IF (ILEN.GT.0) CALL ZDISP(0,LEEL,'KEEL ')
       CALL ZNEW(NEEL,-3,LEEL,'KEEL  ')
       IF (IER.NE.0) GOTO 99999
       CALL SBEEL(KWORK(L(LEDGE)),KWORK(L(LEEL)),0)
       CALL SBEEL(KWORK(L(LEDGE)),KWORK(L(LEEL)),1)
       CALL ZDISP(NEEL*NET,LEEL,'KEEL  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C
      IF (ISA.GT.0) THEN
       CALL ZLEN(LAREA,ILEN)
       IF (ILEN.GT.0) CALL ZDISP(0,LAREA,'KAREA')
       CALL ZNEW(NNAE*NEL,3,LAREA,'KAREA ')
       IF (IER.NE.0) GOTO 99999
C
       IF (ISA.GE.2) THEN
        LNPR1=0
        IF (ISE.GE.2) THEN
         INDA=NVT+NET
        ELSE
         INDA=NVT
        ENDIF
        CALL ZNEW(INDA+8*NAT,3,LNPR1,'KNPR1 ')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY (LNPR,'KNPOLD',LNPR1,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
        CALL ZDISP(0,LNPR,'KNPOLD')
        IF (IER.NE.0) GOTO 99999
        LNPR=LNPR1
        LNPR1=0
       ENDIF
C        
       IF (ISA.GE.3) THEN
        IF (LCORAG.NE.0) THEN
         CALL ZDISP(0,LCORAG,'DCORAG')
         IF (IER.NE.0) GOTO 99999
        ENDIF
        CALL ZNEW(3*8*NAT,1,LCORAG,'DCORAG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
       CALL SBA(KWORK(L(LVERT)),KWORK(L(LNPR)),KWORK(L(LADJ)),
     *          DWORK(L(LCORH)),KWORK(L(LAREA)),DWORK(L(LCORAG)),
     *          ISA,INDA,PARX,PARY,PARZ,SADB)
       IF (IER.NE.0) GOTO 99999
C
       IF (ISA.GE.3) THEN
        CALL ZDISP(3*NAT,LCORAG,'DCORAG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
       IF (ISA.GE.2) THEN
        CALL ZDISP(INDA+NAT,LNPR,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
      ENDIF
C
C
      IF ((ISAEL.GT.0).AND.(ISA.GT.0)) THEN
       CALL ZLEN(LAEL,ILEN)
       IF (ILEN.GT.0) CALL ZDISP(0,LAEL,'KAEL ')
       CALL ZNEW(2*NAT,3,LAEL,'KAEL  ')
       IF (IER.NE.0) GOTO 99999
       CALL SBAEL(KWORK(L(LAREA)),KWORK(L(LAEL)))
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C
      IF (LCORH.NE.0) THEN
       CALL ZDISP(0,LCORH,'DCVOH ')
       IF (IER.NE.0) GOTO 99999
       LCORH =0
      ENDIF
C
C
C
99999 END
