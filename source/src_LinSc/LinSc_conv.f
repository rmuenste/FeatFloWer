************************************************************************
      SUBROUTINE RhoCpConvMat(U1,U2,U3,DA,NU,KCOLA,KLDA,KVERT,KAREA,
     *           KEDGE,DCORVG,KADJ,KVERT1,KAREA1,KEDGE1,
     *           NEL1,NVT1,NET1,NAT1,ELE)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE var_QuadScalar, ONLY : myBoundary,myHeatObjects
      USE Sigma_User, ONLY : myMaterials,mySigma
      USE PP3D_MPI, ONLY:myid
C      
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION U1(*),U2(*),U3(*),DA(*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*),KADJ(NNAE,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KVERT1(NNVE,*),KAREA1(NNAE,*),KEDGE1(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KBFG(NNBAS),KBFL(NNBAS)
      DIMENSION DPHI(NNBAS)
      INTEGER   AdjElem(8)
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
      EXTERNAL EPHI
C
      SAVE
C
      NA=KLDA(NU+1)-1
C
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,4
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      ICUB=7!ICUBN
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL1=1,NEL1
C
      CALL localNDFGL(IEL1,KVERT1,KEDGE1,KAREA1,
     *     NVT1,NET1,NAT1,KBFG,KBFL)
      IF (IER.LT.0) GOTO 99999
C
      AdjElem(1)=IEL1
      AdjElem(2)=KADJ(3,AdjElem(1))
      AdjElem(3)=KADJ(3,AdjElem(2))
      AdjElem(4)=KADJ(3,AdjElem(3))
      AdjElem(5)=KADJ(6,AdjElem(1))
      AdjElem(6)=KADJ(3,AdjElem(5))
      AdjElem(7)=KADJ(3,AdjElem(6))
      AdjElem(8)=KADJ(3,AdjElem(7))

      DO 150 IEL8=1,8

      IEL = AdjElem(IEL8)

      dAvgCoeff = 0d0
      do j=1,8
       ivt = KVERT(j,iel)
       iSeg = myHeatObjects%Segment(ivt)
       if (iSeg.eq.0) then
        iMat= 0
       else
        iMat = mySigma%mySegment(iSeg)%MatInd
       end if
       dAvgCoeff = dAvgCoeff + 0.125d0*
     * (1d0*myMaterials(iMat)%density)*(1d7*myMaterials(iMat)%cp)
      end do
      
      DDD = dAvgCoeff

      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
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
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the bilinear mapping onto the reference element
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
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DU1=0D0
      DU2=0D0
      DU3=0D0
C
      CALL SETXUI(IEL8,XI1,XI2,XI3,XUI1,XUI2,XUI3)
      CALL EPHI(XUI1,XUI2,XUI3,DPHI)
      DO 220 JDFL=1,NNBAS
        HBAS=DPHI(KBFL(JDFL))
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KBFG(JDFL)
         DU1=DU1+U1(JDFG)*HBAS
         DU2=DU2+U2(JDFG)*HBAS
         DU3=DU3+U3(JDFG)*HBAS
        ENDIF
220   CONTINUE
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(1,JDOFEH,1)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2+HBASJ4*DU3
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=HSUMJ*HBASJ1
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         HSUMI=HBASI2*DU1+HBASI3*DU2+HBASI4*DU3
         AH=HSUMI*HBASJ1
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*DDD*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA)-DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
150   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE Conv_LinSc1(U1,U2,U3,DA,NU,KCOLA,
     *           KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,ELE)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION U1(*),U2(*),U3(*),DA(*)
      DIMENSION KCOLA(*),KLDA(*),KINT(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KBFG(NNBAS),KBFL(NNBAS)
      DIMENSION DPHI(NNBAS)
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
      EXTERNAL EPHI
C
      SAVE
C
      NA=KLDA(NU+1)-1

      CALL LCL1(DA,NA)
C
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,4
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      ICUB=7!ICUBN
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      ILINT=0!KINT(IEL)
C
C *** Set zero elements of Jacobian for axiparallel grid
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      CALL NDFGL(IEL,1,13,KVERT,KEDGE,KAREA,KBFG,KBFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
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
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the bilinear mapping onto the reference element
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
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
!       DU1=1D0
!       DU2=0D0
!       DU3=0D0
C
      DU1=0D0
      DU2=0D0
      DU3=0D0
C
      CALL EPHI(XI1,XI2,XI3,DPHI)
      DO 220 JDFL=1,NNBAS
        HBAS=DPHI(KBFL(JDFL))
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KBFG(JDFL)
         DU1=DU1+U1(JDFG)*HBAS
         DU2=DU2+U2(JDFG)*HBAS
         DU3=DU3+U3(JDFG)*HBAS
        ENDIF
220   CONTINUE
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(1,JDOFEH,1)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2+HBASJ4*DU3
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=HSUMJ*HBASJ1
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         HSUMI=HBASI2*DU1+HBASI3*DU2+HBASI4*DU3
         AH=HSUMI*HBASJ1
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA)-DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE Conv_LinSc2(U1,U2,U3,DA,NU,KCOLA,KLDA,KVERT,KAREA,
     *           KEDGE,DCORVG,KADJ,KVERT1,KAREA1,KEDGE1,
     *           NEL1,NVT1,NET1,NAT1,ELE)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION U1(*),U2(*),U3(*),DA(*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*),KADJ(NNAE,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KVERT1(NNVE,*),KAREA1(NNAE,*),KEDGE1(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KBFG(NNBAS),KBFL(NNBAS)
      DIMENSION DPHI(NNBAS)
      INTEGER   AdjElem(8)
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
      EXTERNAL EPHI
C
      SAVE
C
      NA=KLDA(NU+1)-1
C
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,4
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      ICUB=7!ICUBN
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL1=1,NEL1
C
      CALL localNDFGL(IEL1,KVERT1,KEDGE1,KAREA1,
     *     NVT1,NET1,NAT1,KBFG,KBFL)
      IF (IER.LT.0) GOTO 99999
C
      AdjElem(1)=IEL1
      AdjElem(2)=KADJ(3,AdjElem(1))
      AdjElem(3)=KADJ(3,AdjElem(2))
      AdjElem(4)=KADJ(3,AdjElem(3))
      AdjElem(5)=KADJ(6,AdjElem(1))
      AdjElem(6)=KADJ(3,AdjElem(5))
      AdjElem(7)=KADJ(3,AdjElem(6))
      AdjElem(8)=KADJ(3,AdjElem(7))

      DO 150 IEL8=1,8

      IEL = AdjElem(IEL8)

      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
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
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the bilinear mapping onto the reference element
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
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DU1=0D0
      DU2=0D0
      DU3=0D0
C
      CALL SETXUI(IEL8,XI1,XI2,XI3,XUI1,XUI2,XUI3)
      CALL EPHI(XUI1,XUI2,XUI3,DPHI)
      DO 220 JDFL=1,NNBAS
        HBAS=DPHI(KBFL(JDFL))
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KBFG(JDFL)
         DU1=DU1+U1(JDFG)*HBAS
         DU2=DU2+U2(JDFG)*HBAS
         DU3=DU3+U3(JDFG)*HBAS
        ENDIF
220   CONTINUE
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(1,JDOFEH,1)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2+HBASJ4*DU3
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=HSUMJ*HBASJ1
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         HSUMI=HBASI2*DU1+HBASI3*DU2+HBASI4*DU3
         AH=HSUMI*HBASJ1
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA)-DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
150   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
      SUBROUTINE localNDFGL(IEL,KVERT,KEDGE,KAREA,NVT,NET,NAT,KG,KL)
      IMPLICIT NONE
      INTEGER I,IEL,NVT,NET,NAT,KG(*),KL(*)
      INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)

      DO I=1,8
       KG(I) = KVERT(I,IEL)
      END DO

      DO I=1,12
       KG(8+I) = NVT + KEDGE(I,IEL)
      END DO

      DO I=1,6
       KG(20+I) = NVT + NET + KAREA(I,IEL)
      END DO

      KG(27) = NVT + NET + NAT + IEL
      
      DO I=1,27
       KL(I) = I
      END DO

      END 
C
C
C
      SUBROUTINE SETXUI(I,X1,X2,X3,Y1,Y2,Y3)
      IMPLICIT NONE
      REAL*8 X1,X2,X3,Y1,Y2,Y3
      INTEGER i

      IF (I.EQ.1) THEN
       Y1 = X1/2d0 - 0.5d0
       Y2 = X2/2d0 - 0.5d0
       Y3 = X3/2d0 - 0.5d0
       RETURN
      END IF

      IF (I.EQ.2) THEN
       Y1 = X1/2d0 + 0.5d0
       Y2 = X2/2d0 - 0.5d0
       Y3 = X3/2d0 - 0.5d0
       RETURN
      END IF

      IF (I.EQ.3) THEN
       Y1 = X1/2d0 + 0.5d0
       Y2 = X2/2d0 + 0.5d0
       Y3 = X3/2d0 - 0.5d0
       RETURN
      END IF

      IF (I.EQ.4) THEN
       Y1 = X1/2d0 - 0.5d0
       Y2 = X2/2d0 + 0.5d0
       Y3 = X3/2d0 - 0.5d0
       RETURN
      END IF

      IF (I.EQ.5) THEN
       Y1 = X1/2d0 - 0.5d0
       Y2 = X2/2d0 - 0.5d0
       Y3 = X3/2d0 + 0.5d0
       RETURN
      END IF

      IF (I.EQ.6) THEN
       Y1 = X1/2d0 + 0.5d0
       Y2 = X2/2d0 - 0.5d0
       Y3 = X3/2d0 + 0.5d0
       RETURN
      END IF

      IF (I.EQ.7) THEN
       Y1 = X1/2d0 + 0.5d0
       Y2 = X2/2d0 + 0.5d0
       Y3 = X3/2d0 + 0.5d0
       RETURN
      END IF

      IF (I.EQ.8) THEN
       Y1 = X1/2d0 - 0.5d0
       Y2 = X2/2d0 + 0.5d0
       Y3 = X3/2d0 + 0.5d0
       RETURN
      END IF

      END
C
C
C
      SUBROUTINE EPHI(X1,X2,X3,DPHI)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (Q8=0.125D0)
      DIMENSION DPHI(27)
C
      DPHI(1)=
     *-Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0-X3)
      DPHI(2)=
     * Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0-X3)
      DPHI(3)=
     *-Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0-X3)
      DPHI(4)= 
     *Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0-X3)
      DPHI(5)= 
     *Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0+X3)
      DPHI(6)= 
     *-Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0+X3)
      DPHI(7)= 
     *Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0+X3)
      DPHI(8)= 
     *-Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0+X3)
      DPHI(9)= 
     *Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*X3*(1D0-X3)
      DPHI(10)= 
     *-Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
      DPHI(11)= 
     *-Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*X3*(1D0-X3)
      DPHI(12)= 
     *Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
      DPHI(13)= 
     *Q8*X1*(1D0-X1)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
      DPHI(14)= 
     *-Q8*X1*(1D0+X1)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
      DPHI(15)= 
     *Q8*X1*(1D0+X1)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
      DPHI(16)= 
     *-Q8*X1*(1D0-X1)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
      DPHI(17)= 
     *-Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*X3*(1D0+X3)
      DPHI(18)= 
     *Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
      DPHI(19)= 
     *Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*X3*(1D0+X3)
      DPHI(20)= 
     *-Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
      DPHI(21)=
     * -Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
      DPHI(22)=
     * -Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
      DPHI(23)=
     * Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
      DPHI(24)=
     * Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
      DPHI(25)=
     * -Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
      DPHI(26)=
     * Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
      DPHI(27)=
     * Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
C
99999 END
C
C
C
