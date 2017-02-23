************************************************************************
      SUBROUTINE GetForceCyl(U1,U2,U3,P,bALPHA,DVISC,KVERT,KAREA,KEDGE,
     *                     DCORVG,DResForce,ELE)
************************************************************************
*     Discrete convection operator: Q1~ elements (nonparametric)
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      REAL*8  U1(*),U2(*),U3(*),P(*),DVISC(*),DCORVG(NNDIM,*)
      REAL*8  DResForce(3)
      LOGICAL bALPHA(*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
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
      SAVE
C
      IF (myid.eq.0) GOTO 999
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
      ICUB=9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      DResForce(1) = 0D0
      DResForce(2) = 0D0
      DResForce(3) = 0D0
C
C *** Loop over all elements
      nnel = 0
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C-----------------------------------------------------------------
       NJALFA=0
       NIALFA=0
       DO I=1,IDFL
         IG=KDFG(I)
         IF (bALPHA(IG)) THEN
          NJALFA=NJALFA+1
         ENDIF
         IF (.NOT.bALPHA(IG)) THEN
          NIALFA=NIALFA+1
         ENDIF
       ENDDO
C--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
      IF(NJALFA.EQ.27.OR.NIALFA.EQ.27) GOTO 100
C--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
C
      nnel = nnel + 1
      DNY = DVISC(IEL)
C
C *** Evaluation of coordinates of the vertices
      DX0 = 0d0
      DY0 = 0d0
      DZ0 = 0d0
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
      DX0 = DX0 + 0.125d0*DX(IVE)
      DY0 = DY0 + 0.125d0*DY(IVE)
      DZ0 = DZ0 + 0.125d0*DZ(IVE)
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
      CALL ELE(0D0,0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
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
      XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
      YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
      ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
C     Evaluate the solution values and derivatives in the cubature point

       DU1V=0D0     ! U1 value
       DU1X=0D0     ! U1 x deriv
       DU1Y=0D0     ! U1 y deriv
       DU1Z=0D0     ! U1 z deriv
C
       DU2V=0D0     ! U2 value
       DU2X=0D0     ! U2 x deriv
       DU2Y=0D0     ! U2 y deriv
       DU2Z=0D0     ! U2 z deriv
C
       DU3V=0D0     ! U3 value
       DU3X=0D0     ! U3 x deriv
       DU3Y=0D0     ! U3 y deriv
       DU3Z=0D0     ! U3 z deriv
C
       DALV=0D0     ! ALFA value
       DALX=0D0     ! ALFA x deriv
       DALY=0D0     ! ALFA y deriv
       DALZ=0D0     ! ALFA z deriv

       DO I=1,IDFL
         IG=KDFG(I)
         DBI1=DBAS(1,KDFL(I),1)
         DBI2=DBAS(1,KDFL(I),2)
         DBI3=DBAS(1,KDFL(I),3)
         DBI4=DBAS(1,KDFL(I),4)
C---------------FOR U1----------------
         DU1V=DU1V+U1(IG)*DBI1
         DU1X=DU1X+U1(IG)*DBI2
         DU1Y=DU1Y+U1(IG)*DBI3
         DU1Z=DU1Z+U1(IG)*DBI4
C---------------FOR U2----------------
         DU2V=DU2V+U2(IG)*DBI1
         DU2X=DU2X+U2(IG)*DBI2
         DU2Y=DU2Y+U2(IG)*DBI3
         DU2Z=DU2Z+U2(IG)*DBI4
C---------------FOR U3----------------
         DU3V=DU3V+U3(IG)*DBI1
         DU3X=DU3X+U3(IG)*DBI2
         DU3Y=DU3Y+U3(IG)*DBI3
         DU3Z=DU3Z+U3(IG)*DBI4
C---------------FOR ALFA----------------
         IF (bALPHA(IG)) THEN
          DALPHA = 1d0
         ELSE
          DALPHA = 0d0
         END IF
         DALV=DALV+DALPHA*DBI1
         DALX=DALX+DALPHA*DBI2
         DALY=DALY+DALPHA*DBI3
         DALZ=DALZ+DALPHA*DBI4
       ENDDO
C
       JJ = 4*(IEL-1) + 1
       Press =          P(JJ  ) + (XX-DX0)*P(JJ+1) +
     *         (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)
C--------------------------------------------------------
c-----------Form the integrand------------------
       DN1=-DALX
       DN2=-DALY
       DN3=-DALZ
c-------------------Acting force-------------------------
c--------------Deformation calculation-------------
C
       AH1=-Press*DN1 + DNY*(DU1X*DN1 + DU1Y*DN2 + DU1Z*DN3)
       AH2=-Press*DN2 + DNY*(DU2X*DN1 + DU2Y*DN2 + DU2Z*DN3)
       AH3=-Press*DN3 + DNY*(DU3X*DN1 + DU3Y*DN2 + DU3Z*DN3)
C
       DResForce(1) = DResForce(1) + AH1*OM
       DResForce(2) = DResForce(2) + AH2*OM
       DResForce(3) = DResForce(3) + AH3*OM
C
200   CONTINUE
C
100   CONTINUE
C
999   CALL COMM_SUMMN(DResForce,3)
C
99999 CONTINUE

      END
C
C
C
************************************************************************
      SUBROUTINE GetForces(U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,
     *                     DCORVG,ELE)
************************************************************************
*     Discrete convection operator: Q1~ elements (nonparametric)
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      USE var_QuadScalar, ONLY : myFBM
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      REAL*8  U1(*),U2(*),U3(*),P(*),DVISC(*),DCORVG(NNDIM,*)
      INTEGER ALPHA(*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
C
      REAL*8 DResForceX,DResForceY,DResForceZ
      REAL*8 DTrqForceX,DTrqForceY,DTrqForceZ
      REAL*8 Center(3),dForce(6)
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
      SAVE
C
      IF (myid.eq.0) GOTO 999
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
      ICUB=9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      DO IP = 1,myFBM%nParticles

      Center = myFBM%particleNew(IP)%Position
C
      DResForceX = 0D0
      DResForceY = 0D0
      DResForceZ = 0D0
      DTrqForceX = 0d0
      DTrqForceY = 0d0
      DTrqForceZ = 0d0
C
C *** Loop over all elements
      nnel = 0
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C-----------------------------------------------------------------
       NJALFA=0
       NIALFA=0
       DO I=1,IDFL
         IG=KDFG(I)
         IF (ALPHA(IG).EQ.0) THEN
          NJALFA=NJALFA+1
         ENDIF
         IF (ALPHA(IG).EQ.IP) THEN
          NIALFA=NIALFA+1
         ENDIF
       ENDDO
C--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
      IF(NJALFA.EQ.27.OR.NIALFA.EQ.27) GOTO 100
C--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
C
      nnel = nnel + 1
      DNY = DVISC(IEL)
C
C *** Evaluation of coordinates of the vertices
      DX0 = 0d0
      DY0 = 0d0
      DZ0 = 0d0
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
      DX0 = DX0 + 0.125d0*DX(IVE)
      DY0 = DY0 + 0.125d0*DY(IVE)
      DZ0 = DZ0 + 0.125d0*DZ(IVE)
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
      CALL ELE(0D0,0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
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
      XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
      YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
      ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
C     Evaluate the solution values and derivatives in the cubature point

       DU1V=0D0     ! U1 value
       DU1X=0D0     ! U1 x deriv
       DU1Y=0D0     ! U1 y deriv
       DU1Z=0D0     ! U1 z deriv
C
       DU2V=0D0     ! U2 value
       DU2X=0D0     ! U2 x deriv
       DU2Y=0D0     ! U2 y deriv
       DU2Z=0D0     ! U2 z deriv
C
       DU3V=0D0     ! U3 value
       DU3X=0D0     ! U3 x deriv
       DU3Y=0D0     ! U3 y deriv
       DU3Z=0D0     ! U3 z deriv
C
       DALV=0D0     ! ALFA value
       DALX=0D0     ! ALFA x deriv
       DALY=0D0     ! ALFA y deriv
       DALZ=0D0     ! ALFA z deriv

       DO I=1,IDFL
         IG=KDFG(I)
         DBI1=DBAS(1,KDFL(I),1)
         DBI2=DBAS(1,KDFL(I),2)
         DBI3=DBAS(1,KDFL(I),3)
         DBI4=DBAS(1,KDFL(I),4)
C---------------FOR U1----------------
         DU1V=DU1V+U1(IG)*DBI1
         DU1X=DU1X+U1(IG)*DBI2
         DU1Y=DU1Y+U1(IG)*DBI3
         DU1Z=DU1Z+U1(IG)*DBI4
C---------------FOR U2----------------
         DU2V=DU2V+U2(IG)*DBI1
         DU2X=DU2X+U2(IG)*DBI2
         DU2Y=DU2Y+U2(IG)*DBI3
         DU2Z=DU2Z+U2(IG)*DBI4
C---------------FOR U3----------------
         DU3V=DU3V+U3(IG)*DBI1
         DU3X=DU3X+U3(IG)*DBI2
         DU3Y=DU3Y+U3(IG)*DBI3
         DU3Z=DU3Z+U3(IG)*DBI4
C---------------FOR ALFA----------------
         IF (ALPHA(IG).EQ.IP) THEN
          DALPHA = 1d0
         ELSE
          DALPHA = 0d0
         END IF
         DALV=DALV+DALPHA*DBI1
         DALX=DALX+DALPHA*DBI2
         DALY=DALY+DALPHA*DBI3
         DALZ=DALZ+DALPHA*DBI4
       ENDDO
C
       JJ = 4*(IEL-1) + 1
       Press =          P(JJ  ) + (XX-DX0)*P(JJ+1) +
     *         (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)
C--------------------------------------------------------
c-----------Form the integrand------------------
       DN1=-DALX
       DN2=-DALY
       DN3=-DALZ
c-------------------Acting force-------------------------
c--------------Deformation calculation-------------
C
       AH1=-Press*DN1 + DNY*(DU1X*DN1 + DU1Y*DN2 + DU1Z*DN3)
       AH2=-Press*DN2 + DNY*(DU2X*DN1 + DU2Y*DN2 + DU2Z*DN3)
       AH3=-Press*DN3 + DNY*(DU3X*DN1 + DU3Y*DN2 + DU3Z*DN3)
C
!        AH1=-P(IEL)*DN1+DNY*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2 + ! full3D
!      *     (DU1Z+DU3X)*DN3)
!        AH2=-P(IEL)*DN2+DNY*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2 + ! full3D
!      *     (DU2Z+DU3Y)*DN3)
!        AH3=-P(IEL)*DN3+DNY*((DU3X+DU1Z)*DN1+(DU3Y+DU2Z)*DN2 + ! full3D
!      *     (DU3Z+DU3Z)*DN3)
C
c--------------------------------------------------------------
       DResForceX = DResForceX + AH1*OM
       DResForceY = DResForceY + AH2*OM
       DResForceZ = DResForceZ + AH3*OM
c-------------------Torque force------------------------- 
       XTORQUE = XX - Center(1)
       YTORQUE = YY - Center(2)
       ZTORQUE = ZZ - Center(3)
C
       ATQX = YTORQUE*AH3 - ZTORQUE*AH2                        ! full3D
       ATQY = ZTORQUE*AH1 - XTORQUE*AH3                        ! full3D
       ATQZ = XTORQUE*AH2 - YTORQUE*AH1
C
       DTrqForceX = DTrqForceX + ATQX*OM
       DTrqForceY = DTrqForceY + ATQY*OM
       DTrqForceZ = DTrqForceZ + ATQZ*OM
C
C===============================================================
C--------------UPDATE VELOCITY AND POSITION---------------------
C
200   CONTINUE
C
100   CONTINUE
C
      iPointer = 6*(IP-1)
      myFBM%Force(iPointer+1) = DResForceX
      myFBM%Force(iPointer+2) = DResForceY
      myFBM%Force(iPointer+3) = DResForceZ
      myFBM%Force(iPointer+4) = DTrqForceX
      myFBM%Force(iPointer+5) = DTrqForceY
      myFBM%Force(iPointer+6) = DTrqForceZ

      END DO ! nParticles

999   CALL COMM_SUMMN(myFBM%Force,6*myFBM%nParticles)

      DO IP = 1,myFBM%nParticles
       iPointer = 6*(IP-1)+1
       myFBM%ParticleNew(IP)%ResistanceForce = myFBM%Force(iPointer:)
       myFBM%ParticleNew(IP)%TorqueForce = myFBM%Force(iPointer+3:)
      END DO
!       WRITE(*,'(I8,6D12.4)') nnel,dForce

99999 CONTINUE

      END
C
C
C
      SUBROUTINE updateFBM(DensityL,dTime,Gravity,myid)
      USE var_QuadScalar, ONLY : myFBM
      INTEGER myid
      INTEGER IP
      REAL*8 DensityL,dTime,Gravity(3)
      REAL*8 pi,volume,mass,massR,radius,dimomir
      PARAMETER (PI=3.1415926535897931D0)
      REAL*8 RForce(3),dVelocity(3),dOmega(3)

      DO IP = 1,myFBM%nParticles

      radius    = myFBM%particleNew(IP)%sizes(1)
      volume    = (4d0*PI*radius**3d0)/3d0
      mass      = volume*myFBM%particleNew(IP)%density
      massR     = volume*(myFBM%particleNew(IP)%density-DensityL)
      dimomir   = mass*(radius**2d0)*2d0/5d0

      ! Mean resistance forces for the given time step
      RForce    = 0.5D0*(myFBM%particleNew(IP)%ResistanceForce +
     *                   myFBM%particleOld(IP)%ResistanceForce)

      ! Velocity difference for the given time step
      dVelocity = dTime*(RForce+massR*gravity)/mass
!       IF (myid.eq.1) write(*,'(6D12.4)')RForce,dVelocity

       ! Angular velocity difference for the given time step
      dOmega    = dTime*0.5D0*(myFBM%particleNew(IP)%TorqueForce +
     *            myFBM%particleOld(IP)%TorqueForce)/dimomir

!       IF (myid.eq.1) THEN
! !        WRITE(*,*) RForce,dVelocity,dOmega
!        WRITE(*,*) dTime,RForce,massR,gravity,mass
!       END IF
      ! Update the positions
      myFBM%particleNew(IP)%Position = myFBM%particleOld(IP)%Position +
     *dTime*(myFBM%particleOld(IP)%Velocity + 0.5d0*dVelocity)

      ! Update the velocities
      myFBM%particleNew(IP)%Velocity = 
     *myFBM%particleOld(IP)%Velocity + dVelocity

      ! Update the angles
      myFBM%particleNew(IP)%Angle = myFBM%particleOld(IP)%Angle +
     *dTime*(myFBM%particleOld(IP)%AngularVelocity + 0.5d0*dOmega)

      ! Update the angular velocity
      myFBM%particleNew(IP)%AngularVelocity = 
     *myFBM%particleOld(IP)%AngularVelocity + dOmega

      ! Collision with walls
      !----------------------------------------------------------
      IF (myFBM%particleNew(IP)%Position(1).GT.0.4d0-radius) THEN
       myFBM%particleNew(IP)%Position(1) = 0.4d0-radius
       myFBM%particleNew(IP)%Velocity(1) = 
     *-myFBM%particleNew(IP)%Velocity(1)
      END IF
      IF (myFBM%particleNew(IP)%Position(1).LT.radius) THEN
       myFBM%particleNew(IP)%Position(1) = radius
       myFBM%particleNew(IP)%Velocity(1) = 
     *-myFBM%particleNew(IP)%Velocity(1)
      END IF

      IF (myFBM%particleNew(IP)%Position(2).GT.0.4d0-radius) THEN
       myFBM%particleNew(IP)%Position(2) = 0.4d0-radius
       myFBM%particleNew(IP)%Velocity(2) = 
     *-myFBM%particleNew(IP)%Velocity(2)
      END IF
      IF (myFBM%particleNew(IP)%Position(2).LT.radius) THEN
       myFBM%particleNew(IP)%Position(2) = radius
       myFBM%particleNew(IP)%Velocity(2) = 
     *-myFBM%particleNew(IP)%Velocity(2)
      END IF

      IF (myFBM%particleNew(IP)%Position(3).GT.0.4d0-radius) THEN
       myFBM%particleNew(IP)%Position(3) = 0.4d0-radius
       myFBM%particleNew(IP)%Velocity(3) = 
     *-myFBM%particleNew(IP)%Velocity(3)
      END IF
      IF (myFBM%particleNew(IP)%Position(3).LT.radius) THEN
       myFBM%particleNew(IP)%Position(3) = radius
       myFBM%particleNew(IP)%Velocity(3) = 
     *-myFBM%particleNew(IP)%Velocity(3)
      END IF
      !----------------------------------------------------------

      ! Update the forces
      myFBM%particleOld(IP)%ResistanceForce =
     *myFBM%particleNew(IP)%ResistanceForce

      ! Update the torque
      myFBM%particleOld(IP)%TorqueForce =
     *myFBM%particleNew(IP)%TorqueForce

      ! Update the positions
      myFBM%particleOld(IP)%Position =
     *myFBM%particleNew(IP)%Position

      ! Update the velocities
      myFBM%particleOld(IP)%Velocity = 
     *myFBM%particleNew(IP)%Velocity

      ! Update the angles
      myFBM%particleOld(IP)%Angle = 
     *myFBM%particleNew(IP)%Angle

      ! Update the angular velocity
      myFBM%particleOld(IP)%AngularVelocity = 
     *myFBM%particleNew(IP)%AngularVelocity


      IF (myid.eq.1) THEN
      WRITE(*,'(A19,3D12.4)') "ResistanceForce: ",
     *     myFBM%particleNew(IP)%ResistanceForce
      WRITE(*,'(A19,3D12.4)') "TorqueForce: ",
     *     myFBM%particleNew(IP)%TorqueForce
      WRITE(*,'(A19,3D12.4)') "Position: ",
     *     myFBM%particleNew(IP)%Position
      WRITE(*,'(A19,3D12.4)') "Velocity: ",
     *     myFBM%particleNew(IP)%Velocity
      WRITE(*,'(A19,3D12.4)') "Angle: ",
     *     myFBM%particleNew(IP)%Angle
      WRITE(*,'(A19,3D12.4)') "AngularVelocity: ",
     *     myFBM%particleNew(IP)%AngularVelocity
      END IF

      END DO

      END
C
C
C
      SUBROUTINE GetFictKnpr(X,Y,Z,iBndr,inpr,Dist)
      USE var_QuadScalar, ONLY : myFBM
      IMPLICIT NONE
      REAL*8 X,Y,Z,Dist,daux
      REAL*8 PX,PY,PZ,RAD
      INTEGER iBndr,inpr,iP,iaux,ipc,isin

      inpr = 0
      Dist = 0.05d0

      DO iP=1,myFBM%nParticles
       PX = myFBM%particleNew(iP)%Position(1)
       PY = myFBM%particleNew(iP)%Position(2)
       PZ = myFBM%particleNew(iP)%Position(3)
       RAD =myFBM%particleNew(iP)%sizes(1)
       ipc=ip-1
       isin = 0
       if( ((x-px)**2 + (y-py)**2 + (z-pz)**2) .le. rad**2)then
          Dist = 1.0d0
          inpr = IP
       end if
      END DO

      END SUBROUTINE GetFictKnpr
C
C
C
      SUBROUTINE GetVeloFictBCVal(X,Y,Z,ValU,ValV,ValW,IP,t)
      USE var_QuadScalar, ONLY : myFBM
      IMPLICIT NONE
      INTEGER iP
      REAL*8 X,Y,Z,t,ValU,ValV,ValW
      REAL*8 PX,PY,PZ,RAD
      REAL*8 Velo(3),Pos(3),Omega(3)
      REAL*8 DVELZ_X,DVELZ_Y,DVELY_Z,DVELY_X,DVELX_Y,DVELX_Z

      Velo  = myFBM%particleNEW(IP)%Velocity
      Pos   = myFBM%particleNEW(IP)%Position
      Omega = myFBM%particleNEW(IP)%AngularVelocity

      DVELZ_X = -(Y-Pos(2))*Omega(3)
      DVELZ_Y = +(X-Pos(1))*Omega(3)

      DVELY_Z = -(X-Pos(1))*Omega(2)
      DVELY_X = +(Z-Pos(3))*Omega(2)

      DVELX_Y = -(Z-Pos(3))*Omega(1)
      DVELX_Z = +(Y-Pos(2))*Omega(1)

      ValU = Velo(1) + DVELZ_X + DVELY_X
      ValV = Velo(2) + DVELZ_Y + DVELX_Y
      ValW = Velo(3) + DVELX_Z + DVELY_Z

      END SUBROUTINE GetVeloFictBCVal

      subroutine communicateforce(fx,fy,fz,tx,ty,tz)
      use var_QuadScalar, only : myFBM
      implicit none
      real*8 fx(*),fy(*),fz(*),tx(*),ty(*),tz(*)
      integer IP

      do IP = 1,myFBM%nParticles
        fx(IP)=myFBM%particleNew(IP)%ResistanceForce(1)
        fy(IP)=myFBM%particleNew(IP)%ResistanceForce(2)
        fz(IP)=myFBM%particleNew(IP)%ResistanceForce(3)
        tx(IP)=myFBM%particleNew(IP)%TorqueForce(1)
        ty(IP)=myFBM%particleNew(IP)%TorqueForce(2)
        tz(IP)=myFBM%particleNew(IP)%TorqueForce(3)
      end do
      end subroutine







