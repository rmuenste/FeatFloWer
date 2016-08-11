 
************************************************************************
      SUBROUTINE Build_BTMatP1(DAx,DAy,DAz,KLD,KCOL,KVERT,KAREA,KEDGE,
     *           DCORVG,NA,ICUB,ELE)
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
      DIMENSION DAx(*),DAy(*),DAz(*),DBAS1(4)
      DIMENSION DCORVG(NNDIM,*),KLD(*),KCOL(*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG1(NNBAS),KDFL1(NNBAS)
      DIMENSION KDFG2(NNBAS),KDFL2(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRYx(NNBAS,NNBAS)
      DIMENSION DENTRYy(NNBAS,NNBAS),DENTRYz(NNBAS,NNBAS)

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
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I= 1,4
2     BDER(I)=.TRUE.
C
      IELTYP1=12
      IDFL1=NDFL(IELTYP1)
C
      IELTYP2=-1
      CALL ELE(0D0,0D0,0D0,IELTYP2)
      IDFL2=NDFL(IELTYP2)
C
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP1,KVERT,KEDGE,KAREA,KDFG1,KDFL1)
      IF (IER.LT.0) GOTO 99999
C
      CALL NDFGL(IEL,1,IELTYP2,KVERT,KEDGE,KAREA,KDFG2,KDFL2)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
!       WRITE(*,*)IELTYP1,IELTYP2,KDFG1,KDFL1,KDFG2,KDFL2
!       WRITE(*,*)IDFL1,IDFL2
!       pause

      DO 110 JDOFE=1,IDFL1
      JCOL0=KLD(KDFG1(JDOFE))
      DO 110 IDOFE=1,IDFL2
      IDFG=KDFG2(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOL(JCOL).EQ.IDFG) GOTO 111
112   CONTINUE
111   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRYx(JDOFE,IDOFE)=0d0
      DENTRYy(JDOFE,IDOFE)=0d0
      DENTRYz(JDOFE,IDOFE)=0d0
110   CONTINUE
C
!       PAUSE
C *** Evaluation of coordinates of the vertices
      DX0I = 0d0
      DY0I = 0d0
      DZ0I = 0d0
      DO 120 IVE=1,NVE
      IP=KVERT(IVE,IEL)
      KVE(IVE)=IP
      DX(IVE)=DCORVG(1,IP)
      DY(IVE)=DCORVG(2,IP)
      DZ(IVE)=DCORVG(3,IP)
      DX0I = DX0I + 0.125d0*DCORVG(1,IP)
      DY0I = DY0I + 0.125d0*DCORVG(2,IP)
      DZ0I = DZ0I + 0.125d0*DCORVG(3,IP)
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
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
!       WRITE(*,*) detj
!       pause
C
! ----------------------------------------------------------------
!     Computation of the cartesian coordiante of the cubature point
       XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
       YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
       ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
! ----------------------------------------------------------------
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DBAS1(1) = 1d0
      DBAS1(2) = XX-DX0I
      DBAS1(3) = YY-DY0I
      DBAS1(4) = ZZ-DZ0I
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL1
       JDOFEH=KDFL1(JDOFE)
       HBASJ1=DBAS1(JDOFEH)
       HSUMJ=HBASJ1
C
!        write(*,'(I4,4(1xG12.6))') JDOFE,hsumj,hbasj2,du1
       DO 240 IDOFE=1,IDFL2
        IDOFEH=KDFL2(IDOFE)
        HBASI2=DBAS(1,IDOFEH,2)
        HBASI3=DBAS(1,IDOFEH,3)
        HBASI4=DBAS(1,IDOFEH,4)

        AHx=OM*HSUMJ*HBASI2
        AHy=OM*HSUMJ*HBASI3
        AHz=OM*HSUMJ*HBASI4

        DENTRYx(JDOFE,IDOFE)=DENTRYx(JDOFE,IDOFE)+AHx
        DENTRYy(JDOFE,IDOFE)=DENTRYy(JDOFE,IDOFE)+AHy
        DENTRYz(JDOFE,IDOFE)=DENTRYz(JDOFE,IDOFE)+AHz

240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL1
      DO 300 IDOFE=1,IDFL2
        IA    =KENTRY(JDOFE,IDOFE)
        DAx(IA)=DAx(IA)+DENTRYx(JDOFE,IDOFE)
        DAy(IA)=DAy(IA)+DENTRYy(JDOFE,IDOFE)
        DAz(IA)=DAz(IA)+DENTRYz(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
!       IF (iEQ.EQ.1) WRITE(*,*) DMAXDIVU

99999 END
C
C
C
************************************************************************
      SUBROUTINE Build_BMatP1(DAx,DAy,DAz,KLD,KCOL,KVERT,KAREA,KEDGE,
     *           DCORVG,NA,ICUB,ELE)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY : myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DAx(*),DAy(*),DAz(*),DBAS1(4)
      DIMENSION DCORVG(NNDIM,*),KLD(*),KCOL(*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG1(NNBAS),KDFL1(NNBAS)
      DIMENSION KDFG2(NNBAS),KDFL2(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRYx(NNBAS,NNBAS)
      DIMENSION DENTRYy(NNBAS,NNBAS),DENTRYz(NNBAS,NNBAS)

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
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I= 1,4
2     BDER(I)=.TRUE.
C
      IELTYP1=-1
      CALL ELE(0D0,0D0,0D0,IELTYP1)
      IDFL1=NDFL(IELTYP1)
C
      IELTYP2=12
      IDFL2=NDFL(IELTYP2)
C
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP1,KVERT,KEDGE,KAREA,KDFG1,KDFL1)
      IF (IER.LT.0) GOTO 99999
C
      CALL NDFGL(IEL,1,IELTYP2,KVERT,KEDGE,KAREA,KDFG2,KDFL2)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
!       WRITE(*,*)IELTYP1,IELTYP2,KDFG1,KDFL1,KDFG2,KDFL2
!       WRITE(*,*)IDFL1,IDFL2
!       pause

      DO 110 JDOFE=1,IDFL1
      JCOL0=KLD(KDFG1(JDOFE))
      DO 110 IDOFE=1,IDFL2
      IDFG=KDFG2(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOL(JCOL).EQ.IDFG) GOTO 111
112   CONTINUE
111   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRYx(JDOFE,IDOFE)=0d0
      DENTRYy(JDOFE,IDOFE)=0d0
      DENTRYz(JDOFE,IDOFE)=0d0
110   CONTINUE

C
!       PAUSE
C *** Evaluation of coordinates of the vertices
      DX0I = 0d0
      DY0I = 0d0
      DZ0I = 0d0
      DO 120 IVE=1,NVE
      IP=KVERT(IVE,IEL)
      KVE(IVE)=IP
      DX(IVE)=DCORVG(1,IP)
      DY(IVE)=DCORVG(2,IP)
      DZ(IVE)=DCORVG(3,IP)
      DX0I = DX0I + 0.125d0*DCORVG(1,IP)
      DY0I = DY0I + 0.125d0*DCORVG(2,IP)
      DZ0I = DZ0I + 0.125d0*DCORVG(3,IP)
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
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
!       WRITE(*,*) detj
!       pause
C
! ----------------------------------------------------------------
!     Computation of the cartesian coordiante of the cubature point
       XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
       YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
       ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
! ----------------------------------------------------------------
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DBAS1(1) = 1d0
      DBAS1(2) = XX-DX0I
      DBAS1(3) = YY-DY0I
      DBAS1(4) = ZZ-DZ0I
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL1
       JDOFEH=KDFL1(JDOFE)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
C
!        write(*,'(I4,4(1xG12.6))') JDOFE,hsumj,hbasj2,du1
       DO 240 IDOFE=1,IDFL2
        IDOFEH=KDFL2(IDOFE)
        HBASI1=DBAS1(IDOFEH)

        AHx=OM*HBASJ2*HBASI1
        AHy=OM*HBASJ3*HBASI1
        AHz=OM*HBASJ4*HBASI1

        DENTRYx(JDOFE,IDOFE)=DENTRYx(JDOFE,IDOFE)+AHx
        DENTRYy(JDOFE,IDOFE)=DENTRYy(JDOFE,IDOFE)+AHy
        DENTRYz(JDOFE,IDOFE)=DENTRYz(JDOFE,IDOFE)+AHz

240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
           
      DO 300 JDOFE=1,IDFL1
      DO 300 IDOFE=1,IDFL2
        IA    =KENTRY(JDOFE,IDOFE)
        DAx(IA)=DAx(IA)+DENTRYx(JDOFE,IDOFE)
        DAy(IA)=DAy(IA)+DENTRYy(JDOFE,IDOFE)
        DAz(IA)=DAz(IA)+DENTRYz(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
!       IF (iEQ.EQ.1) WRITE(*,*) DMAXDIVU

99999 END

      SUBROUTINE ProlPressure(D1,D2,KADJ2,KVERT2,DCORVG2,NEL)
      IMPLICIT NONE
      INTEGER KVERT2(8,*),KADJ2(6,*),NEL
      REAL*8 DCORVG2(3,*),D1(*),D2(*)
      REAL*8 DAUX(5),DX,DY,DZ,DV,DX0,DY0,DZ0,DXP,DYP,DZP
      INTEGER IEL,JEL(8),I,J,K,IND_1,IND_2,JJEL

      DO IEL=1,NEL

       DX0 = DCORVG2(1,KVERT2(7,IEL))
       DY0 = DCORVG2(2,KVERT2(7,IEL))
       DZ0 = DCORVG2(3,KVERT2(7,IEL))

       JEL(1)  = IEL
       JEL(2)  = KADJ2(3,JEL(1))
       JEL(3)  = KADJ2(3,JEL(2))
       JEL(4)  = KADJ2(3,JEL(3))
       JEL(5)  = KADJ2(6,JEL(1))
       JEL(6)  = KADJ2(6,JEL(2))
       JEL(7)  = KADJ2(6,JEL(3))
       JEL(8)  = KADJ2(6,JEL(4))

       IND_1 = 4*(IEL-1) + 1

       DO J=1,8
        DXP = 0d0
        DYP = 0d0
        DZP = 0d0
        DO K=1,8
         DXP = DXP + 0.125D0*DCORVG2(1,KVERT2(K,JEL(J)))
         DYP = DYP + 0.125D0*DCORVG2(2,KVERT2(K,JEL(J)))
         DZP = DZP + 0.125D0*DCORVG2(3,KVERT2(K,JEL(J)))
        END DO

        IND_2 = 4*(JEL(J)-1)+1

        D2(IND_2  ) = D1(IND_1  )           + D1(IND_1+1)*(DXP-DX0)
     *              + D1(IND_1+2)*(DYP-DY0) + D1(IND_1+3)*(DZP-DZ0)
        D2(IND_2+1) = D1(IND_1+1)
        D2(IND_2+2) = D1(IND_1+2)
        D2(IND_2+3) = D1(IND_1+3)

       END DO

      END DO

      END
