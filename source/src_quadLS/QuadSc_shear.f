************************************************************************
      SUBROUTINE GETNRMSHR(DSTRES,DU1,DU2,DU3,KVERT,KAREA,
     *           KEDGE,DCORVG,ELE)
************************************************************************

*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DSTRES(*),DU1(*),DU2(*),DU3(*)
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
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
      DO 2 I=1,4
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      ICUB=1
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
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices
      DX0=0d0
      DY0=0d0
      DZ0=0d0
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
      DX0=DX0+0.125d0*DX(IVE)
      DY0=DY0+0.125d0*DY(IVE)
      DZ0=DZ0+0.125d0*DZ(IVE)
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
      ICUBP = 1
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
      UH1X=0D0
      UH2X=0D0
      UH3X=0D0
      UH1Y=0D0
      UH2Y=0D0
      UH3Y=0D0
      UH1Z=0D0
      UH2Z=0D0
      UH3Z=0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
C
      UH1X=UH1X+DU1(IEQ)*DBAS(1,ILO,2)
      UH2X=UH2X+DU2(IEQ)*DBAS(1,ILO,2)
      UH3X=UH3X+DU3(IEQ)*DBAS(1,ILO,2)
      UH1Y=UH1Y+DU1(IEQ)*DBAS(1,ILO,3)
      UH2Y=UH2Y+DU2(IEQ)*DBAS(1,ILO,3)
      UH3Y=UH3Y+DU3(IEQ)*DBAS(1,ILO,3)
      UH1Z=UH1Z+DU1(IEQ)*DBAS(1,ILO,4)
      UH2Z=UH2Z+DU2(IEQ)*DBAS(1,ILO,4)
      UH3Z=UH3Z+DU3(IEQ)*DBAS(1,ILO,4)
210   CONTINUE
C
!      IF (DX0.LT.2.25d0) THEN
      DAUX = 0.0002d0**2d0
      DAUX = DAUX + UH1X**2d0       +  UH2Y**2d0       +  UH3Z**2d0 +
     *      (UH1Y+UH2X)**2d0 + (UH1Z+UH3X)**2d0 + (UH2Z+UH3Y)**2d0
      DSTRES(IEL) = DSQRT(DAUX)
!      IF (ISNAN(DSTRES(IEL))) THEN
!       WRITE(*,*) iel,daux
!       WRITE(*,*) UH1X,UH2Y,UH3Z,UH1Y,UH2X,UH1Z,UH3X,UH2Z,UH3Y
!       pause 
!      END IF
 !     ELSE
  !    DSTRES(IEL) = 0d0
   !   END IF
C
C
100   CONTINUE
C
99999 END
C
C
C
      SUBROUTINE GETMGNRMSHR(SHEAR2,SHEAR1,KADJ2,VOL2,NEL1)
      INTEGER KADJ2(6,*)
      REAL*4 VOL2(*)
      REAL*8 SHEAR2(*),SHEAR1(*),DAUX1,DAUX2
      INTEGER IEL,NEL1,JEL(8)

      DO IEL=1,NEL1
       JEL(1)  = IEL
       JEL(2)  = KADJ2(3,JEL(1))
       JEL(3)  = KADJ2(3,JEL(2))
       JEL(4)  = KADJ2(3,JEL(3))
       JEL(5)  = KADJ2(6,JEL(1))
       JEL(6)  = KADJ2(6,JEL(2))
       JEL(7)  = KADJ2(6,JEL(3))
       JEL(8)  = KADJ2(6,JEL(4))

       DAUX1 = 0d0
       DAUX2 = 0d0

       DO I=1,8
        DAUX1 = DAUX1 + DBLE(VOL2(JEL(I)))
        DAUX2 = DAUX2 + DBLE(VOL2(JEL(I)))*SHEAR2(JEL(I))
       END DO

       SHEAR1(IEL) = DAUX2/DAUX1

      END DO

      END

      SUBROUTINE GETSHEARVISC(DSTRS,DVISC,NEL,D1,D2)
      REAL*8 DSTRS(*),DVISC(*),D1,D2
      INTEGER IEL,NEL

      DO IEL=1,NEL
!      IF (DSTRS(IEL).GT.1d-8) THEN
!       DVISC(IEL) = D1*DSTRS(IEL)**(D2-1d0)
       DVISC(IEL) = MIN(D1*DSTRS(IEL)**(D2-1d0),0.5d0)
!       DVISC(IEL) = MAX(D1*DSTRS(IEL)**(D2-1d0),0.00001d0)
 !     ELSE
  !     DVISC(IEL) = D1
   !   END IF
      END DO

      END
