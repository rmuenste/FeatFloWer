************************************************************************
      SUBROUTINE Grav_QuadSc(DG1,DG2,DG3,DD,DG,NU,KVERT,KAREA,KEDGE,
     *           KINT,DCORVG,TSTEP,ELE)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid
      USE var_QuadScalar, ONLY : transform
C     
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DG1(*),DG2(*),DG3(*),DD(*),DG(3)
      DIMENSION KINT(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C     --------------------------- Transformation -------------------------------
      REAL*8    DHELP_Q2(27,4,NNCUBP),DHELP_Q1(8,4,NNCUBP)
      REAL*8    DPP(3)
C     --------------------------------------------------------------------------
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
      ICUB=9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
       CALL E013A(XI1,XI2,XI3,DHELP_Q2,ICUBP)
      END DO
C
      BDER(2:4)=.FALSE.
C     
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
      ILINT = 1
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
       CALL E013A(XI1,XI2,XI3,DHELP_Q2,ICUBP)
      END DO
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
!      ILINT=KINT(IEL)
      DENSITY = DD(IEL)
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
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
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      IF (Transform%ILINT.eq.2) THEN ! Q2
      DO JDOFE=1,27
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q2(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q2(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q2(JDFL,4,ICUBP)
      END DO
      END IF
      IF (Transform%ILINT.eq.1) THEN ! Q1
      DO JDOFE=1,8
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q1(JDFL,4,ICUBP)
      END DO
      END IF
C     
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      DO 230 JDFL=1,IDFL
        IG=KDFG(JDFL)
        IL=KDFL(JDFL)
        HBAS=DBAS(1,IL,1)
        DG1(IG)=DG1(IG) + OM*DENSITY*DG(1)*HBAS*TSTEP
        DG2(IG)=DG2(IG) + OM*DENSITY*DG(2)*HBAS*TSTEP
        DG3(IG)=DG3(IG) + OM*DENSITY*DG(3)*HBAS*TSTEP
230   CONTINUE
C
200   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE IntQ1toP1(DP1,DQ1,KVERT,KAREA,KEDGE,
     *           DCORVG,ELE)
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
      DIMENSION DP1(*),DQ1(*)
      REAL*8 dmyBAS(4),DA(4,4),DIA(4,4),DINTEGRAL(4),DSOL(4)
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
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
      BDER(1)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      ICUB=7
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
      DX0I = 0d0
      DY0I = 0d0
      DZ0I = 0d0
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
      DX0I = DX0I + 0.125d0*DCORVG(1,JP)
      DY0I = DY0I + 0.125d0*DCORVG(2,JP)
      DZ0I = DZ0I + 0.125d0*DCORVG(3,JP)
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
      DINTEGRAL = 0d0
      DA = 0d0
C
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

      dmyBAS(1) = 1d0
      dmyBAS(2) = XX - DX0I
      dmyBAS(3) = YY - DY0I
      dmyBAS(4) = ZZ - DZ0I
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      dVal = 0d0
      DO JDFL=1,IDFL
       IG=KDFG(JDFL)
       IL=KDFL(JDFL)
       dVal = dVal + DBAS(1,IL,1)*DQ1(IG)
      END DO
C
!       WRITE(*,*) icubp,dVal 
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,4
       HBASJ = dmyBAS(JDOFE)
       DO 240 IDOFE=1,4
        HBASI = dmyBAS(IDOFE)
        AH = HBASJ*HBASI
        DA(JDOFE,IDOFE)=DA(JDOFE,IDOFE) + OM*AH
240    CONTINUE
230   CONTINUE
C
      DO 250 JDFL=1,4
        HBAS=dmyBAS(JDFL)
        DINTEGRAL(JDFL) = DINTEGRAL(JDFL) + OM*dVal*HBAS
250   CONTINUE
C
200   CONTINUE
C
      CALL FINDInv(DA, DIA, 4, bflag)
      DSOL = 0d0
      DO JDOFE=1,4
       DO IDOFE=1,4
        DSOL(JDOFE) = DSOL(JDOFE) + DIA(JDOFE,IDOFE)*DINTEGRAL(IDOFE)
       END DO
      END DO
C
      IIII = 4*(IEL-1)
      DP1(IIII+1:IIII+4) = DSOL
C
!       WRITE(*,*) bflag
!       WRITE(*,'(8ES10.2)') dQ1(KVERT(:,iel))
!       WRITE(*,'(20ES10.2)') dsol
!       WRITE(*,'(20ES10.2)') DINTEGRAL
!       WRITE(*,'(20ES10.2)') DA
!       WRITE(*,'(20ES10.2)') DIA
!       pause
C
100   CONTINUE
C
99999 END


************************************************************************
      SUBROUTINE SubIntegrateShearrate(DU1,DU2,DU3,KVERT,KAREA,KEDGE,
     *           DCORVG,dAverageShearRate,ELE)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,COMM_SUMM
      USE var_QuadScalar, ONLY : transform
C     
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C     --------------------------- Transformation -------------------------------
      REAL*8    DHELP_Q2(27,4,NNCUBP),DHELP_Q1(8,4,NNCUBP)
      REAL*8    DPP(3)
C     --------------------------------------------------------------------------
C
      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
      DIMENSION DU3(NNBAS), GRADU3(NNDIM)
      
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
      dVolumeInt = 0d0
      dSQSHRInt  = 0d0
      
      if (myid.ne.0) then
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
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
       CALL E013A(XI1,XI2,XI3,DHELP_Q2,ICUBP)
      END DO
C
      BDER(2:4)=.FALSE.
C     
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
      ILINT = 1
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
       CALL E013A(XI1,XI2,XI3,DHELP_Q2,ICUBP)
      END DO
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C      
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
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
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      IF (Transform%ILINT.eq.2) THEN ! Q2
      DO JDOFE=1,27
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q2(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q2(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q2(JDFL,4,ICUBP)
      END DO
      END IF
      IF (Transform%ILINT.eq.1) THEN ! Q1
      DO JDOFE=1,8
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q1(JDFL,4,ICUBP)
      END DO
      END IF
C     
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C     
C ---=========================---
      GRADU1(1)=0D0!U
      GRADU1(2)=0D0
      GRADU1(3)=0D0
C
      GRADU2(1)=0D0!V
      GRADU2(2)=0D0
      GRADU2(3)=0D0
C
      GRADU3(1)=0D0!W
      GRADU3(2)=0D0
      GRADU3(3)=0D0 
C
      DO 220 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)! local number of basic function
       JDFG=KDFG(JDOFE)! local number of basic function

       GRADU1(1)=GRADU1(1) + DU1(JDFG)*DBAS(1,JDFL,2)!DUX
       GRADU1(2)=GRADU1(2) + DU1(JDFG)*DBAS(1,JDFL,3)!DUY
       GRADU1(3)=GRADU1(3) + DU1(JDFG)*DBAS(1,JDFL,4)!DUZ

       GRADU2(1)=GRADU2(1) + DU2(JDFG)*DBAS(1,JDFL,2)!DVX
       GRADU2(2)=GRADU2(2) + DU2(JDFG)*DBAS(1,JDFL,3)!DVY
       GRADU2(3)=GRADU2(3) + DU2(JDFG)*DBAS(1,JDFL,4)!DVZ

       GRADU3(1)=GRADU3(1) + DU3(JDFG)*DBAS(1,JDFL,2)!DWX
       GRADU3(2)=GRADU3(2) + DU3(JDFG)*DBAS(1,JDFL,3)!DWY
       GRADU3(3)=GRADU3(3) + DU3(JDFG)*DBAS(1,JDFL,4)!DWZ

 220  CONTINUE

C ----=============================================---- 
       dSQSHR = GRADU1(1)**2d0 + GRADU2(2)**2d0 + GRADU3(3)**2d0
     *        + 0.5d0*(GRADU1(2)+GRADU2(1))**2d0
     *        + 0.5d0*(GRADU1(3)+GRADU3(1))**2d0 
     *        + 0.5d0*(GRADU2(3)+GRADU3(2))**2d0
C ----=============================================---- 
C
C *** Summing up over all pairs of multiindices
      DO 230 JDFL=1,IDFL
        IG=KDFG(JDFL)
        IL=KDFL(JDFL)
        HBAS=DBAS(1,IL,1)
        dVolumeInt = dVolumeInt + OM*1d0*HBAS
        dSQSHRInt  = dSQSHRInt  + OM*dSQSHR*HBAS
230   CONTINUE 
C
200   CONTINUE
C
100   CONTINUE
C
      END IF
C      
      CALL COMM_SUMM(dVolumeInt)
      CALL COMM_SUMM(dSQSHRInt)
C      
      dAverageShearRate = (2d0*dSQSHRInt/dVolumeInt)**0.5d0
      
      if (myid.eq.1) then
       write(*,*) 'dAverageShearRate: ',dAverageShearRate,dVolumeInt
      end if
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE HEXProjection(icub,dcbV,DD,DM,KVERT,KAREA,KEDGE,
     *           DCORVG,ELE)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid
      USE var_QuadScalar, ONLY : transform
C     
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION dcbV(3,*),DD(3,*),DM(*)
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C     --------------------------- Transformation -------------------------------
      REAL*8    DHELP_Q2(27,4,NNCUBP),DHELP_Q1(8,4,NNCUBP)
      REAL*8    DPP(3)
C     --------------------------------------------------------------------------
      REAL*8 Velo(3)
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
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
      END DO
C
C *** Loop over all elements
      DO 100 IEL=1,NEL

      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
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
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      DO JDOFE=1,8
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q1(JDFL,4,ICUBP)
      END DO
C     
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      
      Velo = dcbV(:,(iel-1)*ncubp + icubp)
      
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      DO 230 JDFL=1,IDFL
        IG=KDFG(JDFL)
        IL=KDFL(JDFL)
        HBAS=DBAS(1,IL,1)
        DD(:,IG) = DD(:,IG) + OM*Velo(:)*HBAS
        DM(IG) = DM(IG) + OM*HBAS
230   CONTINUE
C
200   CONTINUE
C
100   CONTINUE

!        write(*,*) DD(1:3,1:27)
!        write(*,*) 
!        write(*,*) DM(1:3,1:27)
!        pause
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE SubIntBigU(DU1,DU2,DU3,KVERT,KAREA,KEDGE,
     *           DCORVG,dAverageShearRate,ELE,iComm)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,COMM_SUMM
      USE var_QuadScalar, ONLY : transform
C     
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C     --------------------------- Transformation -------------------------------
      REAL*8    DHELP_Q2(27,4,NNCUBP),DHELP_Q1(8,4,NNCUBP)
      REAL*8    DPP(3)
C     --------------------------------------------------------------------------
C
      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
      DIMENSION DU3(NNBAS), GRADU3(NNDIM)
      
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
      dVolumeInt = 0d0
      dSQSHRInt  = 0d0
      
      if (myid.ne.0) then
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
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
      END DO
C
C *** Loop over all elements

      DO 100 IEL=1,NEL
C      
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
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
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      DO JDOFE=1,8
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q1(JDFL,4,ICUBP)
      END DO
C     
      
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C     
C ---=========================---
      GRADU1(1)=0D0!U
      GRADU1(2)=0D0
      GRADU1(3)=0D0
C
      GRADU2(1)=0D0!V
      GRADU2(2)=0D0
      GRADU2(3)=0D0
C
      GRADU3(1)=0D0!W
      GRADU3(2)=0D0
      GRADU3(3)=0D0 
C
      DO 220 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)! local number of basic function
       JDFG=KDFG(JDOFE)! local number of basic function

       GRADU1(1)=GRADU1(1) + DU1(JDFG)*DBAS(1,JDFL,2)!DUX
       GRADU1(2)=GRADU1(2) + DU1(JDFG)*DBAS(1,JDFL,3)!DUY
       GRADU1(3)=GRADU1(3) + DU1(JDFG)*DBAS(1,JDFL,4)!DUZ

       GRADU2(1)=GRADU2(1) + DU2(JDFG)*DBAS(1,JDFL,2)!DVX
       GRADU2(2)=GRADU2(2) + DU2(JDFG)*DBAS(1,JDFL,3)!DVY
       GRADU2(3)=GRADU2(3) + DU2(JDFG)*DBAS(1,JDFL,4)!DVZ

       GRADU3(1)=GRADU3(1) + DU3(JDFG)*DBAS(1,JDFL,2)!DWX
       GRADU3(2)=GRADU3(2) + DU3(JDFG)*DBAS(1,JDFL,3)!DWY
       GRADU3(3)=GRADU3(3) + DU3(JDFG)*DBAS(1,JDFL,4)!DWZ

 220  CONTINUE

 
C ----=============================================---- 
       dSQSHR = GRADU1(1)**2d0 + GRADU2(2)**2d0 + GRADU3(3)**2d0
     *        + 0.5d0*(GRADU1(2)+GRADU2(1))**2d0
     *        + 0.5d0*(GRADU1(3)+GRADU3(1))**2d0 
     *        + 0.5d0*(GRADU2(3)+GRADU3(2))**2d0
C ----=============================================---- 
C
      !! if sqrt(0.5*(D:D)) then this:
!         dSQSHR = (2d0*dSQSHR)**2d0
      !! otherwise
      dSQSHR = dSQSHR
C      
C *** Summing up over all pairs of multiindices
      DO 230 JDFL=1,IDFL
        IG=KDFG(JDFL)
        IL=KDFL(JDFL)
        HBAS=DBAS(1,IL,1)
        dVolumeInt = dVolumeInt + OM*1d0*HBAS
        dSQSHRInt  = dSQSHRInt  + OM*dSQSHR*HBAS
230   CONTINUE 
C
200   CONTINUE
C
100   CONTINUE
C
      END IF
C      
      if (iComm.eq.1) THEN
       CALL COMM_SUMM(dVolumeInt)
       CALL COMM_SUMM(dSQSHRInt)
      END IF
! C      
!       dAverageShearRate = dSQSHRInt
      
      if (myid.eq.1) then
       write(*,'(A,2ES14.4)') 'Integral D(BigU):D(BigU) =',
     *                         dSQSHRInt,dVolumeInt
      end if
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE SubIntSigmaP(DU1,DU2,DU3, smallU1, smallU2, smallU3,
     *           KVERT,KAREA,KEDGE,
     *           DCORVG,dAverageShearRate,ELE,iComm)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,COMM_SUMM
      USE var_QuadScalar, ONLY : transform
C     
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C     --------------------------- Transformation -------------------------------
      REAL*8    DHELP_Q2(27,4,NNCUBP),DHELP_Q1(8,4,NNCUBP)
      REAL*8    DPP(3)
C     --------------------------------------------------------------------------
C
      ! The velocity BigU on the coarsened 
      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
      DIMENSION DU3(NNBAS), GRADU3(NNDIM)

      ! The velocity on the fine mesh
      DIMENSION smallU1(NNBAS), GRADsmallU1(NNDIM)
      DIMENSION smallU2(NNBAS), GRADsmallU2(NNDIM)
      DIMENSION smallU3(NNBAS), GRADsmallU3(NNDIM)
      
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
      dVolumeInt = 0d0
      dSQSHRInt  = 0d0
      DDuUInt  = 0d0
      
      if (myid.ne.0) then
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
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
      END DO
C
C *** Loop over all elements

      DO 100 IEL=1,NEL
C      
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
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
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      DO JDOFE=1,8
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q1(JDFL,4,ICUBP)
      END DO
C     
      
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C     
C ---=========================---
      GRADU1(1)=0D0!U
      GRADU1(2)=0D0
      GRADU1(3)=0D0
C
      GRADU2(1)=0D0!V
      GRADU2(2)=0D0
      GRADU2(3)=0D0
C
      GRADU3(1)=0D0!W
      GRADU3(2)=0D0
      GRADU3(3)=0D0 
C
      GRADsmallU1(1)=0D0!U
      GRADsmallU1(2)=0D0
      GRADsmallU1(3)=0D0
C
      GRADsmallU2(1)=0D0!V
      GRADsmallU2(2)=0D0
      GRADsmallU2(3)=0D0
C
      GRADsmallU3(1)=0D0!W
      GRADsmallU3(2)=0D0
      GRADsmallU3(3)=0D0 
C
      DO 220 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)! local number of basic function
       JDFG=KDFG(JDOFE)! local number of basic function

       GRADU1(1)=GRADU1(1) + DU1(JDFG)*DBAS(1,JDFL,2)!DUX
       GRADU1(2)=GRADU1(2) + DU1(JDFG)*DBAS(1,JDFL,3)!DUY
       GRADU1(3)=GRADU1(3) + DU1(JDFG)*DBAS(1,JDFL,4)!DUZ

       GRADU2(1)=GRADU2(1) + DU2(JDFG)*DBAS(1,JDFL,2)!DVX
       GRADU2(2)=GRADU2(2) + DU2(JDFG)*DBAS(1,JDFL,3)!DVY
       GRADU2(3)=GRADU2(3) + DU2(JDFG)*DBAS(1,JDFL,4)!DVZ

       GRADU3(1)=GRADU3(1) + DU3(JDFG)*DBAS(1,JDFL,2)!DWX
       GRADU3(2)=GRADU3(2) + DU3(JDFG)*DBAS(1,JDFL,3)!DWY
       GRADU3(3)=GRADU3(3) + DU3(JDFG)*DBAS(1,JDFL,4)!DWZ

       GRADsmallU1(1)=GRADsmallU1(1) + smallU1(JDFG)*DBAS(1,JDFL,2)!DUX
       GRADsmallU1(2)=GRADsmallU1(2) + smallU1(JDFG)*DBAS(1,JDFL,3)!DUY
       GRADsmallU1(3)=GRADsmallU1(3) + smallU1(JDFG)*DBAS(1,JDFL,4)!DUZ

       GRADsmallU2(1)=GRADsmallU2(1) + smallU2(JDFG)*DBAS(1,JDFL,2)!DVX
       GRADsmallU2(2)=GRADsmallU2(2) + smallU2(JDFG)*DBAS(1,JDFL,3)!DVY
       GRADsmallU2(3)=GRADsmallU2(3) + smallU2(JDFG)*DBAS(1,JDFL,4)!DVZ

       GRADsmallU3(1)=GRADsmallU3(1) + smallU3(JDFG)*DBAS(1,JDFL,2)!DWX
       GRADsmallU3(2)=GRADsmallU3(2) + smallU3(JDFG)*DBAS(1,JDFL,3)!DWY
       GRADsmallU3(3)=GRADsmallU3(3) + smallU3(JDFG)*DBAS(1,JDFL,4)!DWZ

       GRADsmallU1(1)=GRADsmallU1(1) - GRADU1(1)
       GRADsmallU1(2)=GRADsmallU1(2) - GRADU1(2)
       GRADsmallU1(3)=GRADsmallU1(3) - GRADU1(3)
                                                
       GRADsmallU2(1)=GRADsmallU2(1) - GRADU2(1)
       GRADsmallU2(2)=GRADsmallU2(2) - GRADU2(2)
       GRADsmallU2(3)=GRADsmallU2(3) - GRADU2(3)
                                                
       GRADsmallU3(1)=GRADsmallU3(1) - GRADU3(1)
       GRADsmallU3(2)=GRADsmallU3(2) - GRADU3(2)
       GRADsmallU3(3)=GRADsmallU3(3) - GRADU3(3)

 220  CONTINUE

C ----=============================================---- 
       DDuU = 1d0 * GRADsmallU1(1) * GradU1(1)
     *       + 1d0 * GRADsmallU2(2) * GradU2(2) 
     *       + 1d0 * GRADsmallU3(3) * GradU3(3)
     *       + 0.5d0 * (GRADsmallU1(2) + GRADsmallU2(1))
     *             * (GRADU1(2) + GRADU2(1)) 
     *       + 0.5d0 * (GRADsmallU1(3) + GRADsmallU3(1)) 
     *             * (GRADU1(3) + GRADU3(1))
     *       + 0.5d0 * (GRADsmallU2(3) + GRADsmallU3(2)) 
     *             * (GRADU2(3) + GRADU3(2)) 
C ----=============================================---- 

 
C ----=============================================---- 
       dSQSHR = GRADU1(1)**2d0 + GRADU2(2)**2d0 + GRADU3(3)**2d0
     *        + 0.5d0*(GRADU1(2)+GRADU2(1))**2d0
     *        + 0.5d0*(GRADU1(3)+GRADU3(1))**2d0 
     *        + 0.5d0*(GRADU2(3)+GRADU3(2))**2d0
C ----=============================================---- 
C
      !! if sqrt(0.5*(D:D)) then this:
!         dSQSHR = (2d0*dSQSHR)**2d0
      !! otherwise
      dSQSHR = dSQSHR
C      
C *** Summing up over all pairs of multiindices
      DO 230 JDFL=1,IDFL
        IG=KDFG(JDFL)
        IL=KDFL(JDFL)
        HBAS=DBAS(1,IL,1)
        dVolumeInt = dVolumeInt + OM*1d0*HBAS
        dSQSHRInt  = dSQSHRInt  + OM*dSQSHR*HBAS
        DDuUInt = DDuUInt + OM * DDuU * HBAS
230   CONTINUE 
C
200   CONTINUE
C
100   CONTINUE
C
      END IF
C      
      if (iComm.eq.1) THEN
       CALL COMM_SUMM(dVolumeInt)
       CALL COMM_SUMM(dSQSHRInt)
       CALL COMM_SUMM(DDuUInt)
      END IF
! C      
!       dAverageShearRate = dSQSHRInt
      
      if (myid.eq.1) then
       write(*,'(A,2ES14.4)') 'Integral D(BigU):D(BigU) =',
     *                         dSQSHRInt,dVolumeInt
       write(*,'(A,1ES14.4)') 'Integral D(u-BigU):D(BigU) =', DDuUInt
      end if
C
99999 END
