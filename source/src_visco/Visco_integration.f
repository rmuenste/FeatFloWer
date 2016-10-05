! ************************************************************************
!       SUBROUTINE ViscoDriver_Mat(U1,U2,U3,T11,T22,T33,T12,T13,T23,
!      *           DS11,DS22,DS33,DS12,DS13,DS23,
!      *           NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,DCORVG,ELE,dLm,dt)
! ************************************************************************
! *
! *-----------------------------------------------------------------------
!       USE PP3D_MPI, ONLY:myid
! C
!       IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
!       CHARACTER SUB*6,FMT*15,CPARAM*120
! C
!       PARAMETER (DEPS=1d-15)
!       PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
!      *           NNDIM=3,NNCOF=10)
!       PARAMETER (Q2=0.5D0,Q8=0.125D0)
! C
!       REAL*8    T11(*),T22(*),T33(*),T12(*),T13(*),T23(*)
!       REAL*8    DS11(*),DS22(*),DS33(*),DS12(*),DS13(*),DS23(*)
!       REAL*8    U1(*),U2(*),U3(*)
! C
!       DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
!       DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
!       DIMENSION KENTRY(NNBAS,NNBAS)
!       REAL*8    S11(NNBAS,NNBAS),S22(NNBAS,NNBAS),S33(NNBAS,NNBAS)
!       REAL*8    S12(NNBAS,NNBAS),S13(NNBAS,NNBAS),S23(NNBAS,NNBAS)
! C
!       REAL*8 AMAT(3,3),QMAT(3,3),EVCT(3),OMAT(3,3),BMAT(3,3)
!       REAL*8 GRADU(3,3),TAU(6)
!      
!       DIMENSION KDFG(NNBAS),KDFL(NNBAS)
!       DIMENSION DU1(NNBAS), GRADU1(NNDIM)
!       DIMENSION DU2(NNBAS), GRADU2(NNDIM)
!       DIMENSION DU3(NNBAS), GRADU3(NNDIM)
! C
!       COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
!       COMMON /ERRCTL/ IER,ICHECK
!       COMMON /CHAR/   SUB,FMT(3),CPARAM
!       COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
!      *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
!      *                IEL,NDIM
!       COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
!      *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
!       COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
!       COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
! 
!       COMMON /COAUX1/ KDFG,KDFL,IDFL
! C
! C *** user COMMON blocks
!       INTEGER  VIPARM 
!       DIMENSION VIPARM(100)
!       EQUIVALENCE (IAUSAV,VIPARM)
!       COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
!      *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
!      *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
!      *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
!      *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
!      *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
! C
!       SAVE
! C
!       DO 1 I= 1,NNDER
! 1     BDER(I)=.FALSE.
! C
!       DO 2 I=1,4
! 2     BDER(I)=.TRUE.
! C
!       IELTYP=-1
!       CALL ELE(0D0,0D0,0D0,IELTYP)
!       IDFL=NDFL(IELTYP)
! C
!       ICUB=9
!       CALL CB3H(ICUB)
!       IF (IER.NE.0) GOTO 99999
! C
! ************************************************************************
! C *** Calculation of the matrix - storage technique 7 or 8
! ************************************************************************
!       ICUBP=ICUB
!       CALL ELE(0D0,0D0,0D0,-2)
! C
! C *** Loop over all elements
!       DO 100 IEL=1,NEL
! C
!       CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
!       IF (IER.LT.0) GOTO 99999
! C
! C *** Determine entry positions in matrix
!       DO 110 JDOFE=1,IDFL
!       ILD=KLDA(KDFG(JDOFE))
!       KENTRY(JDOFE,JDOFE)=ILD
!       S11(JDOFE,JDOFE)=0D0
!       S22(JDOFE,JDOFE)=0D0
!       S33(JDOFE,JDOFE)=0D0
!       S12(JDOFE,JDOFE)=0D0
!       S13(JDOFE,JDOFE)=0D0
!       S23(JDOFE,JDOFE)=0D0
!       JCOL0=ILD
!       DO 111 IDOFE=1,IDFL
!       IF (IDOFE.EQ.JDOFE) GOTO 111
!       IDFG=KDFG(IDOFE)
!       DO 112 JCOL=JCOL0,NA
!       IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
! 112   CONTINUE
! 113   JCOL0=JCOL+1
!       KENTRY(JDOFE,IDOFE)=JCOL
!       S11(JDOFE,IDOFE)=0D0
!       S22(JDOFE,IDOFE)=0D0
!       S33(JDOFE,IDOFE)=0D0
!       S12(JDOFE,IDOFE)=0D0
!       S13(JDOFE,IDOFE)=0D0
!       S23(JDOFE,IDOFE)=0D0
! 111   CONTINUE
! 110   CONTINUE
! C
!       DO 120 IVE=1,NVE
!       JP=KVERT(IVE,IEL)
!       KVE(IVE)=JP
!       DX(IVE)=DCORVG(1,JP)
!       DY(IVE)=DCORVG(2,JP)
!       DZ(IVE)=DCORVG(3,JP)
! 120   CONTINUE
! C
!       DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
!       DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
!       DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
!       DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
!       DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
!       DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
!       DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
!       DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
!       DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
!       DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
!       DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
!       DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
!       DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
!       DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
!       DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
!       DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
!       DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
!       DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
!       DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
!       DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
!       DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
!       DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
!       DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
!       DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
! C
!       DO 130 JDOFE=1,IDFL
!       JDFG=KDFG(JDOFE)
!       JDFL=KDFL(JDOFE)
!       DU1(JDFL) = U1(JDFG) 
!       DU2(JDFL) = U2(JDFG)
!       DU3(JDFL) = U3(JDFG)
!  130  CONTINUE      
! C
! C *** Loop over all cubature points
!       DO 200 ICUBP=1,NCUBP
! C
!       XI1=DXI(ICUBP,1)
!       XI2=DXI(ICUBP,2)
!       XI3=DXI(ICUBP,3)
! C
! C *** Jacobian of the bilinear mapping onto the reference element
!       DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
!       DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
!       DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
!       DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
!       DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
!       DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
!       DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
!       DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
!       DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
!       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
!      *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
!      *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
!       OM=DOMEGA(ICUBP)*ABS(DETJ)
! C
!        XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
!        YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
!        ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
! 
!       CALL ELE(XI1,XI2,XI3,-3)
!       IF (IER.LT.0) GOTO 99999
! C 
!       GRADU1 = 0D0
!       GRADU2 = 0D0
!       GRADU3 = 0D0 
!       TAU    = 0d0
! C
!       DO JDOFE=1,IDFL
!        IG = KDFG(JDOFE)
!        IL = KDFL(JDOFE)
!        GRADU(1,1) = GRADU(1,1) + U1(IG)*DBAS(1,IL,2)
!        GRADU(1,2) = GRADU(1,2) + U1(IG)*DBAS(1,IL,3)
!        GRADU(1,3) = GRADU(1,3) + U1(IG)*DBAS(1,IL,4)
!        GRADU(2,1) = GRADU(2,1) + U2(IG)*DBAS(1,IL,2)
!        GRADU(2,2) = GRADU(2,2) + U2(IG)*DBAS(1,IL,3)
!        GRADU(2,3) = GRADU(2,3) + U2(IG)*DBAS(1,IL,4)
!        GRADU(3,1) = GRADU(3,1) + U3(IG)*DBAS(1,IL,2)
!        GRADU(3,2) = GRADU(3,2) + U3(IG)*DBAS(1,IL,3)
!        GRADU(3,3) = GRADU(3,3) + U3(IG)*DBAS(1,IL,4)
!        TAU(1)    = TAU(1)    + T11(IG)*DBAS(1,IL,1)
!        TAU(2)    = TAU(2)    + T22(IG)*DBAS(1,IL,1)
!        TAU(3)    = TAU(3)    + T33(IG)*DBAS(1,IL,1)
!        TAU(4)    = TAU(4)    + T12(IG)*DBAS(1,IL,1)
!        TAU(5)    = TAU(5)    + T13(IG)*DBAS(1,IL,1)
!        TAU(6)    = TAU(6)    + T23(IG)*DBAS(1,IL,1)
!       END DO
! C
!       TAU(:) = [-4d0,1d0,3d0,1.1d0,-0.2d0,0.01d0]
!       AMAT(1,:) = [TAU(1),TAU(4),TAU(5)]
!       AMAT(2,:) = [TAU(4),TAU(2),TAU(6)]
!       AMAT(3,:) = [TAU(5),TAU(6),TAU(3)]
! C
!       CALL DSYEVD3(AMAT,QMAT,EVCT)
! C
!       CALL GetOtherMatrices(GRADU,QMAT,EVCT,OMAT,BMAT)
! C
!       IF (myid.eq.1) WRITE(*,*) "AMAT"
!       IF (myid.eq.1) WRITE(*,*) AMAT
!       IF (myid.eq.1) WRITE(*,*) "QMAT"
!       IF (myid.eq.1) WRITE(*,*) QMAT
!       IF (myid.eq.1) WRITE(*,*) "EVCT"
!       IF (myid.eq.1) WRITE(*,*) EVCT
!       IF (myid.eq.1) WRITE(*,*) "BMAT"
!       IF (myid.eq.1) WRITE(*,*) BMAT
!       IF (myid.eq.1) WRITE(*,*) "OMAT"
!       IF (myid.eq.1) WRITE(*,*) OMAT
!       IF (DABS(EVCT(1)-EVCT(2)).LT.DEPS) OMAT = 0d0
!       IF (DABS(EVCT(1)-EVCT(3)).LT.DEPS) OMAT = 0d0
!       IF (DABS(EVCT(2)-EVCT(3)).LT.DEPS) OMAT = 0d0
!       IF (myid.eq.1) WRITE(*,*) "OMAT"
!       IF (myid.eq.1) WRITE(*,*) OMAT
!       PAUSE
! 
!       dLam = dLm
! C
!       DO 230 JDOFE=1,IDFL
! 
!        JDOFEH = KDFL(JDOFE)
!        HBASJ1 = DBAS(1,JDOFEH,1)
! 
!        DO 240 IDOFE=1,IDFL
! 
!         IDOFEH = KDFL(IDOFE)
!         HBASI1 = DBAS(1,IDOFEH,1)
! 
!         AH   = HBASI1*HBASJ1
!         AH12 = 2d0*OMAT(1,2)*AH
!         AH13 = ((GRADU1(1) + GRADU3(3)))*AH
!         AH23 = ((GRADU2(2) + GRADU3(3)))*AH
! 
!         S11(JDOFE,IDOFE) = S11(JDOFE,IDOFE) + OM*AH11
!         S22(JDOFE,IDOFE) = S22(JDOFE,IDOFE) + OM*AH22
!         S33(JDOFE,IDOFE) = S33(JDOFE,IDOFE) + OM*AH33
!         S12(JDOFE,IDOFE) = S12(JDOFE,IDOFE) + OM*AH12
!         S13(JDOFE,IDOFE) = S13(JDOFE,IDOFE) + OM*AH13
!         S23(JDOFE,IDOFE) = S23(JDOFE,IDOFE) + OM*AH23
!  240  CONTINUE
!  230  CONTINUE
! C
!  200  CONTINUE
! C
!       DO 300 JDOFE=1,IDFL
!       JDOFEH = KDFL(JDOFE)
!       DO 300 IDOFE=1,IDFL
!         IDOFEH = KDFL(IDOFE)
! !        IA = KENTRY(JDOFEH,IDOFEH) ?????????????????
!         IA = KENTRY(JDOFE,IDOFE)
!         DS11(IA) = DS11(IA) + dt*S11(JDOFE,IDOFE)
!         DS22(IA) = DS22(IA) + dt*S22(JDOFE,IDOFE)
!         DS33(IA) = DS33(IA) + dt*S33(JDOFE,IDOFE)
!         DS12(IA) = DS12(IA) + dt*S12(JDOFE,IDOFE)
!         DS13(IA) = DS13(IA) + dt*S13(JDOFE,IDOFE)
!         DS23(IA) = DS23(IA) + dt*S23(JDOFE,IDOFE)
!  300  CONTINUE
! C
!  100  CONTINUE
! C
! 99999 END
! C
! C
! C
************************************************************************
      SUBROUTINE ViscoDriver_Def(U1,U2,U3,DCOOR,
     *           U11,U22,U33,U12,U13,U23,D11,
     *           D22,D33,D12,D13,D23,KVERT,KAREA,
     *           KEDGE,DCORVG,ELE,dLm,dt)
************************************************************************
*
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,myMPI_Barrier
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (DEPS=1d-13)
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION U1(*),U2(*),U3(*),DCOOR(3,*)
      DIMENSION U11(*),U22(*),U33(*),U12(*),U13(*),U23(*)
      DIMENSION D11(*),D22(*),D33(*),D12(*),D13(*),D23(*)
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
      REAL*8 DU1(NNBAS)
      REAL*8 DU2(NNBAS)
      REAL*8 DU3(NNBAS)

      REAL*8 AMAT(3,3),QMAT(3,3),EVCT(3),OMAT(3,3),BMAT(3,3)
      REAL*8 GRADU(3,3),TAU(6),RMAT(3,3)

      REAL*8 DEF1(NNBAS),DEF4(NNBAS),DUU1(NNBAS),DUU4(NNBAS)
      REAL*8 DEF2(NNBAS),DEF5(NNBAS),DUU2(NNBAS),DUU5(NNBAS)
      REAL*8 DEF3(NNBAS),DEF6(NNBAS),DUU3(NNBAS),DUU6(NNBAS)
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
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS

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
      DXXS = 0d0
      DXXC = 0d0
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
      DO JDOFE=1,IDFL
       JDFG=KDFG(JDOFE)
       DU1(JDOFE)  = U1(JDFG)
       DU2(JDOFE)  = U2(JDFG)
       DU3(JDOFE)  = U3(JDFG)
       DUU1(JDOFE) = U11(JDFG)
       DUU2(JDOFE) = U22(JDFG)
       DUU3(JDOFE) = U33(JDFG)
       DUU4(JDOFE) = U12(JDFG)
       DUU5(JDOFE) = U13(JDFG)
       DUU6(JDOFE) = U23(JDFG)
       DEF1(JDOFE) = 0D0
       DEF2(JDOFE) = 0D0
       DEF3(JDOFE) = 0D0
       DEF4(JDOFE) = 0D0
       DEF5(JDOFE) = 0D0
       DEF6(JDOFE) = 0D0
      END DO
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
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      TAU    = 0d0
      GRADU = 0D0
C
      DO JDOFE=1,IDFL
       IG = KDFG(JDOFE)
       IL = KDFL(JDOFE)
       GRADU(1,1) = GRADU(1,1) + U1(IG)*DBAS(1,IL,2)
       GRADU(1,2) = GRADU(1,2) + U1(IG)*DBAS(1,IL,3)
       GRADU(1,3) = GRADU(1,3) + U1(IG)*DBAS(1,IL,4)
       GRADU(2,1) = GRADU(2,1) + U2(IG)*DBAS(1,IL,2)
       GRADU(2,2) = GRADU(2,2) + U2(IG)*DBAS(1,IL,3)
       GRADU(2,3) = GRADU(2,3) + U2(IG)*DBAS(1,IL,4)
       GRADU(3,1) = GRADU(3,1) + U3(IG)*DBAS(1,IL,2)
       GRADU(3,2) = GRADU(3,2) + U3(IG)*DBAS(1,IL,3)
       GRADU(3,3) = GRADU(3,3) + U3(IG)*DBAS(1,IL,4)
       TAU(1)    = TAU(1)    + U11(IG)*DBAS(1,IL,1)
       TAU(2)    = TAU(2)    + U22(IG)*DBAS(1,IL,1)
       TAU(3)    = TAU(3)    + U33(IG)*DBAS(1,IL,1)
       TAU(4)    = TAU(4)    + U12(IG)*DBAS(1,IL,1)
       TAU(5)    = TAU(5)    + U13(IG)*DBAS(1,IL,1)
       TAU(6)    = TAU(6)    + U23(IG)*DBAS(1,IL,1)
      END DO
C
!      TAU(:) = [-4d0,1d0,3d0,1.1d0,-0.2d0,0.01d0]
C
      AMAT(1,:) = [TAU(1),TAU(4),TAU(5)]
      AMAT(2,:) = [TAU(4),TAU(2),TAU(6)]
      AMAT(3,:) = [TAU(5),TAU(6),TAU(3)]
C
      CALL DSYEVD3(AMAT,QMAT,EVCT)
C
      CALL GetOtherMatrices(GRADU,QMAT,EVCT,OMAT,BMAT,RMAT)
C
!       IF (myid.eq.1) WRITE(*,*) "AMAT"
!       IF (myid.eq.1) WRITE(*,*) AMAT
!       IF (myid.eq.1) WRITE(*,*) "QMAT"
!       IF (myid.eq.1) WRITE(*,*) QMAT
!       IF (myid.eq.1) WRITE(*,*) "EVCT"
!       IF (myid.eq.1) WRITE(*,*) EVCT
!       IF (myid.eq.1) WRITE(*,*) "BMAT"
!       IF (myid.eq.1) WRITE(*,*) BMAT
!       IF (myid.eq.1) WRITE(*,*) "RMAT"
!       IF (myid.eq.1) WRITE(*,*) RMAT
!       IF (myid.eq.1) WRITE(*,*) "OMAT"
!       IF (myid.eq.1) WRITE(*,*) OMAT
      IF (DABS(EVCT(1)-EVCT(2)).LT.DEPS) OMAT = 0d0
      IF (DABS(EVCT(1)-EVCT(3)).LT.DEPS) OMAT = 0d0
      IF (DABS(EVCT(2)-EVCT(3)).LT.DEPS) OMAT = 0d0
!       IF (myid.eq.1) WRITE(*,*) "criterions"
!       IF (myid.eq.1) WRITE(*,*) DABS(EVCT(1)-EVCT(2)),
!      * DABS(EVCT(2)-EVCT(3)),DABS(EVCT(1)-EVCT(3)),DEPS
!       IF (myid.eq.1) WRITE(*,*) "OMAT"
!       IF (myid.eq.1) WRITE(*,*) OMAT
!       pause
      dLam = dLm
C     
      DO 230 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       HBAS = DBAS(1,JDFL,1)

       DEF1(JDOFE) = DEF1(JDOFE) + (
     * 2d0*(OMAT(1,2)*AMAT(2,1) + OMAT(1,3)*AMAT(3,1) +
     * BMAT(1,1)) + (RMAT(1,1) - 1d0)/dLam)*OM*HBAS

       DEF2(JDOFE) = DEF2(JDOFE) + (
     * 2d0*(OMAT(2,1)*AMAT(2,1) + OMAT(2,3)*AMAT(2,3) + 
     * BMAT(2,2)) + (RMAT(2,2) - 1d0)/dLam)*OM*HBAS

       DEF3(JDOFE) = DEF3(JDOFE) + (
     * 2d0*(OMAT(3,1)*AMAT(3,1) + OMAT(3,2)*AMAT(3,2) + 
     * BMAT(3,3)) + (RMAT(3,3) - 1d0)/dLam)*OM*HBAS

       DEF4(JDOFE) = DEF4(JDOFE) + (
     * OMAT(1,2)*(AMAT(2,2) - AMAT(1,1)) +
     * OMAT(1,3)*AMAT(3,2) - OMAT(3,2)*AMAT(1,3) +
     * 2d0*BMAT(1,2) + (RMAT(1,2) - 0d0)/dLam)*OM*HBAS

       DEF5(JDOFE) = DEF5(JDOFE) + (
     * OMAT(1,3)*(AMAT(3,3) - AMAT(1,1)) +
     * OMAT(1,2)*AMAT(2,3) - OMAT(2,3)*AMAT(1,2) +
     * 2d0*BMAT(1,3) + (RMAT(1,3) - 0d0)/dLam)*OM*HBAS

       DEF6(JDOFE) = DEF6(JDOFE) + (
     * OMAT(2,3)*(AMAT(3,3) - AMAT(2,2)) +
     * OMAT(2,1)*AMAT(1,3) - OMAT(1,3)*AMAT(2,1) +
     * 2d0*BMAT(2,3) + (RMAT(2,3) - 0d0)/dLam)*OM*HBAS

230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
       JDFG=KDFG(JDOFE)
       D11(JDFG) = D11(JDFG) + dt*DEF1(JDOFE)
       D22(JDFG) = D22(JDFG) + dt*DEF2(JDOFE)
       D33(JDFG) = D33(JDFG) + dt*DEF3(JDOFE)
       D12(JDFG) = D12(JDFG) + dt*DEF4(JDOFE)
       D13(JDFG) = D13(JDFG) + dt*DEF5(JDOFE)
       D23(JDFG) = D23(JDFG) + dt*DEF6(JDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
      SUBROUTINE GetOtherMatrices(GU,Q,E,O,B,R)
      USE PP3D_MPI, ONLY:myid,myMPI_Barrier
      IMPLICIT NONE
      REAL*8 Q(3,3),E(3),O(3,3),B(3,3),R(3,3),GU(3,3)
      REAL*8 A1(3,3),QT(3,3),M(3,3),EE(3)
      REAL*8 BB(3,3),OO(3,3),RR(3,3)
      REAL*8 DAUX,OMEGA12,OMEGA13,OMEGA23
      INTEGER I,J,K
            
      DO I=1,3
       EE(I) = EXP(E(I))
!      write(*,*)'EE(i)',EE(i) 
!      write(*,*)'EXP(E(I))',EXP(E(I)) 
       DO J=1,3
        QT(I,J) = Q(J,I)
       END DO
      END DO


      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + QT(I,K)*GU(K,J)
        END DO
        A1(I,J) = DAUX
       END DO
      END DO


      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + A1(I,K)*Q(K,J)
        END DO
        M(I,J) = DAUX
       END DO
      END DO

!      write(*,*)'EE(1)',EE(1) 
!      write(*,*)'EE(2)',EE(2) 
!      write(*,*)'EE(3)',EE(3) 


      OMEGA12 = (EE(2)*M(1,2) + EE(1)*M(2,1))/(EE(2) - EE(1))
      OMEGA13 = (EE(3)*M(1,3) + EE(1)*M(3,1))/(EE(3) - EE(1))
      OMEGA23 = (EE(3)*M(2,3) + EE(2)*M(3,2))/(EE(3) - EE(2))

!      write(*,*)'omega12',OMEGA12 
!      write(*,*)'omega13',OMEGA13 
!      write(*,*)'omega23',OMEGA23 
!
!       call myMPI_Barrier()
!       stop

      BB = 0d0
      BB(1,1) = M(1,1)
      BB(2,2) = M(2,2)
      BB(3,3) = M(3,3)

      RR = 0d0
      RR(1,1) = 1d0/EE(1)
      RR(2,2) = 1d0/EE(2)
      RR(3,3) = 1d0/EE(3)

      OO(1,:) = [ 0d0, OMEGA12, OMEGA13]
      OO(2,:) = [-OMEGA12, 0d0, OMEGA23]
      OO(3,:) = [-OMEGA13,-OMEGA23, 0d0]

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + Q(I,K)*BB(K,J)
        END DO
        A1(I,J) = DAUX
       END DO
      END DO

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + A1(I,K)*QT(K,J)
        END DO
        B(I,J) = DAUX
       END DO
      END DO

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + Q(I,K)*OO(K,J)
        END DO
        A1(I,J) = DAUX
       END DO
      END DO

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + A1(I,K)*QT(K,J)
        END DO
        O(I,J) = DAUX
       END DO
      END DO

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + Q(I,K)*RR(K,J)
        END DO
        A1(I,J) = DAUX
       END DO
      END DO

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + A1(I,K)*QT(K,J)
        END DO
        R(I,J) = DAUX
       END DO
      END DO

      END
C
C
C
      SUBROUTINE ConvertPsiToTau(Psi,Tau)
      IMPLICIT NONE
      REAL*8 Psi(6),Tau(6)
      REAL*8 AMAT(3,3),QMAT(3,3),QTMAT(3,3),EVCT(3)
      REAL*8 RMAT(3,3),A1(3,3),EE(3),RR(3,3),DAUX
      INTEGER I,J,K
C      
      AMAT(1,:) = [PSI(1),PSI(4),PSI(5)]
      AMAT(2,:) = [PSI(4),PSI(2),PSI(6)]
      AMAT(3,:) = [PSI(5),PSI(6),PSI(3)]
C
      CALL DSYEVD3(AMAT,QMAT,EVCT)
C
      DO I=1,3
       EE(I) = EXP(EVCT(I))
       DO J=1,3
        QTMAT(I,J) = QMAT(J,I)
       END DO
      END DO

      RR = 0d0
      RR(1,1) = EE(1)
      RR(2,2) = EE(2)
      RR(3,3) = EE(3)

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + QMAT(I,K)*RR(K,J)
        END DO
        A1(I,J) = DAUX
       END DO
      END DO

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + A1(I,K)*QTMAT(K,J)
        END DO
        RMAT(I,J) = DAUX
       END DO
      END DO

      TAU(1) = RMAT(1,1)
      TAU(2) = RMAT(2,2)
      TAU(3) = RMAT(3,3)
      TAU(4) = RMAT(1,2)
      TAU(5) = RMAT(1,3)
      TAU(6) = RMAT(2,3)

      END
C
C
C
      SUBROUTINE ConvertTauToPsi(Tau,Psi)
      IMPLICIT NONE
      REAL*8 Psi(6),Tau(6)
      REAL*8 AMAT(3,3),QMAT(3,3),QTMAT(3,3),EVCT(3)
      REAL*8 RMAT(3,3),A1(3,3),EE(3),RR(3,3),DAUX
      INTEGER I,J,K
C      
      AMAT(1,:) = [TAU(1),TAU(4),TAU(5)]
      AMAT(2,:) = [TAU(4),TAU(2),TAU(6)]
      AMAT(3,:) = [TAU(5),TAU(6),TAU(3)]
C
      CALL DSYEVD3(AMAT,QMAT,EVCT)
C
      DO I=1,3
       EE(I) = LOG(EVCT(I))
       DO J=1,3
        QTMAT(I,J) = QMAT(J,I)
       END DO
      END DO

      RR = 0d0
      RR(1,1) = EE(1)
      RR(2,2) = EE(2)
      RR(3,3) = EE(3)

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + QTMAT(I,K)*RR(K,J)
        END DO
        A1(I,J) = DAUX
       END DO
      END DO

      DO I=1,3
       DO J=1,3
        DAUX = 0d0
        DO K=1,3
        DAUX = DAUX + A1(I,K)*QMAT(K,J)
        END DO
        RMAT(I,J) = DAUX
       END DO
      END DO

      PSI(1) = RMAT(1,1)
      PSI(2) = RMAT(2,2)
      PSI(3) = RMAT(3,3)
      PSI(4) = RMAT(1,2)
      PSI(5) = RMAT(1,3)
      PSI(6) = RMAT(2,3)

      END
