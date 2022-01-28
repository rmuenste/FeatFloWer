      SUBROUTINE GetSurface(KVERT,KAREA,KEDGE,DCORVG,ELE,dArea)
      USE PP3D_MPI, ONLY:myid
      USE var_QuadScalar, ONLY : myBoundary
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      REAL*8    DCORVG(NNDIM,*)
      INTEGER   KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8    DMyOmgP(NNCUBP),DMyCubP(NNCUBP,NNAE,NNDIM)
      REAL*8    dNorm(NNDIM),GRADU(NNDIM,NNDIM),dN(3)
      INTEGER   KENTRY(NNBAS,NNBAS)
      REAL*8    ST(NNBAS,NNBAS),STT(NNBAS,NNBAS)
      REAL*8    dFluidNormal
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM),DMM(9)
      TYPE tF
       REAL*8   DHELP(NNBAS,4,NNCUBP)
      END TYPE
      TYPE(tF) F(6)
      INTEGER   NeighA(4,6),NeighAE(4,6),imap(9)
      DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
      DATA NeighAE/ 9,10,11,12,  9,14,17,13, 10,15,18,14, 
     *             11,16,19,15, 12,13,20,16, 17,18,19,20/
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
      ICUB = 2 
      CALL SetUpMyCub(DMyOmgP,DMyCubP,NCUBP,ICUB)
C
      DO IAT=1,NNAE
      DO ICUBP=1,NCUBP
       XI1=DMyCubP(ICUBP,IAT,1)
       XI2=DMyCubP(ICUBP,IAT,2)
       XI3=DMyCubP(ICUBP,IAT,3)
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
      F(IAT)%DHELP = DHELP
      END DO
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C      
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
C
      IF (IER.LT.0) GOTO 99999
C
      IF (myBoundary%iPhase(NVT+NET+NAT+IEL).eq.1) THEN
       dFluidNormal = +1d0
      ELSE
       dFluidNormal = -1d0
      END IF
C     
      DO 150 IAT=1,6
C
      ivt1 = kvert(NeighA(1,IAT),IEL)
      ivt2 = kvert(NeighA(2,IAT),IEL)
      ivt3 = kvert(NeighA(3,IAT),IEL)
      ivt4 = kvert(NeighA(4,IAT),IEL)
      
      write(*,*)'Bad routine SurfIntegral. File: QuadSc_surftens'
      stop
      if(.true.)then
c      IF (myBoundary%LS_zero(ivt1).and.myBoundary%LS_zero(ivt2).and.
c     *    myBoundary%LS_zero(ivt3).and.myBoundary%LS_zero(ivt4)) THEN
C      
      DHELP = F(IAT)%DHELP
C      
      DO 200 ICUBP=1,NCUBP
C
      XI1=DMyCubP(ICUBP,IAT,1)
      XI2=DMyCubP(ICUBP,IAT,2)
      XI3=DMyCubP(ICUBP,IAT,3)
      OM=DMyOmgP (ICUBP)
C
      DJAC=0d0
      DO JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP(JDFL,4,ICUBP)
      END DO
C
      DETJ = DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
C
      IF (IAT.eq.2.or.iat.eq.4) CALL SurfDet(DJAC,2,dA,dN)
      IF (IAT.eq.1.or.iat.eq.6) CALL SurfDet(DJAC,3,dA,dN)
      IF (IAT.eq.3.or.iat.eq.5) CALL SurfDet(DJAC,1,dA,dN)
C
      CALL ELE(XI1,XI2,XI3,0)
      IF (IER.LT.0) GOTO 99999
C
      dArea = dArea + 0.5d0*dA*OM
C
200   CONTINUE
C
      END IF
150   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE GetVolume(U1,U2,U3,KVERT,KAREA,KEDGE,DCORVG,ELE,
     *           dVol,dCenter,dRise)
************************************************************************
*     Discrete convection operator: Q1~ elements (nonparametric)
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID
      USE var_QuadScalar, ONLY : myBoundary
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      REAL*8  dCenter(3),dVol,dRise(3)
      REAL*8  DCORVG(NNDIM,*),U1(*),U2(*),U3(*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
      REAL*8  DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
      CALL ELE(0D0,0D0,0D0,-2)
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
C      
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      IF (myBoundary%iPhase(NVT+NET+NAT+IEL) .eq. 1) THEN
       dIndicator = 1d0
      ELSE
       dIndicator = 0d0
      END IF
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
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
      DJAC=0d0
      DO JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP(JDFL,4,ICUBP)
      END DO
C
C *** Jacobian of the bilinear mapping onto the reference element
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
       XX = 0d0; YY = 0d0; ZZ = 0d0 
       UX = 0d0; UY = 0d0; UZ = 0d0
C       
       DO I=1,IDFL
        IL = KDFL(I)
        IG = KDFG(I)
        HBAS=DBAS(1,IL,1)
        XX = XX + DCORVG(1,IG)*HBAS
        YY = YY + DCORVG(2,IG)*HBAS
        ZZ = ZZ + DCORVG(3,IG)*HBAS
        UX = UX + U1(IG)*HBAS
        UY = UY + U2(IG)*HBAS
        UZ = UZ + U3(IG)*HBAS
       END DO
C
       dVol = dVol + dIndicator*OM
       dCenter(1) = dCenter(1) + dIndicator*OM*XX
       dCenter(2) = dCenter(2) + dIndicator*OM*YY
       dCenter(3) = dCenter(3) + dIndicator*OM*ZZ
       dRise(1) = dRise(1)     + dIndicator*OM*UX
       dRise(2) = dRise(2)     + dIndicator*OM*UY
       dRise(3) = dRise(3)     + dIndicator*OM*UZ
C
200   CONTINUE
C
100   CONTINUE
C
99999 CONTINUE

      END
!      SUBROUTINE AddSurfaceTension()
!
!      USE def_QuadScalar
!      USE Transport_Q2P1, ONLY : QuadSc,ST_force,Sigma
!      USE PLinScalar, ONLY : dNorm
!      IMPLICIT NONE
!      REAL*8 :: DEPS=0.015d0
!      EXTERNAL E013
!
!      ILEV=NLMAX
!      CALL SETLEV(2)
!!       PAUSE
!      CALL Q1toQ2SurfTens(QuadSc%defU,QuadSc%defV,QuadSc%defW,dNorm,
!     *     ST_force,KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
!     *     KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),TSTEP,Sigma,E013)
!       RETURN
!
!!       CALL SurfTens(QuadSc%defU,QuadSc%defV,QuadSc%defW,normal,
!!      *     LevSet%val(NLMAX)%x,KWORK(L(LVERT)),KWORK(L(LAREA)),
!!      *     KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),
!!      *     VWORK(L(KLVOL(NLMAX))),Sigma,TSTEP,DEPS,E013)
!!       CALL SurfTens_old(QuadSc%defU,QuadSc%defV,QuadSc%defW,normal,
!!      *     LevSet%val(NLMAX)%x,KWORK(L(LVERT)),KWORK(L(LAREA)),
!!      *     KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),
!!      *     Sigma,TSTEP,DEPS,E013)
!!       QuadSc%valU(NLMAX)%x = 0d0
!!       QuadSc%valV(NLMAX)%x = 0d0
!!       QuadSc%valW(NLMAX)%x = 0d0
!!       CALL SurfTens(QuadSc%valU(NLMAX)%x,QuadSc%valV(NLMAX)%x,
!!      *     QuadSc%valW(NLMAX)%x,
!!      *     normal,LevSet%val(NLMAX)%x,KWORK(L(LVERT)),KWORK(L(LAREA)),
!!      *     KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),
!!      *     VWORK(L(KLVOL(NLMAX))),Sigma,TSTEP,DEPS,E013)
!
!      END


************************************************************************
      SUBROUTINE Q1toQ2SurfTens(D1,D2,D3,DNORM,DF,KVERT,KAREA,KEDGE,
     *           KINT,DCORVG,TSTEP,DSIGMA,ELE)
************************************************************************
*     Discrete convection operator: Q1~ elements (nonparametric)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DNORM(3,*),DF(*),D1(*),D2(*),D3(*)
      DIMENSION KINT(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KBFG(NNBAS),KBFL(NNBAS)
      DIMENSION DPHI(NNVE,4)
      REAL*8    TSTEP
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
!       return
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
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      ILINT=KINT(IEL)
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
C--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
!       NJALFA=0
!       NIALFA=0
!       DO I=1,IDFL
!         IG=KDFG(I)
!         IF (DLS(IG).LE.0d0) THEN
!          NJALFA=NJALFA+1
!         ELSE
!          NIALFA=NIALFA+1
!         ENDIF
!       ENDDO
!       IF(NJALFA.EQ.8.OR.NIALFA.EQ.8) GOTO 100
C--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
      IF (IER.LT.0) GOTO 99999
C
      CALL NDFGL(IEL,1,11,KVERT,KEDGE,KAREA,KBFG,KBFL)
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
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      CALL LinElem(XI1,XI2,XI3,DPHI,DETJ,DJAC,.TRUE.)
      IF (IER.LT.0) GOTO 99999
C
      DFORCEx =0D0
      DFORCEy =0D0
      DFORCEz =0D0
C
      DO I=1,NNVE
        IG = KBFG(I)
        IL = KBFL(I)
        DFORCEx = DFORCEx + DSIGMA*DF(IG)*DNORM(1,IG)*DPHI(IL,1)
        DFORCEy = DFORCEy + DSIGMA*DF(IG)*DNORM(2,IG)*DPHI(IL,1)
        DFORCEz = DFORCEz + DSIGMA*DF(IG)*DNORM(3,IG)*DPHI(IL,1)
      ENDDO
C
      DO I=1,IDFL
        IG = KDFG(I)
        IL = KDFL(I)
        D1(IG) = D1(IG) - DFORCEx*DBAS(1,IL,1)*OM*TSTEP
        D2(IG) = D2(IG) - DFORCEy*DBAS(1,IL,1)*OM*TSTEP
        D3(IG) = D3(IG) - DFORCEz*DBAS(1,IL,1)*OM*TSTEP
      ENDDO
C
200   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
      SUBROUTINE LinElem(X1,X2,X3,DBAS,DETJ,DJAC,BDER)
      IMPLICIT NONE
C
      LOGICAL BDER
      REAL*8 Q8,XJ1,X1,X2,X3
      PARAMETER (Q8=0.125D0)
      INTEGER I
      REAL*8 DPHI(8,3),DJAC(3,3),DETJ,DBAS(8,4)
C
      DBAS(1,1)=Q8*(1D0-X1)*(1D0-X2)*(1D0-X3)
      DBAS(2,1)=Q8*(1D0+X1)*(1D0-X2)*(1D0-X3)
      DBAS(3,1)=Q8*(1D0+X1)*(1D0+X2)*(1D0-X3)
      DBAS(4,1)=Q8*(1D0-X1)*(1D0+X2)*(1D0-X3)
      DBAS(5,1)=Q8*(1D0-X1)*(1D0-X2)*(1D0+X3)
      DBAS(6,1)=Q8*(1D0+X1)*(1D0-X2)*(1D0+X3)
      DBAS(7,1)=Q8*(1D0+X1)*(1D0+X2)*(1D0+X3)
      DBAS(8,1)=Q8*(1D0-X1)*(1D0+X2)*(1D0+X3)
C
      IF (BDER) THEN
      DPHI(1,1)=-Q8*(1D0-X2)*(1D0-X3)
      DPHI(2,1)= Q8*(1D0-X2)*(1D0-X3)
      DPHI(3,1)= Q8*(1D0+X2)*(1D0-X3)
      DPHI(4,1)=-Q8*(1D0+X2)*(1D0-X3)
      DPHI(5,1)=-Q8*(1D0-X2)*(1D0+X3)
      DPHI(6,1)= Q8*(1D0-X2)*(1D0+X3)
      DPHI(7,1)= Q8*(1D0+X2)*(1D0+X3)
      DPHI(8,1)=-Q8*(1D0+X2)*(1D0+X3)
C
      DPHI(1,2)=-Q8*(1D0-X1)*(1D0-X3)
      DPHI(2,2)=-Q8*(1D0+X1)*(1D0-X3)
      DPHI(3,2)= Q8*(1D0+X1)*(1D0-X3)
      DPHI(4,2)= Q8*(1D0-X1)*(1D0-X3)
      DPHI(5,2)=-Q8*(1D0-X1)*(1D0+X3)
      DPHI(6,2)=-Q8*(1D0+X1)*(1D0+X3)
      DPHI(7,2)= Q8*(1D0+X1)*(1D0+X3)
      DPHI(8,2)= Q8*(1D0-X1)*(1D0+X3)
C
      DPHI(1,3)=-Q8*(1D0-X1)*(1D0-X2)
      DPHI(2,3)=-Q8*(1D0+X1)*(1D0-X2)
      DPHI(3,3)=-Q8*(1D0+X1)*(1D0+X2)
      DPHI(4,3)=-Q8*(1D0-X1)*(1D0+X2)
      DPHI(5,3)= Q8*(1D0-X1)*(1D0-X2)
      DPHI(6,3)= Q8*(1D0+X1)*(1D0-X2)
      DPHI(7,3)= Q8*(1D0+X1)*(1D0+X2)
      DPHI(8,3)= Q8*(1D0-X1)*(1D0+X2)
C
      XJ1=1D0/DETJ
C
      DO I=1,8
       DBAS(I,2)= XJ1*(
     *   DPHI(I,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *  -DPHI(I,2)*(DJAC(2,1)*DJAC(3,3)-DJAC(3,1)*DJAC(2,3))
     *  +DPHI(I,3)*(DJAC(2,1)*DJAC(3,2)-DJAC(3,1)*DJAC(2,2)))
      END DO
      DO I=1,8
       DBAS(I,3)= XJ1*(
     *   -DPHI(I,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *   +DPHI(I,2)*(DJAC(1,1)*DJAC(3,3)-DJAC(3,1)*DJAC(1,3))
     *   -DPHI(I,3)*(DJAC(1,1)*DJAC(3,2)-DJAC(3,1)*DJAC(1,2)))
      END DO
      DO I=1,8
       DBAS(I,4)= XJ1*(
     *    DPHI(I,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
     *   -DPHI(I,2)*(DJAC(1,1)*DJAC(2,3)-DJAC(2,1)*DJAC(1,3))
     *   +DPHI(I,3)*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2)))
      END DO
C
      END IF
C
99999 END

C
C
C
      SUBROUTINE SurfDet(D,KK,DD,dN)
      REAL*8 D(3,3),DD,dN(3)
      INTEGER kk
      REAL*8 F(2,2)
      INTEGER i,j,ii,jj

      daux = 0d0
      if (kk.eq.1) THEN
       ii = 2; jj = 3
      end if
      if (kk.eq.2) THEN
       ii = 1; jj = 3
      end if
      if (kk.eq.3) THEN
       ii = 1; jj = 2 
      end if

      dN(1) = +D(2,ii)*D(3,jj) - D(2,jj)*D(3,ii)
      dN(2) = -D(1,ii)*D(3,jj) + D(1,jj)*D(3,ii)
      dN(3) = +D(1,ii)*D(2,jj) - D(1,jj)*D(2,ii)
      
      DD = SQRT(dN(1)**2d0 + dN(2)**2d0 + dN(3)**2d0)

      dN = dN/DD
      
      DD = 4d0*DD

      END 
C
C
C
      SUBROUTINE GET_area(P1,P2,P3,P4,dA)
      REAL*8 P1(3),P2(3),P3(3),P4(3),dA
      REAL*8 AX2,AY2,AZ2,AX3,AY3,AZ3,AX4,AY4,AZ4,AX,AY,AZ
      REAL*8 DNX,DNY,DNZ,DNAR1,DNAR2
      
      AX2=P2(1)-P1(1)
      AY2=P2(2)-P1(2)
      AZ2=P2(3)-P1(3)
      AX3=P3(1)-P1(1)
      AY3=P3(2)-P1(2)
      AZ3=P3(3)-P1(3)
      AX4=P4(1)-P1(1)
      AY4=P4(2)-P1(2)
      AZ4=P4(3)-P1(3)
      
      DNX=(AY3*AZ2)-(AZ3*AY2)
      DNY=(AZ3*AX2)-(AX3*AZ2)
      DNZ=(AX3*AY2)-(AY3*AX2)
      DNAR1=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
      !write(*,*) DNX,DNY,DNZ
      
      DNX=(AY4*AZ3)-(AZ4*AY3)
      DNY=(AZ4*AX3)-(AX4*AZ3)
      DNZ=(AX4*AY3)-(AY4*AX3)
      DNAR2=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
      
      dA=0.5D0*(DNAR1+DNAR2)

      END SUBROUTINE GET_area
C
C
C
      SUBROUTINE GET_NormalArea(P1,P2,P3,P4,P5,dn)
      REAL*8 P1(3),P2(3),P3(3),P4(3),P5(3),dn(3)
      REAL*8 AX2,AY2,AZ2,AX3,AY3,AZ3,AX4,AY4,AZ4,AX,AY,AZ
      REAL*8 PA(3),NA(3)
      REAL*8 DNX,DNY,DNZ,DNAR1,DNAR2,DFAC,DHN
      
      AX2=P2(1)-P1(1)
      AY2=P2(2)-P1(2)
      AZ2=P2(3)-P1(3)
      AX3=P3(1)-P1(1)
      AY3=P3(2)-P1(2)
      AZ3=P3(3)-P1(3)
      AX4=P4(1)-P1(1)
      AY4=P4(2)-P1(2)
      AZ4=P4(3)-P1(3)
      
      PA(:)=(P1(:)+P2(:)+P3(:)+P4(:))*0.25D0
      NA(:) =P5(:) - PA(:)

      DN(1)=(AY3*AZ2)-(AZ3*AY2)
      DN(2)=(AZ3*AX2)-(AX3*AZ2)
      DN(3)=(AX3*AY2)-(AY3*AX2)
      DNAR1=SQRT(DN(1)*DN(1)+DN(2)*DN(2)+DN(3)*DN(3))
      
      DN(1)=(AY4*AZ3)-(AZ4*AY3)
      DN(2)=(AZ4*AX3)-(AX4*AZ3)
      DN(3)=(AX4*AY3)-(AY4*AX3)
      DNAR2=SQRT(DN(1)*DN(1)+DN(2)*DN(2)+DN(3)*DN(3))
      
      DFAC=0.5D0*(DNAR1+DNAR2)/DNAR2
      DN(:) =DFAC*DN(:)
C
      DHN = DN(1)*NA(1) + DN(2)*NA(2) + DN(3)*NA(3)
      IF (DHN.LT.0D0) THEN
       DN(:) = -DN(:)
      ENDIF
      
      END SUBROUTINE GET_NormalArea
C
C
C
      SUBROUTINE RecoverNormals(NU,NV,NW,DD,KVERT,
     *           KAREA,KEDGE,DCORVG,ELE)
      USE PP3D_MPI, ONLY:myid
      USE var_QuadScalar, ONLY : myBoundary
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C

      REAL*8 NU(*),NV(*),NW(*)
      REAL*8 DD(*), DCORVG(NNDIM,*)
      INTEGER  KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8 DMyOmgP(NNCUBP),DMyCubP(NNCUBP,NNAE,NNDIM)
      REAL*8    TSTEP,dN(3)
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM),dNorm(3)
      TYPE tF
       REAL*8   DHELP(NNBAS,4,NNCUBP)
      END TYPE
      TYPE(tF) F(6)
      REAL*8    dFluidNormal,PointC(3),PointA(3),dLine(3),dVELO(3)
      INTEGER NeighA(4,6)
      DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

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
      ICUB = 3 
      CALL SetUpMyCub(DMyOmgP,DMyCubP,NCUBP,ICUB)
C
      DO IAT=1,NNAE
      DO ICUBP=1,NCUBP
       XI1=DMyCubP(ICUBP,IAT,1)
       XI2=DMyCubP(ICUBP,IAT,2)
       XI3=DMyCubP(ICUBP,IAT,3)
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
      F(IAT)%DHELP = DHELP
      END DO
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C      
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
C
      IF (IER.LT.0) GOTO 99999
C
      IF (myBoundary%iPhase(NVT+NET+NAT+IEL).ne.0) THEN
       dFluidNormal = +1d0
      ELSE
       dFluidNormal = 0d0
      END IF
C     
      PointC(:) = DCORVG(:,NVT+NET+NAT+IEL)
C      
      DO 150 IAT=1,6
C
      ivt1 = kvert(NeighA(1,IAT),IEL)
      ivt2 = kvert(NeighA(2,IAT),IEL)
      ivt3 = kvert(NeighA(3,IAT),IEL)
      ivt4 = kvert(NeighA(4,IAT),IEL)
      
      IF (myBoundary%LS_zero(ivt1).ne.0.and.
     *    myBoundary%LS_zero(ivt2).ne.0.and.
     *    myBoundary%LS_zero(ivt3).ne.0.and.
     *    myBoundary%LS_zero(ivt4).ne.0) THEN
C     
      PointA(:) = DCORVG(:,NVT+NET+KAREA(IAT,IEL))
      dLine(:)  = PointC(:) - PointA(:)
C      
      DHELP = F(IAT)%DHELP
C      
      DO 200 ICUBP=1,NCUBP
C
      XI1=DMyCubP(ICUBP,IAT,1)
      XI2=DMyCubP(ICUBP,IAT,2)
      XI3=DMyCubP(ICUBP,IAT,3)
      OM=DMyOmgP (ICUBP)
C
      DJAC=0d0
      DO JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP(JDFL,4,ICUBP)
      END DO
C
      DETJ = DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
C
      IF (IAT.eq.2.or.iat.eq.4) CALL SurfDet(DJAC,2,dA,dN)
      IF (IAT.eq.1.or.iat.eq.6) CALL SurfDet(DJAC,3,dA,dN)
      IF (IAT.eq.3.or.iat.eq.5) CALL SurfDet(DJAC,1,dA,dN)
C
      CALL ELE(XI1,XI2,XI3,0)
      IF (IER.LT.0) GOTO 99999
C
      dNorm = dFluidNormal*dN
      dDotProd = dNorm(1)*dLine(1)+dNorm(2)*dLine(2)+dNorm(3)*dLine(3)
      IF (dDotProd.lt.0d0) dNorm = -dNorm
C
      DO I=1,IDFL
        IG = KDFG(I)
        IL = KDFL(I)
        HBAS = DBAS(1,IL,1)
        NU(IG) = NU(IG) - dA*OM*dNorm(1)*HBAS
        NV(IG) = NV(IG) - dA*OM*dNorm(2)*HBAS
        NW(IG) = NW(IG) - dA*OM*dNorm(3)*HBAS
        DD(IG) = DD(IG) + dA*OM*HBAS
      ENDDO

200   CONTINUE
      END IF
150   CONTINUE
100   CONTINUE
99999 END
