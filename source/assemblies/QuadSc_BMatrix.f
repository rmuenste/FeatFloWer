 
***************************************************
      SUBROUTINE Build_BTMatP1(DAx,DAy,DAz,KLD,KCOL,KVERT,KAREA,
     *           KEDGE,DCORVG,NA,ELE)
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
      DIMENSION DAx(*),DAy(*),DAz(*),DBAS1(4)
      DIMENSION DCORVG(NNDIM,*),KLD(*),KCOL(*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG1(NNBAS),KDFL1(NNBAS)
      DIMENSION KDFG2(NNBAS),KDFL2(NNBAS)
      DIMENSION KDFG3(NNBAS),KDFL3(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRYx(NNBAS,NNBAS)
      DIMENSION DENTRYy(NNBAS,NNBAS),DENTRYz(NNBAS,NNBAS)
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
      ICUB = 9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E013A(XI1,XI2,XI3,DHELP_Q2,ICUBP)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
      END DO
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
      DX0I = DCORVG(1,KDFG2(27))
      DY0I = DCORVG(2,KDFG2(27))
      DZ0I = DCORVG(3,KDFG2(27))
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
       JDFL=KDFL2(JDOFE)
       JDFG=KDFG2(JDOFE)
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
       JDFL=KDFL2(JDOFE)
       JDFG=KDFG2(JDOFE)
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
! ----------------------------------------------------------------
!     Computation of the cartesian coordiante of the cubature point
      XX = 0d0; YY = 0d0; ZZ = 0d0
      DO JDOFE=1,IDFL2
       JDFL=KDFL2(JDOFE)
       JDFG=KDFG2(JDOFE)
       HBASI1 = DBAS(1,JDFL,1)
       DPP(:) = DCORVG(:,JDFG)
       XX = XX + DPP(1)*HBASI1
       YY = YY + DPP(2)*HBASI1
       ZZ = ZZ + DPP(3)*HBASI1
      END DO
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
      SUBROUTINE Build_BMatP1(DAx,DAy,DAz,KLD,KCOL,KVERT,KAREA,
     *           KEDGE,DCORVG,NA,ELE)
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
      DIMENSION DAx(*),DAy(*),DAz(*),DBAS1(4)
      DIMENSION DCORVG(NNDIM,*),KLD(*),KCOL(*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG1(NNBAS),KDFL1(NNBAS)
      DIMENSION KDFG2(NNBAS),KDFL2(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRYx(NNBAS,NNBAS)
      DIMENSION DENTRYy(NNBAS,NNBAS),DENTRYz(NNBAS,NNBAS)
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
      ICUB = 9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E013A(XI1,XI2,XI3,DHELP_Q2,ICUBP)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
      END DO
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
C *** Evaluation of coordinates of the vertices
      DX0I = DCORVG(1,KDFG1(27))
      DY0I = DCORVG(2,KDFG1(27))
      DZ0I = DCORVG(3,KDFG1(27))
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
       JDFL=KDFL1(JDOFE)
       JDFG=KDFG1(JDOFE)
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
       JDFL=KDFL1(JDOFE)
       JDFG=KDFG1(JDOFE)
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
! ----------------------------------------------------------------
!     Computation of the cartesian coordiante of the cubature point
      XX = 0d0; YY = 0d0; ZZ = 0d0
      DO JDOFE=1,IDFL1
       JDFL=KDFL1(JDOFE)
       JDFG=KDFG1(JDOFE)
       HBASI1 = DBAS(1,JDFL,1)
       DPP(:) = DCORVG(:,JDFG)
       XX = XX + DPP(1)*HBASI1
       YY = YY + DPP(2)*HBASI1
       ZZ = ZZ + DPP(3)*HBASI1
      END DO
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

99999 END
C
C
C
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
