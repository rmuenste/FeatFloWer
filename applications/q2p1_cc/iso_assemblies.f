
************************************************************************
      SUBROUTINE GetForceCyl_cc_iso(U1,U2,U3,Pc,Po,bALPHA,KVERT,KAREA,
     *                     KEDGE,DCORVG,DResForce,ELE,bNonNewt)
************************************************************************
*     Discrete convection operator: Q1~ elements (nonparametric)
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      USE def_cc, ONLY : Properties
      USE var_QuadScalar,ONLY : theta,itns
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      LOGICAL bNonNewt
      REAL*8  U1(*),U2(*),U3(*),Pc(*),Po(*),DCORVG(NNDIM,*)
      REAL*8  DResForce(7)
      LOGICAL bALPHA(*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM),ViscosityModel
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
      CALL ELE(0D0,0D0,0D0,-2)
C
      DResForce(1) = 0D0
      DResForce(2) = 0D0
      DResForce(3) = 0D0
      DResForce(4) = 0D0
      DResForce(5) = 0D0
      DResForce(6) = 0D0
      DResForce(7) = 0D0
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
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
C
C *** Evaluation of coordinates of the vertices
C *** Evaluation of coordinates of the vertices
      DX0 = DCORVG(1,KDFG(27))
      DY0 = DCORVG(2,KDFG(27))
      DZ0 = DCORVG(3,KDFG(27))
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
      DO JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       HBASI1 = DBAS(1,JDFL,1)
       DPP(:) = DCORVG(:,JDFG)
       XX = XX + DPP(1)*HBASI1
       YY = YY + DPP(2)*HBASI1
       ZZ = ZZ + DPP(3)*HBASI1
      END DO
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
       IF (bNonNewt) THEN
C ----=============================================---- 
       dShearSquare = DU1X**2d0 + DU2Y**2d0 + DU3Z**2d0
     *        + 0.5d0*(DU1Y+DU2X)**2d0
     *        + 0.5d0*(DU1Z+DU3X)**2d0 
     *        + 0.5d0*(DU2Z+DU3Y)**2d0

       dVisc = ViscosityModel(dShearSquare)
       ELSE
       dVisc = Properties%Viscosity(1) 
       END IF
C ----=============================================---- 

       JJ = 4*(IEL-1) + 1

	IF (theta.eq.0.5 .and. itns.gt.1) THEN
       	PressC =          Pc(JJ  ) + (XX-DX0)*Pc(JJ+1) +
     *         (YY-DY0)*Pc(JJ+2) + (ZZ-DZ0)*Pc(JJ+3)
       	PressO =          Po(JJ  ) + (XX-DX0)*Po(JJ+1) +
     *         (YY-DY0)*Po(JJ+2) + (ZZ-DZ0)*Po(JJ+3)

	Press = PressC + 0.5d0 * (PressC-PressO)
	ELSE
       	Press =          Pc(JJ  ) + (XX-DX0)*Pc(JJ+1) +
     *         (YY-DY0)*Pc(JJ+2) + (ZZ-DZ0)*Pc(JJ+3)
	END IF

C--------------------------------------------------------
c-----------Form the integrand------------------
       DN1=-DALX
       DN2=-DALY
       DN3=-DALZ
c-------------------Acting force-------------------------
c--------------Deformation calculation-------------
C
!       AH1=-Press*DN1 + dVisc*(DU1X*DN1 + DU1Y*DN2 + DU1Z*DN3)
!       AH2=-Press*DN2 + dVisc*(DU2X*DN1 + DU2Y*DN2 + DU2Z*DN3)
!       AH3=-Press*DN3 + dVisc*(DU3X*DN1 + DU3Y*DN2 + DU3Z*DN3)

       AH1=-Press*DN1+dVisc*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2 + ! full3D
     *     (DU1Z+DU3X)*DN3)
       AH4=-Press*DN2+dVisc*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2 + ! full3D
     *     (DU2Z+DU3Y)*DN3)
       AH7=-Press*DN3+dVisc*((DU3X+DU1Z)*DN1+(DU3Y+DU2Z)*DN2 + ! full3D
     *     (DU3Z+DU3Z)*DN3)
C
       DResForce(1) = DResForce(1) + AH1*OM
       DResForce(4) = DResForce(4) + AH4*OM
       DResForce(7) = DResForce(7) + AH7*OM

       AH2=dVisc*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2 + ! full3D
     *     (DU1Z+DU3X)*DN3)
       DResForce(2) = DResForce(2) + AH2*OM
       AH5=dVisc*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2 + ! full3D
     *     (DU2Z+DU3Y)*DN3)
       DResForce(5) = DResForce(5) + AH5*OM

       AH3=-Press*DN1
       DResForce(3) = DResForce(3) + AH3*OM
       AH6=-Press*DN2
       DResForce(6) = DResForce(6) + AH6*OM


C
200   CONTINUE
C
100   CONTINUE
C
999   CALL COMM_SUMMN(DResForce,7)
C
99999 CONTINUE

      END 
