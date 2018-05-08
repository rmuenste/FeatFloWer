      FUNCTION Breyer_Carreau(NormShearSquare)
      USE Transport_Q2P1, ONLY : Properties
      IMPLICIT NONE

      real*8 :: Breyer_Carreau
      real*8, intent (in) :: NormShearSquare

      REAL*8 :: dA,dB,dC,dStrs

      dA = 3499.76071*1d1
      dB = 0.08992
      dC = 0.77067
      
       dStrs = DSQRT(NormShearSquare)
       Breyer_Carreau = dA/(1d0+(dStrs/dB)**dC)
      RETURN
      END

      FUNCTION PolyFLOW_Carreau(NormShearSquare)
      USE Transport_Q2P1, ONLY : Properties
      IMPLICIT NONE

      real*8 :: PolyFLOW_Carreau
      real*8, intent (in) :: NormShearSquare

      REAL*8 :: dN

      dN = Properties%PowerLawExp-1d0
      PolyFLOW_Carreau = 
     *Properties%Viscosity(1)*(1d-4 + NormShearSquare)**dN

!      WRITE(*,*) Properties%Viscosity(1), NormShearSquare,dN
      RETURN
      END

!       FUNCTION PolyFLOW_Carreau(NormShearSquare)
!       USE QuadScalar, ONLY : Properties,timens
!       IMPLICIT NONE
! 
!       real*8 :: PolyFLOW_Carreau
!       real*8, intent (in) :: NormShearSquare
!       real*8 dMinLim,dMaxLim,dAlpha,dmu_0,dStrs_0,dStrs,dModel
! 
!       REAL*8 :: dN
! 
!       dMinLim = MAX( 0.1d0,1.6d0 - timens*1d1)
!       dMaxLim = MIN( 3.1d0,1.6d0 + timens*1d1)
! 
!       dStrs = log(NormShearSquare**0.5d0)
! 
!       dModel = 1d0*(dStrs)
! !      dModel = dmu_0/(1d0+(dStrs/dStrs_0)**dAlpha)
! 
!       PolyFLOW_Carreau = MAX(MIN(dModel,dMaxLim),dMinLim)
! 
!       RETURN
!       END 

      FUNCTION HogenPowerlaw(NormShearSquare)
      USE Transport_Q2P1, ONLY : Properties
      IMPLICIT NONE

      real*8 :: HogenPowerlaw
      real*8, intent (in) :: NormShearSquare

      REAL*8 :: dN

      dN = Properties%PowerLawExp-1d0
      HogenPowerlaw = 
     *Properties%Viscosity(1)*(1d-4 + NormShearSquare)**dN

      RETURN
      END 
C
C
C
C$$$************************************************************************
C$$$      SUBROUTINE CUBATURESTRESS(U1,U2,U3,
C$$$     *           DS11,DS22,DS33,DS12,DS13,DS23,DS21,DS31,DS32,
C$$$     *           NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,DCORVG,ELE)
C$$$************************************************************************
C$$$*
C$$$*-----------------------------------------------------------------------
C$$$      USE PP3D_MPI, ONLY:myid

C$$$      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C$$$      CHARACTER SUB*6,FMT*15,CPARAM*120
C$$$C
C$$$      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
C$$$     *           NNDIM=3,NNCOF=10)
C$$$      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C$$$C
C$$$      REAL*8    DS11(*),DS22(*),DS33(*)
C$$$      REAL*8    DS12(*),DS13(*),DS23(*),DS21(*),DS31(*),DS32(*)
C$$$      REAL*8    U1(*),U2(*),U3(*), dETA0,dExp, dVisc, dShear
C$$$C
C$$$      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
C$$$      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C$$$      DIMENSION KENTRY(NNBAS,NNBAS)
C$$$      REAL*8    S11(NNBAS,NNBAS),S22(NNBAS,NNBAS),S33(NNBAS,NNBAS)
C$$$      REAL*8    S12(NNBAS,NNBAS),S13(NNBAS,NNBAS),S23(NNBAS,NNBAS)
C$$$      REAL*8    S21(NNBAS,NNBAS),S31(NNBAS,NNBAS),S32(NNBAS,NNBAS)
C$$$      REAL*8    PolyFLOW_Carreau
C$$$C
C$$$      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C$$$      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
C$$$      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
C$$$      DIMENSION DU3(NNBAS), GRADU3(NNDIM)
C$$$C
C$$$      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C$$$      COMMON /ERRCTL/ IER,ICHECK
C$$$      COMMON /CHAR/   SUB,FMT(3),CPARAM
C$$$      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
C$$$     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
C$$$     *                IEL,NDIM
C$$$      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
C$$$     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
C$$$      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
C$$$      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS

C$$$      COMMON /COAUX1/ KDFG,KDFL,IDFL
C$$$C
C$$$C *** user COMMON blocks
C$$$      INTEGER  VIPARM 
C$$$      DIMENSION VIPARM(100)
C$$$      EQUIVALENCE (IAUSAV,VIPARM)
C$$$      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
C$$$     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
C$$$     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
C$$$     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
C$$$     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
C$$$     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C$$$C
C$$$      SAVE
C$$$C
C$$$      DO 1 I= 1,NNDER
C$$$1     BDER(I)=.FALSE.
C$$$C
C$$$      DO 2 I=1,4
C$$$2     BDER(I)=.TRUE.
C$$$C
C$$$      IELTYP=-1
C$$$      CALL ELE(0D0,0D0,0D0,IELTYP)
C$$$      IDFL=NDFL(IELTYP)
C$$$C
C$$$      ICUB=9
C$$$      CALL CB3H(ICUB)
C$$$      IF (IER.NE.0) GOTO 99999
C$$$C
C$$$************************************************************************
C$$$C *** Calculation of the matrix - storage technique 7 or 8
C$$$************************************************************************
C$$$      ICUBP=ICUB
C$$$      CALL ELE(0D0,0D0,0D0,-2)
C$$$C
C$$$C *** Loop over all elements
C$$$      DO 100 IEL=1,NEL
C$$$C
C$$$C
C$$$      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
C$$$      IF (IER.LT.0) GOTO 99999
C$$$C
C$$$C *** Determine entry positions in matrix
C$$$      DO 110 JDOFE=1,IDFL
C$$$      ILD=KLDA(KDFG(JDOFE))
C$$$      KENTRY(JDOFE,JDOFE)=ILD
C$$$      S11(JDOFE,JDOFE)=0D0
C$$$      S22(JDOFE,JDOFE)=0D0
C$$$      S33(JDOFE,JDOFE)=0D0
C$$$      S12(JDOFE,JDOFE)=0D0
C$$$      S13(JDOFE,JDOFE)=0D0
C$$$      S23(JDOFE,JDOFE)=0D0
C$$$      S21(JDOFE,JDOFE)=0D0
C$$$      S31(JDOFE,JDOFE)=0D0
C$$$      S32(JDOFE,JDOFE)=0D0
C$$$      JCOL0=ILD
C$$$      DO 111 IDOFE=1,IDFL
C$$$      IF (IDOFE.EQ.JDOFE) GOTO 111
C$$$      IDFG=KDFG(IDOFE)
C$$$      DO 112 JCOL=JCOL0,NA
C$$$      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
C$$$112   CONTINUE
C$$$113   JCOL0=JCOL+1
C$$$      KENTRY(JDOFE,IDOFE)=JCOL
C$$$      S11(JDOFE,IDOFE)=0D0
C$$$      S22(JDOFE,IDOFE)=0D0
C$$$      S33(JDOFE,IDOFE)=0D0
C$$$      S12(JDOFE,IDOFE)=0D0
C$$$      S13(JDOFE,IDOFE)=0D0
C$$$      S23(JDOFE,IDOFE)=0D0
C$$$      S21(JDOFE,IDOFE)=0D0
C$$$      S31(JDOFE,IDOFE)=0D0
C$$$      S32(JDOFE,IDOFE)=0D0
C$$$111   CONTINUE
C$$$110   CONTINUE
C$$$C
C$$$C *** Evaluation of coordinates of the vertices
C$$$!      if (myid.eq.1) then
C$$$!      DO idofe=1,idfl
C$$$!       write(*,'(50I10)') myid,iel, idofe, kdfl(idofe),kdfg(idofe)
C$$$!      end do
C$$$!      end if
C$$$!      pause

C$$$      DO 120 IVE=1,NVE
C$$$      JP=KVERT(IVE,IEL)
C$$$      KVE(IVE)=JP
C$$$      DX(IVE)=DCORVG(1,JP)
C$$$      DY(IVE)=DCORVG(2,JP)
C$$$      DZ(IVE)=DCORVG(3,JP)
C$$$120   CONTINUE
C$$$C
C$$$      DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
C$$$      DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
C$$$      DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
C$$$      DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
C$$$      DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
C$$$      DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
C$$$      DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
C$$$      DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
C$$$      DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
C$$$      DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
C$$$      DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
C$$$      DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
C$$$      DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
C$$$      DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
C$$$      DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
C$$$      DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
C$$$      DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
C$$$      DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
C$$$      DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
C$$$      DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
C$$$      DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
C$$$      DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
C$$$      DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
C$$$      DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
C$$$C
C$$$!---============================---
C$$$      DO 130 JDOFE=1,IDFL
C$$$      JDFG=KDFG(JDOFE)
C$$$      JDFL=KDFL(JDOFE)
C$$$!!!   local = global      
C$$$      DU1(JDFL) = U1(JDFG) 
C$$$      DU2(JDFL) = U2(JDFG)
C$$$      DU3(JDFL) = U3(JDFG)

C$$$ 130  CONTINUE      
 
C$$$     dVisc=0.0d0
C$$$! ---===========================---
C$$$C *** Loop1 over all cubature points
C$$$      DO 199 ICUBP=1,NCUBP
C$$$C
C$$$      XI1=DXI(ICUBP,1)
C$$$      XI2=DXI(ICUBP,2)
C$$$      XI3=DXI(ICUBP,3)
C$$$C
C$$$C *** Jacobian of the bilinear mapping onto the reference element
C$$$      DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
C$$$      DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
C$$$      DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
C$$$      DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
C$$$      DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
C$$$      DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
C$$$      DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
C$$$      DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
C$$$      DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
C$$$      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
C$$$     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
C$$$     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
C$$$      OM=DOMEGA(ICUBP)*ABS(DETJ)
C$$$C
C$$$      CALL ELE(XI1,XI2,XI3,-3)
C$$$      IF (IER.LT.0) GOTO 99999
C$$$C 
C$$$C ---=========================---
C$$$      GRADU1(1)=0D0!U
C$$$      GRADU1(2)=0D0
C$$$      GRADU1(3)=0D0

C$$$      GRADU2(1)=0D0!V
C$$$      GRADU2(2)=0D0
C$$$      GRADU2(3)=0D0

C$$$      GRADU3(1)=0D0!W
C$$$      GRADU3(2)=0D0
C$$$      GRADU3(3)=0D0 
C$$$C
C$$$      DO 219 JDOFE=1,IDFL
C$$$       JDFL=KDFL(JDOFE)! local number of basic function
       
C$$$       GRADU1(1)=GRADU1(1) + DU1(JDFL)*DBAS(1,JDFL,2)!DUX
C$$$       GRADU1(2)=GRADU1(2) + DU1(JDFL)*DBAS(1,JDFL,3)!DUY
C$$$       GRADU1(3)=GRADU1(3) + DU1(JDFL)*DBAS(1,JDFL,4)!DUZ

C$$$       GRADU2(1)=GRADU2(1) + DU2(JDFL)*DBAS(1,JDFL,2)!DVX
C$$$       GRADU2(2)=GRADU2(2) + DU2(JDFL)*DBAS(1,JDFL,3)!DVY
C$$$       GRADU2(3)=GRADU2(3) + DU2(JDFL)*DBAS(1,JDFL,4)!DVZ

C$$$       GRADU3(1)=GRADU3(1) + DU3(JDFL)*DBAS(1,JDFL,2)!DWX
C$$$       GRADU3(2)=GRADU3(2) + DU3(JDFL)*DBAS(1,JDFL,3)!DWY
C$$$       GRADU3(3)=GRADU3(3) + DU3(JDFL)*DBAS(1,JDFL,4)!DWZ

C$$$ 219  CONTINUE

C$$$C ----=============================================---- 
C$$$       dShearSquare = GRADU1(1)**2d0 + GRADU2(2)**2d0 
C$$$     *        + GRADU3(3)**2d0 + 0.5d0*(GRADU1(2)+GRADU2(1))**2d0
C$$$     *        + 0.5d0*(GRADU1(3)+GRADU3(1))**2d0 
C$$$     *        + 0.5d0*(GRADU2(3)+GRADU3(2))**2d0

C$$$       dVisc = dVisc + PolyFLOW_Carreau(dShearSquare)

C$$$C ----=============================================---- 
C$$$C
      
C$$$ 199  CONTINUE
C$$$      dVisc = dVisc/real(NCUBP)

C$$$! ---===========================---
C$$$C *** Loop2 over all cubature points
C$$$      DO 200 ICUBP=1,NCUBP
C$$$C
C$$$      XI1=DXI(ICUBP,1)
C$$$      XI2=DXI(ICUBP,2)
C$$$      XI3=DXI(ICUBP,3)
C$$$C
C$$$C *** Jacobian of the bilinear mapping onto the reference element
C$$$      DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
C$$$      DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
C$$$      DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
C$$$      DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
C$$$      DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
C$$$      DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
C$$$      DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
C$$$      DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
C$$$      DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
C$$$      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
C$$$     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
C$$$     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
C$$$      OM=DOMEGA(ICUBP)*ABS(DETJ)
C$$$C
C$$$      CALL ELE(XI1,XI2,XI3,-3)
C$$$      IF (IER.LT.0) GOTO 99999
C$$$C 
C$$$C ---=========================---
C$$$      GRADU1(1)=0D0!U
C$$$      GRADU1(2)=0D0
C$$$      GRADU1(3)=0D0

C$$$      GRADU2(1)=0D0!V
C$$$      GRADU2(2)=0D0
C$$$      GRADU2(3)=0D0

C$$$      GRADU3(1)=0D0!W
C$$$      GRADU3(2)=0D0
C$$$      GRADU3(3)=0D0 
C$$$C
C$$$      DO 220 JDOFE=1,IDFL
C$$$       JDFL=KDFL(JDOFE)! local number of basic function
       
C$$$       GRADU1(1)=GRADU1(1) + DU1(JDFL)*DBAS(1,JDFL,2)!DUX
C$$$       GRADU1(2)=GRADU1(2) + DU1(JDFL)*DBAS(1,JDFL,3)!DUY
C$$$       GRADU1(3)=GRADU1(3) + DU1(JDFL)*DBAS(1,JDFL,4)!DUZ

C$$$       GRADU2(1)=GRADU2(1) + DU2(JDFL)*DBAS(1,JDFL,2)!DVX
C$$$       GRADU2(2)=GRADU2(2) + DU2(JDFL)*DBAS(1,JDFL,3)!DVY
C$$$       GRADU2(3)=GRADU2(3) + DU2(JDFL)*DBAS(1,JDFL,4)!DVZ

C$$$       GRADU3(1)=GRADU3(1) + DU3(JDFL)*DBAS(1,JDFL,2)!DWX
C$$$       GRADU3(2)=GRADU3(2) + DU3(JDFL)*DBAS(1,JDFL,3)!DWY
C$$$       GRADU3(3)=GRADU3(3) + DU3(JDFL)*DBAS(1,JDFL,4)!DWZ

C$$$ 220  CONTINUE

C$$$C
C$$$      DO 230 JDOFE=1,IDFL
C$$$       JDFL=KDFL(JDOFE)! local number of basic function
       
C$$$       HBASJ2=DBAS(1,JDFL,2)
C$$$       HBASJ3=DBAS(1,JDFL,3)
C$$$       HBASJ4=DBAS(1,JDFL,4)

C$$$       DO 240 IDOFE=1,IDFL
C$$$        IF (IDOFE.EQ.JDOFE) THEN
C$$$         AH   = HBASJ2*HBASJ2+HBASJ3*HBASJ3+HBASJ4*HBASJ4
C$$$         AH11 = (AH + HBASJ2*HBASJ2)
C$$$         AH22 = (AH + HBASJ3*HBASJ3)
C$$$         AH33 = (AH + HBASJ4*HBASJ4)
C$$$         AH12 = (     HBASJ2*HBASJ3)
C$$$         AH13 = (     HBASJ2*HBASJ4)
C$$$         AH23 = (     HBASJ3*HBASJ4)
C$$$         AH21 = AH12
C$$$         AH31 = AH13
C$$$         AH32 = AH23
C$$$        ELSE
C$$$         IDOFEH=KDFL(IDOFE)
C$$$         HBASI2=DBAS(1,IDOFEH,2)
C$$$         HBASI3=DBAS(1,IDOFEH,3)
C$$$         HBASI4=DBAS(1,IDOFEH,4)
C$$$         AH   = HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4
C$$$         AH11 = (AH + HBASJ2*HBASI2)
C$$$         AH22 = (AH + HBASJ3*HBASI3)
C$$$         AH33 = (AH + HBASJ4*HBASI4)
C$$$         AH12 = (     HBASI2*HBASJ3)
C$$$         AH13 = (     HBASI2*HBASJ4)
C$$$         AH23 = (     HBASI3*HBASJ4)
C$$$         AH21 = (     HBASI3*HBASJ2)
C$$$         AH31 = (     HBASI4*HBASJ2)
C$$$         AH32 = (     HBASI4*HBASJ3)
C$$$        ENDIF
C$$$        S11(JDOFE,IDOFE)=S11(JDOFE,IDOFE)+dVisc*OM*AH11
C$$$        S22(JDOFE,IDOFE)=S22(JDOFE,IDOFE)+dVisc*OM*AH22
C$$$        S33(JDOFE,IDOFE)=S33(JDOFE,IDOFE)+dVisc*OM*AH33
C$$$        S12(JDOFE,IDOFE)=S12(JDOFE,IDOFE)+dVisc*OM*AH12
C$$$        S13(JDOFE,IDOFE)=S13(JDOFE,IDOFE)+dVisc*OM*AH13
C$$$        S23(JDOFE,IDOFE)=S23(JDOFE,IDOFE)+dVisc*OM*AH23
C$$$        S21(JDOFE,IDOFE)=S21(JDOFE,IDOFE)+dVisc*OM*AH21
C$$$        S31(JDOFE,IDOFE)=S31(JDOFE,IDOFE)+dVisc*OM*AH31
C$$$        S32(JDOFE,IDOFE)=S32(JDOFE,IDOFE)+dVisc*OM*AH32
C$$$ 240  CONTINUE
C$$$ 230  CONTINUE
C$$$C
C$$$ 200  CONTINUE
C$$$C
C$$$      DO 300 JDOFE=1,IDFL
C$$$      DO 300 IDOFE=1,IDFL
C$$$        IA    =KENTRY(JDOFE,IDOFE)
C$$$        DS11(IA)=DS11(IA) + S11(JDOFE,IDOFE)
C$$$        DS22(IA)=DS22(IA) + S22(JDOFE,IDOFE)
C$$$        DS33(IA)=DS33(IA) + S33(JDOFE,IDOFE)
C$$$        DS12(IA)=DS12(IA) + S12(JDOFE,IDOFE)
C$$$        DS13(IA)=DS13(IA) + S13(JDOFE,IDOFE)
C$$$        DS23(IA)=DS23(IA) + S23(JDOFE,IDOFE)
C$$$        DS21(IA)=DS21(IA) + S21(JDOFE,IDOFE)
C$$$        DS31(IA)=DS31(IA) + S31(JDOFE,IDOFE)
C$$$        DS32(IA)=DS32(IA) + S32(JDOFE,IDOFE)
C$$$ 300  CONTINUE
C$$$C
C$$$ 100  CONTINUE
C$$$C
C$$$99999 END
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------


************************************************************************
      SUBROUTINE CUBATURESTRESS(U1,U2,U3,
     *           DS11,DS22,DS33,DS12,DS13,DS23,DS21,DS31,DS32,
     *           NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,DCORVG,ELE)
************************************************************************
*
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
      REAL*8    DS11(*),DS22(*),DS33(*)
      REAL*8    DS12(*),DS13(*),DS23(*),DS21(*),DS31(*),DS32(*)
      REAL*8    U1(*),U2(*),U3(*), dETA0,dExp, dVisc, dShear
C
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS)
      REAL*8    S11(NNBAS,NNBAS),S22(NNBAS,NNBAS),S33(NNBAS,NNBAS)
      REAL*8    S12(NNBAS,NNBAS),S13(NNBAS,NNBAS),S23(NNBAS,NNBAS)
      REAL*8    S21(NNBAS,NNBAS),S31(NNBAS,NNBAS),S32(NNBAS,NNBAS)
      REAL*8    ViscosityModel
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
      DIMENSION DU3(NNBAS), GRADU3(NNDIM)
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
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      S11(JDOFE,JDOFE)=0D0
      S22(JDOFE,JDOFE)=0D0
      S33(JDOFE,JDOFE)=0D0
      S12(JDOFE,JDOFE)=0D0
      S13(JDOFE,JDOFE)=0D0
      S23(JDOFE,JDOFE)=0D0
      S21(JDOFE,JDOFE)=0D0
      S31(JDOFE,JDOFE)=0D0
      S32(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      S11(JDOFE,IDOFE)=0D0
      S22(JDOFE,IDOFE)=0D0
      S33(JDOFE,IDOFE)=0D0
      S12(JDOFE,IDOFE)=0D0
      S13(JDOFE,IDOFE)=0D0
      S23(JDOFE,IDOFE)=0D0
      S21(JDOFE,IDOFE)=0D0
      S31(JDOFE,IDOFE)=0D0
      S32(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
!      if (myid.eq.1) then
!      DO idofe=1,idfl
!       write(*,'(50I10)') myid,iel, idofe, kdfl(idofe),kdfg(idofe)
!      end do
!      end if
!      pause

      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
120   CONTINUE
C
!---============================---
      DO 130 JDOFE=1,IDFL
      JDFG=KDFG(JDOFE)
      JDFL=KDFL(JDOFE)
!!!   local = global      
      DU1(JDFL) = U1(JDFG) 
      DU2(JDFL) = U2(JDFG)
      DU3(JDFL) = U3(JDFG)

 130  CONTINUE      
! ---===========================---
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      IF (Transform%ilint.eq.2) THEN ! Q2
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
      IF (Transform%ilint.eq.1) THEN ! Q1
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

      GRADU2(1)=0D0!V
      GRADU2(2)=0D0
      GRADU2(3)=0D0

      GRADU3(1)=0D0!W
      GRADU3(2)=0D0
      GRADU3(3)=0D0 
C
      DO 220 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)! local number of basic function
       
       GRADU1(1)=GRADU1(1) + DU1(JDFL)*DBAS(1,JDFL,2)!DUX
       GRADU1(2)=GRADU1(2) + DU1(JDFL)*DBAS(1,JDFL,3)!DUY
       GRADU1(3)=GRADU1(3) + DU1(JDFL)*DBAS(1,JDFL,4)!DUZ

       GRADU2(1)=GRADU2(1) + DU2(JDFL)*DBAS(1,JDFL,2)!DVX
       GRADU2(2)=GRADU2(2) + DU2(JDFL)*DBAS(1,JDFL,3)!DVY
       GRADU2(3)=GRADU2(3) + DU2(JDFL)*DBAS(1,JDFL,4)!DVZ

       GRADU3(1)=GRADU3(1) + DU3(JDFL)*DBAS(1,JDFL,2)!DWX
       GRADU3(2)=GRADU3(2) + DU3(JDFL)*DBAS(1,JDFL,3)!DWY
       GRADU3(3)=GRADU3(3) + DU3(JDFL)*DBAS(1,JDFL,4)!DWZ

 220  CONTINUE

C ----=============================================---- 
       dShearSquare = GRADU1(1)**2d0 + GRADU2(2)**2d0 
     *        + GRADU3(3)**2d0 + 0.5d0*(GRADU1(2)+GRADU2(1))**2d0
     *        + 0.5d0*(GRADU1(3)+GRADU3(1))**2d0 
     *        + 0.5d0*(GRADU2(3)+GRADU3(2))**2d0

       dVisc = ViscosityModel(dShearSquare)
!       dVisc = HogenPowerlaw(dShearSquare)
C ----=============================================---- 
C
      DO 230 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)! local number of basic function
       
       HBASJ2=DBAS(1,JDFL,2)
       HBASJ3=DBAS(1,JDFL,3)
       HBASJ4=DBAS(1,JDFL,4)

       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH   = HBASJ2*HBASJ2+HBASJ3*HBASJ3+HBASJ4*HBASJ4
         AH11 = (AH + HBASJ2*HBASJ2)
         AH22 = (AH + HBASJ3*HBASJ3)
         AH33 = (AH + HBASJ4*HBASJ4)
         AH12 = (     HBASJ2*HBASJ3)
         AH13 = (     HBASJ2*HBASJ4)
         AH23 = (     HBASJ3*HBASJ4)
         AH21 = AH12
         AH31 = AH13
         AH32 = AH23
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         AH   = HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4
         AH11 = (AH + HBASJ2*HBASI2)
         AH22 = (AH + HBASJ3*HBASI3)
         AH33 = (AH + HBASJ4*HBASI4)
         AH12 = (     HBASI2*HBASJ3)
         AH13 = (     HBASI2*HBASJ4)
         AH23 = (     HBASI3*HBASJ4)
         AH21 = (     HBASI3*HBASJ2)
         AH31 = (     HBASI4*HBASJ2)
         AH32 = (     HBASI4*HBASJ3)
        ENDIF
        S11(JDOFE,IDOFE)=S11(JDOFE,IDOFE)+dVisc*OM*AH11
        S22(JDOFE,IDOFE)=S22(JDOFE,IDOFE)+dVisc*OM*AH22
        S33(JDOFE,IDOFE)=S33(JDOFE,IDOFE)+dVisc*OM*AH33
        S12(JDOFE,IDOFE)=S12(JDOFE,IDOFE)+dVisc*OM*AH12
        S13(JDOFE,IDOFE)=S13(JDOFE,IDOFE)+dVisc*OM*AH13
        S23(JDOFE,IDOFE)=S23(JDOFE,IDOFE)+dVisc*OM*AH23
        S21(JDOFE,IDOFE)=S21(JDOFE,IDOFE)+dVisc*OM*AH21
        S31(JDOFE,IDOFE)=S31(JDOFE,IDOFE)+dVisc*OM*AH31
        S32(JDOFE,IDOFE)=S32(JDOFE,IDOFE)+dVisc*OM*AH32
 240  CONTINUE
 230  CONTINUE
C
 200  CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DS11(IA)=DS11(IA) + S11(JDOFE,IDOFE)
        DS22(IA)=DS22(IA) + S22(JDOFE,IDOFE)
        DS33(IA)=DS33(IA) + S33(JDOFE,IDOFE)
        DS12(IA)=DS12(IA) + S12(JDOFE,IDOFE)
        DS13(IA)=DS13(IA) + S13(JDOFE,IDOFE)
        DS23(IA)=DS23(IA) + S23(JDOFE,IDOFE)
        DS21(IA)=DS21(IA) + S21(JDOFE,IDOFE)
        DS31(IA)=DS31(IA) + S31(JDOFE,IDOFE)
        DS32(IA)=DS32(IA) + S32(JDOFE,IDOFE)
 300  CONTINUE
C
 100  CONTINUE
C
99999 END
C
C
C





!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!************************************************************************
!      SUBROUTINE CUBATURESTRESS(U1,U2,U3,
!     *           DS11,DS22,DS33,DS12,DS13,DS23,DS21,DS31,DS32,
!     *           NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,DCORVG,ELE)
!************************************************************************
!*
!*-----------------------------------------------------------------------
!      USE PP3D_MPI, ONLY:myid
!
!      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
!      CHARACTER SUB*6,FMT*15,CPARAM*120
!C
!      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
!     *           NNDIM=3,NNCOF=10)
!      PARAMETER (Q2=0.5D0,Q8=0.125D0)
!C
!      REAL*8    DS11(*),DS22(*),DS33(*)
!      REAL*8    DS12(*),DS13(*),DS23(*),DS21(*),DS31(*),DS32(*)
!      REAL*8    U1(*),U2(*),U3(*), dETA0,dExp, dVisc, dShear
!C
!      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
!      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
!      DIMENSION KENTRY(NNBAS,NNBAS)
!      REAL*8    S11(NNBAS,NNBAS),S22(NNBAS,NNBAS),S33(NNBAS,NNBAS)
!      REAL*8    S12(NNBAS,NNBAS),S13(NNBAS,NNBAS),S23(NNBAS,NNBAS)
!      REAL*8    S21(NNBAS,NNBAS),S31(NNBAS,NNBAS),S32(NNBAS,NNBAS)
!      REAL*8    PolyFLOW_Carreau
!C
!      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
!      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
!      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
!      DIMENSION DU3(NNBAS), GRADU3(NNDIM)
!C
!      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
!      COMMON /ERRCTL/ IER,ICHECK
!      COMMON /CHAR/   SUB,FMT(3),CPARAM
!      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
!     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
!     *                IEL,NDIM
!      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
!     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
!      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
!      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
!
!      COMMON /COAUX1/ KDFG,KDFL,IDFL
!C
!C *** user COMMON blocks
!      INTEGER  VIPARM 
!      DIMENSION VIPARM(100)
!      EQUIVALENCE (IAUSAV,VIPARM)
!      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
!     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
!     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
!     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
!     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
!     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
!C
!      SAVE
!C
!      DO 1 I= 1,NNDER
!1     BDER(I)=.FALSE.
!C
!      DO 2 I=1,4
!2     BDER(I)=.TRUE.
!C
!      IELTYP=-1
!      CALL ELE(0D0,0D0,0D0,IELTYP)
!      IDFL=NDFL(IELTYP)
!C
!      ICUB=9
!      CALL CB3H(ICUB)
!      IF (IER.NE.0) GOTO 99999
!C
!************************************************************************
!C *** Calculation of the matrix - storage technique 7 or 8
!************************************************************************
!      ICUBP=ICUB
!      CALL ELE(0D0,0D0,0D0,-2)
!C
!C *** Loop over all elements
!      DO 100 IEL=1,NEL
!C
!C
!      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
!      IF (IER.LT.0) GOTO 99999
!C
!C *** Determine entry positions in matrix
!      DO 110 JDOFE=1,IDFL
!      ILD=KLDA(KDFG(JDOFE))
!      KENTRY(JDOFE,JDOFE)=ILD
!      S11(JDOFE,JDOFE)=0D0
!      S22(JDOFE,JDOFE)=0D0
!      S33(JDOFE,JDOFE)=0D0
!      S12(JDOFE,JDOFE)=0D0
!      S13(JDOFE,JDOFE)=0D0
!      S23(JDOFE,JDOFE)=0D0
!      S21(JDOFE,JDOFE)=0D0
!      S31(JDOFE,JDOFE)=0D0
!      S32(JDOFE,JDOFE)=0D0
!      JCOL0=ILD
!      DO 111 IDOFE=1,IDFL
!      IF (IDOFE.EQ.JDOFE) GOTO 111
!      IDFG=KDFG(IDOFE)
!      DO 112 JCOL=JCOL0,NA
!      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
!112   CONTINUE
!113   JCOL0=JCOL+1
!      KENTRY(JDOFE,IDOFE)=JCOL
!      S11(JDOFE,IDOFE)=0D0
!      S22(JDOFE,IDOFE)=0D0
!      S33(JDOFE,IDOFE)=0D0
!      S12(JDOFE,IDOFE)=0D0
!      S13(JDOFE,IDOFE)=0D0
!      S23(JDOFE,IDOFE)=0D0
!      S21(JDOFE,IDOFE)=0D0
!      S31(JDOFE,IDOFE)=0D0
!      S32(JDOFE,IDOFE)=0D0
!111   CONTINUE
!110   CONTINUE
!C
!C *** Evaluation of coordinates of the vertices
!!      if (myid.eq.1) then
!!      DO idofe=1,idfl
!!       write(*,'(50I10)') myid,iel, idofe, kdfl(idofe),kdfg(idofe)
!!      end do
!!      end if
!!      pause
!
!      DO 120 IVE=1,NVE
!      JP=KVERT(IVE,IEL)
!      KVE(IVE)=JP
!      DX(IVE)=DCORVG(1,JP)
!      DY(IVE)=DCORVG(2,JP)
!      DZ(IVE)=DCORVG(3,JP)
!120   CONTINUE
!C
!      DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
!      DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
!      DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
!      DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
!      DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
!      DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
!      DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
!      DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
!      DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
!      DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
!      DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
!      DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
!      DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
!      DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
!      DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
!      DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
!      DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
!      DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
!      DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
!      DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
!      DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
!      DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
!      DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
!      DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
!C
!!---============================---
!      DO 130 JDOFE=1,IDFL
!      JDFG=KDFG(JDOFE)
!      JDFL=KDFL(JDOFE)
!!!!   local = global      
!      DU1(JDFL) = U1(JDFG) 
!      DU2(JDFL) = U2(JDFG)
!      DU3(JDFL) = U3(JDFG)
!
! 130  CONTINUE      
!! ---===========================---
!C *** Loop over all cubature points
!      DO 200 ICUBP=1,NCUBP
!C
!      XI1=DXI(ICUBP,1)
!      XI2=DXI(ICUBP,2)
!      XI3=DXI(ICUBP,3)
!C
!C *** Jacobian of the bilinear mapping onto the reference element
!      DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
!      DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
!      DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
!      DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
!      DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
!      DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
!      DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
!      DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
!      DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
!      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
!     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
!     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
!      OM=DOMEGA(ICUBP)*ABS(DETJ)
!C
!      CALL ELE(XI1,XI2,XI3,-3)
!      IF (IER.LT.0) GOTO 99999
!C 
!C ---=========================---
!      GRADU1(1)=0D0!U
!      GRADU1(2)=0D0
!      GRADU1(3)=0D0
!
!      GRADU2(1)=0D0!V
!      GRADU2(2)=0D0
!      GRADU2(3)=0D0
!
!      GRADU3(1)=0D0!W
!      GRADU3(2)=0D0
!      GRADU3(3)=0D0 
!C
!      DO 220 JDOFE=1,IDFL
!       JDFL=KDFL(JDOFE)! local number of basic function
!       
!       GRADU1(1)=GRADU1(1) + DU1(JDFL)*DBAS(1,JDFL,2)!DUX
!       GRADU1(2)=GRADU1(2) + DU1(JDFL)*DBAS(1,JDFL,3)!DUY
!       GRADU1(3)=GRADU1(3) + DU1(JDFL)*DBAS(1,JDFL,4)!DUZ
!
!       GRADU2(1)=GRADU2(1) + DU2(JDFL)*DBAS(1,JDFL,2)!DVX
!       GRADU2(2)=GRADU2(2) + DU2(JDFL)*DBAS(1,JDFL,3)!DVY
!       GRADU2(3)=GRADU2(3) + DU2(JDFL)*DBAS(1,JDFL,4)!DVZ
!
!       GRADU3(1)=GRADU3(1) + DU3(JDFL)*DBAS(1,JDFL,2)!DWX
!       GRADU3(2)=GRADU3(2) + DU3(JDFL)*DBAS(1,JDFL,3)!DWY
!       GRADU3(3)=GRADU3(3) + DU3(JDFL)*DBAS(1,JDFL,4)!DWZ
!
! 220  CONTINUE
!
!C ----=============================================---- 
!       dShearSquare = GRADU1(1)**2d0 + GRADU2(2)**2d0 
!     *        + GRADU3(3)**2d0 + 0.5d0*(GRADU1(2)+GRADU2(1))**2d0
!     *        + 0.5d0*(GRADU1(3)+GRADU3(1))**2d0 
!     *        + 0.5d0*(GRADU2(3)+GRADU3(2))**2d0
!
!       dVisc = PolyFLOW_Carreau(dShearSquare)
!!       dVisc = HogenPowerlaw(dShearSquare)
!C ----=============================================---- 
!C
!      DO 230 JDOFE=1,IDFL
!       JDFL=KDFL(JDOFE)! local number of basic function
!       
!       HBASJ2=DBAS(1,JDFL,2)
!       HBASJ3=DBAS(1,JDFL,3)
!       HBASJ4=DBAS(1,JDFL,4)
!
!       DO 240 IDOFE=1,IDFL
!        IF (IDOFE.EQ.JDOFE) THEN
!         AH   = HBASJ2*HBASJ2+HBASJ3*HBASJ3+HBASJ4*HBASJ4
!         AH11 = (AH + HBASJ2*HBASJ2)
!         AH22 = (AH + HBASJ3*HBASJ3)
!         AH33 = (AH + HBASJ4*HBASJ4)
!         AH12 = (     HBASJ2*HBASJ3)
!         AH13 = (     HBASJ2*HBASJ4)
!         AH23 = (     HBASJ3*HBASJ4)
!         AH21 = AH12
!         AH31 = AH13
!         AH32 = AH23
!        ELSE
!         IDOFEH=KDFL(IDOFE)
!         HBASI2=DBAS(1,IDOFEH,2)
!         HBASI3=DBAS(1,IDOFEH,3)
!         HBASI4=DBAS(1,IDOFEH,4)
!         AH   = HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4
!         AH11 = (AH + HBASJ2*HBASI2)
!         AH22 = (AH + HBASJ3*HBASI3)
!         AH33 = (AH + HBASJ4*HBASI4)
!         AH12 = (     HBASI2*HBASJ3)
!         AH13 = (     HBASI2*HBASJ4)
!         AH23 = (     HBASI3*HBASJ4)
!         AH21 = (     HBASI3*HBASJ2)
!         AH31 = (     HBASI4*HBASJ2)
!         AH32 = (     HBASI4*HBASJ3)
!        ENDIF
!        S11(JDOFE,IDOFE)=S11(JDOFE,IDOFE)+dVisc*OM*AH11
!        S22(JDOFE,IDOFE)=S22(JDOFE,IDOFE)+dVisc*OM*AH22
!        S33(JDOFE,IDOFE)=S33(JDOFE,IDOFE)+dVisc*OM*AH33
!        S12(JDOFE,IDOFE)=S12(JDOFE,IDOFE)+dVisc*OM*AH12
!        S13(JDOFE,IDOFE)=S13(JDOFE,IDOFE)+dVisc*OM*AH13
!        S23(JDOFE,IDOFE)=S23(JDOFE,IDOFE)+dVisc*OM*AH23
!        S21(JDOFE,IDOFE)=S21(JDOFE,IDOFE)+dVisc*OM*AH21
!        S31(JDOFE,IDOFE)=S31(JDOFE,IDOFE)+dVisc*OM*AH31
!        S32(JDOFE,IDOFE)=S32(JDOFE,IDOFE)+dVisc*OM*AH32
! 240  CONTINUE
! 230  CONTINUE
!C
! 200  CONTINUE
!C
!      DO 300 JDOFE=1,IDFL
!      DO 300 IDOFE=1,IDFL
!        IA    =KENTRY(JDOFE,IDOFE)
!        DS11(IA)=DS11(IA) + S11(JDOFE,IDOFE)
!        DS22(IA)=DS22(IA) + S22(JDOFE,IDOFE)
!        DS33(IA)=DS33(IA) + S33(JDOFE,IDOFE)
!        DS12(IA)=DS12(IA) + S12(JDOFE,IDOFE)
!        DS13(IA)=DS13(IA) + S13(JDOFE,IDOFE)
!        DS23(IA)=DS23(IA) + S23(JDOFE,IDOFE)
!        DS21(IA)=DS21(IA) + S21(JDOFE,IDOFE)
!        DS31(IA)=DS31(IA) + S31(JDOFE,IDOFE)
!        DS32(IA)=DS32(IA) + S32(JDOFE,IDOFE)
! 300  CONTINUE
!C
! 100  CONTINUE
!C
!99999 END
C
C
C
! ************************************************************************
!       SUBROUTINE STRESSTENSOR(DVISC,DS11,DS22,DS33,DS12,DS13,DS23,NA,
!      *           KCOLA,KLDA,KVERT,KAREA,KEDGE,DCORVG,ELE)
! ************************************************************************
! *
! *-----------------------------------------------------------------------
!       IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
!       CHARACTER SUB*6,FMT*15,CPARAM*120
! C
!       PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
!      *           NNDIM=3,NNCOF=10)
!       PARAMETER (Q2=0.5D0,Q8=0.125D0)
! C
!       REAL*8    DS11(*),DS22(*),DS33(*),DS12(*),DS13(*),DS23(*)
!       DIMENSION KCOLA(*),KLDA(*),DVISC(*),DCORVG(NNDIM,*)
!       DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
!       DIMENSION KENTRY(NNBAS,NNBAS)
!       REAL*8    S11(NNBAS,NNBAS),S22(NNBAS,NNBAS),S33(NNBAS,NNBAS)
!       REAL*8    S12(NNBAS,NNBAS),S13(NNBAS,NNBAS),S23(NNBAS,NNBAS)
! C
!       DIMENSION KDFG(NNBAS),KDFL(NNBAS)
!       DIMENSION DU(NNDIM,NNBAS),DEF(NNDIM,NNBAS),GRADU(NNDIM,NNDIM)
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
!       DVISCOSITY=DVISC(IEL)
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
! C *** Evaluation of coordinates of the vertices
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
!       CALL ELE(XI1,XI2,XI3,-3)
!       IF (IER.LT.0) GOTO 99999
! C
! C *** Summing up over all pairs of multiindices
!       DO 230 JDOFE=1,IDFL
!        JDOFEH=KDFL(JDOFE)
!        HBASJ2=DBAS(1,JDOFEH,2)
!        HBASJ3=DBAS(1,JDOFEH,3)
!        HBASJ4=DBAS(1,JDOFEH,4)
! C
!        DO 240 IDOFE=1,IDFL
!         IF (IDOFE.EQ.JDOFE) THEN
!          AH   = HBASJ2*HBASJ2+HBASJ3*HBASJ3+HBASJ4*HBASJ4
!          AH11 = (AH + HBASJ2*HBASJ2)
!          AH22 = (AH + HBASJ3*HBASJ3)
!          AH33 = (AH + HBASJ4*HBASJ4)
!          AH12 = (     HBASJ2*HBASJ3)
!          AH13 = (     HBASJ2*HBASJ4)
!          AH23 = (     HBASJ3*HBASJ4)
!         ELSE
!          IDOFEH=KDFL(IDOFE)
!          HBASI2=DBAS(1,IDOFEH,2)
!          HBASI3=DBAS(1,IDOFEH,3)
!          HBASI4=DBAS(1,IDOFEH,4)
!          AH   = HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4
!          AH11 = (AH + HBASJ2*HBASI2)
!          AH22 = (AH + HBASJ3*HBASI3)
!          AH33 = (AH + HBASJ4*HBASI4)
!          AH12 = (     HBASJ2*HBASI3)
!          AH13 = (     HBASJ2*HBASI4)
!          AH23 = (     HBASJ3*HBASI4)
!         ENDIF
!         S11(JDOFE,IDOFE)=S11(JDOFE,IDOFE)+OM*AH11
!         S22(JDOFE,IDOFE)=S22(JDOFE,IDOFE)+OM*AH22
!         S33(JDOFE,IDOFE)=S33(JDOFE,IDOFE)+OM*AH33
!         S12(JDOFE,IDOFE)=S12(JDOFE,IDOFE)+OM*AH12
!         S13(JDOFE,IDOFE)=S13(JDOFE,IDOFE)+OM*AH13
!         S23(JDOFE,IDOFE)=S23(JDOFE,IDOFE)+OM*AH23
! 240    CONTINUE
! 230   CONTINUE
! C
! 200   CONTINUE
! C
!       DO 300 JDOFE=1,IDFL
!       DO 300 IDOFE=1,IDFL
!         IA    =KENTRY(JDOFE,IDOFE)
!         DS11(IA)=DS11(IA) + DVISCOSITY*S11(JDOFE,IDOFE)
!         DS22(IA)=DS22(IA) + DVISCOSITY*S22(JDOFE,IDOFE)
!         DS33(IA)=DS33(IA) + DVISCOSITY*S33(JDOFE,IDOFE)
!         DS12(IA)=DS12(IA) + DVISCOSITY*S12(JDOFE,IDOFE)
!         DS13(IA)=DS13(IA) + DVISCOSITY*S13(JDOFE,IDOFE)
!         DS23(IA)=DS23(IA) + DVISCOSITY*S23(JDOFE,IDOFE)
! 300   CONTINUE
! C
! 100   CONTINUE
! C
! 99999 END
C
C
C
************************************************************************
      SUBROUTINE STRESS(U1,U2,U3,D1,D2,D3,DVISCOS,
     *           KVERT,KAREA,KEDGE,DCORVG,ELE)
************************************************************************
*
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
      DIMENSION U1(*),U2(*),U3(*),D1(*),D2(*),D3(*)
      DIMENSION DVISCOS(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DU(NNDIM,NNBAS),DEF(NNDIM,NNBAS),GRADU(NNDIM,NNDIM)

      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
      DIMENSION DU3(NNBAS), GRADU3(NNDIM)
      REAL*8    ViscosityModel
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
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      DVISCOSITY=DVISCOS(IEL)
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
      DO 130 JDOFE=1,IDFL
      JDFG=KDFG(JDOFE)
      DU(1,JDOFE)=U1(JDFG)
      DU(2,JDOFE)=U2(JDFG)
      DU(3,JDOFE)=U3(JDFG)
      DO 140 JDER=1,NNDIM
      DEF(JDER,JDOFE)=0D0
 140  CONTINUE
 130  CONTINUE
C
!---============================---
      DO 150 JDOFE=1,IDFL
      JDFG=KDFG(JDOFE)
      JDFL=KDFL(JDOFE)
!!!   local = global      
      DU1(JDFL) = U1(JDFG) 
      DU2(JDFL) = U2(JDFG)
      DU3(JDFL) = U3(JDFG)

 150  CONTINUE      
! ---===========================---
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      IF (Transform%ilint.eq.2) THEN ! Q2
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
      IF (Transform%ilint.eq.1) THEN ! Q1
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

      GRADU2(1)=0D0!V
      GRADU2(2)=0D0
      GRADU2(3)=0D0

      GRADU3(1)=0D0!W
      GRADU3(2)=0D0
      GRADU3(3)=0D0 
C
      DO 205 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)! local number of basic function
       
       GRADU1(1)=GRADU1(1) + DU1(JDFL)*DBAS(1,JDFL,2)!DUX
       GRADU1(2)=GRADU1(2) + DU1(JDFL)*DBAS(1,JDFL,3)!DUY
       GRADU1(3)=GRADU1(3) + DU1(JDFL)*DBAS(1,JDFL,4)!DUZ

       GRADU2(1)=GRADU2(1) + DU2(JDFL)*DBAS(1,JDFL,2)!DVX
       GRADU2(2)=GRADU2(2) + DU2(JDFL)*DBAS(1,JDFL,3)!DVY
       GRADU2(3)=GRADU2(3) + DU2(JDFL)*DBAS(1,JDFL,4)!DVZ

       GRADU3(1)=GRADU3(1) + DU3(JDFL)*DBAS(1,JDFL,2)!DWX
       GRADU3(2)=GRADU3(2) + DU3(JDFL)*DBAS(1,JDFL,3)!DWY
       GRADU3(3)=GRADU3(3) + DU3(JDFL)*DBAS(1,JDFL,4)!DWZ

 205  CONTINUE

C ----=============================================---- 
       dShearSquare = GRADU1(1)**2d0 + GRADU2(2)**2d0 
     *        + GRADU3(3)**2d0 + 0.5d0*(GRADU1(2)+GRADU2(1))**2d0
     *        + 0.5d0*(GRADU1(3)+GRADU3(1))**2d0 
     *        + 0.5d0*(GRADU2(3)+GRADU3(2))**2d0

       dVisc = ViscosityModel(dShearSquare)
!       dVisc = HogenPowerlaw(dShearSquare)
C ----=============================================---- 
      DO 210 JDER=1,NNDIM
      DO 210 IDER=1,NNDIM
      GRADU(IDER,JDER)=0D0
 210  CONTINUE
C
C
      DO 220 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       DO JDER=1,NNDIM
        DAUX=DU(JDER,JDOFE)
        DO IDER=1,NNDIM
         GRADU(IDER,JDER)=GRADU(IDER,JDER)+DAUX*DBAS(1,JDFL,IDER+1)
        ENDDO
       ENDDO
!        DNUP=DNUP+DNUE(JDOFE)*DBAS(1,JDFL,1)
220   CONTINUE
C
      DO 230 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       DO JDER=1,NNDIM
        DO IDER=1,NNDIM
         DAUX=dVisc*OM*(GRADU(IDER,JDER)+GRADU(JDER,IDER))
         DEF(JDER,JDOFE)=DEF(JDER,JDOFE)+DAUX*DBAS(1,JDFL,IDER+1)
        ENDDO
       ENDDO  
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
       JDFG=KDFG(JDOFE)
       D1(JDFG)=D1(JDFG)-THSTEP*DEF(1,JDOFE)
       D2(JDFG)=D2(JDFG)-THSTEP*DEF(2,JDOFE)
       D3(JDFG)=D3(JDFG)-THSTEP*DEF(3,JDOFE)
!        write(*,*) DEF(:,JDOFE)
300   CONTINUE
!       pause
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE GetGradVelo_rhs_sub(U,D1,D2,D3,
     *           KVERT,KAREA,KEDGE,DCORVG,ELE)
************************************************************************
*
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      REAL*8    U(*),D1(*),D2(*),D3(*)
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DU(NNBAS),DEF(NNDIM,NNBAS),GRADU(NNDIM)
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
      DO 130 JDOFE=1,IDFL
      JDFG=KDFG(JDOFE)
      DU(   JDOFE)=U(JDFG)
      DEF(1,JDOFE)=0D0
      DEF(2,JDOFE)=0D0
      DEF(3,JDOFE)=0D0
 130  CONTINUE
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
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      GRADU(1)=0D0
      GRADU(2)=0D0
      GRADU(3)=0D0
      DO 220 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       DAUX=DU(JDOFE)
       GRADU(1)=GRADU(1)+DAUX*DBAS(1,JDFL,2)
       GRADU(2)=GRADU(2)+DAUX*DBAS(1,JDFL,3)
       GRADU(3)=GRADU(3)+DAUX*DBAS(1,JDFL,4)
220   CONTINUE
C
      DO 230 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       DEF(1,JDOFE)=DEF(1,JDOFE)+OM*GRADU(1)*DBAS(1,JDFL,1)
       DEF(2,JDOFE)=DEF(2,JDOFE)+OM*GRADU(2)*DBAS(1,JDFL,1)
       DEF(3,JDOFE)=DEF(3,JDOFE)+OM*GRADU(3)*DBAS(1,JDFL,1)
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
       JDFG=KDFG(JDOFE)
       D1(JDFG)=D1(JDFG)+DEF(1,JDOFE)
       D2(JDFG)=D2(JDFG)+DEF(2,JDOFE)
       D3(JDFG)=D3(JDFG)+DEF(3,JDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE AssembleViscoStress(D1,D2,D3,bAlpha,U11,U22,U33,U12,
     * U13,U23,KVERT,KAREA,KEDGE,DCORVG,dVis,dt,dLambda,dFRC,ELE)
************************************************************************
*
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      REAL*8 U11(*),U22(*),U33(*),U12(*),U13(*),U23(*)
      REAL*8 DCORVG(NNDIM,*),D1(*),D2(*),D3(*),dFRC(3)
      REAL*8 PSI(6),TAU(6)
      LOGICAL bAlpha(*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
      REAL*8  DU(NNDIM,NNDIM,NNBAS)
      REAL*8  DEF(NNDIM,NNBAS),GRADU(NNDIM,NNDIM),DSS(NNDIM,NNDIM)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     * DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS

      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C*** user COMMON blocks
      INTEGER  VIPARM
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     * ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
       SAVE
C
!       return
C
       dFRC = 0d0
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
       DO 130 JDOFE=1,IDFL
       JDFG=KDFG(JDOFE)
       PSI(1) = U11(JDFG)
       PSI(2) = U22(JDFG)
       PSI(3) = U33(JDFG)
       PSI(4) = U12(JDFG)
       PSI(5) = U13(JDFG)
       PSI(6) = U23(JDFG)

       CALL ConvertPsiToTau(Psi,Tau)

       DU(1,1,JDOFE) = TAU(1)
       DU(2,2,JDOFE) = TAU(2)
       DU(3,3,JDOFE) = TAU(3)
       DU(1,2,JDOFE) = TAU(4)
       DU(2,1,JDOFE) = TAU(4)
       DU(1,3,JDOFE) = TAU(5)
       DU(3,1,JDOFE) = TAU(5)
       DU(2,3,JDOFE) = TAU(6)
       DU(3,2,JDOFE) = TAU(6)
       DO 140 JDER=1,NNDIM
       DEF(JDER,JDOFE)=0D0
  140  CONTINUE
  130  CONTINUE
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
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DALX=0D0     ! ALFA x deriv
      DALY=0D0     ! ALFA y deriv
      DALZ=0D0     ! ALFA z deriv
      DSS =0D0
      DO JDOFE=1,IDFL
       IG = KDFG(JDOFE)
       IL = KDFL(JDOFE)
       DBI1=DBAS(1,IL,1)
       DBI2=DBAS(1,IL,2)
       DBI3=DBAS(1,IL,3)
       DBI4=DBAS(1,IL,4)

       DSS(:,:) = DSS(:,:) + DU(:,:,JDOFE)*DBI1

      IF (bALPHA(IG)) THEN
        DALPHA = 1d0
       ELSE
        DALPHA = 0d0
       END IF
       DALX = DALX + DALPHA*DBI2
       DALY = DALY + DALPHA*DBI3
       DALZ = DALZ + DALPHA*DBI4
      END DO
      DN1=-DALX
      DN2=-DALY
      DN3=-DALZ
C
      GRADU = 0d0
      DO JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       DO JDER=1,NNDIM
        DO IDER=1,NNDIM
         DAUX=DU(JDER,IDER,JDOFE)
         GRADU(IDER,JDER) = GRADU(IDER,JDER) + DAUX*DBAS(1,JDFL,1)
        ENDDO
       ENDDO
      END DO
C
      DNY = dVis/dLambda
      AH1=DNY*((DSS(1,1)-1d0)*DN1+ DSS(1,2)*DN2     + DSS(1,3)*DN3)
      AH2=DNY* (DSS(2,1)*DN1     +(DSS(2,2)-1d0)*DN2+ DSS(2,3)*DN3)
      AH3=DNY* (DSS(3,1)*DN1     + DSS(3,2)*DN2 +(DSS(3,3)-1d0)*DN3)
C
      dFRC(1) = dFRC(1) + AH1*OM
      dFRC(2) = dFRC(2) + AH2*OM
      dFRC(3) = dFRC(3) + AH3*OM
C
      DO 230 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       DO JDER=1,NNDIM
        DO IDER=1,NNDIM
         DAUX= dVis*OM*GRADU(IDER,JDER)/dLambda
         DEF(JDER,JDOFE) = DEF(JDER,JDOFE) + DAUX*DBAS(1,JDFL,IDER+1)
        ENDDO
       ENDDO
230   CONTINUE
C
200   CONTINUE
C
       DO 300 JDOFE=1,IDFL
        JDFG=KDFG(JDOFE)
!        write(*,*) DEF(:,JDOFE),dVis,dLam
        D1(JDFG) = D1(JDFG) - dt*DEF(1,JDOFE)
        D2(JDFG) = D2(JDFG) - dt*DEF(2,JDOFE)
        D3(JDFG) = D3(JDFG) - dt*DEF(3,JDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
