************************************************************************
      SUBROUTINE Build_FluxQP1(U1,U2,U3,DPHI,DNORM,DMID,DF,IPE,DVP,DMP,
     *           KMP,KVERT,KAREA,KEDGE,KINT,DCORVG,KADJ,ICUB,IUPSAM,
     *           ELE1,ELE2,DBC,DTIME,DT,iEQ)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      use PP3D_MPI
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION U1(*),U2(*),U3(*),DPHI(*),DNORM(6,3,*),DF(*)
      DIMENSION DVP(4,*),DMP(3,*),KMP(5,*),DMID(3,*)
      DIMENSION KINT(*),DCORVG(NNDIM,*),IPE(*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*),KADJ(NNAE,*)
C
      DIMENSION KDFG1(NNBAS),KDFL1(NNBAS)
      DIMENSION KDFG2(NNBAS),KDFL2(NNBAS)
C
      INTEGER :: IUPSAM 
      REAL*8 DMyOmgP(NNCUBP),DMyCubP(NNCUBP,NNAE,NNDIM),DBC
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
      IELTYP1=-1
      CALL ELE1(0D0,0D0,0D0,IELTYP1)
      IDFL1=NDFL(IELTYP1)
C
      IELTYP2=-1
      CALL ELE2(0D0,0D0,0D0,IELTYP2)
      IDFL2=NDFL(IELTYP2)
C
!       CALL CB3H(ICUB)
      CALL SetUpMyCub(DMyOmgP,DMyCubP,NCUBP,ICUB)
      IF (IER.NE.0) GOTO 99999
C
      ICUBP=ICUB
      CALL ELE1(0D0,0D0,0D0,-2)
      CALL ELE2(0D0,0D0,0D0,-2)
C
      KPAR = 1
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      IF ((ABS(IPE(IEL)).EQ.100).OR.(iEQ.EQ.2.AND.IPE(IEL).EQ.0)) THEN
       DO IAT=1,6
        IAREA= KAREA(IAT,IEL)
        IF (IAREA.EQ.KMP(5,KPAR)) THEN
         KPAR = KPAR + 1
        END IF
       END DO
       GOTO 100
      END IF
C
      DDS = SIGN(1D0,DBLE(IPE(IEL)))
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
      CALL NDFGL(IEL,1,IELTYP1,KVERT,KEDGE,KAREA,KDFG1,KDFL1)
      IF (IER.LT.0) GOTO 99999
C
      CALL NDFGL(IEL,1,IELTYP2,KVERT,KEDGE,KAREA,KDFG2,KDFL2)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices
      DX0I = DMID(1,IEL)
      DY0I = DMID(2,IEL)
      DZ0I = DMID(3,IEL)
      DO 120 IVE=1,NVE
      IP=KVERT(IVE,IEL)
      KVE(IVE)=IP
      DX(IVE)=DCORVG(1,IP)
      DY(IVE)=DCORVG(2,IP)
      DZ(IVE)=DCORVG(3,IP)
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
C *** Loop over all faces
      DO 150 IAT=1,6
C
      DNX  = -DNORM(IAT,1,IEL) ! FaceArea*n_x
      DNY  = -DNORM(IAT,2,IEL) ! FaceArea*n_y
      DNZ  = -DNORM(IAT,3,IEL) ! FaceArea*n_z
C
      IADJ = KADJ (IAT,IEL)
      IAREA= KAREA(IAT,IEL)
      IF (IAREA.EQ.KMP(5,KPAR)) THEN
       JADJ=1
      ELSE
       JADJ=0
      END IF
C
      DO 200 ICUBP=1,NCUBP
C
      XI1=DMyCubP(ICUBP,IAT,1)
      XI2=DMyCubP(ICUBP,IAT,2)
      XI3=DMyCubP(ICUBP,IAT,3)
      OM=DMyOmgP (ICUBP)! *area, but it is already included in the normal
C
! ----------------------------------------------------------------
!     Computation of the cartesian coordiante of the cubature point
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
      ENDIF
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
! ----------------------------------------------------------------
C
!---- Interpolate the velocity at the cubature point -------------
      CALL ELE2(XI1,XI2,XI3,0)
      DU1=0D0
      DU2=0D0
      DU3=0D0
      DO 220 JDOFE=1,IDFL2
       JDFL=KDFL2(JDOFE)
       JDFG=KDFG2(JDOFE)
       HBASJ1=DBAS(1,JDFL,1)
       DU1=DU1+U1(JDFG)*HBASJ1
       DU2=DU2+U2(JDFG)*HBASJ1
       DU3=DU3+U3(JDFG)*HBASJ1
220   CONTINUE
      IF (iEQ.EQ.2) THEN
       DU1 = DDS*DU1
       DU2 = DDS*DU2
       DU3 = DDS*DU3
      END IF
C
! ----------------------------------------------------------------
C
      IF (IADJ.NE.0) THEN ! ONLY INTERNAL FACES !
C
!---- Midpoint of the neighboring element -------------------------
      DX0J = DMID(1,IADJ)
      DY0J = DMID(2,IADJ)
      DZ0J = DMID(3,IADJ)
C
C---- Interpolate Phi at the cubature point from the neigbor element
      CALL ELE1(XX-DX0J,YY-DY0J,ZZ-DZ0J,0)
      DPHIJ=0D0
      DO 222 JDOFE=1,IDFL1
       JDFL=JDOFE
       JDFG=JDOFE+4*(IADJ-1)
       HBASJ1=DBAS(1,JDFL,1)
       DPHIJ=DPHIJ+DPHI(JDFG)*HBASJ1
222   CONTINUE
C
      END IF
C
      IF (JADJ.NE.0) THEN ! ONLY PARALLEL FACES !
C
!---- Midpoint of the neighboring element -------------------------
      DX0J = DMP(1,KPAR)
      DY0J = DMP(2,KPAR)
      DZ0J = DMP(3,KPAR)
      DD=SQRT(DVP(2,KPAR)**2d0+DVP(3,KPAR)**2d0+DVP(4,KPAR)**2d0)
      IF (DD.LT.1d-4) THEN
!        pause
       JADJ = -1
      END IF
C
C---- Interpolate Phi at the cubature point from the neigbor element
      CALL ELE1(XX-DX0J,YY-DY0J,ZZ-DZ0J,0)
      DPHIJ=0D0
      DO 223 JDOFE=1,IDFL1
       JDFL=JDOFE
       HBASJ1=DBAS(1,JDFL,1)
       DPHIJ=DPHIJ+DVP(JDFL,KPAR)*HBASJ1
223   CONTINUE
C
      END IF
! ----------------------------------------------------------------
C
!---- Interpolate Phi at the cubature point from the current element
      CALL ELE1(XX-DX0I,YY-DY0I,ZZ-DZ0I,0)
      DPHII=0D0
      DO 224 JDOFE=1,IDFL1
       JDFL=JDOFE
       JDFG=JDOFE+4*(IEL -1)
       HBASJ1=DBAS(1,JDFL,1)
       DPHII=DPHII+DPHI(JDFG)*HBASJ1
224   CONTINUE
! ----------------------------------------------------------------
C
      DUN    = DU1*DNX+DU2*DNY+DU3*DNZ
      IF (IADJ.NE.0.AND.JADJ.NE.0) THEN
       WRITE(*,*) "IADJ.NE.0.AND.JADJ.NE.0 in FLUX"
       STOP
      END IF

!       IF (iEQ.EQ.1) THEN
!       IF ((0.3d0-ZZ.LT.0.5d-4).and.(IPE(IEL).GE.0).and.dun.lt.1d-5)THEN
!         write(*,'(I5,20(G12.4))') myid,iel,DUN,du3,xx,yy,zz,icubp,ncubp
!       END IF
!       END IF

      DPHIIJ = DPHII
      IF (DUN.LT.0d0) THEN
       IF (IADJ.NE.0) THEN
!         IF (ABS(IPE(IADJ)).NE.100) THEN
         DPHIIJ = DPHIJ
!         END IF
       ELSE
        IF (0.3d0-ZZ.LT. 0.5d-4) THEN
!          IF (IPE(IEL).GE.0) THEN
          DPHIIJ = DBC(XX,YY,ZZ,T)
!          END IF
        END IF
!         IF (ZZ+0.3d0.LT.0.5d-4) THEN
! !          IF (IPE(IEL).GE.0) THEN
!           IF (IPE(IEL).GT.0) DPHIIJ = MAX(DPHII, 0.0015)
!           IF (IPE(IEL).LT.0) DPHIIJ = MIN(DPHII,-0.0015)
! !          END IF
!         END IF
       END IF
       IF (JADJ.NE.0) THEN
         DPHIIJ = DPHIJ
       END IF
      END IF


!       IF ((IADJ.NE.0.OR.JADJ.NE.0)) THEN ! Internal faces
!        IF (DUN.GT.0d0) THEN
!         DPHIIJ = DPHII
!        ELSE
!         DPHIIJ = DPHIJ
!         !!!! this part here may not be bugfree !!!!
!         IF (JADJ.LT.0) DPHIIJ = DPHII
!         IF (ABS(IPE(IADJ)).EQ.100) DPHIIJ = DPHII
!        ENDIF
!       ELSE                               ! Physical boundary faces
! !        IF (DUN.GT.0d0) THEN
!          DPHIIJ = DPHII                  ! Outflow
! !        ELSE
! !          DPHIIJ = DBC(XX,YY,ZZ,DTIME)    ! Inflow
! !        END IF
! 
!       END IF

C *** Summing up over all pairs of multiindices
!       IF (ABS(DPHIIJ).GT.0.2) THEN
!        write(*,*) DPHIIJ,IADJ,IPE(IADJ)
!          pause
!       END IF

      DO 230 JDOFE=1,IDFL1
       JDFL=KDFL1(JDOFE)
       JDFG=KDFG1(JDOFE)
       HBASJ1=DBAS(1,JDFL,1)
       AH=HBASJ1*DUN*DPHIIJ
       DF(JDFG) = DF(JDFG) - OM*AH*DT
230   CONTINUE
C
200   CONTINUE ! ICUBP
C
      IF (JADJ.NE.0) KPAR = KPAR + 1
C
150   CONTINUE ! IAT
C
100   CONTINUE ! IEL
C
!       IF (iEQ.EQ.1) THEN
!       pause
!       end if

99999 END
C
C
C
      SUBROUTINE SetUpMyCub(DMyOmgP,DMyCubP,MyNCubP,ICUB)
      IMPLICIT NONE
      INTEGER NNCUBP,NNAE,NNDIM
      PARAMETER (NNCUBP=36,NNAE=6,NNDIM=3)
      INTEGER ICUB,MyNCubP
      REAL*8 DMyOmgP(NNCUBP),DMyCubP(NNCUBP,NNAE,NNDIM)
C
      IF (ICUB.EQ.1) THEN
C
       MyNCubP      =  1
       DMyOmgP(1)   =  1d0

       DMyCubP(1,1,1) =  0d0
       DMyCubP(1,1,2) =  0d0
       DMyCubP(1,1,3) = -1d0

       DMyCubP(1,2,1) =  0d0
       DMyCubP(1,2,2) = -1d0
       DMyCubP(1,2,3) =  0d0

       DMyCubP(1,3,1) =  1d0
       DMyCubP(1,3,2) =  0d0
       DMyCubP(1,3,3) =  0d0

       DMyCubP(1,4,1) =  0d0
       DMyCubP(1,4,2) =  1d0
       DMyCubP(1,4,3) =  0d0

       DMyCubP(1,5,1) = -1d0
       DMyCubP(1,5,2) =  0d0
       DMyCubP(1,5,3) =  0d0

       DMyCubP(1,6,1) =  0d0
       DMyCubP(1,6,2) =  0d0
       DMyCubP(1,6,3) =  1d0

       GOTO 9999
      ENDIF
C
      IF (ICUB.EQ.2) THEN

       MyNCubP      =  4
       DMyOmgP(1)   =  0.25d0
       DMyOmgP(2)   =  0.25d0
       DMyOmgP(3)   =  0.25d0
       DMyOmgP(4)   =  0.25d0

       DMyCubP(1,1,1) =  0.577350269189626D0
       DMyCubP(1,1,2) =  0.577350269189626D0
       DMyCubP(1,1,3) = -1d0
       DMyCubP(2,1,1) = -0.577350269189626D0
       DMyCubP(2,1,2) =  0.577350269189626D0
       DMyCubP(2,1,3) = -1d0
       DMyCubP(3,1,1) = -0.577350269189626D0
       DMyCubP(3,1,2) = -0.577350269189626D0
       DMyCubP(3,1,3) = -1d0
       DMyCubP(4,1,1) =  0.577350269189626D0
       DMyCubP(4,1,2) = -0.577350269189626D0
       DMyCubP(4,1,3) = -1d0

       DMyCubP(1,2,1) =  0.577350269189626D0
       DMyCubP(1,2,2) = -1d0
       DMyCubP(1,2,3) =  0.577350269189626D0
       DMyCubP(2,2,1) = -0.577350269189626D0
       DMyCubP(2,2,2) = -1d0
       DMyCubP(2,2,3) =  0.577350269189626D0
       DMyCubP(3,2,1) = -0.577350269189626D0
       DMyCubP(3,2,2) = -1d0
       DMyCubP(3,2,3) = -0.577350269189626D0
       DMyCubP(4,2,1) =  0.577350269189626D0
       DMyCubP(4,2,2) = -1d0
       DMyCubP(4,2,3) = -0.577350269189626D0

       DMyCubP(1,3,1) =  1d0
       DMyCubP(1,3,2) =  0.577350269189626D0
       DMyCubP(1,3,3) =  0.577350269189626D0
       DMyCubP(2,3,1) =  1d0
       DMyCubP(2,3,2) = -0.577350269189626D0
       DMyCubP(2,3,3) =  0.577350269189626D0
       DMyCubP(3,3,1) =  1d0
       DMyCubP(3,3,2) = -0.577350269189626D0
       DMyCubP(3,3,3) = -0.577350269189626D0
       DMyCubP(4,3,1) =  1d0
       DMyCubP(4,3,2) =  0.577350269189626D0
       DMyCubP(4,3,3) = -0.577350269189626D0

       DMyCubP(1,4,1) =  0.577350269189626D0
       DMyCubP(1,4,2) =  1d0
       DMyCubP(1,4,3) =  0.577350269189626D0
       DMyCubP(2,4,1) = -0.577350269189626D0
       DMyCubP(2,4,2) =  1d0
       DMyCubP(2,4,3) =  0.577350269189626D0
       DMyCubP(3,4,1) = -0.577350269189626D0
       DMyCubP(3,4,2) =  1d0
       DMyCubP(3,4,3) = -0.577350269189626D0
       DMyCubP(4,4,1) =  0.577350269189626D0
       DMyCubP(4,4,2) =  1d0
       DMyCubP(4,4,3) = -0.577350269189626D0

       DMyCubP(1,5,1) = -1d0
       DMyCubP(1,5,2) =  0.577350269189626D0
       DMyCubP(1,5,3) =  0.577350269189626D0
       DMyCubP(2,5,1) = -1d0
       DMyCubP(2,5,2) = -0.577350269189626D0
       DMyCubP(2,5,3) =  0.577350269189626D0
       DMyCubP(3,5,1) = -1d0
       DMyCubP(3,5,2) = -0.577350269189626D0
       DMyCubP(3,5,3) = -0.577350269189626D0
       DMyCubP(4,5,1) = -1d0
       DMyCubP(4,5,2) =  0.577350269189626D0
       DMyCubP(4,5,3) = -0.577350269189626D0

       DMyCubP(1,6,1) =  0.577350269189626D0
       DMyCubP(1,6,2) =  0.577350269189626D0
       DMyCubP(1,6,3) =  1d0
       DMyCubP(2,6,1) = -0.577350269189626D0
       DMyCubP(2,6,2) =  0.577350269189626D0
       DMyCubP(2,6,3) =  1d0
       DMyCubP(3,6,1) = -0.577350269189626D0
       DMyCubP(3,6,2) = -0.577350269189626D0
       DMyCubP(3,6,3) =  1d0
       DMyCubP(4,6,1) =  0.577350269189626D0
       DMyCubP(4,6,2) = -0.577350269189626D0
       DMyCubP(4,6,3) =  1d0

       GOTO 9999
      END IF
C
      IF (ICUB.EQ.3) THEN
       MyNCubP      =  9
       DMyOmgP(1)   =  0.077160493827161d0
       DMyOmgP(2)   =  0.077160493827161d0
       DMyOmgP(3)   =  0.077160493827161d0
       DMyOmgP(4)   =  0.077160493827161d0
       DMyOmgP(5)   =  0.123456790123457D0
       DMyOmgP(6)   =  0.123456790123457D0
       DMyOmgP(7)   =  0.123456790123457D0
       DMyOmgP(8)   =  0.123456790123457D0
       DMyOmgP(9)   =  0.197530864197531D0

       ! IAT = 1 !
       DMyCubP(1,1,1) =  0.774596669241483D0
       DMyCubP(1,1,2) =  0.774596669241483D0
       DMyCubP(1,1,3) = -1d0
       DMyCubP(2,1,1) = -0.774596669241483D0
       DMyCubP(2,1,2) =  0.774596669241483D0
       DMyCubP(2,1,3) = -1d0
       DMyCubP(3,1,1) =  0.774596669241483D0
       DMyCubP(3,1,2) = -0.774596669241483D0
       DMyCubP(3,1,3) = -1d0
       DMyCubP(4,1,1) = -0.774596669241483D0
       DMyCubP(4,1,2) = -0.774596669241483D0
       DMyCubP(4,1,3) = -1d0
       DMyCubP(5,1,1) =  0.774596669241483D0
       DMyCubP(5,1,2) =  0D0
       DMyCubP(5,1,3) = -1d0
       DMyCubP(6,1,1) = -0.774596669241483D0
       DMyCubP(6,1,2) =  0D0
       DMyCubP(6,1,3) = -1d0
       DMyCubP(7,1,1) =  0D0
       DMyCubP(7,1,2) =  0.774596669241483D0
       DMyCubP(7,1,3) = -1d0
       DMyCubP(8,1,1) =  0D0
       DMyCubP(8,1,2) = -0.774596669241483D0
       DMyCubP(8,1,3) = -1d0
       DMyCubP(9,1,1) =  0D0
       DMyCubP(9,1,2) =  0D0
       DMyCubP(9,1,3) = -1d0

       ! IAT = 2 !
       DMyCubP(1,2,1) =  0.774596669241483D0
       DMyCubP(1,2,2) = -1d0
       DMyCubP(1,2,3) =  0.774596669241483D0
       DMyCubP(2,2,1) = -0.774596669241483D0
       DMyCubP(2,2,2) = -1d0
       DMyCubP(2,2,3) =  0.774596669241483D0
       DMyCubP(3,2,1) =  0.774596669241483D0
       DMyCubP(3,2,2) = -1d0
       DMyCubP(3,2,3) = -0.774596669241483D0
       DMyCubP(4,2,1) = -0.774596669241483D0
       DMyCubP(4,2,2) = -1d0
       DMyCubP(4,2,3) = -0.774596669241483D0
       DMyCubP(5,2,1) =  0.774596669241483D0
       DMyCubP(5,2,2) = -1d0
       DMyCubP(5,2,3) =  0D0
       DMyCubP(6,2,1) = -0.774596669241483D0
       DMyCubP(6,2,2) = -1d0
       DMyCubP(6,2,3) =  0D0
       DMyCubP(7,2,1) =  0D0
       DMyCubP(7,2,2) = -1d0
       DMyCubP(7,2,3) =  0.774596669241483D0
       DMyCubP(8,2,1) =  0D0
       DMyCubP(8,2,2) = -1d0
       DMyCubP(8,2,3) = -0.774596669241483D0
       DMyCubP(9,2,1) =  0D0
       DMyCubP(9,2,2) = -1d0
       DMyCubP(9,2,3) =  0D0

       ! IAT = 3 !
       DMyCubP(1,3,1) =  1d0
       DMyCubP(1,3,2) =  0.774596669241483D0
       DMyCubP(1,3,3) =  0.774596669241483D0
       DMyCubP(2,3,1) =  1d0
       DMyCubP(2,3,2) =  0.774596669241483D0
       DMyCubP(2,3,3) = -0.774596669241483D0
       DMyCubP(3,3,1) =  1d0
       DMyCubP(3,3,2) = -0.774596669241483D0
       DMyCubP(3,3,3) =  0.774596669241483D0
       DMyCubP(4,3,1) =  1d0
       DMyCubP(4,3,2) = -0.774596669241483D0
       DMyCubP(4,3,3) = -0.774596669241483D0
       DMyCubP(5,3,1) =  1d0
       DMyCubP(5,3,2) =  0D0
       DMyCubP(5,3,3) =  0.774596669241483D0
       DMyCubP(6,3,1) =  1d0
       DMyCubP(6,3,2) =  0D0
       DMyCubP(6,3,3) = -0.774596669241483D0
       DMyCubP(7,3,1) =  1d0
       DMyCubP(7,3,2) =  0.774596669241483D0
       DMyCubP(7,3,3) =  0D0
       DMyCubP(8,3,1) =  1d0
       DMyCubP(8,3,2) = -0.774596669241483D0
       DMyCubP(8,3,3) =  0D0
       DMyCubP(9,3,1) =  1D0
       DMyCubP(9,3,2) =  0D0
       DMyCubP(9,3,3) =  0d0

       ! IAT = 4 !
       DMyCubP(1,4,1) =  0.774596669241483D0
       DMyCubP(1,4,2) =  1d0
       DMyCubP(1,4,3) =  0.774596669241483D0
       DMyCubP(2,4,1) = -0.774596669241483D0
       DMyCubP(2,4,2) =  1d0
       DMyCubP(2,4,3) =  0.774596669241483D0
       DMyCubP(3,4,1) =  0.774596669241483D0
       DMyCubP(3,4,2) =  1d0
       DMyCubP(3,4,3) = -0.774596669241483D0
       DMyCubP(4,4,1) = -0.774596669241483D0
       DMyCubP(4,4,2) =  1d0
       DMyCubP(4,4,3) = -0.774596669241483D0
       DMyCubP(5,4,1) =  0.774596669241483D0
       DMyCubP(5,4,2) =  1d0
       DMyCubP(5,4,3) =  0D0
       DMyCubP(6,4,1) = -0.774596669241483D0
       DMyCubP(6,4,2) =  1d0
       DMyCubP(6,4,3) =  0D0
       DMyCubP(7,4,1) =  0D0
       DMyCubP(7,4,2) =  1d0
       DMyCubP(7,4,3) =  0.774596669241483D0
       DMyCubP(8,4,1) =  0D0
       DMyCubP(8,4,2) =  1d0
       DMyCubP(8,4,3) = -0.774596669241483D0
       DMyCubP(9,4,1) =  0D0
       DMyCubP(9,4,2) =  1d0
       DMyCubP(9,4,3) =  0D0

       ! IAT = 5 !
       DMyCubP(1,5,1) = -1d0
       DMyCubP(1,5,2) =  0.774596669241483D0
       DMyCubP(1,5,3) =  0.774596669241483D0
       DMyCubP(2,5,1) = -1d0
       DMyCubP(2,5,2) =  0.774596669241483D0
       DMyCubP(2,5,3) = -0.774596669241483D0
       DMyCubP(3,5,1) = -1d0
       DMyCubP(3,5,2) = -0.774596669241483D0
       DMyCubP(3,5,3) =  0.774596669241483D0
       DMyCubP(4,5,1) = -1d0
       DMyCubP(4,5,2) = -0.774596669241483D0
       DMyCubP(4,5,3) = -0.774596669241483D0
       DMyCubP(5,5,1) = -1d0
       DMyCubP(5,5,2) =  0D0
       DMyCubP(5,5,3) =  0.774596669241483D0
       DMyCubP(6,5,1) = -1d0
       DMyCubP(6,5,2) =  0D0
       DMyCubP(6,5,3) = -0.774596669241483D0
       DMyCubP(7,5,1) = -1d0
       DMyCubP(7,5,2) =  0.774596669241483D0
       DMyCubP(7,5,3) =  0D0
       DMyCubP(8,5,1) = -1d0
       DMyCubP(8,5,2) = -0.774596669241483D0
       DMyCubP(8,5,3) =  0D0
       DMyCubP(9,5,1) = -1D0
       DMyCubP(9,5,2) =  0D0
       DMyCubP(9,5,3) =  0d0

       ! IAT = 6 !
       DMyCubP(1,6,1) =  0.774596669241483D0
       DMyCubP(1,6,2) =  0.774596669241483D0
       DMyCubP(1,6,3) =  1d0
       DMyCubP(2,6,1) = -0.774596669241483D0
       DMyCubP(2,6,2) =  0.774596669241483D0
       DMyCubP(2,6,3) =  1d0
       DMyCubP(3,6,1) =  0.774596669241483D0
       DMyCubP(3,6,2) = -0.774596669241483D0
       DMyCubP(3,6,3) =  1d0
       DMyCubP(4,6,1) = -0.774596669241483D0
       DMyCubP(4,6,2) = -0.774596669241483D0
       DMyCubP(4,6,3) =  1d0
       DMyCubP(5,6,1) =  0.774596669241483D0
       DMyCubP(5,6,2) =  0D0
       DMyCubP(5,6,3) =  1d0
       DMyCubP(6,6,1) = -0.774596669241483D0
       DMyCubP(6,6,2) =  0D0
       DMyCubP(6,6,3) =  1d0
       DMyCubP(7,6,1) =  0D0
       DMyCubP(7,6,2) =  0.774596669241483D0
       DMyCubP(7,6,3) =  1d0
       DMyCubP(8,6,1) =  0D0
       DMyCubP(8,6,2) = -0.774596669241483D0
       DMyCubP(8,6,3) =  1d0
       DMyCubP(9,6,1) =  0D0
       DMyCubP(9,6,2) =  0D0
       DMyCubP(9,6,3) =  1d0

       GOTO 9999
      END IF

      WRITE(*,*) "Not defined cubature formula ..."
      STOP
9999  CONTINUE
      END
C
C
C
      SUBROUTINE BuildRHS1(DF,IPE,DMID,DK,DM,DU,NEL,DT,IDEF)
      IMPLICIT NONE
      REAL*8 DU(*),DF(*),DK(4,4,*),DM(4,4,*),DT
      INTEGER IPE(*),NEL,IDEF
      REAL*8 DSUM,DS
      INTEGER I,J,IEL,II,JJ
      REAL*8 DMID(3,*)
      LOGICAL bBC


      DS = DBLE(IDEF)

      DO IEL=1,NEL

      bBC = .NOT.((DMID(3,IEL).GT.0.285).AND.(IPE(IEL).EQ.0))
!       IF (ABS(IPE(IEL)).NE.100.AND.bBC) THEN 
      IF (ABS(IPE(IEL)).NE.100) THEN 

       DO I=1,4
        DSUM = 0d0
        DO J=1,4
         JJ = 4*(IEL-1)+J
         DSUM = DSUM + (DS*DM(I,J,IEL) + DT*DK(I,J,IEL))*DU(JJ) !
        END DO
        II = 4*(IEL-1)+I
        DF(II) = DF(II) + DSUM
       END DO

      ELSE

      DO I=1,4
       II = 4*(IEL-1)+I
       DF(II) = 0d0
      END DO

      END IF

      END DO

      END
C
C
C
      SUBROUTINE BuildRHS2(DF,IPE,DK,DM,DU,NEL,DT,IDEF)
      IMPLICIT NONE
      REAL*8 DU(*),DF(*),DK(4,4,*),DM(4,4,*),DT
      INTEGER IPE(*),NEL,IDEF
      REAL*8 DSUM,DS
      INTEGER I,J,IEL,II,JJ


      DS = DBLE(IDEF)

      DO IEL=1,NEL

      IF (ABS(IPE(IEL)).NE.100.AND.IPE(IEL).NE.0) THEN 

       DO I=1,4
        DSUM = 0d0
        DO J=1,4
         JJ = 4*(IEL-1)+J
         DSUM = DSUM + (DS*DM(I,J,IEL) + DT*DK(I,J,IEL))*DU(JJ) !
        END DO
        II = 4*(IEL-1)+I
        DF(II) = DF(II) + DSUM
       END DO

      ELSE

      DO I=1,4
       II = 4*(IEL-1)+I
       DF(II) = 0d0
      END DO

      END IF

      END DO

      END
C
C
C
      SUBROUTINE imMatMultF1(DU,DF,IPE,DMID,DIM,NEL)
      IMPLICIT NONE
      REAL*8 DU(*),DF(*),DIM(4,4,*)
      INTEGER IPE(*),NEL
      REAL*8 DSUM
      INTEGER I,J,IEL,II,JJ
      REAL*8 DMID(3,*)
      LOGICAL bBC

      DO IEL=1,NEL

      bBC = .NOT.((DMID(3,IEL).GT.0.285).AND.(IPE(IEL).EQ.0))
!       IF (ABS(IPE(IEL)).NE.100.AND.bBC) THEN 
      IF (ABS(IPE(IEL)).NE.100) THEN 
       DO I=1,4
        DSUM = 0d0
        DO J=1,4
         JJ = 4*(IEL-1)+J
         DSUM = DSUM + DIM(I,J,IEL)*DF(JJ)
        END DO
        II = 4*(IEL-1)+I
        DU(II) = DSUM
       END DO

      ELSE

      DO I=1,4
       II = 4*(IEL-1)+I
       DU(II) = 0d0
      END DO

      END IF

      END DO

      END
C
C
C
      SUBROUTINE imMatMultF2(DU,DF,IPE,DIM,NEL)
      IMPLICIT NONE
      REAL*8 DU(*),DF(*),DIM(4,4,*)
      INTEGER IPE(*),NEL
      REAL*8 DSUM
      INTEGER I,J,IEL,II,JJ

      DO IEL=1,NEL

      IF (ABS(IPE(IEL)).NE.100.AND.IPE(IEL).NE.0) THEN 
       DO I=1,4
        DSUM = 0d0
        DO J=1,4
         JJ = 4*(IEL-1)+J
         DSUM = DSUM + DIM(I,J,IEL)*DF(JJ)
        END DO
        II = 4*(IEL-1)+I
        DU(II) = DSUM
       END DO

      ELSE

      DO I=1,4
       II = 4*(IEL-1)+I
       DU(II) = 0d0
      END DO

      END IF

      END DO

      END
C
C
C
      SUBROUTINE imMatMultF3(DU,DF,IPE,DIM,NEL)
      IMPLICIT NONE
      REAL*8 DU(*),DF(*),DIM(4,4,*)
      INTEGER IPE(*),NEL
      REAL*8 DSUM
      INTEGER I,J,IEL,II,JJ

      DO IEL=1,NEL

!       IF (ABS(IPE(IEL)).LT.100) THEN
       DO I=1,4
        DSUM = 0d0
        DO J=1,4
         JJ = 4*(IEL-1)+J
         DSUM = DSUM + DIM(I,J,IEL)*DF(JJ)
        END DO
        II = 4*(IEL-1)+I
        DU(II) = DSUM
       END DO

!       ELSE

!       DO I=1,4
!        II = 4*(IEL-1)+I
!        DU(II) = 0d0
!       END DO

!       END IF

      END DO

      END

