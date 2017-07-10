************************************************************************
      SUBROUTINE Build_P1Mass(DA,KVERT,KAREA,KEDGE,DCORVG,ICUB,ELE)
************************************************************************
*     Discrete diffusion operator: Q1 elements
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DA(4,4,*)
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
      BDER(1)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
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
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
C
      IF (IER.LT.0) GOTO 99999
C
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
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
! ----------------------------------------------------------------
!     Computation of the cartesian coordiante of the cubature point
      XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
      YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
      ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
! ----------------------------------------------------------------
      CALL ELE(XX-DX0I,YY-DY0I,ZZ-DZ0I,0)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ=DBAS(1,JDOFEH,1)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=HBASJ*HBASJ
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI=DBAS(1,IDOFEH,1)
         AH=HBASJ*HBASI
        ENDIF
        DA(JDOFE,IDOFE,IEL)=DA(JDOFE,IDOFE,IEL)+OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
100   CONTINUE
C
99999 END
      SUBROUTINE CC_GetDefect_sub(DEFU,DEFV,DEFW,DEFP,
     *           VALU,VALV,VALW,VALP,NDOF,NEL)

      USE var_QuadScalar,ONLY : CC_EMat,mg_A11mat,mg_A22mat,mg_A33mat,
     *    mg_A12mat,mg_A13mat,mg_A23mat,mg_A21mat,mg_A31mat,mg_A32mat,
     *    mg_qMat,mg_qlMat,mg_BXmat,mg_BYmat,mg_BZmat,
     *    mg_lqMat,mg_BTXmat,mg_BTYmat,mg_BTZmat,NLMAX,TSTEP,ILEV
      USE PP3D_MPI, ONLY: MGE013,myid,master

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*8  DEFU(*),DEFV(*),DEFW(*),DEFP(*)
      REAL*8  VALU(*),VALV(*),VALW(*),VALP(*)
      INTEGER NDOF

      IF (myid.eq.0) GOTO 1

      dt = TSTEP

      DO i=1,NDOF

!        DEFU(i) = 0d0
!        DEFV(i) = 0d0
!        DEFW(i) = 0d0

 
!        DO j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
!         jVelo = mg_qMat(ILEV)%ColA(j)
!         DEFU(i) = DEFU(i) + mg_A11mat(ILEV)%a(j)*VALU(jVelo)
!         DEFV(i) = DEFV(i) + mg_A22mat(ILEV)%a(j)*VALV(jVelo)
!         DEFW(i) = DEFW(i) + mg_A33mat(ILEV)%a(j)*VALW(jVelo)
!        END DO

       DO j=mg_qlMat(ILEV)%LdA(i),mg_qlMat(ILEV)%LdA(i+1)-1
        jPres = mg_qlMat(ILEV)%ColA(j)
        dP    = VALP(jPres)
        DEFU(i) = DEFU(i) - dt*mg_BXmat(ILEV)%a(j)*dP
        DEFV(i) = DEFV(i) - dt*mg_BYmat(ILEV)%a(j)*dP
        DEFW(i) = DEFW(i) - dt*mg_BZmat(ILEV)%a(j)*dP
       END DO
      END DO

      DO i=1,4*NEL
       DEFP(i) = 0d0
       DO j=mg_lqMat(ILEV)%LdA(i),mg_lqMat(ILEV)%LdA(i+1)-1
        jVelo = mg_lqMat(ILEV)%ColA(j)
        DEFP(i) = DEFP(i) + mg_BTXmat(ILEV)%a(j)*VALU(jVelo)
     *                    + mg_BTYmat(ILEV)%a(j)*VALV(jVelo)    
     *                    + mg_BTZmat(ILEV)%a(j)*VALW(jVelo)    
       END DO
      END DO

1     CONTINUE

      END
      !
      !
      !
      SUBROUTINE CC_Extraction_sub(KVERT,KAREA,KEDGE,NEL,KNU,KNV,KNW)
      USE var_QuadScalar,ONLY : CC_EMat,mg_A11mat,mg_A22mat,mg_A33mat,
     *    mg_A12mat,mg_A13mat,mg_A23mat,mg_A21mat,mg_A31mat,mg_A32mat,
     *    mg_qMat,mg_qlMat,mg_BXmat,mg_BYmat,mg_BZmat,
     *    mg_lqMat,mg_BTXmat,mg_BTYmat,mg_BTZmat,NLMAX,TSTEP,ILEV
      USE PP3D_MPI, ONLY: MGE013,myid,master
      USE UMFPackSolver_CC, ONLY : myUmfPack_CCFactorize

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      INTEGER KNU(*),KNV(*),KNW(*)
      INTEGER KVERT(8,*),KAREA(6,*),KEDGE(12,*),NEL
      INTEGER KDFG1(27),KDFL1(27),KDFG2(27),KDFL2(27)
      INTEGER LDA(86),COLA(85,85)

      IF (myid.eq.0) GOTO 1

      dt = TSTEP

      DO i = 1, 86
       LdA(i) =  ((i-1)*85 + 1) - 1
      END DO
      DO i = 1, 85
       DO j = 1, 85
        ColA (i,j) = (i)-1
       END DO
      END DO

      DO IEL = 1,NEL

       CALL NDFGL(IEL,1,13,KVERT,KEDGE,KAREA,KDFG1,KDFL1)

       CC_EMat(ILEV)%E(IEL)%a = 0d0

       DO idof = 1,27
        IL = KDFL1(idof)
        IG = KDFG1(idof)
        DO LG = mg_qMat(ILEV)%LdA(IG),mg_qMat(ILEV)%LdA(IG+1)-1
         KG = mg_qMat(ILEV)%ColA(LG)
         DO KL=1,27
          IF (KDFG1(KL).EQ.KG) THEN
           JL = KDFL1(KL)

!            IF (IL.EQ.JL) THEN
! !             CC_EMat(ILEV)%E(IEL)%a( 0 + IL, 0 + JL) = 
! !      *      MGE013(ILEV)%UE11(IG)
! !             CC_EMat(ILEV)%E(IEL)%a(27 + IL,27 + JL) = 
! !      *      MGE013(ILEV)%UE22(IG)
! !             CC_EMat(ILEV)%E(IEL)%a(54 + IL,54 + JL) = 
! !      *      MGE013(ILEV)%UE33(IG)
! 
!             CC_EMat(ILEV)%E(IEL)%a( 0 + IL, 0 + JL) = 
!      *      MGE013(ILEV)%UE11(IG)
!             CC_EMat(ILEV)%E(IEL)%a( 0 + IL,27 + JL) = 
!      *      mg_A12mat(ILEV)%a(LG)
!             CC_EMat(ILEV)%E(IEL)%a( 0 + IL,54 + JL) = 
!      *      mg_A13mat(ILEV)%a(LG)
! 
!             CC_EMat(ILEV)%E(IEL)%a(27 + IL, 0 + JL) = 
!      *      mg_A21mat(ILEV)%a(LG)
!             CC_EMat(ILEV)%E(IEL)%a(27 + IL,27 + JL) = 
!      *      MGE013(ILEV)%UE22(IG)
!             CC_EMat(ILEV)%E(IEL)%a(27 + IL,54 + JL) = 
!      *      mg_A23mat(ILEV)%a(LG)
! 
!             CC_EMat(ILEV)%E(IEL)%a(54 + IL, 0 + JL) = 
!      *      mg_A31mat(ILEV)%a(LG)
!             CC_EMat(ILEV)%E(IEL)%a(54 + IL,27 + JL) = 
!      *      mg_A32mat(ILEV)%a(LG)
!             CC_EMat(ILEV)%E(IEL)%a(54 + IL,54 + JL) = 
!      *      MGE013(ILEV)%UE33(IG)
!      
!             ELSE
!             CC_EMat(ILEV)%E(IEL)%a( 0 + IL, 0 + JL) = 
!      *      mg_A11mat(ILEV)%a(LG)
!             CC_EMat(ILEV)%E(IEL)%a(27 + IL,27 + JL) = 
!      *      mg_A22mat(ILEV)%a(LG)
!             CC_EMat(ILEV)%E(IEL)%a(54 + IL,54 + JL) = 
!      *      mg_A33mat(ILEV)%a(LG)
! 
            CC_EMat(ILEV)%E(IEL)%a( 0 + IL, 0 + JL) = 
     *      mg_A11mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a( 0 + IL,27 + JL) = 
     *      mg_A12mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a( 0 + IL,54 + JL) = 
     *      mg_A13mat(ILEV)%a(LG)

            CC_EMat(ILEV)%E(IEL)%a(27 + IL, 0 + JL) = 
     *      mg_A21mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a(27 + IL,27 + JL) = 
     *      mg_A22mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a(27 + IL,54 + JL) = 
     *      mg_A23mat(ILEV)%a(LG)

            CC_EMat(ILEV)%E(IEL)%a(54 + IL, 0 + JL) = 
     *      mg_A31mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a(54 + IL,27 + JL) = 
     *      mg_A32mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a(54 + IL,54 + JL) = 
     *      mg_A33mat(ILEV)%a(LG)
!            END IF
          END IF
         END DO
        END DO

        DO LG = mg_qlMat(ILEV)%LdA(IG),mg_qlMat(ILEV)%LdA(IG+1)-1,4
         KG = mg_qlMat(ILEV)%ColA(LG)
         IF (KG.EQ.4*(IEL-1)+1) THEN
!           IF (myid.eq.1) WRITE(*,*) "yes!",IEL,KG
          CC_EMat(ILEV)%E(IEL)%a( 0 + IL, 82:85) =
     *    -mg_BXmat(ILEV)%a(LG:LG+3)*dt
          CC_EMat(ILEV)%E(IEL)%a(27 + IL, 82:85) =
     *    -mg_BYmat(ILEV)%a(LG:LG+3)*dt
          CC_EMat(ILEV)%E(IEL)%a(54 + IL, 82:85) =
     *    -mg_BZmat(ILEV)%a(LG:LG+3)*dt
         END IF
        END DO

        IF (KNU(IG).EQ.1) THEN
         CC_EMat(ILEV)%E(IEL)%a( 0 + IL, 1:85  ) = 0d0
         CC_EMat(ILEV)%E(IEL)%a( 0 + IL, 0 + IL) = 1d0
        END IF
        IF (KNV(IG).EQ.1) THEN
         CC_EMat(ILEV)%E(IEL)%a(27 + IL, 1:85  ) = 0d0
         CC_EMat(ILEV)%E(IEL)%a(27 + IL,27 + IL) = 1d0
        END IF
        IF (KNW(IG).EQ.1) THEN
         CC_EMat(ILEV)%E(IEL)%a(54 + IL, 1:85  ) = 0d0
         CC_EMat(ILEV)%E(IEL)%a(54 + IL,54 + IL) = 1d0
        END IF

       END DO

       DO IL = 1,4
!         IL = 1
        IG = 4*(IEL-1) + IL
        DO LG = mg_lqMat(ILEV)%LdA(IG),mg_lqMat(ILEV)%LdA(IG+1)-1        
         KG = mg_lqMat(ILEV)%ColA(LG)
         DO KL=1,27
          IF (KDFG1(KL).EQ.KG) THEN
           JL = KDFL1(KL)
           CC_EMat(ILEV)%E(IEL)%a(81 + IL, 0 + JL) = 
     *     mg_BTXmat(ILEV)%a(LG)
           CC_EMat(ILEV)%E(IEL)%a(81 + IL,27 + JL) =
     *     mg_BTYmat(ILEV)%a(LG)
           CC_EMat(ILEV)%E(IEL)%a(81 + IL,54 + JL) = 
     *     mg_BTZmat(ILEV)%a(LG)
          END IF
         END DO
        END DO
       END DO

!       if (myid.eq.1) then
!       open (file='tratrb.txt',unit=221)
!       do il=1,85
!       write(221,'(85ES12.4)') CC_EMat(ILEV)%E(IEL)%a(il,:)
!       end do
!       close(221)
!       end if
!       pause


       CALL myUmfPack_CCFactorize(CC_EMat(ILEV)%E(IEL)%a,LdA,ColA,
     *      CC_EMat(ILEV)%E(IEL)%H(1),CC_EMat(ILEV)%E(IEL)%H(2),85)

       IF (CC_EMat(ILEV)%E(IEL)%H(1).EQ.-1)
     *      then
        WRITE(*,*) "Sym: PROBLEM with the handle number ... ", iel
        pause
       END IF

       IF (CC_EMat(ILEV)%E(IEL)%H(2).EQ.-1)
     *      then
        WRITE(*,*) "Num: PROBLEM with the handle number ... ", iel
        pause
       END IF

      END DO

1     CONTINUE

      END
