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
      USE var_QuadScalar_newton, ONLY:zeitstep,itns,tsm

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*8  DEFU(*),DEFV(*),DEFW(*),DEFP(*)
      REAL*8  VALU(*),VALV(*),VALW(*),VALP(*)
      INTEGER NDOF
      REAL*8 dt

	IF (tsm.EQ.0 .OR. itns.EQ.1) THEN
		dt = tstep
	ELSE 
		dt = zeitstep
	END IF

      IF (myid.eq.0) GOTO 1


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
      SUBROUTINE CC_Extraction_sub(KVERT,KAREA,KEDGE,NEL,KNU,KNV,KNW,
     *				iStringPos,iString)
      USE var_QuadScalar,ONLY : CC_EMat,
     *    mg_qMat,mg_qlMat,mg_BXmat,mg_BYmat,mg_BZmat,
     *    mg_lqMat,mg_BTXmat,mg_BTYmat,mg_BTZmat,NLMAX,TSTEP,ILEV

      USE var_QuadScalar_newton, ONLY : 
     *    mg_AA11mat,mg_AA22mat,mg_AA33mat,
     *    mg_AA12mat,mg_AA13mat,mg_AA23mat,
     *    mg_AA21mat,mg_AA31mat,mg_AA32mat

      USE PP3D_MPI, ONLY: MGE013,myid,master
      USE UMFPackSolver_CC, ONLY : myUmfPack_CCFactorize
      USE var_QuadScalar_newton, ONLY:zeitstep,itns,tsm

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      INTEGER KNU(*),KNV(*),KNW(*),iStringPos,iString
      INTEGER KVERT(8,*),KAREA(6,*),KEDGE(12,*),NEL
      INTEGER KDFG1(27),KDFL1(27),KDFG2(27),KDFL2(27)
      INTEGER LDA(86),COLA(85,85)
      REAL*8 dt

	IF (tsm.EQ.0 .OR. itns.EQ.1) THEN
		dt = tstep
	ELSE 
		dt = zeitstep
	END IF

      IF (myid.eq.0) GOTO 1


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
     *      mg_AA11mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a( 0 + IL,27 + JL) = 
     *      mg_AA12mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a( 0 + IL,54 + JL) = 
     *      mg_AA13mat(ILEV)%a(LG)

            CC_EMat(ILEV)%E(IEL)%a(27 + IL, 0 + JL) = 
     *      mg_AA21mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a(27 + IL,27 + JL) = 
     *      mg_AA22mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a(27 + IL,54 + JL) = 
     *      mg_AA23mat(ILEV)%a(LG)

            CC_EMat(ILEV)%E(IEL)%a(54 + IL, 0 + JL) = 
     *      mg_AA31mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a(54 + IL,27 + JL) = 
     *      mg_AA32mat(ILEV)%a(LG)
            CC_EMat(ILEV)%E(IEL)%a(54 + IL,54 + JL) = 
     *      mg_AA33mat(ILEV)%a(LG)
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


!!!!!!!!!!!! Visual OUTPUT !!!!!!!!!!!!!!!!!!!!
	   if (20*iel/nel.ge.iString) THEN
	    DO i=iStringPos+1,20*iel/nel
	     IF (myid.eq.1) write(*,'(A$)') '%'
	    END DO
	    iString=20*iel/nel+1
	    iStringPos=iString-1
	   END IF
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

************************************************************************
      SUBROUTINE IntQ2toP1_New(DEF,DU,KVERT,KAREA,KEDGE,
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
      REAL*8  DEF(4,*),DU(*),DP1BAS(4),DCORVG(NNDIM,*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
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
      DMX = 0d0 
      DMY = 0d0 
      DMZ = 0d0 
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
      DMX = DMX + 0.125d0*DX(IVE)
      DMY = DMY + 0.125d0*DY(IVE)
      DMZ = DMZ + 0.125d0*DZ(IVE)
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
      XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
      YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
      ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2

      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DVAL = 0d0
      DO JDFL=1,IDFL
       IL = KDFL(JDFL)
       IG = KDFG(JDFL)
       DVAL = DVAL + DBAS(1,IL,1)*DU(IG)
      END DO
C
      DP1BAS(1) = 1d0
      DP1BAS(2) = (XX-DMX)
      DP1BAS(3) = (YY-DMY)
      DP1BAS(4) = (ZZ-DMZ)
C
C *** Summing up over all pairs of multiindices
      DO 230 JDFL=1,4
        DEF(JDFL,IEL) = DEF(JDFL,IEL) + OM*DP1BAS(JDFL)*DVAL
230   CONTINUE
C
200   CONTINUE
C
100   CONTINUE
C
99999 END
************************************************************************
      SUBROUTINE IntP1toQ2_cc(DEF,DML,DP,KVERT,KAREA,KEDGE,
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
      DIMENSION DEF(*),DML(*),DP(4,*)
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
      DMX = 0d0 
      DMY = 0d0 
      DMZ = 0d0 
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
      DMX = DMX + 0.125d0*DX(IVE)
      DMY = DMY + 0.125d0*DY(IVE)
      DMZ = DMZ + 0.125d0*DZ(IVE)
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
      XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
      YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
      ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2

      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DVAL = DP(1,IEL) + DP(2,IEL)*(XX-DMX) + DP(3,IEL)*(YY-DMY)
     *                 + DP(3,IEL)*(ZZ-DMZ)  

C *** Summing up over all pairs of multiindices
      DO 230 JDFL=1,IDFL
        IG=KDFG(JDFL)
        IL=KDFL(JDFL)
        HBAS=DBAS(1,IL,1)
        DEF(IG) = DEF(IG) + OM*HBAS*DVAL
        DML(IG) = DML(IG) + OM*HBAS
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
      SUBROUTINE GetForceCyl_cc(U1,U2,U3,P,bALPHA,KVERT,KAREA,KEDGE,
     *                     DCORVG,DResForce,ELE)
************************************************************************
*     Discrete convection operator: Q1~ elements (nonparametric)
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      REAL*8  U1(*),U2(*),U3(*),P(*),DCORVG(NNDIM,*)
      REAL*8  DResForce(3)
      LOGICAL bALPHA(*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
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
     * ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
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
      DResForce(1) = 0D0
      DResForce(2) = 0D0
      DResForce(3) = 0D0
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
      DX0 = 0d0
      DY0 = 0d0
      DZ0 = 0d0
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
      DX0 = DX0 + 0.125d0*DX(IVE)
      DY0 = DY0 + 0.125d0*DY(IVE)
      DZ0 = DZ0 + 0.125d0*DZ(IVE)
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
C ----=============================================----
       dShearSquare = DU1X**2d0 + DU2Y**2d0 + DU3Z**2d0
     *        + 0.5d0*(DU1Y+DU2X)**2d0
     *        + 0.5d0*(DU1Z+DU3X)**2d0
     *        + 0.5d0*(DU2Z+DU3Y)**2d0

       dVisc = HogenPowerlaw(dShearSquare)
C ----=============================================----

       JJ = 4*(IEL-1) + 1
       Press =          P(JJ  ) + (XX-DX0)*P(JJ+1) +
     *         (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)
C--------------------------------------------------------
c-----------Form the integrand------------------
       DN1=-DALX
       DN2=-DALY
       DN3=-DALZ
c-------------------Acting force-------------------------
c--------------Deformation calculation-------------
C
       AH1=-Press*DN1 + dVisc*(DU1X*DN1 + DU1Y*DN2 + DU1Z*DN3)
       AH2=-Press*DN2 + dVisc*(DU2X*DN1 + DU2Y*DN2 + DU2Z*DN3)
       AH3=-Press*DN3 + dVisc*(DU3X*DN1 + DU3Y*DN2 + DU3Z*DN3)
C
       DResForce(1) = DResForce(1) + AH1*OM
       DResForce(2) = DResForce(2) + AH2*OM
       DResForce(3) = DResForce(3) + AH3*OM
C
200   CONTINUE
C
100   CONTINUE
C
999   CALL COMM_SUMMN(DResForce,3)
C
99999 CONTINUE

      END 



