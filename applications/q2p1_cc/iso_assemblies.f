!
! ----------------------------------------------
! MMat
! ----------------------------------------------
!
************************************************************************
      SUBROUTINE BuildMRhoMat_iso(DENS,DA,NA,KCOLA,KLDA,KVERT,KAREA,
     *                  KEDGE,DCORVG,ICUB,ELE)
************************************************************************
*     Discrete diffusion operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DENS(*),DA(*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
      BDER(1:4)=.TRUE.
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
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
C      
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      DENSITY = DENS(IEL)
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the mapping onto the reference element
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
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ=DBAS(1,JDOFEH,1)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=DENSITY*(HBASJ*HBASJ)
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI=DBAS(1,IDOFEH,1)
         AH=DENSITY*(HBASJ*HBASI)
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA)+DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
!
! ----------------------------------------------
! barMMat
! ----------------------------------------------
!!  CALL Build_barMMat_iso(mgDensity(ILEV)%x,qMat%na,qMat%ColA,qMat%LdA,&
!   mg_mesh%level(ILEV)%kvert,mg_mesh%level(ILEV)%karea,&
!   mg_mesh%level(ILEV)%kedge,mg_mesh%level(ILEV)%dcorvg,9,E013,&
!   mg_barM11mat(ILEV)%a,mg_barM12mat(ILEV)%a,mg_barM13mat(ILEV)%a,&
!   mg_barM21mat(ILEV)%a,mg_barM22mat(ILEV)%a,mg_barM23mat(ILEV)%a,&
!   mg_barM31mat(ILEV)%a,mg_barM32mat(ILEV)%a,mg_barM33mat(ILEV)%a,&
!   myScalar%valU,myScalar%valV,myScalar%valW)

************************************************************************
      SUBROUTINE Build_barMMat_iso(DENS,NA,KCOLA,KLDA,KVERT,KAREA,
     *                  KEDGE,DCORVG,ICUB,ELE,
     *                  DA11,DA12,DA13,DA21,DA22,DA23,DA31,DA32,DA33,
     *                  U1,U2,U3)
************************************************************************
*     Discrete diffusion operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DENS(*)

      DIMENSION DA11(*),DA12(*),DA13(*)
      DIMENSION DA21(*),DA22(*),DA23(*)
      DIMENSION DA31(*),DA32(*),DA33(*)

      REAL*8    M11(NNBAS,NNBAS),M22(NNBAS,NNBAS),M33(NNBAS,NNBAS)
      REAL*8    M12(NNBAS,NNBAS),M13(NNBAS,NNBAS),M23(NNBAS,NNBAS)
      REAL*8    M21(NNBAS,NNBAS),M31(NNBAS,NNBAS),M32(NNBAS,NNBAS)

      REAL*8    U1(*),U2(*),U3(*)
      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
      DIMENSION DU3(NNBAS), GRADU3(NNDIM)

      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
      BDER(1:4)=.TRUE.
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
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
C      
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      DENSITY = DENS(IEL)
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      M11(JDOFE,JDOFE)=0D0
      M22(JDOFE,JDOFE)=0D0
      M33(JDOFE,JDOFE)=0D0
      M12(JDOFE,JDOFE)=0D0
      M13(JDOFE,JDOFE)=0D0
      M23(JDOFE,JDOFE)=0D0
      M21(JDOFE,JDOFE)=0D0
      M31(JDOFE,JDOFE)=0D0
      M32(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      M11(JDOFE,IDOFE)=0D0
      M22(JDOFE,IDOFE)=0D0
      M33(JDOFE,IDOFE)=0D0
      M12(JDOFE,IDOFE)=0D0
      M13(JDOFE,IDOFE)=0D0
      M23(JDOFE,IDOFE)=0D0
      M21(JDOFE,IDOFE)=0D0
      M31(JDOFE,IDOFE)=0D0
      M32(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
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
C *** Jacobian of the mapping onto the reference element
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
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999

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
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ=DBAS(1,JDOFEH,1)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH11=DENSITY*(HBASJ*HBASJ)*GRADU1(1)
         AH12=DENSITY*(HBASJ*HBASJ)*GRADU1(2)
         AH13=DENSITY*(HBASJ*HBASJ)*GRADU1(3)
         AH21=DENSITY*(HBASJ*HBASJ)*GRADU2(1)
         AH22=DENSITY*(HBASJ*HBASJ)*GRADU2(2)
         AH23=DENSITY*(HBASJ*HBASJ)*GRADU2(3)
         AH31=DENSITY*(HBASJ*HBASJ)*GRADU3(1)
         AH32=DENSITY*(HBASJ*HBASJ)*GRADU3(2)
         AH33=DENSITY*(HBASJ*HBASJ)*GRADU3(3)
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI=DBAS(1,IDOFEH,1)
         AH11=DENSITY*(HBASJ*HBASI)*GRADU1(1)
         AH12=DENSITY*(HBASJ*HBASI)*GRADU1(2)
         AH13=DENSITY*(HBASJ*HBASI)*GRADU1(3)
         AH21=DENSITY*(HBASJ*HBASI)*GRADU2(1)
         AH22=DENSITY*(HBASJ*HBASI)*GRADU2(2)
         AH23=DENSITY*(HBASJ*HBASI)*GRADU2(3)
         AH31=DENSITY*(HBASJ*HBASI)*GRADU3(1)
         AH32=DENSITY*(HBASJ*HBASI)*GRADU3(2)
         AH33=DENSITY*(HBASJ*HBASI)*GRADU3(3)
        ENDIF
        M11(JDOFE,IDOFE)=M11(JDOFE,IDOFE)+OM*AH11
        M12(JDOFE,IDOFE)=M12(JDOFE,IDOFE)+OM*AH12
        M13(JDOFE,IDOFE)=M13(JDOFE,IDOFE)+OM*AH13
        M21(JDOFE,IDOFE)=M21(JDOFE,IDOFE)+OM*AH21
        M22(JDOFE,IDOFE)=M22(JDOFE,IDOFE)+OM*AH22
        M23(JDOFE,IDOFE)=M23(JDOFE,IDOFE)+OM*AH23
        M31(JDOFE,IDOFE)=M31(JDOFE,IDOFE)+OM*AH31
        M32(JDOFE,IDOFE)=M32(JDOFE,IDOFE)+OM*AH32
        M33(JDOFE,IDOFE)=M33(JDOFE,IDOFE)+OM*AH33
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA11(IA)=DA11(IA)+M11(JDOFE,IDOFE)
        DA12(IA)=DA12(IA)+M12(JDOFE,IDOFE)
        DA13(IA)=DA13(IA)+M13(JDOFE,IDOFE)
        DA21(IA)=DA21(IA)+M21(JDOFE,IDOFE)
        DA22(IA)=DA22(IA)+M22(JDOFE,IDOFE)
        DA23(IA)=DA23(IA)+M23(JDOFE,IDOFE)
        DA31(IA)=DA31(IA)+M31(JDOFE,IDOFE)
        DA32(IA)=DA32(IA)+M32(JDOFE,IDOFE)
        DA33(IA)=DA33(IA)+M33(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
!
! ----------------------------------------------
! KMat
! ----------------------------------------------
!
************************************************************************
      SUBROUTINE CONVQ2_iso(DENS,U1,U2,U3,dMeshVelo,DA,NU,KCOLA,
     *           KLDA,KVERT,KAREA,KEDGE,DCORVG,ELE)
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
      DIMENSION DENS(*),U1(*),U2(*),U3(*),DA(*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
      NA=KLDA(NU+1)-1
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
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      DENSITY=DENS(IEL)
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the mapping onto the reference element
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
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DU1=0D0
      DU2=0D0
      DU3=0D0
      DO 220 JDOFE=1,IDFL
        JDOFEH=KDFL(JDOFE)
        HBAS=DBAS(1,JDOFEH,1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDOFE)
         DU1=DU1+U1(JDFG)*HBAS
         DU2=DU2+U2(JDFG)*HBAS
         DU3=DU3+U3(JDFG)*HBAS
        ENDIF
220   CONTINUE
C
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(1,JDOFEH,1)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2+HBASJ4*DU3
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=HSUMJ*HBASJ1
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         HSUMI=HBASI2*DU1+HBASI3*DU2+HBASI4*DU3
         AH=HSUMI*HBASJ1
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+DENSITY*OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA)+DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
!
! ----------------------------------------------
! SMat
! ----------------------------------------------
!
      FUNCTION PolyFLOW_Carreau_iso(NormShearSquare)
      USE Transport_Q2P1, ONLY : Properties
      IMPLICIT NONE

      real*8 :: PolyFLOW_Carreau_iso
      real*8, intent (in) :: NormShearSquare

      REAL*8 :: dN

      dN = Properties%PowerLawExp-1d0
      PolyFLOW_Carreau_iso = 
     *Properties%Viscosity(1)*(1d-4 + NormShearSquare)**dN

!      WRITE(*,*) Properties%Viscosity(1), NormShearSquare,dN
      RETURN
      END
************************************************************************
      SUBROUTINE CUBATURESTRESS_iso(U1,U2,U3,
     *           DS11,DS22,DS33,DS12,DS13,DS23,DS21,DS31,DS32,
     *           NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,DCORVG,ELE)
************************************************************************
*
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid

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
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
      DIMENSION DU3(NNBAS), GRADU3(NNDIM)

      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
C *** Jacobian of the mapping onto the reference element
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
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
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

       dVisc = PolyFLOW_Carreau_iso(dShearSquare)
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
************************************************************************
      SUBROUTINE STRESSTENSOR_iso(DVISC,DS11,DS22,DS33,DS12,DS13,DS23,
     *           NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,DCORVG,ELE)
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
      REAL*8    DS11(*),DS22(*),DS33(*),DS12(*),DS13(*),DS23(*)
      DIMENSION KCOLA(*),KLDA(*),DVISC(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS)
      REAL*8    S11(NNBAS,NNBAS),S22(NNBAS,NNBAS),S33(NNBAS,NNBAS)
      REAL*8    S12(NNBAS,NNBAS),S13(NNBAS,NNBAS),S23(NNBAS,NNBAS)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DU(NNDIM,NNBAS),DEF(NNDIM,NNBAS),GRADU(NNDIM,NNDIM)

      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
      DVISCOSITY=DVISC(IEL)
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
111   CONTINUE
110   CONTINUE
C

C

C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the mapping onto the reference element
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
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH   = HBASJ2*HBASJ2+HBASJ3*HBASJ3+HBASJ4*HBASJ4
         AH11 = (AH + HBASJ2*HBASJ2)
         AH22 = (AH + HBASJ3*HBASJ3)
         AH33 = (AH + HBASJ4*HBASJ4)
         AH12 = (     HBASJ2*HBASJ3)
         AH13 = (     HBASJ2*HBASJ4)
         AH23 = (     HBASJ3*HBASJ4)
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         AH   = HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4
         AH11 = (AH + HBASJ2*HBASI2)
         AH22 = (AH + HBASJ3*HBASI3)
         AH33 = (AH + HBASJ4*HBASI4)
         AH12 = (     HBASJ2*HBASI3)
         AH13 = (     HBASJ2*HBASI4)
         AH23 = (     HBASJ3*HBASI4)
        ENDIF
        S11(JDOFE,IDOFE)=S11(JDOFE,IDOFE)+OM*AH11
        S22(JDOFE,IDOFE)=S22(JDOFE,IDOFE)+OM*AH22
        S33(JDOFE,IDOFE)=S33(JDOFE,IDOFE)+OM*AH33
        S12(JDOFE,IDOFE)=S12(JDOFE,IDOFE)+OM*AH12
        S13(JDOFE,IDOFE)=S13(JDOFE,IDOFE)+OM*AH13
        S23(JDOFE,IDOFE)=S23(JDOFE,IDOFE)+OM*AH23
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DS11(IA)=DS11(IA) + DVISCOSITY*S11(JDOFE,IDOFE)
        DS22(IA)=DS22(IA) + DVISCOSITY*S22(JDOFE,IDOFE)
        DS33(IA)=DS33(IA) + DVISCOSITY*S33(JDOFE,IDOFE)
        DS12(IA)=DS12(IA) + DVISCOSITY*S12(JDOFE,IDOFE)
        DS13(IA)=DS13(IA) + DVISCOSITY*S13(JDOFE,IDOFE)
        DS23(IA)=DS23(IA) + DVISCOSITY*S23(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE STRESS_iso(U1,U2,U3,D1,D2,D3,DVISCOS,
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
      DIMENSION U1(*),U2(*),U3(*),D1(*),D2(*),D3(*)
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DU(NNDIM,NNBAS),DEF(NNDIM,NNBAS),GRADU(NNDIM,NNDIM)

      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices

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
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      DO 210 JDER=1,NNDIM
      DO 210 IDER=1,NNDIM
      GRADU(IDER,JDER)=0D0
 210  CONTINUE
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
C ----=============================================---- 
       dSS = GRADU(1,1)**2d0 + GRADU(2,2)**2d0 + GRADU(3,3)**2d0
     *     + 0.5d0*(GRADU(1,2)+GRADU(2,1))**2d0
     *     + 0.5d0*(GRADU(1,3)+GRADU(3,1))**2d0 
     *     + 0.5d0*(GRADU(2,3)+GRADU(3,2))**2d0

       DVISCOSITY = PolyFLOW_Carreau_iso(dSS)
C ----=============================================---- 

      DO 230 JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       DO JDER=1,NNDIM
        DO IDER=1,NNDIM
         DAUX=DVISCOSITY*OM*(GRADU(IDER,JDER)+GRADU(JDER,IDER))
         DEF(JDER,JDOFE)=DEF(JDER,JDOFE)+DAUX*DBAS(1,JDFL,IDER+1)
        ENDDO
       ENDDO  
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
       JDFG=KDFG(JDOFE)
       D1(JDFG)=D1(JDFG)+THSTEP*DEF(1,JDOFE)
       D2(JDFG)=D2(JDFG)+THSTEP*DEF(2,JDOFE)
       D3(JDFG)=D3(JDFG)+THSTEP*DEF(3,JDOFE)
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
      SUBROUTINE GetGradVelo_rhs_sub_iso(U,D1,D2,D3,
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

      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C


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
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
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
      SUBROUTINE AssembleViscoStress_iso(D1,D2,D3,bAlpha,U11,U22,U33,
     * U12,U13,U23,KVERT,KAREA,KEDGE,DCORVG,dVis,dt,dLambda,dFRC,ELE)
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
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)

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
!
! ----------------------------------------------
! DiffMat
! ----------------------------------------------
!
************************************************************************
      SUBROUTINE DIFFQ2_NNEWT_iso(U1,U2,U3,DA,NA,KCOLA,KLDA,KVERT,KAREA,
     *                  KEDGE,DCORVG,ELE)
************************************************************************
*     Discrete diffusion operator: Q2 elements ---PREPARED !!
*-----------------------------------------------------------------------
      USE var_QuadScalar, ONLY : bNonNewtonian
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DA(*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8    U1(*),U2(*),U3(*), dVisc, dShear
      DIMENSION DU1(NNBAS), GRADU1(NNDIM)
      DIMENSION DU2(NNBAS), GRADU2(NNDIM)
      DIMENSION DU3(NNBAS), GRADU3(NNDIM)

      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
!       CALL LCL1(DA,NA)
! C
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
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C

C
C *** Loop over all cubature points
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
C
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the mapping onto the reference element
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
      IF (bNonNewtonian) THEN
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

       DVisco = PolyFLOW_Carreau_iso(dShearSquare)
C ----=============================================---- 
      END IF

C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=(HBASJ2*HBASJ2+HBASJ3*HBASJ3+HBASJ4*HBASJ4)
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         AH=(HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4)
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+DVisco*OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA)+DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
 
************************************************************************
      SUBROUTINE DIFFQ2_NEWT_iso(DA,NA,KCOLA,KLDA,KVERT,KAREA,
     *                       KEDGE,DCORVG,ELE)
************************************************************************
*     Discrete diffusion operator: Q2 elements ---PREPARED !!
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DA(*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)

      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
!       CALL LCL1(DA,NA)
! C
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
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C

C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the mapping onto the reference element
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
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=(HBASJ2*HBASJ2+HBASJ3*HBASJ3+HBASJ4*HBASJ4)
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         AH=(HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4)
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA) + DVisco*DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
 
!
! ----------------------------------------------
! BMat
! ----------------------------------------------
!
***************************************************
      SUBROUTINE Build_BTMatP1_iso(DAx,DAy,DAz,KLD,KCOL,KVERT,KAREA,
     *           KEDGE,DCORVG,NA,ICUB,ELE)
************************************************************************
*     Discrete convection operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid
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
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
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
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC=0d0
      DO JDOFE=1,IDFL2
       JDFL=KDFL2(JDOFE)
       JDFG=KDFG2(JDOFE)
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
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
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
      SUBROUTINE Build_BMatP1_iso(DAx,DAy,DAz,KLD,KCOL,KVERT,KAREA,
     *           KEDGE,DCORVG,NA,ICUB,ELE)
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
      DIMENSION DAx(*),DAy(*),DAz(*),DBAS1(4)
      DIMENSION DCORVG(NNDIM,*),KLD(*),KCOL(*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
C
      DIMENSION KDFG1(NNBAS),KDFL1(NNBAS)
      DIMENSION KDFG2(NNBAS),KDFL2(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRYx(NNBAS,NNBAS)
      DIMENSION DENTRYy(NNBAS,NNBAS),DENTRYz(NNBAS,NNBAS)
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
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
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC=0d0
      DO JDOFE=1,IDFL1
       JDFL=KDFL1(JDOFE)
       JDFG=KDFG1(JDOFE)
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
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
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

************************************************************************
      SUBROUTINE GetForceCyl_cc_iso(U1,U2,U3,P,bALPHA,KVERT,KAREA,KEDGE,
     *                     DCORVG,DResForce,ELE,bNonNewt)
************************************************************************
*     Discrete convection operator: Q1~ elements (nonparametric)
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      USE def_cc, ONLY : Properties
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      LOGICAL bNonNewt
      REAL*8  U1(*),U2(*),U3(*),P(*),DCORVG(NNDIM,*)
      REAL*8  DResForce(7)
      LOGICAL bALPHA(*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM)
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

       dVisc = PolyFLOW_Carreau_iso(dShearSquare)
       ELSE
       dVisc = Properties%Viscosity(1) 
       END IF
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
