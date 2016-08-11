MODULE def_SCALAR

USE def_FEAT
USE PP3D_MPI, ONLY:CommSum,CommMat,myid,showID

IMPLICIT NONE

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------

TYPE scalar
 INTEGER :: ndof,na
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: knpr
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: aux,rhs,def,val,val_old
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: src,snk
 REAL*4  , DIMENSION(:)  , ALLOCATABLE :: Amat
END TYPE scalar

TYPE mg_scalar
 TYPE(scalar), DIMENSION(:), ALLOCATABLE :: l
 REAL*8  defCrit,epsCrit,MinDef
 INTEGER NLmin,NLmax
 INTEGER SolvIter,SolvType
 CHARACTER cName*7
END TYPE mg_scalar

TYPE vector
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: val
END TYPE vector

TYPE mg_vector
 TYPE(vector), DIMENSION(:), ALLOCATABLE :: l
END TYPE mg_vector

INTEGER I,J,K
INTEGER NEDGE

! -------------- Subroutines -------------------
CONTAINS
! ----------------------------------------------
SUBROUTINE Init_General_Scalar(myScalar,Knpr_Scalar)

TYPE(mg_scalar), INTENT(INOUT) :: myScalar
EXTERNAL Knpr_Scalar

ALLOCATE (myScalar%l(NLMIN:NLMAX))

DO ILEV=NLMIN,NLMAX

 nat=KNAT(ILEV); myScalar%l(ILEV)%ndof=KNAT(ILEV)
 na =KNA(ILEV) ; myScalar%l(ILEV)%na  =KNA(ILEV)

 ALLOCATE (myScalar%l(ILEV)%Amat(na))

 ALLOCATE (myScalar%l(ILEV)%knpr(nat))
 ALLOCATE (myScalar%l(ILEV)%src(nat))
 ALLOCATE (myScalar%l(ILEV)%snk(nat))
 ALLOCATE (myScalar%l(ILEV)%val(nat))

END DO

ALLOCATE (myScalar%l(NLMAX)%aux(nat))
ALLOCATE (myScalar%l(NLMAX)%val_old(nat))
ALLOCATE (myScalar%l(NLMAX)%rhs(nat))
ALLOCATE (myScalar%l(NLMAX)%def(nat))

DO ILEV=NLMIN,NLMAX
 CALL Knpr_Scalar(myScalar%l(ILEV)%knpr,DWORK(L(KLCAG(ILEV))),&
      myScalar%l(ILEV)%ndof)
END DO

END SUBROUTINE Init_General_Scalar
!
! ----------------------------------------------
!
SUBROUTINE Matdef_General_Scalar(myScalar,idef,imat)
INTEGER :: idef,imat
TYPE(mg_scalar), INTENT(INOUT) :: myScalar

! IF (IUPW.EQ.0) THEN
 CALL MatdefAFC_Scalar(myScalar%l(ILEV)%val,&
 myScalar%l(ILEV)%def,myScalar%l(ILEV)%Amat, &
 DWORK(L(KLDK(ILEV))),DWORK(KST1),KWORK(KCOLA),&
 KWORK(KLDA),myScalar%l(ILEV)%ndof,&
 DWORK(L(LNUT)),DWORK(KM1),KWORK(L(LINOD)),&
 KWORK(L(LJNOD)),DWORK(L(LAEDGE)),idef,imat)
! ELSE
!  CALL Matdef_Scalar(myScalar%l(ILEV)%val,&
!  myScalar%l(ILEV)%def,myScalar%l(ILEV)%Amat, &
!  DWORK(L(KLDK(ILEV))),DWORK(KST1),KWORK(KCOLA),&
!  KWORK(KLDA),myScalar%l(ILEV)%ndof,&
!  DWORK(L(LNUT)),DWORK(KM1),idef)
! END IF


END SUBROUTINE Matdef_General_Scalar
!
! ----------------------------------------------
!
SUBROUTINE Matdef_Scalar(U,D,A,DK,DS,KCOLA,KLDA,NU,DNUT,DML,IDEF)

REAL*8  U(*),D(*),DK(*),DS(*),DNUT(*),DML(*)
REAL*4  A(*)
INTEGER KCOLA(*),KLDA(*),NU,IDEF
!---------------------------------------------------------------
REAL*8  DAUX,SIGML,DD
INTEGER II_LOC,IJ_LOC,ILOC

!     Perform global matrix assembly
IF (IDEF.LE.0) THEN
 DO I=1,NU

  II_LOC=KLDA(I)
  DAUX=DK(II_LOC)+DS(II_LOC)
  A(II_LOC)=REAL(DML(I))-REAL(THSTEP*DAUX)

  DO IJ_LOC=KLDA(I)+1,KLDA(I+1)-1
   DAUX=DK(IJ_LOC)+DS(IJ_LOC)
   A(IJ_LOC)=-REAL(THSTEP*DAUX)
  ENDDO
 ENDDO
END IF

!     Augment the defect if required
IF (IDEF.NE.0) THEN
 SIGML=DBLE(SIGN(1,IDEF))
 DO I=1,NU
  DAUX=SIGML*DML(I)
  DD=(DK(KLDA(I))+DS(KLDA(I)))*U(I)
  DO ILOC=KLDA(I)+1,KLDA(I+1)-1
   J=KCOLA(ILOC)
   DD=DD+(DK(ILOC)+DS(ILOC))*U(J)
  ENDDO
  D(I)=D(I)+DAUX*U(I)+THSTEP*DD
 ENDDO
ENDIF


END SUBROUTINE Matdef_Scalar
!
! ----------------------------------------------
!
SUBROUTINE InitializeAFC_General_Scalar()

! Construct the convection operator
CALL GetConvection_scalar(DWORK(KU1),DWORK(KU2),DWORK(KU3),&
     DWORK(L(KLDK(ILEV))),KWORK(KCOLA),KWORK(KLDA),&
     KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),&
     KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)))

CALL GetAFCStuctures_Scalar(DWORK(L(KLDK(ILEV))),KWORK(KCOLA),&
     KWORK(KLDA),NU,KWORK(L(LISEP)),KWORK(L(LIAUX)),&
     KWORK(L(LINOD)),KWORK(L(LJNOD)),DWORK(L(LAEDGE)),NEDGE)

END SUBROUTINE InitializeAFC_General_Scalar
!
! ----------------------------------------------
!
SUBROUTINE GetConvection_scalar(U1,U2,U3,DK,KCOLA,KLDA,KVERT,&
           KAREA,KEDGE,KINT,DCORVG)
REAL*8  U1(*),U2(*),U3(*),DK(*),DCORVG(*)
INTEGER KCOLA(*),KLDA(*),KVERT(*),KAREA(*),KEDGE(*),KINT(*)
INTEGER NA

NA=KLDA(NU+1)-1
IF (IELT.EQ.0) CALL ST_CONVDG(U1,U2,U3,&
    DK,NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,E031)
IF (IELT.EQ.1) CALL ST_CONVDG(U1,U2,U3,&
    DK,NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,E030)
IF (IELT.EQ.2) CALL ST_CONVNP(U1,U2,U3,&
    DK,NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,EM31)
IF (IELT.EQ.3) CALL ST_CONVNP(U1,U2,U3,&
    DK,NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,EM30)

END SUBROUTINE GetConvection_scalar
!
! ----------------------------------------------
!
SUBROUTINE GetAFCStuctures_Scalar(DK,KCOLA,KLDA,&
           NU,ISEP,IAUX,INOD,JNOD,AEDGE,IEDGE)

! Matrix assembly
REAL*8  DK(*)
INTEGER KCOLA(*),KLDA(*),NU
!AFC
INTEGER ISEP(*),IAUX(*),INOD(*),JNOD(*),IEDGE
REAL*8  AEDGE(*)

!---------------------------------------------------------------
REAL*8  DK_IJ,DK_JI,D_IJ
INTEGER II_LOC,IJ_LOC,JI_LOC,JJ_LOC,ILOC
INTEGER NA

EXTERNAL E030,E031,EM31,EM30

IEDGE = 0

DO I=1,NU
  ILOC=KLDA(I)
1 ILOC=ILOC+1 
  IF (KCOLA(ILOC).LT.I.AND.ILOC.LT.KLDA(I+1)) GOTO 1
  ISEP(I)=ILOC-1
END DO

CALL LCP3(ISEP,IAUX,NU)

DO I=1,NU
  II_LOC=KLDA(I)
  DO IJ_LOC=KLDA(I)+1,ISEP(I)
!   Position of entries in global matrices
    J=KCOLA(IJ_LOC)
    IAUX(J)=IAUX(J)+1
    JI_LOC=IAUX(J)
    JJ_LOC=KLDA(J)

!   Discrete diffusion coefficients
    DK_IJ=DK(IJ_LOC)
    DK_JI=DK(JI_LOC)
    D_IJ=MAX(-DK_IJ,0D0,-DK_JI)

!   Elimination of negative entries
    DK(II_LOC)=DK(II_LOC)-D_IJ
    DK(IJ_LOC)=DK(IJ_LOC)+D_IJ
    DK(JJ_LOC)=DK(JJ_LOC)-D_IJ
    DK(JI_LOC)=DK(JI_LOC)+D_IJ

!   Determine which node is located upwind 
    IEDGE=IEDGE+1
    IF (DK_JI.GT.DK_IJ) THEN
      INOD(IEDGE)=I
      JNOD(IEDGE)=J
      AEDGE(IEDGE)=MIN(D_IJ,DK_JI+D_IJ)
    ELSE
      INOD(IEDGE)=J
      JNOD(IEDGE)=I
      AEDGE(IEDGE)=MIN(D_IJ,DK_IJ+D_IJ)
    ENDIF

  ENDDO

ENDDO

END SUBROUTINE GetAFCStuctures_Scalar
!
! ----------------------------------------------
!
SUBROUTINE MatdefAFC_Scalar(U,D,A,DK,DS,KCOLA,KLDA,NU,DNUT,DML,&
           INOD,JNOD,AEDGE,IDEF,IMAT)

REAL*8  U(*),D(*),DK(*),DS(*),DNUT(*),DML(*)
REAL*4  A(*)
INTEGER KCOLA(*),KLDA(*),NU,IDEF,IMAT
!AFC
INTEGER INOD(*),JNOD(*)
REAL*8  AEDGE(*)

!---------------------------------------------------------------
REAL*8  DAUX,SIGML,DK_IJ,DK_JI,D_IJ
INTEGER II_LOC,IJ_LOC,JI_LOC,JJ_LOC,ILOC
INTEGER NA

! Add the contribution of the mass matrix to the defect
IF (IDEF.NE.0) THEN
 SIGML=DBLE(SIGN(1,IDEF))
 DO I=1,NU
  DAUX=SIGML*DML(I)
  D(I)=D(I)+DAUX*U(I)
 ENDDO
ENDIF

NA=KLDA(NU+1)-1

! Perform global matrix assembly
IF (IDEF.LE.0.AND.IMAT.EQ.1) THEN

 DO I=1,NU

  II_LOC=KLDA(I)
  DAUX=DK(II_LOC)!+DS(II_LOC)
  A(II_LOC)=REAL(DML(I))-REAL(THSTEP*DAUX)

  DO IJ_LOC=KLDA(I)+1,KLDA(I+1)-1
   DAUX=DK(IJ_LOC)!+DS(IJ_LOC)
   A(IJ_LOC)=-REAL(THSTEP*DAUX)
  ENDDO
 ENDDO

END IF

! Augment the defect if required
IF (IDEF.NE.0) THEN
 CALL DefTVD_Scalar(DK,KCOLA,KLDA,NU,U,D,&
 DS,INOD,JNOD,AEDGE)
ENDIF

END SUBROUTINE MatdefAFC_Scalar
!
! ----------------------------------------------
!
SUBROUTINE DefTVD_Scalar(DK,KCOLA,KLDA,NU,U,D,DS,&
           INOD,JNOD,AEDGE)

REAL*8 DEPS
REAL*8 DK(*),DS(*),U(*),D(*),AEDGE(*)
INTEGER INOD(*),JNOD(*),KCOLA(*),KLDA(*),NU

REAL*8 PP(NU),PM(NU),QP(NU),QM(NU)
PARAMETER (DEPS=1D-15)
INTEGER IEDGE,ILOC
REAL*8  DD,DAUX

DO I=1,NU

  PP(I)=0D0
  PM(I)=0D0
  QP(I)=0D0
  QM(I)=0D0

  DD=(DK(KLDA(I)))*U(I)!+DS(KLDA(I))
  DO ILOC=KLDA(I)+1,KLDA(I+1)-1
    J=KCOLA(ILOC)
    DD=DD+(DK(ILOC))*U(J)!+DS(ILOC)
  ENDDO
  D(I)=D(I)+THSTEP*DD
ENDDO

!RETURN

DO IEDGE=1,NEDGE

  I=INOD(IEDGE)
  J=JNOD(IEDGE)

  DAUX=AEDGE(IEDGE)*(U(I)-U(J))

  PP(I)=PP(I)+MAX(0D0,DAUX)
  PM(I)=PM(I)+MIN(0D0,DAUX)

  QP(I)=QP(I)+MAX(0D0,-DAUX)
  QM(I)=QM(I)+MIN(0D0,-DAUX)

  QP(J)=QP(J)+MAX(0D0,DAUX)
  QM(J)=QM(J)+MIN(0D0,DAUX)

ENDDO

CALL CommSum(PP,ILEV)          ! PARALLEL
CALL CommSum(PM,ILEV)          ! PARALLEL
CALL CommSum(QP,ILEV)          ! PARALLEL
CALL CommSum(QM,ILEV)          ! PARALLEL

DO I=1,NU
  IF (PP(I).GT.QP(I)+DEPS) THEN
    PP(I)=QP(I)/PP(I)
  ELSE
    PP(I)=1D0
  ENDIF
  IF (PM(I).LT.QM(I)-DEPS) THEN
    PM(I)=QM(I)/PM(I)
  ELSE
    PM(I)=1D0
  ENDIF
ENDDO

DO IEDGE=1,NEDGE

  I=INOD(IEDGE)
  J=JNOD(IEDGE)

  DAUX=THSTEP*AEDGE(IEDGE)*(U(I)-U(J))

  IF (DAUX.GT.0) THEN
    DAUX=PP(I)*DAUX
  ELSE
    DAUX=PM(I)*DAUX
  ENDIF

 D(I)=D(I)+DAUX
 D(J)=D(J)-DAUX

ENDDO

END SUBROUTINE DefTVD_Scalar
!
! ----------------------------------------------
!
SUBROUTINE Resdfk_General_Scalar(myScalar,resScalar,defScalar,rhsScalar)
TYPE(mg_scalar), INTENT(INOUT) :: myScalar
REAL*8  resScalar,defScalar,rhsScalar,RESF,RESU

CALL LL21 (myScalar%l(NLMAX)%rhs,myScalar%l(NLMAX)%ndof,RESF)
RESF=MAX(1D-15,RESF)

CALL LL21 (myScalar%l(NLMAX)%def,myScalar%l(NLMAX)%ndof,RESU)

resScalar = RESU/RESF
defScalar = RESU
rhsScalar = RESF
!write(*,*) RESU,RESF

END SUBROUTINE Resdfk_General_Scalar
!
! ----------------------------------------------
!
SUBROUTINE Boundary_General_Scalar_mat(myScalar)
TYPE(mg_scalar), INTENT(INOUT) :: myScalar

CALL Boundary_Scalar_mat(myScalar%l(ILEV)%Amat,&
KWORK(L(KLLDA(ILEV))),myScalar%l(ILEV)%knpr,myScalar%l(ILEV)%ndof)

END SUBROUTINE Boundary_General_Scalar_mat
!
! ----------------------------------------------
!
SUBROUTINE Boundary_Scalar_mat(VA,KLD,KNPR,NU)
REAL*4  VA(*)
INTEGER KLD(*),KNPR(*),NU,ICOL

DO I=1,NU
 IF (KNPR(I).eq.1) THEN
   VA(KLD(I))=1E0
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    VA(ICOL)=0E0
   END DO
 END IF
END DO

END SUBROUTINE Boundary_Scalar_mat
!
! ----------------------------------------------
!
SUBROUTINE RESTR_general_scalar(myScalar)
TYPE(mg_scalar), INTENT(INOUT) :: myScalar
INTEGER I1,I2

I1 = ILEV-1
I2 = ILEV
CALL RESTR_scalar(myScalar%l(I1)%val,myScalar%l(I2)%val, &
     KWORK(L(KLVERT(I1))),KWORK(L(KLAREA (I1))),         &
     KWORK(L(KLADJ(I1))),KNU(I1),KNP(I1),KNVT(I1),       &
     KWORK(L(KLVERT(I2))),KWORK(L(KLAREA (I2))),         &
     KWORK(L(KLADJ(I2))),KNU(I2),KNP(I2),KNVT(I2),       &
     VWORK(L(KLVOL(I1))),1)

END SUBROUTINE RESTR_general_scalar
!
! ----------------------------------------------
!
SUBROUTINE RESTR_scalar(DU1,DU2,              &
           KVERT1,KAREA1,KADJ1,NEQ1,NEL1,NVT1,&
           KVERT2,KAREA2,KADJ2,NEQ2,NEL2,NVT2,&
           AVOL,IVEL)

IMPLICIT DOUBLE PRECISION (A,D,R)
IMPLICIT INTEGER (I)

PARAMETER (A1=0.125D0,A2=0.25D0,A3=0.08333333D0)
PARAMETER (R1=0.25D0,R2=0.5D0,R3=0.16666666D0)
REAL*8 DU1(*),DU2(*)
INTEGER KVERT1(8,*),KAREA1(6,*),KADJ1(6,*)
INTEGER KVERT2(8,*),KAREA2(6,*),KADJ2(6,*)
INTEGER NEQ1,NEL1,NVT1,NEQ2,NEL2,NVT2,IVEL
REAL*4  AVOL(*)

      DO 10 IEL1=1,NEL1

      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)

      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(6,IELH2)
      IELH7=KADJ2(6,IELH3)
      IELH8=KADJ2(6,IELH4)

      I1=KAREA2(1,IELH1)
      I2=KAREA2(2,IELH1)
      I3=KAREA2(3,IELH1)
      I4=KAREA2(4,IELH1)
      I5=KAREA2(5,IELH1)
      I6=KAREA2(6,IELH1)
      I7=KAREA2(1,IELH2)
      I8=KAREA2(2,IELH2)
      I9=KAREA2(3,IELH2)
      I10=KAREA2(5,IELH2)
      I11=KAREA2(6,IELH2)
      I12=KAREA2(1,IELH3)
      I13=KAREA2(2,IELH3)
      I14=KAREA2(3,IELH3)
      I15=KAREA2(5,IELH3)
      I16=KAREA2(6,IELH3)
      I17=KAREA2(1,IELH4)
      I18=KAREA2(2,IELH4)
      I19=KAREA2(5,IELH4)
      I20=KAREA2(6,IELH4)
      I21=KAREA2(1,IELH5)
      I22=KAREA2(2,IELH5)
      I23=KAREA2(3,IELH5)
      I24=KAREA2(4,IELH5)
      I25=KAREA2(5,IELH5)
      I26=KAREA2(1,IELH6)
      I27=KAREA2(2,IELH6)
      I28=KAREA2(3,IELH6)
      I29=KAREA2(5,IELH6)
      I30=KAREA2(1,IELH7)
      I31=KAREA2(2,IELH7)
      I32=KAREA2(3,IELH7)
      I33=KAREA2(5,IELH7)
      I34=KAREA2(1,IELH8)
      I35=KAREA2(2,IELH8)
      I36=KAREA2(5,IELH8)

      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
      DUH13=DU2(I13)
      DUH14=DU2(I14)
      DUH15=DU2(I15)
      DUH16=DU2(I16)
      DUH17=DU2(I17)
      DUH18=DU2(I18)
      DUH19=DU2(I19)
      DUH20=DU2(I20)
      DUH21=DU2(I21)
      DUH22=DU2(I22)
      DUH23=DU2(I23)
      DUH24=DU2(I24)
      DUH25=DU2(I25)
      DUH26=DU2(I26)
      DUH27=DU2(I27)
      DUH28=DU2(I28)
      DUH29=DU2(I29)
      DUH30=DU2(I30)
      DUH31=DU2(I31)
      DUH32=DU2(I32)
      DUH33=DU2(I33)
      DUH34=DU2(I34)
      DUH35=DU2(I35)
      DUH36=DU2(I36)

      IF (IVEL.EQ.1) THEN
      AVOL11=DBLE(AVOL(IEL1))
      ENDIF

      IF (IVEL.EQ.0) THEN
      AVOL11=1D0
      ENDIF

! *** The area IA1

      IF (KADJ1(1,IEL1).NE.0) THEN
!     case of an inner area
       IF (KADJ1(1,IEL1).GT.IEL1) THEN
        DU1(IA1)=(A1*(DUH1+DUH7+DUH12+DUH17)+&
          A2*(DUH3+DUH14+DUH4+DUH9)-&
          A3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+&
              DUH5+DUH6+DUH11+&
              DUH16+DUH20))*AVOL11
      ELSE
        DU1(IA1)=DU1(IA1)+&
          (A1*(DUH1+DUH7+DUH12+DUH17)+A2*(DUH3+DUH14+DUH4+DUH9)-&
          A3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+&
              DUH5+DUH6+DUH11+&
              DUH16+DUH20))*AVOL11
      ENDIF
      ELSE
!     case of a boundary edge
        DU1(IA1)=R1*(DUH1+DUH7+DUH12+DUH17)+&
           R2*(DUH3+DUH14+DUH4+DUH9)-&
          R3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+&
              DUH5+DUH6+DUH11+&
              DUH16+DUH20)
      ENDIF

! *** The area IA2

      IF (KADJ1(2,IEL1).NE.0) THEN
!     case of an inner area
       IF (KADJ1(2,IEL1).GT.IEL1) THEN
        DU1(IA2)=(A1*(DUH2+DUH10+DUH29+DUH22)+&
           A2*(DUH3+DUH23+DUH6+DUH11)-&
           A3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+&
               DUH5+DUH4+DUH9+&
               DUH24+DUH28))*AVOL11  
      ELSE
        DU1(IA2)=DU1(IA2)+&
           (A1*(DUH2+DUH10+DUH29+DUH22)+A2*(DUH3+DUH23+DUH6+DUH11)-&
           A3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+&
               DUH5+DUH4+DUH9+&
               DUH24+DUH28))*AVOL11  
       ENDIF
       ELSE  
!     case of a boundary area
        DU1(IA2)=R1*(DUH2+DUH10+DUH29+DUH22)+&
           R2*(DUH3+DUH23+DUH6+DUH11)-&
           R3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+&
               DUH5+DUH4+DUH9+&
               DUH24+DUH28)  
      ENDIF

! *** The edge IA3

      IF (KADJ1(3,IEL1).NE.0) THEN
!     case of an inner area
       IF (KADJ1(3,IEL1).GT.IEL1) THEN
        DU1(IA3)=(A1*(DUH8+DUH15+DUH33+DUH27)+&
           A2*(DUH11+DUH16+DUH9+DUH28)-&
           A3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+&
               DUH31+DUH30+DUH26+&
                DUH29+DUH10))*AVOL11  
      ELSE
        DU1(IA3)= DU1(IA3)+&
           (A1*(DUH8+DUH15+DUH33+DUH27)+A2*(DUH11+DUH16+DUH9+DUH28)-&
           A3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+&
               DUH31+DUH30+DUH26+&
                DUH29+DUH10))*AVOL11
        ENDIF
        ELSE  
!     case of a boundary area 
        DU1(IA3)=R1*(DUH8+DUH15+DUH33+DUH27)+&
           R2*(DUH11+DUH16+DUH9+DUH28)-&
           R3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+&
               DUH31+DUH30+DUH26+&
                DUH29+DUH10)  
      ENDIF

! *** The area IA4

      IF (KADJ1(4,IEL1).NE.0) THEN 
!     case of an inner area
       IF (KADJ1(4,IEL1).GT.IEL1) THEN
        DU1(IA4)=(A1*(DUH13+DUH31+DUH36+DUH19)+&
         A2*(DUH32+DUH14+DUH16+DUH20)-&
         A3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+&
             DUH18+DUH4+DUH9+&
             DUH28+DUH24))*AVOL11
      ELSE
        DU1(IA4)=DU1(IA4)+&
         (A1*(DUH13+DUH31+DUH36+DUH19)+A2*(DUH32+DUH14+DUH16+DUH20)-&
         A3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+&
             DUH18+DUH4+DUH9+&
             DUH28+DUH24))*AVOL11  
        ENDIF
        ELSE  
!     case of a boundary area
        DU1(IA4)=R1*(DUH13+DUH31+DUH36+DUH19)+&
         R2*(DUH32+DUH14+DUH16+DUH20)-&
         R3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+&
             DUH18+DUH4+DUH9+&
             DUH28+DUH24)  
      ENDIF

! *** The area IA5

      IF (KADJ1(5,IEL1).NE.0) THEN 
!     case of an inner area
       IF (KADJ1(5,IEL1).GT.IEL1) THEN
        DU1(IA5)=(A1*(DUH5+DUH18+DUH35+DUH25)+&
          A2*(DUH6+DUH20+DUH4+DUH24)-&
          A3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+&
              DUH2+DUH3+DUH14+&
              DUH32+DUH23))*AVOL11  
      ELSE
        DU1(IA5)=DU1(IA5)+&
          (A1*(DUH5+DUH18+DUH35+DUH25)+A2*(DUH6+DUH20+DUH4+DUH24)-&
          A3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+&
              DUH2+DUH3+DUH14+&
              DUH32+DUH23))*AVOL11  
       ENDIF
       ELSE  
!     case of a boundary area
        DU1(IA5)=R1*(DUH5+DUH18+DUH35+DUH25)+&
          R2*(DUH6+DUH20+DUH4+DUH24)-&
          R3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+&
              DUH2+DUH3+DUH14+&
              DUH32+DUH23)  
      ENDIF

! *** The area IA6

      IF (KADJ1(6,IEL1).NE.0) THEN 
!     case of an inner area
       IF (KADJ1(6,IEL1).GT.IEL1) THEN
        DU1(IA6)=(A1*(DUH26+DUH30+DUH34+DUH21)+&
          A2*(DUH23+DUH32+DUH24+DUH28)-&
          A3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+&
              DUH25+DUH6+DUH11+&
              DUH16+DUH20))*AVOL11  
      ELSE
        DU1(IA6)= DU1(IA6)+&
          (A1*(DUH26+DUH30+DUH34+DUH21)+A2*(DUH23+DUH32+DUH24+DUH28)-&
          A3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+&
              DUH25+DUH6+DUH11+&
              DUH16+DUH20))*AVOL11  
        ENDIF
        ELSE 
!     case of a boundary area
        DU1(IA6)=R1*(DUH26+DUH30+DUH34+DUH21)+&
          R2*(DUH23+DUH32+DUH24+DUH28)-&
          R3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+&
              DUH25+DUH6+DUH11+&
              DUH16+DUH20)  
      ENDIF

10    CONTINUE

      IF (IVEL.EQ.1) THEN
      DO 555 IEL1=1,NEL1
      DO 666 I=1,6
      IADJ=KADJ1(I,IEL1)
      IA=KAREA1(I,IEL1)
      IF (IADJ.EQ.0) THEN
      GOTO 666
      ENDIF
      IF (IADJ.GT.IEL1) THEN
      AVOL1=DBLE(AVOL(IEL1))
      AVOL2=DBLE(AVOL(IADJ))
      DU1(IA)=(1D0/SQRT(AVOL1+AVOL2))*DU1(IA)
      ELSE
      AVOL1=DBLE(AVOL(IEL1))
      AVOL2=DBLE(AVOL(IADJ))
      DU1(IA)=(2D0/SQRT(AVOL1+AVOL2))*DU1(IA)
      ENDIF
666   CONTINUE
555   CONTINUE
      ENDIF

END SUBROUTINE RESTR_scalar
!
! ----------------------------------------------
!
SUBROUTINE Protocol_General_Scalar(mfile,myScalar,nINL,&
           ResScalar,DefScalar,RhsScalar)
TYPE(mg_scalar), INTENT(INOUT) :: myScalar
INTEGER nINL,mfile
INTEGER length
REAL*8 ResScalar,DefScalar,RhsScalar
CHARACTER C1*14,C2*14,C3*14


IF (myid.eq.showID) THEN
length =  LEN(myScalar%cName)

C1='              '
C2='              '
C3='              '
WRITE(C1(1:3+length),'(A3,A7)') 'Res',myScalar%cName
WRITE(C2(1:3+length),'(A3,A7)') 'Def',myScalar%cName
WRITE(C3(1:7+length),'(A7,A7)') 'GlobDef',myScalar%cName

IF (nINL.EQ.0) THEN
 WRITE(MTERM,5)
 WRITE(MFILE,5)
 WRITE(MTERM,'(A8,5(2X,A14))') "INL",TRIM(C1),TRIM(C2),TRIM(C3)
 WRITE(MFILE,'(A8,5(2X,A14))') "INL",TRIM(C1),TRIM(C2),TRIM(C3)
 WRITE(MTERM,5)
 WRITE(MFILE,5)
 WRITE(MTERM,'(A8,6XA10,5(6X,D10.4))') "Criteria"," ",DefScalar*myScalar%defCrit,RhsScalar
 WRITE(MFILE,'(A8,6XA10,5(6X,D10.4))') "Criteria"," ",DefScalar*myScalar%defCrit,RhsScalar
 WRITE(MTERM,5)
 WRITE(MFILE,5)
 WRITE(MTERM,'(I8,5(6X,D10.4))') 0,ResScalar,DefScalar
 WRITE(MFILE,'(I8,5(6X,D10.4))') 0,ResScalar,DefScalar
ELSE
 WRITE(MTERM,'(I8,5(6X,D10.4))') nINL,ResScalar,DefScalar,RhsScalar
 WRITE(MFILE,'(I8,5(6X,D10.4))') nINL,ResScalar,DefScalar,RhsScalar
END IF

END IF

5  FORMAT(80('-'))

END SUBROUTINE Protocol_General_Scalar
!
! ----------------------------------------------
!
SUBROUTINE Solve_General_Scalar(myScalar,Bndry_Val)
TYPE(mg_scalar), INTENT(INOUT) :: myScalar
EXTERNAL Bndry_Val

ILEV=NLMAX
CALL SETLEV(2)
CALL Boundary_General_Scalar_mat(myScalar)
CALL CommMat(myScalar%l(ILEV)%Amat,&                   ! PARALLEL
     KWORK(L(KLLDA(ILEV))),KNU(ILEV),ILEV)                ! PARALLEL

CALL LCL1 (myScalar%l(ILEV)%val,myScalar%l(ILEV)%ndof)
CALL ParIA237(myScalar%l(ILEV)%Amat,KWORK(L(KLCOLA(ILEV))),&
     KWORK(L(KLLDA(ILEV))),myScalar%l(ILEV)%val,&
     myScalar%l(ILEV)%def,myScalar%l(ILEV)%aux,&
     myScalar%l(ILEV)%ndof,myScalar%SolvIter,RLXSMU,ILEV)

! Update the solution
CALL LLC1(myScalar%l(ILEV)%val_old,myScalar%l(ILEV)%val,&
     myScalar%l(ILEV)%ndof,1D0,1D0)

! Set dirichlet boundary conditions on the solution
CALL Bndry_Val(DWORK(L(KLCAG(ILEV))))

END SUBROUTINE Solve_General_Scalar

END MODULE

