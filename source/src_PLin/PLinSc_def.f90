MODULE def_PLinScalar

USE def_FEAT
USE PP3D_MPI, ONLY:myid,showID,E012_CommSetUp,E012_CommVal
use types

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



TYPE TMatrix
 INTEGER :: nu,na
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: ColA,LdA
END TYPE
TYPE(TMatrix) :: lMat,lMatQ0

REAL*8  , DIMENSION(:,:,:)  , ALLOCATABLE :: mMat,imMat,kMat
REAL*8  , DIMENSION(:,:,:)  , ALLOCATABLE :: fNorm
REAL*8  , DIMENSION(:,:)    , ALLOCATABLE :: dNorm
REAL*8  , DIMENSION(:)      , ALLOCATABLE :: midElem
REAL*8  , DIMENSION(:)      , ALLOCATABLE :: FracFieldQ0,FracFieldQ1
INTEGER , DIMENSION(:)      , ALLOCATABLE :: IntPhaseElem
REAL*8,   DIMENSION(:,:)    , ALLOCATABLE :: ElemMid

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE InitCommP1(myPLinSc)
TYPE(lScalar) myPLinSc
INTEGER I,J
EXTERNAL E012

ILEV=NLMAX
CALL SETLEV(2)

CALL E012_CommSetUp(1,myPLinSc%iParFace,myPLinSc%iNumFace)
ALLOCATE (myPLinSc%iParFace(5,myPLinSc%iNumFace))
myPLinSc%iParFace=0d0
CALL E012_CommSetUp(2,myPLinSc%iParFace,myPLinSc%iNumFace)
ALLOCATE (myPLinSc%dParMidC(3,myPLinSc%iNumFace))
myPLinSc%dParMidC = 0d0
ALLOCATE (myPLinSc%dParFace(4,myPLinSc%iNumFace))
myPLinSc%dParFace = 0d0

CALL GetParMidValues(DWORK(L(LCORVG)),myPLinSc%dParFace,&
     KWORK(L(LAREA)),KWORK(L(LVERT)),myPLinSc%iParFace,NEL)

CALL E012_CommVal(myPLinSc%dParFace,myPLinSc%iParFace)

myPLinSc%dParMidC(1,:) = myPLinSc%dParFace(1,:)
myPLinSc%dParMidC(2,:) = myPLinSc%dParFace(2,:)
myPLinSc%dParMidC(3,:) = myPLinSc%dParFace(3,:)

END SUBROUTINE InitCommP1
!
! ----------------------------------------------
!
SUBROUTINE PLinScP1toQ1(myPLinSc)
TYPE(lScalar) myPLinSc
EXTERNAL E012

ILEV=NLMAX
CALL SETLEV(2)

myPLinSc%Q1=0d0
myPLinSc%aux=0d0
CALL IntP1toQ1(myPLinSc%Q1,myPLinSc%val(NLMAX)%x,&
     myPLinSc%aux,KWORK(L(LVERT)),DWORK(L(LCORVG)),&
     VWORK(L(KLVOL(ILEV))),E012)

END SUBROUTINE PLinScP1toQ1
!
! ----------------------------------------------
!
SUBROUTINE PLinSc_Solve(iEquation,myPLinSc)
TYPE(lScalar) myPLinSc
INTEGER I,J,K,iEquation

ILEV=NLMAX
CALL SETLEV(2)

myPLinSc%val_old(:) = myPLinSc%val(NLMAX)%x(:)
myPLinSc%val(NLMAX)%x(:) = 0d0
IF     (iEquation.EQ.1) THEN
 CALL imMatMultF1(myPLinSc%val(NLMAX)%x,myPLinSc%def,&
      IntPhaseElem,ElemMid,imMat,NEL)
ELSEIF (iEquation.EQ.2) THEN
 CALL imMatMultF2(myPLinSc%val(NLMAX)%x,myPLinSc%def,&
      IntPhaseElem,imMat,NEL)
ELSE
 CALL imMatMultF3(myPLinSc%val(NLMAX)%x,myPLinSc%def,&
      IntPhaseElem,imMat,NEL)
END IF

END SUBROUTINE PLinSc_Solve
!
! ----------------------------------------------
!
SUBROUTINE PLinSc_RHS(iEquation,myPLinSc,U,V,W,dt,idef)
REAL*8 U(*),V(*),W(*),dt,dtimelevel
TYPE(lScalar) myPLinSc
INTEGER idef
INTEGER I,J,K,iEquation
EXTERNAL E011,E012,E013,GetBCVal

ILEV=NLMAX
CALL SETLEV(2)

CALL GetParValues(myPLinSc%val(NLMAX)%x,myPLinSc%dParFace,&
     KWORK(L(LAREA)),myPLinSc%iParFace,NEL)

CALL E012_CommVal(myPLinSc%dParFace,myPLinSc%iParFace)

IF (idef.eq. 1) dtimelevel = timens-tstep
IF (idef.eq.-1) dtimelevel = timens

IF (iEquation.EQ.1) THEN
 CALL Build_FluxQP1(U,V,W,myPLinSc%val(NLMAX)%x,fNorm,ElemMid,&
      myPLinSc%def,IntPhaseElem,myPLinSc%dParFace,&
      myPLinSc%dParMidC,myPLinSc%iParFace,KWORK(L(LVERT)),&
      KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),&
      DWORK(L(LCORVG)),KWORK(L(LADJ)),myPLinSc%prm%SrfCubF,&
      myPLinSc%prm%StabScheme,E012,E013,GetBCVal,dtimelevel,dt,1)

 CALL BuildRHS1(myPLinSc%def,IntPhaseElem,ElemMid,kMat,mMat,&
      myPLinSc%val(NLMAX)%x,NEL,dt,idef)
ELSE
 CALL Build_FluxQP1(U,V,W,myPLinSc%val(NLMAX)%x,fNorm,ElemMid,&
      myPLinSc%def,IntPhaseElem,myPLinSc%dParFace,&
      myPLinSc%dParMidC,myPLinSc%iParFace,KWORK(L(LVERT)),&
      KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),&
      DWORK(L(LCORVG)),KWORK(L(LADJ)),myPLinSc%prm%SrfCubF,&
      myPLinSc%prm%StabScheme,E012,E011,GetBCVal,dtimelevel,dt,2)

 CALL BuildRHS2(myPLinSc%def,IntPhaseElem,kMat,mMat,&
      myPLinSc%val(NLMAX)%x,NEL,dt,idef)
END IF


END SUBROUTINE PLinSc_RHS
!
! ----------------------------------------------
!
SUBROUTINE InitFNormMatrix()
INTEGER I,J,K


 ILEV=NLMAX
 CALL SETLEV(2)
 ALLOCATE( fNorm(6,3,NEL))
 CALL Build_NormQP1(fNorm,KWORK(L(LVERT)),DWORK(L(LCORVG)),NEL)

END SUBROUTINE InitFNormMatrix
!
! ----------------------------------------------
!
SUBROUTINE InitConvMatrix(iEquation,myPLinSc,U,V,W)
TYPE(lScalar) myPLinSc
REAL*8  U(*),V(*),W(*)
INTEGER I,J,K,iEquation
EXTERNAL E011,E012,E013

 ILEV=NLMAX
 CALL SETLEV(2)
 IF (.NOT.ALLOCATED(kMat)) THEN
  ALLOCATE(kMat(4,4,NEL))
 END IF
 kMat = 0d0
 IF (iEquation.EQ.1) THEN
  CALL Build_ConvQP1(U,V,W,kMat,IntPhaseElem,KWORK(L(LVERT)),&
       KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(NLMAX))),&
       DWORK(L(LCORVG)),myPLinSc%prm%VolCubF,E012,E013,1)
 ELSE
  CALL Build_ConvQP1(U,V,W,kMat,IntPhaseElem,KWORK(L(LVERT)),&
       KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(NLMAX))),&
       DWORK(L(LCORVG)),myPLinSc%prm%VolCubF,E012,E011,2)
 END IF

! if (myid.eq.1) OPEN(987,FILE='KMAT1.txt')
! if (myid.eq.2) OPEN(987,FILE='KMAT2.txt')
! if (myid.eq.3) OPEN(987,FILE='KMAT3.txt')
! 
! DO I=1,NEL
! WRITE(987,*) "---",I
!  DO J=1,4
!   WRITE(987,*) (kMat(K,J,I),K=1,4)
!  END DO
! END DO
! CLOSE(987)


END SUBROUTINE InitConvMatrix
!
! ----------------------------------------------
!
SUBROUTINE RISource(myPLinSc,dt)
TYPE(lScalar) myPLinSc
REAL*8 dt
EXTERNAL E012

CALL Build_RISrcQP1(myPLinSc%def,IntPhaseElem,KWORK(L(LVERT)),&
     KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),&
     DWORK(L(LCORVG)),myPLinSc%prm%VolCubF,E012,dt)

END SUBROUTINE RISource
!
! ----------------------------------------------
!
SUBROUTINE InitMassMatrix(myPLinSc)
TYPE(lScalar) myPLinSc
INTEGER I,J,K
EXTERNAL E012
integer :: errorflag = 0


 ILEV=NLMAX
 CALL SETLEV(2)
 ALLOCATE( mMat(4,4,NEL))
 mMat=0d0
 CALL Build_MassQP1(mMat,KWORK(L(LVERT)),KWORK(L(LAREA)),&
      KWORK(L(LEDGE)),KWORK(L(KLINT(NLMAX))),DWORK(L(LCORVG)),&
      myPLinSc%prm%VolCubF,E012)

 ALLOCATE(imMat(4,4,NEL))

 DO I=1,NEL
  CALL FINDInv(mMat(:,:,I), imMat(:,:,I), 4, errorflag)
 END DO

END SUBROUTINE InitMassMatrix
!
! ----------------------------------------------
!
SUBROUTINE Create_MatStruct()
INTEGER iSymm,nERow,nMat
INTEGER , DIMENSION(:)  , ALLOCATABLE :: TempColA
INTEGER I,J

ILEV=NLMAX
CALL SETLEV(2)

END SUBROUTINE Create_MatStruct
!
! ----------------------------------------------
!
SUBROUTINE InitField(myPLinSc)
TYPE(lScalar) myPLinSc

ILEV=NLMAX
CALL SETLEV(2)

! Initialization of the levelset field
myPLinSc%ndof = 4*NEL
myPLinSc%nel = NEL
ALLOCATE(myPLinSc%val(NLMIN:NLMAX))
ALLOCATE(myPLinSc%val(NLMAX)%x(myPLinSc%ndof))
ALLOCATE(myPLinSc%def(myPLinSc%ndof))
ALLOCATE(myPLinSc%rhs(myPLinSc%ndof))
ALLOCATE(myPLinSc%aux(myPLinSc%ndof))
ALLOCATE(myPLinSc%val_old(myPLinSc%ndof))
ALLOCATE(myPLinSc%val_ad(myPLinSc%ndof))
ALLOCATE(myPLinSc%Q1(NVT))


END SUBROUTINE InitField
!
! ----------------------------------------------
!
END MODULE def_PLinScalar