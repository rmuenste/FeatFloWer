SUBROUTINE Configure_AFCStruct(KCOLA,KLDA,NU,ISEP)
USE def_LinScalar, ONLY : AFC
IMPLICIT NONE
INTEGER KCOLA(*),KLDA(*),NU
INTEGER ISEP(*)
!---------------------------------------------------------------
INTEGER II_LOC,IJ_LOC,JI_LOC,JJ_LOC,ILOC
INTEGER IEDGE,JEDGE
INTEGER I,J,iSEP_OLD

IEDGE = 0
JEDGE = 0

DO I=1,NU
 ISEP_OLD = ISEP(I)
 II_LOC=KLDA(I)
 J=KCOLA(II_LOC)
 ISEP(I) = II_LOC
 
 DO IJ_LOC=KLDA(I)+1,KLDA(I+1)-1
  J=KCOLA(IJ_LOC)
  IF (J.gt.I) IEDGE = IEDGE + 1
  IF (J.lt.I) JEDGE = JEDGE + 1
  
  IF (J.lt.I) ISEP(I)=IJ_LOC
  
 END DO
 
END DO

IF ((iedge.ne.AFC%nedge).or.(jedge .ne. AFC%nedge)) THEN
  write(*,*) "Ultimate problem in the AFC stab constructor "
  STOP
END IF

END SUBROUTINE Configure_AFCStruct
!
! ----------------------------------------------
!
SUBROUTINE AFC_LinScalar(DK,KCOLA,KLDA,&
           NU,ISEP,IAUX,INOD,JNOD,AEDGE)
USE def_LinScalar, ONLY : AFC
IMPLICIT NONE
! Matrix assembly
REAL*8  DK(*)
INTEGER KCOLA(*),KLDA(*),NU
!AFC
INTEGER ISEP(*),IAUX(*),INOD(*),JNOD(*),IEDGE
REAL*8  AEDGE(*)
!---------------------------------------------------------------
REAL*8  DK_IJ,DK_JI,D_IJ
INTEGER II_LOC,IJ_LOC,JI_LOC,JJ_LOC,ILOC
INTEGER I,J

CALL LCP3(ISEP,IAUX,NU)

IEDGE = 0

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

IF (iedge .ne. AFC%nedge) THEN
  write(*,*) "Ultimate problem in the AFC stab constructor "
  STOP
END IF

END
!
! ----------------------------------------------
!
SUBROUTINE DefTVD_linScalar(U,D,dt)

USE PP3D_MPI, ONLY:E011Sum,myid
USE def_LinScalar, ONLY : AFC
IMPLICIT NONE

REAL*8 dt
REAL*8 U(*),D(*)

REAL*8 DEPS,DAUX,F_IJ
PARAMETER (DEPS=1D-15)
INTEGER IEDGE,ILOC,I,J

CALL LCL1 (AFC%PP,AFC%nu)
CALL LCL1 (AFC%PM,AFC%nu)
CALL LCL1 (AFC%QP,AFC%nu)
CALL LCL1 (AFC%QM,AFC%nu)

DO IEDGE=1,AFC%nedge

  I=AFC%inod(IEDGE)
  J=AFC%jnod(IEDGE)

  DAUX=AFC%aedge(IEDGE)*(U(I)-U(J))

  AFC%PP(I)=AFC%PP(I)+MAX(0D0,DAUX)
  AFC%PM(I)=AFC%PM(I)+MIN(0D0,DAUX)

  AFC%QP(I)=AFC%QP(I)+MAX(0D0,-DAUX)
  AFC%QM(I)=AFC%QM(I)+MIN(0D0,-DAUX)

  AFC%QP(J)=AFC%QP(J)+MAX(0D0,DAUX)
  AFC%QM(J)=AFC%QM(J)+MIN(0D0,DAUX)

ENDDO

CALL E011Sum(AFC%PP)          ! PARALLEL
CALL E011Sum(AFC%PM)          ! PARALLEL
CALL E011Sum(AFC%QP)          ! PARALLEL
CALL E011Sum(AFC%QM)          ! PARALLEL

DO I=1,AFC%nu
  IF (AFC%PP(I).GT.AFC%QP(I)+DEPS) THEN
    AFC%PP(I)=AFC%QP(I)/AFC%PP(I)
  ELSE
    AFC%PP(I)=1D0
  ENDIF
  IF (AFC%PM(I).LT.AFC%QM(I)-DEPS) THEN
    AFC%PM(I)=AFC%QM(I)/AFC%PM(I)
  ELSE
    AFC%PM(I)=1D0
  ENDIF
ENDDO

DO IEDGE=1,AFC%nedge

  I=AFC%inod(IEDGE)
  J=AFC%jnod(IEDGE)

  F_IJ=dt*AFC%aedge(IEDGE)*(U(I)-U(J))

  IF (F_IJ.GT.0d0) THEN
    DAUX=AFC%PP(I)*F_IJ
  ELSE
    DAUX=AFC%PM(I)*F_IJ
  ENDIF
  
  D(I)=D(I)+DAUX
  D(J)=D(J)-DAUX

!  FLUX(IEDGE) = F_IJ - DAUX
!  IF (J.GT.I) FLUX(IEDGE) = -FLUX(IEDGE)

ENDDO

END
