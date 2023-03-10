MODULE def_cc

USE PP3D_MPI, ONLY:E011Sum,E011DMat,myid,showID,MGE013,&
                   COMM_Maximum,COMM_SUMM,COMM_NLComplete,SORT2D
USE var_QuadScalar_newton
USE def_QuadScalar, only:outputmatrix


USE UMFPackSolver_CC, ONLY : myUmfPack_CCFree, myUmfPack_CCFactorizeLocalMat

use mg_cc, only: MG_Solver_cc
USE mg_QuadScalar, ONLY : mgProlRestInit

IMPLICIT NONE

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE Matdef_general_QuadScalar_cc(myScalar,idef,alpha)
EXTERNAL E013
INTEGER :: idef
INTEGER I,J
TYPE(TQuadScalar) myScalar
REAL*8 daux,tttx1,tttx0,alpha




! Build up the matrix
 IF (idef.eq.-1) THEN
  DO ILEV=NLMIN,NLMAX

    IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
     A11Mat     => mg_A11Mat(ILEV)%a
     A22Mat     => mg_A22Mat(ILEV)%a
     A33Mat     => mg_A33Mat(ILEV)%a
     A12Mat     => mg_A12Mat(ILEV)%a
     A13Mat     => mg_A13Mat(ILEV)%a
     A23Mat     => mg_A23Mat(ILEV)%a
     A21Mat     => mg_A21Mat(ILEV)%a
     A31Mat     => mg_A31Mat(ILEV)%a
     A32Mat     => mg_A32Mat(ILEV)%a
     AA11Mat     => mg_AA11Mat(ILEV)%a
     AA22Mat     => mg_AA22Mat(ILEV)%a
     AA33Mat     => mg_AA33Mat(ILEV)%a
     AA12Mat     => mg_AA12Mat(ILEV)%a
     AA13Mat     => mg_AA13Mat(ILEV)%a
     AA23Mat     => mg_AA23Mat(ILEV)%a
     AA21Mat     => mg_AA21Mat(ILEV)%a
     AA31Mat     => mg_AA31Mat(ILEV)%a
     AA32Mat     => mg_AA32Mat(ILEV)%a
    ELSE
     A11Mat     => mg_A11Mat(ILEV)%a
     A22Mat     => mg_A22Mat(ILEV)%a
     A33Mat     => mg_A33Mat(ILEV)%a
     AA11Mat     => mg_AA11Mat(ILEV)%a
     AA22Mat     => mg_AA22Mat(ILEV)%a
     AA33Mat     => mg_AA33Mat(ILEV)%a
    END IF

    IF (myMatrixRenewal%S.GE.1) THEN
     S11Mat   => mg_S11Mat(ILEV)%a
     S22Mat   => mg_S22Mat(ILEV)%a
     S33Mat   => mg_S33Mat(ILEV)%a
     S12Mat   => mg_S12Mat(ILEV)%a
     S13Mat   => mg_S13Mat(ILEV)%a
     S23Mat   => mg_S23Mat(ILEV)%a
     S21Mat   => mg_S21Mat(ILEV)%a
     S31Mat   => mg_S31Mat(ILEV)%a
     S32Mat   => mg_S32Mat(ILEV)%a
    END IF

    IF (myMatrixRenewal%D.GE.1) THEN
     DMat     => mg_DMat(ILEV)%a
    END IF

    IF (myMatrixRenewal%K.GE.1) THEN
     KMat     => mg_KMat(ILEV)%a
     IF(bNonNewtonian) THEN
      barM11Mat     => mg_barM11Mat(ILEV)%a
      barM22Mat     => mg_barM22Mat(ILEV)%a
      barM33Mat     => mg_barM33Mat(ILEV)%a
      barM12Mat     => mg_barM12Mat(ILEV)%a
      barM13Mat     => mg_barM13Mat(ILEV)%a
      barM23Mat     => mg_barM23Mat(ILEV)%a
      barM21Mat     => mg_barM21Mat(ILEV)%a
      barM31Mat     => mg_barM31Mat(ILEV)%a
      barM32Mat     => mg_barM32Mat(ILEV)%a
     ELSE
      barM11Mat     => mg_barM11Mat(ILEV)%a
      barM22Mat     => mg_barM22Mat(ILEV)%a
      barM33Mat     => mg_barM33Mat(ILEV)%a
     END IF
    END IF

    IF (myMatrixRenewal%M.GE.1) THEN
     MMat     => mg_MMat(ILEV)%a
     MlRhoMat => mg_MlRhoMat(ILEV)%a
    END IF

    qMat     => mg_qMat(ILEV)
    !!-------------------    POINTER Setup  -------------------!!

    !!-------------------  MATRIX Assembly -------------------!!

    IF (bNonNewtonian) THEN ! Non-Newtonian
     IF (myMatrixRenewal%S.EQ.0) THEN
!      IF (myMatrixRenewal%K.GE.1) THEN ! Non-Newtonian Navier-Stokes with D
!       DO I=1,qMat%nu
!        DO J=qMat%LdA(I),qMat%LdA(I+1)-1
!         daux = MMat(J) + zeitstep*(2d0*DMat(J)+KMat(J))! 
!         daux2 = MMat(J) + zeitstep*(2d0*DMat(J)+KMat(J)+alpha*KnewMat(J))! 
!         A11mat(J) =  daux
!         A22mat(J) =  daux
!         A33mat(J) =  daux
!         AA11mat(J) =  daux2
!         AA22mat(J) =  daux2
!         AA33mat(J) =  daux2
!        END DO
!       END DO
!      ELSE !Non-Newtonian Stokes with D
!       DO I=1,qMat%nu
!        DO J=qMat%LdA(I),qMat%LdA(I+1)-1
!         daux =  MMat(J) + zeitstep*(2d0*DMat(J))!
!         A11mat(J) =  daux
!         A22mat(J) =  daux
!         A33mat(J) =  daux
!        END DO
!       END DO 
!      END IF     
     ELSE
      IF (myMatrixRenewal%K.GE.1) THEN ! Non-Newtonian Navier-Stokes with S
       DO I=1,qMat%nu
        DO J=qMat%LdA(I),qMat%LdA(I+1)-1
         ! non-linear defect
         A11mat(J) = MMat(J) + zeitstep*(S11Mat(J) + KMat(J))
         A22mat(J) = MMat(J) + zeitstep*(S22Mat(J) + KMat(J))
         A33mat(J) = MMat(J) + zeitstep*(S33Mat(J) + KMat(J))
         A12Mat(J) = zeitstep*S12Mat(J)
         A13Mat(J) = zeitstep*S13Mat(J)
         A23Mat(J) = zeitstep*S23Mat(J)
         A21Mat(J) = zeitstep*S21Mat(J)
         A31Mat(J) = zeitstep*S31Mat(J)
         A32Mat(J) = zeitstep*S32Mat(J)
         ! for multigrid
         AA11mat(J) = MMat(J) + zeitstep*(S11Mat(J) + KMat(J) + alpha * barM11Mat(J))
         AA22mat(J) = MMat(J) + zeitstep*(S22Mat(J) + KMat(J) + alpha * barM22Mat(J))
         AA33mat(J) = MMat(J) + zeitstep*(S33Mat(J) + KMat(J) + alpha * barM33Mat(J))
         AA12Mat(J) = zeitstep*(S12Mat(J) + alpha * barM12Mat(J))
         AA13Mat(J) = zeitstep*(S13Mat(J) + alpha * barM13Mat(J))
         AA23Mat(J) = zeitstep*(S23Mat(J) + alpha * barM23Mat(J))
         AA21Mat(J) = zeitstep*(S21Mat(J) + alpha * barM21Mat(J))
         AA31Mat(J) = zeitstep*(S31Mat(J) + alpha * barM31Mat(J))
         AA32Mat(J) = zeitstep*(S32Mat(J) + alpha * barM32Mat(J))
        END DO
       END DO
      ELSE ! Non-Newtonian Stokes with S
       DO I=1,qMat%nu
        DO J=qMat%LdA(I),qMat%LdA(I+1)-1
         A11mat(J) = MMat(J) + zeitstep*S11Mat(J)
         A22mat(J) = MMat(J) + zeitstep*S22Mat(J)
         A33mat(J) = MMat(J) + zeitstep*S33Mat(J)
         A12Mat(J) = zeitstep*S12Mat(J)
         A13Mat(J) = zeitstep*S13Mat(J)
         A23Mat(J) = zeitstep*S23Mat(J)
         A21Mat(J) = zeitstep*S21Mat(J)
         A31Mat(J) = zeitstep*S31Mat(J)
         A32Mat(J) = zeitstep*S32Mat(J)
         AA11mat(J) = A11mat(J)
         AA22mat(J) = A22mat(J)
         AA33mat(J) = A33mat(J)
         AA12Mat(J) = A12mat(J)
         AA13Mat(J) = A13mat(J)
         AA23Mat(J) = A23mat(J)
         AA21Mat(J) = A21mat(J)
         AA31Mat(J) = A31mat(J)
         AA32Mat(J) = A32mat(J)
        END DO
       END DO
      END IF
     END IF

    ELSE 

! Newtonian Navier-Stokes
!  ---> can be done with Non-newtonian problem

!     IF (myMatrixRenewal%K.GE.1) THEN
!      DO I=1,qMat%nu
!       DO J=qMat%LdA(I),qMat%LdA(I+1)-1
!        daux = MMat(J) + zeitstep*(DMat(J)+KMat(J))
!        daux2 = MMat(J) + zeitstep*(DMat(J)+KMat(J)+alpha*KnewMat(J))
!        A11mat(J) =  daux
!        A22mat(J) =  daux
!        A33mat(J) =  daux
!        AA11mat(J) =  daux2
!        AA22mat(J) =  daux2
!        AA33mat(J) =  daux2
!       END DO
!      END DO
!     ELSE ! Newtonian Stokes
!      DO I=1,qMat%nu
!       J = qMat%LdA(I)
!       DO J=qMat%LdA(I),qMat%LdA(I+1)-1
!        daux = MMat(J) + zeitstep*(DMat(J))
!        A11mat(J) =  daux
!        A22mat(J) =  daux
!        A33mat(J) =  daux
!       END DO
!      END DO
!     END IF
    END IF

  END DO
 END IF


    !!-------------------  MATRIX Assembly -------------------!!

    !!-------------------    POINTER Setup  -------------------!!
 ILEV=NLMAX
 CALL SETLEV(2)

 IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
  A11Mat     => mg_A11Mat(ILEV)%a
  A22Mat     => mg_A22Mat(ILEV)%a
  A33Mat     => mg_A33Mat(ILEV)%a
  A12Mat     => mg_A12Mat(ILEV)%a
  A13Mat     => mg_A13Mat(ILEV)%a
  A23Mat     => mg_A23Mat(ILEV)%a
  A21Mat     => mg_A21Mat(ILEV)%a
  A31Mat     => mg_A31Mat(ILEV)%a
  A32Mat     => mg_A32Mat(ILEV)%a
  AA11Mat     => mg_AA11Mat(ILEV)%a
  AA22Mat     => mg_AA22Mat(ILEV)%a
  AA33Mat     => mg_AA33Mat(ILEV)%a
  AA12Mat     => mg_AA12Mat(ILEV)%a
  AA13Mat     => mg_AA13Mat(ILEV)%a
  AA23Mat     => mg_AA23Mat(ILEV)%a
  AA21Mat     => mg_AA21Mat(ILEV)%a
  AA31Mat     => mg_AA31Mat(ILEV)%a
  AA32Mat     => mg_AA32Mat(ILEV)%a
 ELSE
  A11Mat     => mg_A11Mat(ILEV)%a
  A22Mat     => mg_A22Mat(ILEV)%a
  A33Mat     => mg_A33Mat(ILEV)%a
  AA11Mat     => mg_AA11Mat(ILEV)%a
  AA22Mat     => mg_AA22Mat(ILEV)%a
  AA33Mat     => mg_AA33Mat(ILEV)%a
 END IF

 IF (myMatrixRenewal%S.GE.1) THEN
  S11Mat   => mg_S11Mat(ILEV)%a
  S22Mat   => mg_S22Mat(ILEV)%a
  S33Mat   => mg_S33Mat(ILEV)%a
  S12Mat   => mg_S12Mat(ILEV)%a
  S13Mat   => mg_S13Mat(ILEV)%a
  S23Mat   => mg_S23Mat(ILEV)%a
  S21Mat   => mg_S21Mat(ILEV)%a
  S31Mat   => mg_S31Mat(ILEV)%a
  S32Mat   => mg_S32Mat(ILEV)%a
 END IF

 IF (myMatrixRenewal%D.GE.1) THEN
  DMat     => mg_DMat(ILEV)%a
 END IF

 IF (myMatrixRenewal%K.GE.1) THEN
  KMat     => mg_KMat(ILEV)%a
     IF(bNonNewtonian) THEN
      barM11Mat     => mg_barM11Mat(ILEV)%a
      barM22Mat     => mg_barM22Mat(ILEV)%a
      barM33Mat     => mg_barM33Mat(ILEV)%a
      barM12Mat     => mg_barM12Mat(ILEV)%a
      barM13Mat     => mg_barM13Mat(ILEV)%a
      barM23Mat     => mg_barM23Mat(ILEV)%a
      barM21Mat     => mg_barM21Mat(ILEV)%a
      barM31Mat     => mg_barM31Mat(ILEV)%a
      barM32Mat     => mg_barM32Mat(ILEV)%a
     ELSE
      barM11Mat     => mg_barM11Mat(ILEV)%a
      barM22Mat     => mg_barM22Mat(ILEV)%a
      barM33Mat     => mg_barM33Mat(ILEV)%a
     END IF
 END IF

 IF (myMatrixRenewal%M.GE.1) THEN
  MMat     => mg_MMat(ILEV)%a
  MlRhoMat => mg_MlRhoMat(ILEV)%a
 END IF

 qMat     => mg_qMat(ILEV)
    !!-------------------    POINTER Setup  -------------------!!

 ! Build up the rhs
 IF (idef.eq. 1) THEN

   myScalar%defU = 0d0
   myScalar%defV = 0d0
   myScalar%defW = 0d0


!! BE, CN or starting !!
   IF (tsm.EQ.0 .OR. itns.EQ.1) THEN

   IF (myMatrixRenewal%M.GE.1) THEN
! DAX := A1*DA*DX + A2*DAX 
! SUBROUTINE LAX17(DA,KCOL,KLD,NEQ,DX,DAX,A1,A2)
     !M u_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defU,1d0,1d0)
     !M v_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defV,1d0,1d0)
     !M w_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defW,1d0,1d0)
   END IF

   IF (myMatrixRenewal%S.GE.1) THEN
     !zeitstep S11 u_n
     CALL LAX17(S11Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defU,-zeitstep,1d0)
     !zeitstep S12 v_n
     CALL LAX17(S12Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defU,-zeitstep,1d0)
     !zeitstep S13 w_n
     CALL LAX17(S13Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defU,-zeitstep,1d0)

     !zeitstep S21 u_n
     CALL LAX17(S21Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defV,-zeitstep,1d0)
     !zeitstep S22 v_n
     CALL LAX17(S22Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defV,-zeitstep,1d0)
     !zeitstep S23 w_n
     CALL LAX17(S23Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defV,-zeitstep,1d0)

     !zeitstep S31 u_n
     CALL LAX17(S31Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defW,-zeitstep,1d0)
     !zeitstep S32 v_n
     CALL LAX17(S32Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defW,-zeitstep,1d0)
     !zeitstep S33 w_n
     CALL LAX17(S33Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defW,-zeitstep,1d0)
   END IF

   IF (myMatrixRenewal%K.GE.1) THEN
     !zeitstep K u_n
     CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defU,-zeitstep,1d0)
     !zeitstep K v_n
     CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defV,-zeitstep,1d0)
     !zeitstep K w_n
     CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defW,-zeitstep,1d0)
   END IF


!! BDF(2) !!
   ELSE IF (tsm.EQ.2 .OR. itns.EQ.2) THEN

   IF (myMatrixRenewal%M.GE.1) THEN
! DAX := A1*DA*DX + A2*DAX 
! SUBROUTINE LAX17(DA,KCOL,KLD,NEQ,DX,DAX,A1,A2)
     !4/3 M u_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defU,4d0/3d0,1d0)
     !4/3 M v_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defV,4d0/3d0,1d0)
     !4/3 M w_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defW,4d0/3d0,1d0)

     !-1/3 M u_(n-1)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU_old1,myScalar%defU,-1d0/3d0,1d0)
     !-1/3 M v_(n-1)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV_old1,myScalar%defV,-1d0/3d0,1d0)
     !-1/3 M w_(n-1)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW_old1,myScalar%defW,-1d0/3d0,1d0)
   END IF

!! BDF(3) !!
   ELSE IF (tsm.EQ.3) THEN

   IF (myMatrixRenewal%M.GE.1) THEN
! DAX := A1*DA*DX + A2*DAX 
! SUBROUTINE LAX17(DA,KCOL,KLD,NEQ,DX,DAX,A1,A2)
     !18/11 M u_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defU,18d0/11d0,1d0)
     !18/11 M v_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defV,18d0/11d0,1d0)
     !18/11 M w_n
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defW,18d0/11d0,1d0)

     !-9/11 M u_(n-1)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU_old1,myScalar%defU,-9d0/11d0,1d0)
     !-9/11 M v_(n-1)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV_old1,myScalar%defV,-9d0/11d0,1d0)
     !-9/11 M w_(n-1)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW_old1,myScalar%defW,-9d0/11d0,1d0)

     !2/11 M u_(n-2)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU_old2,myScalar%defU,2d0/11d0,1d0)
     !2/11 M v_(n-2)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV_old2,myScalar%defV,2d0/11d0,1d0)
     !2/11 M w_(n-2)
     CALL LAX17(Mmat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW_old2,myScalar%defW,2d0/11d0,1d0)
   END IF
   END IF
   

 END IF

END SUBROUTINE Matdef_general_QuadScalar_cc
!
! ----------------------------------------------
!
subroutine output_mata()
implicit none

CALL OutputMatrix("MatG",qMat,A11mat,3)

end subroutine output_mata

subroutine output_mmat()
implicit none

CALL OutputMatrix("MatM",qMat,MMat,3)

end subroutine output_mmat
!
! ----------------------------------------------
!
SUBROUTINE Create_P1MMat()
EXTERNAL E012

 IF (.not.ALLOCATED(mg_P1MMat))      ALLOCATE(mg_P1MMat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_P1iMMat))      ALLOCATE(mg_P1iMMat(NLMIN:NLMAX))

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') "P1 Mass matrix and its inverse "

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  IF (.not.ALLOCATED(mg_P1MMat(ILEV)%a)) ALLOCATE(mg_P1MMat(ILEV)%a(4*4*KNEL(ILEV)))
  IF (.not.ALLOCATED(mg_P1iMMat(ILEV)%a)) ALLOCATE(mg_P1iMMat(ILEV)%a(4*4*KNEL(ILEV)))

  mg_P1MMat(ILEV)%a = 0d0

  CALL Build_P1Mass(mg_P1MMat(ILEV)%a,&
                    mg_mesh%level(ilev)%kvert,&
                    mg_mesh%level(ilev)%karea,&
                    mg_mesh%level(ilev)%kedge,&
                    mg_mesh%level(ilev)%dcorvg,&
                    7,E012)

  CALL InversionOfP1Mass(mg_P1MMat(ILEV)%a,mg_P1iMMat(ILEV)%a,KNEL(ILEV))

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') "created!"

 CONTAINS

 SUBROUTINE InversionOfP1Mass(mMat,imMat,N)
 INTEGER I,N
 REAL*8 mMat(4,4,*),imMat(4,4,*)
 LOGICAL BFLAG

 DO I=1,N
  CALL FINDInv(mMat(:,:,I), imMat(:,:,I), 4, bflag)
 END DO

 END SUBROUTINE InversionOfP1Mass

END SUBROUTINE Create_P1MMat
!
! ----------------------------------------------
!
SUBROUTINE Create_Special_CCStructures()
INTEGER I,J,istat
CHARACTER cFile*20

IF (myid.ne.0) THEN
 
 ALLOCATE(my_mg_CCPiece(1:3))
 DO ILEV=NLMIN,NLMAX
   CALL SETLEV(2)
   ALLOCATE(my_mg_CCPiece(ILEV)%a(3,(2**(ILEV)+1)**3))
   WRITE(Cfile,'(A,I2.2,A)')'_data/struct_',ilev,'.dat'

   open (unit=541,file=ADJUSTL(TRIM(cFile)),action="read",iostat=istat)
   if(istat .ne. 0)then
     write(*,*)"Could not open file for writing in Create_Special_CCStructures."
   stop          
   end if    

   DO i=1,(2**(ILEV)+1)**3
    READ(541,*) my_mg_CCPiece(ILEV)%a(:,i)
   END DO
   READ(541,*) 

   my_mg_CCPiece(ILEV)%nu_qq = (2**(ILEV)+1)**3
   ALLOCATE(my_mg_CCPiece(ILEV)%LdA_qq(my_mg_CCPiece(ILEV)%nu_qq+1))
   DO i=1,my_mg_CCPiece(ILEV)%nu_qq+1
    READ(541,*) my_mg_CCPiece(ILEV)%LdA_qq(i)
   END DO   
   READ(541,*) 
   my_mg_CCPiece(ILEV)%na_qq = my_mg_CCPiece(ILEV)%LdA_qq(my_mg_CCPiece(ILEV)%nu_qq+1)-1
   ALLOCATE(my_mg_CCPiece(ILEV)%ColA_qq(my_mg_CCPiece(ILEV)%na_qq))
   DO i=1,my_mg_CCPiece(ILEV)%nu_qq
    READ(541,*) (my_mg_CCPiece(ILEV)%ColA_qq(j),j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1)
   END DO   

   my_mg_CCPiece(ILEV)%nu_lq = 4*((2**(ILEV-1))**3)
   READ(541,*) 
   ALLOCATE(my_mg_CCPiece(ILEV)%LdA_lq(my_mg_CCPiece(ILEV)%nu_lq+1))
   DO i=1,my_mg_CCPiece(ILEV)%nu_lq+1
    READ(541,*) my_mg_CCPiece(ILEV)%LdA_lq(i)
   END DO   
   READ(541,*) 
   my_mg_CCPiece(ILEV)%na_lq = my_mg_CCPiece(ILEV)%LdA_lq(my_mg_CCPiece(ILEV)%nu_lq+1)-1
!    WRITE(*,*) my_mg_CCPiece(ILEV)%nu_lq,my_mg_CCPiece(ILEV)%na_lq
   ALLOCATE(my_mg_CCPiece(ILEV)%ColA_lq(my_mg_CCPiece(ILEV)%na_lq))
   DO i=1,my_mg_CCPiece(ILEV)%nu_lq
    READ(541,*) (my_mg_CCPiece(ILEV)%ColA_lq(j),j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1)
   END DO   

   my_mg_CCPiece(ILEV)%nu_ql = (2**(ILEV)+1)**3
   READ(541,*) 
   ALLOCATE(my_mg_CCPiece(ILEV)%LdA_ql(my_mg_CCPiece(ILEV)%nu_ql+1))
   DO i=1,my_mg_CCPiece(ILEV)%nu_ql+1
    READ(541,*) my_mg_CCPiece(ILEV)%LdA_ql(i)
   END DO   
   READ(541,*) 
   my_mg_CCPiece(ILEV)%na_ql = my_mg_CCPiece(ILEV)%LdA_ql(my_mg_CCPiece(ILEV)%nu_ql+1)-1
!    WRITE(*,*) my_mg_CCPiece(ILEV)%nu_ql,my_mg_CCPiece(ILEV)%na_ql
   ALLOCATE(my_mg_CCPiece(ILEV)%ColA_ql(my_mg_CCPiece(ILEV)%na_ql))
   DO i=1,my_mg_CCPiece(ILEV)%nu_ql
    READ(541,*) (my_mg_CCPiece(ILEV)%ColA_ql(j),j=my_mg_CCPiece(ILEV)%LdA_ql(i),my_mg_CCPiece(ILEV)%LdA_ql(i+1)-1)
   END DO   

   CLOSE(541)
 END DO 

!  ILEV = 3
!  CALL SETLEV(2)
!  
!  call WRITE_TRI_MESH(KWORK(L(LVERT)),DWORK(L(KLCVG(ILEV+1))),nel,nvt)

!  WRITE(*,*) '   '

 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)
  CALL CC_CreateStruct()
  ALLOCATE(my_mg_CCPiece(ILEV)%E(KNEL(1)))
  DO i=1,KNEL(1)
   ALLOCATE(my_mg_CCPiece(ILEV)%E(i)%E_qq(my_mg_CCPiece(ILEV)%na_qq))
   ALLOCATE(my_mg_CCPiece(ILEV)%E(i)%E_ql(my_mg_CCPiece(ILEV)%na_ql))
   ALLOCATE(my_mg_CCPiece(ILEV)%E(i)%E_lq(my_mg_CCPiece(ILEV)%na_lq))
  END DO
 END DO

 ILEV = 1
 CALL SETLEV(2)
 DO i=1,KNEL(1)

  call GetIndices1(i,&
                   mg_mesh%level(ilev)%kvert,&
                   mg_mesh%level(ilev)%kedge,&
                   mg_mesh%level(ilev)%karea,&
                   mg_mesh%level(ilev)%nel,&
                   mg_mesh%level(ilev)%nvt,&
                   mg_mesh%level(ilev)%net,&
                   mg_mesh%level(ilev)%nat)

 END DO

 ILEV = 2
 CALL SETLEV(2)
 DO i=1,KNEL(1)
  call GetIndices2(i,&
                   mg_mesh%level(ilev)%kadj,&
                   mg_mesh%level(ilev)%kvert,&
                   mg_mesh%level(ilev)%kedge,&
                   mg_mesh%level(ilev)%karea,&
                   mg_mesh%level(ilev)%nel,&
                   mg_mesh%level(ilev)%nvt,&
                   mg_mesh%level(ilev)%net,&
                   mg_mesh%level(ilev)%nat)
 END DO
 
 IF (NLMAX.GE.3) THEN
  ILEV = 3
  CALL SETLEV(2)
  DO i=1,KNEL(1)
   call GetIndices3(i,&
                    mg_mesh%level(ilev-1)%kadj,&
                    mg_mesh%level(ilev)%kadj,&
                    mg_mesh%level(ilev)%kvert,&
                    mg_mesh%level(ilev)%kedge,&
                    mg_mesh%level(ilev)%karea,&
                    mg_mesh%level(ilev)%nel,&
                    mg_mesh%level(ilev)%nvt,&
                    mg_mesh%level(ilev)%net,&
                    mg_mesh%level(ilev)%nat)
  END DO
 END IF

END IF

END SUBROUTINE Create_Special_CCStructures
!
! ----------------------------------------------
!
SUBROUTINE GetIndices3(myel,kadj2,kadj3,kvert,kedge,karea,nel,nvt,net,nat)
INTEGER kadj2(6,*),kadj3(6,*),kvert(8,*),kedge(12,*),karea(6,*),nel,nvt,net,nat
INTEGER myel,i,ivt,iel,JEL(64),ind,nEntry
INTEGER iComp,iaux, j,jj,kk,iLoc,jLoc,iGlob,jGlob,myEntry,nLoc,iPos,iMod,jAux
INTEGER, ALLOCATABLE :: pairV(:,:),pairE(:,:)

JEL(1)  = myel
JEL(2)  = KADJ2(3,JEL(1))
JEL(3)  = KADJ2(3,JEL(2))
JEL(4)  = KADJ2(3,JEL(3))
JEL(5)  = KADJ2(6,JEL(1))
JEL(6)  = KADJ2(3,JEL(5))
JEL(7)  = KADJ2(3,JEL(6))
JEL(8)  = KADJ2(3,JEL(7))

j = 8
DO i=1,8
 JEL(j+1)  = KADJ3(3,JEL(i)); j = j+1
 JEL(j+1)  = KADJ3(3,JEL(j)); j = j+1
 JEL(j+1)  = KADJ3(3,JEL(j)); j = j+1
 JEL(j+1)  = KADJ3(6,JEL(i)); j = j+1
 JEL(j+1)  = KADJ3(3,JEL(j)); j = j+1
 JEL(j+1)  = KADJ3(3,JEL(j)); j = j+1
 JEL(j+1)  = KADJ3(3,JEL(j)); j = j+1
END DO

nEntry = SIZE(my_mg_CCPiece(ILEV)%a(1,:))
ALLOCATE(pairV(4,nEntry))
ALLOCATE(pairE(4,64))

DO i=1,64
 pairE(:,i) = [i,JEL(i),i,JEL(i)]
END DO

DO i=1,nEntry
 
 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.1) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = kvert(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.2) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+kedge(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.3) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+net+karea(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.4) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+net+nat+jel(iel)
 END IF

 pairV(:,i) = [i,ind,i,ind]

END DO

CALL SORT2D(pairV(4,1:nEntry),pairV(3,1:nEntry),nEntry)

DO i=1,my_mg_CCPiece(ILEV)%nu_qq

 iLoc  = i; iGlob = pairV(2,iLoc)
 iaux = my_mg_CCPiece(ILEV)%nu_qq
        
 DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1

   jLoc  = my_mg_CCPiece(ILEV)%ColA_qq(j)
   jGlob = pairV(2,jLoc)
   nLoc = mg_qMat(ILEV)%LdA(iGlob+1)-mg_qMat(ILEV)%LdA(iGlob)
   iPos = mg_qMat(ILEV)%LdA(iGlob)
   CALL FindEntry(mg_qMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
   my_mg_CCPiece(ILEV)%E(myel)%E_qq(j) = iPos + myEntry
   
 END DO

 DO j=my_mg_CCPiece(ILEV)%LdA_ql(i),my_mg_CCPiece(ILEV)%LdA_ql(i+1)-1

   iMod = MOD(my_mg_CCPiece(ILEV)%ColA_ql(j)-1,4) + 1
   jLoc = (my_mg_CCPiece(ILEV)%ColA_ql(j)-iMod)/4 + 1
   jAux = pairE(2,jLoc)
   jGlob = 4*(jAux-1)+iMod
   nLoc = mg_qlMat(ILEV)%LdA(iGlob+1)-mg_qlMat(ILEV)%LdA(iGlob)
   iPos = mg_qlMat(ILEV)%LdA(iGlob)
   CALL FindEntry(mg_qlMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
   my_mg_CCPiece(ILEV)%E(myel)%E_ql(j) = iPos + myEntry

 END DO

END DO

DO i=1,my_mg_CCPiece(ILEV)%nu_lq

 iMod = MOD(i-1,4) + 1
 iLoc = (i-iMod)/4 + 1
 jAux = pairE(2,iLoc)
 iGlob = 4*(jAux-1)+iMod
   
 DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1

  jLoc  = my_mg_CCPiece(ILEV)%ColA_lq(j)
  jGlob = pairV(2,jLoc)
  nLoc = mg_lqMat(ILEV)%LdA(iGlob+1)-mg_lqMat(ILEV)%LdA(iGlob)
  iPos = mg_lqMat(ILEV)%LdA(iGlob)
  CALL FindEntry(mg_lqMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
  my_mg_CCPiece(ILEV)%E(myel)%E_lq(j) = iPos + myEntry

 END DO
        
END DO

ALLOCATE (my_mg_CCPiece(ILEV)%E(myel)%pairV(4,nEntry))
ALLOCATE (my_mg_CCPiece(ILEV)%E(myel)%pairE(4,64))

my_mg_CCPiece(ILEV)%E(myel)%pairV = pairV
my_mg_CCPiece(ILEV)%E(myel)%pairE = pairE

DEALLOCATE(pairV,pairE)

END SUBROUTINE GetIndices3
!
! ----------------------------------------------
!
SUBROUTINE GetIndices2(myel,kadj,kvert,kedge,karea,nel,nvt,net,nat)
INTEGER kadj(6,*),kvert(8,*),kedge(12,*),karea(6,*),nel,nvt,net,nat
INTEGER myel,i,ivt,iel,JEL(8),ind,nEntry
INTEGER iComp,iaux, j,jj,kk,iLoc,jLoc,iGlob,jGlob,myEntry,nLoc,iPos,iMod,jAux
INTEGER, ALLOCATABLE :: pairV(:,:),pairE(:,:)

JEL(1)  = myel
JEL(2)  = KADJ(3,JEL(1))
JEL(3)  = KADJ(3,JEL(2))
JEL(4)  = KADJ(3,JEL(3))
JEL(5)  = KADJ(6,JEL(1))
JEL(6)  = KADJ(3,JEL(5))
JEL(7)  = KADJ(3,JEL(6))
JEL(8)  = KADJ(3,JEL(7))

nEntry = SIZE(my_mg_CCPiece(ILEV)%a(1,:))
ALLOCATE(pairV(4,nEntry))
ALLOCATE(pairE(4,8))

DO i=1,8
 pairE(:,i) = [i,JEL(i),i,JEL(i)]
END DO

DO i=1,nEntry
 
 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.1) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = kvert(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.2) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+kedge(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.3) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+net+karea(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.4) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+net+nat+jel(iel)
 END IF

 pairV(:,i) = [i,ind,i,ind]

END DO

CALL SORT2D(pairV(4,1:nEntry),pairV(3,1:nEntry),nEntry)

DO i=1,my_mg_CCPiece(ILEV)%nu_qq

 iLoc  = i; iGlob = pairV(2,iLoc)
 iaux = my_mg_CCPiece(ILEV)%nu_qq
        
 DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1

   jLoc  = my_mg_CCPiece(ILEV)%ColA_qq(j)
   jGlob = pairV(2,jLoc)
   nLoc = mg_qMat(ILEV)%LdA(iGlob+1)-mg_qMat(ILEV)%LdA(iGlob)
   iPos = mg_qMat(ILEV)%LdA(iGlob)
   CALL FindEntry(mg_qMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
   my_mg_CCPiece(ILEV)%E(myel)%E_qq(j) = iPos + myEntry
   
 END DO

 DO j=my_mg_CCPiece(ILEV)%LdA_ql(i),my_mg_CCPiece(ILEV)%LdA_ql(i+1)-1

   iMod = MOD(my_mg_CCPiece(ILEV)%ColA_ql(j)-1,4) + 1
   jLoc = (my_mg_CCPiece(ILEV)%ColA_ql(j)-iMod)/4 + 1
   jAux = pairE(2,jLoc)
   jGlob = 4*(jAux-1)+iMod
   nLoc = mg_qlMat(ILEV)%LdA(iGlob+1)-mg_qlMat(ILEV)%LdA(iGlob)
   iPos = mg_qlMat(ILEV)%LdA(iGlob)
   CALL FindEntry(mg_qlMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
   my_mg_CCPiece(ILEV)%E(myel)%E_ql(j) = iPos + myEntry

 END DO

END DO

DO i=1,my_mg_CCPiece(ILEV)%nu_lq

 iMod = MOD(i-1,4) + 1
 iLoc = (i-iMod)/4 + 1
 jAux = pairE(2,iLoc)
 iGlob = 4*(jAux-1)+iMod
   
 DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1

  jLoc  = my_mg_CCPiece(ILEV)%ColA_lq(j)
  jGlob = pairV(2,jLoc)
  nLoc = mg_lqMat(ILEV)%LdA(iGlob+1)-mg_lqMat(ILEV)%LdA(iGlob)
  iPos = mg_lqMat(ILEV)%LdA(iGlob)
  CALL FindEntry(mg_lqMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
  my_mg_CCPiece(ILEV)%E(myel)%E_lq(j) = iPos + myEntry

 END DO
        
END DO

ALLOCATE (my_mg_CCPiece(ILEV)%E(myel)%pairV(4,nEntry))
ALLOCATE (my_mg_CCPiece(ILEV)%E(myel)%pairE(4,8))

my_mg_CCPiece(ILEV)%E(myel)%pairV = pairV
my_mg_CCPiece(ILEV)%E(myel)%pairE = pairE

DEALLOCATE(pairV,pairE)

END SUBROUTINE GetIndices2
!
! ----------------------------------------------
!
SUBROUTINE GetIndices1(myel,kvert,kedge,karea,nel,nvt,net,nat)
INTEGER kvert(8,*),kedge(12,*),karea(6,*),nel,nvt,net,nat
INTEGER myel,i,ivt,iel,JEL(1),ind,nEntry
INTEGER iComp,iaux, j,jj,kk,iLoc,jLoc,iGlob,jGlob,myEntry,nLoc,iPos,iMod,jAux
INTEGER, ALLOCATABLE :: pairV(:,:),pairE(:,:)

JEL(1)  = myel

nEntry = SIZE(my_mg_CCPiece(ILEV)%a(1,:))
ALLOCATE(pairV(4,nEntry))
ALLOCATE(pairE(4,1))

DO i=1,1
 pairE(:,i) = [i,JEL(i),i,JEL(i)]
END DO

DO i=1,nEntry
 
 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.1) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = kvert(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.2) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+kedge(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.3) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+net+karea(ivt,jel(iel))
 END IF

 IF (my_mg_CCPiece(ILEV)%a(1,i).eq.4) THEN
  iel = my_mg_CCPiece(ILEV)%a(2,i)
  ivt = my_mg_CCPiece(ILEV)%a(3,i)
  ind = nvt+net+nat+jel(iel)
 END IF

 pairV(:,i) = [i,ind,i,ind]

END DO

CALL SORT2D(pairV(4,1:nEntry),pairV(3,1:nEntry),nEntry)

! IF (myid.eq.1.and.myel.eq.9)THEN
! WRITE(*,'(27I5)') pairV(2,:)
! WRITE(*,*) pairE(2,:)
! WRITE(*,*) 
! END IF

DO i=1,my_mg_CCPiece(ILEV)%nu_qq

 iLoc  = i; iGlob = pairV(2,iLoc)
 iaux = my_mg_CCPiece(ILEV)%nu_qq
        
 DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1

   jLoc  = my_mg_CCPiece(ILEV)%ColA_qq(j)
   jGlob = pairV(2,jLoc)
   nLoc = mg_qMat(ILEV)%LdA(iGlob+1)-mg_qMat(ILEV)%LdA(iGlob)
   iPos = mg_qMat(ILEV)%LdA(iGlob)
   CALL FindEntry(mg_qMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
!    IF (myid.eq.1.and.myel.eq.9)THEN
!     write(*,*) jGlob, mg_qMat(ILEV)%ColA(iPos + myEntry)
!    END IF
   
   my_mg_CCPiece(ILEV)%E(myel)%E_qq(j) = iPos + myEntry
   
 END DO

!  IF (myid.eq.1.and.myel.eq.9)THEN
!  WRITE(*,'(A1,27I5)') '+',  mg_qMat(ILEV)%ColA(my_mg_CCPiece(ILEV)%E(myel)%E_qq(my_mg_CCPiece(ILEV)%LdA_qq(i):my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1))
!  WRITE(*,'(A1,27I5)') '-',  pairV(2,my_mg_CCPiece(ILEV)%ColA_qq(my_mg_CCPiece(ILEV)%LdA_qq(i):my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1))
!  END IF

 DO j=my_mg_CCPiece(ILEV)%LdA_ql(i),my_mg_CCPiece(ILEV)%LdA_ql(i+1)-1

   iMod = MOD(my_mg_CCPiece(ILEV)%ColA_ql(j)-1,4) + 1
   jLoc = (my_mg_CCPiece(ILEV)%ColA_ql(j)-iMod)/4 + 1
   jAux = pairE(2,jLoc)
   jGlob = 4*(jAux-1)+iMod
   nLoc = mg_qlMat(ILEV)%LdA(iGlob+1)-mg_qlMat(ILEV)%LdA(iGlob)
   iPos = mg_qlMat(ILEV)%LdA(iGlob)
   CALL FindEntry(mg_qlMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
   my_mg_CCPiece(ILEV)%E(myel)%E_ql(j) = iPos + myEntry

 END DO

END DO

! IF (myid.eq.1.and.myel.eq.9)THEN
!  pause
! END IF

 DO i=1,my_mg_CCPiece(ILEV)%nu_lq

 iMod = MOD(i-1,4) + 1
 iLoc = (i-iMod)/4 + 1
 jAux = pairE(2,iLoc)
 iGlob = 4*(jAux-1)+iMod
   
 DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1

  jLoc  = my_mg_CCPiece(ILEV)%ColA_lq(j)
  jGlob = pairV(2,jLoc)
  nLoc = mg_lqMat(ILEV)%LdA(iGlob+1)-mg_lqMat(ILEV)%LdA(iGlob)
  iPos = mg_lqMat(ILEV)%LdA(iGlob)
  CALL FindEntry(mg_lqMat(ILEV)%ColA(iPos:),jGlob,myEntry,nLoc)
  my_mg_CCPiece(ILEV)%E(myel)%E_lq(j) = iPos + myEntry

 END DO
        
END DO

ALLOCATE (my_mg_CCPiece(ILEV)%E(myel)%pairV(4,nEntry))
ALLOCATE (my_mg_CCPiece(ILEV)%E(myel)%pairE(4,1))

my_mg_CCPiece(ILEV)%E(myel)%pairV = pairV
my_mg_CCPiece(ILEV)%E(myel)%pairE = pairE

DEALLOCATE(pairV,pairE)

END SUBROUTINE GetIndices1
!
! ----------------------------------------------
!
SUBROUTINE FindEntry(List,ii,jj,nList)
INTEGEr ii,jj,nList
INTEGER iMid,iMAx,iMin
INTEGER List(nList)

! WRITE(*,*) " new search for :",  ii
! WRITE(*,*) " in :",  List(:)
! WRITE(*,*) 
! WRITE(*,*) 

IF (ii.eq.List(1)) THEN
 jj = 1 - 1 
 RETURN
END IF

iMin = 2
iMid = nList/2
iMax = nList

IF (ii.eq.List(iMax)) THEN
 jj = iMax - 1
 RETURN
END IF

IF (ii.eq.List(iMin)) THEN
 jj = iMin - 1
 RETURN
END IF

DO
 IF (ii.eq.List(iMid)) THEN
  jj = iMid-1
  RETURN
 ELSE
  IF (ii.lt.List(iMid)) THEN
   iMax = iMid
   iMid = (iMax + iMin)/2
!    IF (ii.eq.List(List(2,iMax)) THEN
!     jj = ii
!     EXIT
!    END IF
  ELSE
   iMin = iMid
   iMid = (iMax + iMin)/2
!    IF (ii.eq.List(List(2,iMin)) THEN
!     jj = ii
!     EXIT
!    END IF
  END IF
 END IF
!  WRITE(*,*) iMin,iMid,iMax
END DO

! pause

END SUBROUTINE FindEntry
!
! ----------------------------------------------
!
SUBROUTINE CC_CreateStruct()
INTEGER i,j,iComp,iaux,jj,kk
CHARACTER cFile*14

my_mg_CCPiece(ILEV)%MPatch%nu = 3*my_mg_CCPiece(ILEV)%nu_qq + 1*my_mg_CCPiece(ILEV)%nu_lq
my_mg_CCPiece(ILEV)%MPatch%na = 9*my_mg_CCPiece(ILEV)%na_qq + 3*my_mg_CCPiece(ILEV)%na_lq + 3*my_mg_CCPiece(ILEV)%na_ql

ALLOCATE (my_mg_CCPiece(ILEV)%MPatch%LdA (my_mg_CCPiece(ILEV)%MPatch%nu+1))
ALLOCATE (my_mg_CCPiece(ILEV)%MPatch%ColA(my_mg_CCPiece(ILEV)%MPatch%na  ))
ALLOCATE (my_mg_CCPiece(ILEV)%MPatchCopy%LdA (my_mg_CCPiece(ILEV)%MPatch%nu+1))
ALLOCATE (my_mg_CCPiece(ILEV)%MPatchCopy%ColA(my_mg_CCPiece(ILEV)%MPatch%na  ))

my_mg_CCPiece(ILEV)%MPatch%LdA(1) = 1

jj = 1
DO iComp =1,3
 DO i=1,my_mg_CCPiece(ILEV)%nu_qq
  iaux = 3*(my_mg_CCPiece(ILEV)%LdA_qq(i+1)-my_mg_CCPiece(ILEV)%LdA_qq(i)) + &
         1*(my_mg_CCPiece(ILEV)%LdA_ql(i+1)-my_mg_CCPiece(ILEV)%LdA_ql(i))
  my_mg_CCPiece(ILEV)%MPatch%LdA(jj+1) = my_mg_CCPiece(ILEV)%MPatch%LdA(jj) + iaux
  jj = jj + 1
 END DO
END DO

! WRITE(*,*) '-+-',9*my_mg_CCPiece(ILEV)%na_qq + 3*my_mg_CCPiece(ILEV)%na_lq,my_mg_CCPiece(ILEV)%MPatch%LdA(jj)-1

DO i=1,my_mg_CCPiece(ILEV)%nu_lq
 iaux = 3*(my_mg_CCPiece(ILEV)%LdA_lq(i+1)-my_mg_CCPiece(ILEV)%LdA_lq(i))
 my_mg_CCPiece(ILEV)%MPatch%LdA(jj+1) = my_mg_CCPiece(ILEV)%MPatch%LdA(jj) + iaux
 jj = jj + 1
END DO

! WRITE(*,*) '-+-',9*my_mg_CCPiece(ILEV)%na_qq + 6*my_mg_CCPiece(ILEV)%na_lq,my_mg_CCPiece(ILEV)%MPatch%LdA(jj)-1

jj = 1
kk = 1
DO iComp =1,3
 DO i=1,my_mg_CCPiece(ILEV)%nu_qq

 iaux = my_mg_CCPiece(ILEV)%nu_qq
        
 DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
   my_mg_CCPiece(ILEV)%MPatch%ColA(kk) = 0*iaux + my_mg_CCPiece(ILEV)%ColA_qq(j)
   kk = kk + 1
  END DO
  DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
   my_mg_CCPiece(ILEV)%MPatch%ColA(kk) = 1*iaux + my_mg_CCPiece(ILEV)%ColA_qq(j)
   kk = kk + 1
  END DO
  DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
   my_mg_CCPiece(ILEV)%MPatch%ColA(kk) = 2*iaux + my_mg_CCPiece(ILEV)%ColA_qq(j)
   kk = kk + 1
  END DO
  DO j=my_mg_CCPiece(ILEV)%LdA_ql(i),my_mg_CCPiece(ILEV)%LdA_ql(i+1)-1
   my_mg_CCPiece(ILEV)%MPatch%ColA(kk) = 3*iaux + my_mg_CCPiece(ILEV)%ColA_ql(j)
   kk = kk + 1
  END DO
  jj = jj + 1
 END DO
END DO

! WRITE(*,*) ilev,my_mg_CCPiece(ILEV)%MPatch%nu,my_mg_CCPiece(ILEV)%MPatch%na,jj-1,kk-1

DO i=1,my_mg_CCPiece(ILEV)%nu_lq

 iaux = my_mg_CCPiece(ILEV)%nu_qq
        
 DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1
  my_mg_CCPiece(ILEV)%MPatch%ColA(kk) = 0*iaux + my_mg_CCPiece(ILEV)%ColA_lq(j)
  kk = kk + 1
 END DO
 DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1
  my_mg_CCPiece(ILEV)%MPatch%ColA(kk) = 1*iaux + my_mg_CCPiece(ILEV)%ColA_lq(j)
  kk = kk + 1
 END DO
 DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1
  my_mg_CCPiece(ILEV)%MPatch%ColA(kk) = 2*iaux + my_mg_CCPiece(ILEV)%ColA_lq(j)
  kk = kk + 1
 END DO
 jj = jj + 1
END DO


! WRITE(*,*) ilev,my_mg_CCPiece(ILEV)%MPatch%nu,my_mg_CCPiece(ILEV)%MPatch%na,jj-1,kk-1

! IF (myid.eq.1) THEN
!  WRITE(Cfile,'(A,I2.2,A)')'sparse_',ilev,'.txt'
!  OPEN(FILE=Cfile,UNIT=142)
!  DO i=1,my_mg_CCPiece(ILEV)%MPatch%nu
!   DO j=my_mg_CCPiece(ILEV)%MPatch%LdA(i),my_mg_CCPiece(ILEV)%MPatch%LdA(i+1)-1
!    WRITE(142,*) i,my_mg_CCPiece(ILEV)%MPatch%nu-my_mg_CCPiece(ILEV)%MPatch%ColA(j)+1
!   END DO
!  END DO
!  CLOSE(142)
! END IF

END SUBROUTINE CC_CreateStruct
!
! ----------------------------------------------
!
SUBROUTINE CC_GetDefect(qScalar,lScalar)
TYPE(TQuadScalar) qScalar
TYPE(TLinScalar)  lScalar
INTEGER ndof
EXTERNAL E013

 ILEV = NLMAX
 CALL SETLEV(2)

 IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
  A11Mat     => mg_A11Mat(ILEV)%a
  A22Mat     => mg_A22Mat(ILEV)%a
  A33Mat     => mg_A33Mat(ILEV)%a
  A12Mat     => mg_A12Mat(ILEV)%a
  A13Mat     => mg_A13Mat(ILEV)%a
  A23Mat     => mg_A23Mat(ILEV)%a
  A21Mat     => mg_A21Mat(ILEV)%a
  A31Mat     => mg_A31Mat(ILEV)%a
  A32Mat     => mg_A32Mat(ILEV)%a
 ELSE
  A11Mat     => mg_A11Mat(ILEV)%a
  A22Mat     => mg_A22Mat(ILEV)%a
  A33Mat     => mg_A33Mat(ILEV)%a
 END IF

 IF (myMatrixRenewal%S.GE.1) THEN
  S11Mat   => mg_S11Mat(ILEV)%a
  S22Mat   => mg_S22Mat(ILEV)%a
  S33Mat   => mg_S33Mat(ILEV)%a
  S12Mat   => mg_S12Mat(ILEV)%a
  S13Mat   => mg_S13Mat(ILEV)%a
  S23Mat   => mg_S23Mat(ILEV)%a
  S21Mat   => mg_S21Mat(ILEV)%a
  S31Mat   => mg_S31Mat(ILEV)%a
  S32Mat   => mg_S32Mat(ILEV)%a
 END IF

 IF (myMatrixRenewal%D.GE.1) THEN
  DMat     => mg_DMat(ILEV)%a
 END IF

 IF (myMatrixRenewal%K.GE.1) THEN
  KMat     => mg_KMat(ILEV)%a
 END IF

 IF (myMatrixRenewal%M.GE.1) THEN
  MMat     => mg_MMat(ILEV)%a
  MlRhoMat => mg_MlRhoMat(ILEV)%a
 END IF

 qMat     => mg_qMat(ILEV)
 
 lScalar%rhsP(ILEV)%x = 0d0
! qScalar%rhsU = 0d0
! qScalar%rhsV = 0d0
! qScalar%rhsW = 0d0

   IF (bNonNewtonian) THEN
    IF(myMatrixRenewal%S.GE.1) THEN

     qScalar%defU = 0d0
     qScalar%defV = 0d0
     qScalar%defW = 0d0
! DAX := A1*DA*DX + A2*DAX 
! SUBROUTINE LAX17(DA,KCOL,KLD,NEQ,DX,DAX,A1,A2)
     CALL LAX17(A11mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valU,qScalar%defU,1d0,1d0)

     CALL LAX17(A12Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valV,qScalar%defU,1d0,1d0)

     CALL LAX17(A13Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valW,qScalar%defU,1d0,1d0)

     CALL LAX17(A21Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valU,qScalar%defV,1d0,1d0)

     CALL LAX17(A22mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valV,qScalar%defV,1d0,1d0)

     CALL LAX17(A23Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valW,qScalar%defV,1d0,1d0)

     CALL LAX17(A31Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valU,qScalar%defW,1d0,1d0)

     CALL LAX17(A32Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valV,qScalar%defW,1d0,1d0)

     CALL LAX17(A33mat,qMat%ColA,qMat%LdA,qMat%nu,&
     qScalar%valW,qScalar%defW,1d0,1d0)

    ELSE
!     qScalar%defU = 0d0
!     qScalar%defV = 0d0
!     qScalar%defW = 0d0

!     IF (myMatrixRenewal%M.GE.1) THEN
!      CALL LAX17(MMat,qMat%ColA,qMat%LdA,qMat%nu,&
!      qScalar%valU,qScalar%defU,zeitstep,1d0)
!      CALL LAX17(MMat,qMat%ColA,qMat%LdA,qMat%nu,&
!      qScalar%valV,qScalar%defV,zeitstep,1d0)
!      CALL LAX17(MMat,qMat%ColA,qMat%LdA,qMat%nu,&
!      qScalar%valW,qScalar%defW,zeitstep,1d0)
!     END IF

!     IF (myMatrixRenewal%K.GE.1) THEN
!      CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
!      qScalar%valU,qScalar%defU,zeitstep,1d0)
!      CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
!      qScalar%valV,qScalar%defV,zeitstep,1d0)
!      CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
!      qScalar%valW,qScalar%defW,zeitstep,1d0)
!     END IF

!!      IF (myMatrixRenewal%D.GE.1) THEN
!!       CALL LAX17(DMat,qMat%ColA,qMat%LdA,qMat%nu,&
!!       qScalar%valU,qScalar%defU,2d0*zeitstep,1d0)
!!       CALL LAX17(DMat,qMat%ColA,qMat%LdA,qMat%nu,&
!!       qScalar%valV,qScalar%defV,2d0*zeitstep,1d0)
!!       CALL LAX17(DMat,qMat%ColA,qMat%LdA,qMat%nu,&
!!       qScalar%valW,qScalar%defW,2d0*zeitstep,1d0)
!!      END IF

!!      CALL ZTIME(tttx0)
!!      ILEV = NLMAX
!!      CALL SETLEV(2)
!     write(*,*)'stress called, which is not adjusted'
!     pause
!     CALL STRESS(qScalar%valU,qScalar%valV,qScalar%valW,&
!     qScalar%defU, qScalar%defV, qScalar%defW,&
!     KWORK(L(LVERT)),KWORK(L(LAREA)),&
!     KWORK(L(LEDGE)),DWORK(L(LCORVG)),E013 ) ! S*u
!!      CALL ZTIME(tttx1)
!!      myStat%tSMat = myStat%tSMat + (tttx1-tttx0)
    END IF   
   ELSE

    CALL LAX17(A11mat,qMat%ColA,qMat%LdA,qMat%nu,&
    qScalar%valU,qScalar%defU,1d0,0d0)
    CALL LAX17(A22mat,qMat%ColA,qMat%LdA,qMat%nu,&
    qScalar%valV,qScalar%defV,1d0,0d0)
    CALL LAX17(A33mat,qMat%ColA,qMat%LdA,qMat%nu,&
    qScalar%valW,qScalar%defW,1d0,0d0)

   END IF

 CALL CC_GetDefect_sub(qScalar%defU,qScalar%defV,qScalar%defW,lScalar%defP(NLMAX)%x,&
      qScalar%ValU,qScalar%ValV,qScalar%ValW,&
      lScalar%ValP(NLMAX)%x,qScalar%ndof,NEL)

 qScalar%defU = qScalar%rhsU - qScalar%defU
 qScalar%defV = qScalar%rhsV - qScalar%defV
 qScalar%defW = qScalar%rhsW - qScalar%defW
 lScalar%defP(ILEV)%x = lScalar%rhsP(ILEV)%x - lScalar%defP(ILEV)%x

END SUBROUTINE CC_GetDefect
!
! ----------------------------------------------
!
SUBROUTINE GetDefNorms(qScalar,lScalar,DefNorm)
TYPE(TQuadScalar) qScalar
TYPE(TLinScalar)  lScalar
REAL*8 DefNorm(4)
INTEGER ndof

IF (myid.ne.0) THEN

 ILEV = NLMAX
 CALL SETLEV(2)

 qScalar%auxU = qScalar%defU
 qScalar%auxV = qScalar%defV
 qScalar%auxW = qScalar%defW
 CALL E013SUM(qScalar%auxU)
 CALL E013SUM(qScalar%auxV)
 CALL E013SUM(qScalar%auxW)

 ndof = qScalar%ndof
 CALL LL21(qScalar%auxU,ndof ,DefNorm(1))
 CALL LL21(qScalar%auxV,ndof ,DefNorm(2))
 CALL LL21(qScalar%auxW,ndof ,DefNorm(3))
 CALL LL21(lScalar%defP(ILEV)%x(      1:),4*nel,DefNorm(4))
END IF

CALL COMM_Maximum(DefNorm(1))
CALL COMM_Maximum(DefNorm(2))
CALL COMM_Maximum(DefNorm(3))
CALL COMM_Maximum(DefNorm(4))

IF (myid.eq.showid) WRITE(*,'(A,4ES12.3)') "non-linear Defect",DefNorm

END SUBROUTINE GetDefNorms
!
! ----------------------------------------------
!
SUBROUTINE CC_Extraction(qScalar)
TYPE(TQuadScalar) qScalar
INTEGER iStringPos,iString

IF (.NOT.(ALLOCATED(CC_EMat))) ALLOCATE(CC_EMat(NLMIN:NLMAX))

DO ILEV = NLMIN,NLMAX
 IF (myid.eq.1) write(*,'(A,I1.1,A$)') " Lev",ILEV,": "
 iStringPos = 0
 iString = 1

 CALL SETLEV(2)

 IF (.NOT.(ALLOCATED(CC_EMat(ILEV)%E)))then
   ALLOCATE(CC_EMat(ILEV)%E(mg_mesh%level(ilev)%nel))
 endif

 IF (myid.ne.0)then
   CALL E013UVWMAT(mg_AA11mat(ILEV)%a,mg_AA22mat(ILEV)%a,&
     mg_AA33mat(ILEV)%a,mg_qMat(ILEV)%LdA,mg_qMat(ILEV)%nu)

   CALL E013BCMAT(qScalar%knprU(ILEV)%x,&
     qScalar%knprV(ILEV)%x,qScalar%knprW(ILEV)%x)

   CALL E013_PUTMAT(mg_AA11mat(ILEV)%a,&
     mg_AA22mat(ILEV)%a,mg_AA33mat(ILEV)%a,mg_qMat(ILEV)%LdA,mg_qMat(ILEV)%nu)
 end if

 CALL CC_Extraction_sub(mg_mesh%level(ilev)%kvert,&
                        mg_mesh%level(ilev)%karea,&
                        mg_mesh%level(ilev)%kedge,&
                        mg_mesh%level(ilev)%nel,&
                        qScalar%knprU(ILEV)%x,&
                        qScalar%knprV(ILEV)%x,&
                        qScalar%knprW(ILEV)%x,&
			iStringPos,iString)

 IF (myid.ne.0)then
   CALL E013_GETMAT(mg_AA11mat(ILEV)%a,mg_AA22mat(ILEV)%a,&
                    mg_AA33mat(ILEV)%a,mg_qMat(ILEV)%LdA,&
                    mg_qMat(ILEV)%nu)
 endif
 IF (myid.eq.1) write(*,*)
END DO

END SUBROUTINE CC_Extraction
!
! ----------------------------------------------
!
SUBROUTINE E013BCMAT(kU,kV,kW)
integer kU(*),kV(*),kW(*)
integer i,j,ndof,ii,jj
real*8 daux

ndof = nvt+net+nat+nel 
DO i = 1,ndof
 IF (kU(i).eq.1) THEN
  ii = mg_qMat(ILEV)%LdA(i)
  daux = mg_AA11Mat(ILEV)%a(ii)
  DO j = ii,mg_qMat(ILEV)%LdA(i+1)-1
   mg_AA11Mat(ILEV)%a(j) = 0d0
   mg_AA12Mat(ILEV)%a(j) = 0d0
   mg_AA13Mat(ILEV)%a(j) = 0d0
  END DO
  mg_AA11Mat(ILEV)%a(ii) = daux
  ii = mg_qlMat(ILEV)%LdA(i)
  DO j = ii,mg_qlMat(ILEV)%LdA(i+1)-1
   mg_BXMat(ILEV)%a(j) = 0d0
  END DO
 END IF

 IF (kV(i).eq.1) THEN
  ii = mg_qMat(ILEV)%LdA(i)
  daux = mg_AA22Mat(ILEV)%a(ii)
  DO j = ii,mg_qMat(ILEV)%LdA(i+1)-1
   mg_AA21Mat(ILEV)%a(j) = 0d0
   mg_AA22Mat(ILEV)%a(j) = 0d0
   mg_AA23Mat(ILEV)%a(j) = 0d0
  END DO
  mg_AA22Mat(ILEV)%a(ii) = daux
  ii = mg_qlMat(ILEV)%LdA(i)
  DO j = ii,mg_qlMat(ILEV)%LdA(i+1)-1
   mg_BYMat(ILEV)%a(j) = 0d0
  END DO
 END IF

  IF (kW(i).eq.1) THEN
  ii = mg_qMat(ILEV)%LdA(i)
  daux = mg_AA33Mat(ILEV)%a(ii)
  DO j = ii,mg_qMat(ILEV)%LdA(i+1)-1
   mg_AA31Mat(ILEV)%a(j) = 0d0
   mg_AA32Mat(ILEV)%a(j) = 0d0
   mg_AA33Mat(ILEV)%a(j) = 0d0
  END DO
  mg_AA33Mat(ILEV)%a(ii) = daux
  ii = mg_qlMat(ILEV)%LdA(i)
  DO j = ii,mg_qlMat(ILEV)%LdA(i+1)-1
   mg_BZMat(ILEV)%a(j) = 0d0
  END DO
 END IF

END DO

! DO i = 1,4*nel
!  ii = mg_lqMat(ILEV)%LdA(i)
!  DO j = ii,mg_lqMat(ILEV)%LdA(i+1)-1
!   jj = mg_lqMat(ILEV)%ColA(j)
!   IF (kU(jj).eq.1) mg_BTXMat(ILEV)%a(j) = 0d0
!   IF (kV(jj).eq.1) mg_BTYMat(ILEV)%a(j) = 0d0
!   IF (kW(jj).eq.1) mg_BTZMat(ILEV)%a(j) = 0d0
!  END DO
! END DO

END SUBROUTINE E013BCMAT
!
! ----------------------------------------------
!
SUBROUTINE Special_CC_Coarse(qScalar)
TYPE(TQuadScalar), INTENT(INOUT), TARGET :: qScalar
INTEGER iLenV,iLenP,jEq,iEq
INTEGER i,j,ii,jj,k
CHARACTER cFile*(20)
REAL*8 dBC,dt


	IF (tsm.EQ.0 .OR. itns.EQ.1) THEN
		dt = tstep
	ELSE 
		dt = zeitstep
	END IF




IF (myid.ne.0) THEN

 ILEV = NLMIN
 CALL SETLEV(2)

 qScalar%auxU=DBLE(myid)
 CALL E013Min(qScalar%auxU)
!  if (myid.eq.4) WRITE(*,*) myid,INT(qScalar%auxU(1:knvt(2)))

 iLenV = nel+nvt+nat+net
 iLenP = 4*nel 

 myCrsMat%na = 9*mg_qMat(ILEV)%na + 6*mg_lqMat(ILEV)%nA
 IF (.not.ALLOCATED(myCrsMat%Row)) ALLOCATE(myCrsMat%Row(myCrsMat%na))
 IF (.not.ALLOCATED(myCrsMat%Col)) ALLOCATE(myCrsMat%Col(myCrsMat%na))
 IF (.not.ALLOCATED(myCrsMat%A))   ALLOCATE(myCrsMat%A  (myCrsMat%na))
 
!  WRITE(336,*)  9*mg_qMat(ILEV)%na + 6*mg_lqMat(ILEV)%nA
 
 k = 0
 do i=1,mg_qMat(ILEV)%nu
  ii = 0*iLenV + i
  iEq = i
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 0*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA11Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 1*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA12Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 2*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA13Mat(ILEV)%a(j)
  end do
  do j=mg_qlMat(ILEV)%LdA(i),mg_qlMat(ILEV)%LdA(i+1)-1
   jj = 3*iLenV + mg_qlMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = -mg_BXMat(ILEV)%a(j)*dt
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(jj)
   myCrsMat%Col(k) = GlobalNumbering(ii)
   myCrsMat%A  (k) = -mg_BXMat(ILEV)%a(j)
  end do
 end do

 do i=1,mg_qMat(ILEV)%nu
  iEq = i
  ii = 1*iLenV+i
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 0*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA21Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 1*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA22Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 2*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA23Mat(ILEV)%a(j)
  end do
  do j=mg_qlMat(ILEV)%LdA(i),mg_qlMat(ILEV)%LdA(i+1)-1
   jj = 3*iLenV + mg_qlMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = -mg_BYMat(ILEV)%a(j)*dt
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(jj)
   myCrsMat%Col(k) = GlobalNumbering(ii)
   myCrsMat%A  (k) = -mg_BYMat(ILEV)%a(j)
  end do
 end do

 do i=1,mg_qMat(ILEV)%nu
  iEq = i
  ii = 2*iLenV+i
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 0*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA31Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 1*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA32Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 2*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = mg_AA33Mat(ILEV)%a(j)
  end do
  do j=mg_qlMat(ILEV)%LdA(i),mg_qlMat(ILEV)%LdA(i+1)-1
   jj = 3*iLenV + mg_qlMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(ii)
   myCrsMat%Col(k) = GlobalNumbering(jj)
   myCrsMat%A  (k) = -mg_BZMat(ILEV)%a(j)*dt
   k = k + 1
   myCrsMat%Row(k) = GlobalNumbering(jj)
   myCrsMat%Col(k) = GlobalNumbering(ii)
   myCrsMat%A  (k) = -mg_BZMat(ILEV)%a(j)
  end do
 end do

 !  do i=1,mg_lqMat(ILEV)%nu
!   ii = 3*iLenV+i
!   do j=mg_lqMat(ILEV)%LdA(i),mg_lqMat(ILEV)%LdA(i+1)-1
!    jj = 0*iLenV + mg_lqMat(ILEV)%ColA(j)
!    k = k + 1
!    myCrsMat%Row(k) = GlobalNumbering(ii)
!    myCrsMat%Col(k) = GlobalNumbering(jj)
!    myCrsMat%A  (k) = mg_BTXMat(ILEV)%a(j)
!   end do
!   do j=mg_lqMat(ILEV)%LdA(i),mg_lqMat(ILEV)%LdA(i+1)-1
!    jj = 1*iLenV + mg_lqMat(ILEV)%ColA(j)
!    k = k + 1
!    myCrsMat%Row(k) = GlobalNumbering(ii)
!    myCrsMat%Col(k) = GlobalNumbering(jj)
!    myCrsMat%A  (k) = mg_BTYMat(ILEV)%a(j)
!   end do
!   do j=mg_lqMat(ILEV)%LdA(i),mg_lqMat(ILEV)%LdA(i+1)-1
!    jj = 2*iLenV + mg_lqMat(ILEV)%ColA(j)
!    k = k + 1
!    myCrsMat%Row(k) = GlobalNumbering(ii)
!    myCrsMat%Col(k) = GlobalNumbering(jj)
!    myCrsMat%A  (k) = mg_BTZMat(ILEV)%a(j)
!   end do
!  end do

else

 ILEV = NLMIN
 CALL SETLEV(2)

 myCrsMat%na = 9*mg_qMat(ILEV)%na + 6*mg_lqMat(ILEV)%nA
 myCrsMat%nu = 3*mg_qMat(ILEV)%nu + 1*mg_lqMat(ILEV)%nu
 
 IF (.not.ALLOCATED(myCrsMat%Row)) ALLOCATE(myCrsMat%Row(myCrsMat%na))
 IF (.not.ALLOCATED(myCrsMat%Col)) ALLOCATE(myCrsMat%Col(myCrsMat%na))
 IF (.not.ALLOCATED(myCrsMat%A))   ALLOCATE(myCrsMat%A  (myCrsMat%na))
 IF (.not.ALLOCATED(myCrsMat%D))   ALLOCATE(myCrsMat%D(  myCrsMat%nu))

 iLenV = nel+nvt+nat+net
 iLenP = 4*nel 

!  iLenV = 1*((2**ilev) + 1)**3
 iLenP = 4*(2**(ilev-1)+0)**3
 k = 0
 
 do i=1,mg_qMat(ILEV)%nu
  ii = 0*iLenV + i
  dBc = 1d0
  IF (qScalar%knprU(ILEV)%x(i).ne.0) dBC = 0d0
  
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jEq = mg_qMat(ILEV)%ColA(j)
   jj = 0*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA11Mat(ILEV)%a(j)
   IF (dBC.eq.0d0.and.ii.eq.jj) myCrsMat%A  (k) = 1d0
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jj = 1*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA12Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jj = 2*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA13Mat(ILEV)%a(j)
  end do
  do j=mg_qlMat(ILEV)%LdA(i),mg_qlMat(ILEV)%LdA(i+1)-1
   jj = 3*iLenV + mg_qlMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = -dBC*mg_BXMat(ILEV)%a(j)*dt
   k = k + 1
   myCrsMat%Row(k) = jj
   myCrsMat%Col(k) = ii
   myCrsMat%A  (k) = -mg_BXMat(ILEV)%a(j)
  end do
 end do

 do i=1,mg_qMat(ILEV)%nu
  dBc = 1d0
  IF (qScalar%knprV(ILEV)%x(i).ne.0) dBC = 0d0
  ii = 1*iLenV+i
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jj = 0*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA21Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jj = 1*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA22Mat(ILEV)%a(j)
   IF (dBC.eq.0d0.and.ii.eq.jj) myCrsMat%A  (k) = 1d0
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jj = 2*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA23Mat(ILEV)%a(j)
  end do
  do j=mg_qlMat(ILEV)%LdA(i),mg_qlMat(ILEV)%LdA(i+1)-1
   jj = 3*iLenV + mg_qlMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = -dBC*mg_BYMat(ILEV)%a(j)*dt
   k = k + 1
   myCrsMat%Row(k) = jj
   myCrsMat%Col(k) = ii
   myCrsMat%A  (k) = -mg_BYMat(ILEV)%a(j)
  end do
 end do
 do i=1,mg_qMat(ILEV)%nu
  dBc = 1d0
  IF (qScalar%knprW(ILEV)%x(i).ne.0) dBC = 0d0
  ii = 2*iLenV+i
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jj = 0*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA31Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jj = 1*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA32Mat(ILEV)%a(j)
  end do
  do j=mg_qMat(ILEV)%LdA(i),mg_qMat(ILEV)%LdA(i+1)-1
   jj = 2*iLenV + mg_qMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = dBC*mg_AA33Mat(ILEV)%a(j)
   IF (dBC.eq.0d0.and.ii.eq.jj) myCrsMat%A  (k) = 1d0
  end do
  do j=mg_qlMat(ILEV)%LdA(i),mg_qlMat(ILEV)%LdA(i+1)-1
   jj = 3*iLenV + mg_qlMat(ILEV)%ColA(j)
   k = k + 1
   myCrsMat%Row(k) = ii
   myCrsMat%Col(k) = jj
   myCrsMat%A  (k) = -dBC*mg_BZMat(ILEV)%a(j)*dt
   k = k + 1
   myCrsMat%Row(k) = jj
   myCrsMat%Col(k) = ii
   myCrsMat%A  (k) = -mg_BZMat(ILEV)%a(j)
  end do
 end do

!  WRITE(*,*)  myCrsMat%nu,myCrsMat%na

end if

!  WRITE(Cfile,'(A,I2.2,A)')'input_',myid,'.txt'
!  OPEN(FILE=cFile,unit=336)
!  DO i =1,k
!   write(336,*) myCrsMat%Row(i),myCrsMat%Col(i)
!  END DO
!  CLOSE(336)

!  CALL COMM_WriteMatrix()

!WRITE(*,*) 'i like it@!!  ', myid,k
!  pause
   !  CLOSE(336)

END SUBROUTINE Special_CC_Coarse
!
! ----------------------------------------------
!
SUBROUTINE CC_MemFree()
INTEGER IEL

IF (myid.eq.0) THEN
!   WRITE(*,*) "handles before: ", CC_H(1),CC_H(2)
!   CALL myUmfPack_CCcrsFree(CC_H(2))
!   WRITE(*,*) "handles after: ", CC_H(1),CC_H(2)
END IF


 DO ILEV = NLMIN,NLMAX

  CALL SETLEV(2)

  DO IEL=1,NEL
   CALL myUmfPack_CCFree(CC_EMat(ILEV)%E(IEL)%H(1),CC_EMat(ILEV)%E(IEL)%H(2))
  END DO

 END DO
! END IF

ILEV = NLMAX
CALL SETLEV(2)

END SUBROUTINE CC_MemFree
!
! ----------------------------------------------
!
SUBROUTINE CC_mgSolve(qScalar,lScalar,mfile,iter,digitcriterion)
TYPE(TQuadScalar), INTENT(INOUT), TARGET :: qScalar
TYPE(TLinScalar), INTENT(INOUT), TARGET  :: lScalar
REAL*8 daux,iter,digitcriterion
INTEGER mfile,i,j
INTEGER ndof_u,ndof_p
REAL*8 dmaxx(3)

 IF (myid.EQ.0) GOTO 1
 
 ILEV = NLMAX
 CALL SETLEV(2)
 ndof_u = nel+net+nat+nvt
 ndof_p = 4*nel

 qScalar%sol(NLMAX)%x(0*ndof_u+1:1*ndof_u) = qScalar%ValU
 qScalar%sol(NLMAX)%x(1*ndof_u+1:2*ndof_u) = qScalar%ValV
 qScalar%sol(NLMAX)%x(2*ndof_u+1:3*ndof_u) = qScalar%ValW
 
 qScalar%rhs(NLMAX)%x(0*ndof_u+1:1*ndof_u) = qScalar%defU
 qScalar%rhs(NLMAX)%x(1*ndof_u+1:2*ndof_u) = qScalar%defV
 qScalar%rhs(NLMAX)%x(2*ndof_u+1:3*ndof_u) = qScalar%defW
 lScalar%rhsP(NLMAX)%x                     = lScalar%defP(NLMAX)%x !0d0             !Div U = 0d0

! Matrices
 IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
  MyMG%A11  => mg_AA11Mat
  MyMG%A22  => mg_AA22Mat
  MyMG%A33  => mg_AA33Mat
  MyMG%A12  => mg_AA12Mat
  MyMG%A13  => mg_AA13Mat
  MyMG%A23  => mg_AA23Mat
  MyMG%A21  => mg_AA21Mat
  MyMG%A31  => mg_AA31Mat
  MyMG%A32  => mg_AA32Mat
 ELSE
  MyMG%A11    => mg_AA11Mat
  MyMG%A22    => mg_AA22Mat
  MyMG%A33    => mg_AA33Mat
 END IF
 MyMG%BX    => mg_BXMat
 MyMG%BY    => mg_BYMat
 MyMG%BZ    => mg_BZMat
 MyMG%BTX   => mg_BTXMat
 MyMG%BTY   => mg_BTYMat
 MyMG%BTZ   => mg_BTZMat

! Matrix structures 
 MyMG%Lq    => mg_qMat
 MyMG%Lql   => mg_qlMat
 MyMG%Llq   => mg_lqMat

 ! Velociy related vectors
 MyMG%X_u   => qScalar%sol
 MyMG%B_u   => qScalar%rhs
 MyMG%dX_u  => qScalar%dsol
 MyMG%D_u   => qScalar%def
 MyMG%A_u   => qScalar%aux

 ! Pressure related vectors
 MyMG%X_p   => lScalar%ValP
 MyMG%B_p   => lScalar%rhsP
 MyMG%dX_p  => lScalar%dValP
 MyMG%D_p   => lScalar%defP
 MyMG%A_p   => lScalar%auxP

 ! Dirichlet filters
 MyMG%KNPRU => qScalar%knprU(NLMAX)%x
 MyMG%KNPRV => qScalar%knprV(NLMAX)%x
 MyMG%KNPRW => qScalar%knprW(NLMAX)%x

 MyMG%cVariable = "Coupled "

1 CONTINUE 

 MyMG%MinIterCycle       = ccParams%MinIterCycle
 MyMG%MaxIterCycle       = ccParams%MaxIterCycle
 MyMG%nSmootherSteps     = ccParams%nSmootherSteps
 MyMG%VANKA              = ccParams%VANKA
 MyMG%RLX                = ccParams%RLX
 MyMG%CycleType          = ccParams%CycleType
 MyMG%Criterion1         = digitcriterion
 MyMG%Criterion2         = ccParams%Criterion
 MyMG%MinLev		 = ccParams%MinLev
 MyMG%MedLev             = ccParams%MedLev

 daux = DBLE(NLMAX)
 CALL COMM_Maximum(daux)
 MyMG%MaxLev             = NINT(daux)


 CALL MG_Solver_cc(mfile,mterm,iter)

 IF (myid.NE.0) THEN
  qScalar%ValU           = qScalar%dsol(NLMAX)%x(0*ndof_u+1:1*ndof_u)
  qScalar%ValV           = qScalar%dsol(NLMAX)%x(1*ndof_u+1:2*ndof_u)
  qScalar%ValW           = qScalar%dsol(NLMAX)%x(2*ndof_u+1:3*ndof_u)
  lScalar%ValP(NLMAX)%x  = lScalar%dValP(NLMAX)%x
 END IF

END SUBROUTINE CC_mgSolve
!
! ----------------------------------------------
!
SUBROUTINE Special_CCExtraction(qScalar)
TYPE(TQuadScalar) qScalar
INTEGER iel,i,ii,j,jj,ijEntry,iPos,jVert,jel
INTEGER my3by3_LocPatch(3,3),my1by3_LocPatch(1,3)
integer ilenp,ilenv,iO1,iO2
REAL*8,aLLOCATABLE ::  A(:,:)
CHARACTER cFile*20
REAL*8 :: dt
INTEGER iStringPos,iString

	IF (tsm.EQ.0 .OR. itns.EQ.1) THEN
		dt = tstep
	ELSE 
		dt = zeitstep
	END IF

IF (myid.eq.0) RETURN

DO ILEV = NLMIN,NLMAX

 IF (myid.eq.1) write(*,'(A,I1.1,A$)') "Lev",ILEV,": "

 iStringPos = 0
 iString = 1
 
 CALL SETLEV(2)

 CALL E013UVWMAT(mg_AA11mat(ILEV)%a,mg_AA22mat(ILEV)%a,&
                 mg_AA33mat(ILEV)%a,mg_qMat(ILEV)%LdA,&
                 mg_qMat(ILEV)%nu)

 CALL E013BCMAT(qScalar%knprU(ILEV)%x,&
                qScalar%knprV(ILEV)%x,&
                qScalar%knprW(ILEV)%x)

 CALL E013_PUTMAT(mg_AA11mat(ILEV)%a,mg_AA22mat(ILEV)%a,&
                  mg_AA33mat(ILEV)%a,mg_qMat(ILEV)%LdA,&
                  mg_qMat(ILEV)%nu)
 
 DO iel = 1,KNEL(1)

  IF (.not.ALLOCATED(my_mg_CCPiece(ILEV)%E(iel)%A)) ALLOCATE (my_mg_CCPiece(ILEV)%E(iel)%A(my_mg_CCPiece(ILEV)%MPatch%na))
  my_mg_CCPiece(ILEV)%E(iel)%A = 0d0
  
  iPos = 1
  DO i=1,my_mg_CCPiece(ILEV)%nu_qq

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A11Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A12Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A13Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO
   
   DO j=my_mg_CCPiece(ILEV)%LdA_ql(i),my_mg_CCPiece(ILEV)%LdA_ql(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_ql(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = -mg_BXMat(ILEV)%A(ijEntry)*dt
    iPos = iPos + 1
   END DO

  END DO
  
  DO i=1,my_mg_CCPiece(ILEV)%nu_qq

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A21Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A22Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A23Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO
   
   DO j=my_mg_CCPiece(ILEV)%LdA_ql(i),my_mg_CCPiece(ILEV)%LdA_ql(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_ql(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = -mg_BYMat(ILEV)%A(ijEntry)*dt
    iPos = iPos + 1
   END DO

  END DO
  
  DO i=1,my_mg_CCPiece(ILEV)%nu_qq

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A31Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A32Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO

   DO j=my_mg_CCPiece(ILEV)%LdA_qq(i),my_mg_CCPiece(ILEV)%LdA_qq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_qq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_A33Mat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO
   
   DO j=my_mg_CCPiece(ILEV)%LdA_ql(i),my_mg_CCPiece(ILEV)%LdA_ql(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_ql(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = -mg_BZMat(ILEV)%A(ijEntry)*dt
    iPos = iPos + 1
   END DO

  END DO
  
  DO i=1,my_mg_CCPiece(ILEV)%nu_lq

   DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_lq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_BTXMat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO

   DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_lq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_BTYMat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO

   DO j=my_mg_CCPiece(ILEV)%LdA_lq(i),my_mg_CCPiece(ILEV)%LdA_lq(i+1)-1
    ijEntry = my_mg_CCPiece(ILEV)%E(iel)%E_lq(j)
    my_mg_CCPiece(ILEV)%E(iel)%A(iPos) = mg_BTZMat(ILEV)%A(ijEntry)
    iPos = iPos + 1
   END DO
  
  END DO

  my_mg_CCPiece(ILEV)%MPatchCopy%nu  = my_mg_CCPiece(ILEV)%MPatch%nu
  my_mg_CCPiece(ILEV)%MPatchCopy%na  = my_mg_CCPiece(ILEV)%MPatch%na
  my_mg_CCPiece(ILEV)%MPatchCopy%LdA  = my_mg_CCPiece(ILEV)%MPatch%LdA
  my_mg_CCPiece(ILEV)%MPatchCopy%ColA = my_mg_CCPiece(ILEV)%MPatch%ColA

  CALL myUmfPack_CCFactorizeLocalMat(my_mg_CCPiece(ILEV)%E(iel)%A,&
                                     my_mg_CCPiece(ILEV)%MPatchCopy,&
                                     my_mg_CCPiece(ILEV)%E(iel)%sym,&
                                     my_mg_CCPiece(ILEV)%E(iel)%num)


  IF (my_mg_CCPiece(ILEV)%E(iel)%sym.EQ.-1.OR.my_mg_CCPiece(ILEV)%E(iel)%num.EQ.-1) then
   IF (myid.eq.1) write(*,'(A$)') '\'
  ELSE
   if (20*iel/knel(1).ge.iString) THEN
    DO i=iStringPos+1,20*iel/knel(1)
     IF (myid.eq.1) write(*,'(A$)') '%'
    END DO
    iString=20*iel/knel(1)+1
    iStringPos=iString-1
   END IF
  END IF

 
 END DO

 IF (myid.eq.1) write(*,*) 

 CALL E013_GETMAT(mg_A11mat(ILEV)%a,&
                  mg_A22mat(ILEV)%a,&
                  mg_A33mat(ILEV)%a,&
                  mg_qMat(ILEV)%LdA,&
                  mg_qMat(ILEV)%nu)

 END DO

END SUBROUTINE Special_CCExtraction
!
! ----------------------------------------------
!
SUBROUTINE Special_CCMemFree()
INTEGER IEL

IF (myid.ne.0) then

 DO ILEV = NLMIN+1,NLMAX
  DO IEL=1,KNEL(1)
   CALL myUmfPack_CCFree(my_mg_CCPiece(ILEV)%E(iel)%sym,my_mg_CCPiece(ILEV)%E(iel)%num)
  END DO
 END DO

END IF

ILEV = NLMAX
CALL SETLEV(2)

END SUBROUTINE Special_CCMemFree
!
! ----------------------------------------------
!
SUBROUTINE QuadScP1toQ2_cc(lSc,qSc)
TYPE(TLinScalar) lSc
TYPE(TQuadScalar) qSc
INTEGER i
EXTERNAL E013

IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)

 qSc%defU = 0d0
 qSc%defV = 0d0

 CALL IntP1toQ2_cc(qSc%defU,&
                   qSc%defV,&
                   lSc%valP(ILEV)%x,&
                   mg_mesh%level(NLMAX)%kvert,&
                   mg_mesh%level(NLMAX)%karea,&
                   mg_mesh%level(NLMAX)%kedge,&
                   mg_mesh%level(NLMAX)%dcorvg,&
                   E013)

 CALL E013Sum(qSc%defU)
 CALL E013Sum(qSc%defV)

 DO i=1,qSc%ndof
  lSc%Q2(i) = qSc%defU(i)/qSc%defV(i)
 END DO
END IF

END SUBROUTINE QuadScP1toQ2_cc
!
! ----------------------------------------------
! barMMat
! ----------------------------------------------
!
SUBROUTINE Create_barMMat_iso(myScalar)
TYPE(TQuadScalar) myScalar
EXTERNAL E013

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_barM11mat)) ALLOCATE(mg_barM11mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_barM22mat)) ALLOCATE(mg_barM22mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_barM33mat)) ALLOCATE(mg_barM33mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_barM12mat)) ALLOCATE(mg_barM12mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_barM13mat)) ALLOCATE(mg_barM13mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_barM23mat)) ALLOCATE(mg_barM23mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_barM21mat)) ALLOCATE(mg_barM21mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_barM31mat)) ALLOCATE(mg_barM31mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_barM32mat)) ALLOCATE(mg_barM32mat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_barM11mat(ILEV)%a)) ALLOCATE(mg_barM11mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_barM22mat(ILEV)%a)) ALLOCATE(mg_barM22mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_barM33mat(ILEV)%a)) ALLOCATE(mg_barM33mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_barM12mat(ILEV)%a)) ALLOCATE(mg_barM12mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_barM13mat(ILEV)%a)) ALLOCATE(mg_barM13mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_barM23mat(ILEV)%a)) ALLOCATE(mg_barM23mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_barM21mat(ILEV)%a)) ALLOCATE(mg_barM21mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_barM31mat(ILEV)%a)) ALLOCATE(mg_barM31mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_barM32mat(ILEV)%a)) ALLOCATE(mg_barM32mat(ILEV)%a(qMat%na))

  mg_barM11mat(ILEV)%a=0d0
  mg_barM22mat(ILEV)%a=0d0
  mg_barM33mat(ILEV)%a=0d0
  mg_barM12mat(ILEV)%a=0d0
  mg_barM13mat(ILEV)%a=0d0
  mg_barM23mat(ILEV)%a=0d0
  mg_barM21mat(ILEV)%a=0d0
  mg_barM31mat(ILEV)%a=0d0
  mg_barM32mat(ILEV)%a=0d0


  CALL Build_barMMat(mgDensity(ILEV)%x,qMat%na,qMat%ColA,qMat%LdA,&
   mg_mesh%level(ILEV)%kvert,mg_mesh%level(ILEV)%karea,&
   mg_mesh%level(ILEV)%kedge,mg_mesh%level(ILEV)%dcorvg,E013,&
   mg_barM11mat(ILEV)%a,mg_barM12mat(ILEV)%a,mg_barM13mat(ILEV)%a,&
   mg_barM21mat(ILEV)%a,mg_barM22mat(ILEV)%a,mg_barM23mat(ILEV)%a,&
   mg_barM31mat(ILEV)%a,mg_barM32mat(ILEV)%a,mg_barM33mat(ILEV)%a,&
   myScalar%valU,myScalar%valV,myScalar%valW)


 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat  => mg_qMat(NLMAX)
 barM11Mat => mg_barM11Mat(NLMAX)%a
 barM22Mat => mg_barM22Mat(NLMAX)%a
 barM33Mat => mg_barM33Mat(NLMAX)%a
 barM12Mat => mg_barM12Mat(NLMAX)%a
 barM13Mat => mg_barM13Mat(NLMAX)%a
 barM23Mat => mg_barM23Mat(NLMAX)%a
 barM21Mat => mg_barM21Mat(NLMAX)%a
 barM31Mat => mg_barM31Mat(NLMAX)%a
 barM32Mat => mg_barM32Mat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tCMat = myStat%tCMat + (myStat%t1-myStat%t0)


END SUBROUTINE Create_barMMat_iso
!
! ----------------------------------------------
!
!
! ----------------------------------------------
! BMat
! ----------------------------------------------
!
SUBROUTINE Create_BMat_mod() !(B)
INTEGER nERow,pNEL
INTEGER I,J
real*8 ddx,ddy,ddz
CHARACTER*10 myFile
EXTERNAL E011,E013

 IF (.NOT.ALLOCATED(mg_BXMat)) ALLOCATE(mg_BXMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BYMat)) ALLOCATE(mg_BYMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BZMat)) ALLOCATE(mg_BZMat(NLMIN:NLMAX))
! Helping Force calculation
 IF (.NOT.ALLOCATED(mg_BXMat_new)) ALLOCATE(mg_BXMat_new(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BYMat_new)) ALLOCATE(mg_BYMat_new(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BZMat_new)) ALLOCATE(mg_BZMat_new(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)

  qlMat => mg_qlMat(ILEV)
  IF (.NOT.ALLOCATED(mg_BXMat(ILEV)%a)) ALLOCATE(mg_BXMat(ILEV)%a(qlMat%na))
  IF (.NOT.ALLOCATED(mg_BYMat(ILEV)%a)) ALLOCATE(mg_BYMat(ILEV)%a(qlMat%na))
  IF (.NOT.ALLOCATED(mg_BZMat(ILEV)%a)) ALLOCATE(mg_BZMat(ILEV)%a(qlMat%na))
  mg_BXMat(ILEV)%a=0d0
  mg_BYMat(ILEV)%a=0d0
  mg_BZMat(ILEV)%a=0d0
! Helping Force calculation
  IF (.NOT.ALLOCATED(mg_BXMat_new(ILEV)%a)) ALLOCATE(mg_BXMat_new(ILEV)%a(qlMat%na))
  IF (.NOT.ALLOCATED(mg_BYMat_new(ILEV)%a)) ALLOCATE(mg_BYMat_new(ILEV)%a(qlMat%na))
  IF (.NOT.ALLOCATED(mg_BZMat_new(ILEV)%a)) ALLOCATE(mg_BZMat_new(ILEV)%a(qlMat%na))
  mg_BXMat_new(ILEV)%a=0d0
  mg_BYMat_new(ILEV)%a=0d0
  mg_BZMat_new(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [B]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
  CALL Build_BMatP1(mg_BXMat(ILEV)%a,mg_BYMat(ILEV)%a,&
       mg_BZMat(ILEV)%a,qlMat%LdA,qlMat%ColA,&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       qlMat%na,E013)
! Helping Force calculation
  mg_BXMat_new(ILEV)%a = mg_BXMat(ILEV)%a
  mg_BYMat_new(ILEV)%a = mg_BYMat(ILEV)%a
  mg_BZMat_new(ILEV)%a = mg_BZMat(ILEV)%a


 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

 IF (.NOT.ALLOCATED(mg_BTXMat)) ALLOCATE(mg_BTXMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BTYMat)) ALLOCATE(mg_BTYMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BTZMat)) ALLOCATE(mg_BTZMat(NLMIN:NLMAX))

DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)

  lqMat => mg_lqMat(ILEV)
  IF (.NOT.ALLOCATED(mg_BTXMat(ILEV)%a)) ALLOCATE(mg_BTXMat(ILEV)%a(lqMat%na))
  IF (.NOT.ALLOCATED(mg_BTYMat(ILEV)%a)) ALLOCATE(mg_BTYMat(ILEV)%a(lqMat%na))
  IF (.NOT.ALLOCATED(mg_BTZMat(ILEV)%a)) ALLOCATE(mg_BTZMat(ILEV)%a(lqMat%na))
  mg_BTXMat(ILEV)%a=0d0
  mg_BTYMat(ILEV)%a=0d0
  mg_BTZMat(ILEV)%a=0d0


  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [B{T}]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
  CALL Build_BTMatP1(mg_BTXMat(ILEV)%a,mg_BTYMat(ILEV)%a,&
       mg_BTZMat(ILEV)%a,lqMat%LdA,lqMat%ColA,&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       lqMat%na,E013)

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

55 CONTINUE

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat  => mg_qlMat(NLMAX)
 BXMat => mg_BXMat(NLMAX)%a
 BYMat => mg_BYMat(NLMAX)%a
 BZMat => mg_BZMat(NLMAX)%a
! Helping Force calculation
 BXMat_new => mg_BXMat_new(NLMAX)%a
 BYMat_new => mg_BYMat_new(NLMAX)%a
 BZMat_new => mg_BZMat_new(NLMAX)%a

 lqMat  => mg_lqMat(NLMAX)
 BTXMat => mg_BTXMat(NLMAX)%a
 BTYMat => mg_BTYMat(NLMAX)%a
 BTZMat => mg_BTZMat(NLMAX)%a

END SUBROUTINE Create_BMat_mod
!
! 
SUBROUTINE InitializeProlRest_cc(Param)

TYPE(tParamCC) Param
LOGICAL, TARGET :: PR

 PR = .FALSE.
 MyMG%MinLev             = Param%MinLev
 MyMG%MedLev             = Param%MedLev
 myMG%MaxLev=NLMAX
 IF(myid.eq.showid) WRITE(*,*) "Initialization of velocity prolongation matrix"
 myMG%bProlRest => PR
 MyMG%cVariable = "Velocity"
 CALL mgProlRestInit()

 PR = .FALSE.
 IF(myid.eq.showid) WRITE(*,*) "Initialization of pressure prolongation matrix"
 myMG%bProlRest => PR
 MyMG%cVariable = "Pressure"
 CALL mgProlRestInit()

END SUBROUTINE InitializeProlRest_cc
!
! ----------------------------------------------
!
SUBROUTINE Create_AMat_new()
INTEGER NA,complete

IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
  ALLOCATE(mg_A11mat(NLMIN:NLMAX))
  ALLOCATE(mg_A22mat(NLMIN:NLMAX))
  ALLOCATE(mg_A33mat(NLMIN:NLMAX))
  ALLOCATE(mg_A12mat(NLMIN:NLMAX))
  ALLOCATE(mg_A13mat(NLMIN:NLMAX))
  ALLOCATE(mg_A23mat(NLMIN:NLMAX))
  ALLOCATE(mg_A21mat(NLMIN:NLMAX))
  ALLOCATE(mg_A31mat(NLMIN:NLMAX))
  ALLOCATE(mg_A32mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA11mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA22mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA33mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA12mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA13mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA23mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA21mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA31mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA32mat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

   IF (myid.eq.showID) THEN
    IF (ILEV.EQ.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') "Allocation of Aii Matrix on Level [", ILEV,"]"
    END IF
    IF (ILEV.EQ.NLMAX) THEN
     WRITE(MTERM,'(A,I1,A)', advance='yes') ", [",ILEV,"]"
    END IF
    IF (ILEV.NE.NLMAX.AND.ILEV.NE.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
    END IF
   END IF

   NA = mg_qMat(ILEV)%na
   ALLOCATE(mg_A11mat(ILEV)%a(NA))
   ALLOCATE(mg_A22mat(ILEV)%a(NA))
   ALLOCATE(mg_A33mat(ILEV)%a(NA))
   ALLOCATE(mg_A12mat(ILEV)%a(NA))
   ALLOCATE(mg_A13mat(ILEV)%a(NA))
   ALLOCATE(mg_A23mat(ILEV)%a(NA))
   ALLOCATE(mg_A21mat(ILEV)%a(NA))
   ALLOCATE(mg_A31mat(ILEV)%a(NA))
   ALLOCATE(mg_A32mat(ILEV)%a(NA))
   ALLOCATE(mg_AA11mat(ILEV)%a(NA))
   ALLOCATE(mg_AA22mat(ILEV)%a(NA))
   ALLOCATE(mg_AA33mat(ILEV)%a(NA))
   ALLOCATE(mg_AA12mat(ILEV)%a(NA))
   ALLOCATE(mg_AA13mat(ILEV)%a(NA))
   ALLOCATE(mg_AA23mat(ILEV)%a(NA))
   ALLOCATE(mg_AA21mat(ILEV)%a(NA))
   ALLOCATE(mg_AA31mat(ILEV)%a(NA))
   ALLOCATE(mg_AA32mat(ILEV)%a(NA))
  END DO
  A11Mat => mg_A11mat(NLMAX)%a
  A22Mat => mg_A22mat(NLMAX)%a
  A33Mat => mg_A33mat(NLMAX)%a
  A12Mat => mg_A12mat(NLMAX)%a
  A13Mat => mg_A13mat(NLMAX)%a
  A23Mat => mg_A23mat(NLMAX)%a
  A21Mat => mg_A21mat(NLMAX)%a
  A31Mat => mg_A31mat(NLMAX)%a
  A32Mat => mg_A32mat(NLMAX)%a
  AA11Mat => mg_AA11mat(NLMAX)%a
  AA22Mat => mg_AA22mat(NLMAX)%a
  AA33Mat => mg_AA33mat(NLMAX)%a
  AA12Mat => mg_AA12mat(NLMAX)%a
  AA13Mat => mg_AA13mat(NLMAX)%a
  AA23Mat => mg_AA23mat(NLMAX)%a
  AA21Mat => mg_AA21mat(NLMAX)%a
  AA31Mat => mg_AA31mat(NLMAX)%a
  AA32Mat => mg_AA32mat(NLMAX)%a
 ELSE
  ALLOCATE(mg_A11mat(NLMIN:NLMAX))
  ALLOCATE(mg_A22mat(NLMIN:NLMAX))
  ALLOCATE(mg_A33mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA11mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA22mat(NLMIN:NLMAX))
  ALLOCATE(mg_AA33mat(NLMIN:NLMAX))
  DO ILEV=NLMIN,NLMAX

   IF (myid.eq.showID) THEN
    IF (ILEV.EQ.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') "Allocation of A Matrix on Level [", ILEV,"]"
    END IF
    IF (ILEV.EQ.NLMAX) THEN
     WRITE(MTERM,'(A,I1,A)', advance='yes') ", [",ILEV,"]"
    END IF
    IF (ILEV.NE.NLMAX.AND.ILEV.NE.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
    END IF
   END IF

   NA = mg_qMat(ILEV)%na
   ALLOCATE(mg_A11mat(ILEV)%a(NA))
   ALLOCATE(mg_A22mat(ILEV)%a(NA))
   ALLOCATE(mg_A33mat(ILEV)%a(NA))
   ALLOCATE(mg_AA11mat(ILEV)%a(NA))
   ALLOCATE(mg_AA22mat(ILEV)%a(NA))
   ALLOCATE(mg_AA33mat(ILEV)%a(NA))
  END DO
  A11Mat => mg_A11mat(NLMAX)%a
  A22Mat => mg_A22mat(NLMAX)%a
  A33Mat => mg_A33mat(NLMAX)%a
  AA11Mat => mg_AA11mat(NLMAX)%a
  AA22Mat => mg_AA22mat(NLMAX)%a
  AA33Mat => mg_AA33mat(NLMAX)%a
 END IF

END SUBROUTINE Create_AMat_new
!
! ----------------------------------------------
!
SUBROUTINE myOutputMatrix(cFile,myMat,Mat)
TYPE(TMatrix) myMat
REAL*8 Mat(*),DD
CHARACTER*4 cFile
CHARACTER*12 myFile
INTEGER I,J,II,K


 IF (myid.EQ.0) THEN

  WRITE(myFile(1:12),'(A3,A4)') cFile,".txt"

!  WRITE(*,*) myid,myFile
 OPEN(987,FILE=myFile)

 DO I=1,myMat%nu
   DO J=myMat%LdA(I),myMat%LdA(I+1)-1
     WRITE(987,'(I10,I10,E16.8)') I,myMat%ColA(J),Mat(J)
   END DO
 END DO
 CLOSE(987)

 END IF

END SUBROUTINE myOutputMatrix
END MODULE def_cc
