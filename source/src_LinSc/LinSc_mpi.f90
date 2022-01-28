! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
!SUBROUTINE E011Sum(FX) !ok
!USE PP3D_MPI
!USE def_feat, ONLY: ILEV,NLMIN,NLMAX
!USE var_QuadScalar, only: knvt
!IMPLICIT NONE
!
!REAL*8  FX(*)
!INTEGER I,pID,pJD,nSIZE,nEIGH,NU
!
!IF (myid.ne.MASTER) THEN
!
! DO pID=1,subnodes
!  IF (myid.NE.pID) THEN
!   DO pJD=1,subnodes
!    IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
!     nSIZE = E011ST(pJD)%Num
!
!     DO I=1,nSIZE
!      E011ST(pJD)%SDVect(I) = FX(E011ST(pJD)%VertLink(1,I))
!     END DO
!
!     CALL SENDD_myMPI(E011ST(pJD)%SDVect,nSIZE,pID)
!
!    END IF
!   END DO
!  ELSE
!   DO pJD=1,subnodes
!     IF (E011ST(pJD)%Num.GT.0) THEN
!      nSIZE = E011ST(pJD)%Num
!
!      CALL RECVD_myMPI(E011ST(pJD)%RDVect,nSIZE,pJD)
!
!     END IF
!   END DO
!  END IF
! END DO
!
! DO pJD=1,subnodes
!
!   NU = KNVT(ILEV)
!   IF (ILEV.lt.NLMIN.OR.ILEV.GT.NLMAX+1) THEN
!    WRITE(*,*) myid,ILEV,'problem'
!    pause
!   END IF
!   IF (E011ST(pJD)%Num.GT.0) THEN
!     nSIZE = E011ST(pJD)%Num
!
!     DO I=1,nSIZE
!      IF (E011ST(pJD)%VertLink(2,I).le.NU) THEN
!       FX(E011ST(pJD)%VertLink(2,I)) = &
!       FX(E011ST(pJD)%VertLink(2,I)) +  E011ST(pJD)%RDVect(I)
!      END IF
!     END DO
!
!   END IF
! END DO
!
!END IF
!
!if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
!
!END SUBROUTINE E011Sum
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
!SUBROUTINE E011MAT(A,KLDA,NU) !ok
!USE PP3D_MPI
!IMPLICIT NONE
!
!REAL*8 A(*)
!INTEGER KLDA(*),NU
!INTEGER I,J
!INTEGER pID,pJD,nSIZE
!
!IF (myid.ne.MASTER) THEN
!
! DO I=1,NU
!  E011_UE(I)=A(KLDA(I))
! ENDDO
!
! DO pID=1,subnodes
!  IF (myid.NE.pID) THEN
!   DO pJD=1,subnodes
!    IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
!     nSIZE = E011ST(pJD)%Num
!
!     DO I=1,nSIZE
!      E011ST(pJD)%SVVect(I) = E011_UE(E011ST(pJD)%VertLink(1,I))
!     END DO
!
!     CALL SENDV_myMPI(E011ST(pJD)%SVVect,nSIZE,pID)
!
!    END IF
!   END DO
!  ELSE
!   DO pJD=1,subnodes
!    IF (E011ST(pJD)%Num.GT.0) THEN
!
!     nSIZE = E011ST(pJD)%Num
!     CALL RECVV_myMPI(E011ST(pJD)%RVVect,nSIZE,pJD)
!
!    END IF
!
!   END DO
!  END IF
! END DO
!
! DO pJD=1,subnodes
!
!   IF (E011ST(pJD)%Num.GT.0) THEN
!
!     nSIZE = E011ST(pJD)%Num
!
!     DO I=1,nSIZE
!       E011_UE(E011ST(pJD)%VertLink(2,I)) = &
!       E011_UE(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RVVect(I)
!     END DO
!
!   END IF
!
! END DO
!
!END IF
!
!if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
!
!END SUBROUTINE E011MAT
! ----------------------------------------------
SUBROUTINE E011_GenLinSc_Q1_UMAT(A,KLDA,NU,iC) !ok
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE PP3D_MPI
IMPLICIT NONE

REAL*8 A(*)
INTEGER KLDA(*),NU
INTEGER I,J,iC
INTEGER pID,pJD,nSIZE

IF (myid.ne.MASTER) THEN

 DO I=1,NU
  MGE011(ILEV)%UE(iC)%x(I)=A(KLDA(I))
 ENDDO

 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
     nSIZE = E011ST(pJD)%Num

!      WRITE(*,*) myid, nSIZE, SIZE(E011ST(pJD)%SDVect)
     DO I=1,nSIZE
        IF (E011ST(pJD)%VertLink(1,I).le.nu) then
         E011ST(pJD)%SDVect(I) = MGE011(ILEV)%UE(iC)%x(E011ST(pJD)%VertLink(1,I))
        END IF
     END DO

     CALL SENDD_myMPI(E011ST(pJD)%SDVect,nSIZE,pID)

    END IF
   END DO
  ELSE
   DO pJD=1,subnodes
    IF (E011ST(pJD)%Num.GT.0) THEN

     nSIZE = E011ST(pJD)%Num
     CALL RECVD_myMPI(E011ST(pJD)%RDVect,nSIZE,pJD)

    END IF

   END DO
  END IF
 END DO

 DO pJD=1,subnodes

   IF (E011ST(pJD)%Num.GT.0) THEN

     nSIZE = E011ST(pJD)%Num

     DO I=1,nSIZE
       IF (E011ST(pJD)%VertLink(2,I).le.NU) THEN
        MGE011(ILEV)%UE(iC)%x(E011ST(pJD)%VertLink(2,I)) = &
        MGE011(ILEV)%UE(iC)%x(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RDVect(I)
       END IF
     END DO

   END IF

 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011_GenLinSc_Q1_UMAT
! ----------------------------------------------
SUBROUTINE E011_UMAT(A,KLDA,NU,iC) !ok
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE PP3D_MPI
IMPLICIT NONE

REAL*8 A(*)
INTEGER KLDA(*),NU
INTEGER I,J,iC
INTEGER pID,pJD,nSIZE

IF (myid.ne.MASTER) THEN

 DO I=1,NU
  IF (iC.eq.1) MGE011(ILEV)%UE11(I)=A(KLDA(I))
  IF (iC.eq.2) MGE011(ILEV)%UE22(I)=A(KLDA(I))
  IF (iC.eq.3) MGE011(ILEV)%UE33(I)=A(KLDA(I))
 ENDDO

 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
     nSIZE = E011ST(pJD)%Num

!      WRITE(*,*) myid, nSIZE, SIZE(E011ST(pJD)%SDVect)
     DO I=1,nSIZE
        IF (E011ST(pJD)%VertLink(1,I).le.nu) then
         IF (iC.eq.1) E011ST(pJD)%SDVect(I) = MGE011(ILEV)%UE11(E011ST(pJD)%VertLink(1,I))
         IF (iC.eq.2) E011ST(pJD)%SDVect(I) = MGE011(ILEV)%UE22(E011ST(pJD)%VertLink(1,I))
         IF (iC.eq.3) E011ST(pJD)%SDVect(I) = MGE011(ILEV)%UE33(E011ST(pJD)%VertLink(1,I))
        END IF
     END DO

     CALL SENDD_myMPI(E011ST(pJD)%SDVect,nSIZE,pID)

    END IF
   END DO
  ELSE
   DO pJD=1,subnodes
    IF (E011ST(pJD)%Num.GT.0) THEN

     nSIZE = E011ST(pJD)%Num
     CALL RECVD_myMPI(E011ST(pJD)%RDVect,nSIZE,pJD)

    END IF

   END DO
  END IF
 END DO

 DO pJD=1,subnodes

   IF (E011ST(pJD)%Num.GT.0) THEN

     nSIZE = E011ST(pJD)%Num

     DO I=1,nSIZE
       IF (E011ST(pJD)%VertLink(2,I).le.NU) THEN
        IF (iC.eq.1)  MGE011(ILEV)%UE11(E011ST(pJD)%VertLink(2,I)) = &
        MGE011(ILEV)%UE11(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RDVect(I)
        IF (iC.eq.2)  MGE011(ILEV)%UE22(E011ST(pJD)%VertLink(2,I)) = &
        MGE011(ILEV)%UE22(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RDVect(I)
        IF (iC.eq.3)  MGE011(ILEV)%UE33(E011ST(pJD)%VertLink(2,I)) = &
        MGE011(ILEV)%UE33(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RDVect(I)
       END IF
     END DO

   END IF

 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011_UMAT
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E011XYZMAT(A11,A22,A33,KLDA,NU) !ok
USE def_feat, ONLY: ILEV
USE PP3D_MPI
IMPLICIT NONE

REAL*8 A11(*),A22(*),A33(*)
INTEGER KLDA(*),NU
INTEGER I,J
INTEGER pID,pJD,nSIZE
INTEGER MEQ,MEQ1,MEQ2,MEQ3

IF (myid.ne.MASTER) THEN

 DO I=1,NU
  MGE011(ILEV)%UE11(I)=A11(KLDA(I))
  MGE011(ILEV)%UE22(I)=A22(KLDA(I))
  MGE011(ILEV)%UE33(I)=A33(KLDA(I))
 ENDDO

 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
     nSIZE = E011ST(pJD)%Num

     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     DO I=1,nSIZE
      E011ST(pJD)%SDVect(MEQ1+I) = MGE011(ILEV)%UE11(E011ST(pJD)%VertLink(1,I))
      E011ST(pJD)%SDVect(MEQ2+I) = MGE011(ILEV)%UE22(E011ST(pJD)%VertLink(1,I))
      E011ST(pJD)%SDVect(MEQ3+I) = MGE011(ILEV)%UE33(E011ST(pJD)%VertLink(1,I))
     END DO

!      DO I=1,nSIZE
!       E011ST(pJD)%SVVect(I) = E011_UE(E011ST(pJD)%VertLink(1,I))
!      END DO

     CALL SENDD_myMPI(E011ST(pJD)%SDVect,3*nSIZE,pID)

    END IF
   END DO
  ELSE
   DO pJD=1,subnodes
    IF (E011ST(pJD)%Num.GT.0) THEN

     nSIZE = E011ST(pJD)%Num
     CALL RECVD_myMPI(E011ST(pJD)%RDVect,3*nSIZE,pJD)

    END IF

   END DO
  END IF
 END DO

 DO pJD=1,subnodes

   IF (E011ST(pJD)%Num.GT.0) THEN

     nSIZE = E011ST(pJD)%Num
     MEQ = nSize
     MEQ1 = 0
     MEQ2 = MEQ
     MEQ3 = 2*MEQ

     DO I=1,nSIZE
      MGE011(ILEV)%UE11(E011ST(pJD)%VertLink(2,I)) = &
      MGE011(ILEV)%UE11(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RDVect(MEQ1+I)
      MGE011(ILEV)%UE22(E011ST(pJD)%VertLink(2,I)) = &
      MGE011(ILEV)%UE22(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RDVect(MEQ2+I)
      MGE011(ILEV)%UE33(E011ST(pJD)%VertLink(2,I)) = &
      MGE011(ILEV)%UE33(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RDVect(MEQ3+I)
     END DO

   END IF

 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011XYZMAT
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
!SUBROUTINE E011_CreateComm(NVT)
!USE PP3D_MPI
!USE def_feat, ONLY: ILEV,NLMIN,NLMAX
!IMPLICIT NONE
!
!INTEGER NVT,pID,jAux,i
!
!IF (myid.EQ.0) RETURN
!
!ALLOCATE(E011ST(subnodes))
!
!ILEV = NLMAX
!
!DO pID=1,subnodes
!  E011ST(pID)%Num = MGE013(ILEV)%ST(pID)%Num
!END DO
!
!ALLOCATE(E011_UE(NVT))
!
!DO pID=1,subnodes
! IF (E011ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
!  jAux = E011ST(pID)%Num
!  ALLOCATE(E011ST(pID)%VertLink(2,E011ST(pID)%Num))
!  ALLOCATE(E011ST(pID)%SVVect  (  E011ST(pID)%Num))
!  ALLOCATE(E011ST(pID)%RVVect  (  E011ST(pID)%Num))
!  ALLOCATE(E011ST(pID)%SDVect  (  E011ST(pID)%Num))
!  ALLOCATE(E011ST(pID)%RDVect  (  E011ST(pID)%Num))
!  ALLOCATE(E011ST(pID)%SBVect  (  E011ST(pID)%Num))
!  ALLOCATE(E011ST(pID)%RBVect  (  E011ST(pID)%Num))
!  DO i=1,jAux
!   E011ST(pID)%VertLink(1,i) =  MGE013(ILEV)%ST(pID)%VertLink(1,i)
!   E011ST(pID)%VertLink(2,i) =  MGE013(ILEV)%ST(pID)%VertLink(2,i)
!  END DO
! END IF
!END DO
!
!END


