      SUBROUTINE SSORSmoother_GenLinSc_Q1(AXX,KCOL,KLD,
     *           DX,DB,DD,KNPR,NEQ,NIT,NFLD,OMEGA)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE011,myid,E011SUM
      USE var_QuadScalar, ONLY : tMGFldMatrix
      
      IMPLICIT NONE
      REAL*8 DX(*),DB(*),DD(*)
      INTEGER KCOL(*),KLD(*),KNPR(*)
      INTEGER NEQ,NIT,NFLD
      REAL*8 OMEGA
      TYPE(tMGFldMatrix) AXX(*)
      INTEGER iMEQ,IEQ,J,iFLD,ITE,ICOL
      REAL*8 AUX
!
      DO iFLD =1,NFLD

      iMEQ = (iFLD-1)*NEQ
!
      DO 1 ITE=1,NIT
!
!     The real SOR part is here
!     ----------------------------------------------
      DO IEQ=1,NEQ
       IF (KNPR(IEQ).EQ.0) THEN
        AUX=DB(iMEQ+IEQ)
        DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
         J=KCOL(ICOL)
         AUX = AUX - AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
        END DO
        DX(iMEQ+IEQ) = DX(iMEQ+IEQ) + 
     *  OMEGA*(AUX/AXX(iFLD)%FLD(ILEV)%a(KLD(IEQ))-DX(iMEQ+IEQ))
       END IF
      END DO
!     ----------------------------------------------
!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO IEQ=1,NEQ
       DD(iMEQ+IEQ) = DB(iMEQ+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(iMEQ+IEQ) = DD(iMEQ+IEQ) -
     *  AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
       END DO
      END DO
!
      CALL E011Sum(DD(iMEQ+1))
!      
      DO IEQ=1,NEQ
       DX(iMEQ+IEQ)=(1D0-OMEGA)*DX(iMEQ+IEQ)+
     * OMEGA*DD(iMEQ+IEQ)/MGE011(ILEV)%UE(iFLD)%x(IEQ)
      END DO
!     ----------------------------------------------
!     The real SOR part is here (going backwards)
!     ----------------------------------------------
      DO IEQ=NEQ,1,-1
       IF (KNPR(IEQ).EQ.0) THEN
        AUX=DB(iMEQ+IEQ)
        DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
         J=KCOL(ICOL)
         AUX = AUX - AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
        END DO
        DX(iMEQ+IEQ) = DX(iMEQ+IEQ) + 
     *  OMEGA*(AUX/AXX(iFLD)%FLD(ILEV)%a(KLD(IEQ))-DX(iMEQ+IEQ))
       END IF
      END DO
!     ----------------------------------------------
!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO IEQ=1,NEQ
       DD(iMEQ+IEQ) = DB(iMEQ+IEQ)
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(iMEQ+IEQ) = DD(iMEQ+IEQ) -
     *  AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
       END DO
      END DO
!
      CALL E011Sum(DD(iMEQ+1))
!
      DO IEQ=1,NEQ
       DX(iMEQ+IEQ)=(1D0-OMEGA)*DX(iMEQ+IEQ)+
     * OMEGA*DD(iMEQ+IEQ)/MGE011(ILEV)%UE(iFLD)%x(IEQ)
      END DO
!     ----------------------------------------------
C
1     CONTINUE
C
      END DO ! iFLD

      END
C
C
C
      SUBROUTINE SSORSolver_GenLinSc_Q1(AXX,KCOL,KLD,
     *           DX,DB,DD,KNPR,NEQ,NIT,NFLD,OMEGA,DEF)
C
      USE def_feat, ONLY: ILEV
      USE PP3D_MPI, ONLY : MGE011,myid,E011SUM
      USE var_QuadScalar, ONLY : tMGFldMatrix
      REAL*8 DX(*),DB(*),DD(*)
      INTEGER KCOL(*),KLD(*),KNPR(*)
      INTEGER NEQ,NIT,NFLD
      REAL*8 OMEGA,DEF
      TYPE(tMGFldMatrix) AXX(*)
      INTEGER iMEQ,IEQ,J,iFLD,ITE,ICOL
      REAL*8 AUX
!


      DO iFLD =1,NFLD

      iMEQ = (iFLD-1)*NEQ
!
      DO ITE=1,NIT
!     ----------------------------------------------
!     The real SOR part is here
!     ----------------------------------------------
      DO IEQ=1,NEQ
       IF (KNPR(IEQ).EQ.0) THEN
        AUX=DB(iMEQ+IEQ)
        DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
         J=KCOL(ICOL)
         AUX = AUX - AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
        END DO
        DX(iMEQ+IEQ) = DX(iMEQ+IEQ) + 
     *  OMEGA*(AUX/AXX(iFLD)%FLD(ILEV)%a(KLD(IEQ))-DX(iMEQ+IEQ))
       END IF
      END DO
!     ----------------------------------------------
!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO IEQ=1,NEQ
       DD(iMEQ+IEQ) = DB(iMEQ+IEQ)
!        DD(iMEQ+IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(iMEQ+IEQ) = DD(iMEQ+IEQ) -
     *  AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
       END DO
      END DO
!
      CALL E011Sum(DD(iMEQ+1))
!      
      DO IEQ=1,NEQ
       DX(iMEQ+IEQ)=(1D0-OMEGA)*DX(iMEQ+IEQ)+
     * OMEGA*DD(iMEQ+IEQ)/MGE011(ILEV)%UE(iFLD)%x(IEQ)
!      * OMEGA*(DD(iMEQ+IEQ) + 
!      * DB(iMEQ+IEQ))/MGE011(ILEV)%UE(iFLD)%x(IEQ)
      END DO
!     ----------------------------------------------
! !     The real SOR part is here (going backwards)
! !     ----------------------------------------------
      DO IEQ=NEQ,1,-1
       IF (KNPR(IEQ).EQ.0) THEN
        AUX=DB(iMEQ+IEQ)
        DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
         J=KCOL(ICOL)
         AUX = AUX - AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
        END DO
        DX(iMEQ+IEQ) = DX(iMEQ+IEQ) + 
     *  OMEGA*(AUX/AXX(iFLD)%FLD(ILEV)%a(KLD(IEQ))-DX(iMEQ+IEQ))
       END IF
      END DO
!     ----------------------------------------------
!     The parallel nodes are handeld by Jacobi ...
!     ----------------------------------------------
      DO IEQ=1,NEQ
       DD(iMEQ+IEQ) = DB(iMEQ+IEQ)
!        DD(iMEQ+IEQ) = 0d0
       DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
        J = KCOL(ICOL)
        DD(iMEQ+IEQ) = DD(iMEQ+IEQ) -
     *  AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
       END DO
      END DO
!
      CALL E011Sum(DD(iMEQ+1))
!      
      DO IEQ=1,NEQ
       DX(iMEQ+IEQ)=(1D0-OMEGA)*DX(iMEQ+IEQ)+
     * OMEGA*DD(iMEQ+IEQ)/MGE011(ILEV)%UE(iFLD)%x(IEQ)
!      * OMEGA*(DD(iMEQ+IEQ) + 
!      * DB(iMEQ+IEQ))/MGE011(ILEV)%UE(iFLD)%x(IEQ)
      END DO
!     ----------------------------------------------

      END DO ! ITE

      DO IEQ=1,NEQ
       DD(iMEQ+IEQ) = DB(iMEQ+IEQ)
       DO ICOL=KLD(IEQ),KLD(IEQ+1)-1
        J=KCOL(ICOL)
        DD(iMEQ+IEQ)=DD(iMEQ+IEQ)-
     *  AXX(iFLD)%FLD(ILEV)%a(ICOL)*DX(iMEQ+J)
       END DO
      END DO
C
      CALL E011Sum(DD(iMEQ+1))
      
      END DO !iFLD
     
      CALL LL21(DD,NFLD*NEQ,DEF)

      END
C
C
C
