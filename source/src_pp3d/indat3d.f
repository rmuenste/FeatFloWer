************************************************************************
      DOUBLE PRECISION FUNCTION FDATIN(ITYP,IBLOC,X,Y,Z,TIMENS,RE)
*
*     Prescribed data for files coeff.f and bndry.f
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (PI=3.1415926535897931D0)
C
C
C
      FDATIN=0D0
C
C
C
C=======================================================================
C *** Case 1: Velocity boundary values and/or exact solution
C=======================================================================
C
      IF (ITYP.EQ.1) THEN

       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        IF (Z.LT.1d-4) THEN
         DIST = SQRT(X**2d0+Y**2D0)
         IF ((DIST-0.25d0).LT.1d-4) THEN
          FDATIN=(0.25d0-DIST)*(DIST+0.25d0)*16d0*4/3
         ELSE
          FDATIN=(0.50d0-DIST)*(DIST-0.25d0)*64d0*4/3
         END IF
        END IF
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 2: Velocity x-derivative of exact solution
C=======================================================================
C
      IF (ITYP.EQ.2) THEN
C
       IF (IBLOC.EQ.1) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 3: Velocity y-derivative of exact solution
C=======================================================================
C
      IF (ITYP.EQ.3) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 4: Velocity z-derivative of exact solution
C=======================================================================
C
      IF (ITYP.EQ.4) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 5: Exact pressure solution
C=======================================================================
C
      IF (ITYP.EQ.5) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C=======================================================================
C *** Case 6: Right hand side for momentum equation
C=======================================================================
C
      IF (ITYP.EQ.6) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 7: Right hand side for continuity equation
C=======================================================================
C
      IF (ITYP.EQ.7) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C=======================================================================
C *** Case 8: Mean pressure values
C=======================================================================
C
      IF (ITYP.EQ.8) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE NEUDAT(PX,PY,PZ,TIMENS,IFLAG,JVERT)
*
*     Neumann-boundary part
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C     
C=======================================================================
C *** Set Neumann-boundary parts
C         IFLAG =  1   Dirichlet BC for all components
C         IFLAG =  0   Neumann BC for all components
C         IFLAG = -1   No-penetration boundary condition
C         IFLAG = -2   Symmetry boundary condition
C         IFLAG =  3   Periodic boundary condition
C=======================================================================
C
!     Dirichlet everywhere
      IFLAG=1
C
!     Outflow
      IF ((1.0d0-PZ).LT.1d-4) THEN
       DIST = 0.5d0-SQRT(PX**2d0+PY**2D0)
       IF (.NOT.(DIST.LT.1d-4)) THEN
         IFLAG=0
       END IF
      END IF
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE PERIODa(I,PI,PJ)
************************************************************************
*     Purpose: Associate the midpoint pairs for the PERIODIC BC
*-----------------------------------------------------------------------
C
      INTEGER I
      REAL*8 PI(3),PJ(3),tol
      DATA tol/5d-5/
C
      I = 0
      IF ((PI(2).GT.(PJ(2)-tol).AND.PI(2).LT.(PJ(2)+tol)).AND.
     *    (PI(3).GT.(PJ(3)-tol).AND.PI(3).LT.(PJ(3)+tol))) THEN
       IF (PI(1).GT.0.99999D0) I = -1
       IF (PJ(1).GT.0.99999D0) I = +1
      END IF

      IF ((PI(1).GT.(PJ(1)-tol).AND.PI(1).LT.(PJ(1)+tol)).AND.
     *    (PI(3).GT.(PJ(3)-tol).AND.PI(3).LT.(PJ(3)+tol))) THEN
       IF (PI(2).LT.0.00001D0) I = +1
       IF (PJ(2).LT.0.00001D0) I = -1
       RETURN
      END IF

      END
C
C
C
************************************************************************
      SUBROUTINE BDPDAT(IEL,INPR,PX,PY,PZ,TIMENS,IFLAG1,IFLAG2)
*
*     Pressure integral boundary part
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C=======================================================================
C *** Set pressure integral boundary parts
C=======================================================================
C
C      IF (((PZ.GT.0.00D0).AND.(PZ.LT.0.41D0)).AND.
C     *    ((PX.GE.0.45D0).AND.(PX.LE.0.55D0)).AND.
C     *    ((PY.GE.0.15D0).AND.(PY.LE.0.25D0))) THEN
C       IFLAG1=1      
C      ENDIF
C
C      IF (((PZ.GT.0.00D0).AND.(PZ.LT.0.41D0)).AND.
C     *    ((PX.GE.0.45D0).AND.(PX.LE.0.50D0)).AND.
C     *    ((PY.GE.0.15D0).AND.(PY.LE.0.25D0))) THEN
C       IFLAG2=1
C      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE BDFDAT(IEL,INPR,PX,PY,PZ,TIMENS,DNU,IFLAG,DPF1,DPF2)
*
*     lift and drag data
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C     
C=======================================================================
C *** Parameters for lift (DFW) and drag (DAW) in bdforc.f (INPR=2)
C ***
C *** dfw=2 int_s [dpf(1) dut/dn n_y - p n_x] ds / dpf(2)
C *** daw=2 int_s [dpf(1) dut/dn n_x + p n_y] ds / dpf(2)
C ***
C=======================================================================
C
C      RHO  =1.0D0
C      DIST =0.041D0
C      UMEAN=1.0D0
C
C      DPF1=RHO*DNU
C      DPF2=RHO*DIST*UMEAN**2
C
C=======================================================================
C
C      IF (((PZ.GT.0.00D0).AND.(PZ.LT.0.41D0)).AND.
C     *    ((PX.GE.0.45D0).AND.(PX.LE.0.55D0)).AND.
C     *    ((PY.GE.0.15D0).AND.(PY.LE.0.25D0))) THEN
C       IFLAG=1      
C      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE PTSDAT(TIMENS,DNU)
*
*     Data for Point-output (for fpost)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /NSPTS/  KPU(2),KPP(4)
      SAVE
C
C     
C
C=======================================================================
C *** Points for velocity and pressure
C=======================================================================
C
C      KPU(1)=3421
C      KPU(2)=677
C
C      KPP(1)=687
C      KPP(2)=690
C      KPP(3)=3498
C      KPP(4)=677
C
C
C
99999 END
