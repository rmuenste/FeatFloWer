************************************************************************
      SUBROUTINE XOPTCN(DU1,DU2,DU3,DU1OLD,DU2OLD,DU3OLD,NU,
     *                  DELU1,DELU2,DELU3,OMEGA)
************************************************************************
*    Purpose: - Computes for given vectors U and UOLD the optimal 
*               weighted correction OMEGA*V with V:=U-UOLD
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      DIMENSION DU1(*),DU2(*),DU3(*),DU1OLD(*),DU2OLD(*),DU3OLD(*)
C
      SAVE 
C
C=======================================================================
C *** Calculate the optimal correction  U:=OMEGA*(U-UOLD)
C=======================================================================
C
      CALL LSC1(DU1,NU,OMEGA)
      CALL LSC1(DU2,NU,OMEGA)
      CALL LSC1(DU3,NU,OMEGA)
C
C=======================================================================
C *** Calculate maximum changes   DELU
C=======================================================================
C
      CALL LLI1(DU1,NU,DELU1,INDU1)
      CALL LLI1(DU2,NU,DELU2,INDU2)
      CALL LLI1(DU3,NU,DELU3,INDU3)
C
C=======================================================================
C *** Update the solution   U:=UOLD+U
C=======================================================================
C
      CALL LLC1(DU1OLD,DU1,NU,1D0,1D0)
      CALL LLC1(DU2OLD,DU2,NU,1D0,1D0)
      CALL LLC1(DU3OLD,DU3,NU,1D0,1D0)
C
C=======================================================================
C *** relative maximum changes   DELU
C=======================================================================
C
      CALL LLI1(DU1,NU,DELT1,INDT1)
      CALL LLI1(DU2,NU,DELT2,INDT2)
      CALL LLI1(DU3,NU,DELT3,INDT3)
C
      DELT=MAX(1D-12,DELT1,DELT2,DELT3)
      DELU1=DELU1/DELT
      DELU2=DELU2/DELT
      DELU3=DELU3/DELT
C
C
99999 END
