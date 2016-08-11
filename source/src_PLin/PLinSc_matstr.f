************************************************************************
      SUBROUTINE   QP1_STRC_NUM(KADJ,NEL,NC)
************************************************************************
*    Purpose:  calculates the number of matrix entries
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      DIMENSION KADJ(NNAE,*)
C
C
C
      NC=0
C
      DO 10 IEL=1,NEL
      NC=NC+1
C
      DO 20 IAT=1,6
      IADJ=KADJ(IAT,IEL)
      IF (IADJ.EQ.0) GOTO 20
C
      NC=NC+1
20    CONTINUE
C
10    CONTINUE
C
C
C
      END
************************************************************************
      SUBROUTINE   QP1_STRC_ENT(KCOLCW,KLDCW,KCOLC,KLDC,KADJ,NEL,NC)
************************************************************************
*    Purpose:  calculates the matrix structure for Quadrilateral P1 el
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      DIMENSION KCOLC(*),KLDC(*),KCOLCW(*),KLDCW(*),KADJ(NNAE,*)
C
C
C
      NC=0
      KLDC(1)=1
C
      DO 10 IEL=1,NEL
      NC=NC+1
      KCOLC(NC)=IEL
C
      DO 20 IAT=1,6
      IADJ=KADJ(IAT,IEL)
      IF (IADJ.EQ.0) GOTO 20
C
      NC=NC+1
      KCOLC(NC)=IADJ
20    CONTINUE
C
      KLDC(IEL+1)=NC+1
C
10    CONTINUE
C
C
      DO 30 IEQ=1,NEL
C
31    BSORT=.TRUE.
      DO 32 ICOL=KLDC(IEQ)+1,KLDC(IEQ+1)-2
      IF (KCOLC(ICOL).GT.KCOLC(ICOL+1)) THEN
       IHELP=KCOLC(ICOL)
       KCOLC(ICOL)=KCOLC(ICOL+1)
       KCOLC(ICOL+1)=IHELP
       BSORT=.FALSE.
      ENDIF
32    CONTINUE
      IF (.NOT.BSORT) GOTO 31
C
30    CONTINUE
C
C
      ICOLW=0
      DO 40 IEQ=1,NEL
      DO 41 ICOL=KLDC(IEQ),KLDC(IEQ+1)-1
       KCOLCW(ICOLW+1) = KCOLC(ICOL)
       KCOLCW(ICOLW+2) = KCOLC(ICOL) +   NEL
       KCOLCW(ICOLW+3) = KCOLC(ICOL) + 2*NEL
       KCOLCW(ICOLW+4) = KCOLC(ICOL) + 3*NEL
       ICOLW = ICOLW + 4
41    CONTINUE
40    CONTINUE
C
      KLDCW(1)=1
      DO 50 IEQ=1,NEL
      IN = (KLDC(IEQ+1) - KLDC(IEQ))*4
      KLDCW(IEQ+1) = KLDCW(IEQ) + IN
50    CONTINUE
C
C
C
      END
