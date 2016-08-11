      SUBROUTINE GetParValues(DVAL,DBP,KAREA,KBP,NEL)
      IMPLICIT NONE
      REAL*8 DVAL(*),DBP(4,*)
      INTEGER KAREA(6,*),KBP(5,*),NEL
      INTEGER IEL,IAT,IAREA,K,L

      K=1
      DO IEL=1,NEL
      DO IAT=1,6
       IAREA = KAREA(IAT,IEL)
       IF (IAREA.EQ.KBP(5,K)) THEN
        L = 4*(IEL-1)
        DBP(1,K) = DVAL(L+1)
        DBP(2,K) = DVAL(L+2)
        DBP(3,K) = DVAL(L+3)
        DBP(4,K) = DVAL(L+4)
        K = K + 1
       END IF
      END DO
      END DO

      END
C
      SUBROUTINE GetParInts(IVAL,IBP,KAREA,KBP,NEL)
      IMPLICIT NONE
      INTEGER IVAL(*),IBP(*)
      INTEGER KAREA(6,*),KBP(5,*),NEL
      INTEGER IEL,IAT,IAREA,K,L

      K=1
      DO IEL=1,NEL
      DO IAT=1,6
       IAREA = KAREA(IAT,IEL)
       IF (IAREA.EQ.KBP(5,K)) THEN
        IBP(K) = IVAL(IEL)
        K = K + 1
       END IF
      END DO
      END DO

      END
C
      SUBROUTINE GetParMidValuesO(DCORAG,DMP,KAREA,KVERT,KBP,NEL)
      IMPLICIT NONE
      REAL*8 DCORAG(3,*),DMP(4,*)
      INTEGER KAREA(6,*),KVERT(8,*),KBP(5,*),NEL
      INTEGER IEL,IVT,IAT,IAREA,IVERT,K

      K=1
      DO IEL=1,NEL
      DO IAT=1,6
       IAREA = KAREA(IAT,IEL)
       IF (IAREA.EQ.KBP(5,K)) THEN
        DMP(1,K) = DCORAG(1,IAREA)
        DMP(2,K) = DCORAG(2,IAREA)
        DMP(3,K) = DCORAG(3,IAREA)
        DMP(4,K) = 0d0
        K = K + 1
       END IF
      END DO
      END DO

      END
C
      SUBROUTINE GetParMidValues(DCORVG,DMP,KAREA,KVERT,KBP,NEL)
      IMPLICIT NONE
      REAL*8 DCORVG(3,*),DMP(4,*)
      INTEGER KAREA(6,*),KVERT(8,*),KBP(5,*),NEL
      INTEGER IEL,IVT,IAT,IAREA,IVERT,K

      K=1
      DO IEL=1,NEL
      DO IAT=1,6
       IAREA = KAREA(IAT,IEL)
       IF (IAREA.EQ.KBP(5,K)) THEN
        DMP(1,K) = 0d0
        DMP(2,K) = 0d0
        DMP(3,K) = 0d0
        DMP(4,K) = 0d0
        DO IVT=1,8
         IVERT = KVERT(IVT,IEL)
         DMP(1,K) = DMP(1,K) + 0.125d0*DCORVG(1,IVERT)
         DMP(2,K) = DMP(2,K) + 0.125d0*DCORVG(2,IVERT)
         DMP(3,K) = DMP(3,K) + 0.125d0*DCORVG(3,IVERT)
        END DO
        K = K + 1
       END IF
      END DO
      END DO

      END
