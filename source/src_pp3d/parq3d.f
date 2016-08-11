      DOUBLE PRECISION FUNCTION PARX(T1,T2,T3,IBCT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARX=T1
99999 END
C
C
C
      DOUBLE PRECISION FUNCTION PARY(T1,T2,T3,IBCT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARY=T2
99999 END
C
C
C
      DOUBLE PRECISION FUNCTION PARZ(T1,T2,T3,IBCT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARZ=T3
      GOTO 99999
99999 END
C
C
C
************************************************************************
      SUBROUTINE   EXTRUDER(IOUT,DCORVG,KNPR,KVEL,DCORV2,KNPR2,BD,BD_Z,
     *             iBD_Z,NVT,NVEL,ILEV,NLEV)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNLEV=9)
      CHARACTER CFILE*(60),CTEXT*(8) 
      REAL*8 DCORVG(3,*),DCORV2(2,*),BD(2,*),BD_Z(2,*),MINDIST
      INTEGER KNPR2(*),KVEL(NVEL,*),KNPR(*),iBD_Z(2,*),MFILE
      CHARACTER CTRIB(NNLEV)*15
C
      DATA MFILE /91/ 
      DATA CTRIB/ '#ext/par2D1.ext','#ext/par2D2.ext','#ext/par2D3.ext',
     *            '#ext/par2D4.ext','#ext/par2D5.ext','#ext/par2D6.ext',
     *            '#ext/par2D7.ext','#ext/par2D8.ext','#ext/par2D9.ext'/

      REAL*8 MIN_Z,MAX_Z
C-----------------------------------------------------------------------
      RETURN
      IF ((IOUT.EQ.1).AND.(ILEV.EQ.1)) WRITE(*,102)
C
      MIN_Z = 1D30
      MAX_Z = -1D30
      DO IVT=1,NVT
       IF (MIN_Z.GT.DCORVG(3,IVT)) MIN_Z = DCORVG(3,IVT)
       IF (MAX_Z.LT.DCORVG(3,IVT)) MAX_Z = DCORVG(3,IVT)
      END DO
C-----------------------------------------------------------------------
      OPEN(UNIT=MFILE,FILE=CTRIB(ILEV))
      READ(MFILE,*) 
      READ(MFILE,*) NVT2
      READ(MFILE,*) 
      READ(MFILE,*) (DCORV2(1,IVT),IVT=1,NVT2)
      READ(MFILE,*) (DCORV2(2,IVT),IVT=1,NVT2)
      READ(MFILE,*) 
      iCOUNT = 0
      DO IVT=1,NVT2
       READ(MFILE,*) KNPR2(IVT)
       IF (KNPR2(IVT).NE.0) iCOUNT = iCOUNT + 1
      END DO
      DO IVT=1,NVT2
       iBD_Z(1,IVT) = 0
       BD_Z(1,IVT) = 1D30
       iBD_Z(2,IVT) = 0
       BD_Z(2,IVT) = 1D30
      END DO
      iCOUNT = 0
      DO IVT=1,NVT2
       IF (KNPR2(IVT).NE.0) THEN
        iCOUNT = iCOUNT + 1
        BD(1,iCOUNT) = DCORV2(1,IVT)
        BD(2,iCOUNT) = DCORV2(2,IVT)
       END IF
      END DO
      CLOSE(MFILE)
      KBD = 0
C-----------------------------------------------------------------------
      DO 10 IVT=1,NVT
       INPR = KNPR(IVT)
       IF (INPR.EQ.0) GOTO 10
       X=DCORVG(1,IVT)
       Y=DCORVG(2,IVT)
       Z=DCORVG(3,IVT)
       MINDIST = 1D30
       DO JVT=1,iCOUNT
        DIST = DSQRT((DCORVG(1,IVT)-BD(1,JVT))**2+
     *               (DCORVG(2,IVT)-BD(2,JVT))**2)
        IF (DIST.LT.MINDIST) THEN 
         MINDIST = DIST
         JBD = JVT
        END IF
       END DO
       IF ((Z.EQ.MIN_Z).OR.(Z.EQ.MAX_Z)) THEN
        IF (KVEL(3,IVT).NE.0) GOTO 10
        IF (Z.EQ.MIN_Z) iPOS = 1
        IF (Z.EQ.MAX_Z) iPOS = 2
        IF (BD_Z(iPOS,JBD).GT.MINDIST) THEN
         BD_Z(iPOS,JBD) = MINDIST
         iBD_Z(iPOS,JBD) = IVT
        END IF
       ELSE
        DCORVG(1,IVT) = BD(1,JBD)
        DCORVG(2,IVT) = BD(2,JBD)
        KBD = KBD + 1
       END IF
10    CONTINUE
C
      DO JBD=1,iCOUNT
       IF ((iBD_Z(1,JBD).LE.NVT).AND.(iBD_Z(1,JBD).GE.1).AND.
     *     (iBD_Z(2,JBD).LE.NVT).AND.(iBD_Z(2,JBD).GE.1)) THEN
        DCORVG(1,iBD_Z(1,JBD)) = BD(1,JBD)
        DCORVG(2,iBD_Z(1,JBD)) = BD(2,JBD)
        KBD = KBD + 1
	DCORVG(1,iBD_Z(2,JBD)) = BD(1,JBD)
	DCORVG(2,iBD_Z(2,JBD)) = BD(2,JBD)
 	KBD = KBD + 1
       END IF
      END DO 
C-----------------------------------------------------------------------
      IF (IOUT.EQ.1) WRITE(*,101) iCOUNT,KBD
C-----------------------------------------------------------------------
      IF ((IOUT.EQ.1).AND.(ILEV.EQ.NLEV)) WRITE(*,103)
C
101   FORMAT('2D BND POINTS: ',I7,'  3D BND POINTS: ',I7)
102   FORMAT(17('-'),(' NUMBERS OF MANIPULATED VERTICES IN 2D AND 3D '),
     *       17('-'))
103   FORMAT(80('-'))
      END

      SUBROUTINE MEM_EXT(iLOCATE,NVT2,LCOR2,LKPR2,LDBD,LDBDZ,LiDBZ)
      INTEGER iLOCATE,NVT2,LCOR2,LKPR2,LDBD,LDBDZ,LiDBZ
C
      IF (iLOCATE.EQ.0) THEN
       CALL ZNEW(2*NVT2,1,LCOR2,'COR2  ')
       CALL ZNEW(  NVT2,3,LKPR2,'KPR2  ')
       CALL ZNEW(2*NVT2,1,LDBD, 'DBD   ')
       CALL ZNEW(2*NVT2,1,LDBDZ,'DBDZ  ')
       CALL ZNEW(2*NVT2,3,LiDBZ,'iDBZ  ')
      ELSE
       CALL ZDISP(0,LCOR2,'COR2  ')
       CALL ZDISP(0,LKPR2,'KPR2  ')
       CALL ZDISP(0,LDBD, 'DBD   ')
       CALL ZDISP(0,LDBDZ,'DBDZ  ')
       CALL ZDISP(0,LiDBZ,'iDBZ  ')
      END IF
C
      END 
C
C
C
      SUBROUTINE DEFORM(DCORVG,NVT,II)
!      USE PP3D_MPI, ONLY : myid
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*8 DCORVG(3,*), DX, DY, DZ
      REAL*8 RX,RY,RADMAX,RADMIN,RATIO
      DATA RX/0.5D0/,RY/0.2D0/,RAD/0.05d0/
      DATA RADMIN/0.011D0/,RADMAX/0.046D0/

      !return
      IF (II.EQ.2) DCRIT = 0.001d0
      IF (II.EQ.3) DCRIT = 0.0006d0
      IF (II.EQ.4) DCRIT = 0.00035d0

      DO IVT=1,NVT
          DX=DCORVG(1,IVT)-RX
          DY=DCORVG(2,IVT)-RY
          DRAD = SQRT(DX*DX+DY*DY)
          IF (ABS(DRAD-RAD).LT.DCRIT) THEN
           DFACT = RAD/DRAD
           DCORVG(1,IVT) = RX + DFACT*DX
           DCORVG(2,IVT) = RY + DFACT*DY
        END IF

      END DO





!!!!!!!!!!!!!!!!!!!!!!!!original!!!!!!!!!!!!!!!!!!!!!
!      IF (II.EQ.2) DCRIT = 0.0075d0
!      IF (II.EQ.3) DCRIT = 0.0025d0
!      IF (II.EQ.4) DCRIT = 0.0010d0

!      DO IVT=1,NVT
!        DX=DCORVG(1,IVT)
!        DY=DCORVG(2,IVT)
!        IF (DY.LT.-5d-4) THEN
!         DIST = SQRT((DX-0d0)**2d0+(DY+0.125d0)**2d0)
!         IF (DIST.LT.0.075d0) THEN
!          dnx = dx-0.000d0
!          dny = dy+0.125d0
!          DCORVG(1,IVT) = dnx*0.075d0/DIST + 0.000d0
!          DCORVG(2,IVT) = dny*0.075d0/DIST - 0.125d0
!         END IF
!         IF (0.15d0-DIST.LT.DCRIT) THEN
!          dnx = dx-0.000d0
!          dny = dy+0.125d0
!          DCORVG(1,IVT) = dnx*0.15d0/DIST + 0.000d0
!          DCORVG(2,IVT) = dny*0.15d0/DIST - 0.125d0
!         END IF
!        END IF
!        IF (DY.GT.+5d-4) THEN
!         DIST = SQRT((DX-0d0)**2d0+(DY-0.125d0)**2d0)
!         IF (DIST.LT.0.075d0) THEN
!          dnx = dx-0.000d0
!          dny = dy-0.125d0
!          DCORVG(1,IVT) = dnx*0.075d0/DIST + 0.000d0
!          DCORVG(2,IVT) = dny*0.075d0/DIST + 0.125d0
!         END IF
!         IF (0.15d0-DIST.LT.DCRIT) THEN
!          dnx = dx-0.000d0
!          dny = dy-0.125d0
!          DCORVG(1,IVT) = dnx*0.15d0/DIST + 0.000d0
!          DCORVG(2,IVT) = dny*0.15d0/DIST + 0.125d0
!         END IF
!        END IF

!      END DO

      RETURN

      END
C
C
C
      SUBROUTINE INPRSET(DCORVG,KNPR,NVT,C)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNLEV=9)
      CHARACTER CFILE*(60),CTEXT*(8) 
      REAL*8 DCORVG(3,*),RX,RY
      INTEGER KNPR(*)
      CHARACTER C*1
      DATA RX/0.1D0/,RY/0.1D0/, DISTREF/0.011D0/

      RETURN
      J = 0
      DO I = 1,NVT
!       WRITE(*,*) DCORVG(1,I),DCORVG(2,I),DCORVG(3,I)
       IF ((DCORVG(3,I).GT.0.0D0).AND.(DCORVG(3,I).LT.0.0625D0).AND. !0.51D-2
     *     (DCORVG(1,I).GT.1D-5   ).AND.(DCORVG(1,I).LT.1-1D-5 ).AND.
     *     (DCORVG(2,I).GT.1D-5   ).AND.(DCORVG(2,I).LT.1-1D-5 )) THEN
           J = J + 1
           KNPR(I) = 0
       END IF
      END DO

      WRITE(*,*) "number of changes due to INPRSET:", J," for ",C
      END
C
C
C

      SUBROUTINE MIDDEF(DCORVG,KVERT,KAREA,NEL,DCORAG)
      INTEGER KAREA(6,*),KVERT(8,*)
      REAL*8 DCORVG(3,*),DCORAG(3,*)
      REAL*8 P1X,P2X,P3X,P4X,P1Y,P2Y,P3Y,P4Y,P1Z,P2Z,P3Z,P4Z
      REAL*8 PX,PY,PZ
      INTEGER IEL,IAT,NEL,II
      SAVE

      !RETURN
      II = 0
      DO 1 IEL=1,NEL
C
       DO 2 IAT=1,6
        IAR =KAREA(IAT,IEL)
C
        IF (IAT.EQ.1) THEN
         IVT1=KVERT(1,IEL)
         IVT2=KVERT(2,IEL)
         IVT3=KVERT(3,IEL)
         IVT4=KVERT(4,IEL) 
         GOTO 3
        ENDIF
C
        IF (IAT.EQ.2) THEN
         IVT1=KVERT(1,IEL)
         IVT2=KVERT(2,IEL)
         IVT3=KVERT(6,IEL)
         IVT4=KVERT(5,IEL)
         GOTO 3
        ENDIF
C
        IF (IAT.EQ.3) THEN
         IVT1=KVERT(2,IEL)
         IVT2=KVERT(3,IEL)
         IVT3=KVERT(7,IEL)
         IVT4=KVERT(6,IEL)
         GOTO 3
        ENDIF
C
        IF (IAT.EQ.4) THEN
         IVT1=KVERT(3,IEL)
         IVT2=KVERT(4,IEL)
         IVT3=KVERT(8,IEL)
         IVT4=KVERT(7,IEL)
         GOTO 3
        ENDIF
C
        IF (IAT.EQ.5) THEN
         IVT1=KVERT(4,IEL)
         IVT2=KVERT(1,IEL)
         IVT3=KVERT(5,IEL)
         IVT4=KVERT(8,IEL)
         GOTO 3
        ENDIF
C
        IF (IAT.EQ.6) THEN
         IVT1=KVERT(5,IEL)
         IVT2=KVERT(6,IEL)
         IVT3=KVERT(7,IEL)
         IVT4=KVERT(8,IEL)
         GOTO 3
        ENDIF
C
3       P1X=DCORVG(1,IVT1)
        P1Y=DCORVG(2,IVT1)
        P1Z=DCORVG(3,IVT1)
        P2X=DCORVG(1,IVT2)
        P2Y=DCORVG(2,IVT2)
        P2Z=DCORVG(3,IVT2)
        P3X=DCORVG(1,IVT3)
        P3Y=DCORVG(2,IVT3)
        P3Z=DCORVG(3,IVT3)
        P4X=DCORVG(1,IVT4)
        P4Y=DCORVG(2,IVT4)
        P4Z=DCORVG(3,IVT4)
C
        PX=(P1X+P2X+P3X+P4X)*0.25D0
        PY=(P1Y+P2Y+P3Y+P4Y)*0.25D0
        PZ=(P1Z+P2Z+P3Z+P4Z)*0.25D0
C
        DCORAG(1,IAR) = PX
        DCORAG(2,IAR) = PY
        DCORAG(3,IAR) = PZ
C
2      CONTINUE
1     CONTINUE
C
!       WRITE(*,*)
!      *"             Midpoint translation due to 'DISTORSION' ..."
!       WRITE(*,103)
C
C
C
C
103   FORMAT(80('-'))
      END
