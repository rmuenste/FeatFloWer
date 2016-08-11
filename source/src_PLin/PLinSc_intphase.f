      SUBROUTINE GetSurf(IntPhase,DVAL,DMID,DCORVG,KVERT,KADJ,IO,NEL)
      USE PP3D_MPI, ONLY:myid
      USE DTP, ONLY: GetPolygon
      IMPLICIT NONE
      REAL*8  DCORVG(3,*),DVAL(*),DMID(3,*)
      INTEGER KVERT(8,*),NEL,IntPhase(*),KADJ(6,*),IO
      INTEGER IEL,IAT,JEL,IVT,IVERT,IP,JP,ILayer,NP
      REAL*8  DV(4),DX(3,8),DM(3),DG(3),DD,DS,DP(3,12)
      INTEGER I,J,K
      CHARACTER cp*(21)

      cp="#gmv/poly_          "

      IF(myid.lt.10) THEN
       WRITE(cp(11:13),'(I1,I1,A1)') 0,myid,'_'
      ELSE
        WRITE(cp(11:13),'(I2,A1)') myid,'_'
      END IF

      IF     ((iO.GE.0   ).AND.(iO.LT.10   )) THEN
       WRITE(cp(14:21),'(A3,I1,A4)') "000",iO,".gmv"
      ELSEIF ((iO.GE.10  ).AND.(iO.LT.100  )) THEN
       WRITE(cp(14:21),'(A2,I2,A4)') "00",iO,".gmv"
      ELSEIF ((iO.GE.100 ).AND.(iO.LT.1000 )) THEN
       WRITE(cp(14:21),'(A1,I3,A4)') "0",iO,".gmv"
      ELSEIF ((iO.GE.1000).AND.(iO.LT.10000)) THEN
       WRITE(cp(14:21),'(   I4,A4)') iO,".gmv"
      ELSEIF (iO.GE.10000)                    THEN
       STOP
      END IF

      OPEN (987,FILE=cp)
      WRITE(987,'(A)')'gmvinput ascii'
      WRITE(987,*)  'polygons'

      k=0
      DO IEL=1,NEL
       IF (IntPhase(IEL).EQ.0) THEN
        IP = 4*(IEL-1)
        DS=SQRT(DVAL(IP+2)**2d0+DVAL(IP+3)**2d0+DVAL(IP+4)**2d0)
        DV(1) = DVAL(IP+1)/DS
        DV(2) = DVAL(IP+2)/DS
        DV(3) = DVAL(IP+3)/DS
        DV(4) = DVAL(IP+4)/DS
        DM(1) = DMID(1,IEL)
        DM(2) = DMID(2,IEL)
        DM(3) = DMID(3,IEL)
        DO IVT=1,8
         IVERT = KVERT(IVT,IEL)
         DX(1,IVT) = DCORVG(1,IVERT)
         DX(2,IVT) = DCORVG(2,IVERT)
         DX(3,IVT) = DCORVG(3,IVERT)
        END DO

        CALL GetPolygon(DX,DV,DM,DP,NP)
        IF (NP.NE.0) THEN
         k=k+1
         WRITE(987,'(2i6,36G12.4)') 2,NP,((DP(I,J),J=1,NP),I=1,3)
        END IF

       END IF
      END DO ! IEL

      WRITE(987,*)  'endpoly'
      WRITE(987,*)  'endgmv'
      CLOSE(987)

      END
C
C
C
      SUBROUTINE SliceSetUp(SliceMap,DMID,KADJ,NEL,nS)
      USE PP3D_MPI, ONLY:myid
      IMPLICIT NONE
      INTEGER KADJ(6,*),SliceMap(*)
      INTEGER iS,jS,IAT,IEL,JEL,KEL,NEL,nS,IADJ
      REAL*8  DMZ,DEPS,DRAD,DMIN
      REAL*8  DMID(3,*)
      PARAMETER (DEPS = 0.0001d0)

!       if (myid.ne.0) WRITE(*,*) myid,DMID(3,1)

      iS = 0
!
      DMIN = 1d30
      DO IEL=1,NEL
       DRAD = DSQRT(DMID(1,IEL)**2d0+DMID(2,IEL)**2d0)
       IF (DRAD.LT.0.044d0.AND.DMID(3,IEL).LT.DMIN) THEN
        JEL = IEL
        DMIN = DMID(3,IEL)
       END IF
      END DO

      DO
       iS = iS + 1
       jS = 0
       DO KEL=1,NEL
        DRAD = DSQRT(DMID(1,KEL)**2d0+DMID(2,KEL)**2d0)
        IF (DRAD.LT.0.044d0.AND.DABS(DMID(3,JEL)-
     *                               DMID(3,KEL)).LT.DEPS) THEN
         jS = jS + 1
         SliceMap(KEL) = iS
        END IF
       END DO
       IEL = JEL
!        WRITE(*,*) myid,"|",iS,jS,IEL,DMID(3,IEL)
!        if (myid.eq.5) WRITE(*,*) iS,jS,IEL,DMID(3,IEL)
       DO IAT=1,6
        IADJ = KADJ(IAT,IEL)
        IF (IADJ.NE.0) THEN
         IF (DMID(3,IADJ)-DMID(3,IEL).GT.DEPS) THEN
!          IF (myid.eq.5) write(*,*) DMID(3,KADJ(IAT,IEL))-DMID(3,IEL)
         GOTO 1
         END IF
        END IF
       END DO
       GOTO 2

1      JEL = IADJ
      END DO

2     CONTINUE

      nS = iS
!       WRITE(*,*) myid,"|",nS
!       if (myid.eq.4) then
! !        DO IEL=1,NEL
!         WRITE(*,*) SliceMap(1:1000)
! !        END DO
!       end if

      END
C
C
C
      SUBROUTINE CorrectIntPhaseII(DVAL,AVOL,SliceMap,IPE,nS,NEL,BDS)
      USE PP3D_MPI, ONLY:myid,COMM_NLComplete,E012_CommSlice
      IMPLICIT NONE
      REAL*8 AVOL(*)
      REAL*8 DVAL(4,*),DAUX
      INTEGER SliceMap(*),IPE(*),nS,NEL
      INTEGER IEL,iS
      INTEGER IB1(nS+4),IB2(nS+4)
      LOGICAL BBU2,BBU1,BBM,BBL1,BBL2,BDS
      INTEGER IDS

!       RETURN
      IDS = 1
      IF (myid.ne.0) THEN

      IB1(1:nS+4) = 0
      IB2(1:nS+4) = 0
!       DO iS=1,nS
!        IB1(iS+2) = (8-myid)*100+iS
!        IB2(iS+2) = (8-myid)*100+iS
!       END DO

      DO IEL=1,NEL
       IF (SliceMap(IEL).NE.0) THEN
        iS = SliceMap(IEL)+2
        IF (IPE(IEL).LE. 0) IB1(iS) = IB1(iS) + 1
        IF (IPE(IEL).LE.-2) IB2(iS) = IB2(iS) + 1
       END IF
      END DO

      CALL E012_CommSlice(IB1,nS)
      CALL E012_CommSlice(IB2,nS)
!       write (*,'(I4,A5,80I4)') myid,"|1|",IB1
!       write (*,'(I4,A5,80I4)') myid,"|2|",IB2

      DO iS=3,nS+2
!      IF (myid.eq.1) WRITE(*,*) iS,iB1(iS),iB2(iS)
!        IF (iS.GT.4.AND.iS.LT.nS+1) THEN
        BBU2 = .FALSE.
        BBU1 = .FALSE.
        BBM  = .FALSE.
        BBL1 = .FALSE.
        BBL2 = .FALSE.
        IF (IB1(iS-2).GT.0.AND.IB2(iS-2).EQ.0) BBU2 = .TRUE.
        IF (IB1(iS-1).GT.0.AND.IB2(iS-1).EQ.0) BBU1 = .TRUE.
        IF (IB1(iS  ).GT.0.AND.IB2(iS  ).EQ.0) BBM  = .TRUE.
        IF (IB1(iS+1).GT.0.AND.IB2(iS+1).EQ.0) BBL1 = .TRUE.
        IF (IB1(iS+2).GT.0.AND.IB2(iS+2).EQ.0) BBL2 = .TRUE.
        IF (BBU2.AND.BBU1.AND.BBM.AND.BBL1.AND.BBL2) THEN
         IDS = 0
         DO IEL=1,NEL
          IF (SliceMap(IEL).EQ.iS-2) THEN
           IF (IPE(IEL).LE.0) THEN
            DAUX = DBLE(AVOL(IEL))**0.33333333d0
            DVAL(1,IEL) = DAUX
            DVAL(2,IEL) = 0d0
            DVAL(3,IEL) = 0d0
            DVAL(4,IEL) = 0d0
           END IF
          END IF
         END DO
        END IF
!        END IF
      END DO

      END IF

      CALL COMM_NLComplete(IDS)
      BDS=.FALSE.
      IF (IDS.EQ.0) BDS=.TRUE.

      END
C
C
C
      SUBROUTINE CorrectIntPhaseI(DF1,IPE,KVERT,DMID,NEL,DEPS)
      USE PP3D_MPI, ONLY:myid,COMM_NLComplete
      IMPLICIT NONE
      INTEGER IPE(*),KVERT(8,*),IEL,NEL,II,JJ,IVT,IVERT
      LOGICAL BB1,BB2
      INTEGER iCleanUp
      REAL*8 DMID(3,*),DF1(*),DEPS

      IF (myid.ne.0) THEN
!       II = 0
!       JJ = 0
!       DO IEL=1,NEL
!        IF (IPE(IEL).EQ.0) THEN
!         BB1 = .FALSE.
!         BB2 = .FALSE.
!         DO IVT=1,8
!           IVERT = KVERT(IVT,IEL)
!           IF (    DF1(IVERT).GT.DEPS) BB1 = .TRUE.
!           IF (1D0-DF1(IVERT).GT.DEPS) BB2 = .TRUE.
!         END DO
!         IF (.NOT.BB1) THEN
!          IPE(IEL) =  1
!          II = II + 1
!         END IF
!         IF (.NOT.BB2) THEN
!          IPE(IEL) = -1
!          JJ = JJ + 1
!         END IF
!        END IF
!       END DO

      iCleanUp = 1
      DO IEL=1,NEL
       IF (DMID(3,IEL).LT.-0.42.AND.DMID(3,IEL).GT.-0.425) THEN
        IF (ABS(IPE(IEL)).LE.4) THEN
         iCleanUp=0
!          write(*,*) iel,DMID(3,IEL)
        END IF
       END IF
      END DO

      END IF

      CALL COMM_NLComplete(iCleanUp)

      IF (myid.ne.0) THEN

      IF (iCleanUp) THEN
!        WRITE(*,*) myid,"Cleaning up the outflow subdomain ..."
       DO IEL=1,NEL
        IF (DMID(3,IEL).LT.-0.42) THEN
         IPE(IEL) = 100
        END IF
       END DO
      END IF

      END IF

      IF (myid.ne.0) THEN
!        IF(II.GE.1.OR.JJ.GE.1)WRITE(*,*)II,JJ,
!      *   " elements were reset on ",myid
      END IF

      END
C
C
C
      SUBROUTINE IntPhaseFinder(IntPhase,DVAL,DMID,DCORVG,
     *           KVERT,KAREA,KADJ,KVEL,NVEL,KMP,NEL)
      USE PP3D_MPI, ONLY:E012_CommInt,myid
      IMPLICIT NONE
      REAL*8  DCORVG(3,*),DVAL(*),DMID(3,*)
      INTEGER NEL,NVEL
      INTEGER KVERT(8,*),KAREA(6,*),KVEL(NVEL,*),IntPhase(*)
      INTEGER KADJ(6,*),KMP(5,*)
      INTEGER IEL,JEL,IAT,IVT,IVERT,IP,JP,IINT,JINT,KPAR,IAREA,IADJ
      REAL*8  DV,DX(3),DM(3),DMINDV,dReSet
      INTEGER iLayer,nLayer,IntphaseB(NEL),IntPhase_old(NEL)
      INTEGER iCount,jCount,jjCount,IVE,JVE,iFound,myElems(50),nElem
      PARAMETER (nLayer = 4)
      LOGICAL BB,BB1,BB2,bFound

      IntPhase_old(1:nel) = IntPhase(1:nel)

      DO IEL=1,NEL
       IP = 4*(IEL-1)
       DM(1)  = DMID(1,IEL)
       DM(2)  = DMID(2,IEL)
       DM(3)  = DMID(3,IEL)
       IINT = 0
       JINT = 0
       DO IVT=1,8
        IVERT = KVERT(IVT,IEL)
        DX(1) = DM(1)-DCORVG(1,IVERT)
        DX(2) = DM(2)-DCORVG(2,IVERT)
        DX(3) = DM(3)-DCORVG(3,IVERT)
        DV    = DVAL(IP+1)       + DVAL(IP+2)*DX(1) +
     *          DVAL(IP+3)*DX(2) + DVAL(IP+4)*DX(3)
        IF (DV.LE. 0.00075d0) IINT = IINT + 1
        IF (DV.GE.-0.00075d0) JINT = JINT + 1
       END DO
       IntPhase(IEL)=0
       IF (IINT.EQ.8.OR.JINT.EQ.8) THEN
        IF (IINT.EQ.8.AND.JINT.EQ.0) THEN
         IntPhase(IEL)=-100
        END IF
        IF (JINT.EQ.8.AND.IINT.EQ.0) THEN
         IntPhase(IEL)= 100
        END IF
       END IF
      END DO

      DO iLayer = 1,nLayer

      CALL GetParInts(IntPhase,IntPhaseB,KAREA,KMP,NEL)
      CALL E012_CommInt(IntPhaseB,KMP)

      KPAR = 1
      DO IEL=1,NEL
       DO IAT = 1,6
        IAREA= KAREA(IAT,IEL)
        IF (IAREA.EQ.KMP(5,KPAR)) THEN
         IF (ABS(IntPhaseB(KPAR)).EQ.iLayer-1) THEN
!          WRITE(*,*) myid,IntPhaseB(KPAR)
          IF (IntPhaseB(KPAR).EQ.0) THEN
            IF (IntPhase(IEL).LT.0) THEN
             IntPhase(IEL) = max(-1,IntPhase(IEL))
            END IF
            IF (IntPhase(IEL).GT.0) THEN
             IntPhase(IEL) = min(1,IntPhase(IEL))
            END IF
          END IF
          IF (IntPhaseB(KPAR).GT.0) THEN
           IntPhase(IEL) = min(IntPhaseB(KPAR)+1,IntPhase(IEL))
          END IF
          IF (IntPhaseB(KPAR).LT.0) THEN
           IntPhase(IEL) = max(IntPhaseB(KPAR)-1,IntPhase(IEL))
          END IF
         END IF
         KPAR = KPAR + 1
        END IF
       END DO
      END DO

      DO IEL=1,NEL
       IF (ABS(IntPhase(IEL)).EQ.iLayer-1) THEN
        DO IAT = 1,6
         JEL = KADJ(IAT,IEL)
         IF (JEL.NE.0.AND.ABS(IntPhase(JEL)).GT.iLayer) THEN
          IF (IntPhase(JEL).LT.0) IntPhase(JEL) = -iLayer
          IF (IntPhase(JEL).GT.0) IntPhase(JEL) = +iLayer
         END IF
        END DO
       END IF
      END DO

      END DO
C
      END
C
C
C
      SUBROUTINE AddElem(jel,elems,nel)
      INTEGER iel,jel,nel,elems(*)

      DO iel=1,nel
       IF (elems(iel).EQ.jel) GOTO 1
      END DO

      nel = nel + 1
      elems(nel) = jel

1     CONTINUE

      END
C
C
C
      SUBROUTINE ReinitFarField(IntPhase,dMid,DVAL,NEL)
      IMPLICIT NONE
      REAL*8  DVAL(4,*),dMid(3,*)
      INTEGER NEL,IntPhase(*)
      REAL*8  DV(4),DAUX,DS,dFar,dmz
      INTEGER IEL,IP

      DO IEL=1,NEL
       DMZ = DMID(3,IEL)
       dFar = 0.04d0
       IF (ABS(IntPhase(IEL)).EQ.100) THEN
        DS = DSIGN(1d0,DBLE(IntPhase(IEL)))
        DVAL(1,IEL) = DS*dFar
        DVAL(2,IEL) = 0d0
        DVAL(3,IEL) = 0d0
        DVAL(4,IEL) = 0d0
       END IF
      END DO

      END
C
C
C
      SUBROUTINE ReinitCloseField(IntPhase,DVAL,NEL)
      IMPLICIT NONE
      REAL*8  DVAL(4,*)
      INTEGER NEL,IntPhase(*)
      REAL*8  DV(4),DAUX,DS
      INTEGER IEL,IP

!       RETURN
      DO IEL=1,NEL
       IF (IntPhase(IEL).EQ.0) THEN
        DV(1) = DVAL(1,IEL)
        DV(2) = DVAL(2,IEL)
        DV(3) = DVAL(3,IEL)
        DV(4) = DVAL(4,IEL)
        DAUX = SQRT(DV(2)*DV(2)+DV(3)*DV(3)+DV(4)*DV(4))
        DVAL(1,IEL) = DV(1)/DAUX
        DVAL(2,IEL) = DV(2)/DAUX
        DVAL(3,IEL) = DV(3)/DAUX
        DVAL(4,IEL) = DV(4)/DAUX
       END IF
      END DO

      END
C
C
C
      SUBROUTINE ReinitWholeField(IntPhase,DVAL,NEL)
      IMPLICIT NONE
      REAL*8  DVAL(4,*)
      INTEGER NEL,IntPhase(*)
      REAL*8  DV(4),DAUX,DS
      INTEGER IEL,IP

!       RETURN
      DO IEL=1,NEL
       IF (ABS(IntPhase(IEL)).NE.100) THEN
        DV(1) = DVAL(1,IEL)
        DV(2) = DVAL(2,IEL)
        DV(3) = DVAL(3,IEL)
        DV(4) = DVAL(4,IEL)
        DAUX = SQRT(DV(2)*DV(2)+DV(3)*DV(3)+DV(4)*DV(4))
        IF (DAUX.GT.1d-4) THEN
         DVAL(1,IEL) = DV(1)/DAUX
         DVAL(2,IEL) = DV(2)/DAUX
         DVAL(3,IEL) = DV(3)/DAUX
         DVAL(4,IEL) = DV(4)/DAUX
        END IF
       END IF
      END DO

      END
C
C
C
      SUBROUTINE ReinitInitField(IntPhase,DVAL,NEL)
      IMPLICIT NONE
      REAL*8  DVAL(4,*)
      INTEGER NEL,IntPhase(*)
      REAL*8  DV(4),DAUX,DS
      INTEGER IEL,IP

!       RETURN
      DO IEL=1,NEL
       IF (IntPhase(IEL).NE.0) THEN
         DVAL(1,IEL) = 0d0
         DVAL(2,IEL) = 0d0
         DVAL(3,IEL) = 0d0
         DVAL(4,IEL) = 0d0
       END IF
      END DO

      END
C
C
C
      SUBROUTINE GetElementMids(DMP,DCORVG,KVERT,NEL)
      IMPLICIT NONE
      REAL*8 DCORVG(3,*),DMP(3,*)
      INTEGER KVERT(8,*),NEL
      INTEGER IEL,IVT,IVERT

      DO IEL=1,NEL
       DMP(1,IEL) = 0d0
       DMP(2,IEL) = 0d0
       DMP(3,IEL) = 0d0
       DO IVT=1,8
        IVERT = KVERT(IVT,IEL)
        DMP(1,IEL) = DMP(1,IEL) + 0.125d0*DCORVG(1,IVERT)
        DMP(2,IEL) = DMP(2,IEL) + 0.125d0*DCORVG(2,IVERT)
        DMP(3,IEL) = DMP(3,IEL) + 0.125d0*DCORVG(3,IVERT)
       END DO
      END DO

      END
C
C
C
      SUBROUTINE GetAverageGradPhi(DVAL,NEL,IPE,DGRAD)
      USE PP3D_MPI, ONLY: Comm_Summ,myid
      IMPLICIT NONE
      REAL*8 DVAL(4,*),DGRAD
      INTEGER NEL,IPE(*)
      REAL*8 DD,DAUX
      INTEGER IEL,JEL

      IF (myid.ne.0) THEN
       JEL = 0
       DGRAD = 0d0
       DAUX=0d0
       DO IEL=1,NEL
        IF (ABS(IPE(IEL)).LT.4) THEN
         JEL = JEL + 1
         DD = DSQRT(DVAL(2,IEL)**2d0+DVAL(3,IEL)**2d0+DVAL(4,IEL)**2d0)
         DD = DABS(1d0-DD)
         DGRAD = DGRAD + DD
        ELSE
         IF (ABS(IPE(IEL)).LT.100) THEN
         DD = DSQRT(DVAL(2,IEL)**2d0+DVAL(3,IEL)**2d0+DVAL(4,IEL)**2d0)
         DD = DABS(1d0-DD)
         DAUX = MAX(DAUX,DD)
         END IF
        END IF
       END DO
      END IF

      DD = DBLE(JEL)
!       IF (DAUX.GT.2d0) DGRAD = 1D33
      CALL Comm_Summ(DGRAD)
      CALL Comm_Summ(DD)

      DGRAD = DGRAD/DD

      END SUBROUTINE GetAverageGradPhi

      SUBROUTINE CheckGradPhiQualitySub(DVAL,IPE,NEL,DQ)
      USE PP3D_MPI, ONLY: Comm_Summ,myid
      INTEGER IPE(*),NEL
      REAL*8 DVAL(4,*)
      INTEGER IEL
      REAL*8 DD,DQ

      DQ = 0d0
      DO IEL=1,NEL
       IF (IPE(IEL).GT.0.AND.IPE(IEL).LT.4) THEN
        DD = DSQRT(DVAL(2,IEL)**2d0+DVAL(3,IEL)**2d0+DVAL(4,IEL)**2d0)
        IF (ABS(DD-1d0).GT.0.025d0) DQ= DQ + 1d0
       END IF
      END DO
      write(*,*) dq

      END

