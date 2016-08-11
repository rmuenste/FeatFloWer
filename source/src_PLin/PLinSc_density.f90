      SUBROUTINE IntegrateMass(DMASS,DRHO_ELEM,VVOL,NEL)
      REAL*4  VVOL(*)
      REAL*8  DRHO_ELEM(*),DMASS
      INTEGER NEL
      INTEGER IEL

      DMASS = 0d0
      DO IEL=1,NEL
       DMASS = DMASS + DRHO_ELEM(IEL)*DBLE(VVOL(IEL))
      END DO

      END
!
!
!
      SUBROUTINE GetP1Density(DRHO_ELEM,IntPhase,DVAL,DMID,DFRAC,&
                 DCORVG,KVERT,VVOL,DENS1,DENS2,DFRAC1,DFRAC2,NEL)

      USE PP3D_MPI, ONLY:myid
      IMPLICIT NONE
      REAL*4 VVOL(*)
      INTEGER KVERT(8,*),IntPhase(*),NEL
      REAL*8  DCORVG(3,*),DVAL(*),DRHO_ELEM(*),DMID(3,*),DENS1,DENS2
      REAL*8  DFRAC(*)
      INTEGER IEL,IP,IVT,IVERT,ITET
      INTEGER ITETRA1,ITETRA2,ITETRA3,ITETRA4
      REAL*8  DX(3,8),DV(8),DXTETRA(3,4),DVTETRA(4)
      REAL*8  DVOL1,DVOL2,DV1,DV2,DFRAC1,DFRAC2,DAUX

      DFRAC1 = 0d0
      DFRAC2 = 0d0

      DO 100 IEL=1,NEL

!--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
      IF(IntPhase(IEL).NE.0) THEN
       IF(IntPhase(IEL).GT.0) THEN
        DFRAC1 = DFRAC1 + DBLE(VVOL(IEL))
        DRHO_ELEM(IEL) = DENS1
        DFRAC(IEL) = 0d0
       END IF
       IF(IntPhase(IEL).LT.0) THEN
        DFRAC2 = DFRAC2 + DBLE(VVOL(IEL))
        DRHO_ELEM(IEL) = DENS2
        DFRAC(IEL) = 1d0
       END IF
       GOTO 100
      END IF
!-------------------------------------------------------------------

      IP = 4*(IEL-1)
      DO IVT=1,8
       IVERT = KVERT(IVT,IEL)
       DX(1,IVT) = DMID(1,IEL)-DCORVG(1,IVERT)
       DX(2,IVT) = DMID(2,IEL)-DCORVG(2,IVERT)
       DX(3,IVT) = DMID(3,IEL)-DCORVG(3,IVERT)
       DAUX = DSQRT(DX(1,IVT)**2d0+DX(2,IVT)**2d0+DX(3,IVT)**2d0)
       DV(IVT)   = (DVAL(IP+1)           + DVAL(IP+2)*DX(1,IVT) + &
                    DVAL(IP+3)*DX(2,IVT) + DVAL(IP+4)*DX(3,IVT))/DAUX
      END DO

      DVOL1 = 0d0
      DVOL2 = 0d0

      DO ITET=1,6
       IF (ITET.EQ.1) THEN
        ITETRA1 = 1
        ITETRA2 = 2
        ITETRA3 = 4
        ITETRA4 = 6
        GOTO 110
       END IF

       IF (ITET.EQ.2) THEN
        ITETRA1 = 2
        ITETRA2 = 3
        ITETRA3 = 4
        ITETRA4 = 7
        GOTO 110
       END IF

       IF (ITET.EQ.3) THEN
        ITETRA1 = 2
        ITETRA2 = 4
        ITETRA3 = 6
        ITETRA4 = 7
        GOTO 110
       END IF

       IF (ITET.EQ.4) THEN
        ITETRA1 = 1
        ITETRA2 = 6
        ITETRA3 = 4
        ITETRA4 = 5
        GOTO 110
       END IF

       IF (ITET.EQ.5) THEN
        ITETRA1 = 4
        ITETRA2 = 6
        ITETRA3 = 7
        ITETRA4 = 8
        GOTO 110
       END IF

       IF (ITET.EQ.6) THEN
        ITETRA1 = 4
        ITETRA2 = 5
        ITETRA3 = 6
        ITETRA4 = 8
        GOTO 110
       END IF

110    CONTINUE

       DXTETRA(:,1) = DX(:,ITETRA1)
       DXTETRA(:,2) = DX(:,ITETRA2)
       DXTETRA(:,3) = DX(:,ITETRA3)
       DXTETRA(:,4) = DX(:,ITETRA4)

       DVTETRA(1) = DV(ITETRA1)
       DVTETRA(2) = DV(ITETRA2)
       DVTETRA(3) = DV(ITETRA3)
       DVTETRA(4) = DV(ITETRA4)

       CALL SubTetra(DV1,DV2,DXTETRA,DVTETRA)
       DVOL1 = DVOL1 + DV1
       DVOL2 = DVOL2 + DV2
       DFRAC1 = DFRAC1 + DV1
       DFRAC2 = DFRAC2 + DV2

      END DO

      DRHO_ELEM(IEL) = (DVOL1*DENS2 + DVOL2*DENS1)/(DVOL1+DVOL2)
      DFRAC(IEL) = DVOL1/(DVOL1+DVOL2)

100   CONTINUE

      END
