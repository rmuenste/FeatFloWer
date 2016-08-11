      SUBROUTINE SubTetra(DV1,DV2,DTETRA,DIST)
      IMPLICIT NONE
      REAL*8 DV,DV1,DV2,DTETRA(3,4),DIST(4),DSTETRA(3,4)
      REAL*8 PAR,DSV,DSPV,DSPRISM(3,6)
      INTEGER NJALFA,NIALFA,I,K,IP,JP(3),KP1(2),KP2(2),k1,k2

      NJALFA=0
      NIALFA=0
      DO I=1,4
        IF (DIST(I).GE.0d0) THEN
         NJALFA=NJALFA+1
        ENDIF
        IF (DIST(I).LT.0d0) THEN
         NIALFA=NIALFA+1
        ENDIF
      ENDDO

      CALL GetTetraVolume(DTETRA,DV)

      IF(NJALFA.EQ.4.OR.NIALFA.EQ.4) THEN
       IF(NIALFA.EQ.4) THEN
        DV1 = DV
        DV2 = 0d0
       END IF
       IF(NJALFA.EQ.4) THEN
        DV1 = 0d0
        DV2 = DV
       END IF
       GOTO 100
      ELSE
       IF(NJALFA.EQ.1.OR.NIALFA.EQ.1) THEN
        DO I=1,4
         IF (DIST(I).GE.0d0.AND.NJALFA.EQ.1) IP=I
         IF (DIST(I).LT.0d0.AND.NIALFA.EQ.1) IP=I
        ENDDO
        k=1
        DO I=1,4
         IF (IP.NE.I) THEN
          JP(k)=I
          k=k+1
         END IF
        ENDDO
        DSTETRA(:,1) = DTETRA(:,IP)
        DO I=1,3
         PAR = -DIST(IP)/(DIST(JP(I))-DIST(IP))
         DSTETRA(:,I+1)=DTETRA(:,IP)+PAR*(DTETRA(:,JP(I))-DTETRA(:,IP))
        END DO
        CALL GetTetraVolume(DSTETRA,DSV)
        IF(NIALFA.EQ.1) THEN
         DV1 = DSV
         DV2 = DV-DSV
        END IF
        IF(NJALFA.EQ.1) THEN
         DV1 = DV-DSV
         DV2 = DSV
        END IF
       ELSE

        k1=1
        k2=1
        DO I=1,4
         IF (DIST(I).GE.0d0) THEN
          KP1(k1)=I
          k1 = k1 + 1
         END IF
         IF (DIST(I).LT.0d0) THEN
          KP2(k2)=I
          k2 = k2 + 1
         END IF
        ENDDO

!         WRITE(*,*) KP1
!         WRITE(*,*) KP2

        DSPRISM(:,1) = DTETRA(:,KP1(1))
        DSPRISM(:,4) = DTETRA(:,KP1(2))

        DO I=1,2
         PAR = -DIST(KP1(1))/(DIST(KP2(I))-DIST(KP1(1)))
         DSPRISM(:,I+1)=DTETRA(:,KP1(1))+PAR*(DTETRA(:,KP2(I))-DTETRA(:,KP1(1)))
        END DO

        DO I=1,2
         PAR = -DIST(KP1(2))/(DIST(KP2(I))-DIST(KP1(2)))
         DSPRISM(:,I+4)=DTETRA(:,KP1(2))+PAR*(DTETRA(:,KP2(I))-DTETRA(:,KP1(2)))
        END DO

        DSTETRA(:,1) = DSPRISM(:,1)
        DSTETRA(:,2) = DSPRISM(:,2)
        DSTETRA(:,3) = DSPRISM(:,3)
        DSTETRA(:,4) = DSPRISM(:,5)
        CALL GetTetraVolume(DSTETRA,DSPV)
        DSV = 0d0 + DSPV

        DSTETRA(:,1) = DSPRISM(:,4)
        DSTETRA(:,2) = DSPRISM(:,5)
        DSTETRA(:,3) = DSPRISM(:,6)
        DSTETRA(:,4) = DSPRISM(:,3)
        CALL GetTetraVolume(DSTETRA,DSPV)
        DSV = DSV + DSPV

        DSTETRA(:,1) = DSPRISM(:,1)
        DSTETRA(:,2) = DSPRISM(:,5)
        DSTETRA(:,3) = DSPRISM(:,3)
        DSTETRA(:,4) = DSPRISM(:,4)
        CALL GetTetraVolume(DSTETRA,DSPV)
        DSV = DSV + DSPV

        DV1 = DV-DSV
        DV2 = DSV

!         do i=1,6
!         WRITE(*,*) DSPRISM(:,i)
!         end do

       END IF
      END IF

100   CONTINUE
!       WRITE(*,*) DV1,DV2

      END

      SUBROUTINE GetTetraVolume(DT,DV)
      IMPLICIT NONE
      REAL*8 DV,DT(3,4),DA(3),DN(3)
      REAL*8 D1(3),D2(3),D3(3)

       D1(1) = (DT(1,2)-DT(1,1))
       D1(2) = (DT(2,2)-DT(2,1))
       D1(3) = (DT(3,2)-DT(3,1))

       D2(1) = (DT(1,3)-DT(1,1))
       D2(2) = (DT(2,3)-DT(2,1))
       D2(3) = (DT(3,3)-DT(3,1))

       DA(1) = +(D1(2)*D2(3)-D1(3)*D2(2))
       DA(2) = -(D1(1)*D2(3)-D1(3)*D2(1))
       DA(3) = +(D1(1)*D2(2)-D1(2)*D2(1))

       D3(1) = (DT(1,4)-DT(1,1))
       D3(2) = (DT(2,4)-DT(2,1))
       D3(3) = (DT(3,4)-DT(3,1))

       DV    = ABS(DA(1)*D3(1)+DA(2)*D3(2)+DA(3)*D3(3))/6d0

      END

