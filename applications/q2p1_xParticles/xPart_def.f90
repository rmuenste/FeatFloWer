MODULE xPart_def

USE var_QuadScalar, only : myParticleParam

REAL*8 dist_CGAL,cdx,cdy,cdz
REAL*8 distance, dNormal(3),shellE(27)

real*8 DJ(8,3),q8
PARAMETER (Q8=0.125D0)
REAL*8 DJACI(3,3),DJAC(3,3),DBAS(27),DerBAS(27,3),DV0(3,27),DV1(3,27),DETJ
REAL*8 XI1,XI2,XI3,XX,YY,ZZ
REAL*8 :: dPI = 4d0*dATAN(1d0)
REAL*8  DV_Loc(3)
REAL*8 :: dParticleVelo(3) = 0d0
REAL*8  P8(3,8), P(3),PR(3),P_new(3)
LOGICAL OK_FLAG
REAL*8 diff1,diff2,diff3
INTEGER iter

REAL*8 :: DeltaT

CHARACTER cMeshFile*128
CHARACTER :: cParamFile*(128) = '_data/ParticleTracerInput.dat'
INTEGER :: nIter = 1000
real*8  :: xFactor = 1.0
integer :: xChunks = 5

REAL*8 :: d_CorrDist = 0.25d0,dTimeStep=1.00d1, minDist = -0.01d0
INTEGER :: nTime = 10

 CONTAINS
 
!========================================================================================
!                            Sub: GetPointFromElement
!========================================================================================
SUBROUTINE GetPointFromElement(P,iStart,BF,PR,iP)
LOGICAL BF

REAL*8  P(3),PR(3)

DJ(1,1)=( P8(1,1)+P8(1,2)+P8(1,3)+P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(1,2)=( P8(2,1)+P8(2,2)+P8(2,3)+P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(1,3)=( P8(3,1)+P8(3,2)+P8(3,3)+P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(2,1)=(-P8(1,1)+P8(1,2)+P8(1,3)-P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(2,2)=(-P8(2,1)+P8(2,2)+P8(2,3)-P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(2,3)=(-P8(3,1)+P8(3,2)+P8(3,3)-P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(3,1)=(-P8(1,1)-P8(1,2)+P8(1,3)+P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(3,2)=(-P8(2,1)-P8(2,2)+P8(2,3)+P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(3,3)=(-P8(3,1)-P8(3,2)+P8(3,3)+P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(4,1)=(-P8(1,1)-P8(1,2)-P8(1,3)-P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(4,2)=(-P8(2,1)-P8(2,2)-P8(2,3)-P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(4,3)=(-P8(3,1)-P8(3,2)-P8(3,3)-P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(5,1)=( P8(1,1)-P8(1,2)+P8(1,3)-P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(5,2)=( P8(2,1)-P8(2,2)+P8(2,3)-P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(5,3)=( P8(3,1)-P8(3,2)+P8(3,3)-P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(6,1)=( P8(1,1)-P8(1,2)-P8(1,3)+P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(6,2)=( P8(2,1)-P8(2,2)-P8(2,3)+P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(6,3)=( P8(3,1)-P8(3,2)-P8(3,3)+P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(7,1)=( P8(1,1)+P8(1,2)-P8(1,3)-P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(7,2)=( P8(2,1)+P8(2,2)-P8(2,3)-P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(7,3)=( P8(3,1)+P8(3,2)-P8(3,3)-P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(8,1)=(-P8(1,1)+P8(1,2)-P8(1,3)+P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(8,2)=(-P8(2,1)+P8(2,2)-P8(2,3)+P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(8,3)=(-P8(3,1)+P8(3,2)-P8(3,3)+P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8

IF (istart.eq.1) THEN
 XI1 = -1d0; XI2 = -1d0; XI3 = -1d0
END IF
IF (istart.eq.2) THEN
 XI1 = +1d0; XI2 = -1d0; XI3 = -1d0
END IF
IF (istart.eq.3) THEN
 XI1 = +1d0; XI2 = +1d0; XI3 = -1d0
END IF
IF (istart.eq.4) THEN
 XI1 = -1d0; XI2 = +1d0; XI3 = -1d0
END IF
IF (istart.eq.5) THEN
 XI1 = -1d0; XI2 = -1d0; XI3 = 1d0
END IF
IF (istart.eq.6) THEN
 XI1 = +1d0; XI2 = -1d0; XI3 = 1d0
END IF
IF (istart.eq.7) THEN
 XI1 = +1d0; XI2 = +1d0; XI3 = 1d0
END IF
IF (istart.eq.8) THEN
 XI1 = -1d0; XI2 = +1d0; XI3 = 1d0
END IF

DO iter = 1,80

  DJAC(1,1)=DJ(2,1)+DJ(5,1)*XI2+DJ(6,1)*XI3+DJ(8,1)*XI2*XI3
  DJAC(1,2)=DJ(3,1)+DJ(5,1)*XI1+DJ(7,1)*XI3+DJ(8,1)*XI1*XI3
  DJAC(1,3)=DJ(4,1)+DJ(6,1)*XI1+DJ(7,1)*XI2+DJ(8,1)*XI1*XI2
  DJAC(2,1)=DJ(2,2)+DJ(5,2)*XI2+DJ(6,2)*XI3+DJ(8,2)*XI2*XI3
  DJAC(2,2)=DJ(3,2)+DJ(5,2)*XI1+DJ(7,2)*XI3+DJ(8,2)*XI1*XI3
  DJAC(2,3)=DJ(4,2)+DJ(6,2)*XI1+DJ(7,2)*XI2+DJ(8,2)*XI1*XI2
  DJAC(3,1)=DJ(2,3)+DJ(5,3)*XI2+DJ(6,3)*XI3+DJ(8,3)*XI2*XI3
  DJAC(3,2)=DJ(3,3)+DJ(5,3)*XI1+DJ(7,3)*XI3+DJ(8,3)*XI1*XI3
  DJAC(3,3)=DJ(4,3)+DJ(6,3)*XI1+DJ(7,3)*XI2+DJ(8,3)*XI1*XI2
!   DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
!        -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
!        +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))

  XX = DJ(1,1) + DJAC(1,1)*XI1 + DJ(3,1)*XI2 + DJ(4,1)*XI3 + DJ(7,1)*XI2*XI3
  YY = DJ(1,2) + DJ(2,2)*XI1 + DJAC(2,2)*XI2 + DJ(4,2)*XI3 + DJ(6,2)*XI1*XI3
  ZZ = DJ(1,3) + DJ(2,3)*XI1 + DJ(3,3)*XI2 + DJAC(3,3)*XI3 + DJ(5,3)*XI1*XI2

  CALL M33INV (DJAC, DJACI, OK_FLAG)

  dFact = -0.33d0
  diff1 = -dFact*(DJACI(1,1)*(P(1)-XX) + DJACI(1,2)*(P(2)-YY) + DJACI(1,3)*(P(3)-ZZ))
  diff2 = -dFact*(DJACI(2,1)*(P(1)-XX) + DJACI(2,2)*(P(2)-YY) + DJACI(2,3)*(P(3)-ZZ))
  diff3 = -dFact*(DJACI(3,1)*(P(1)-XX) + DJACI(3,2)*(P(2)-YY) + DJACI(3,3)*(P(3)-ZZ))
  XI1 = XI1 + diff1
  XI2 = XI2 + diff2
  XI3 = XI3 + diff3
  daux = diff1*diff1 + diff2*diff2 + diff3*diff3

!   WRITE(*,'(I3,9ES12.2)') iter,XI1,XI2,XI3, XX,YY,ZZ, diff1,diff2,diff3
  IF (iter.gt.4.and.(ABS(XI1).GT.1.5d0.OR.ABS(XI2).GT.1.5d0.OR.ABS(XI3).GT.1.5d0)) GOTO 2
  IF (iter.gt.4.and.daux.lt.0.00001d0) GOTO 1

END DO

1 CONTINUE
IF (ABS(XI1).GT.1.01d0.OR.ABS(XI2).GT.1.01d0.OR.ABS(XI3).GT.1.01d0) GOTO 2

BF = .TRUE.
PR = [XI1,XI2,XI3]
RETURN

2 CONTINUE

END SUBROUTINE GetPointFromElement
!========================================================================================
!                           Sub: Search_Particle
!========================================================================================
SUBROUTINE MovePointThroughElement(tLevel,tEnd,tStart,iP,iE)
USE var_QuadScalar, ONLY : myParticleParam
USE PP3D_MPI, ONLY : myid
implicit none
LOGICAL :: BOUT=.FALSE.
REAL*8  tLevel,tEnd,tStart
REAL*8 :: dAlpha,XB,YB,ZB
REAL*8 :: RK_Velo(3,4),dVelo,dElemSize,dFracTime,dist,dfrac
integer imove,ip,ie
logical :: ok_flag
REAL*8 :: rho_p,rho_l,d_r,d_A,d_V,d_g(3),d_Up(3),d_U(3),dauxU(3),C_D,magU,d_Mu,dRE
REAL*8 :: tau_p, f_p(3)

iMove = 0
dVelo = 0d0

XI1 = PR(1)
XI2 = PR(2)
XI3 = PR(3)

! WRITE(*,'(A,7E14.4)') 'my original position: ', P,tLevel
IF (bOut) WRITE(*,'(A,7E14.4)') 'my original position: ', XI1,XI2,XI3,P,tLevel

CALL RETURN_Velo()

DJ(1,1)=( P8(1,1)+P8(1,2)+P8(1,3)+P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(1,2)=( P8(2,1)+P8(2,2)+P8(2,3)+P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(1,3)=( P8(3,1)+P8(3,2)+P8(3,3)+P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(2,1)=(-P8(1,1)+P8(1,2)+P8(1,3)-P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(2,2)=(-P8(2,1)+P8(2,2)+P8(2,3)-P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(2,3)=(-P8(3,1)+P8(3,2)+P8(3,3)-P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(3,1)=(-P8(1,1)-P8(1,2)+P8(1,3)+P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(3,2)=(-P8(2,1)-P8(2,2)+P8(2,3)+P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(3,3)=(-P8(3,1)-P8(3,2)+P8(3,3)+P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(4,1)=(-P8(1,1)-P8(1,2)-P8(1,3)-P8(1,4)+P8(1,5)+P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(4,2)=(-P8(2,1)-P8(2,2)-P8(2,3)-P8(2,4)+P8(2,5)+P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(4,3)=(-P8(3,1)-P8(3,2)-P8(3,3)-P8(3,4)+P8(3,5)+P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(5,1)=( P8(1,1)-P8(1,2)+P8(1,3)-P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(5,2)=( P8(2,1)-P8(2,2)+P8(2,3)-P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(5,3)=( P8(3,1)-P8(3,2)+P8(3,3)-P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(6,1)=( P8(1,1)-P8(1,2)-P8(1,3)+P8(1,4)-P8(1,5)+P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(6,2)=( P8(2,1)-P8(2,2)-P8(2,3)+P8(2,4)-P8(2,5)+P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(6,3)=( P8(3,1)-P8(3,2)-P8(3,3)+P8(3,4)-P8(3,5)+P8(3,6)+P8(3,7)-P8(3,8))*Q8
DJ(7,1)=( P8(1,1)+P8(1,2)-P8(1,3)-P8(1,4)-P8(1,5)-P8(1,6)+P8(1,7)+P8(1,8))*Q8
DJ(7,2)=( P8(2,1)+P8(2,2)-P8(2,3)-P8(2,4)-P8(2,5)-P8(2,6)+P8(2,7)+P8(2,8))*Q8
DJ(7,3)=( P8(3,1)+P8(3,2)-P8(3,3)-P8(3,4)-P8(3,5)-P8(3,6)+P8(3,7)+P8(3,8))*Q8
DJ(8,1)=(-P8(1,1)+P8(1,2)-P8(1,3)+P8(1,4)+P8(1,5)-P8(1,6)+P8(1,7)-P8(1,8))*Q8
DJ(8,2)=(-P8(2,1)+P8(2,2)-P8(2,3)+P8(2,4)+P8(2,5)-P8(2,6)+P8(2,7)-P8(2,8))*Q8
DJ(8,3)=(-P8(3,1)+P8(3,2)-P8(3,3)+P8(3,4)+P8(3,5)-P8(3,6)+P8(3,7)-P8(3,8))*Q8

DJAC(1,1)=DJ(2,1)
DJAC(1,2)=DJ(3,1)
DJAC(1,3)=DJ(4,1)
DJAC(2,1)=DJ(2,2)
DJAC(2,2)=DJ(3,2)
DJAC(2,3)=DJ(4,2)
DJAC(3,1)=DJ(2,3)
DJAC(3,2)=DJ(3,3)
DJAC(3,3)=DJ(4,3)
DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))

dElemSize = (ABS(DETJ)**0.3333d0)

dParticleVelo = DV_Loc
 
dVelo = SQRT(dParticleVelo(1)*dParticleVelo(1) + dParticleVelo(2)*dParticleVelo(2) + dParticleVelo(3)*dParticleVelo(3))
DeltaT = myParticleParam%hSize*dElemSize/dVelo

IF ((tLevel + DeltaT).GT.tEnd) THEN
 DeltaT = tEnd - tLevel
END IF
dFracTime = (tLevel + 0.5d0*DeltaT - tStart)/(tEnd - tStart)

55 CONTINUE

 CALL FindPoint(1d0,dParticleVelo)  ! --> XI1,XI2,XI3

 CALL RETURN_Velo()                     ! --> DV_loc
 
 dParticleVelo = DV_Loc

 
 IF (isnan(magU)) pause
 
!   write(*,'(a,80es12.4)') 'aasaa : ',P_new-P
  
  P(1) = P_New(1)
  P(2) = P_New(2)
  P(3) = P_New(3)

  tLevel = tLevel + DeltaT
  iMove  = iMove + 1

  IF ((ABS(XI1).GT.1.01d0.OR.ABS(XI2).GT.1.01d0.OR.ABS(XI3).GT.1.01d0).OR.(tLevel.ge.tEnd)) THEN
  ! Reached the end of the element -> leave the routine
  GOTO 5
 ELSE

! We are still in the element - go back and continue marching in the element
  GOTO 55

 END IF

5 CONTINUE


end subroutine MovePointThroughElement
!========================================================================================
!                           Sub: RETURN_Velo
!========================================================================================
SUBROUTINE RETURN_Velo()
implicit none
!!!!!!!!!!!!!!!!!!!!!!
integer i

! Loop not neccesary because nothing depends on i?
! DO i=1,27
 DBAS( 1)=-Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
 DBAS( 2)= Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
 DBAS( 3)=-Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
 DBAS( 4)=Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
 DBAS( 5)=Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
 DBAS( 6)=-Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
 DBAS( 7)=Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
 DBAS( 8)=-Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
 DBAS( 9)=Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
 DBAS(10)=-Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
 DBAS(11)=-Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
 DBAS(12)=Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
 DBAS(13)=Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
 DBAS(14)=-Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
 DBAS(15)=Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
 DBAS(16)=-Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
 DBAS(17)=-Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
 DBAS(18)=Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
 DBAS(19)=Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
 DBAS(20)=-Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
 DBAS(21)= -Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
 DBAS(22)= -Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
 DBAS(23)= Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
 DBAS(24)= Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
 DBAS(25)= -Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
 DBAS(26)= Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
 DBAS(27)= Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
! END DO

DV_loc = 0d0
DO i=1,27
 DV_loc(1) = DV_loc(1) + DBAS(I)*DV0(1,I)
 DV_loc(2) = DV_loc(2) + DBAS(I)*DV0(2,I)
 DV_loc(3) = DV_loc(3) + DBAS(I)*DV0(3,I)
END DO

END SUBROUTINE RETURN_Velo
!========================================================================================
!                           Sub: RETURN_Distance
!========================================================================================
SUBROUTINE RETURN_Distance()
implicit none
!!!!!!!!!!!!!!!!!!!!!!
integer i
REAL*8 DerBas1,DerBas2,DerBas3,DETI

DBAS( 1)=-Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
DBAS( 2)= Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
DBAS( 3)=-Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
DBAS( 4)=Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
DBAS( 5)=Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
DBAS( 6)=-Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
DBAS( 7)=Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
DBAS( 8)=-Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
DBAS( 9)=Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
DBAS(10)=-Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
DBAS(11)=-Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
DBAS(12)=Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
DBAS(13)=Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
DBAS(14)=-Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
DBAS(15)=Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
DBAS(16)=-Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
DBAS(17)=-Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
DBAS(18)=Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
DBAS(19)=Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
DBAS(20)=-Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
DBAS(21)= -Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
DBAS(22)= -Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
DBAS(23)= Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
DBAS(24)= Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
DBAS(25)= -Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
DBAS(26)= Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
DBAS(27)= Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)

Distance = 0d0
DO i=1,27
 Distance  = Distance  + DBAS(I)*shellE(I)
END DO

END SUBROUTINE RETURN_Distance
!========================================================================================
!                           Sub: RETURN_Normal
!========================================================================================
SUBROUTINE RETURN_Normal()
implicit none
!!!!!!!!!!!!!!!!!!!!!!
integer i
REAL*8 DerBas1,DerBas2,DerBas3,DETI

DerBas(1,1)= Q8*(-1D0+2D0*XI1)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
DerBas(2,1)= Q8*(1D0+2D0*XI1)*XI2*(1D0-XI2)*XI3*(1D0-XI3)
DerBas(3,1)= -Q8*(1D0+2D0*XI1)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
DerBas(4,1)= -Q8*(-1D0+2D0*XI1)*XI2*(1D0+XI2)*XI3*(1D0-XI3)
DerBas(5,1)= -Q8*(-1D0+2D0*XI1)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
DerBas(6,1)= -Q8*(1D0+2D0*XI1)*XI2*(1D0-XI2)*XI3*(1D0+XI3)
DerBas(7,1)= Q8*(1D0+2D0*XI1)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
DerBas(8,1)= Q8*(-1D0+2D0*XI1)*XI2*(1D0+XI2)*XI3*(1D0+XI3)
DerBas(9,1)= -4D0*Q8*XI1*XI2*(1D0-XI2)*XI3*(1D0-XI3)
DerBas(10,1)= -Q8*(1D0+2D0*XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
DerBas(11,1)= 4D0*Q8*XI1*XI2*(1D0+XI2)*XI3*(1D0-XI3)
DerBas(12,1)=-Q8*(-1D0+2D0*XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
DerBas(13,1)=-Q8*(-1D0+2D0*XI1)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
DerBas(14,1)=-Q8*(1D0+2D0*XI1)*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
DerBas(15,1)=Q8*(1D0+2D0*XI1)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
DerBas(16,1)=Q8*(-1D0+2D0*XI1)*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
DerBas(17,1)=4D0*Q8*XI1*XI2*(1D0-XI2)*XI3*(1D0+XI3)
DerBas(18,1)=Q8*(1D0+2D0*XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
DerBas(19,1)=-4D0*Q8*XI1*XI2*(1D0+XI2)*XI3*(1D0+XI3)
DerBas(20,1)=Q8*(-1D0+2D0*XI1)*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
DerBas(21,1)=4D0*Q8*XI1*(2D0-2D0*XI2**2D0)*XI3*(1D0-XI3)
DerBas(22,1)=4D0*Q8*XI1*XI2*(1D0-XI2)*(2D0-2D0*XI3**2D0)
DerBas(23,1)=Q8*(1D0+2D0*XI1)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
DerBas(24,1)=-4D0*Q8*XI1*XI2*(1D0+XI2)*(2D0-2D0*XI3**2D0)
DerBas(25,1)=Q8*(-1D0+2D0*XI1)*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)
DerBas(26,1)=-4D0*Q8*XI1*(2D0-2D0*XI2**2D0)*XI3*(1D0+XI3)
DerBas(27,1)=-4D0*Q8*XI1*(2D0-2D0*XI2**2D0)*(2D0-2D0*XI3**2D0)

DerBas(1,2)=Q8*XI1*(1D0-XI1)*(-1D0+2D0*XI2)*XI3*(1D0-XI3)
DerBas(2,2)=-Q8*XI1*(1D0+XI1)*(-1D0+2D0*XI2)*XI3*(1D0-XI3)
DerBas(3,2)=-Q8*XI1*(1D0+XI1)*(1D0+2D0*XI2)*XI3*(1D0-XI3)
DerBas(4,2)=Q8*XI1*(1D0-XI1)*(1D0+2D0*XI2)*XI3*(1D0-XI3)
DerBas(5,2)=-Q8*XI1*(1D0-XI1)*(-1D0+2D0*XI2)*XI3*(1D0+XI3)
DerBas(6,2)=Q8*XI1*(1D0+XI1)*(-1D0+2D0*XI2)*XI3*(1D0+XI3)
DerBas(7,2)=Q8*XI1*(1D0+XI1)*(1D0+2D0*XI2)*XI3*(1D0+XI3)
DerBas(8,2)=-Q8*XI1*(1D0-XI1)*(1D0+2D0*XI2)*XI3*(1D0+XI3)
DerBas(9,2)= -Q8*(2D0-2D0*XI1**2D0)*(-1D0+2D0*XI2)*XI3*(1D0-XI3)
DerBas(10,2)=4D0*Q8*XI1*(1D0+XI1)*XI2*XI3*(1D0-XI3)
DerBas(11,2)=-Q8*(2D0-2D0*XI1**2D0)*(1D0+2D0*XI2)*XI3*(1D0-XI3)
DerBas(12,2)=-4D0*Q8*XI1*(1D0-XI1)*XI2*XI3*(1D0-XI3)
DerBas(13,2)=-Q8*XI1*(1D0-XI1)*(-1D0+2D0*XI2)*(2D0-2D0*XI3**2D0)
DerBas(14,2)=Q8*XI1*(1D0+XI1)*(-1D0+2D0*XI2)*(2D0-2D0*XI3**2D0)
DerBas(15,2)=Q8*XI1*(1D0+XI1)*(1D0+2D0*XI2)*(2D0-2D0*XI3**2D0)
DerBas(16,2)=-Q8*XI1*(1D0-XI1)*(1D0+2D0*XI2)*(2D0-2D0*XI3**2D0)
DerBas(17,2)=Q8*(2D0-2D0*XI1**2D0)*(-1D0+2D0*XI2)*XI3*(1D0+XI3)
DerBas(18,2)=-4D0*Q8*XI1*(1D0+XI1)*XI2*XI3*(1D0+XI3)
DerBas(19,2)=Q8*(2D0-2D0*XI1**2D0)*(1D0+2D0*XI2)*XI3*(1D0+XI3)
DerBas(20,2)=4D0*Q8*XI1*(1D0-XI1)*XI2*XI3*(1D0+XI3)
DerBas(21,2)=4D0*Q8*(2D0-2D0*XI1**2D0)*XI2*XI3*(1D0-XI3)
DerBas(22,2)=Q8*(2D0-2D0*XI1**2D0)*(-1D0+2D0*XI2)*(2D0-2D0*XI3**2D0)
DerBas(23,2)=-4D0*Q8*XI1*(1D0+XI1)*XI2*(2D0-2D0*XI3**2D0)
DerBas(24,2)=Q8*(2D0-2D0*XI1**2D0)*(1D0+2D0*XI2)*(2D0-2D0*XI3**2D0)
DerBas(25,2)=4D0*Q8*XI1*(1D0-XI1)*XI2*(2D0-2D0*XI3**2D0)
DerBas(26,2)=-4D0*Q8*(2D0-2D0*XI1**2D0)*XI2*XI3*(1D0+XI3)
DerBas(27,2)=-4D0*Q8*(2D0-2D0*XI1**2D0)*XI2*(2D0-2D0*XI3**2D0)

DerBas(1,3)=Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*(-1D0+2D0*XI3)
DerBas(2,3)=-Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*(-1D0+2D0*XI3)
DerBas(3,3)=Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*(-1D0+2D0*XI3)
DerBas(4,3)=-Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*(-1D0+2D0*XI3)
DerBas(5,3)=Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*(1D0+2D0*XI3)
DerBas(6,3)=-Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*(1D0+2D0*XI3)
DerBas(7,3)=Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*(1D0+2D0*XI3)
DerBas(8,3)=-Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*(1D0+2D0*XI3)
DerBas(9,3)=-Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*(-1D0+2D0*XI3)
DerBas(10,3)=Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*(-1D0+2D0*XI3)
DerBas(11,3)=Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*(-1D0+2D0*XI3)
DerBas(12,3)=Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*(-1D0+2D0*XI3)
DerBas(13,3)=-4D0*Q8*XI1*(1D0-XI1)*XI2*(1D0-XI2)*XI3
DerBas(14,3)=4D0*Q8*XI1*(1D0+XI1)*XI2*(1D0-XI2)*XI3
DerBas(15,3)=-4D0*Q8*XI1*(1D0+XI1)*XI2*(1D0+XI2)*XI3
DerBas(16,3)=4D0*Q8*XI1*(1D0-XI1)*XI2*(1D0+XI2)*XI3
DerBas(17,3)=-Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*(1D0+2D0*XI3)
DerBas(18,3)=Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*(1D0+2D0*XI3)
DerBas(19,3)=Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*(1D0+2D0*XI3)
DerBas(20,3)=-Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*(1D0+2D0*XI3)
DerBas(21,3)=Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*(-1D0+2D0*XI3)
DerBas(22,3)=4D0*Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0-XI2)*XI3
DerBas(23,3)=-4D0*Q8*XI1*(1D0+XI1)*(2D0-2D0*XI2**2D0)*XI3
DerBas(24,3)=-4D0*Q8*(2D0-2D0*XI1**2D0)*XI2*(1D0+XI2)*XI3
DerBas(25,3)=4D0*Q8*XI1*(1D0-XI1)*(2D0-2D0*XI2**2D0)*XI3
DerBas(26,3)=Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*(1D0+2D0*XI3)
DerBas(27,3)=-4D0*Q8*(2D0-2D0*XI1**2D0)*(2D0-2D0*XI2**2D0)*XI3

DETI=1D0/DETJ

DO i=1,27
 DerBas1=DerBas(i,1)
 DerBas2=DerBas(i,2)
 DerBas3=DerBas(i,3)

 DerBas(I,1)= DETI*( &
 +DerBas1*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
 -DerBas2*(DJAC(2,1)*DJAC(3,3)-DJAC(3,1)*DJAC(2,3)) &
 +DerBas3*(DJAC(2,1)*DJAC(3,2)-DJAC(3,1)*DJAC(2,2)))

 DerBas(I,2)= DETI*( &
 -DerBas1*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
 +DerBas2*(DJAC(1,1)*DJAC(3,3)-DJAC(3,1)*DJAC(1,3)) &
 -DerBas3*(DJAC(1,1)*DJAC(3,2)-DJAC(3,1)*DJAC(1,2)))

 DerBas(I,3)= DETI*( &
 +DerBas1*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3)) &
 -DerBas2*(DJAC(1,1)*DJAC(2,3)-DJAC(2,1)*DJAC(1,3)) &
 +DerBas3*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2)))
 
END DO

dNormal  = 0d0
DO i=1,27
 dNormal(1) = dNormal(1) + DerBas(I,1)*shellE(I)
 dNormal(2) = dNormal(2) + DerBas(I,2)*shellE(I)
 dNormal(3) = dNormal(3) + DerBas(I,3)*shellE(I)
END DO

END SUBROUTINE RETURN_Normal
!========================================================================================
!                           Sub: FindPoint
!========================================================================================
SUBROUTINE FindPoint(fFact,dLocalVelo)
implicit none
REAL*8 fFact,dLocalVelo(3)
!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 :: DV_Rot(3),dVV(3)
REAL*8 :: dAlpha,XB,YB,ZB,dtheta
REAL*8 :: dRR,dU_r,dU_theta,dU_z,dRho,dRho0,dRho1,dZ,daux,dFact
!LOGICAL :: bRotationalMovement=.false.

P_New = P + fFact*DeltaT * dParticleVelo ! dParticleVelo

DO iter = 1,80

  DJAC(1,1)=DJ(2,1)+DJ(5,1)*XI2+DJ(6,1)*XI3+DJ(8,1)*XI2*XI3
  DJAC(1,2)=DJ(3,1)+DJ(5,1)*XI1+DJ(7,1)*XI3+DJ(8,1)*XI1*XI3
  DJAC(1,3)=DJ(4,1)+DJ(6,1)*XI1+DJ(7,1)*XI2+DJ(8,1)*XI1*XI2
  DJAC(2,1)=DJ(2,2)+DJ(5,2)*XI2+DJ(6,2)*XI3+DJ(8,2)*XI2*XI3
  DJAC(2,2)=DJ(3,2)+DJ(5,2)*XI1+DJ(7,2)*XI3+DJ(8,2)*XI1*XI3
  DJAC(2,3)=DJ(4,2)+DJ(6,2)*XI1+DJ(7,2)*XI2+DJ(8,2)*XI1*XI2
  DJAC(3,1)=DJ(2,3)+DJ(5,3)*XI2+DJ(6,3)*XI3+DJ(8,3)*XI2*XI3
  DJAC(3,2)=DJ(3,3)+DJ(5,3)*XI1+DJ(7,3)*XI3+DJ(8,3)*XI1*XI3
  DJAC(3,3)=DJ(4,3)+DJ(6,3)*XI1+DJ(7,3)*XI2+DJ(8,3)*XI1*XI2

  XX = DJ(1,1) + DJAC(1,1)*XI1 + DJ(3,1)*XI2 + DJ(4,1)*XI3 + DJ(7,1)*XI2*XI3
  YY = DJ(1,2) + DJ(2,2)*XI1 + DJAC(2,2)*XI2 + DJ(4,2)*XI3 + DJ(6,2)*XI1*XI3
  ZZ = DJ(1,3) + DJ(2,3)*XI1 + DJ(3,3)*XI2 + DJAC(3,3)*XI3 + DJ(5,3)*XI1*XI2

  CALL M33INV (DJAC, DJACI, OK_FLAG)

  dFact = -0.33d0
  diff1 = -dFact*(DJACI(1,1)*(P_New(1)-XX) + DJACI(1,2)*(P_New(2)-YY) + DJACI(1,3)*(P_New(3)-ZZ))
  diff2 = -dFact*(DJACI(2,1)*(P_New(1)-XX) + DJACI(2,2)*(P_New(2)-YY) + DJACI(2,3)*(P_New(3)-ZZ))
  diff3 = -dFact*(DJACI(3,1)*(P_New(1)-XX) + DJACI(3,2)*(P_New(2)-YY) + DJACI(3,3)*(P_New(3)-ZZ))
  XI1 = XI1 + diff1
  XI2 = XI2 + diff2
  XI3 = XI3 + diff3
  daux = diff1*diff1 + diff2*diff2 + diff3*diff3

  IF (iter.gt.3.and.daux.lt.0.00001d0) GOTO 1

END DO

1 CONTINUE

end subroutine FindPoint
!========================================================================================
!                                    Sub: M33Inv
!========================================================================================
subroutine M33INV (A, AINV, OK_FLAG)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
LOGICAL, INTENT(OUT) :: OK_FLAG

DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
DOUBLE PRECISION :: DET
DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


DET =   A(1,1)*A(2,2)*A(3,3)  &
      - A(1,1)*A(2,3)*A(3,2)  &
      - A(1,2)*A(2,1)*A(3,3)  &
      + A(1,2)*A(2,3)*A(3,1)  &
      + A(1,3)*A(2,1)*A(3,2)  &
      - A(1,3)*A(2,2)*A(3,1)

 IF (ABS(DET) .LE. EPS) THEN
    AINV = 0.0D0
    OK_FLAG = .FALSE.
    RETURN
 END IF

 COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
 COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
 COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
 COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
 COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
 COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
 COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
 COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
 COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

 AINV = TRANSPOSE(COFACTOR) / DET

 OK_FLAG = .TRUE.

 RETURN

end subroutine M33INV

END MODULE xPart_def