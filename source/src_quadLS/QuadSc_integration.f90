!=========================================================================
! QuadSc_integration.f90
!
! Integration and analysis utility functions
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE Integrate_DIE_Flowrate(dcorvg,karea,kvert,nel,dVolFlow,iPar)
REAL*8 dcorvg(3,*),dVolFlow
INTEGER karea(6,*),kvert(8,*),nel,iPar
!---------------------------------
INTEGER NeighA(4,6)
REAL*8 P(3),dAN(3),dV
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4

dVolFlow = 0d0

if (iPar.eq.0) then
 k=1
 DO i=1,nel
  DO j=1,6
   IF (k.eq.karea(j,i)) THEN
    IF (myBoundary%iInflow(nvt+net+k).ne.0) THEN
     ivt1 = kvert(NeighA(1,j),i)
     ivt2 = kvert(NeighA(2,j),i)
     ivt3 = kvert(NeighA(3,j),i)
     ivt4 = kvert(NeighA(4,j),i)
     CALL GET_NormalArea(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dcorvg(1:3,nvt+net+nat+i),dAN)
     dV = dAN(1)*QuadSc%ValU(nvt+net+k) + dAN(2)*QuadSc%ValV(nvt+net+k) + dAN(3)*QuadSc%ValW(nvt+net+k)
!      write(*,'(4ES12.4)') dAN
     dVolFlow = dVolFlow + dV
    END IF
    k = k + 1
   END IF
  END DO
 END DO
end if

if (iPar.eq.1) then
 k=1
 DO i=1,nel
  DO j=1,6
   IF (k.eq.karea(j,i)) THEN
    IF (myBoundary%bOutflow(nvt+net+k)) THEN
     ivt1 = kvert(NeighA(1,j),i)
     ivt2 = kvert(NeighA(2,j),i)
     ivt3 = kvert(NeighA(3,j),i)
     ivt4 = kvert(NeighA(4,j),i)
     CALL GET_NormalArea(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dcorvg(1:3,nvt+net+nat+i),dAN)
     dAN = - dAN
     dV = dAN(1)*QuadSc%ValU(nvt+net+k) + dAN(2)*QuadSc%ValV(nvt+net+k) + dAN(3)*QuadSc%ValW(nvt+net+k)
     dVolFlow = dVolFlow + dV
    END IF
    k = k + 1
   END IF
  END DO
 END DO
end if

END SUBROUTINE Integrate_DIE_Flowrate
!=========================================================================
!
!=========================================================================
SUBROUTINE IntegrateFlowrate(dcorvg,karea,kvert,nel,dVolFlow,dPar)
REAL*8 dcorvg(3,*),dVolFlow
INTEGER karea(6,*),kvert(8,*),nel
REAL*8 dPar
!---------------------------------
INTEGER NeighA(4,6)
REAL*8 P(3),dA
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4

dVolFlow = 0d0

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   IF (abs(dcorvg(3,ivt1)-dPar).lt.1d-4.and.abs(dcorvg(3,ivt2)-dPar).lt.1d-4.and. &
       abs(dcorvg(3,ivt3)-dPar).lt.1d-4.and.abs(dcorvg(3,ivt4)-dPar).lt.1d-4) then
       CALL GET_area(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dA)
       dVolFlow = dVolFlow + dA*(QuadSc%ValW(nvt+net+k))
   END IF
   k = k + 1
  END IF
 END DO
END DO

END SUBROUTINE IntegrateFlowrate
!=========================================================================
!
!=========================================================================
SUBROUTINE IntegratePressure(dcorvg,karea,kvert,nel,dIntPres,dArea,dPar)
REAL*8 dcorvg(3,*),dIntPres
INTEGER karea(6,*),kvert(8,*),nel
REAL*8 dPar,dArea
!---------------------------------
INTEGER NeighA(4,6)
REAL*8 P(3),dA
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4

dIntPres = 0d0
dArea   = 0d0

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   IF (abs(dcorvg(3,ivt1)-dPar).lt.1d-4.and.abs(dcorvg(3,ivt2)-dPar).lt.1d-4.and. &
       abs(dcorvg(3,ivt3)-dPar).lt.1d-4.and.abs(dcorvg(3,ivt4)-dPar).lt.1d-4) then
       CALL GET_area(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dA)
       dIntPres = dIntPres + dA*(LinSc%ValP(NLMAX)%x(4*(i-1)+1))
       dArea = dArea + dA
   END IF
   k = k + 1
  END IF
 END DO
END DO

END SUBROUTINE IntegratePressure
!=========================================================================
!
!=========================================================================
SUBROUTINE IntegratePressureAtInflow(dcorvg,karea,kvert,nel,dIntPres,dArea)
REAL*8 dcorvg(3,*),dIntPres
INTEGER karea(6,*),kvert(8,*),nel
REAL*8 dArea
!---------------------------------
INTEGER NeighA(4,6)
REAL*8 P(3),dA
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4

dIntPres = 0d0
dArea   = 0d0

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   IF (myBoundary%iInflow(ivt1).lt.0.and.myBoundary%iInflow(ivt2).lt.0.and.myBoundary%iInflow(ivt3).lt.0.and.myBoundary%iInflow(ivt4).lt.0) then
       CALL GET_area(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dA)
       dIntPres = dIntPres + dA*(LinSc%ValP(NLMAX)%x(4*(i-1)+1))
       dArea = dArea + dA
   END IF
   k = k + 1
  END IF
 END DO
END DO

END SUBROUTINE IntegratePressureAtInflow
!=========================================================================
!
!=========================================================================
subroutine IntegrateQuantities(mfile)
  INTEGER mfile
  REAL*8 dArray(8)
  REAL*8 :: dR=0.25d0, myPI=4d0*ATAN(1d0), dWidth=0.05d0
  EXTERNAL E013

  write(*,*)'New untested subroutine'
  stop

  IF (myid.ne.0) THEN
    dArray = 0d0
    ILEV = NLMAX
    CALL SETLEV(2)
    CALL GetSurface(KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),&
      myQ2coor,E013,dArray(1))
    CALL GetVolume(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),&
      myQ2coor,E013,dArray(2),dArray(3:5),dArray(6:8))
  END IF

  CALL Comm_SummN(dArray,8)

  IF (myid.eq.1) THEN
    WRITE(MTERM,'(A)') "Time Circularity Mass Center RiseVelo RefFrame"
    WRITE(MTERM,'(A,ES12.4,2ES14.6,10ES12.4)') "Stats: ", timens,&
      (dWidth*(2*myPI*dR))/dArray(1),(dWidth*(myPI*dR*dR))/dArray(2),dArray(3:5)/dArray(2),dArray(6:8)/dArray(2),myALE%dFrameVelocity(2)
    WRITE(MFILE,'(A)') "Time Circularity Mass Center RiseVelo RefFrame"
    WRITE(MFILE,'(A,ES12.4,2ES14.6,10ES12.4)') "Stats: ", timens,&
      (dWidth*(2*myPI*dR))/dArray(1),(dWidth*(myPI*dR*dR))/dArray(2),dArray(3:5)/dArray(2),dArray(6:8)/dArray(2),myALE%dFrameVelocity(2)
  END IF

  IF (myALE%bUseFrameVelocity) THEN
    myALE%dFrameVelocityChange = dArray(6:8)/dArray(2)
    myALE%dFrameVelocity = myALE%dFrameVelocity + 1d0*myALE%dFrameVelocityChange
  ELSE
    myALE%dFrameVelocityChange = 0d0
    myALE%dFrameVelocity       = 0d0
  END IF

  IF (myid.eq.1) THEN
    WRITE(MTERM,'(A,3ES12.4)') "ReferenceFrame: ", myALE%dFrameVelocity(2), myALE%dFrameVelocityChange(2)/TSTEP
    WRITE(MFILE,'(A,3ES12.4)') "ReferenceFrame: ", myALE%dFrameVelocity(2), myALE%dFrameVelocityChange(2)/TSTEP
  END IF
end subroutine IntegrateQuantities
!=========================================================================
!
!=========================================================================
SUBROUTINE ExtractElemSizeDistribution()
real*8 x(8),y(8),z(8),dVol
real*8, allocatable :: daux(:),daux2(:)
integer iel,i,nQ2_dof

if (myid.ne.0) then

 ILEV = NLMAX
 CALL SETLEV(2)

 nQ2_dof = mg_mesh%level(ilev)%nvt + &
           mg_mesh%level(ilev)%net + &
           mg_mesh%level(ilev)%nat + &
           mg_mesh%level(ilev)%nel

 allocate(daux(nQ2_dof))
 allocate(daux2(nQ2_dof))
 IF (.not.allocated(ElemSizeDist)) allocate(ElemSizeDist(mg_mesh%level(ilev)%nvt))

 daux2 = 0d0
 daux  = 0d0

 do iel=1,mg_mesh%level(ilev)%nel

  x(:) = mg_mesh%level(ilev)%dcorvg(1,mg_mesh%level(ilev)%kvert(:,iel))
  y(:) = mg_mesh%level(ilev)%dcorvg(2,mg_mesh%level(ilev)%kvert(:,iel))
  z(:) = mg_mesh%level(ilev)%dcorvg(3,mg_mesh%level(ilev)%kvert(:,iel))
  CALL GetElemVol(x,y,z,dVol)

  daux2(mg_mesh%level(ilev)%kvert(:,iel)) = daux2(mg_mesh%level(ilev)%kvert(:,iel)) + 0.125d0*(dVol**0.333d0)
  daux(mg_mesh%level(ilev)%kvert(:,iel))  = daux(mg_mesh%level(ilev)%kvert(:,iel)) + 0.125d0

 end do

 CALL E013Sum(daux2)
 CALL E013Sum(daux)

 do i=1,mg_mesh%level(ilev)%nvt
  ElemSizeDist(i) = daux2(i)/daux(i)
 end do

 deallocate(daux,daux2)

end if

END SUBROUTINE ExtractElemSizeDistribution
!=========================================================================
!
!=========================================================================
SUBROUTINE ExtractVeloGradients()

  ILEV = NLMAX
  CALL SETLEV(2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValU)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,1)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValV)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValW)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,3)

END SUBROUTINE ExtractVeloGradients
!=========================================================================
!
!=========================================================================
SUBROUTINE GetPressureSample(dcorvg,nvt)
  INTEGER NVT
  REAL*8 dcorvg(3,*)
  INTEGER I,J1,J2
  REAL*8 :: P1X=0.55d0,P1Y=0.20d0,P1Z=0.205d0
  REAL*8 :: P2X=0.45d0,P2Y=0.20d0,P2Z=0.205d0
  REAL*8 DIST1,DIST2,MINDIST1,MINDIST2,PX,PY,PZ

  IF (myid.ne.0) THEN
    MINDIST1 = 1d30
    MINDIST2 = 1d30
    DO i=1,nvt
    PX = dcorvg(1,I)
    PY = dcorvg(2,I)
    PZ = dcorvg(3,I)
    DIST1 = SQRT((PX-P1X)**2d0 + (PY-P1Y)**2d0 + (PZ-P1Z)**2d0)
    DIST2 = SQRT((PX-P2X)**2d0 + (PY-P2Y)**2d0 + (PZ-P2Z)**2d0)
    IF (DIST1.LT.MINDIST1) THEN
      MINDIST1 = DIST1
      J1 = I
    END IF
    IF (DIST2.LT.MINDIST2) THEN
      MINDIST2 = DIST2
      J2 = I
    END IF
    END DO
    MINDIST1 = -MINDIST1
    MINDIST2 = -MINDIST2
    DIST1 = MINDIST1
    DIST2 = MINDIST2
    !  WRITE(*,*) myid,J1,J2,MINDIST1,MINDIST2!PressureSample
  END IF

  CALL COMM_Maximum(MINDIST1)
  CALL COMM_Maximum(MINDIST2)

  IF (myid.ne.0) THEN
    !  WRITE(*,*) myid,J1,J2,MINDIST1,MINDIST2!PressureSample
    PressureSample = 0
    IF (DIST1.EQ.MINDIST1) THEN
      PressureSample(1) = J1
    END IF
    IF (DIST2.EQ.MINDIST2) THEN
      PressureSample(2) = J2
    END IF
    !  WRITE(*,*) myid,MINDIST1,MINDIST2,PressureSample
  END IF

END SUBROUTINE GetPressureSample
