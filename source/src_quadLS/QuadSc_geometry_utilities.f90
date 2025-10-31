!=========================================================================
! QuadSc_geometry_utilities.f90
!
! Geometry and FBM (Fictitious Boundary Method) operations
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE updateFBMGeometry_Wangen()
use cinterface, only: calculateFBM

 if (calculateFBM()) then
  if (myid.eq.showid) write(*,*) '> FBM computation step'

  ILEV=NLMAX
  CALL SETLEV(2)

  CALL QuadScalar_FictKnpr_Wangen(mg_mesh%level(ilev)%dcorvg,&
    mg_mesh%level(ilev)%dcorag,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%karea)

  CALL E013Max_SUPER(FictKNPR)

 end if

END SUBROUTINE  updateFBMGeometry_Wangen
!=========================================================================
!
!=========================================================================
SUBROUTINE updateFBMGeometry()
use cinterface, only: calculateFBM

 if (calculateFBM()) then
  if (myid.eq.showid) write(*,*) '> FBM computation step'

  ILEV=NLMAX
  CALL SETLEV(2)

  CALL QuadScalar_FictKnpr(mg_mesh%level(ilev)%dcorvg,&
    mg_mesh%level(ilev)%dcorag,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%karea)

  ! Write a warning:
  CALL E013Max_SUPER(FictKNPR)

 else
  if (myid.eq.showid) write(*,*) '> FBM disabled'
 end if

END SUBROUTINE  updateFBMGeometry
!=========================================================================
!
!=========================================================================
SUBROUTINE updateMixerGeometry(mfile)
use geometry_processing, only : calcDistanceFunction_sse, QuadScalar_MixerKnpr,&
    calcDistanceFunction_netzsch,dEpsDist

integer, intent(in) :: mfile

integer :: i

REAL :: tttt0,tttt1

CALL myMPI_Barrier()
CALL ZTIME(tttt0)

ILEV=NLMAX
CALL SETLEV(2)
QuadSc%AuxU = dEpsDist
QuadSc%AuxV = dEpsDist

MixerKNPR(:) = 0

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 CALL calcDistanceFunction_sse(mg_mesh%level(ilev)%dcorvg,&
                           mg_mesh%level(ilev)%kvert,&
                           mg_mesh%level(ilev)%kedge,&
                           mg_mesh%level(ilev)%karea,&
                           mg_mesh%level(ilev)%nel,&
                           mg_mesh%level(ilev)%nvt,&
                           mg_mesh%level(ilev)%nat,&
                           mg_mesh%level(ilev)%net,&
                           QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."NETZSCH") THEN
 CALL calcDistanceFunction_netzsch(mg_mesh%level(ilev)%dcorvg,&
                           mg_mesh%level(ilev)%kvert,&
                           mg_mesh%level(ilev)%kedge,&
                           mg_mesh%level(ilev)%karea,&
                           mg_mesh%level(ilev)%nel,&
                           mg_mesh%level(ilev)%nvt,&
                           mg_mesh%level(ilev)%nat,&
                           mg_mesh%level(ilev)%net,&
                           QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
 CALL QuadScalar_MixerKnpr(mg_mesh%level(ilev)%dcorvg,&
                           mg_mesh%level(ilev)%kvert,&
                           mg_mesh%level(ilev)%kedge,&
                           mg_mesh%level(ilev)%karea,&
                           mg_mesh%level(ilev)%nel,&
                           mg_mesh%level(ilev)%nvt,&
                           mg_mesh%level(ilev)%nat,&
                           mg_mesh%level(ilev)%net,&
                           QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW)
END IF


CALL myMPI_Barrier()
CALL ZTIME(tttt1)
IF (myid.eq.1) WRITE(mterm,"(A,F6.1,A)") "Time used for FINE mesh distance estimation was: ", tttt1-tttt0, "s!"
IF (myid.eq.1) WRITE(mfile,"(A,F6.1,A)") "Time used for FINE mesh distance estimation was: ", tttt1-tttt0, "s!"

END SUBROUTINE  updateMixerGeometry
!=========================================================================
!
!=========================================================================
SUBROUTINE MoveInterfacePoints(dcoor,MFILE)
  INTEGER mfile
  REAL*8 dcoor(3,*)
  REAL*8 Velo(3),Displacement(3),dMaxVelo,daux,dArea
  INTEGER i,iInterface

  write(*,*)'New untested subroutine'
  stop

  IF (myid.ne.0) then
    DO i=1,QuadSc%ndof
    Velo = [QuadSc%ValU(i),QuadSc%ValV(i),QuadSc%ValW(i)]
    Displacement = 1d0*tstep*Velo
    myQ2Coor(:,i) = myQ2Coor(:,i) + Displacement
    END DO
  END IF

  IF (.NOT.ALLOCATED(myTSurf)) ALLOCATE(myTSurf(Properties%nInterface))

  IF (myid.ne.0) then
    ILEV=NLMAX
    CALL SETLEV(2)
    DO i=1,QuadSc%ndof
    dcoor(:,i) = myQ2Coor(:,i)
    END DO

    DO iInterface=1,Properties%nInterface
    CALL BuildUpTriangulation(KWORK(L(LVERT)),KWORK(L(LEDGE)),KWORK(L(LAREA)),myQ2Coor,iInterface)
    END DO
  END IF

  CALL CommunicateSurface()

  IF (myid.ne.0) then
    ILEV=NLMAX
    CALL SETLEV(2)
    DO i=1,QuadSc%ndof
    dcoor(:,i) = myALE%Q2coor_old(:,i)
    myQ2Coor(:,i) = myALE%Q2coor_old(:,i)
    END DO
  END IF

  ! IF (myid.eq.1) THEN
  !  CALL GetCompleteArea(dArea,1)
  !  WRITE(MFILE,'(A,3ES14.6)') "CompleteSurfaceAreaAndCircularity: ", TIMENS,dArea, (0.05*(2*(4d0*ATAN(1d0))*0.25d0))/dArea
  !  WRITE(MTERM,'(A,3ES14.6)') "CompleteSurfaceAreaAndCircularity: ", TIMENS,dArea, (0.05*(2*(4d0*ATAN(1d0))*0.25d0))/dArea
  ! END IF

END SUBROUTINE MoveInterfacePoints
!=========================================================================
!
!=========================================================================
SUBROUTINE GetCompleteArea(DCompleteArea,iIF)
  REAL*8 DCompleteArea
  REAL*8 DA(3,3),dArea
  INTEGER i,j,IP1,IP2,IP3,iIF

  write(*,*)'New untested subroutine'
  stop

  DCompleteArea = 0d0
  DO i=1,myTSurf(iIF)%nT
  DO j=1,8
  IP1 = j
  IP2 = MOD(j,8)+1
  IP3 = 9
  DA(1,2)=myTSurf(iIF)%T(i)%C(1,IP2) - myTSurf(iIF)%T(i)%C(1,IP1) !P2X-P1X
  DA(2,2)=myTSurf(iIF)%T(i)%C(2,IP2) - myTSurf(iIF)%T(i)%C(2,IP1) !P2Y-P1Y
  DA(3,2)=myTSurf(iIF)%T(i)%C(3,IP2) - myTSurf(iIF)%T(i)%C(3,IP1) !P2Z-P1Z
  DA(1,3)=myTSurf(iIF)%T(i)%C(1,IP3) - myTSurf(iIF)%T(i)%C(1,IP1) !P3X-P1X
  DA(2,3)=myTSurf(iIF)%T(i)%C(2,IP3) - myTSurf(iIF)%T(i)%C(2,IP1) !P3Y-P1Y
  DA(3,3)=myTSurf(iIF)%T(i)%C(3,IP3) - myTSurf(iIF)%T(i)%C(3,IP1) !P3Z-P1Z

  DA(1,1) = DA(2,3)*DA(3,2) - DA(3,3)*DA(2,2)
  DA(2,1) = DA(3,3)*DA(1,2) - DA(1,3)*DA(3,2)
  DA(3,1) = DA(1,3)*DA(2,2) - DA(2,3)*DA(1,2)
  dArea = 0.5*SQRT(DA(1,1)*DA(1,1) + DA(2,1)*DA(2,1) + DA(3,1)*DA(3,1))
  DCompleteArea = DCompleteArea + dArea
  END DO
  END DO

END SUBROUTINE GetCompleteArea
!=========================================================================
!
!=========================================================================
subroutine BuildUpTriangulation(kvert,kedge,karea,dcorvg,iIF)
  REAL*8 dcorvg(3,*)
  INTEGER karea(6,*),kvert(8,*),kedge(12,*),iIF
  INTEGER iel,i,j,k,ivt1,ivt2,ivt3,ivt4,ivt5,iT
  INTEGER NeighA(4,6),NeighU(4,6)
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
  DATA NeighU/1,2,3,4, 1,6,9,5, 2,7,10,6, 3,8,11,7, 4,5,12,8, 9,10,11,12/

  write(*,*)'New untested subroutine'
  stop

  iT = 0
  k=1
  DO i=1,nel
  DO j=1,6
  IF (k.eq.karea(j,i)) THEN
    IF (myBoundary%LS_zero(nvt+net+k).eq.iIF) THEN
      iT = iT + 1
    END IF
    k = k + 1
  END IF
  END DO
  END DO

  IF (ALLOCATED(myTSurf(iIF)%T)) DEALLOCATE(myTSurf(iIF)%T)

  myTSurf(iIF)%nT = iT
  ALLOCATE(myTSurf(iIF)%T(myTSurf(iIF)%nT))

  iT = 0
  k=1
  DO i=1,nel
  DO j=1,6
  IF (k.eq.karea(j,i)) THEN
    IF (myBoundary%LS_zero(nvt+net+k).eq.iIF) THEN
      iT = iT + 1
      ivt1 = kvert(NeighA(1,j),i)
      ivt2 = kvert(NeighA(2,j),i)
      ivt3 = kvert(NeighA(3,j),i)
      ivt4 = kvert(NeighA(4,j),i)
      myTSurf(iIF)%T(iT)%C(:,1) = dcorvg(:,ivt1)
      myTSurf(iIF)%T(iT)%C(:,3) = dcorvg(:,ivt2)
      myTSurf(iIF)%T(iT)%C(:,5) = dcorvg(:,ivt3)
      myTSurf(iIF)%T(iT)%C(:,7) = dcorvg(:,ivt4)
      ivt1 = kedge(NeighU(1,j),i)
      ivt2 = kedge(NeighU(2,j),i)
      ivt3 = kedge(NeighU(3,j),i)
      ivt4 = kedge(NeighU(4,j),i)
      myTSurf(iIF)%T(iT)%C(:,2) = dcorvg(:,nvt+ivt1)
      myTSurf(iIF)%T(iT)%C(:,4) = dcorvg(:,nvt+ivt2)
      myTSurf(iIF)%T(iT)%C(:,6) = dcorvg(:,nvt+ivt3)
      myTSurf(iIF)%T(iT)%C(:,8) = dcorvg(:,nvt+ivt4)
      myTSurf(iIF)%T(iT)%C(:,9) = dcorvg(:,nvt+net+k)
    END IF
    k = k + 1
  END IF
  END DO
  END DO

end subroutine BuildUpTriangulation
