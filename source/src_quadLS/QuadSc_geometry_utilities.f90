!=========================================================================
! QuadSc_geometry_utilities.f90
!
! Geometry and FBM (Fictitious Boundary Method) operations
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
! updateFBMGeometry_Wangen - Update FBM geometry using Wangen method
! Computes fictitious boundary method geometry for Wangen approach
!=========================================================================
subroutine updateFBMGeometry_Wangen()
  use cinterface, only: calculateFBM
  implicit none

  if (calculateFBM()) then
    if (myid == showid) write (*, *) '> FBM computation step'

    ILEV = NLMAX
    call SETLEV(2)

    call QuadScalar_FictKnpr_Wangen(mg_mesh%level(ilev)%dcorvg, &
                                    mg_mesh%level(ilev)%dcorag, &
                                    mg_mesh%level(ilev)%kvert, &
                                    mg_mesh%level(ilev)%kedge, &
                                    mg_mesh%level(ilev)%karea)

    call E013Max_SUPER(FictKNPR)

  end if

end subroutine updateFBMGeometry_Wangen

!=========================================================================
! updateFBMGeometry - Update FBM geometry (standard method)
! Computes fictitious boundary method geometry for standard approach
!=========================================================================
subroutine updateFBMGeometry()
  use cinterface, only: calculateFBM
  implicit none

  if (calculateFBM()) then
    if (myid == showid) write (*, *) '> FBM computation step'

    ILEV = NLMAX
    call SETLEV(2)

    call QuadScalar_FictKnpr(mg_mesh%level(ilev)%dcorvg, &
                             mg_mesh%level(ilev)%dcorag, &
                             mg_mesh%level(ilev)%kvert, &
                             mg_mesh%level(ilev)%kedge, &
                             mg_mesh%level(ilev)%karea)

    ! Write a warning:
    call E013Max_SUPER(FictKNPR)

  else
    if (myid == showid) write (*, *) '> FBM disabled'
  end if

end subroutine updateFBMGeometry
!=========================================================================
! updateMixerGeometry - Update mixer/die geometry with distance functions
! Computes signed distance functions for various mixer types:
! SSE (Single Screw Extruder), DIE, NETZSCH, TSE (Twin Screw Extruder)
!=========================================================================
subroutine updateMixerGeometry(mfile)
  use geometry_processing, only: calcDistanceFunction_sse, QuadScalar_MixerKnpr, &
                                 calcDistanceFunction_netzsch, dEpsDist
  implicit none
  integer, intent(in) :: mfile
  integer :: i
  real :: tttt0, tttt1

  call myMPI_Barrier()
  call ZTIME(tttt0)

  ILEV = NLMAX
  call SETLEV(2)
  QuadSc%AuxU = dEpsDist
  QuadSc%AuxV = dEpsDist

  MixerKNPR(:) = 0

  ! Calculate distance function based on simulation type
  if (adjustl(trim(mySigma%cType)) == "SSE" .or. adjustl(trim(mySigma%cType)) == "DIE") then
    call calcDistanceFunction_sse(mg_mesh%level(ilev)%dcorvg, &
                                  mg_mesh%level(ilev)%kvert, &
                                  mg_mesh%level(ilev)%kedge, &
                                  mg_mesh%level(ilev)%karea, &
                                  mg_mesh%level(ilev)%nel, &
                                  mg_mesh%level(ilev)%nvt, &
                                  mg_mesh%level(ilev)%nat, &
                                  mg_mesh%level(ilev)%net, &
                                  QuadSc%AuxU, QuadSc%AuxV, QuadSc%AuxW)
  end if

  if (adjustl(trim(mySigma%cType)) == "NETZSCH") then
    call calcDistanceFunction_netzsch(mg_mesh%level(ilev)%dcorvg, &
                                      mg_mesh%level(ilev)%kvert, &
                                      mg_mesh%level(ilev)%kedge, &
                                      mg_mesh%level(ilev)%karea, &
                                      mg_mesh%level(ilev)%nel, &
                                      mg_mesh%level(ilev)%nvt, &
                                      mg_mesh%level(ilev)%nat, &
                                      mg_mesh%level(ilev)%net, &
                                      QuadSc%AuxU, QuadSc%AuxV, QuadSc%AuxW)
  end if

  if (adjustl(trim(mySigma%cType)) == "TSE") then
    call QuadScalar_MixerKnpr(mg_mesh%level(ilev)%dcorvg, &
                              mg_mesh%level(ilev)%kvert, &
                              mg_mesh%level(ilev)%kedge, &
                              mg_mesh%level(ilev)%karea, &
                              mg_mesh%level(ilev)%nel, &
                              mg_mesh%level(ilev)%nvt, &
                              mg_mesh%level(ilev)%nat, &
                              mg_mesh%level(ilev)%net, &
                              QuadSc%AuxU, QuadSc%AuxV, QuadSc%AuxW)
  end if

  call myMPI_Barrier()
  call ZTIME(tttt1)

  ! Report timing information
  if (myid == 1) then
    write (mterm, "(A,F6.1,A)") "Time used for FINE mesh distance estimation was: ", &
                                tttt1 - tttt0, "s!"
    write (mfile, "(A,F6.1,A)") "Time used for FINE mesh distance estimation was: ", &
                                tttt1 - tttt0, "s!"
  end if

end subroutine updateMixerGeometry
!=========================================================================
! MoveInterfacePoints - Move interface points based on velocity field
! WARNING: This subroutine is untested and currently disabled
!=========================================================================
subroutine MoveInterfacePoints(dcoor, MFILE)
  implicit none
  real(8), intent(inout) :: dcoor(3, *)
  integer, intent(in) :: mfile
  real(8) :: Velo(3), Displacement(3), dMaxVelo, daux, dArea
  integer :: i, iInterface

  ! WARNING: Untested subroutine - disabled for safety
  write (*, *) 'ERROR: MoveInterfacePoints is untested and currently disabled'
  stop

  if (myid /= 0) then
    do i = 1, QuadSc%ndof
      Velo = [QuadSc%ValU(i), QuadSc%ValV(i), QuadSc%ValW(i)]
      Displacement = 1.0d0 * tstep * Velo
      myQ2Coor(:, i) = myQ2Coor(:, i) + Displacement
    end do
  end if

  if (.not. allocated(myTSurf)) allocate (myTSurf(Properties%nInterface))

  if (myid /= 0) then
    ILEV = NLMAX
    call SETLEV(2)
    do i = 1, QuadSc%ndof
      dcoor(:, i) = myQ2Coor(:, i)
    end do

    do iInterface = 1, Properties%nInterface
      call BuildUpTriangulation(KWORK(L(LVERT)), KWORK(L(LEDGE)), KWORK(L(LAREA)), &
                                myQ2Coor, iInterface)
    end do
  end if

  call CommunicateSurface()

  if (myid /= 0) then
    ILEV = NLMAX
    call SETLEV(2)
    do i = 1, QuadSc%ndof
      dcoor(:, i) = myALE%Q2coor_old(:, i)
      myQ2Coor(:, i) = myALE%Q2coor_old(:, i)
    end do
  end if

  ! IF (myid.eq.1) THEN
  !  CALL GetCompleteArea(dArea,1)
  !  WRITE(MFILE,'(A,3ES14.6)') "CompleteSurfaceAreaAndCircularity: ", TIMENS,dArea, (0.05*(2*(4d0*ATAN(1d0))*0.25d0))/dArea
  !  WRITE(MTERM,'(A,3ES14.6)') "CompleteSurfaceAreaAndCircularity: ", TIMENS,dArea, (0.05*(2*(4d0*ATAN(1d0))*0.25d0))/dArea
  ! END IF

end subroutine MoveInterfacePoints

!=========================================================================
! GetCompleteArea - Calculate complete surface area of an interface
! WARNING: This subroutine is untested and currently disabled
!=========================================================================
subroutine GetCompleteArea(DCompleteArea, iIF)
  implicit none
  real(8), intent(out) :: DCompleteArea
  integer, intent(in) :: iIF
  real(8) :: DA(3, 3), dArea
  integer :: i, j, IP1, IP2, IP3

  ! WARNING: Untested subroutine - disabled for safety
  write (*, *) 'ERROR: GetCompleteArea is untested and currently disabled'
  stop

  DCompleteArea = 0.0d0
  do i = 1, myTSurf(iIF)%nT
    do j = 1, 8
      IP1 = j
      IP2 = mod(j, 8) + 1
      IP3 = 9

      ! Compute edge vectors for triangle (P2-P1, P3-P1)
      DA(1, 2) = myTSurf(iIF)%T(i)%C(1, IP2) - myTSurf(iIF)%T(i)%C(1, IP1) ! P2X-P1X
      DA(2, 2) = myTSurf(iIF)%T(i)%C(2, IP2) - myTSurf(iIF)%T(i)%C(2, IP1) ! P2Y-P1Y
      DA(3, 2) = myTSurf(iIF)%T(i)%C(3, IP2) - myTSurf(iIF)%T(i)%C(3, IP1) ! P2Z-P1Z
      DA(1, 3) = myTSurf(iIF)%T(i)%C(1, IP3) - myTSurf(iIF)%T(i)%C(1, IP1) ! P3X-P1X
      DA(2, 3) = myTSurf(iIF)%T(i)%C(2, IP3) - myTSurf(iIF)%T(i)%C(2, IP1) ! P3Y-P1Y
      DA(3, 3) = myTSurf(iIF)%T(i)%C(3, IP3) - myTSurf(iIF)%T(i)%C(3, IP1) ! P3Z-P1Z

      ! Compute cross product to get normal vector
      DA(1, 1) = DA(2, 3) * DA(3, 2) - DA(3, 3) * DA(2, 2)
      DA(2, 1) = DA(3, 3) * DA(1, 2) - DA(1, 3) * DA(3, 2)
      DA(3, 1) = DA(1, 3) * DA(2, 2) - DA(2, 3) * DA(1, 2)

      ! Area = 0.5 * |cross product|
      dArea = 0.5d0 * sqrt(DA(1, 1) * DA(1, 1) + DA(2, 1) * DA(2, 1) + DA(3, 1) * DA(3, 1))
      DCompleteArea = DCompleteArea + dArea
    end do
  end do

end subroutine GetCompleteArea

!=========================================================================
! BuildUpTriangulation - Build Q2 surface triangulation for an interface
! Creates triangulated surface representation from Q2 finite element mesh
! WARNING: This subroutine is untested and currently disabled
!=========================================================================
subroutine BuildUpTriangulation(kvert, kedge, karea, dcorvg, iIF)
  implicit none
  integer, intent(in) :: kvert(8, *), kedge(12, *), karea(6, *), iIF
  real(8), intent(in) :: dcorvg(3, *)
  integer :: iel, i, j, k, ivt1, ivt2, ivt3, ivt4, ivt5, iT

  ! Element connectivity arrays for hexahedral Q2 elements
  integer :: NeighA(4, 6), NeighU(4, 6)
  data NeighA/1, 2, 3, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 4, 8, 7, 4, 1, 5, 8, 5, 6, 7, 8/
  data NeighU/1, 2, 3, 4, 1, 6, 9, 5, 2, 7, 10, 6, 3, 8, 11, 7, 4, 5, 12, 8, 9, 10, 11, 12/

  ! WARNING: Untested subroutine - disabled for safety
  write (*, *) 'ERROR: BuildUpTriangulation is untested and currently disabled'
  stop

  ! Count interface triangles
  iT = 0
  k = 1
  do i = 1, nel
    do j = 1, 6
      if (k == karea(j, i)) then
        if (myBoundary%LS_zero(nvt + net + k) == iIF) then
          iT = iT + 1
        end if
        k = k + 1
      end if
    end do
  end do

  ! Allocate triangulation structure
  if (allocated(myTSurf(iIF)%T)) deallocate (myTSurf(iIF)%T)

  myTSurf(iIF)%nT = iT
  allocate (myTSurf(iIF)%T(myTSurf(iIF)%nT))

  ! Build triangulation with Q2 coordinates
  iT = 0
  k = 1
  do i = 1, nel
    do j = 1, 6
      if (k == karea(j, i)) then
        if (myBoundary%LS_zero(nvt + net + k) == iIF) then
          iT = iT + 1

          ! Store corner vertices (Q2 nodes 1,3,5,7)
          ivt1 = kvert(NeighA(1, j), i)
          ivt2 = kvert(NeighA(2, j), i)
          ivt3 = kvert(NeighA(3, j), i)
          ivt4 = kvert(NeighA(4, j), i)
          myTSurf(iIF)%T(iT)%C(:, 1) = dcorvg(:, ivt1)
          myTSurf(iIF)%T(iT)%C(:, 3) = dcorvg(:, ivt2)
          myTSurf(iIF)%T(iT)%C(:, 5) = dcorvg(:, ivt3)
          myTSurf(iIF)%T(iT)%C(:, 7) = dcorvg(:, ivt4)

          ! Store edge midpoint nodes (Q2 nodes 2,4,6,8)
          ivt1 = kedge(NeighU(1, j), i)
          ivt2 = kedge(NeighU(2, j), i)
          ivt3 = kedge(NeighU(3, j), i)
          ivt4 = kedge(NeighU(4, j), i)
          myTSurf(iIF)%T(iT)%C(:, 2) = dcorvg(:, nvt + ivt1)
          myTSurf(iIF)%T(iT)%C(:, 4) = dcorvg(:, nvt + ivt2)
          myTSurf(iIF)%T(iT)%C(:, 6) = dcorvg(:, nvt + ivt3)
          myTSurf(iIF)%T(iT)%C(:, 8) = dcorvg(:, nvt + ivt4)

          ! Store face center node (Q2 node 9)
          myTSurf(iIF)%T(iT)%C(:, 9) = dcorvg(:, nvt + net + k)
        end if
        k = k + 1
      end if
    end do
  end do

end subroutine BuildUpTriangulation
