!=========================================================================
! QuadSc_integration.f90
!
! Integration and analysis utility functions
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================

!=========================================================================
! Integrate_DIE_Flowrate - Integrate volumetric flowrate for DIE simulation
! Computes flowrate through inflow (iPar=0) or outflow (iPar=1) boundaries
! using face-based integration with normal area vectors
!=========================================================================
subroutine Integrate_DIE_Flowrate(dcorvg, karea, kvert, nel, dVolFlow, iPar)
  implicit none
  real(8), intent(in) :: dcorvg(3, *)
  integer, intent(in) :: karea(6, *), kvert(8, *), nel, iPar
  real(8), intent(out) :: dVolFlow
  integer :: i, j, k, ivt1, ivt2, ivt3, ivt4
  real(8) :: P(3), dAN(3), dV

  ! Face connectivity for hexahedral elements
  integer :: NeighA(4, 6)
  data NeighA/1, 2, 3, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 4, 8, 7, 4, 1, 5, 8, 5, 6, 7, 8/

  dVolFlow = 0.0d0

  ! Integrate over inflow boundary
  if (iPar == 0) then
    k = 1
    do i = 1, nel
      do j = 1, 6
        if (k == karea(j, i)) then
          if (myBoundary%iInflow(nvt + net + k) /= 0) then
            ivt1 = kvert(NeighA(1, j), i)
            ivt2 = kvert(NeighA(2, j), i)
            ivt3 = kvert(NeighA(3, j), i)
            ivt4 = kvert(NeighA(4, j), i)

            ! Compute normal area vector for this face
            call GET_NormalArea(dcorvg(1:3, ivt1), dcorvg(1:3, ivt2), &
                                dcorvg(1:3, ivt3), dcorvg(1:3, ivt4), &
                                dcorvg(1:3, nvt + net + nat + i), dAN)

            ! Compute flux: Q = A·v (dot product of area vector and velocity)
            dV = dAN(1) * QuadSc%ValU(nvt + net + k) + &
                 dAN(2) * QuadSc%ValV(nvt + net + k) + &
                 dAN(3) * QuadSc%ValW(nvt + net + k)

            dVolFlow = dVolFlow + dV
          end if
          k = k + 1
        end if
      end do
    end do
  end if

  ! Integrate over outflow boundary
  if (iPar == 1) then
    k = 1
    do i = 1, nel
      do j = 1, 6
        if (k == karea(j, i)) then
          if (myBoundary%bOutflow(nvt + net + k)) then
            ivt1 = kvert(NeighA(1, j), i)
            ivt2 = kvert(NeighA(2, j), i)
            ivt3 = kvert(NeighA(3, j), i)
            ivt4 = kvert(NeighA(4, j), i)

            ! Compute normal area vector (flip direction for outflow)
            call GET_NormalArea(dcorvg(1:3, ivt1), dcorvg(1:3, ivt2), &
                                dcorvg(1:3, ivt3), dcorvg(1:3, ivt4), &
                                dcorvg(1:3, nvt + net + nat + i), dAN)
            dAN = -dAN

            ! Compute flux: Q = A·v
            dV = dAN(1) * QuadSc%ValU(nvt + net + k) + &
                 dAN(2) * QuadSc%ValV(nvt + net + k) + &
                 dAN(3) * QuadSc%ValW(nvt + net + k)

            dVolFlow = dVolFlow + dV
          end if
          k = k + 1
        end if
      end do
    end do
  end if

end subroutine Integrate_DIE_Flowrate

!=========================================================================
! IntegrateFlowrate - Integrate flowrate through a z-plane
! Computes volumetric flowrate through faces at specified z-coordinate
!=========================================================================
subroutine IntegrateFlowrate(dcorvg, karea, kvert, nel, dVolFlow, dPar)
  implicit none
  real(8), intent(in) :: dcorvg(3, *), dPar
  integer, intent(in) :: karea(6, *), kvert(8, *), nel
  real(8), intent(out) :: dVolFlow
  integer :: i, j, k, ivt1, ivt2, ivt3, ivt4
  real(8) :: P(3), dA

  ! Face connectivity for hexahedral elements
  integer :: NeighA(4, 6)
  data NeighA/1, 2, 3, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 4, 8, 7, 4, 1, 5, 8, 5, 6, 7, 8/

  dVolFlow = 0.0d0

  k = 1
  do i = 1, nel
    do j = 1, 6
      if (k == karea(j, i)) then
        ivt1 = kvert(NeighA(1, j), i)
        ivt2 = kvert(NeighA(2, j), i)
        ivt3 = kvert(NeighA(3, j), i)
        ivt4 = kvert(NeighA(4, j), i)

        ! Check if all four vertices lie on the z-plane (within tolerance)
        if (abs(dcorvg(3, ivt1) - dPar) < 1.0d-4 .and. &
            abs(dcorvg(3, ivt2) - dPar) < 1.0d-4 .and. &
            abs(dcorvg(3, ivt3) - dPar) < 1.0d-4 .and. &
            abs(dcorvg(3, ivt4) - dPar) < 1.0d-4) then

          ! Compute face area and integrate w-velocity
          call GET_area(dcorvg(1:3, ivt1), dcorvg(1:3, ivt2), &
                        dcorvg(1:3, ivt3), dcorvg(1:3, ivt4), dA)
          dVolFlow = dVolFlow + dA * QuadSc%ValW(nvt + net + k)
        end if
        k = k + 1
      end if
    end do
  end do

end subroutine IntegrateFlowrate

!=========================================================================
! IntegratePressure - Integrate pressure through a z-plane
! Computes area-weighted average pressure at specified z-coordinate
!=========================================================================
subroutine IntegratePressure(dcorvg, karea, kvert, nel, dIntPres, dArea, dPar)
  implicit none
  real(8), intent(in) :: dcorvg(3, *), dPar
  integer, intent(in) :: karea(6, *), kvert(8, *), nel
  real(8), intent(out) :: dIntPres, dArea
  integer :: i, j, k, ivt1, ivt2, ivt3, ivt4
  real(8) :: P(3), dA

  ! Face connectivity for hexahedral elements
  integer :: NeighA(4, 6)
  data NeighA/1, 2, 3, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 4, 8, 7, 4, 1, 5, 8, 5, 6, 7, 8/

  dIntPres = 0.0d0
  dArea = 0.0d0

  k = 1
  do i = 1, nel
    do j = 1, 6
      if (k == karea(j, i)) then
        ivt1 = kvert(NeighA(1, j), i)
        ivt2 = kvert(NeighA(2, j), i)
        ivt3 = kvert(NeighA(3, j), i)
        ivt4 = kvert(NeighA(4, j), i)

        ! Check if all four vertices lie on the z-plane (within tolerance)
        if (abs(dcorvg(3, ivt1) - dPar) < 1.0d-4 .and. &
            abs(dcorvg(3, ivt2) - dPar) < 1.0d-4 .and. &
            abs(dcorvg(3, ivt3) - dPar) < 1.0d-4 .and. &
            abs(dcorvg(3, ivt4) - dPar) < 1.0d-4) then

          ! Compute face area and integrate pressure (P1 element value)
          call GET_area(dcorvg(1:3, ivt1), dcorvg(1:3, ivt2), &
                        dcorvg(1:3, ivt3), dcorvg(1:3, ivt4), dA)
          dIntPres = dIntPres + dA * LinSc%ValP(NLMAX)%x(4 * (i - 1) + 1)
          dArea = dArea + dA
        end if
        k = k + 1
      end if
    end do
  end do

end subroutine IntegratePressure

!=========================================================================
! IntegratePressureAtInflow - Integrate pressure at inflow boundary
! Computes area-weighted average pressure at inflow faces
!=========================================================================
subroutine IntegratePressureAtInflow(dcorvg, karea, kvert, nel, dIntPres, dArea)
  implicit none
  real(8), intent(in) :: dcorvg(3, *)
  integer, intent(in) :: karea(6, *), kvert(8, *), nel
  real(8), intent(out) :: dIntPres, dArea
  integer :: i, j, k, ivt1, ivt2, ivt3, ivt4
  real(8) :: P(3), dA

  ! Face connectivity for hexahedral elements
  integer :: NeighA(4, 6)
  data NeighA/1, 2, 3, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 4, 8, 7, 4, 1, 5, 8, 5, 6, 7, 8/

  dIntPres = 0.0d0
  dArea = 0.0d0

  k = 1
  do i = 1, nel
    do j = 1, 6
      if (k == karea(j, i)) then
        ivt1 = kvert(NeighA(1, j), i)
        ivt2 = kvert(NeighA(2, j), i)
        ivt3 = kvert(NeighA(3, j), i)
        ivt4 = kvert(NeighA(4, j), i)

        ! Check if all four vertices are on inflow boundary
        if (myBoundary%iInflow(ivt1) < 0 .and. myBoundary%iInflow(ivt2) < 0 .and. &
            myBoundary%iInflow(ivt3) < 0 .and. myBoundary%iInflow(ivt4) < 0) then

          ! Compute face area and integrate pressure
          call GET_area(dcorvg(1:3, ivt1), dcorvg(1:3, ivt2), &
                        dcorvg(1:3, ivt3), dcorvg(1:3, ivt4), dA)
          dIntPres = dIntPres + dA * LinSc%ValP(NLMAX)%x(4 * (i - 1) + 1)
          dArea = dArea + dA
        end if
        k = k + 1
      end if
    end do
  end do

end subroutine IntegratePressureAtInflow

!=========================================================================
! IntegrateQuantities - Integrate surface, volume, and mass properties
! Computes surface area, volume, center of mass, and velocities for ALE
! WARNING: This subroutine is untested and currently disabled
!=========================================================================
subroutine IntegrateQuantities(mfile)
  implicit none
  integer, intent(in) :: mfile
  real(8) :: dArray(8)
  real(8) :: dR = 0.25d0, myPI = 4.0d0 * atan(1.0d0), dWidth = 0.05d0
  external E013

  ! WARNING: Untested subroutine - disabled for safety
  write (*, *) 'ERROR: IntegrateQuantities is untested and currently disabled'
  stop

  if (myid /= 0) then
    dArray = 0.0d0
    ILEV = NLMAX
    call SETLEV(2)
    call GetSurface(KWORK(L(LVERT)), KWORK(L(LAREA)), KWORK(L(LEDGE)), &
                    myQ2coor, E013, dArray(1))
    call GetVolume(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                   KWORK(L(LVERT)), KWORK(L(LAREA)), KWORK(L(LEDGE)), &
                   myQ2coor, E013, dArray(2), dArray(3:5), dArray(6:8))
  end if

  call Comm_SummN(dArray, 8)

  if (myid == 1) then
    write (MTERM, '(A)') "Time Circularity Mass Center RiseVelo RefFrame"
    write (MTERM, '(A,ES12.4,2ES14.6,10ES12.4)') "Stats: ", timens, &
      (dWidth * (2.0d0 * myPI * dR)) / dArray(1), &
      (dWidth * (myPI * dR * dR)) / dArray(2), &
      dArray(3:5) / dArray(2), dArray(6:8) / dArray(2), myALE%dFrameVelocity(2)
    write (MFILE, '(A)') "Time Circularity Mass Center RiseVelo RefFrame"
    write (MFILE, '(A,ES12.4,2ES14.6,10ES12.4)') "Stats: ", timens, &
      (dWidth * (2.0d0 * myPI * dR)) / dArray(1), &
      (dWidth * (myPI * dR * dR)) / dArray(2), &
      dArray(3:5) / dArray(2), dArray(6:8) / dArray(2), myALE%dFrameVelocity(2)
  end if

  ! Update ALE reference frame velocity
  if (myALE%bUseFrameVelocity) then
    myALE%dFrameVelocityChange = dArray(6:8) / dArray(2)
    myALE%dFrameVelocity = myALE%dFrameVelocity + 1.0d0 * myALE%dFrameVelocityChange
  else
    myALE%dFrameVelocityChange = 0.0d0
    myALE%dFrameVelocity = 0.0d0
  end if

  if (myid == 1) then
    write (MTERM, '(A,3ES12.4)') "ReferenceFrame: ", myALE%dFrameVelocity(2), &
                                  myALE%dFrameVelocityChange(2) / TSTEP
    write (MFILE, '(A,3ES12.4)') "ReferenceFrame: ", myALE%dFrameVelocity(2), &
                                  myALE%dFrameVelocityChange(2) / TSTEP
  end if

end subroutine IntegrateQuantities

!=========================================================================
! ExtractElemSizeDistribution - Compute characteristic element size at vertices
! Computes h = V^(1/3) averaged over elements adjacent to each vertex
! Used for error estimation and adaptive mesh refinement
!=========================================================================
subroutine ExtractElemSizeDistribution()
  implicit none
  real(8) :: x(8), y(8), z(8), dVol
  real(8), allocatable :: daux(:), daux2(:)
  integer :: iel, i, nQ2_dof

  if (myid /= 0) then

    ILEV = NLMAX
    call SETLEV(2)

    ! Total Q2 degrees of freedom
    nQ2_dof = mg_mesh%level(ilev)%nvt + &
              mg_mesh%level(ilev)%net + &
              mg_mesh%level(ilev)%nat + &
              mg_mesh%level(ilev)%nel

    allocate (daux(nQ2_dof))
    allocate (daux2(nQ2_dof))
    if (.not. allocated(ElemSizeDist)) allocate (ElemSizeDist(mg_mesh%level(ilev)%nvt))

    daux2 = 0.0d0
    daux = 0.0d0

    ! Accumulate element size contributions
    do iel = 1, mg_mesh%level(ilev)%nel

      x(:) = mg_mesh%level(ilev)%dcorvg(1, mg_mesh%level(ilev)%kvert(:, iel))
      y(:) = mg_mesh%level(ilev)%dcorvg(2, mg_mesh%level(ilev)%kvert(:, iel))
      z(:) = mg_mesh%level(ilev)%dcorvg(3, mg_mesh%level(ilev)%kvert(:, iel))
      call GetElemVol(x, y, z, dVol)

      ! Accumulate h = V^(1/3) and count for averaging (0.125 = 1/8 vertices per element)
      daux2(mg_mesh%level(ilev)%kvert(:, iel)) = &
        daux2(mg_mesh%level(ilev)%kvert(:, iel)) + 0.125d0 * (dVol**0.333d0)
      daux(mg_mesh%level(ilev)%kvert(:, iel)) = &
        daux(mg_mesh%level(ilev)%kvert(:, iel)) + 0.125d0

    end do

    ! Sum contributions from all MPI processes
    call E013Sum(daux2)
    call E013Sum(daux)

    ! Compute average element size at each vertex
    do i = 1, mg_mesh%level(ilev)%nvt
      ElemSizeDist(i) = daux2(i) / daux(i)
    end do

    deallocate (daux, daux2)

  end if

end subroutine ExtractElemSizeDistribution

!=========================================================================
! ExtractVeloGradients - Compute velocity gradients at all Q2 nodes
! Computes ∂u/∂x, ∂u/∂y, ∂u/∂z, ∂v/∂x, ∂v/∂y, ∂v/∂z, ∂w/∂x, ∂w/∂y, ∂w/∂z
! Stored in QuadSc%ValUx, ValUy, ValUz, ValVx, ValVy, ValVz, ValWx, ValWy, ValWz
!=========================================================================
subroutine ExtractVeloGradients()
  implicit none

  ILEV = NLMAX
  call SETLEV(2)

  ! Compute gradients of U component
  call GetGradVelo_rhs(QuadSc, QuadSc%ValU)
  call E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)
  call GetGradVelo_val(QuadSc, 1)

  ! Compute gradients of V component
  call GetGradVelo_rhs(QuadSc, QuadSc%ValV)
  call E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)
  call GetGradVelo_val(QuadSc, 2)

  ! Compute gradients of W component
  call GetGradVelo_rhs(QuadSc, QuadSc%ValW)
  call E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)
  call GetGradVelo_val(QuadSc, 3)

end subroutine ExtractVeloGradients

!=========================================================================
! GetPressureSample - Find nearest vertices to two sampling points
! Locates vertices closest to specified probe points for pressure monitoring
! Uses MPI reduction to find global minimum distance across all processes
!=========================================================================
subroutine GetPressureSample(dcorvg, nvt)
  implicit none
  real(8), intent(in) :: dcorvg(3, *)
  integer, intent(in) :: nvt
  integer :: i, J1, J2
  real(8) :: DIST1, DIST2, MINDIST1, MINDIST2, PX, PY, PZ

  ! Sampling point coordinates
  real(8), parameter :: P1X = 0.55d0, P1Y = 0.20d0, P1Z = 0.205d0
  real(8), parameter :: P2X = 0.45d0, P2Y = 0.20d0, P2Z = 0.205d0

  if (myid /= 0) then
    MINDIST1 = 1.0d30
    MINDIST2 = 1.0d30

    ! Find closest vertices to each sampling point
    do i = 1, nvt
      PX = dcorvg(1, i)
      PY = dcorvg(2, i)
      PZ = dcorvg(3, i)

      ! Distance to first sampling point
      DIST1 = sqrt((PX - P1X)**2.0d0 + (PY - P1Y)**2.0d0 + (PZ - P1Z)**2.0d0)
      if (DIST1 < MINDIST1) then
        MINDIST1 = DIST1
        J1 = i
      end if

      ! Distance to second sampling point
      DIST2 = sqrt((PX - P2X)**2.0d0 + (PY - P2Y)**2.0d0 + (PZ - P2Z)**2.0d0)
      if (DIST2 < MINDIST2) then
        MINDIST2 = DIST2
        J2 = i
      end if
    end do

    ! Negate for MPI maximum reduction (finds minimum)
    MINDIST1 = -MINDIST1
    MINDIST2 = -MINDIST2
    DIST1 = MINDIST1
    DIST2 = MINDIST2
  end if

  ! Find global minimum distances across all processes
  call COMM_Maximum(MINDIST1)
  call COMM_Maximum(MINDIST2)

  if (myid /= 0) then
    PressureSample = 0

    ! Assign vertex indices if this process owns the closest vertices
    if (DIST1 == MINDIST1) then
      PressureSample(1) = J1
    end if
    if (DIST2 == MINDIST2) then
      PressureSample(2) = J2
    end if
  end if

end subroutine GetPressureSample
