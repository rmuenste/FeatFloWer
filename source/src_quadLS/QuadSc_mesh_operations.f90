!=========================================================================
! QuadSc_mesh_operations.f90
!
! Mesh deformation and ALE (Arbitrary Lagrangian-Eulerian) operations
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================

!=========================================================================
! STORE_OLD_MESH - Store old mesh coordinates for ALE
! Saves current mesh coordinates for computing mesh velocity
!=========================================================================
subroutine STORE_OLD_MESH(dcoor)
  implicit none
  real(8), intent(in) :: dcoor(3, *)
  integer :: i

  do i = 1, QuadSc%ndof
    myALE%OldCoor(:, i) = dcoor(:, i)
  end do

end subroutine STORE_OLD_MESH

!=========================================================================
! STORE_NEW_MESH - Store new mesh coordinates for ALE
! Saves updated mesh coordinates after deformation
!=========================================================================
subroutine STORE_NEW_MESH(dcoor)
  implicit none
  real(8), intent(in) :: dcoor(3, *)
  integer :: i

  do i = 1, QuadSc%ndof
    myALE%NewCoor(:, i) = dcoor(:, i)
  end do

end subroutine STORE_NEW_MESH

!=========================================================================
! GET_MESH_VELO - Compute mesh velocity from coordinate differences
! Calculates mesh velocity for ALE formulation: v_mesh = (x_new - x_old)/dt
!=========================================================================
subroutine GET_MESH_VELO()
  implicit none
  integer :: i
  real(8) :: dmax, daux

  dmax = 0.0d0
  do i = 1, QuadSc%ndof
    myALE%MeshVelo(:, i) = (myALE%NewCoor(:, i) - myALE%OldCoor(:, i)) / tstep

    ! Compute magnitude of mesh velocity
    daux = sqrt(myALE%MeshVelo(1, i)**2 + myALE%MeshVelo(2, i)**2 + &
                myALE%MeshVelo(3, i)**2)
    if (daux > dmax) dmax = daux
  end do

  call COMM_Maximum(dmax)
  if (myid == showid) write (*, *) 'Maximum mesh velocity: ', dmax

end subroutine GET_MESH_VELO

!=========================================================================
! RotateMyMesh - Rotate mesh for rotating geometry simulations
! Applies rotation transformation for time-dependent rotation (e.g., mixers)
!=========================================================================
subroutine RotateMyMesh(dcoor)
  implicit none
  real(8), intent(inout) :: dcoor(3, *)
  real(8) :: x, y, z, r, a, dAlpha
  real(8), parameter :: PI = 3.141592654d0
  integer :: i, nnn

  ! Compute rotation angle: 250 RPM
  dAlpha = 1.0d0 * (250.0d0 / 60.0d0) * 2.0d0 * timens * PI

  if (myid == 0) then
    nnn = KNVT(NLMAX)
  else
    nnn = KNVT(NLMAX + 1)
  end if

  ! Apply rotation transformation
  do i = 1, nnn
    x = myALE%OrigCoor(1, i)
    y = myALE%OrigCoor(2, i)
    r = sqrt(x * x + y * y)
    a = atan(y / x)
    if (x < 0.0d0) a = a + PI
    a = a + dAlpha
    dcoor(1, i) = r * cos(a)
    dcoor(2, i) = r * sin(a)
  end do

  ! Update coordinates on coarse level
  ILEV = NLMIN
  call SETLEV(2)
  call ExchangeNodeValuesOnCoarseLevel(DWORK(L(LCORVG)), KWORK(L(LVERT)), NVT, NEL)

  ILEV = NLMAX
  call SETLEV(2)

end subroutine RotateMyMesh

!=========================================================================
! StoreOrigCoor - Store original mesh coordinates
! Saves reference coordinates for rotation or deformation
!=========================================================================
subroutine StoreOrigCoor(dcoor)
  implicit none
  real(8), intent(in) :: dcoor(3, *)
  integer :: i, nnn

  if (myid == 0) then
    nnn = KNVT(NLMAX)
  else
    nnn = KNVT(NLMAX + 1)
  end if

  do i = 1, nnn
    myALE%OrigCoor(:, i) = dcoor(:, i)
  end do

end subroutine StoreOrigCoor

!=========================================================================
! GetMeshVelocity2 - Alternative mesh velocity computation
! Computes mesh velocity from Q2 coordinate differences
!=========================================================================
subroutine GetMeshVelocity2(mfile)
  implicit none
  integer, intent(in) :: mfile
  real(8) :: dMaxVelo, daux
  integer :: i

  dMaxVelo = 0.0d0
  if (myid /= 0) then
    do i = 1, QuadSc%ndof
      myALE%MeshVelo(:, i) = (myQ2Coor(:, i) - myALE%Q2Coor_old(:, i)) / tstep
      daux = myALE%MeshVelo(1, i)**2 + myALE%MeshVelo(2, i)**2 + myALE%MeshVelo(3, i)**2
      if (dMaxVelo < daux) dMaxVelo = daux
    end do
  end if

  call COMM_Maximum(dMaxVelo)

  if (myid == 1) then
    write (mfile, *) "Maximum Mesh Velocity: ", sqrt(dMaxVelo)
    write (mterm, *) "Maximum Mesh Velocity: ", sqrt(dMaxVelo)
  end if

end subroutine GetMeshVelocity2

!=========================================================================
! StaticMeshAdaptation - Apply static mesh adaptation
! Loads pre-computed adapted mesh from file
!=========================================================================
subroutine StaticMeshAdaptation()
  implicit none
  integer :: iAdaptMeshLevel, idL

  if (.not. bMeshAdaptation) return

  ! Read adaptation level from file
  open (474, file=adjustl(trim(cAdaptedMeshFile))//"/level.prf")
  read (474, *) iAdaptMeshLevel
  close (474)

  idL = iAdaptMeshLevel - NLMAX

  call CreateDumpStructures(idL)
  call LoadSmartAdaptedMeshFile(DWORK(L(KLCVG(1))), cAdaptedMeshFile, idL)

  ! Refresh coordinates on all levels
  do ILEV = iAdaptMeshLevel, NLMAX
    call SETLEV(2)
    write (*, *) 'Mesh levels:', ilev, iAdaptMeshLevel, NLMAX
    call RefreshCoordinates(DWORK(L(KLCVG(ILEV + 1))), DWORK(L(KLCAG(ILEV))), &
                            KWORK(L(KLVERT(ILEV))), KWORK(L(KLEDGE(ILEV))), KWORK(L(KLAREA(ILEV))))
  end do

end subroutine StaticMeshAdaptation

!=========================================================================
! CoorWriter - Write coordinates to file
! Utility for writing mesh coordinates for debugging/visualization
!=========================================================================
subroutine CoorWriter(dcorvg, nvt, cF)
  implicit none
  character(len=*), intent(in) :: cF
  real(8), intent(in) :: dcorvg(3, *)
  integer, intent(in) :: nvt
  integer :: i

  open (547, file=trim(adjustl(cF)))
  do i = 1, nvt
    write (547, *) dcorvg(:, i)
  end do
  close (547)

end subroutine CoorWriter

!=========================================================================
! RefreshCoordinates - Recompute edge, face, and element coordinates
! Updates coordinates of edge midpoints, face centers, and element centers
! from vertex coordinates for Q2 finite element discretization
!=========================================================================
subroutine RefreshCoordinates(dcorvg, dcorag, kvert, kedge, karea)
  implicit none
  real(8), intent(inout) :: dcorvg(3, *)
  real(8), intent(inout) :: dcorag(3, *)
  integer, intent(in) :: kvert(8, *), kedge(12, *), karea(6, *)

  real(8) :: px, py, pz, dist
  integer :: i, j, k, ivt1, ivt2, ivt3, ivt4

  ! Element connectivity arrays
  integer :: NeighE(2, 12), NeighA(4, 6)
  data NeighE/1, 2, 2, 3, 3, 4, 4, 1, 1, 5, 2, 6, 3, 7, 4, 8, 5, 6, 6, 7, 7, 8, 8, 5/
  data NeighA/1, 2, 3, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 4, 8, 7, 4, 1, 5, 8, 5, 6, 7, 8/

  ! Compute edge midpoint coordinates
  k = 1
  do i = 1, nel
    do j = 1, 12
      if (k == kedge(j, i)) then
        ivt1 = kvert(NeighE(1, j), i)
        ivt2 = kvert(NeighE(2, j), i)
        px = 0.5d0 * (dcorvg(1, ivt1) + dcorvg(1, ivt2))
        py = 0.5d0 * (dcorvg(2, ivt1) + dcorvg(2, ivt2))
        pz = 0.5d0 * (dcorvg(3, ivt1) + dcorvg(3, ivt2))
        dcorvg(1, nvt + k) = px
        dcorvg(2, nvt + k) = py
        dcorvg(3, nvt + k) = pz
        k = k + 1
      end if
    end do
  end do

  ! Compute face center coordinates
  k = 1
  do i = 1, nel
    do j = 1, 6
      if (k == karea(j, i)) then
        ivt1 = kvert(NeighA(1, j), i)
        ivt2 = kvert(NeighA(2, j), i)
        ivt3 = kvert(NeighA(3, j), i)
        ivt4 = kvert(NeighA(4, j), i)
        px = 0.25d0 * (dcorvg(1, ivt1) + dcorvg(1, ivt2) + dcorvg(1, ivt3) + dcorvg(1, ivt4))
        py = 0.25d0 * (dcorvg(2, ivt1) + dcorvg(2, ivt2) + dcorvg(2, ivt3) + dcorvg(2, ivt4))
        pz = 0.25d0 * (dcorvg(3, ivt1) + dcorvg(3, ivt2) + dcorvg(3, ivt3) + dcorvg(3, ivt4))
        dcorag(1, k) = px
        dcorag(2, k) = py
        dcorag(3, k) = pz
        dcorvg(1, nvt + net + k) = px
        dcorvg(2, nvt + net + k) = py
        dcorvg(3, nvt + net + k) = pz
        k = k + 1
      end if
    end do
  end do

  ! Compute element center coordinates
  do i = 1, nel
    px = 0.0d0
    py = 0.0d0
    pz = 0.0d0
    do j = 1, 8
      px = px + 0.125d0 * dcorvg(1, kvert(j, i))
      py = py + 0.125d0 * dcorvg(2, kvert(j, i))
      pz = pz + 0.125d0 * dcorvg(3, kvert(j, i))
    end do
    dcorvg(1, nvt + net + nat + i) = px
    dcorvg(2, nvt + net + nat + i) = py
    dcorvg(3, nvt + net + nat + i) = pz
  end do

end subroutine RefreshCoordinates
