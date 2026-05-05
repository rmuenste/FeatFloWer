module attraction_smoother
!-------------------------------------------------------------------------------------------------
! Shared node-attraction utilities for mesh smoothing.
!-------------------------------------------------------------------------------------------------
  use var_QuadScalar, only : tMultiMesh

  implicit none

  private

  type, public :: tAttractionConfig
    logical :: enabled = .false.
    logical :: umbrellaWeightEnabled = .false.
    character(len=32) :: mode = "node_attraction"
    character(len=32) :: distanceProvider = "cylinder"
    integer :: nAttractSteps = 8
    integer :: nPostUmbrellaSteps = 2
    integer :: applyOnLevel = 0
    real*8 :: omega = 0.2d0
    real*8 :: bandOuterStrong = 0.04d0
    real*8 :: bandOuterMax = 0.15d0
    real*8 :: bandInner = 0.04d0
    real*8 :: innerMaxAlpha = 0.45d0
    real*8 :: outerMinAlpha = 0.6d0
    real*8 :: outerMaxAlpha = 1.0d0
    real*8 :: outerFarExponent = 2.0d0
    real*8 :: umbrellaWeightCap = 25d0
    real*8 :: umbrellaWeightPower = 2.3d0
    real*8 :: cylinderCenter(3) = (/0.5d0, 0.2d0, 0.205d0/)
    real*8 :: cylinderRadius = 0.05d0
    real*8 :: cylinderLength = 0.4105d0
    logical :: cylinderUseCaps = .false.
  end type tAttractionConfig

  public :: AttractionSmoother_GetConfig
  public :: AttractionSmoother_Apply
  public :: AttractionSmoother_ComputeDistance
  public :: AttractionSmoother_ComputeWeight
  public :: AttractionSmoother_UseUmbrellaDistanceWeight
  public :: AttractionSmoother_GetUmbrellaWeight

contains

  subroutine AttractionSmoother_GetConfig(cfg)
    use var_QuadScalar, only : bAttractionSmootherEnable, bUmbrellaDistanceWeightEnable, &
      cAttractionMode, cAttractionDistanceProvider, nAttractionSteps, &
      nPostAttractionUmbrellaSteps, iAttractionApplyOnLevel, dAttractionOmega, &
      dAttractionBandOuterStrong, dAttractionBandOuterMax, dAttractionBandInner, &
      dAttractionInnerMaxAlpha, dAttractionOuterMinAlpha, dAttractionOuterMaxAlpha, &
      dAttractionOuterFarExponent, dUmbrellaWeightCap, dUmbrellaWeightPower, &
      bAttractionCylinderUseCaps, bFAC3D_CylUmbrellaWeight, dFAC3D_CylCenter, &
      dFAC3D_CylRadius, dFAC3D_CylLength
    implicit none
    type(tAttractionConfig), intent(out) :: cfg

    cfg%enabled = bAttractionSmootherEnable .OR. bFAC3D_CylUmbrellaWeight
    cfg%umbrellaWeightEnabled = bUmbrellaDistanceWeightEnable .OR. bFAC3D_CylUmbrellaWeight
    cfg%mode = cAttractionMode
    cfg%distanceProvider = cAttractionDistanceProvider
    cfg%nAttractSteps = nAttractionSteps
    cfg%nPostUmbrellaSteps = nPostAttractionUmbrellaSteps
    cfg%applyOnLevel = iAttractionApplyOnLevel
    cfg%omega = dAttractionOmega
    cfg%bandOuterStrong = dAttractionBandOuterStrong
    cfg%bandOuterMax = dAttractionBandOuterMax
    cfg%bandInner = dAttractionBandInner
    cfg%innerMaxAlpha = dAttractionInnerMaxAlpha
    cfg%outerMinAlpha = dAttractionOuterMinAlpha
    cfg%outerMaxAlpha = dAttractionOuterMaxAlpha
    cfg%outerFarExponent = dAttractionOuterFarExponent
    cfg%umbrellaWeightCap = dUmbrellaWeightCap
    cfg%umbrellaWeightPower = dUmbrellaWeightPower
    cfg%cylinderCenter = dFAC3D_CylCenter
    cfg%cylinderRadius = dFAC3D_CylRadius
    cfg%cylinderLength = dFAC3D_CylLength
    cfg%cylinderUseCaps = bAttractionCylinderUseCaps

  end subroutine AttractionSmoother_GetConfig

  logical function AttractionSmoother_UseUmbrellaDistanceWeight()
    implicit none
    type(tAttractionConfig) :: cfg

    call AttractionSmoother_GetConfig(cfg)
    AttractionSmoother_UseUmbrellaDistanceWeight = cfg%umbrellaWeightEnabled

  end function AttractionSmoother_UseUmbrellaDistanceWeight

  subroutine AttractionSmoother_Apply(mgMesh, ilev, cfg)
    use Parametrization, only : ParametrizeBndryPoints_STRCT
    implicit none
    type(tMultiMesh), intent(inout) :: mgMesh
    integer, intent(in) :: ilev
    type(tAttractionConfig), intent(in) :: cfg

    integer :: i, iStep, nVtx
    real*8, allocatable :: dist(:)
    real*8, allocatable :: witness(:,:)
    logical, allocatable :: attractable(:)
    real*8 :: alpha

    if (.not. cfg%enabled) return
    if (cfg%nAttractSteps .le. 0) return
    if (trim(adjustl(cfg%distanceProvider)) .ne. "cylinder") return

    nVtx = mgMesh%level(ilev)%nvt
    allocate(dist(nVtx))
    allocate(witness(3,nVtx))
    allocate(attractable(nVtx))

    do iStep = 1, cfg%nAttractSteps
      call AttractionSmoother_ComputeDistance(nVtx, mgMesh%level(ilev)%dcorvg, &
        dist, witness, attractable, cfg)

      do i = 1, nVtx
        if (.not. attractable(i)) cycle
        call compute_alpha(dist(i), cfg, alpha)
        if (alpha .eq. 0d0) cycle

        mgMesh%level(ilev)%dcorvg(1,i) = mgMesh%level(ilev)%dcorvg(1,i) + &
          cfg%omega*alpha*(witness(1,i) - mgMesh%level(ilev)%dcorvg(1,i))
        mgMesh%level(ilev)%dcorvg(2,i) = mgMesh%level(ilev)%dcorvg(2,i) + &
          cfg%omega*alpha*(witness(2,i) - mgMesh%level(ilev)%dcorvg(2,i))
        mgMesh%level(ilev)%dcorvg(3,i) = mgMesh%level(ilev)%dcorvg(3,i) + &
          cfg%omega*alpha*(witness(3,i) - mgMesh%level(ilev)%dcorvg(3,i))
      end do

      call ParametrizeBndryPoints_STRCT(mgMesh, ilev)
    end do

    deallocate(dist, witness, attractable)

  end subroutine AttractionSmoother_Apply

  subroutine AttractionSmoother_ComputeDistance(nVtx, dcorvg, dist, witness, attractable, cfgIn)
    implicit none
    integer, intent(in) :: nVtx
    real*8, dimension(:,:), intent(in) :: dcorvg
    real*8, dimension(:), intent(out) :: dist
    real*8, dimension(:,:), intent(out), optional :: witness
    logical, dimension(:), intent(out), optional :: attractable
    type(tAttractionConfig), intent(in), optional :: cfgIn

    type(tAttractionConfig) :: cfg

    if (present(cfgIn)) then
      cfg = cfgIn
    else
      call AttractionSmoother_GetConfig(cfg)
    end if

    select case (trim(adjustl(cfg%distanceProvider)))
    case ("cylinder")
      call compute_cylinder_witness(nVtx, dcorvg, dist, witness, attractable, cfg)
    case default
      dist(1:nVtx) = 0d0
      if (present(witness)) witness(:,1:nVtx) = dcorvg(:,1:nVtx)
      if (present(attractable)) attractable(1:nVtx) = .false.
    end select

  end subroutine AttractionSmoother_ComputeDistance

  subroutine AttractionSmoother_ComputeWeight(nVtx, dist, weight, cfgIn)
    implicit none
    integer, intent(in) :: nVtx
    real*8, intent(in) :: dist(*)
    real*8, intent(out) :: weight(*)
    type(tAttractionConfig), intent(in), optional :: cfgIn

    integer :: i
    real*8 :: alpha
    type(tAttractionConfig) :: cfg

    if (present(cfgIn)) then
      cfg = cfgIn
    else
      call AttractionSmoother_GetConfig(cfg)
    end if

    do i = 1, nVtx
      call compute_alpha(dist(i), cfg, alpha)
      if (dist(i) .lt. 0d0) then
        weight(i) = -alpha
      else
        weight(i) = alpha
      end if
    end do

  end subroutine AttractionSmoother_ComputeWeight

  subroutine AttractionSmoother_GetUmbrellaWeight(dist, weight, cfgIn)
    implicit none
    real*8, intent(in) :: dist
    real*8, intent(out) :: weight
    type(tAttractionConfig), intent(in), optional :: cfgIn

    type(tAttractionConfig) :: cfg
    real*8 :: kernel

    if (present(cfgIn)) then
      cfg = cfgIn
    else
      call AttractionSmoother_GetConfig(cfg)
    end if

    call KernelFunction(dist, kernel)
    weight = min(kernel, cfg%umbrellaWeightCap)
    weight = weight**cfg%umbrellaWeightPower

  end subroutine AttractionSmoother_GetUmbrellaWeight

  subroutine compute_cylinder_witness(nVtx, dcorvg, dist, witness, attractable, cfg)
    implicit none
    integer, intent(in) :: nVtx
    real*8, dimension(:,:), intent(in) :: dcorvg
    real*8, dimension(:), intent(out) :: dist
    real*8, dimension(:,:), intent(out), optional :: witness
    logical, dimension(:), intent(out), optional :: attractable
    type(tAttractionConfig), intent(in) :: cfg

    integer :: i
    real*8 :: dx, dy, dz, dzAbs, rxy, zGap, radialGap, halfLen
    real*8 :: nearX, nearY, nearZ, invR

    halfLen = 0.5d0*cfg%cylinderLength

    do i = 1, nVtx
      dx = dcorvg(1,i) - cfg%cylinderCenter(1)
      dy = dcorvg(2,i) - cfg%cylinderCenter(2)
      dz = dcorvg(3,i) - cfg%cylinderCenter(3)
      dzAbs = abs(dz)
      rxy = sqrt(dx*dx + dy*dy)
      zGap = dzAbs - halfLen
      radialGap = rxy - cfg%cylinderRadius
      dist(i) = max(radialGap, zGap)

      if (rxy .le. 1d-14) then
        if (present(witness)) witness(:,i) = dcorvg(:,i)
        if (present(attractable)) attractable(i) = .false.
        cycle
      end if

      invR = 1d0/rxy
      nearX = cfg%cylinderCenter(1) + cfg%cylinderRadius*dx*invR
      nearY = cfg%cylinderCenter(2) + cfg%cylinderRadius*dy*invR
      nearZ = dcorvg(3,i)

      if (zGap .gt. radialGap) then
        if (.not. cfg%cylinderUseCaps) then
          if (present(attractable)) attractable(i) = .false.
          if (present(witness)) witness(:,i) = (/nearX, nearY, nearZ/)
          cycle
        end if
        nearZ = cfg%cylinderCenter(3) + sign(halfLen, dz)
      end if

      if (present(witness)) witness(:,i) = (/nearX, nearY, nearZ/)
      if (present(attractable)) attractable(i) = .true.
    end do

  end subroutine compute_cylinder_witness

  subroutine compute_alpha(dist, cfg, alpha)
    implicit none
    real*8, intent(in) :: dist
    type(tAttractionConfig), intent(in) :: cfg
    real*8, intent(out) :: alpha

    real*8 :: absDist, t

    alpha = 0d0

    if (dist .gt. 0d0 .and. dist .le. cfg%bandOuterMax) then
      if (dist .lt. cfg%bandOuterStrong) then
        t = dist/cfg%bandOuterStrong
        alpha = cfg%outerMinAlpha + (cfg%outerMaxAlpha - cfg%outerMinAlpha)*t
      else
        t = (cfg%bandOuterMax - dist)/(cfg%bandOuterMax - cfg%bandOuterStrong)
        alpha = t**cfg%outerFarExponent
      end if
    else if (dist .lt. 0d0) then
      absDist = abs(dist)
      if (absDist .le. cfg%bandInner) then
        alpha = cfg%innerMaxAlpha*(1d0 - absDist/cfg%bandInner)
      end if
    end if

  end subroutine compute_alpha

  subroutine KernelFunction(dIn, f)
    implicit none
    real*8, intent(in) :: dIn
    real*8, intent(out) :: f
    real*8 :: d, daux

    d = dIn
    if (d .lt. 0d0) d = 2.5d0*abs(d)

    if (d .lt. 3.0d0) then
      daux = 2d0 + 1.0d0*(3.0d0-d)
    else
      daux = 2d0 - 0.2d0*(d-3.0d0)
    end if

    f = max(daux, 0.8d0)

  end subroutine KernelFunction

end module attraction_smoother
