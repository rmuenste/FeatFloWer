!=========================================================================
! QuadSc_material_properties.f90
!
! Material property and density update functions
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================

!=========================================================================
! GetNonNewtViscosity - Compute non-Newtonian viscosity from velocity field
! Uses L2 projection to compute shear rate and viscosity at Q2 nodes
!=========================================================================
subroutine GetNonNewtViscosity()
  implicit none
  integer :: i, ilevel
  real(8) :: ViscosityModel
  external E013

  ilev = mg_mesh%nlmax
  ilevel = mg_mesh%nlmax

  if (myid /= 0) then
    QuadSc%defU = 0.0d0
    QuadSc%defV = 0.0d0
    QuadSc%defW = 0.0d0
    call L2ProjVisco(QuadSc%ValU, QuadSc%ValV, QuadSc%ValW, &
                     Temperature, MaterialDistribution(ilev)%x, &
                     QuadSc%defU, QuadSc%defV, QuadSc%defW, &
                     mg_mesh%level(ilevel)%kvert, &
                     mg_mesh%level(ilevel)%karea, &
                     mg_mesh%level(ilevel)%kedge, &
                     mg_mesh%level(ilevel)%dcorvg, E013)

    call E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)

    ! Allocate shear rate array if needed
    if (.not. allocated(Shearrate)) allocate (Shearrate(QuadSc%ndof))

    ! Compute shear rate and viscosity at each node
    do i = 1, QuadSc%ndof
      Shearrate(i) = QuadSc%defV(i) / QuadSc%defW(i)
      Viscosity(i) = QuadSc%defU(i) / QuadSc%defW(i)
    end do

  end if

end subroutine GetNonNewtViscosity

!=========================================================================
! FilterColdElements - Mark cold elements as solidified
! Identifies and marks elements below crystallization temperature (80°C)
!=========================================================================
subroutine FilterColdElements(mfile)
  implicit none
  integer, intent(in) :: mfile
  integer :: i
  real(8) :: dj

  if (myid /= 0) then
    dj = 0.0d0
    do i = 1, QuadSc%ndof
      ! Mark elements below 80°C as solidified
      if (MixerKNPR(i) == 0 .and. Temperature(i) < 80.0d0) then
        MixerKNPR(i) = 105
        Shell(i) = -1.0d0
        dj = dj + 1.0d0
      end if
    end do
  end if

  call COMM_SUMM(dj)

  if (myid == 1) then
    write (mfile, *) 'Number of solidified dofs: ', nint(dj)
    write (mterm, *) 'Number of solidified dofs: ', nint(dj)
  end if

end subroutine FilterColdElements
!=========================================================================
! UpdateDensityDistribution_XSE - Update density based on temperature
! Computes element-wise density from temperature using material models:
! - DENSITY model: ρ = ρ₀ - T*α (thermal expansion)
! - SPECVOLUME model: ρ = 1/(v₀ + T*β) (specific volume)
!=========================================================================
subroutine UpdateDensityDistribution_XSE(mfile)
  implicit none
  integer, intent(in) :: mfile
  integer :: i, iel, ifld, iMat
  real(8) :: daux, taux, dAlpha, dMaxMat, dWSFactor
  real(8) :: AlphaViscosityMatModel, WallSlip

  if (.not. bMasterTurnedOn) return

  if (myid == 1) then
    write (MTERM, *) "Update of the density distribution!"
    write (MFILE, *) "Update of the density distribution!"
  end if

  do ILEV = NLMIN, NLMAX

    do iel = 1, mg_mesh%level(ilev)%nel

      i = mg_mesh%level(ilev)%nvt + &
          mg_mesh%level(ilev)%net + &
          mg_mesh%level(ilev)%nat + &
          iel

      taux = Temperature(i)

      ! Determine material type for multi-material simulations
      if (myMultiMat%nOfMaterials > 1) then
        iMat = myMultiMat%InitMaterial
        dMaxMat = 1.0d-5

        ! Find dominant material based on volume fraction
        do iFld = 2, GenLinScalar%nOfFields
          if (GenLinScalar%Fld(iFld)%val(i) > dMaxMat) then
            iMat = iFld - 1
            dMaxMat = GenLinScalar%Fld(iFld)%val(i)
          end if
        end do
      else
        iMat = 1
      end if

      ! Apply density model based on material thermodynamic properties
      if (adjustl(trim(myMultiMat%Mat(iMat)%Thermodyn%DensityModel)) == "DENSITY") then
        ! Linear thermal expansion: ρ = ρ₀ - T*α
        mgDensity(ILEV)%x(iel) = myMultiMat%Mat(iMat)%Thermodyn%densityT0 - &
                                 taux * myMultiMat%Mat(iMat)%Thermodyn%densitySteig
      end if

      if (adjustl(trim(myMultiMat%Mat(iMat)%Thermodyn%DensityModel)) == "SPECVOLUME") then
        ! Specific volume model: ρ = 1/(v₀ + T*β)
        mgDensity(ILEV)%x(iel) = 1.0d0 / (myMultiMat%Mat(iMat)%Thermodyn%densityT0 + &
                                 taux * myMultiMat%Mat(iMat)%Thermodyn%densitySteig)
      end if

    end do

  end do

  ! Send density data to master process for coarse-level solvers
  ILEV = LinSc%prm%MGprmIn%MedLev

  if (LinSc%prm%MGprmIn%MedLev >= 1 .and. LinSc%prm%MGprmIn%CrsSolverType <= 4) then
    call E010GATHR_L1(mgDensity(1)%x, mg_mesh%level(1)%nel)
  end if

  if (LinSc%prm%MGprmIn%MedLev >= 2 .and. LinSc%prm%MGprmIn%CrsSolverType <= 4) then
    call E010GATHR_L2(mgDensity(2)%x, mg_mesh%level(2)%nel)
  end if

  if (LinSc%prm%MGprmIn%MedLev >= 3 .and. LinSc%prm%MGprmIn%CrsSolverType <= 4) then
    call E010GATHR_L3(mgDensity(3)%x, mg_mesh%level(3)%nel)
  end if

  ILEV = NLMAX

end subroutine UpdateDensityDistribution_XSE

!=========================================================================
! GetAlphaNonNewtViscosity_sse - Compute non-Newtonian viscosity for SSE
! Calculates shear rate and viscosity with:
! - Temperature-dependent rheology
! - Multi-material support
! - Wall slip effects
!=========================================================================
subroutine GetAlphaNonNewtViscosity_sse()
  implicit none
  integer :: i, ifld, iMat
  real(8) :: daux, taux, dAlpha, dMaxMat, dWSFactor
  real(8) :: AlphaViscosityMatModel, WallSlip

  if (.not. bMasterTurnedOn) return

  ILEV = NLMAX
  call SETLEV(2)

  ! Compute velocity gradients for all three components
  call GetGradVelo_rhs(QuadSc, QuadSc%ValU)
  call E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)
  call GetGradVelo_val(QuadSc, 1)

  call GetGradVelo_rhs(QuadSc, QuadSc%ValV)
  call E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)
  call GetGradVelo_val(QuadSc, 2)

  call GetGradVelo_rhs(QuadSc, QuadSc%ValW)
  call E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)
  call GetGradVelo_val(QuadSc, 3)

  ! Compute shear rate and viscosity at each node
  do i = 1, size(QuadSc%ValU)

    ! Compute second invariant of strain rate tensor: II = tr(D²)
    daux = QuadSc%ValUx(i)**2.0d0 + QuadSc%ValVy(i)**2.0d0 + QuadSc%ValWz(i)**2.0d0 + &
           0.5d0 * (QuadSc%ValUy(i) + QuadSc%ValVx(i))**2.0d0 + &
           0.5d0 * (QuadSc%ValUz(i) + QuadSc%ValWx(i))**2.0d0 + &
           0.5d0 * (QuadSc%ValVz(i) + QuadSc%ValWy(i))**2.0d0

    ! Get temperature (sampled from Temperature field)
    taux = Temperature(i)

    ! Compute viscosity based on material type
    if (myMultiMat%nOfMaterials > 1) then
      ! Multi-material: determine dominant material
      iMat = myMultiMat%InitMaterial
      dMaxMat = 1.0d-5
      do iFld = 2, GenLinScalar%nOfFields
        if (GenLinScalar%Fld(iFld)%val(i) > dMaxMat) then
          iMat = iFld - 1
          dMaxMat = GenLinScalar%Fld(iFld)%val(i)
        end if
      end do

      Shearrate(i) = sqrt(2.0d0 * daux)
      Viscosity(i) = AlphaViscosityMatModel(daux, iMat, taux)

      ! Apply wall slip correction if enabled
      if (myMultiMat%Mat(iMat)%Rheology%bWallSlip) then
        dWSFactor = WallSlip(shell(i), screw(i), iMat, Viscosity(i) * Shearrate(i))
        Viscosity(i) = dWSFactor * Viscosity(i)
      end if
    else
      ! Single material
      Shearrate(i) = sqrt(2.0d0 * daux)
      Viscosity(i) = AlphaViscosityMatModel(daux, 1, taux)

      if (myMultiMat%Mat(1)%Rheology%bWallSlip) then
        dWSFactor = WallSlip(shell(i), screw(i), 1, Viscosity(i) * Shearrate(i))
        Viscosity(i) = dWSFactor * Viscosity(i)
      end if
    end if

  end do

end subroutine GetAlphaNonNewtViscosity_sse

!=========================================================================
! UpdateMaterialProperties - Update material distribution on all grid levels
! Propagates material information from finest to coarser multigrid levels
! using volume fraction averaging for multi-material flows
!=========================================================================
subroutine UpdateMaterialProperties
  implicit none
  integer :: ii, iMat, jel(8), ndof_nel, iMaxFrac
  real(8) :: dMaxFrac
  real(8), allocatable :: dFrac(:)

  if (myid == 1) write (*, *) 'Updating Material Properties Distribution... '

  if (myid /= 0) then

    if (allocated(MaterialDistribution)) then

      ! Initialize finest level from coarser level on restart
      if (istart == 2) then
        do iel = 1, mg_mesh%level(nlmax - 1)%nel
          ! Get 8 child elements on finer level
          jel(1) = iel
          jel(2) = mg_mesh%level(ilev + 1)%kadj(3, jel(1))
          jel(3) = mg_mesh%level(ilev + 1)%kadj(3, jel(2))
          jel(4) = mg_mesh%level(ilev + 1)%kadj(3, jel(3))
          jel(5) = mg_mesh%level(ilev + 1)%kadj(6, jel(1))
          jel(6) = mg_mesh%level(ilev + 1)%kadj(6, jel(2))
          jel(7) = mg_mesh%level(ilev + 1)%kadj(6, jel(3))
          jel(8) = mg_mesh%level(ilev + 1)%kadj(6, jel(4))

          ! Assign parent material to all children
          iMat = MaterialDistribution(nlmax - 1)%x(jel(1))
          do ii = 1, 8
            MaterialDistribution(nlmax)%x(jel(1)) = iMat
          end do
        end do
      end if

      ndof_nel = (knvt(NLMAX) + knat(NLMAX) + knet(NLMAX))
      allocate (dFrac(myMultiMat%nOfMaterials))

      ! Propagate material distribution from fine to coarse levels
      do ilev = NLMAX - 1, NLMIN, -1
        ndof_nel = (knvt(ilev + 1) + knat(ilev + 1) + knet(ilev + 1))

        do iel = 1, mg_mesh%level(ilev)%nel
          ! Get 8 child elements on finer level
          jel(1) = iel
          jel(2) = mg_mesh%level(ilev + 1)%kadj(3, jel(1))
          jel(3) = mg_mesh%level(ilev + 1)%kadj(3, jel(2))
          jel(4) = mg_mesh%level(ilev + 1)%kadj(3, jel(3))
          jel(5) = mg_mesh%level(ilev + 1)%kadj(6, jel(1))
          jel(6) = mg_mesh%level(ilev + 1)%kadj(6, jel(2))
          jel(7) = mg_mesh%level(ilev + 1)%kadj(6, jel(3))
          jel(8) = mg_mesh%level(ilev + 1)%kadj(6, jel(4))

          ! Compute volume fraction for each material
          dFrac = 0.0d0
          do ii = 1, 8
            iMat = MaterialDistribution(ilev + 1)%x(jel(ii))
            dFrac(iMat) = dFrac(iMat) + 1.0d0
          end do

          ! Assign dominant material to parent element
          dMaxFrac = 0.0d0
          do ii = 1, myMultiMat%nOfMaterials
            if (dFrac(ii) > dMaxFrac) then
              dMaxFrac = dFrac(ii)
              iMaxFrac = ii
            end if
          end do

          MaterialDistribution(ilev)%x(jel(1)) = iMaxFrac

        end do
      end do

      deallocate (dFrac)
    end if

  else
    ! Master process: initialize all levels with default material
    do ilev = NLMAX, NLMIN, -1
      MaterialDistribution(ilev)%x = myMultiMat%InitMaterial
    end do
  end if

end subroutine UpdateMaterialProperties
