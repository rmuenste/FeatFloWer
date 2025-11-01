!=========================================================================
! QuadSc_corrections.f90
!
! Solution correction and RHS assembly functions for QuadScalar module
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================

!=========================================================================
! Velocity_Correction - Apply velocity correction step
! Updates velocity field: U = U~ - k M^-1 B P
!=========================================================================
subroutine Velocity_Correction()
  implicit none
  integer :: i

  ! Update of U = U~ - k M^-1 B P
  ILEV = NLMAX
  qlMat => mg_qlMat(ILEV)
  BXMat => mg_BXMat(ILEV)%a
  BYMat => mg_BYMat(ILEV)%a
  BZMat => mg_BZMat(ILEV)%a
  MlRhoPmat => mg_MlRhoPmat(ILEV)%a

  call B_Mul_U(qlMat%ColA, qlMAt%LdA, BXMat, BYMat, BZMat, LinSc%ValP(NLMAX)%x, &
               QuadSc%defU, QuadSc%defV, QuadSc%defW, QuadSc%ndof, -TSTEP, 0.0d0)

  call E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)

  call Boundary_QuadScalar_Def()

  do i = 1, QuadSc%ndof
    QuadSc%valU(i) = QuadSc%valU(i) - QuadSc%defU(i) / MlRhoPmat(i)
    QuadSc%valV(i) = QuadSc%valV(i) - QuadSc%defV(i) / MlRhoPmat(i)
    QuadSc%valW(i) = QuadSc%valW(i) - QuadSc%defW(i) / MlRhoPmat(i)
  end do

end subroutine Velocity_Correction

!=========================================================================
! Pressure_Correction - Apply pressure correction step
! Updates pressure field with diffusion and extrapolation
!=========================================================================
subroutine Pressure_Correction()
  implicit none
  integer :: i
  real(8) :: dR, daux

  if (GAMMA > 0.0d0) then
    call SETLEV(2)
    P1iMMat => mg_P1iMMat(NLMAX)%a

    daux = THSTEP * (GAMMA + Properties%Viscosity(1))
    call AddDiffPrec(P1iMMat, LinSc%rhsP(NLMAX)%x, LinSc%ValP(NLMAX)%x, knel(nlmax), daux)
  end if

  do i = 1, lMat%nu
    LinSc%valP(NLMAX)%x(i) = LinSc%valP(NLMAX)%x(i) + LinSc%valP_old(i)
    LinSc%P_new(i) = 1.5d0 * LinSc%valP(NLMAX)%x(i) - 0.5d0 * LinSc%valP_old(i)
  end do

end subroutine Pressure_Correction

!=========================================================================
! AddPressureGradient - Add pressure gradient to velocity defect
!=========================================================================
subroutine AddPressureGradient()
  implicit none
  integer :: i, j, iel
  real(8) :: ddx, ddy, ddz, ddp

  ILEV = NLMAX
  qlMat => mg_qlMat(ILEV)
  BXMat => mg_BXMat(ILEV)%a
  BYMat => mg_BYMat(ILEV)%a
  BZMat => mg_BZMat(ILEV)%a

  call B_Mul_U(qlMat%ColA, qlMAt%LdA, BXMat, BYMat, BZMat, LinSc%valP(NLMAX)%x, &
               QuadSc%defU, QuadSc%defV, QuadSc%defW, QuadSc%ndof, TSTEP, 1.0d0)

end subroutine AddPressureGradient

!=========================================================================
! AddPressureGradientWithJump - Add pressure gradient with periodic jump
! Handles pressure discontinuities at periodic boundaries
!=========================================================================
subroutine AddPressureGradientWithJump()
  implicit none
  integer :: i, j, iel, ivt
  real(8) :: ddx, ddy, ddz, ddp, P(3)

  if (.not. allocated(dPeriodicVector)) allocate (dPeriodicVector(QuadSc%ndof))
  dPeriodicVector = 0.0d0

  ! Set up auxiliary pressure field
  do i = 1, size(LinSc%AuxP(NLMAX)%x)
    if (mod(i, 4) == 1) then
      LinSc%AuxP(NLMAX)%x(i) = 1.0d1
    else
      LinSc%AuxP(NLMAX)%x(i) = 0.0d0
    end if
  end do

  call B_Mul_U(qlMat%ColA, qlMAt%LdA, BXMat, BYMat, BZMat, LinSc%AuxP(NLMAX)%x, &
               QuadSc%AuxU, QuadSc%AuxV, QuadSc%AuxW, QuadSc%ndof, +1.0d0, 0.0d0)

  ! Identify periodic boundary points
  do i = 1, QuadSc%ndof
    if ((abs(myQ2Coor(1, i)) < +1.0d-1) .and. (myQ2Coor(2, i) < +0.0d0) .and. myid == 8) then
      dPeriodicVector(i) = QuadSc%auxU(i)
    else
      dPeriodicVector(i) = 0.0d0
    end if
  end do

  call B_Mul_U(qlMat%ColA, qlMAt%LdA, BXMat, BYMat, BZMat, LinSc%valP(NLMAX)%x, &
               QuadSc%defU, QuadSc%defV, QuadSc%defW, QuadSc%ndof, TSTEP, 1.0d0)

  do i = 1, QuadSc%ndof
    QuadSc%defU(i) = QuadSc%defU(i) + TSTEP * dPeriodicVector(i)
  end do

end subroutine AddPressureGradientWithJump

!=========================================================================
! AddPeriodicPressureGradient - Add periodic pressure gradient
! Adds pressure gradient contribution for periodic boundary conditions
!=========================================================================
subroutine AddPeriodicPressureGradient()
  implicit none
  integer :: i, j, iel
  real(8) :: ddx, ddy, ddz, ddp

  call B_Mul_U(qlMat%ColA, qlMAt%LdA, BXMat, BYMat, BZMat, LinSc%valP(NLMAX)%x, &
               QuadSc%defU, QuadSc%defV, QuadSc%defW, QuadSc%ndof, TSTEP, 1.0d0)

  do i = 1, QuadSc%ndof
    QuadSc%defW(i) = QuadSc%defW(i) + TSTEP * dPeriodicVector(i)
  end do

end subroutine AddPeriodicPressureGradient

!=========================================================================
! AddGravForce - Add gravitational body force
! Assembles gravitational force contribution to RHS
!=========================================================================
subroutine AddGravForce
  implicit none
  integer :: i, lint
  real(8) :: daux
  external E013

  ILEV = NLMAX
  call SETLEV(2)

  ! Compute magnitude of gravity vector
  daux = sqrt(Properties%Gravity(1)**2 + Properties%Gravity(2)**2 + Properties%Gravity(3)**2)

  if (daux > 0.0d0) then
    call Grav_QuadSc(QuadSc%defU, QuadSc%defV, QuadSc%defW, mgDensity(ILEV)%x, &
                     Properties%Gravity, QuadSc%ndof, &
                     mg_mesh%level(ilev)%kvert, &
                     mg_mesh%level(ilev)%karea, &
                     mg_mesh%level(ilev)%kedge, &
                     KWORK(L(KLINT(NLMAX))), &
                     mg_mesh%level(ilev)%dcorvg, &
                     tstep, E013)
  end if

end subroutine AddGravForce

!=========================================================================
! OperatorRegenaration - Regenerate system matrices
! Rebuilds various system matrices based on type parameter
!=========================================================================
subroutine OperatorRegenaration(iType)
  implicit none
  integer, intent(in) :: iType
  logical :: bHit

  bHit = .false.

  ! Diffusion matrix
  if (iType == myMatrixRenewal%D) then
    call Create_DiffMat(QuadSc)
    bHit = .true.
  end if

  ! Convection matrix
  if (iType == myMatrixRenewal%K) then
    call Create_KMat(QuadSc)
    bHit = .true.
  end if

  ! Mass matrix
  if (iType == myMatrixRenewal%M) then
    call Create_MRhoMat()
    bHit = .true.
  end if

  ! Stabilization matrix
  if (iType == myMatrixRenewal%S) then
    call Create_SMat(QuadSc)
    bHit = .true.
  end if

  ! Pressure coupling matrix
  if (iType == myMatrixRenewal%C) then
    call Create_BMat()

    if (myid /= master) then
      call Fill_QuadLinParMat()
    end if

    call SetSlipOnBandBT()

    call Create_CMat(QuadSc%knprU, QuadSc%knprV, QuadSc%knprW, LinSc%knprP, &
                     LinSc%prm%MGprmIn%MinLev, LinSc%prm%MGprmIn%CrsSolverType)

    if (myid /= master) then
      call Create_ParCMat(QuadSc%knprU, QuadSc%knprV, QuadSc%knprW)
    end if
    bHit = .true.
  end if

  ! Navier-Stokes stabilization
  if (iType == 1 .and. bNS_Stabilization) then
    call Create_hDiffMat()
    bHit = .true.
  end if

  ! Newton linearization for Burgers
  if (iType == 3 .and. (NewtonForBurgers /= 0.0d0)) then
    call Create_barMMat_iso(QuadSc)
    bHit = .true.
  end if

  if (myid == ShowID .and. bHit) write (MTERM, '(A)', advance='yes') " "

end subroutine OperatorRegenaration

!=========================================================================
! AddViscoStress - Add viscoelastic stress contribution
! Assembles viscoelastic stress tensor contributions to RHS
!=========================================================================
subroutine AddViscoStress()
  implicit none
  external E013

  ILEV = NLMAX
  call SETLEV(2)

  call AssembleViscoStress(QuadSc%defU, QuadSc%defV, QuadSc%defW, BndrForce, &
                           ViscoSc%Val11, ViscoSc%Val22, ViscoSc%Val33, &
                           ViscoSc%Val12, ViscoSc%Val13, ViscoSc%Val23, &
                           mg_mesh%level(ilev)%kvert, &
                           mg_mesh%level(ilev)%karea, &
                           mg_mesh%level(ilev)%kedge, &
                           mg_mesh%level(ilev)%dcorvg, &
                           Properties%Viscosity(2), tstep, Properties%ViscoLambda, &
                           ViscoElasticForce, E013)

  call Comm_SummN(ViscoElasticForce, 3)

end subroutine AddViscoStress
