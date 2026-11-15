!=========================================================================
! QuadSc_el_mm.f90
!
! Euler-Lagrange momentum models: fluid-side adapter. Included in Transport_Q2P1 (host
! association to mgDensity, mg_mesh, Properties, Create_MRhoMat, myMatrixRenewal, ...).
!
! The model policy lives in EL_MOMENTUM_MODEL; this file does the mesh, MPI and operator work
! the policy asks for. It is entered from Transport_q2p1_UxyzP_el only, and only when
! SimPar@ELMomentumModel selects a model: every other application and every 'none' run never
! reaches it, and no shared assembly routine has a branch for it.
!
! eps-weighted inertia (all models). The solver takes its inertia coefficient from the
! element-wise field mgDensity: lumped mass (Create_MRhoMat), convection (Create_KMat/CONVQ2),
! body forces (Grav_QuadSc) and, through MlRhoPmat, the velocity correction of the projection.
! EL_MM_Update_Fluid_Coefficients writes mgDensity = rho*eps_f on the finest level, restricts
! it volume-weighted to the coarser levels, hands the coarsest level to the master and
! rebuilds the mass matrices. With the porosity of the new particle positions (the pass after
! the particle advance) the momentum equation of the coming solve is
!     eps^{n+1} rho (u^{n+1} - u^n)/dt + eps^{n+1} rho (u . grad) u = ...
! i.e. the non-conservative form eps rho Du/Dt of cfdemSolverPiso. CFDEM's backward-Euler
! assembly ddt(eps,U) - Sp(ddt(eps),U) evaluates to eps^n (U^{n+1} - U^n)/dt, the same term
! with the porosity one step older; the difference is O(dt d eps/dt du/dt) and deliberate
! (one coefficient field for time derivative and convection).
!
! Not rebuilt: the pressure matrix C0 = B^T (rho Ml)^-1 B of start-up (MatrixRenewal C1). It
! is the preconditioner of the conservative continuity iteration
! (EL_MOMENTUM_MODEL:EL_MM_RESIDUAL_DIVISOR); the coefficients stored in mg_CMat do not depend
! on later mass rebuilds.
!
! Lifecycle of the mass arrays under a model: mg_MlRhoMat / mg_MlRhoPmat hold the weighted
! lumped mass from the first fluid solve on (also for the P1->Q2 pressure output projection
! and the histogram/profile post-processing, whose nodal weights are then eps-weighted);
! mg_MMat is overwritten by every rebuild and must not be used as a geometric-mass cache.
!
! Collective: the coarse-level hand-over is joined by every rank, so the caller guards on the
! deck key only (EL_MM_WEIGHTED_INERTIA).
!=========================================================================
SUBROUTINE EL_MM_Update_Fluid_Coefficients()

  USE EL_FIELDS, ONLY: el_field_data
  USE PP3D_MPI, ONLY: myid, master

  INTEGER :: ilev_save, jlev, iel, j, nel_c, child
  REAL*8 :: rho, vol, vsum, wsum
  LOGICAL, SAVE :: checked = .FALSE.

  IF (.NOT.checked) THEN
    CALL EL_MM_Check_Solver_Settings()
    checked = .TRUE.
  END IF

  rho = Properties%Density(1)
  ilev_save = ILEV

  IF (myid.NE.master) THEN
    IF (.NOT.ALLOCATED(el_field_data%epsilon_f)) THEN
      WRITE(*,*) 'EL_MM_Update_Fluid_Coefficients: porosity field not allocated on rank ', myid
      STOP 1
    END IF
    DO iel = 1, mg_mesh%level(NLMAX)%nel
      mgDensity(NLMAX)%x(iel) = rho*el_field_data%epsilon_f(iel)
    END DO
    ! volume-weighted restriction: a coarse element and its seven further children
    ! (refinement numbering: child 1 keeps the parent index, children 2..8 follow at
    ! nel_c + 7*(iel-1) + 1..7)
    DO jlev = NLMAX-1, NLMIN, -1
      nel_c = mg_mesh%level(jlev)%nel
      DO iel = 1, nel_c
        vol = mg_mesh%level(jlev+1)%dvol(iel)
        vsum = vol
        wsum = vol*mgDensity(jlev+1)%x(iel)
        DO j = 1, 7
          child = nel_c + 7*(iel-1) + j
          vol = mg_mesh%level(jlev+1)%dvol(child)
          vsum = vsum + vol
          wsum = wsum + vol*mgDensity(jlev+1)%x(child)
        END DO
        mgDensity(jlev)%x(iel) = wsum/vsum
      END DO
    END DO
  END IF

  ! coarsest level to the master (workers send, the master receives; every rank takes part)
  CALL E010_CollectCoarseVector(mgDensity(NLMIN)%x, mg_mesh%level(NLMIN)%nel)

  ! lumped mass and its parallel sum on all levels from the new coefficient
  CALL Create_MRhoMat()

  ILEV = ilev_save
  CALL SETLEV(2)

END SUBROUTINE EL_MM_Update_Fluid_Coefficients

!=========================================================================
! Solver settings a momentum model relies on, checked once at the first call (same deck
! values on every rank, hence a coordinated stop):
!   MatrixRenewal  M >= 1 (the mass arrays exist), K = 2 or 3 (convection rebuilt every step
!                  from the weighted coefficient), C = 1 (the start-up pressure matrix is the
!                  preconditioner; a per-step rebuild would take the weighted mass),
!   GAMMA = 0      (the rotational pressure update assumes the unweighted divergence source),
!   convective form of the convection term (the divergence form lacks the compensating
!                  -u div(eps u) of the non-conservative weighted inertia),
!   Newtonian fluid with one material.
!=========================================================================
SUBROUTINE EL_MM_Check_Solver_Settings()

  USE var_QuadScalar, ONLY: iConvectionForm, GAMMA, bNonNewtonian

  LOGICAL :: ok

  ok = .TRUE.
  IF (myMatrixRenewal%M.LT.1 .OR. (myMatrixRenewal%K.NE.2 .AND. myMatrixRenewal%K.NE.3) .OR. &
      myMatrixRenewal%C.NE.1) THEN
    WRITE(*,'(A,3(A,I0))') 'ELMomentumModel needs SimPar@MatrixRenewal with M >= 1, K = 2 or 3, C = 1;', &
      ' got M', myMatrixRenewal%M, ' K', myMatrixRenewal%K, ' C', myMatrixRenewal%C
    ok = .FALSE.
  END IF
  IF (GAMMA.NE.0.0d0) THEN
    WRITE(*,'(A)') 'ELMomentumModel needs GAMMA = 0.'
    ok = .FALSE.
  END IF
  IF (iConvectionForm.NE.0) THEN
    WRITE(*,'(A)') 'ELMomentumModel needs ELConvectionForm = convective.'
    ok = .FALSE.
  END IF
  IF (bNonNewtonian) THEN
    WRITE(*,'(A)') 'ELMomentumModel is implemented for the Newtonian fluid path only.'
    ok = .FALSE.
  END IF
  IF (.NOT.ok) STOP 1

END SUBROUTINE EL_MM_Check_Solver_Settings

!=========================================================================
! Dirichlet-boundary mask of the velocity dofs for the particle viscous force (EL_VISCOUS):
! a dof is marked when all three velocity components carry a Dirichlet condition (no-slip
! walls, inlets). Built once from the boundary flags of the finest level; workers only.
!=========================================================================
SUBROUTINE EL_MM_Prepare_Wall_Mask()

  USE EL_FIELDS, ONLY: el_field_data
  USE PP3D_MPI, ONLY: myid, master

  INTEGER :: i, ndof

  IF (myid.EQ.master) RETURN
  IF (ALLOCATED(el_field_data%wall_dof)) RETURN
  ndof = QuadSc%ndof
  ALLOCATE(el_field_data%wall_dof(ndof))
  DO i = 1, ndof
    el_field_data%wall_dof(i) = (QuadSc%knprU(NLMAX)%x(i).NE.0 .AND. &
                                 QuadSc%knprV(NLMAX)%x(i).NE.0 .AND. &
                                 QuadSc%knprW(NLMAX)%x(i).NE.0)
  END DO

END SUBROUTINE EL_MM_Prepare_Wall_Mask

!=========================================================================
! Model A, operator form: solid-fraction corrections of the unweighted operators
! (EL_MM_OPERATOR). Both routines are called from the fluid core behind the guard
! ".NOT.enable_fbm .AND. EL_MM_OPERATOR_FORM()" (deck-only), next to the unweighted operation
! they complete; the shared routines AddPressureGradient and Velocity_Correction are unchanged.
!
! EL_MM_Add_Operator_Terms - momentum right-hand side, after AddPressureGradient:
!     def += dt B p              (AddPressureGradient: int p div N, the unweighted -grad p)
!     def -= dt fp(p)            -> int p div(eps_h N), i.e. -eps grad p
!     def += dt c fv(u^n)        the solid-fraction part of the viscous term, explicit: together
!                                with the implicit unweighted diffusion matrix the fluid sees
!                                div(eps tau). The explicit share phi < 1 of an implicit term
!                                is unconditionally stable ((1 + phi a)/(1 + a) < 1).
! c is the coefficient of the fluid's diffusion operator (EL_FORCES:
! EL_FLUID_OPERATOR_VISCOSITY), the one the particle viscous force uses.
!=========================================================================
SUBROUTINE EL_MM_Add_Operator_Terms()

  USE EL_FIELDS, ONLY: el_field_data
  USE EL_MM_OPERATOR, ONLY: EL_MM_PHI_OPERATORS
  USE EL_FORCES, ONLY: EL_FLUID_OPERATOR_VISCOSITY
  USE PP3D_MPI, ONLY: myid, master

  REAL*8, ALLOCATABLE :: fp(:,:), fv(:,:)
  REAL*8 :: c
  INTEGER :: ndof, i

  IF (myid.EQ.master) RETURN
  ndof = QuadSc%ndof
  ALLOCATE(fp(3,ndof), fv(3,ndof))
  CALL EL_MM_PHI_OPERATORS(mg_mesh, NLMAX, mg_mesh%level(NLMAX)%dcorvg, &
    el_field_data%epsilon_transport_q2, LinSc%valP(NLMAX)%x, QuadSc%valU, QuadSc%valV, &
    QuadSc%valW, .TRUE., fp, fv)
  c = EL_FLUID_OPERATOR_VISCOSITY(Properties%Density(1), Properties%Viscosity(1))
  DO i = 1, ndof
    QuadSc%defU(i) = QuadSc%defU(i) + TSTEP*(c*fv(1,i) - fp(1,i))
    QuadSc%defV(i) = QuadSc%defV(i) + TSTEP*(c*fv(2,i) - fp(2,i))
    QuadSc%defW(i) = QuadSc%defW(i) + TSTEP*(c*fv(3,i) - fp(3,i))
  END DO
  DEALLOCATE(fp, fv)

END SUBROUTINE EL_MM_Add_Operator_Terms

!=========================================================================
! EL_MM_Correct_Velocity_Operator - projection, directly after Velocity_Correction and before
! Pressure_Correction (LinSc%valP(NLMAX) still holds the pressure increment q):
!     Velocity_Correction:   u += dt M^-1 B q        (unweighted gradient of the increment)
!     here:                  u -= dt M^-1 fp(q)      -> u += dt M^-1 int q div(eps_h N)
! with the same parallel sum, Dirichlet treatment and weighted lumped mass.
!=========================================================================
SUBROUTINE EL_MM_Correct_Velocity_Operator()

  USE EL_FIELDS, ONLY: el_field_data
  USE EL_MM_OPERATOR, ONLY: EL_MM_PHI_OPERATORS
  USE PP3D_MPI, ONLY: myid, master

  REAL*8, ALLOCATABLE :: fp(:,:), fv(:,:)
  INTEGER :: ndof, i

  IF (myid.EQ.master) RETURN
  ndof = QuadSc%ndof
  ALLOCATE(fp(3,ndof), fv(3,ndof))
  CALL EL_MM_PHI_OPERATORS(mg_mesh, NLMAX, mg_mesh%level(NLMAX)%dcorvg, &
    el_field_data%epsilon_transport_q2, LinSc%valP(NLMAX)%x, QuadSc%valU, QuadSc%valV, &
    QuadSc%valW, .FALSE., fp, fv)
  DO i = 1, ndof
    QuadSc%defU(i) = TSTEP*fp(1,i)
    QuadSc%defV(i) = TSTEP*fp(2,i)
    QuadSc%defW(i) = TSTEP*fp(3,i)
  END DO
  DEALLOCATE(fp, fv)

  ILEV = NLMAX
  MlRhoPmat => mg_MlRhoPmat(ILEV)%a
  CALL E013Sum3(QuadSc%defU, QuadSc%defV, QuadSc%defW)
  CALL Boundary_QuadScalar_Def()
  DO i = 1, ndof
    QuadSc%valU(i) = QuadSc%valU(i) - QuadSc%defU(i)/MlRhoPmat(i)
    QuadSc%valV(i) = QuadSc%valV(i) - QuadSc%defV(i)/MlRhoPmat(i)
    QuadSc%valW(i) = QuadSc%valW(i) - QuadSc%defW(i)/MlRhoPmat(i)
  END DO

END SUBROUTINE EL_MM_Correct_Velocity_Operator
