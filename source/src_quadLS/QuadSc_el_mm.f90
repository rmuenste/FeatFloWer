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
! Not rebuilt by default: the pressure matrix C0 = B^T (rho Ml)^-1 B of start-up
! (MatrixRenewal C1). It is the preconditioner of the conservative continuity iteration
! (EL_MOMENTUM_MODEL:EL_MM_RESIDUAL_DIVISOR); the coefficients stored in mg_CMat do not depend
! on later mass rebuilds. SimPar@ELMMPressureOperator = weighted replaces it every step by
! the operator of the model (EL_MM_Pressure_Operator below).
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
! Pressure matrix of the operator form (model A) under SimPar@ELMMPressureOperator = weighted.
!
! With the weighted lumped mass W = lump(eps rho) the velocity correction of model A is
!     du = dt W^-1 (B - fp) q = dt W^-1 D_eps^T q        (EL_MM_Correct_Velocity_Operator)
! and the conservative continuity residual changes by D_eps du, D_eps the conservative
! divergence: the projection operator is D_eps W^-1 D_eps^T, a discrete div((eps/rho) grad).
! The start-up matrix scaled by the element porosity (divisor eps) misses grad eps . grad q:
! at the edge of a dense bed the flux correction through the bed comes out too large by up to
! 1/eps, the residual iteration alternates and needs several correctors (case 05, moving bed).
!
! Here the matrix is rebuilt with the existing assembly B^T m^-1 B and the nodal mass
!     m = U0^2 / W,     U0 = rho Ml (the unweighted lumped mass, mg_MlPMat of Create_MMat),
! i.e. W divided by the square of the nodal porosity W/U0 (D_eps ~ diag(eps) B^T): a nodal
! approximation of the operator with the sparsity of C0, symmetric, from data each rank
! already holds (no new exchange). The residual divisor is then 1 (EL_MM_RESIDUAL_DIVISOR).
! It stays a preconditioner: the fixed point of the residual iteration does not depend on it.
!
! Called after EL_MM_Update_Fluid_Coefficients (mg_MlRhoPmat holds W): builds C on the levels
! of this rank from m, the parallel matrix on the workers and the coarse solver (the tail of
! the MatrixRenewal C branch of the solver; B, its parallel part and the slip treatment do
! not depend on the mass and are not rebuilt), then puts W back. On return mg_MlRhoPmat holds
! W again and ILEV / MlRhoPmat point to NLMAX as after Create_MRhoMat.
!
! Coarse solver: the previous UMFPACK factorisation (types 2, 4) is released here before the
! new one is made; under MatrixRenewal C1 nobody else does (PARDISO, type 9, releases its own).
! The master receives the porosity on level NLMIN only, so the pressure multigrid must have
! its coarse and its master level there, and the hypre solvers keep a matrix copy of their
! own: EL_MM_Check_Solver_Settings refuses the other settings.
! Collective (every rank, deck-only guard EL_MM_PRESSURE_REBUILD at the call site). The
! rebuild is done every step, also while the porosity does not change (fixed bed).
!=========================================================================
SUBROUTINE EL_MM_Pressure_Operator()

  USE PP3D_MPI, ONLY: myid, master
  USE UMFPackSolver, ONLY: myUmfPack_Free

  TYPE(mg_Matrix), ALLOCATABLE :: weighted(:)
  INTEGER :: jlev, ilev_save, i, n, crs
  REAL*8 :: rho, u0

  IF (.NOT.ALLOCATED(mg_MlPMat) .OR. .NOT.ALLOCATED(mg_MlRhoPmat)) THEN
    WRITE(*,*) 'EL_MM_Pressure_Operator: lumped mass arrays not built on rank ', myid
    STOP 1
  END IF
  rho = Properties%Density(1)
  ilev_save = ILEV

  ALLOCATE(weighted(NLMIN:NLMAX))
  DO jlev = NLMIN, NLMAX
    n = SIZE(mg_MlRhoPmat(jlev)%a)
    IF (SIZE(mg_MlPMat(jlev)%a).NE.n) THEN
      WRITE(*,*) 'EL_MM_Pressure_Operator: mass arrays of different size on level ', jlev
      STOP 1
    END IF
    ALLOCATE(weighted(jlev)%a(n))
    weighted(jlev)%a = mg_MlRhoPmat(jlev)%a
    DO i = 1, n
      u0 = rho*mg_MlPMat(jlev)%a(i)
      mg_MlRhoPmat(jlev)%a(i) = u0*u0/weighted(jlev)%a(i)
    END DO
  END DO

  crs = LinSc%prm%MGprmIn%CrsSolverType
  IF (myid.EQ.master .AND. (crs.EQ.2 .OR. crs.EQ.4)) CALL myUmfPack_Free()
  CALL Create_CMat(QuadSc%knprU, QuadSc%knprV, QuadSc%knprW, LinSc%knprP, &
                   LinSc%prm%MGprmIn%MinLev, crs)
  IF (myid.NE.master) CALL Create_ParCMat(QuadSc%knprU, QuadSc%knprV, QuadSc%knprW)

  DO jlev = NLMIN, NLMAX
    mg_MlRhoPmat(jlev)%a = weighted(jlev)%a
  END DO
  DEALLOCATE(weighted)

  ILEV = ilev_save
  CALL SETLEV(2)
  MlRhoPmat => mg_MlRhoPmat(NLMAX)%a

END SUBROUTINE EL_MM_Pressure_Operator

!=========================================================================
! Solver settings a momentum model relies on, checked once at the start of the first
! Euler-Lagrange step, before any pass (EL_MM_Begin_Step; same deck values on every rank,
! hence a coordinated stop):
!   MatrixRenewal  M >= 1 (the mass arrays exist), K = 2 or 3 (convection rebuilt every step
!                  from the weighted coefficient), C = 1 (the start-up pressure matrix is the
!                  preconditioner; a per-step rebuild would take the weighted mass - the
!                  rebuild of ELMMPressureOperator = weighted is done by the adapter itself
!                  and needs Pres@MGCrsSolverType 1-4 or 9 and Pres@MGMinLev = MGMedLev =
!                  coarsest level, see EL_MM_Pressure_Operator),
!   GAMMA = 0      (the rotational pressure update assumes the unweighted divergence source),
!   convective form of the convection term (the divergence form lacks the compensating
!                  -u div(eps u) of the non-conservative weighted inertia),
!   Newtonian fluid path without viscoelastic stress.
!=========================================================================
SUBROUTINE EL_MM_Check_Solver_Settings()

  USE var_QuadScalar, ONLY: iConvectionForm, GAMMA, bNonNewtonian, bViscoElastic
  USE PP3D_MPI, ONLY: myid, showid, master

  LOGICAL :: ok, tell

  ok = .TRUE.
  ! one line per message on the terminal: the master stops first, so it reports as well
  tell = (myid.EQ.showid .OR. myid.EQ.master)
  IF (myMatrixRenewal%M.LT.1 .OR. (myMatrixRenewal%K.NE.2 .AND. myMatrixRenewal%K.NE.3) .OR. &
      myMatrixRenewal%C.NE.1) THEN
    IF (tell) WRITE(*,'(A,3(A,I0))') 'ELMomentumModel needs SimPar@MatrixRenewal with M >= 1, K = 2 or 3, C = 1;', &
      ' got M', myMatrixRenewal%M, ' K', myMatrixRenewal%K, ' C', myMatrixRenewal%C
    ok = .FALSE.
  END IF
  IF (EL_MM_PRESSURE_REBUILD()) THEN
    ! the per-step rebuild of the pressure matrix (EL_MM_Pressure_Operator)
    IF (ALL(LinSc%prm%MGprmIn%CrsSolverType.NE.(/ 1, 2, 3, 4, 9 /))) THEN
      IF (tell) WRITE(*,'(A,I0)') 'ELMMPressureOperator = weighted needs Pres@MGCrsSolverType '// &
        '1, 2, 3, 4 or 9 (the hypre solvers keep their own copy of the matrix); got ', &
        LinSc%prm%MGprmIn%CrsSolverType
      ok = .FALSE.
    END IF
    IF (LinSc%prm%MGprmIn%MedLev.NE.LinSc%prm%MGprmIn%MinLev .OR. &
        LinSc%prm%MGprmIn%MinLev.NE.NLMIN) THEN
      IF (tell) WRITE(*,'(A,3(A,I0))') 'ELMMPressureOperator = weighted needs Pres@MGMinLev = '// &
        'Pres@MGMedLev = coarsest mesh level (the master receives the porosity there only);', &
        ' got MinLev ', LinSc%prm%MGprmIn%MinLev, ' MedLev ', LinSc%prm%MGprmIn%MedLev, &
        ' coarsest level ', NLMIN
      ok = .FALSE.
    END IF
  END IF
  IF (GAMMA.NE.0.0d0) THEN
    IF (tell) WRITE(*,'(A)') 'ELMomentumModel needs GAMMA = 0.'
    ok = .FALSE.
  END IF
  IF (iConvectionForm.NE.0) THEN
    IF (tell) WRITE(*,'(A)') 'ELMomentumModel needs ELConvectionForm = convective.'
    ok = .FALSE.
  END IF
  IF (bNonNewtonian .OR. bViscoElastic) THEN
    IF (tell) WRITE(*,'(A)') 'ELMomentumModel is implemented for the Newtonian fluid path only (FlowType = Newtonian, no viscoelastic stress).'
    ok = .FALSE.
  END IF
  IF (.NOT.ok) STOP 1

END SUBROUTINE EL_MM_Check_Solver_Settings

!=========================================================================
! EL_MM_Begin_Step - start of every Euler-Lagrange step under a momentum model, before the
! first particle-mesh pass (called from Transport_q2p1_UxyzP_el on every rank, deck-only
! guard). Once per run: the solver-setting check. Whenever the fields were (re)allocated: the
! boundary-face mask of the viscous force. Once per run: the porosity priming.
!
! Priming: a fresh run starts with the porosity field at its initialisation value 1, and the
! pre-advance pass samples the forces before it deposits. Without priming the first particle
! advance would use eps = 1 in the drag closure, in the 1/eps scale of model B and in the
! inertia weight. One deposit-only pass (step label 0: no PE arming, no history advance,
! nothing captured for the fluid) builds the porosity of the initial positions first; the
! pre-advance pass that follows then finds eps_old = eps and a zero rate, which is the state
! of a suspension at rest. Never entered for 'none', whose start-up sequence is unchanged.
!=========================================================================
SUBROUTINE EL_MM_Begin_Step(mfile)

  USE EL_TRANSFER, ONLY: EL_PARTICLE_MESH_PASS
  USE EL_MOMENTUM_MODEL, ONLY: EL_MM_VISCOUS_FORCE
  USE PP3D_MPI, ONLY: myid, master

  INTEGER, INTENT(IN) :: mfile
  LOGICAL, SAVE :: checked = .FALSE., primed = .FALSE.
  REAL*8 :: dummy_velocity(1)

  IF (.NOT.checked) THEN
    CALL EL_MM_Check_Solver_Settings()
    checked = .TRUE.
  END IF

  IF (EL_MM_VISCOUS_FORCE()) CALL EL_MM_Prepare_Boundary_Mask()

  IF (.NOT.primed) THEN
    dummy_velocity = 0.0d0
    IF (myid.EQ.master) THEN
      CALL EL_PARTICLE_MESH_PASS(mg_mesh, NLMAX, dummy_velocity, &
        dummy_velocity, dummy_velocity, dummy_velocity, Properties%Density(1), &
        EL_CLOSURE_VISCOSITY(Properties%Viscosity(1)), Properties%Gravity, tstep, mfile, 0)
    ELSE
      CALL EL_PARTICLE_MESH_PASS(mg_mesh, NLMAX, QuadSc%valU, QuadSc%valV, &
        QuadSc%valW, LinSc%valP(NLMAX)%x, Properties%Density(1), &
        EL_CLOSURE_VISCOSITY(Properties%Viscosity(1)), &
        Properties%Gravity, tstep, mfile, 0)
    END IF
    primed = .TRUE.
  END IF

END SUBROUTINE EL_MM_Begin_Step

!=========================================================================
! Boundary-face mask for the particle viscous force (EL_VISCOUS): the face-centre dof of an
! element face is marked when the face has exactly one adjacent element over all ranks, i.e.
! lies on the physical boundary of the domain (interior faces and subdomain interfaces have
! two). Independent of the boundary condition: walls, inlets, outlets and slip planes all
! cut the kernel support. Built when the mask is not allocated or has the wrong size (first
! step, changed mesh); every worker joins the parallel sum, the master holds no fields.
!=========================================================================
SUBROUTINE EL_MM_Prepare_Boundary_Mask()

  USE EL_FIELDS, ONLY: el_field_data
  USE EL_VA_FLUX, ONLY: EL_VA_ELEMENT_DOFS, NQ2, NFACE
  USE EL_VISCOUS, ONLY: EL_VISCOUS_FACE_CENTRE_DOF
  USE PP3D_MPI, ONLY: myid, master, MPI_COMM_SUBS
  USE MPI, ONLY: MPI_LOGICAL, MPI_LOR

  INTEGER :: i, ndof, iel, f, dofs(NQ2), ilev_save, ierr
  REAL*8, ALLOCATABLE :: count(:)
  LOGICAL :: need_local, need_any

  ! the master holds no fields and is not a member of the workers' communicator
  IF (myid.EQ.master) RETURN
  ! worker-uniform decision (the parallel sum below is collective among the workers)
  need_local = .NOT.ALLOCATED(el_field_data%boundary_face_dof)
  IF (.NOT.need_local) need_local = (SIZE(el_field_data%boundary_face_dof).NE.QuadSc%ndof)
  CALL MPI_Allreduce(need_local, need_any, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_SUBS, ierr)
  IF (.NOT.need_any) RETURN

  ndof = QuadSc%ndof
  ALLOCATE(count(ndof))
  count = 0.0d0
  DO iel = 1, mg_mesh%level(NLMAX)%nel
    CALL EL_VA_ELEMENT_DOFS(mg_mesh, NLMAX, iel, dofs)
    DO f = 1, NFACE
      i = dofs(EL_VISCOUS_FACE_CENTRE_DOF(f))
      count(i) = count(i) + 1.0d0
    END DO
  END DO
  ! parallel sum on the finest level (the global level pinned as in EL_TRANSFER)
  ilev_save = ILEV
  ILEV = NLMAX
  CALL E013Sum(count)
  ILEV = ilev_save

  IF (ALLOCATED(el_field_data%boundary_face_dof)) DEALLOCATE(el_field_data%boundary_face_dof)
  ALLOCATE(el_field_data%boundary_face_dof(ndof))
  DO i = 1, ndof
    el_field_data%boundary_face_dof(i) = (ABS(count(i) - 1.0d0).LT.0.5d0)
  END DO
  DEALLOCATE(count)

END SUBROUTINE EL_MM_Prepare_Boundary_Mask

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
  IF (.NOT.ALLOCATED(el_field_data%epsilon_transport_q2)) THEN
    WRITE(*,*) 'Momentum model A: transport porosity not built on rank ', myid
    STOP 1
  END IF
  ndof = QuadSc%ndof
  ALLOCATE(fp(3,ndof), fv(3,ndof))
  CALL EL_MM_PHI_OPERATORS(mg_mesh, NLMAX, mg_mesh%level(NLMAX)%dcorvg, &
    el_field_data%epsilon_transport_q2, LinSc%valP(NLMAX)%x, QuadSc%valU, QuadSc%valV, &
    QuadSc%valW, .TRUE., fp, fv)
  c = EL_FLUID_OPERATOR_VISCOSITY(Properties%Density(1), &
                                  EL_CLOSURE_VISCOSITY(Properties%Viscosity(1)))
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
! with the same parallel sum, Dirichlet treatment and weighted lumped mass: the tail below
! repeats the one of Velocity_Correction on purpose, so that the shared routine stays as it is.
!=========================================================================
SUBROUTINE EL_MM_Correct_Velocity_Operator()

  USE EL_FIELDS, ONLY: el_field_data
  USE EL_MM_OPERATOR, ONLY: EL_MM_PHI_OPERATORS
  USE PP3D_MPI, ONLY: myid, master

  REAL*8, ALLOCATABLE :: fp(:,:), fv(:,:)
  INTEGER :: ndof, i

  IF (myid.EQ.master) RETURN
  IF (.NOT.ALLOCATED(el_field_data%epsilon_transport_q2)) THEN
    WRITE(*,*) 'Momentum model A: transport porosity not built on rank ', myid
    STOP 1
  END IF
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
