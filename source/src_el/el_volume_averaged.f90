MODULE EL_VA_CONTINUITY

  ! Volume-averaged continuity for the unresolved Euler-Lagrange coupling.
  !
  ! The single-phase solver enforces div u = 0. In a particle-laden fluid the
  ! volume-averaged continuity equation is
  !
  !     d eps/dt + div(eps u) = 0        (eps = fluid volume fraction)
  !
  ! i.e.  div u = g,   g = -( d eps/dt + u . grad eps ) / eps.
  !
  ! Stage 'continuity' (SimPar@ELVolumeAveraged = continuity) keeps the
  ! single-phase momentum operator and imposes the source g through the
  ! pressure Poisson right-hand side: the discrete constraint is B^T u = G
  ! with G_q = int_K q g dV, and the projection solves C dp = -(1/k) B^T u~ + G/k,
  ! so G/k is added to defP after Matdef_General_LinScalar(...,1) and before
  ! it is copied to rhsP (see Transport_q2p1_UxyzP_fluid_core). Only the
  ! constant P1 test function receives the source (the three gradient entries
  ! vanish for the axis-aligned brick cells of every EL mesh).
  !
  ! Inputs are the element fields of EL_FIELDS on level NLMAX, fresh from the
  ! post-advance transfer pass: eps, d eps/dt = (eps^{n+1}-eps^n)/dt and the
  ! cell-averaged kernel gradient grad(alpha_p) (grad eps = -grad alpha_p).
  ! The element velocity is the mean of its eight vertex values.
  !
  ! Compatibility: in a closed or periodic box the pressure matrix has a
  ! constant nullspace, so the source must have zero volume integral. The
  ! continuous form div(eps u) has that property exactly (constant particle
  ! volume), the divided form g loses it at O(phi); its volume integral is
  ! therefore removed through the configured sink distribution (uniform mean
  ! by default, see below). Reductions run over the worker communicator; the
  ! CFD master never calls this routine.
  !
  ! Time integration of the rate (SimPar@ELVASourceTimeScheme, EL_VA_CONTINUITY_RATE):
  ! 'be' uses the backward difference d eps/dt directly, 'bdf2' the endpoint BDF2 value
  ! (1+a) (d eps/dt)^{n+1} - a (d eps/dt)^n, a = dt/(dt + dt_prev) (= 1/2 for equal steps:
  ! 1.5/-0.5), i.e. the two backward differences centred at t^{n+1/2} and t^{n-1/2}
  ! extrapolated to t^{n+1}. The previous rate and its interval are kept in EL_FIELDS
  ! (deps_f_dt_prev, rate_dt_prev, saved at the history advance); BE whenever either rate
  ! is not a genuine backward difference (first step, restart without history).
  !
  ! Compatibility (SimPar@ELVACompatibility, EL_VA_SINK): the volume integral int_g of the raw
  ! source is removed through a sink distribution s_e with sum_e s_e |e| = int_g,
  !   'uniform'  s_e = gbar = int_g / V              (the campaign code, default)
  !   'alpha'    s_e = int_g alpha_p(e) / W, W = sum alpha_p |e|   (sink where the particles are;
  !              uniform fallback when W < 1e-30)
  ! so that the moment of the compensating sink sits at the cloud instead of the box centre
  ! (s04 RUNBOOK section 10: closes the level-3 return-flow deficit). The deck choice is the
  ! requested mode; the effective mode may fall back to uniform at run time (no particles).
  ! Every non-uniform mode is a weighted sink s_e = int_g w_e/W with W = sum w_e |e|: a further
  ! mode is one entry in el_va_compat_names (EL_CONFIG; deck value = label, the deck is
  ! validated against that table), its id below, and one CASE in EL_VA_SINK_WEIGHT. A further time scheme touches EL_VA_USE_*
  ! (history eligibility), the weight, EL_VA_CONTINUITY_RATE and the scheme label.
  !
  ! Source/projection correctors (SimPar@ELVASourceCorrectors) re-enter this routine with
  ! icorr = 1..n from EL_VA_Projection_Correctors (Transport_Q2P1); the routine itself is
  ! pass-agnostic, icorr only labels the log lines and int_g_out returns the volume integral
  ! of the raw source for the residual log.
  !
  ! Continuity form (SimPar@ELVAContinuityForm, EL_VA_FORM_ID of EL_CONFIG): the third,
  ! orthogonal axis.
  !   'advective'    the per-element loop above (campaign code, kept verbatim).
  !   'conservative' div(eps u) = -d eps/dt in the weak form of all four P1 test functions
  !                  q_1 = 1, q_m = x_m - xc_m, with the shared face fluxes of the continuous
  !                  transport porosity eps_h (EL_VA_FLUX; the P0 eps stays the storage
  !                  variable). Per element and row m, with D_m = int_e q_m div u^k read from
  !                  the assembled defp (the quadrature of B itself), the unsinked residual
  !                      r_m = F_m(u^k) + d_cont Q_m,   Q_1 = |e|, Q_m = |e| int_e q_m / |e|_geom
  !                  and the preconditioning coefficient eps_e = MAX(eps, 1e-3), the source is
  !                      gV_m = D_m - r_m/eps_e,  defp_m := -r_m/(eps_e dt)  (row 1: - s_e|e|/dt)
  !                  i.e. defp_m + D_m/dt = gV_m/dt. Row 1 enters the sink bookkeeping as
  !                  g = gV_1/|e| (same slot meaning as the advective form; the "advective"
  !                  term of the audit is then adv = -(eps_e gV_1 + d_cont|e|)/|e|). Iterated
  !                  through the corrector passes it is form C at the fixed point:
  !                  F_m = -d_cont Q_m - delta_m1 eps_e |e| s_e, and summing row 1 gives
  !                  sum_e eps_e |e| s_e = -(S + Phi) with S = sum_e d_cont |e| (storage rate)
  !                  and Phi the boundary flux: the sink vanishes exactly when storage and
  !                  boundary flux are compatible (closed box, no clipping, no relaxation).
  !                  The sink stays on the constant row only (an algebraic compatibility
  !                  correction, not the four-row weak form of a spatial sink). The audit lines
  !                  EL_VA_FLUX / EL_VA_FLUXMOM report S, Phi, the eps-weighted sink integral,
  !                  the unsinked residual norms, the coordinate-moment identities
  !                  sum_e [xc_m F_1 + F_1+m] + int eps_h u_m = 0 (closed box) and both the
  !                  reconstructed-eps and the P0-eps superficial velocities.
  !
  ! Collectives, in order: the campaign block (core, every pass); the temporal block (audit
  ! steps); the weight block (weighted sink modes); the flux block and its max (conservative
  ! form on audit steps); the corrected-source integral (weighted modes on audit steps); the
  ! eps-weighted sink integral (conservative form on audit steps). The blocks are reduced
  ! separately so that the collective pattern - and hence the floating-point result - of the
  ! default configuration is unchanged by any switch. Which collectives are reached depends
  ! only on the deck and the step number, never on rank-local data; the audit section below
  ! the reductions contains no collective.

  USE MPI
  USE PP3D_MPI, ONLY: myid, master, showid, MPI_COMM_SUBS
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_volume_averaged, el_momentum_audit_freq, &
                       el_va_source_time_scheme, el_va_source_correctors, el_va_compatibility, &
                       el_va_compat_names, EL_VA_FORM_ID
  USE EL_FIELDS, ONLY: el_field_data
  USE EL_VA_FLUX, ONLY: EL_VA_FLUX_MOMENTS, EL_VA_CONSERVATIVE_ROWS
  USE EL_GEOMETRY, ONLY: EL_DOMAIN_Z_CENTER

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: EL_ADD_CONTINUITY_SOURCE, EL_VOLUME_AVERAGED_ACTIVE
  PUBLIC :: EL_VA_LOG_CONSTRAINT_RESIDUAL, EL_VA_AUDIT_STEP

  ! Slots of the worker reductions (EL_ADD_CONTINUITY_SOURCE).
  ! core (every pass; this is the campaign reduction and must not change):
  INTEGER, PARAMETER :: N_CORE = 8
  INTEGER, PARAMETER :: C_INT_G = 1    ! sum g |e|                     raw source integral
  INTEGER, PARAMETER :: C_VOL = 2      ! sum |e|                       V
  INTEGER, PARAMETER :: C_L1_DEDT = 3  ! sum |d eps/dt / eps| |e|      L1 of the backward difference
  INTEGER, PARAMETER :: C_L1_ADV = 4   ! sum |adv / eps| |e|           L1 of the advective term
  INTEGER, PARAMETER :: C_M1RAW = 5    ! sum z_c g |e|                 first moment of the raw source
  INTEGER, PARAMETER :: C_Q2 = 6       ! sum dt/2 d_cont^2 / eps |e|   O(dt) temporal defect (T4)
  INTEGER, PARAMETER :: C_ADV_NET = 7  ! sum adv / eps |e|             net advective source
  INTEGER, PARAMETER :: C_S2 = 8       ! sum alpha^2 / eps |e|         O(alpha^2) weight
  ! temporal block (audit steps only):
  INTEGER, PARAMETER :: N_TSRC = 2
  INTEGER, PARAMETER :: T_DEDT_NET = 1 ! sum d_cont / eps |e|          signed temporal source
  INTEGER, PARAMETER :: T_L1_DCONT = 2 ! sum |d_cont / eps| |e|
  ! weight block (weighted sink modes only):
  INTEGER, PARAMETER :: N_WSUM = 2
  INTEGER, PARAMETER :: W_SUM = 1      ! sum w_e |e|                   W
  INTEGER, PARAMETER :: W_ZSUM = 2     ! sum z_c w_e |e|               W z_w

  ! flux block (conservative form, audit steps only; EL_VA_FLUX / EL_VA_FLUXMOM):
  INTEGER, PARAMETER :: N_FLUX = 15
  INTEGER, PARAMETER :: F_S = 1        ! sum d_cont |e|               storage rate S
  INTEGER, PARAMETER :: F_PHI = 2      ! sum F_1                      boundary flux Phi
  INTEGER, PARAMETER :: F_RES1_L1 = 3  ! sum |r_1|                    unsinked residual, row 1
  INTEGER, PARAMETER :: F_RESL_L1 = 4  ! sum (|r_2|+|r_3|+|r_4|)/h_e  linear rows per length
  INTEGER, PARAMETER :: F_ID_X = 5     ! sum xc F_1 + F_2 + int eps_h u_x   moment identities
  INTEGER, PARAMETER :: F_ID_Y = 6
  INTEGER, PARAMETER :: F_ID_Z = 7
  INTEGER, PARAMETER :: F_EPSU_Z = 8   ! sum int_e eps_h u_z          reconstructed-eps flow
  INTEGER, PARAMETER :: F_P0U_Z = 9    ! sum eps_e int_e u_z          P0-eps flow (EL_MEAN_SLIP)
  INTEGER, PARAMETER :: F_M1RATE = 10  ! sum d_cont (xc_z |e| + Q_z)  storage first-moment rate
  INTEGER, PARAMETER :: F_VGEOM = 11   ! sum |e|_geom                 quadrature volume
  INTEGER, PARAMETER :: F_RES1_SQ = 12 ! sum r_1^2
  INTEGER, PARAMETER :: F_L1_F = 13    ! sum |F_1|
  INTEGER, PARAMETER :: F_L1_EPSD = 14 ! sum |eps_e D_1|
  INTEGER, PARAMETER :: F_DIVCHK = 15  ! sum |sum_f int_f u.n - D_1|  face flux vs B divergence

  ! Sink modes: id = index into el_va_compat_names (EL_CONFIG; deck value = audit label, the
  ! deck is validated against the same table); EL_VA_SINK_WEIGHT gives the per-element weight.
  INTEGER, PARAMETER :: SINK_UNIFORM = 0, SINK_ALPHA = 1
  ! Continuity forms: id = index into el_va_form_names (EL_CONFIG, EL_VA_FORM_ID).
  INTEGER, PARAMETER :: FORM_ADVECTIVE = 0, FORM_CONSERVATIVE = 1

  ! longest audit line: EL_VA_FLUXMOM with ' corr= -1', about 275 characters (an internal
  ! WRITE past LOG_LEN is a runtime error, not a silent truncation)
  INTEGER, PARAMETER :: LOG_LEN = 320

CONTAINS

  LOGICAL FUNCTION EL_VOLUME_AVERAGED_ACTIVE()
    EL_VOLUME_AVERAGED_ACTIVE = (TRIM(el_volume_averaged).EQ.'continuity') .AND. &
      ALLOCATED(el_field_data%grad_alpha_p) .AND. ALLOCATED(el_field_data%deps_f_dt)
  END FUNCTION EL_VOLUME_AVERAGED_ACTIVE

  ! .TRUE. on the steps where the EL_VA_* audit lines are written (ELMomentumAuditFreq).
  LOGICAL FUNCTION EL_VA_AUDIT_STEP(istep)
    INTEGER, INTENT(IN) :: istep
    EL_VA_AUDIT_STEP = .FALSE.
    IF (el_momentum_audit_freq.GT.0) &
      EL_VA_AUDIT_STEP = (MOD(ABS(istep), el_momentum_audit_freq).EQ.0)
  END FUNCTION EL_VA_AUDIT_STEP

  ! .TRUE. when the BDF2 rate can be formed: scheme requested and both the current and
  ! the previous rate are genuine backward differences over positive intervals.
  LOGICAL FUNCTION EL_VA_USE_BDF2()
    EL_VA_USE_BDF2 = (TRIM(el_va_source_time_scheme).EQ.'bdf2') .AND. &
      el_field_data%rate_valid .AND. el_field_data%prev_rate_valid .AND. &
      el_field_data%rate_dt.GT.0.0d0 .AND. el_field_data%rate_dt_prev.GT.0.0d0 .AND. &
      ALLOCATED(el_field_data%deps_f_dt_prev)
  END FUNCTION EL_VA_USE_BDF2

  ! BDF2 extrapolation weight a = dt/(dt + dt_prev) (exactly 0.5 for equal steps).
  REAL*8 FUNCTION EL_VA_BDF2_WEIGHT()
    EL_VA_BDF2_WEIGHT = el_field_data%rate_dt/(el_field_data%rate_dt + el_field_data%rate_dt_prev)
  END FUNCTION EL_VA_BDF2_WEIGHT

  ! Rate entering the source of element iel: BE, or BDF2 with weight a.
  REAL*8 FUNCTION EL_VA_CONTINUITY_RATE(iel, use_bdf2, a)
    INTEGER, INTENT(IN) :: iel
    LOGICAL, INTENT(IN) :: use_bdf2
    REAL*8, INTENT(IN) :: a
    EL_VA_CONTINUITY_RATE = el_field_data%deps_f_dt(iel)
    IF (use_bdf2) EL_VA_CONTINUITY_RATE = (1.0d0+a)*el_field_data%deps_f_dt(iel) &
                                        - a*el_field_data%deps_f_dt_prev(iel)
  END FUNCTION EL_VA_CONTINUITY_RATE

  ! Compatibility mode requested by the deck (the effective mode may still fall back to
  ! uniform at run time, see EL_VA_EFFECTIVE_SINK_MODE).
  INTEGER FUNCTION EL_VA_SINK_MODE()
    INTEGER :: m
    EL_VA_SINK_MODE = SINK_UNIFORM
    DO m = 0, UBOUND(el_va_compat_names,1)
      IF (TRIM(el_va_compatibility).EQ.TRIM(el_va_compat_names(m))) EL_VA_SINK_MODE = m
    END DO
  END FUNCTION EL_VA_SINK_MODE

  ! A weighted mode without weight (W < 1e-30: no particles) is uniform. W is the reduced
  ! total, identical on every worker, so the fallback is taken by all workers alike and never
  ! changes which collectives are reached.
  PURE INTEGER FUNCTION EL_VA_EFFECTIVE_SINK_MODE(requested_mode, w_total)
    INTEGER, INTENT(IN) :: requested_mode
    REAL*8, INTENT(IN) :: w_total
    EL_VA_EFFECTIVE_SINK_MODE = requested_mode
    IF (requested_mode.NE.SINK_UNIFORM .AND. w_total.LT.1.0d-30) &
      EL_VA_EFFECTIVE_SINK_MODE = SINK_UNIFORM
  END FUNCTION EL_VA_EFFECTIVE_SINK_MODE

  ! Per-element weight w_e of a weighted sink mode (a new mode is one CASE here).
  PURE REAL*8 FUNCTION EL_VA_SINK_WEIGHT(mode, alpha_e)
    INTEGER, INTENT(IN) :: mode
    REAL*8, INTENT(IN) :: alpha_e
    SELECT CASE (mode)
    CASE (SINK_ALPHA)
      EL_VA_SINK_WEIGHT = alpha_e
    CASE DEFAULT
      EL_VA_SINK_WEIGHT = 0.0d0
    END SELECT
  END FUNCTION EL_VA_SINK_WEIGHT

  ! Sink of an element: the amount removed from its raw source so that the corrected source
  ! has zero volume integral. uniform: gbar (campaign code); weighted: int_g w_e/W (alpha: T9).
  ! The expressions are kept verbatim from the validated runs - do not reorder the operations.
  PURE REAL*8 FUNCTION EL_VA_SINK(mode, gbar, int_g, w_e, w_total)
    INTEGER, INTENT(IN) :: mode
    REAL*8, INTENT(IN) :: gbar, int_g, w_e, w_total
    SELECT CASE (mode)
    CASE (SINK_UNIFORM)
      EL_VA_SINK = gbar
    CASE DEFAULT
      EL_VA_SINK = int_g*w_e/w_total
    END SELECT
  END FUNCTION EL_VA_SINK

  ! One audit line to the protocol file and to the terminal (rank showid only).
  SUBROUTINE EL_VA_WRITE_LINE(mfile, line)
    INTEGER, INTENT(IN) :: mfile
    CHARACTER(LEN=*), INTENT(IN) :: line
    IF (myid.NE.showid) RETURN
    WRITE(mfile,'(A)') TRIM(line)
    WRITE(*,'(A)') TRIM(line)
  END SUBROUTINE EL_VA_WRITE_LINE

  SUBROUTINE EL_ADD_CONTINUITY_SOURCE(defp, val_u, val_v, val_w, mesh, ilev, &
                                      dt, time, istep, mfile, icorr, int_g_out)

    REAL*8, INTENT(INOUT) :: defp(:)
    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, istep, mfile
    REAL*8, INTENT(IN) :: dt, time
    ! pass label: 0 predictor, 1..n corrector, -1 final (diagnostic) assembly; the label is
    ! only written to the log lines when correctors are enabled (default logs unchanged)
    INTEGER, INTENT(IN), OPTIONAL :: icorr
    ! volume integral of the raw source g (before the sink removal)
    REAL*8, INTENT(OUT), OPTIONAL :: int_g_out

    INTEGER :: iel, nel, i, ig, ierr, m
    ! requested_mode: the deck choice, decides which collectives and audit lines exist (the
    ! same on every worker); effective_mode: what is applied, may fall back to uniform
    INTEGER :: requested_mode, effective_mode
    LOGICAL :: weighted
    ! conservative: the continuity form of the deck (the same on every worker)
    LOGICAL :: conservative
    REAL*8 :: fsum(N_FLUX), fsum_g(N_FLUX), fmax, fmax_g, s_e, sink_eps, sink_eps_g
    REAL*8, ALLOCATABLE :: res(:,:), epsk(:)
    REAL*8 :: u(3), eps, deps_dt, d_cont, adv, g, g_corr, vol, zc, alpha, a_bdf2, w_e
    REAL*8 :: core(N_CORE), core_g(N_CORE), tsrc(N_TSRC), tsrc_g(N_TSRC)
    REAL*8 :: wsum(N_WSUM), wsum_g(N_WSUM)
    REAL*8 :: gbar, vol_total, w_total, z_w, mshift, int_gc, int_gc_g
    LOGICAL :: use_bdf2, audit
    CHARACTER(LEN=16) :: corr_tag
    CHARACTER(LEN=LOG_LEN) :: line
    REAL*8, ALLOCATABLE :: gk(:), wk(:)

    IF (PRESENT(int_g_out)) int_g_out = 0.0d0
    IF (myid.EQ.master) RETURN
    ! The conservative form validates the readiness of every worker collectively before the
    ! rank-local returns below (inherited from the campaign code, uniform in practice), so
    ! that an unready worker stops the run coherently instead of leaving a collective.
    conservative = (EL_VA_FORM_ID().EQ.FORM_CONSERVATIVE)
    IF (conservative) CALL EL_VA_CONSERVATIVE_READINESS(defp, mesh, ilev, dt, mfile)
    IF (.NOT.EL_VOLUME_AVERAGED_ACTIVE()) RETURN
    IF (dt.LE.0.0d0) RETURN

    nel = mesh%level(ilev)%nel
    IF (nel.LE.0) RETURN
    IF (SIZE(el_field_data%deps_f_dt).LT.nel) RETURN
    IF (SIZE(defp).LT.4*nel) RETURN

    requested_mode = EL_VA_SINK_MODE()
    effective_mode = requested_mode
    weighted = (requested_mode.NE.SINK_UNIFORM)
    use_bdf2 = EL_VA_USE_BDF2()
    a_bdf2 = 0.0d0
    IF (use_bdf2) a_bdf2 = EL_VA_BDF2_WEIGHT()
    audit = EL_VA_AUDIT_STEP(istep)
    corr_tag = ''
    IF (PRESENT(icorr) .AND. el_va_source_correctors.GT.0) WRITE(corr_tag,'(A,I0)') ' corr= ', icorr

    ALLOCATE(gk(nel))
    IF (weighted) ALLOCATE(wk(nel))
    core = 0.0d0
    tsrc = 0.0d0
    wsum = 0.0d0
    fsum = 0.0d0
    fmax = 0.0d0

    ! ---- raw source per element and the local sums --------------------------------------
    IF (conservative) THEN
      ALLOCATE(res(4,nel), epsk(nel))
      CALL EL_VA_CONSERVATIVE_RAW_SOURCE(defp, val_u, val_v, val_w, mesh, ilev, dt, &
        use_bdf2, a_bdf2, requested_mode, weighted, gk, wk, epsk, res, core, tsrc, wsum, &
        fsum, fmax)
    ELSE
    ! campaign loop, verbatim; indentation kept so that the diff against the parent is empty
    DO iel = 1, nel
      u = 0.0d0
      zc = 0.0d0
      DO i = 1, 8
        ig = mesh%level(ilev)%kvert(i,iel)
        u(1) = u(1) + val_u(ig)
        u(2) = u(2) + val_v(ig)
        u(3) = u(3) + val_w(ig)
        zc = zc + mesh%level(ilev)%dcorvg(3,ig)
      END DO
      u = 0.125d0*u
      zc = 0.125d0*zc
      eps = MAX(el_field_data%epsilon_f(iel), 1.0d-3)
      deps_dt = el_field_data%deps_f_dt(iel)
      d_cont = EL_VA_CONTINUITY_RATE(iel, use_bdf2, a_bdf2)
      ! grad eps = - grad alpha_p  (cell-averaged kernel gradient)
      adv = -DOT_PRODUCT(u, el_field_data%grad_alpha_p(:,iel))
      g = -(d_cont + adv)/eps
      gk(iel) = g
      vol = mesh%level(ilev)%dvol(iel)
      alpha = 1.0d0 - eps
      IF (ALLOCATED(el_field_data%alpha_p)) THEN
        IF (SIZE(el_field_data%alpha_p).GE.nel) alpha = el_field_data%alpha_p(iel)
      END IF
      core(C_INT_G) = core(C_INT_G) + g*vol
      core(C_VOL) = core(C_VOL) + vol
      core(C_L1_DEDT) = core(C_L1_DEDT) + ABS(deps_dt/eps)*vol
      core(C_L1_ADV) = core(C_L1_ADV) + ABS(adv/eps)*vol
      core(C_M1RAW) = core(C_M1RAW) + zc*g*vol
      core(C_Q2) = core(C_Q2) + 0.5d0*dt*d_cont*d_cont/eps*vol
      core(C_ADV_NET) = core(C_ADV_NET) + adv/eps*vol
      core(C_S2) = core(C_S2) + alpha*alpha/eps*vol
      tsrc(T_DEDT_NET) = tsrc(T_DEDT_NET) + d_cont/eps*vol
      tsrc(T_L1_DCONT) = tsrc(T_L1_DCONT) + ABS(d_cont/eps)*vol
      IF (weighted) THEN
        wk(iel) = EL_VA_SINK_WEIGHT(requested_mode, alpha)
        wsum(W_SUM) = wsum(W_SUM) + wk(iel)*vol
        wsum(W_ZSUM) = wsum(W_ZSUM) + zc*wk(iel)*vol
      END IF
    END DO
    END IF

    ! ---- worker reductions ----------------------------------------------------------------
    CALL MPI_Allreduce(core, core_g, N_CORE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    tsrc_g = 0.0d0
    IF (audit) CALL MPI_Allreduce(tsrc, tsrc_g, N_TSRC, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    wsum_g = 0.0d0
    IF (weighted) CALL MPI_Allreduce(wsum, wsum_g, N_WSUM, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    fsum_g = 0.0d0
    fmax_g = 0.0d0
    IF (conservative .AND. audit) THEN
      CALL MPI_Allreduce(fsum, fsum_g, N_FLUX, MPI_DOUBLE_PRECISION, MPI_SUM, &
        MPI_COMM_SUBS, ierr)
      CALL MPI_Allreduce(fmax, fmax_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_SUBS, ierr)
    END IF

    vol_total = MAX(core_g(C_VOL), 1.0d-300)
    gbar = core_g(C_INT_G)/vol_total
    IF (PRESENT(int_g_out)) int_g_out = core_g(C_INT_G)
    w_total = wsum_g(W_SUM)
    effective_mode = EL_VA_EFFECTIVE_SINK_MODE(requested_mode, w_total)
    z_w = 0.0d0
    IF (effective_mode.NE.SINK_UNIFORM) z_w = wsum_g(W_ZSUM)/w_total

    ! ---- corrected source into the P1 dofs -------------------------------------------------
    int_gc = 0.0d0
    sink_eps = 0.0d0
    IF (conservative) THEN
      ! all four rows as the direct residual right-hand side: defp_m := -r_m/(eps_e dt),
      ! the sink on the constant row only (see the module header)
      DO iel = 1, nel
        vol = mesh%level(ilev)%dvol(iel)
        w_e = 0.0d0
        IF (weighted) w_e = wk(iel)
        s_e = EL_VA_SINK(effective_mode, gbar, core_g(C_INT_G), w_e, w_total)
        g_corr = gk(iel) - s_e
        defp(4*(iel-1)+1) = -res(1,iel)/(epsk(iel)*dt) - s_e*vol/dt
        DO m = 2, 4
          defp(4*(iel-1)+m) = -res(m,iel)/(epsk(iel)*dt)
        END DO
        int_gc = int_gc + g_corr*vol
        sink_eps = sink_eps + epsk(iel)*s_e*vol
      END DO
    ELSE
    ! campaign update, verbatim; indentation kept so that the diff against the parent is empty
    DO iel = 1, nel
      vol = mesh%level(ilev)%dvol(iel)
      w_e = 0.0d0
      IF (weighted) w_e = wk(iel)
      g_corr = gk(iel) - EL_VA_SINK(effective_mode, gbar, core_g(C_INT_G), w_e, w_total)
      defp(4*(iel-1)+1) = defp(4*(iel-1)+1) + g_corr*vol/dt
      int_gc = int_gc + g_corr*vol
    END DO
    END IF
    ! integral of the corrected source: a roundoff check (about 1e-20) for the weighted modes,
    ! reduced on audit steps only
    int_gc_g = 0.0d0
    IF (audit .AND. weighted) CALL MPI_Allreduce(int_gc, int_gc_g, 1, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    ! eps-weighted sink integral sum eps_e s_e |e| (= -(S + Phi) at the conservative fixed
    ! point: zero for compatible storage, nonzero under clipping or relaxation); the last
    ! collective of this routine
    sink_eps_g = 0.0d0
    IF (conservative .AND. audit) CALL MPI_Allreduce(sink_eps, sink_eps_g, 1, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)

    ! ---- audit lines (rank showid; no collective below this point) --------------------------
    IF (audit) THEN
      WRITE(line,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A)') &
        'EL_VA_SOURCE t= ', time, ' step= ', istep, &
        ' int_g= ', core_g(C_INT_G), ' gbar= ', gbar, &
        ' L1_dedt= ', core_g(C_L1_DEDT), ' L1_adv= ', core_g(C_L1_ADV), TRIM(corr_tag)
      CALL EL_VA_WRITE_LINE(mfile, line)
      WRITE(line,'(A,ES14.6,A,I0,A,ES14.7,A,ES13.6,A,ES13.6,A,ES13.6,A)') &
        'EL_VA_MOMENT t= ', time, ' step= ', istep, &
        ' M1raw= ', core_g(C_M1RAW), ' Q2= ', core_g(C_Q2), &
        ' adv_net= ', core_g(C_ADV_NET), ' S2= ', core_g(C_S2), TRIM(corr_tag)
      CALL EL_VA_WRITE_LINE(mfile, line)
      WRITE(line,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,A,A)') &
        'EL_VA_TSRC t= ', time, ' step= ', istep, &
        ' dedt_net= ', tsrc_g(T_DEDT_NET), ' L1_dcont= ', tsrc_g(T_L1_DCONT), &
        ' scheme= ', TRIM(MERGE('BDF2', 'BE  ', use_bdf2)), TRIM(corr_tag)
      CALL EL_VA_WRITE_LINE(mfile, line)
      IF (weighted) THEN
        ! mshift: the moment the sink moves relative to the uniform removal (zero under the
        ! no-particle fallback); the field names keep the wording of the T9 record (z_alpha)
        mshift = 0.0d0
        IF (effective_mode.NE.SINK_UNIFORM) mshift = (z_w-EL_DOMAIN_Z_CENTER())*core_g(C_INT_G)
        WRITE(line,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,A,A)') &
          'EL_VA_COMPAT t= ', time, ' step= ', istep, &
          ' int_g_raw= ', core_g(C_INT_G), ' int_g_corr= ', int_gc_g, ' W= ', w_total, &
          ' z_alpha= ', z_w, ' mshift= ', mshift, &
          ' mode= ', TRIM(el_va_compat_names(effective_mode)), TRIM(corr_tag)
        CALL EL_VA_WRITE_LINE(mfile, line)
      END IF
      IF (conservative) THEN
        ! conservation: storage rate, boundary flux, eps-weighted sink, unsinked residuals
        WRITE(line,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A)') &
          'EL_VA_FLUX t= ', time, ' step= ', istep, &
          ' S= ', fsum_g(F_S), ' Phi= ', fsum_g(F_PHI), ' sink_eps= ', sink_eps_g, &
          ' res1_L1= ', fsum_g(F_RES1_L1), ' res1_L2= ', SQRT(fsum_g(F_RES1_SQ)), &
          ' res1_max= ', fmax_g, ' resL_L1= ', fsum_g(F_RESL_L1), TRIM(corr_tag)
        CALL EL_VA_WRITE_LINE(mfile, line)
        ! moments: the three coordinate identities, both superficial velocities, the storage
        ! first-moment rate, the storage-vs-quadrature volume gap
        WRITE(line,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A)') &
          'EL_VA_FLUXMOM t= ', time, ' step= ', istep, &
          ' id_x= ', fsum_g(F_ID_X), ' id_y= ', fsum_g(F_ID_Y), ' id_z= ', fsum_g(F_ID_Z), &
          ' uf_eps_z= ', fsum_g(F_EPSU_Z)/vol_total, ' uf_p0_z= ', fsum_g(F_P0U_Z)/vol_total, &
          ' m1rate= ', fsum_g(F_M1RATE), ' dvol_gap= ', (core_g(C_VOL)-fsum_g(F_VGEOM))/vol_total, &
          ' L1_F= ', fsum_g(F_L1_F), ' L1_epsD= ', fsum_g(F_L1_EPSD), ' divchk= ', fsum_g(F_DIVCHK), &
          TRIM(corr_tag)
        CALL EL_VA_WRITE_LINE(mfile, line)
      END IF
    END IF

    DEALLOCATE(gk)
    IF (ALLOCATED(wk)) DEALLOCATE(wk)
    IF (ALLOCATED(res)) DEALLOCATE(res)
    IF (ALLOCATED(epsk)) DEALLOCATE(epsk)

  END SUBROUTINE EL_ADD_CONTINUITY_SOURCE

  ! Readiness of a worker for the conservative form, reduced over the workers: the EL fields,
  ! the transport porosity, the rate array, defp and dt as the evaluation needs them. A
  ! worker without elements is ready (contributes zero). On any unready worker every worker
  ! stops the run (message on showid). Collective; called by every worker on every pass.
  SUBROUTINE EL_VA_CONSERVATIVE_READINESS(defp, mesh, ilev, dt, mfile)

    REAL*8, INTENT(IN) :: defp(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, mfile
    REAL*8, INTENT(IN) :: dt

    INTEGER :: nel, ndof, ierr
    LOGICAL :: ready, all_ready

    nel = mesh%level(ilev)%nel
    ndof = mesh%level(ilev)%nvt + mesh%level(ilev)%net + mesh%level(ilev)%nat + nel
    ready = EL_VOLUME_AVERAGED_ACTIVE() .AND. dt.GT.0.0d0 .AND. &
            ALLOCATED(el_field_data%epsilon_transport_q2)
    IF (ready) ready = SIZE(el_field_data%epsilon_transport_q2).GE.ndof .AND. &
                       SIZE(el_field_data%deps_f_dt).GE.nel .AND. &
                       SIZE(el_field_data%epsilon_f).GE.nel .AND. SIZE(defp).GE.4*nel .AND. &
                       SIZE(mesh%level(ilev)%dcorvg,2).GE.ndof
    CALL MPI_Allreduce(ready, all_ready, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_SUBS, ierr)
    IF (.NOT.all_ready) THEN
      IF (myid.EQ.showid) THEN
        WRITE(*,'(A)') 'Fatal: a worker is not ready for ELVAContinuityForm = conservative (EL fields, transport porosity, rate or pressure arrays missing or undersized).'
        WRITE(mfile,'(A)') 'Fatal: a worker is not ready for ELVAContinuityForm = conservative.'
      END IF
      CALL MPI_Abort(MPI_COMM_SUBS, 908, ierr)
    END IF

  END SUBROUTINE EL_VA_CONSERVATIVE_READINESS

  ! Raw source of the conservative form (module header): the flux moments of the current
  ! velocity, the unsinked residual r_m and the row-1 source gV_1 = D_1 - r_1/eps_e per
  ! element (EL_VA_CONSERVATIVE_ROWS); fills the same campaign slots as the advective loop
  ! (gk = gV_1/|e|, adv from it), the sink weights and the flux block. Two element centres are
  ! used on purpose: zc, the vertex mean of the campaign slots (M1raw, weights), and xc, the
  ! centre node of the P1 basis of B, for the moments and identities; they coincide on the
  ! straight-edged EL meshes. Local sums only, no collective, no early return.
  SUBROUTINE EL_VA_CONSERVATIVE_RAW_SOURCE(defp, val_u, val_v, val_w, mesh, ilev, dt, &
      use_bdf2, a_bdf2, requested_mode, weighted, gk, wk, epsk, res, core, tsrc, wsum, &
      fsum, fmax)

    REAL*8, INTENT(IN) :: defp(:), val_u(:), val_v(:), val_w(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, requested_mode
    REAL*8, INTENT(IN) :: dt, a_bdf2
    LOGICAL, INTENT(IN) :: use_bdf2, weighted
    REAL*8, INTENT(OUT) :: gk(:), epsk(:), res(:,:)
    REAL*8, ALLOCATABLE, INTENT(INOUT) :: wk(:)
    REAL*8, INTENT(INOUT) :: core(N_CORE), tsrc(N_TSRC), wsum(N_WSUM), fsum(N_FLUX), fmax

    INTEGER :: iel, nel, i, ig, m, ioff, ierr
    REAL*8 :: eps, deps_dt, d_cont, adv, g, gv, vol, vgeom, zc, alpha, h_e
    REAL*8 :: xc(3), q(4), d(4), r(4)
    REAL*8, ALLOCATABLE :: flux(:,:), eps_u(:,:), vol_u(:,:), geom(:,:), div_u(:)

    nel = mesh%level(ilev)%nel
    ioff = mesh%level(ilev)%nvt + mesh%level(ilev)%net + mesh%level(ilev)%nat
    ALLOCATE(flux(4,nel), eps_u(3,nel), vol_u(3,nel), geom(0:3,nel), div_u(nel))
    CALL EL_VA_FLUX_MOMENTS(mesh, ilev, mesh%level(ilev)%dcorvg, &
      el_field_data%epsilon_transport_q2, val_u, val_v, val_w, flux, eps_u, vol_u, geom, &
      div_u, ierr)
    ! unreachable after EL_VA_CONSERVATIVE_READINESS; a programming error, stop the job
    IF (ierr.NE.0) CALL MPI_Abort(MPI_COMM_SUBS, 909, ierr)

    DO iel = 1, nel
      zc = 0.0d0
      DO i = 1, 8
        ig = mesh%level(ilev)%kvert(i,iel)
        zc = zc + mesh%level(ilev)%dcorvg(3,ig)
      END DO
      zc = 0.125d0*zc
      xc = mesh%level(ilev)%dcorvg(:,ioff+iel)
      eps = MAX(el_field_data%epsilon_f(iel), 1.0d-3)
      epsk(iel) = eps
      deps_dt = el_field_data%deps_f_dt(iel)
      d_cont = EL_VA_CONTINUITY_RATE(iel, use_bdf2, a_bdf2)
      vol = mesh%level(ilev)%dvol(iel)
      vgeom = geom(0,iel)
      ! storage moments with the deposition volume |e| = dvol
      q(1) = vol
      q(2:4) = vol*geom(1:3,iel)/vgeom
      DO m = 1, 4
        d(m) = -dt*defp(4*(iel-1)+m)
      END DO
      CALL EL_VA_CONSERVATIVE_ROWS(flux(:,iel), d_cont, q, d, eps, r, gv)
      res(:,iel) = r
      g = gv/vol
      gk(iel) = g
      adv = -(eps*gv + d_cont*vol)/vol
      alpha = 1.0d0 - eps
      IF (ALLOCATED(el_field_data%alpha_p)) THEN
        IF (SIZE(el_field_data%alpha_p).GE.nel) alpha = el_field_data%alpha_p(iel)
      END IF
      core(C_INT_G) = core(C_INT_G) + g*vol
      core(C_VOL) = core(C_VOL) + vol
      core(C_L1_DEDT) = core(C_L1_DEDT) + ABS(deps_dt/eps)*vol
      core(C_L1_ADV) = core(C_L1_ADV) + ABS(adv/eps)*vol
      core(C_M1RAW) = core(C_M1RAW) + zc*g*vol
      core(C_Q2) = core(C_Q2) + 0.5d0*dt*d_cont*d_cont/eps*vol
      core(C_ADV_NET) = core(C_ADV_NET) + adv/eps*vol
      core(C_S2) = core(C_S2) + alpha*alpha/eps*vol
      tsrc(T_DEDT_NET) = tsrc(T_DEDT_NET) + d_cont/eps*vol
      tsrc(T_L1_DCONT) = tsrc(T_L1_DCONT) + ABS(d_cont/eps)*vol
      IF (weighted) THEN
        wk(iel) = EL_VA_SINK_WEIGHT(requested_mode, alpha)
        wsum(W_SUM) = wsum(W_SUM) + wk(iel)*vol
        wsum(W_ZSUM) = wsum(W_ZSUM) + zc*wk(iel)*vol
      END IF
      h_e = vgeom**(1.0d0/3.0d0)
      fsum(F_S) = fsum(F_S) + d_cont*vol
      fsum(F_PHI) = fsum(F_PHI) + flux(1,iel)
      fsum(F_RES1_L1) = fsum(F_RES1_L1) + ABS(r(1))
      fsum(F_RES1_SQ) = fsum(F_RES1_SQ) + r(1)*r(1)
      fmax = MAX(fmax, ABS(r(1)))
      fsum(F_RESL_L1) = fsum(F_RESL_L1) + (ABS(r(2))+ABS(r(3))+ABS(r(4)))/h_e
      DO m = 1, 3
        fsum(F_ID_X-1+m) = fsum(F_ID_X-1+m) + xc(m)*flux(1,iel) + flux(1+m,iel) + eps_u(m,iel)
      END DO
      fsum(F_EPSU_Z) = fsum(F_EPSU_Z) + eps_u(3,iel)
      fsum(F_P0U_Z) = fsum(F_P0U_Z) + el_field_data%epsilon_f(iel)*vol_u(3,iel)
      fsum(F_M1RATE) = fsum(F_M1RATE) + d_cont*(xc(3)*vol + q(4))
      fsum(F_VGEOM) = fsum(F_VGEOM) + vgeom
      fsum(F_L1_F) = fsum(F_L1_F) + ABS(flux(1,iel))
      fsum(F_L1_EPSD) = fsum(F_L1_EPSD) + ABS(eps*d(1))
      fsum(F_DIVCHK) = fsum(F_DIVCHK) + ABS(div_u(iel) - d(1))
    END DO

    DEALLOCATE(flux, eps_u, vol_u, geom, div_u)

  END SUBROUTINE EL_VA_CONSERVATIVE_RAW_SOURCE

  ! Constraint residual of the assembled pressure right-hand side (constant P1 dof of every
  ! element, scaled by dt): r_e = dt*defp_e = -integral_e div u dV + (g_e - s_e)|e| in m^3/s
  ! (defp = -(1/dt) B^T u + source; under the conservative form -r_1/eps_e - s_e|e|).
  ! Written as EL_VA_RESIDUAL on the audit steps only (two Allreduces over the workers).
  SUBROUTINE EL_VA_LOG_CONSTRAINT_RESIDUAL(defp, nel, dt, time, istep, icorr, int_g, mfile)

    REAL*8, INTENT(IN) :: defp(:)
    INTEGER, INTENT(IN) :: nel, istep, icorr, mfile
    REAL*8, INTENT(IN) :: dt, time, int_g

    INTEGER :: iel, ierr
    REAL*8 :: r, l1_local, l1_global, mx_local, mx_global
    CHARACTER(LEN=LOG_LEN) :: line

    IF (myid.EQ.master) RETURN
    IF (.NOT.EL_VA_AUDIT_STEP(istep)) RETURN
    ! no rank-local return beyond this point: the Allreduces below are collective over the
    ! workers (an empty or undersized partition contributes zero)
    l1_local = 0.0d0
    mx_local = 0.0d0
    IF (SIZE(defp).GE.4*nel) THEN
      DO iel = 1, nel
        r = dt*defp(4*(iel-1)+1)
        l1_local = l1_local + ABS(r)
        mx_local = MAX(mx_local, ABS(r))
      END DO
    END IF
    CALL MPI_Allreduce(l1_local, l1_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(mx_local, mx_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      MPI_COMM_SUBS, ierr)
    WRITE(line,'(A,ES14.6,A,I0,A,I0,A,ES13.6,A,ES13.6,A,ES13.6)') &
      'EL_VA_RESIDUAL t= ', time, ' step= ', istep, ' corr= ', icorr, &
      ' L1= ', l1_global, ' max= ', mx_global, ' int_g= ', int_g
    CALL EL_VA_WRITE_LINE(mfile, line)

  END SUBROUTINE EL_VA_LOG_CONSTRAINT_RESIDUAL

END MODULE EL_VA_CONTINUITY
