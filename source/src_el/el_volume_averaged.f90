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
  ! Collectives, in order: the campaign block (core, every pass); the temporal block (audit
  ! steps); the weight block (weighted sink modes); the corrected-source integral (weighted
  ! modes on audit steps). The blocks are reduced separately so that the collective pattern -
  ! and hence the floating-point result - of the default configuration is unchanged by any
  ! switch. Which collectives are reached depends only on the deck and the step number, never
  ! on rank-local data; the audit section below the reductions contains no collective.

  USE MPI
  USE PP3D_MPI, ONLY: myid, master, showid, MPI_COMM_SUBS
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_volume_averaged, el_momentum_audit_freq, &
                       el_va_source_time_scheme, el_va_source_correctors, el_va_compatibility, &
                       el_va_compat_names
  USE EL_FIELDS, ONLY: el_field_data
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

  ! Sink modes: id = index into el_va_compat_names (EL_CONFIG; deck value = audit label, the
  ! deck is validated against the same table); EL_VA_SINK_WEIGHT gives the per-element weight.
  INTEGER, PARAMETER :: SINK_UNIFORM = 0, SINK_ALPHA = 1

  ! longest audit line: EL_VA_COMPAT with ' corr= -1', about 190 characters
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

    INTEGER :: iel, nel, i, ig, ierr
    ! requested_mode: the deck choice, decides which collectives and audit lines exist (the
    ! same on every worker); effective_mode: what is applied, may fall back to uniform
    INTEGER :: requested_mode, effective_mode
    LOGICAL :: weighted
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

    ! ---- raw source per element and the local sums --------------------------------------
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

    ! ---- worker reductions ----------------------------------------------------------------
    CALL MPI_Allreduce(core, core_g, N_CORE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    tsrc_g = 0.0d0
    IF (audit) CALL MPI_Allreduce(tsrc, tsrc_g, N_TSRC, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    wsum_g = 0.0d0
    IF (weighted) CALL MPI_Allreduce(wsum, wsum_g, N_WSUM, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)

    vol_total = MAX(core_g(C_VOL), 1.0d-300)
    gbar = core_g(C_INT_G)/vol_total
    IF (PRESENT(int_g_out)) int_g_out = core_g(C_INT_G)
    w_total = wsum_g(W_SUM)
    effective_mode = EL_VA_EFFECTIVE_SINK_MODE(requested_mode, w_total)
    z_w = 0.0d0
    IF (effective_mode.NE.SINK_UNIFORM) z_w = wsum_g(W_ZSUM)/w_total

    ! ---- corrected source into the constant P1 dof -----------------------------------------
    int_gc = 0.0d0
    DO iel = 1, nel
      vol = mesh%level(ilev)%dvol(iel)
      w_e = 0.0d0
      IF (weighted) w_e = wk(iel)
      g_corr = gk(iel) - EL_VA_SINK(effective_mode, gbar, core_g(C_INT_G), w_e, w_total)
      defp(4*(iel-1)+1) = defp(4*(iel-1)+1) + g_corr*vol/dt
      int_gc = int_gc + g_corr*vol
    END DO
    ! integral of the corrected source: a roundoff check (about 1e-20) for the weighted modes,
    ! reduced on audit steps only; the last collective of this routine
    int_gc_g = 0.0d0
    IF (audit .AND. weighted) CALL MPI_Allreduce(int_gc, int_gc_g, 1, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)

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
    END IF

    DEALLOCATE(gk)
    IF (ALLOCATED(wk)) DEALLOCATE(wk)

  END SUBROUTINE EL_ADD_CONTINUITY_SOURCE

  ! Constraint residual of the assembled pressure right-hand side (constant P1 dof of every
  ! element, scaled by dt): r_e = dt*defp_e = integral_e div u dV + (g_e - s_e)|e| in m^3/s.
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
