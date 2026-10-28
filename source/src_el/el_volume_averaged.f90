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
  ! volume), the divided form g loses it at O(phi); the volume mean of g is
  ! therefore removed (one Allreduce over the worker communicator; the CFD
  ! master never calls this routine).
  !
  ! Time integration of the rate (SimPar@ELVASourceTimeScheme, EL_VA_CONTINUITY_RATE):
  ! 'be' uses the backward difference d eps/dt directly, 'bdf2' the endpoint BDF2 value
  ! (1+a) (d eps/dt)^{n+1} - a (d eps/dt)^n, a = dt/(dt + dt_prev) (= 1/2 for equal steps:
  ! 1.5/-0.5), i.e. the two backward differences centred at t^{n+1/2} and t^{n-1/2}
  ! extrapolated to t^{n+1}. The previous rate and its interval are kept in EL_FIELDS
  ! (deps_f_dt_prev, rate_dt_prev, saved at the history advance); BE whenever either rate
  ! is not a genuine backward difference (first step, restart without history).
  ! Source/projection correctors (SimPar@ELVASourceCorrectors) re-enter this routine with
  ! icorr = 1..n from EL_VA_Projection_Correctors (Transport_Q2P1); the routine itself is
  ! pass-agnostic, icorr only labels the log lines and int_g_out returns the volume integral
  ! of the raw source for the residual log.

  USE MPI
  USE PP3D_MPI, ONLY: myid, master, showid, MPI_COMM_SUBS
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_volume_averaged, el_momentum_audit_freq, &
                       el_va_source_time_scheme, el_va_source_correctors
  USE EL_FIELDS, ONLY: el_field_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: EL_ADD_CONTINUITY_SOURCE, EL_VOLUME_AVERAGED_ACTIVE
  PUBLIC :: EL_VA_LOG_CONSTRAINT_RESIDUAL, EL_VA_AUDIT_STEP

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
    ! volume integral of the raw source g (before the mean removal)
    REAL*8, INTENT(OUT), OPTIONAL :: int_g_out

    INTEGER :: iel, nel, i, ig, ierr
    REAL*8 :: u(3), eps, deps_dt, d_cont, adv, g, vol, zc, alpha, a_bdf2
    LOGICAL :: use_bdf2, audit
    ! local(1:8)/global(1:8): the campaign reduction (source, volume, L1 terms, T4 moments);
    ! local(9:10): temporal part of the source (signed, L1), reduced on audit steps only so
    ! that the collective pattern of the default scheme is unchanged
    REAL*8 :: local(10), global(10), gbar, vol_total
    CHARACTER(LEN=16) :: corr_tag
    REAL*8, ALLOCATABLE :: gk(:)

    IF (PRESENT(int_g_out)) int_g_out = 0.0d0
    IF (myid.EQ.master) RETURN
    IF (.NOT.EL_VOLUME_AVERAGED_ACTIVE()) RETURN
    IF (dt.LE.0.0d0) RETURN

    nel = mesh%level(ilev)%nel
    IF (nel.LE.0) RETURN
    IF (SIZE(el_field_data%deps_f_dt).LT.nel) RETURN
    IF (SIZE(defp).LT.4*nel) RETURN

    ALLOCATE(gk(nel))
    local = 0.0d0
    use_bdf2 = EL_VA_USE_BDF2()
    a_bdf2 = 0.0d0
    IF (use_bdf2) a_bdf2 = EL_VA_BDF2_WEIGHT()
    audit = EL_VA_AUDIT_STEP(istep)
    corr_tag = ''
    IF (PRESENT(icorr) .AND. el_va_source_correctors.GT.0) WRITE(corr_tag,'(A,I0)') ' corr= ', icorr

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
      local(1) = local(1) + g*vol
      local(2) = local(2) + vol
      local(3) = local(3) + ABS(deps_dt/eps)*vol
      local(4) = local(4) + ABS(adv/eps)*vol
      ! T4 diagnostics (SPEC 04 section 12.4): first moment of the raw source, the O(dt)
      ! backward-difference term, the net advection source and the O(alpha^2) weight
      alpha = 1.0d0 - eps
      IF (ALLOCATED(el_field_data%alpha_p)) THEN
        IF (SIZE(el_field_data%alpha_p).GE.nel) alpha = el_field_data%alpha_p(iel)
      END IF
      local(5) = local(5) + zc*g*vol
      local(6) = local(6) + 0.5d0*dt*d_cont*d_cont/eps*vol
      local(7) = local(7) + adv/eps*vol
      local(8) = local(8) + alpha*alpha/eps*vol
      ! temporal part of the source (signed and L1), scheme-dependent (EL_VA_TSRC)
      local(9) = local(9) + d_cont/eps*vol
      local(10) = local(10) + ABS(d_cont/eps)*vol
    END DO

    CALL MPI_Allreduce(local(1:8), global(1:8), 8, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    global(9:10) = 0.0d0
    IF (audit) CALL MPI_Allreduce(local(9:10), global(9:10), 2, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_SUBS, ierr)
    vol_total = MAX(global(2), 1.0d-300)
    gbar = global(1)/vol_total
    IF (PRESENT(int_g_out)) int_g_out = global(1)

    DO iel = 1, nel
      vol = mesh%level(ilev)%dvol(iel)
      defp(4*(iel-1)+1) = defp(4*(iel-1)+1) + (gk(iel)-gbar)*vol/dt
    END DO

    IF (audit .AND. myid.EQ.showid) THEN
      WRITE(mfile,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A)') &
        'EL_VA_SOURCE t= ', time, ' step= ', istep, &
        ' int_g= ', global(1), ' gbar= ', gbar, &
        ' L1_dedt= ', global(3), ' L1_adv= ', global(4), TRIM(corr_tag)
      WRITE(*,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6,A)') &
        'EL_VA_SOURCE t= ', time, ' step= ', istep, &
        ' int_g= ', global(1), ' gbar= ', gbar, &
        ' L1_dedt= ', global(3), ' L1_adv= ', global(4), TRIM(corr_tag)
      WRITE(mfile,'(A,ES14.6,A,I0,A,ES14.7,A,ES13.6,A,ES13.6,A,ES13.6,A)') &
        'EL_VA_MOMENT t= ', time, ' step= ', istep, &
        ' M1raw= ', global(5), ' Q2= ', global(6), &
        ' adv_net= ', global(7), ' S2= ', global(8), TRIM(corr_tag)
      WRITE(*,'(A,ES14.6,A,I0,A,ES14.7,A,ES13.6,A,ES13.6,A,ES13.6,A)') &
        'EL_VA_MOMENT t= ', time, ' step= ', istep, &
        ' M1raw= ', global(5), ' Q2= ', global(6), &
        ' adv_net= ', global(7), ' S2= ', global(8), TRIM(corr_tag)
      WRITE(mfile,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,A,A)') &
        'EL_VA_TSRC t= ', time, ' step= ', istep, &
        ' dedt_net= ', global(9), ' L1_dcont= ', global(10), &
        ' scheme= ', TRIM(MERGE('BDF2', 'BE  ', use_bdf2)), TRIM(corr_tag)
      WRITE(*,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,A,A)') &
        'EL_VA_TSRC t= ', time, ' step= ', istep, &
        ' dedt_net= ', global(9), ' L1_dcont= ', global(10), &
        ' scheme= ', TRIM(MERGE('BDF2', 'BE  ', use_bdf2)), TRIM(corr_tag)
    END IF

    DEALLOCATE(gk)

  END SUBROUTINE EL_ADD_CONTINUITY_SOURCE

  ! Constraint residual of the assembled pressure right-hand side (constant P1 dof of every
  ! element, scaled by dt): r_e = dt*defp_e = integral_e div u dV + (g_e - gbar)|e| in m^3/s.
  ! Written as EL_VA_RESIDUAL on the audit steps only (two Allreduces over the workers).
  SUBROUTINE EL_VA_LOG_CONSTRAINT_RESIDUAL(defp, nel, dt, time, istep, icorr, int_g, mfile)

    REAL*8, INTENT(IN) :: defp(:)
    INTEGER, INTENT(IN) :: nel, istep, icorr, mfile
    REAL*8, INTENT(IN) :: dt, time, int_g

    INTEGER :: iel, ierr
    REAL*8 :: r, l1_local, l1_global, mx_local, mx_global

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
    IF (myid.EQ.showid) THEN
      WRITE(mfile,'(A,ES14.6,A,I0,A,I0,A,ES13.6,A,ES13.6,A,ES13.6)') &
        'EL_VA_RESIDUAL t= ', time, ' step= ', istep, ' corr= ', icorr, &
        ' L1= ', l1_global, ' max= ', mx_global, ' int_g= ', int_g
      WRITE(*,'(A,ES14.6,A,I0,A,I0,A,ES13.6,A,ES13.6,A,ES13.6)') &
        'EL_VA_RESIDUAL t= ', time, ' step= ', istep, ' corr= ', icorr, &
        ' L1= ', l1_global, ' max= ', mx_global, ' int_g= ', int_g
    END IF

  END SUBROUTINE EL_VA_LOG_CONSTRAINT_RESIDUAL

END MODULE EL_VA_CONTINUITY
