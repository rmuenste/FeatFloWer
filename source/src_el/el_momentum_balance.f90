!=========================================================================
! EL_MOMENTUM_BALANCE - momentum balance of particles and gas over every macro step (diagnostic)
!
! SimPar@ELMomentumBalanceAudit = Yes (EL_CONFIG; default No). Read-only: every value is formed
! from solver arrays and their products in scratch storage; no operand, operation order or
! communication of the solve changes. Lines written by rank showid after the fluid solve of the
! step, on the terminal and the protocol. Accumulators and snapshots live here; the fluid-side
! extraction is the adapter QuadSc_el_mombal.f90 (Transport_Q2P1), the particle side records from
! EL_TRANSFER (pre-advance pass) and EL_DIAGNOSTICS (EL_NEWTON_PAIR_END).
!
! Conventions. Every force-like value is the MEAN FORCE over the macro step, impulse / dt, in N,
! as a vector (x y z); positive along +z is upward (the column axis). "On the gas" means: the
! value the term adds to the right-hand side of the discrete gas momentum equation, summed over
! the velocity dofs, divided by dt. Momentum values (mom_*) are in kg m/s.
!
! Discrete gas equation of the step (backward Euler, the only scheme the application accepts with
! the key, app_init), per velocity dof i, lumped mass M_i (rho eps-weighted under a momentum
! model), diffusion matrix D, convection matrix K(u~), drag diagonal B_i (drag_B_source):
!   predictor  M_i (u~_i - u^n_i) = dt [ (B p^n)_i - fp_i(p^n) + c fv_i(u^n) + f_i + B_i u^n_i
!                                        - B_i u~_i - (D u~)_i - (K u~)_i ] + other_i - r_i
!   projection M_i (u^{n+1}_i - u~_i) = dt [ (B q)_i - fp_i(q) ]   on the free rows, 0 on Dirichlet
! with p^{n+1} = p^n + q (all corrector passes together), f the explicit particle reaction
! (fluid_feedback_source), B_i u^n_i the explicit partner of ELDragImplicitWeights = hrz (zero
! without it), fp / fv the solid-fraction operators of momentum model A (EL_MM_OPERATOR; zero
! otherwise), r the residual of the nonlinear momentum solve (the defect at u~) and 'other'
! whatever else the right-hand side carries (body forces; zero in case 11). On Dirichlet rows the
! equation is replaced by the prescribed value; the impulse the constraint exerts there is
!   J_D = sum_{D} [ (S u~)_i - rhs_i ] / dt,     S = M + dt (D + K + B).
!
! Global sums. Every rank-local array of an element assembly (defect, rhs, M_i, D u, K u, B p, fp,
! fv, f, B_i) holds PARTIAL values at dofs shared between ranks (and periodic partners): the
! communication sum adds them. The velocity holds the full value on every copy. Hence the sum over
! the local dofs of a partial array, or of a partial array times the velocity, summed over the
! ranks, counts every dof exactly once; nothing here sums a communicated array. The same rule
! holds for the Dirichlet-row sums (every copy of a Dirichlet dof carries the mask). Verified by
! running a case on two partitions (rank-count independence of the first-step sums).
!
! Terms (gas, index MB_*); each is printed as the sum over ALL rows (line GA) and over the
! Dirichlet rows only (line GD); the free-row value is GA - GD:
!   inertia_pred  sum M (u~ - u^n)/dt           inertia_proj  sum M (u^{n+1} - u~)/dt
!   pgrad_old     sum B p^n                     mmp_old       -sum fp(p^n)
!   mmv_old       c sum fv(u^n)                 feedback      sum f
!   partner       sum B_i u^n_i                 rhs_other     sum [(rhs - M u^n)/dt] - the five before
!   dragimp_pred  -sum B_i u~_i                 visc_pred     -sum (D u~)
!   conv_pred     -sum (K u~)                   defect        sum r/dt (rows before zeroing)
!   rhs_full      sum rhs/dt (rows before zeroing; used for J_D)
!   pgrad_q       sum B q                       mmp_q         -sum fp(q)
!   dragimp_end   -sum B_i u^{n+1}_i            pgrad_end     sum B p^{n+1}
!   mmp_end       -sum fp(p^{n+1})
!   hrz_pred = partner + dragimp_pred = -sum B_i (u~_i - u^n_i)        (derived, printed)
!   hrz_end  = partner + dragimp_end  = -sum B_i (u^{n+1}_i - u^n_i)   (derived, printed)
!
! Particles (line P), from the transfer pass before the particle advance and the particle momentum
! balance over the advance (all ranks summed there):
!   drag_expl  closure drag at the start of the step (EL_FORCE_BUDGET convention)
!   drag_impl  drag impulse/dt of the sub-cycled semi-implicit pe update as the CFD side computes it
!              for the fluid (drag_force_eff, issue-D fix); this is the drag the particles are
!              ARMED with. pe integrates the drag sub-steps on a free-flight reference velocity
!              (HardContactEulerLagrange.h), so contacts inside the sub-steps do not change the
!              drag impulse a particle receives; pe does not report it per body
!   pressure, viscous, lift    the other hydrodynamic forces armed in pe (model A: -V grad p, V div tau)
!   weight     the buoyant weight (rho_p - rho_f) V g armed in pe
!   hydro      drag_impl + pressure + viscous + lift  (H_p of the F1 budget, RUNBOOK section 15)
!   fb_used    sum of the explicit feedback the particles deposit (hydro reaction actually deposited)
!   deposit    sum of the deposited nodal source (= -fb_used up to the kernel normalisation)
!   dp_dt      measured change of the total particle momentum over the pe advance / dt
!   floor_bal  dp_dt - (hydro + weight): the force of all walls from the particle balance (bed audit)
!   wall_contact_raw, wall_lub_raw   -(impulse ON the z walls)/dt as pe reports it; the contact
!              impulse is the unrelaxed PGS value, pe applies PE_WALL_RELAXATION times it in a
!              multi-process run (pe HardContactEulerLagrange.h, restitution notes)
!
! Boundaries and storage (line B), after the fluid solve of the step (u^{n+1}, p^{n+1}):
!   p_in_out, area_in_out   area-mean P1 pressure and area of the bottom / top plane (EL_BED_AUDIT)
!   dpA        (p_in - p_out) area_in
!   flux_in/out   sum over the plane faces of rho_eps_e u (u.n) dS, n outward, rho_eps_e the
!              element coefficient of the convection operator: the momentum the flow carries OUT
!              of the domain per unit time (the gas receives -flux)
!   trac_in/out   sum over the plane faces of c (grad u) n dS (unweighted Laplacian form of D,
!              c the operator viscosity): the viscous force the outside exerts on the gas
!   eps_in_out area mean of the transport porosity eps_h on the plane (the Q1 field the model-A
!              operator fp uses; 1 without that field): separates rec_pB from rec_pfp, since
!              -sum_all fp(p) ~ -P_p + p_in (1 - eps_in) A in the continuous limit
!   fbD_class  z part of feedback_D (the particle reaction on Dirichlet rows) on vertex, edge and
!              face dofs (N; centre dofs never lie on a plane)
! The momentum values mom_start / mom_pred / mom_end (sum M u^n, sum M u~, sum M u^{n+1}, this
! lumped mass of the step) are in lines GA / GD.
!
! Closures (line C):
!   res_gas_F  = inertia_F - [pgrad_end + mmp_end + mmv_old + visc_pred + conv_pred + feedback
!                + partner + dragimp_pred + rhs_other]_F, F = free rows. Holds to the SOLVER
!                TOLERANCE: it equals -defect_F (the residual of the nonlinear momentum solve).
!   gas_roundoff = res_gas_F + defect_F: holds to ROUND-OFF (projection division by the lumped
!                mass, p^{n+1} - p^n = q, the split of the rhs).
!   proj_F     = inertia_proj_F - (pgrad_q + mmp_q)_F: ROUND-OFF. proj_D = inertia_proj_D: exactly 0.
!   jd         J_D = -(defect_D + rhs_full_D); jd_check = jd - {inertia_pred_D - [pgrad_old +
!                mmp_old + mmv_old + feedback + partner + dragimp_pred + visc_pred + conv_pred
!                + rhs_other]_D}: ROUND-OFF.
!   res_p      = dp_dt - (hydro + weight + PE_WALL_RELAXATION wall_contact_raw + wall_lub_raw):
!                the particle balance against the own wall impulse of pe. Holds only to the accuracy of
!                that bookkeeping (the factor is the fold of a multi-process pe run; the application
!                refuses the audit with a single worker, where pe applies the impulse in full);
!                res_p_raw: the same with factor 1.
!   Reconstruction of hydro - dpA (the F1 residual) from the free-row gas balance; the terms add
!   up to it to round-off plus the solver residual (which is itself one of the terms):
!     rec_pB     = pgrad_end_F - dpA              unweighted pressure term vs the plane pressures
!     rec_pfp    = pressure + mmp_end_F           particle pressure force vs model-A solid share
!     rec_visc   = viscous + (mmv_old + visc_pred)_F
!     rec_conv   = conv_pred_F
!     rec_fbD    = -feedback_D                    particle reaction deposited on Dirichlet rows
!     rec_hrz    = (partner + dragimp_pred)_F     implicit drag increment, not mirrored
!     rec_other  = rhs_other_F
!     rec_solver = -defect_F
!     rec_inert  = -(inertia_pred + inertia_proj)_F
!     rec_dep    = drag_impl + lift + feedback_A  deposit vs the armed drag + lift
!     rec_sum    = the sum of the ten;  hydro_dpA = hydro - dpA
!
! Collectives: EL_MOMBAL_WRITE reduces over MPI_COMM_SUBS; it is called by every worker under
! the deck value only (never by the master, which holds no fields); the record routines are local.
!=========================================================================
MODULE EL_MOMENTUM_BALANCE

  USE MPI
  USE PP3D_MPI, ONLY: myid, showid, MPI_COMM_SUBS
  USE TYPES, ONLY: tMultiMesh
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
  USE EL_CONFIG, ONLY: el_momentum_balance_audit
  USE EL_VA_FLUX, ONLY: tElVaFluxTables, EL_VA_FLUX_TABLES, EL_VA_ELEMENT_DOFS, &
                        EL_VA_FACE_GEOMETRY, EL_VA_INVERT3, NQ2, NFACE, NGPF
#ifdef HAVE_PE
  USE DEM_QUERY, ONLY: getElWallContactImpulse, getElWallLubImpulse
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EL_MOMBAL_ACTIVE, EL_MOMBAL_BEGIN_STEP, EL_MOMBAL_SET_PREDICTOR, EL_MOMBAL_RECORD, &
            EL_MOMBAL_RECORD_PLANES, EL_MOMBAL_RECORD_PARTICLES, EL_MOMBAL_RECORD_PARTICLE_DP, &
            EL_MOMBAL_WRITE, EL_MOMBAL_PLANE_LOCAL, EL_MOMBAL_SUMS, EL_MOMBAL_RECORD_FBD_CLASS
  ! snapshots of the step (read by the adapter, written here only)
  PUBLIC :: mb_u_old, mb_u_pred, mb_p_old

  ! gas terms (see the header for the meaning and sign of each)
  INTEGER, PARAMETER, PUBLIC :: MB_MOM_START = 1, MB_MOM_PRED = 2, MB_MOM_END = 3, &
    MB_INERTIA_PRED = 4, MB_INERTIA_PROJ = 5, MB_PGRAD_OLD = 6, MB_MMP_OLD = 7, MB_MMV_OLD = 8, &
    MB_FEEDBACK = 9, MB_PARTNER = 10, MB_RHS_OTHER = 11, MB_DRAGIMP_PRED = 12, &
    MB_VISC_PRED = 13, MB_CONV_PRED = 14, MB_DEFECT = 15, MB_RHS_FULL = 16, MB_PGRAD_Q = 17, &
    MB_MMP_Q = 18, MB_DRAGIMP_END = 19, MB_PGRAD_END = 20, MB_MMP_END = 21
  INTEGER, PARAMETER :: MB_NG = 21
  CHARACTER(LEN=13), PARAMETER :: gas_names(MB_NG) = (/ &
    'mom_start    ', 'mom_pred     ', 'mom_end      ', 'inertia_pred ', 'inertia_proj ', &
    'pgrad_old    ', 'mmp_old      ', 'mmv_old      ', 'feedback     ', 'partner      ', &
    'rhs_other    ', 'dragimp_pred ', 'visc_pred    ', 'conv_pred    ', 'defect       ', &
    'rhs_full     ', 'pgrad_q      ', 'mmp_q        ', 'dragimp_end  ', 'pgrad_end    ', &
    'mmp_end      ' /)

  ! planes: 1 area_in, 2 int p in, 3 area_out, 4 int p out, 5:7 flux_in, 8:10 flux_out,
  ! 11:13 trac_in, 14:16 trac_out, 17 int eps_h in, 18 int eps_h out
  INTEGER, PARAMETER, PUBLIC :: MB_NPLANE = 18

  ! pe folds every contact correction with this relaxation in a multi-process run; its wall
  ! contact diagnostics report the unrelaxed impulse (libs/pe HardContactEulerLagrange.h)
  REAL*8, PARAMETER :: PE_WALL_RELAXATION = 0.9d0

  REAL*8, ALLOCATABLE, SAVE :: mb_u_old(:,:), mb_u_pred(:,:), mb_p_old(:)
  LOGICAL, ALLOCATABLE, SAVE :: dmask(:,:)
  ! rank-local partial sums: (component, 1 all rows / 2 Dirichlet rows, term)
  REAL*8, SAVE :: gas_sum(3, 2, MB_NG) = 0.0d0
  LOGICAL, SAVE :: gas_set(MB_NG) = .FALSE.
  REAL*8, SAVE :: plane_loc(MB_NPLANE) = 0.0d0
  LOGICAL, SAVE :: planes_set = .FALSE.
  ! z of the particle reaction on Dirichlet rows by dof class (vertex, edge, face), rank-local
  REAL*8, SAVE :: fbd_class(3) = 0.0d0

  ! particle values, already global (summed over the ranks by the caller of the recording routine)
  INTEGER, SAVE :: part_step = -1, dp_step = -1
  REAL*8, SAVE :: p_drag_expl(3) = 0.0d0, p_drag_impl(3) = 0.0d0, p_pressure(3) = 0.0d0
  REAL*8, SAVE :: p_viscous(3) = 0.0d0, p_lift(3) = 0.0d0, p_weight(3) = 0.0d0
  REAL*8, SAVE :: p_fb_used(3) = 0.0d0, p_deposit(3) = 0.0d0
  REAL*8, SAVE :: p_dp(3) = 0.0d0, p_expected(3) = 0.0d0

  TYPE(tElVaFluxTables), SAVE :: tables

CONTAINS

  LOGICAL FUNCTION EL_MOMBAL_ACTIVE()
    EL_MOMBAL_ACTIVE = el_momentum_balance_audit
  END FUNCTION EL_MOMBAL_ACTIVE

  ! Start of the fluid step: Dirichlet mask (per component), u^n and p^n; clears the gas terms.
  SUBROUTINE EL_MOMBAL_BEGIN_STEP(mask, u, v, w, p)
    LOGICAL, INTENT(IN) :: mask(:,:)
    REAL*8, INTENT(IN) :: u(:), v(:), w(:), p(:)
    INTEGER :: n

    n = SIZE(mask, 2)
    IF (ALLOCATED(dmask)) THEN
      IF (SIZE(dmask, 2).NE.n) DEALLOCATE(dmask, mb_u_old, mb_u_pred)
    END IF
    IF (.NOT.ALLOCATED(dmask)) ALLOCATE(dmask(3, n), mb_u_old(3, n), mb_u_pred(3, n))
    IF (ALLOCATED(mb_p_old)) THEN
      IF (SIZE(mb_p_old).NE.SIZE(p)) DEALLOCATE(mb_p_old)
    END IF
    IF (.NOT.ALLOCATED(mb_p_old)) ALLOCATE(mb_p_old(SIZE(p)))
    dmask = mask
    mb_u_old(1,:) = u(1:n)
    mb_u_old(2,:) = v(1:n)
    mb_u_old(3,:) = w(1:n)
    mb_u_pred = mb_u_old
    mb_p_old = p
    gas_sum = 0.0d0
    gas_set = .FALSE.
    plane_loc = 0.0d0
    planes_set = .FALSE.
    fbd_class = 0.0d0
  END SUBROUTINE EL_MOMBAL_BEGIN_STEP

  ! The converged predictor u~ (the adapter calls it with every nonlinear iterate; the last wins).
  SUBROUTINE EL_MOMBAL_SET_PREDICTOR(u, v, w)
    REAL*8, INTENT(IN) :: u(:), v(:), w(:)
    INTEGER :: n
    IF (.NOT.ALLOCATED(mb_u_pred)) RETURN
    n = SIZE(mb_u_pred, 2)
    mb_u_pred(1,:) = u(1:n)
    mb_u_pred(2,:) = v(1:n)
    mb_u_pred(3,:) = w(1:n)
  END SUBROUTINE EL_MOMBAL_SET_PREDICTOR

  ! Sums over all rows and over the Dirichlet rows of scale * f (f(3,n), partial values).
  PURE SUBROUTINE EL_MOMBAL_SUMS(f, mask, scale, sum_all, sum_dir)
    REAL*8, INTENT(IN) :: f(:,:), scale
    LOGICAL, INTENT(IN) :: mask(:,:)
    REAL*8, INTENT(OUT) :: sum_all(3), sum_dir(3)
    INTEGER :: i, c
    sum_all = 0.0d0
    sum_dir = 0.0d0
    DO i = 1, SIZE(f, 2)
      DO c = 1, 3
        sum_all(c) = sum_all(c) + f(c, i)
        IF (mask(c, i)) sum_dir(c) = sum_dir(c) + f(c, i)
      END DO
    END DO
    sum_all = scale*sum_all
    sum_dir = scale*sum_dir
  END SUBROUTINE EL_MOMBAL_SUMS

  ! Gas term `term` := scale * (sums of f) (set, not added: a term is recorded once per step, the
  ! defect of the nonlinear loop once per iteration, the last one wins).
  SUBROUTINE EL_MOMBAL_RECORD(term, f, scale)
    INTEGER, INTENT(IN) :: term
    REAL*8, INTENT(IN) :: f(:,:), scale
    IF (.NOT.ALLOCATED(dmask)) RETURN
    CALL EL_MOMBAL_SUMS(f, dmask, scale, gas_sum(:,1,term), gas_sum(:,2,term))
    gas_set(term) = .TRUE.
  END SUBROUTINE EL_MOMBAL_RECORD

  SUBROUTINE EL_MOMBAL_RECORD_FBD_CLASS(vals)
    REAL*8, INTENT(IN) :: vals(3)
    fbd_class = vals
  END SUBROUTINE EL_MOMBAL_RECORD_FBD_CLASS

  SUBROUTINE EL_MOMBAL_RECORD_PLANES(vals)
    REAL*8, INTENT(IN) :: vals(MB_NPLANE)
    plane_loc = vals
    planes_set = .TRUE.
  END SUBROUTINE EL_MOMBAL_RECORD_PLANES

  ! Particle side, transfer pass before the advance of step istep (values summed over the ranks).
  SUBROUTINE EL_MOMBAL_RECORD_PARTICLES(istep, drag_expl, drag_impl, pressure, viscous, lift, &
                                        weight, fb_used, deposit)
    INTEGER, INTENT(IN) :: istep
    REAL*8, INTENT(IN) :: drag_expl(3), drag_impl(3), pressure(3), viscous(3), lift(3), &
                          weight(3), fb_used(3), deposit(3)
    IF (.NOT.el_momentum_balance_audit) RETURN
    part_step = istep
    p_drag_expl = drag_expl
    p_drag_impl = drag_impl
    p_pressure = pressure
    p_viscous = viscous
    p_lift = lift
    p_weight = weight
    p_fb_used = fb_used
    p_deposit = deposit
  END SUBROUTINE EL_MOMBAL_RECORD_PARTICLES

  ! Particle momentum change over the pe advance of step istep and the impulse armed (N s).
  SUBROUTINE EL_MOMBAL_RECORD_PARTICLE_DP(istep, dp, expected)
    INTEGER, INTENT(IN) :: istep
    REAL*8, INTENT(IN) :: dp(3), expected(3)
    IF (.NOT.el_momentum_balance_audit) RETURN
    dp_step = istep
    p_dp = dp
    p_expected = expected
  END SUBROUTINE EL_MOMBAL_RECORD_PARTICLE_DP

  ! Momentum flux and viscous traction through the element faces of this rank that lie in the
  ! plane z = zplane (all 9 face points within tol): area, sum rho_eps_e u (u.n) dS and
  ! sum c (grad u) n dS, n the outward normal of the element (the outward normal of the domain on a
  ! boundary plane). coor(3,*) the Q2 dof coordinates, coef_e(nel) the element coefficient of the
  ! convection operator, visc the operator viscosity c. Isoparametric Q2 map, 3 x 3 Gauss points.
  ! eps_q2 / eps_int (optional, together): the integral over the faces of the nodal field eps_q2.
  SUBROUTINE EL_MOMBAL_PLANE_LOCAL(mesh, ilev, coor, val_u, val_v, val_w, coef_e, visc, zplane, &
                                   tol, area, flux, traction, eps_q2, eps_int)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: coor(:,:), val_u(:), val_v(:), val_w(:), coef_e(:), visc, zplane, tol
    REAL*8, INTENT(INOUT) :: area, flux(3), traction(3)
    REAL*8, INTENT(IN), OPTIONAL :: eps_q2(:)
    REAL*8, INTENT(INOUT), OPTIONAL :: eps_int

    INTEGER :: iel, f, p, l, k, dofs(NQ2)
    REAL*8 :: x(3,NQ2), u(3,NQ2), jac(3,3), inv(3,3), detj, det2, nds(3), uq(3), gu(3,3), dn(3)
    REAL*8 :: e(NQ2)
    LOGICAL :: on_plane, with_eps

    IF (.NOT.tables%ready) CALL EL_VA_FLUX_TABLES(tables)
    with_eps = PRESENT(eps_q2) .AND. PRESENT(eps_int)
    DO iel = 1, mesh%level(ilev)%nel
      CALL EL_VA_ELEMENT_DOFS(mesh, ilev, iel, dofs)
      DO l = 1, NQ2
        x(:,l) = coor(:, dofs(l))
      END DO
      DO f = 1, NFACE
        on_plane = .TRUE.
        DO p = 1, NGPF
          IF (ABS(DOT_PRODUCT(x(3,:), tables%fval(:,p,f)) - zplane).GE.tol) on_plane = .FALSE.
        END DO
        IF (.NOT.on_plane) CYCLE
        DO l = 1, NQ2
          u(1,l) = val_u(dofs(l))
          u(2,l) = val_v(dofs(l))
          u(3,l) = val_w(dofs(l))
          IF (with_eps) e(l) = eps_q2(dofs(l))
        END DO
        DO p = 1, NGPF
          CALL EL_VA_FACE_GEOMETRY(tables, f, p, x, jac, detj, nds)
          CALL EL_VA_INVERT3(jac, inv, det2)
          uq = 0.0d0
          gu = 0.0d0
          DO l = 1, NQ2
            uq = uq + u(:,l)*tables%fval(l,p,f)
            ! d N_l / d x_k = sum_m dref(m,l) inv(m,k)
            DO k = 1, 3
              dn(k) = tables%fder(1,l,p,f)*inv(1,k) + tables%fder(2,l,p,f)*inv(2,k) + &
                      tables%fder(3,l,p,f)*inv(3,k)
            END DO
            DO k = 1, 3
              gu(:,k) = gu(:,k) + u(:,l)*dn(k)
            END DO
          END DO
          area = area + SQRT(SUM(nds*nds))
          IF (with_eps) eps_int = eps_int + DOT_PRODUCT(e, tables%fval(:,p,f))*SQRT(SUM(nds*nds))
          flux = flux + coef_e(iel)*uq*DOT_PRODUCT(uq, nds)
          traction = traction + visc*MATMUL(gu, nds)
        END DO
      END DO
    END DO
  END SUBROUTINE EL_MOMBAL_PLANE_LOCAL

  ! Reduce and write the lines of step istep (collective over the workers, see the header).
  SUBROUTINE EL_MOMBAL_WRITE(time, istep, dt, mfile)
    REAL*8, INTENT(IN) :: time, dt
    INTEGER, INTENT(IN) :: istep, mfile

    INTEGER, PARAMETER :: NBUF = 6*MB_NG + MB_NPLANE + 12 + 3
    REAL*8 :: loc(NBUF), glob(NBUF), wall(3,4), wtmp(3), nan
    REAL*8 :: ga(3, MB_NG), gd(3, MB_NG), gf(3, MB_NG), pl(MB_NPLANE)
    REAL*8 :: pp(3,8), hydro(3), floor_bal(3), dp_dt(3), expected_dt(3), wc(3), wl(3), res_p(3), &
              res_p_raw(3)
    REAL*8 :: p_in, p_out, eps_in, eps_out, fbdc(3), dpa(3), terms_f(3), inertia_f(3), res_gas(3), proj_f(3), jd(3), jdc(3)
    REAL*8 :: rec(3, 10), hydro_dpa(3), hrz_pred(3, 2), hrz_end(3, 2)
    INTEGER :: ierr, k, iw
    LOGICAL :: have_p, have_dp

    IF (.NOT.el_momentum_balance_audit) RETURN
    nan = IEEE_VALUE(1.0d0, IEEE_QUIET_NAN)

    wall = 0.0d0
#ifdef HAVE_PE
    DO iw = 0, 1
      CALL getElWallContactImpulse(iw, wtmp)
      wall(:, 1+iw) = wtmp
      CALL getElWallLubImpulse(iw, wtmp)
      wall(:, 3+iw) = wtmp
    END DO
#endif
    loc(1:3*MB_NG) = RESHAPE(gas_sum(:,1,:), (/3*MB_NG/))
    loc(3*MB_NG+1:6*MB_NG) = RESHAPE(gas_sum(:,2,:), (/3*MB_NG/))
    loc(6*MB_NG+1:6*MB_NG+MB_NPLANE) = plane_loc
    loc(6*MB_NG+MB_NPLANE+1:6*MB_NG+MB_NPLANE+12) = RESHAPE(wall, (/12/))
    loc(NBUF-2:NBUF) = fbd_class
    CALL MPI_Allreduce(loc, glob, NBUF, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    ga = RESHAPE(glob(1:3*MB_NG), (/3, MB_NG/))
    gd = RESHAPE(glob(3*MB_NG+1:6*MB_NG), (/3, MB_NG/))
    pl = glob(6*MB_NG+1:6*MB_NG+MB_NPLANE)
    wall = RESHAPE(glob(6*MB_NG+MB_NPLANE+1:6*MB_NG+MB_NPLANE+12), (/3, 4/))
    fbdc = glob(NBUF-2:NBUF)
    IF (myid.NE.showid) RETURN

    ! a term not recorded this step is printed as NaN (never as 0)
    DO k = 1, MB_NG
      IF (.NOT.gas_set(k)) THEN
        ga(:,k) = nan
        gd(:,k) = nan
      END IF
    END DO
    IF (.NOT.planes_set) pl = nan
    gf = ga - gd
    hrz_pred(:,1) = ga(:,MB_PARTNER) + ga(:,MB_DRAGIMP_PRED)
    hrz_pred(:,2) = gd(:,MB_PARTNER) + gd(:,MB_DRAGIMP_PRED)
    hrz_end(:,1) = ga(:,MB_PARTNER) + ga(:,MB_DRAGIMP_END)
    hrz_end(:,2) = gd(:,MB_PARTNER) + gd(:,MB_DRAGIMP_END)

    ! particles (values of another step are printed as NaN)
    have_p = (part_step.EQ.istep)
    have_dp = (dp_step.EQ.istep) .AND. dt.GT.0.0d0
    pp(:,1) = p_drag_expl
    pp(:,2) = p_drag_impl
    pp(:,3) = p_pressure
    pp(:,4) = p_viscous
    pp(:,5) = p_lift
    pp(:,6) = p_weight
    pp(:,7) = p_fb_used
    pp(:,8) = p_deposit
    IF (.NOT.have_p) pp = nan
    hydro = pp(:,2) + pp(:,3) + pp(:,4) + pp(:,5)
    dp_dt = nan
    expected_dt = nan
    IF (have_dp) THEN
      dp_dt = p_dp/dt
      expected_dt = p_expected/dt
    END IF
    floor_bal = dp_dt - (hydro + pp(:,6))
    wc = -(wall(:,1) + wall(:,2))/dt
    wl = -(wall(:,3) + wall(:,4))/dt
    res_p_raw = dp_dt - (hydro + pp(:,6) + wc + wl)
    res_p = dp_dt - (hydro + pp(:,6) + PE_WALL_RELAXATION*wc + wl)

    ! boundaries
    p_in = nan
    p_out = nan
    IF (pl(1).GT.0.0d0) p_in = pl(2)/pl(1)
    IF (pl(3).GT.0.0d0) p_out = pl(4)/pl(3)
    eps_in = nan
    eps_out = nan
    IF (pl(1).GT.0.0d0) eps_in = pl(17)/pl(1)
    IF (pl(3).GT.0.0d0) eps_out = pl(18)/pl(3)
    IF (.NOT.gas_set(MB_FEEDBACK)) fbdc = nan
    dpa = 0.0d0
    dpa(3) = (p_in - p_out)*pl(1)

    ! gas closures (free rows F = all - Dirichlet)
    inertia_f = gf(:,MB_INERTIA_PRED) + gf(:,MB_INERTIA_PROJ)
    terms_f = gf(:,MB_PGRAD_END) + gf(:,MB_MMP_END) + gf(:,MB_MMV_OLD) + gf(:,MB_VISC_PRED) + &
              gf(:,MB_CONV_PRED) + gf(:,MB_FEEDBACK) + gf(:,MB_PARTNER) + gf(:,MB_DRAGIMP_PRED) + &
              gf(:,MB_RHS_OTHER)
    res_gas = inertia_f - terms_f
    proj_f = gf(:,MB_INERTIA_PROJ) - (gf(:,MB_PGRAD_Q) + gf(:,MB_MMP_Q))
    jd = -(gd(:,MB_DEFECT) + gd(:,MB_RHS_FULL))
    jdc = jd - (gd(:,MB_INERTIA_PRED) - (gd(:,MB_PGRAD_OLD) + gd(:,MB_MMP_OLD) + gd(:,MB_MMV_OLD) + &
          gd(:,MB_FEEDBACK) + gd(:,MB_PARTNER) + gd(:,MB_DRAGIMP_PRED) + gd(:,MB_VISC_PRED) + &
          gd(:,MB_CONV_PRED) + gd(:,MB_RHS_OTHER)))

    ! reconstruction of hydro - dpA
    rec(:,1) = gf(:,MB_PGRAD_END) - dpa
    rec(:,2) = pp(:,3) + gf(:,MB_MMP_END)
    rec(:,3) = pp(:,4) + gf(:,MB_MMV_OLD) + gf(:,MB_VISC_PRED)
    rec(:,4) = gf(:,MB_CONV_PRED)
    rec(:,5) = -gd(:,MB_FEEDBACK)
    rec(:,6) = gf(:,MB_PARTNER) + gf(:,MB_DRAGIMP_PRED)
    rec(:,7) = gf(:,MB_RHS_OTHER)
    rec(:,8) = -gf(:,MB_DEFECT)
    rec(:,9) = -inertia_f
    rec(:,10) = pp(:,2) + pp(:,5) + ga(:,MB_FEEDBACK)
    hydro_dpa = hydro - dpa

    CALL WRITE_LINES(6)
    CALL WRITE_LINES(mfile)

  CONTAINS

    SUBROUTINE WRITE_LINES(unit)
      INTEGER, INTENT(IN) :: unit
      CHARACTER(LEN=*), PARAMETER :: V = '3ES24.15'
      CHARACTER(LEN=64) :: head
      INTEGER :: j

      WRITE(head,'(A,ES24.16,A,I0)') ' t= ', time, ' step= ', istep
      WRITE(unit,'(A,A,A,16(A,'//V//'),A,ES12.4)') &
        'EL_MOM_BALANCE_P', TRIM(head), ' unit= N', &
        ' drag_expl= ', pp(:,1), ' drag_impl= ', pp(:,2), ' pressure= ', pp(:,3), &
        ' viscous= ', pp(:,4), ' lift= ', pp(:,5), ' weight= ', pp(:,6), ' hydro= ', hydro, &
        ' fb_used= ', pp(:,7), ' deposit= ', pp(:,8), ' dp_dt= ', dp_dt, &
        ' expected_dt= ', expected_dt, ' floor_bal= ', floor_bal, &
        ' wall_contact_raw= ', wc, ' wall_lub_raw= ', wl, ' res_p_raw= ', res_p_raw, &
        ' res_p= ', res_p, ' relax= ', PE_WALL_RELAXATION
      WRITE(unit,'(A,A,A,23(A,'//V//'))') 'EL_MOM_BALANCE_GA', TRIM(head), &
        ' unit= N(mom_*:kg_m/s)', (' '//TRIM(gas_names(j))//'= ', ga(:,j), j = 1, MB_NG), &
        ' hrz_pred= ', hrz_pred(:,1), ' hrz_end= ', hrz_end(:,1)
      WRITE(unit,'(A,A,A,23(A,'//V//'))') 'EL_MOM_BALANCE_GD', TRIM(head), &
        ' unit= N(mom_*:kg_m/s)', (' '//TRIM(gas_names(j))//'= ', gd(:,j), j = 1, MB_NG), &
        ' hrz_pred= ', hrz_pred(:,2), ' hrz_end= ', hrz_end(:,2)
      WRITE(unit,'(A,A,A,2ES24.15,A,2ES24.15,5(A,'//V//'),A,2ES24.15,A,'//V//')') &
        'EL_MOM_BALANCE_B', TRIM(head), &
        ' p_in_out= ', p_in, p_out, ' area_in_out= ', pl(1), pl(3), ' dpA= ', dpa, &
        ' flux_in= ', pl(5:7), ' flux_out= ', pl(8:10), ' trac_in= ', pl(11:13), &
        ' trac_out= ', pl(14:16), ' eps_in_out= ', eps_in, eps_out, ' fbD_class= ', fbdc
      WRITE(unit,'(A,A,19(A,'//V//'))') 'EL_MOM_BALANCE_C', TRIM(head), &
        ' res_gas_F= ', res_gas, ' defect_F= ', gf(:,MB_DEFECT), ' gas_roundoff= ', &
        res_gas + gf(:,MB_DEFECT), ' proj_F= ', proj_f, ' proj_D= ', gd(:,MB_INERTIA_PROJ), &
        ' jd= ', jd, ' jd_check= ', jdc, ' hydro_dpA= ', hydro_dpa, ' rec_sum= ', SUM(rec, DIM=2), &
        ' rec_pB= ', rec(:,1), ' rec_pfp= ', rec(:,2), ' rec_visc= ', rec(:,3), &
        ' rec_conv= ', rec(:,4), ' rec_fbD= ', rec(:,5), ' rec_hrz= ', rec(:,6), &
        ' rec_other= ', rec(:,7), ' rec_solver= ', rec(:,8), ' rec_inert= ', rec(:,9), &
        ' rec_dep= ', rec(:,10)
    END SUBROUTINE WRITE_LINES

  END SUBROUTINE EL_MOMBAL_WRITE

END MODULE EL_MOMENTUM_BALANCE
