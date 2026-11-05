# Conservative continuity form: design consultation (GPT-6 Astra high, 2026-11-04)

Consultation run before any code was written (plan mode), on the branch head fc4fa9c2. The brief
and the answer are reproduced verbatim; the design note el_va_conservative_form_2026-11-05.md
states which recommendations were adopted.

## Brief

```
# Design consultation: conservative form of the EL volume-averaged continuity constraint ("form C")

You are consulted on a DESIGN, before any code is written. Read the code under the working
directory (FeatFloWer, branch feature/el-volume-averaged, HEAD fc4fa9c2) where you need to verify a
claim. Do not modify any files. Fortran 90/2008 (gfortran 14.3), CMake/ninja, MPI (OpenMPI 4.1.6).
Your answer is read by the owner, who decides what gets implemented.

## 1. Problem

FeatFloWer is a Q2/P1 (triquadratic velocity, discontinuous P1 pressure: per element the constant
mode plus physical x-xc, y-yc, z-zc, xc = the Q2 element-centre node) Navier-Stokes solver with MPI
domain decomposition (rank 0 = master without mesh; workers own element subdomains, Q2 nodes on
partition interfaces are duplicated and summed with E013Sum*). The Euler-Lagrange (EL) branch
deposits a particle volume fraction alpha_p per element (kernel, source/src_el/el_quadrature.f90:489-605),
eps = 1 - alpha_p per element only (no vertex/face values; no ghost elements), its backward rate
deps_f_dt per element, and the cell-averaged ANALYTIC kernel gradient grad_alpha_p(3,nel).

The volume-averaged continuity constraint is enforced in the projection as div u = g per element:
source/src_el/el_volume_averaged.f90, EL_ADD_CONTINUITY_SOURCE (lines 204-367):
  u8 = mean of the 8 corner values of the Q2 velocity (lines 262-270)
  eps = MAX(epsilon_f(iel),1d-3); d_cont = BE or BDF2 rate (EL_VA_CONTINUITY_RATE, 139-146)
  adv = -DOT_PRODUCT(u8, grad_alpha_p(:,iel))          ! = u . grad eps, lines 274-275
  g   = -(d_cont + adv)/eps                             ! line 276
  ... reductions (one 8-slot MPI_Allreduce on MPI_COMM_SUBS, always; a 2-slot one on audit steps;
      a 2-slot one in weighted sink modes) ...
  g_corr = gk(iel) - EL_VA_SINK(mode, gbar, int_g, w_e, w_total)   ! uniform: gbar; alpha: int_g*alpha_e/W
  defp(4*(iel-1)+1) = defp(4*(iel-1)+1) + g_corr*vol/dt            ! constant P1 row only, vol = dvol(iel)
Only the constant mode receives the source; the linear rows keep the unweighted
int_e (x_m - xc_m) div u = 0. defp enters the routine as assembled by Matdef_General_LinScalar
(source/src_quadLS/QuadSc_def.f90:1912-1943 -> BT_Mul_U_mod, QuadSc_proj.f:727-739):
  defP = -(1/dt) * B^T u   (B assembled in source/assemblies/QuadSc_BMatrix.f with CB3H(9) = 3x3x3 Gauss
  and the isoparametric 27-node Q2 Jacobian, Transform%ILINT = 2), so at entry
  -dt*defp(4*(iel-1)+m) = int_e q_m div u_h  for the four P1 test functions q_m.
The pressure solve then yields u^{k+1} with int_e q_m div u^{k+1} = -dt*defp_m (to solver tolerance),
Velocity_Correction/Pressure_Correction in source/src_quadLS/QuadSc_corrections.f90.

Deck-switchable extensions already on the branch (all default off, bitwise at defaults, reviewed):
  SimPar@ELVASourceTimeScheme = BE|BDF2 ; SimPar@ELVASourceCorrectors = 0..20 ;
  SimPar@ELVACompatibility = uniform|alpha   (el_config.f90:89-105, name table el_va_compat_names,
  validation 288-308, echo 349-351/383-385; parser cases applications/q2p1_el_pipeflow/app_init.f90:803-808)
The corrector loop EL_VA_Projection_Correctors (source/src_quadLS/QuadSc_el_va.f90:50-130, an include of
Transport_Q2P1 in QuadSc_main.f90) re-assembles defP with the corrected velocity, calls the source
with icorr=1..n, re-solves the pressure, corrects u and p; a final assembly-only pass (-1) on audit
steps logs the residual (EL_VA_LOG_CONSTRAINT_RESIDUAL, el_volume_averaged.f90:372-404). Predictor call:
QuadSc_main.f90:800-804 with icorr=0; corrector call site 832-836 guarded by
EL_VA_SOURCE_CORRECTORS_ENABLED() (deck-only, so the master joins every collective solve).

Measured defect (s04 voidfraction twin against CFDEM, VB return-flow case, closed box, cloud of
particles settling): with BDF2 + 1 corrector + alpha sink, the eps-weighted return flow G5 =
int eps u_z / (N V_p |w|) is 1.0055 at level 2 and 1.0021 at level 3, dt-independent. The analysis
(the peer's spec and an earlier Codex report) attributes the remainder to the discrete product-rule
defect: the net raw source int_g = sum_e eps_e... precisely R = sum_e eps_e int_e div u_h +
sum_e |e| u8_e . grad eps_e is the discrete form of int div(eps u) dV, which is exactly 0 in a closed
box but not discretely (u8 . grad eps_e is not a conservative flux). The "conservative form C" was
ranked as the durable solution: enforce div(eps u) = -d eps/dt with a SHARED face flux so that
D_eps(1) = 0 (telescoping) and the z-moment identity D_eps(z) = -int eps u_z hold exactly, and the
suggested staged route was: keep the existing projection as the preconditioner and iterate the
conservative residual G(U^k) = E^{-1}(f - L_eps U^k), L_eps = D_eps - E D, E = diag(eps_e) on the
four pressure rows, through the existing corrector loop; algebraically form C at convergence.

## 2. Proposed design (what I want you to critique)

### 2.1 Discrete operator
Nodal eps: build eps_h in the Q2 space, eps_q2(ndof), nodal value = volume-weighted mean of the eps
of the elements adjacent to that Q2 node (all 27 dofs per element: vertices, edges, faces, centre):
  num(kdfg(i)) += eps_e*V_e ; den(kdfg(i)) += V_e  over local elements, then E013Sum2(num,den)
  (the certified cross-rank Q2 nodal sum used for the velocity, needs ILEV = NLMAX pinned),
  eps_q2 = num/den.  Bounded (convex combination), continuous across ranks by construction.
  Alternative considered: the Q1 trilinear interpolant of vertex averages expressed in Q2 dofs
  (edge = mean of 2 vertices etc.), which is lower degree (Q1xQ2xcofactor = degree 5 -> 3-point
  Gauss exact on trilinear cells).
Element flux moments (face quadrature, template ExtractBoundaryNormals_Q2 in
source/src_quadLS/QuadSc_force.f90:620: SetUpMyCub ICUB=3 (3x3 Gauss per face), E013A, 27-node
Jacobian, SurfDet, normal oriented outward by the element centre node):
  F_e^(1)(u) = sum_{6 faces} int_f eps_h u_h . n dS
  F_e^(m)(u) = sum_f int_f (x_m - xc_m) eps_h u_h . n dS - int_e eps_h u_{h,m} dV   (m = 2..4, the
               conservative weak form for the linear test functions; the volume term by CB3H(9))
  plus for diagnostics int_e eps_h u_h dV (the eps_h-weighted velocity integral) and
  Q_m = int_e (x_m - xc_m) dV.
Both elements adjacent to an interior face evaluate the SAME nodal data at the SAME 9 face points
with opposite normals, so sum_e F_e^(1) = boundary flux exactly (to rounding) on ANY mesh, with no
neighbour-element lookup and no face exchange (kadj is 0 on both physical and partition faces here
and is not needed). On no-slip walls all 9 face nodes are Dirichlet, so the wall flux is exactly 0.

### 2.2 Fixed point through the existing projection
Per source evaluation with the current velocity u^k (predictor: the momentum-predictor velocity;
correctors: the corrected velocity), per element and pressure row m:
  D_e^(m)(u^k) = -dt*defp_pre(4*(iel-1)+m)        ! B's own quadrature, read at routine entry
  J_e^(m) = F_e^(m)(u^k) - eps_e * D_e^(m)(u^k)   ! the conservative "advective" part
  target: F_e^(m)(u^{k+1}) = -d_e * Q_m            ! Q_1 = V_e (= dvol, as now)
  linearised with D_e^(m)(u^{k+1}) - D_e^(m)(u^k) treated implicitly, J explicit:
  g_e^(m) := -(d_e*Q_m + J_e^(m))/(eps_e*V_e)   ;  defp_m += g_e^(m)*V_e/dt  (m=1 additionally
  minus the sink s_e, exactly as now).
For m = 1 this is literally the current formula g = -(d_cont + adv)/eps with
  adv := J_e^(1)/V_e = (1/V_e)[ F_e^(1)(u^k) - eps_e int_e div u^k ]
i.e. the advective term in Green-Gauss (conservative) form instead of u8 . grad eps_e. At the fixed
point F_e^(m) = -d_e Q_m - eps_e V_e s_e, summing the m=1 rows gives sum_e eps_e V_e s_e = 0, so the
compatibility sink vanishes at convergence (uniform: gbar = 0; alpha: int_g = 0) and the local
conservative balance holds exactly on any mesh because D_e is read from B's own quadrature (no
mismatch between the evaluator's and B's volume rule). Contraction factor ~ ||E^{-1} L_eps|| ~
h|grad eps|/eps: small for dilute clouds (1-2 correctors), larger for dense beds.
With the linear rows included, sum_e [zc F_e^(1) + F_e^(z)] = -int eps_h u_z exactly, hence at the
fixed point int eps_h u_z = sum_e zc V_e d_e + sum_e d_e Q_z = the deposited first-moment rate
(= N V_p |w| by the G3 identity; Q_z = 0 on bricks).

### 2.3 Code structure (sustainability is an explicit owner requirement)
- NEW module source/src_el/el_va_flux.f90 (EL_VA_FLUX): (a) EL_VA_ACCUMULATE_NODAL_EPS(num, den,
  mesh, ilev, epsilon_f) pure local accumulation, the caller does the E013Sum2 and the division (keeps
  the module MPI-free and unit-testable); (b) EL_VA_ELEMENT_FLUX_MOMENTS(iel, ...) -> flux(4),
  ueps(3), qmom(3): pure per-element face+volume quadrature; (c) a driver over all local elements
  that does the E013A/cubature setup once. No COMMON-block state beyond what the existing EL
  quadrature routines already use (/ELEM/ /CUB/ /TRIAD/ ...).
- el_fields.f90: epsilon_q2(:) (ndof), allocated only under the conservative form, rebuilt in the
  post-advance transfer pass right after EL_FINALIZE_VOID_FRACTION (eps^{n+1} final there, ILEV
  pinned at that call site), not in the restart file.
- el_config.f90: key SimPar@ELVAContinuityForm = advective | conservative (default advective),
  name table el_va_form_names like el_va_compat_names; validation: requires ELVolumeAveraged =
  continuity and ELVASourceCorrectors >= 1 (without a corrector it is not form C); echo lines.
  Possibly a second table entry 'conservative_p0' = constant row only (staging/diagnosis)?
- el_volume_averaged.f90: form id resolved once per call; in the element loop a SELECT CASE on the
  form for adv (default branch keeps the current statements verbatim -> bitwise at the default);
  the linear-row targets and the extra reductions only under the conservative form (no new
  collective at the default; collectives keyed on deck values only, never on rank-local state);
  a new audit line EL_VA_FLUX (sum_e F_e^(1) telescoping check, L1/max of F_e^(1)+V_e d_e,
  z-identity residual, int eps_h u_z / V as the like-quantity superficial velocity) with the
  existing corr= tag; EL_VA_LOG_CONSTRAINT_RESIDUAL extended to all four rows under the form.
- Call sites: ILEV = NLMAX; CALL SETLEV(2) before the source call under the conservative form
  (Matdef sets ILEV but not the /TRIAD/ common that NDFGL/E013 read); the predictor call in
  QuadSc_main.f90 and the corrector include.
- Unit test source/src_el/tests/test_el_va_flux.f90 on el_test_helpers' BUILD_STRUCTURED_MESH
  (unit cube, kvert/kedge/karea/dcorvg/dvol, /TRIAD/ set; no kadj, none needed): telescoping to
  1e-14 for a manufactured velocity with u.n = 0 on the walls and a random nodal eps; constant eps
  and linear velocity -> F_e^(1) = eps*V_e*div u exactly; the z-identity to 1e-14; orientation.
- Validation ladder on s04 VB (level 2 and 3, dt and corrector-count series), gates: every
  certified configuration bitwise at the default (tier-2 twin, VB_L2, T9a keys, s01, semi-implicit,
  P4/P5 controls) via tools/el_stream_compare.py; then the usual two reviews.

## 3. Questions (answer each; say [V] verified in code, [I] inferred, [P] proposal)
1. Is the fixed-point formulation in 2.2 correct, including the use of -dt*defp_pre as D_e, the
   sign/scale conventions of defp and the sink interaction? Any term missing (e.g. the pressure
   extrapolation P_new, GAMMA>0 rotational branch, Boundary_LinScalar_Def on rhsP)?
2. Nodal eps representation: Q2-nodal volume-weighted average vs Q1 trilinear interpolant vs P0
   face mean (would need the per-face exchange mg_mpi(ILEV)%parST / CommSumHalf in
   source/src_mpi/pp3d_mpi.f90:849-1195). Which do you recommend for conservation, consistency with
   the P0 rate d_e in the storage term, boundedness, and the diagnostics ("like quantities")?
3. Face form vs volume form int_e div(eps_h u_h) dV (same 3x3x3 rule as B): the volume form is
   simpler (template EL_INTEGRATE_FLUID_MOMENTUM, el_quadrature.f90:158-242) but telescopes exactly
   only on affine cells. Is the face form worth its bookkeeping?
4. Linear rows: include them from the start (full four-row residual) or stage the constant row
   first? Is there an inf-sup / solvability concern in modifying the linear-row targets through the
   existing C = B^T M_l^-1 B solve (unsymmetric solvers are used: SOR smoother, UMFPACK coarse)?
5. Convergence of the outer iteration: is a residual-based stopping key (e.g.
   ELVASourceCorrectorTol) needed now, or is the fixed corrector count enough for the ladder?
6. MPI safety: any collective or rank-local early return in the proposed flow that could desync
   the master/workers or the workers among themselves? (Existing early returns at
   el_volume_averaged.f90:232-240 precede the collectives.)
7. Code structure: is the module split, the key/table design, and the SELECT CASE strategy slot
   clean, maintainable and extensible (e.g. for a later full-operator replacement B_eps in the
   pressure matrix, or a fourth time scheme)? Name anything that would make a future maintainer
   trip. Propose better names if ours are poor.
8. Performance: cost per source evaluation (6 faces x 9 points x 27 basis functions per element
   plus the 27-point volume rule) relative to the deposition; anything to precompute per step?
9. What would you test first to falsify the design cheaply (before the cluster ladder)?
10. Anything else: risks, the dvol single-precision storage (SETARE stores REAL(AAA)), the
    el_eps_f_min clip and eps_f_relax interplay, open boundaries (inflow/outflow, periodic E1 bed),
    restart.

Output: findings ordered by severity (file:line where applicable, what is wrong, why it matters,
the fix), then your answers to the ten questions, then a one-paragraph verdict with a recommended
staging (what to implement first). Do not modify any files.
```

## Answer (Codex gpt-6-astra, reasoning effort high)

The fixed-point idea is sound for the closed-box problem, but several guarantees in the proposal are too strong. I recommend proceeding with **shared face quadrature, Q1 porosity reconstructed into Q2 DOFs, and all four pressure rows**, with explicit residual and compatibility diagnostics.

[V] Reviewed HEAD `fc4fa9c2958423061f751b5348a7ad1978cb299c`. No files were modified or simulations run. Below, **[V]** means verified in this checkout, **[I]** means derived or inferred, and **[P]** means recommended design.

Findings, ordered by severity:

1. **High — a vanishing sink requires compatible storage, not merely convergence.**

   [I] Write \(S=\sum_e d_eV_e\), and let \(\Phi\) be the boundary flux of \(\epsilon_hu_h\). Your fixed-point equation implies
   \[
   \sum_e\epsilon_eV_es_e=-(S+\Phi).
   \]
   Therefore the uniform/alpha sink vanishes only when \(S+\Phi=0\), with positive effective porosity and a nondegenerate sink weight.

   [V] Clipping and relaxation modify the stored porosity and its rate in [el_transfer.f90:523](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_transfer.f90:523). Clipping can destroy constant total fluid volume even when particle volume is constant.

   [P] Audit \(S\), \(\Phi\), the weighted sink integral, and the **unsinked** conservative residual separately. A small sink-corrected projection residual must not certify conservation. For initial closed-box certification, require negligible clipping and negligible \(S\).

2. **High — the exact moment identity controls a different observable from the existing G5.**

   [V] The current momentum integration multiplies velocity by **elementwise** `epsilon_f` and uses `EL_Q1_MAP`: [el_quadrature.f90:223](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_quadrature.f90:223). The proposed identity controls the new quadrature integral of \(\epsilon_hu_z\).

   [I] Even at perfect convergence,
   \[
   \int \epsilon_hu_z
   \ne \sum_e\epsilon_e\int_eu_z
   \]
   in general. A new G5 equal to one does not prove that the old G5 defect disappeared.

   Nor is the deposited first-moment rate identically \(NV_p|w|\): the [s04 runbook:25](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/applications/q2p1_el_pipeflow/validation_cases/s04_voidfraction_cfdem_twin/RUNBOOK.md:25) reports measurable G3 scatter.

   [P] Preserve the old diagnostic and add the reconstructed-porosity diagnostic, their difference, and the actual BE/BDF2 storage-moment rate. Separate projection error, reconstruction difference, and deposition/time-discretization error.

3. **High — Q2 nodal averaging does not produce a bounded Q2 field.**

   [I] The nodal averages are convex combinations; their Q2 interpolant is not. A realizable one-dimensional example has cell porosities \(0.4,1\). In the first cell, endpoint/centre/endpoint values can be \(0.4,0.4,0.7\). Their quadratic interpolant reaches \(0.3625\), below the prescribed minimum. Reversing the contrast produces an overshoot above one.

   [P] Prefer the Q1 vertex reconstruction for the initial implementation. It preserves bounds throughout the reference cell. Clipping Q2 values at quadrature points would introduce another operator definition and invalidate the polynomial exactness arguments.

4. **High — the proposed transfer insertion point does not have global `ILEV` pinned. Existing local returns can also deadlock.**

   [V] At [el_transfer.f90:337](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_transfer.f90:337), global `ILEV` is pinned for `E013Sum3`, then restored **before** `EL_FINALIZE_VOID_FRACTION`. `/TRIAD/` remains pinned there, but the communication level does not.

   [V] The returns at [el_volume_averaged.f90:232](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:232) depend on worker-local allocation and sizes before worker collectives.

   [P] Explicitly pin and restore global `ILEV` around the new nodal sum. In conservative mode, validate worker readiness collectively; empty workers contribute zero, and inconsistent allocation/size failures terminate coherently. Do not add another allocation-dependent return.

5. **High for general deployment — pressure boundary rows and periodic moments invalidate the unconditional claims.**

   [V] The pressure solver applies `Boundary_LinScalar_Def` to `rhsP` at [QuadSc_def.f90:2006](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_def.f90:2006). That routine zeros **all four rows** of elements with `KNPRP=1`: [QuadSc_boundary.f90:688](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_boundary.f90:688). Those constraints are not enforced by your source update.

   [I] For open boundaries, the moment identity contains \(\int_{\partial\Omega}z\epsilon_hu\cdot n\). In a vertically periodic domain, \(z\) is not a periodic test function: opposite-face mass fluxes cancel, but their \(z\)-weighted fluxes do not.

   [P] First certify the closed, impermeable box. Treat open-boundary compatibility and periodic first moments explicitly before claiming support for them.

6. **Medium — “any mesh” conflates algebraic cancellation with geometry and quadrature consistency.**

   [V] `SETARE` stores a double-precision array entry after `REAL(AAA)` conversion and computes geometry from eight vertices: [mesh_refine.f90:3348](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_mesh/mesh_refine.f90:3348), [mesh_refine.f90:3420](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_mesh/mesh_refine.f90:3420). Deposition divides by this `dvol`, while using a Q1 map.

   [I] Consequently, `Q1=dvol` and quadrature-derived linear `Q_m` need not describe one consistent integration measure. Reading \(D\) from the matrix avoids a preconditioner mismatch; it does not repair the storage geometry.

   [P] Distinguish `storage_volume` from `geometry_volume`. Initially retain the deposition-compatible storage volume and quantify their discrepancy. Resolve the storage moments consistently before extending certification to curved geometry.

7. **Medium — the face template contains sign and factor-of-four traps.**

   [V] In [QuadSc_force.f90:725](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_force.f90:725), `dLine` points from face centre toward element centre. The sign test makes `dNorm` point inward; the later subtraction produces the outward contribution.

   [V] `SurfDet` returns four times the tangent-cross-product magnitude at [QuadSc_surftens.f:951](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/assemblies/QuadSc_surftens.f:951), paired with face weights normalized to sum to one in [QuadSc_stress.f:1443](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/assemblies/QuadSc_stress.f:1443).

   [P] Prefer an explicitly signed reference-face cofactor vector and one documented weight convention. Use the centre test as a validation check, rather than the primary orientation rule for curved cells.

8. **Medium — the proposed harness is insufficient as written.**

   [V] `BUILD_STRUCTURED_MESH` allocates `dcorvg(3,nvt)`, not coordinates for all Q2 DOFs: [el_test_helpers.f90:89](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/tests/el_test_helpers.f90:89). `BUILD_DOF_COORDINATES` supplies the missing coordinates separately.

   [V] The existing transfer-test documentation explicitly excludes real shared-DOF communication coverage: [test_el_transfer.f90:20](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/tests/test_el_transfer.f90:20).

   [P] Extend the geometry fixture and add a tiny production-initialized partition test. Running the synthetic test on several ranks does not test `E013Sum2`.

Answers to the ten questions:

**1. Fixed-point formulation, signs, pressure terms, and sink**

[V] The entry interpretation is correct for the assembled divergence operator:

- `Build_BTMatP1` assembles positive \(\int q\,\partial_j\phi_i\).
- `BT_Mul_U_mod` multiplies by \(-1/\Delta t\).
- `Velocity_Correction` consequently gives \(D\Delta u=\Delta t\,\mathrm{rhsP}\), assuming the pressure solve and boundary restrictions are consistent.

See [QuadSc_BMatrix.f:193](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/assemblies/QuadSc_BMatrix.f:193), [QuadSc_proj.f:737](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_proj.f:737), and [QuadSc_corrections.f90:24](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_corrections.f90:24).

[I] Define
\[
r_{em}=F_{em}(u^k)+d_eQ_{em}.
\]
Your update simplifies exactly, algebraically, to
\[
\begin{aligned}
\mathrm{defp}^{post}_{em}
&=-\frac{D_{em}}{\Delta t}
-\frac{d_eQ_{em}+F_{em}-\epsilon_eD_{em}}{\epsilon_e\Delta t}
-\delta_{m1}\frac{V_es_e}{\Delta t}\\
&=-\frac{r_{em}}{\epsilon_e\Delta t}
-\delta_{m1}\frac{V_es_e}{\Delta t}.
\end{aligned}
\]
Thus the sign and scaling are correct. At convergence, \(r_{em}=-\delta_{m1}\epsilon_eV_es_e\).

[P] Implement the conservative branch as this direct residual RHS. Preserve the pre-source divergence where needed for diagnostics and the existing sink calculation, but avoid forming \(J\), dividing by \(V\), multiplying by \(V\), and cancelling \(D\) numerically.

[I] The statement that the corrected divergence equals `-dt*defp_post` is incorrect: that array represents the correction RHS. The corrected divergence is \(D(u^k)+\Delta t\,\mathrm{defp}^{post}\).

[V] `P_new` is rebuilt from final pressure and restored \(P^n\) in [QuadSc_el_va.f90:122](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_el_va.f90:122). No extra extrapolated-pressure term belongs in the continuity residual.

[V/I] `GAMMA>0` modifies pressure **after** velocity correction, using `rhsP`, in [QuadSc_corrections.f90:48](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_corrections.f90:48). It does not alter that pass’s velocity constraint. Its pressure interpretation under repeated conservative corrections needs a separate test; begin certification with `GAMMA=0`.

[I] The constant-row sink is an algebraic compatibility correction. On non-centroidal cells it is not the four-row weak representation of a spatially constant scalar sink. Keeping it constant-row-only is reasonable for this iteration, provided that distinction is documented.

**2. Porosity representation**

[P] My preference is:

| Representation | Assessment |
|---|---|
| Q1 vertex averages, embedded in Q2 DOFs | Best initial choice: continuous, bounded, simple communication, favorable quadrature |
| Independent averages at all Q2 DOFs | Continuous but not pointwise bounded; extra variation without demonstrated benefit |
| Shared P0 face mean | Good conservative mass flux and close connection to P0 storage; requires explicit face communication and a separate volume-flux definition |

[I] None of the nodal reconstructions generally preserves every deposited cell average. That is acceptable when P0 porosity is the conserved storage variable and \(\epsilon_h\) is a reconstruction used for transport. It should be described that way, rather than as a single continuous porosity field satisfying both storage and flux equations.

[P] With a face-only porosity, the four-row formula still needs a specified volume term. For example, using \(\epsilon_eu\) there would make the moment identity control the P0-weighted integral. That is a legitimate alternative design, but different from a single reconstructed flux field.

**3. Face versus volume form**

[P] Use the face form for the durable implementation. Its main benefit is cancellation of shared interface contributions independently of integration accuracy.

[I] The “volume form telescopes only on affine cells” claim is too restrictive:

- With Q1 porosity, Q2 velocity, and trilinear geometry, three-point tensor Gauss quadrature is sufficient for all four weak moments.
- With Q2 porosity and trilinear geometry, the constant-row divergence integral is also exactly integrated by that rule; the linear moments need not be.
- General Q2 geometry requires higher orders for exact integration.

These degree counts use the standard \(2n-1\) exactness of \(n\)-point Gaussian quadrature; see [NIST DLMF §3.5(v)](https://dlmf.nist.gov/3.5#v).

[I] With face quadrature, telescoping and the combined moment identity remain algebraically valid for conforming shared geometry, traces, point sets, weights, and opposite orientations—even when the integrals are underintegrated. Call the result a **quadrature identity** on curved cells. It is not automatically an exact physical integral.

**4. Linear rows and solvability**

[P] Include all four rows in the first solver implementation. A constant-row prototype is useful inside tests, but I would not add a permanent `conservative_p0` deck option yet. Constant-row conservation alone cannot establish the requested return-flow identity.

[I] Changing the RHS does not change the existing matrix’s inf-sup properties. But convergence requires more than solvability of that matrix.

Let \(H\) be the boundary-constrained velocity correction map, so \(C=DH\). The effective conservative pressure operator is
\[
A_\epsilon=F H.
\]
Its rank and conditioning on the pressure quotient space matter. It need not be symmetric, and using SOR/UMFPACK for the existing \(C\) solve does not establish convergence for \(A_\epsilon\).

[P] On a small mesh, inspect the nullity and smallest nonzero singular values of both operators. A very small no-slip mesh can have insufficient velocity freedom to support all pressure modes, so do not infer general stability from a one-element projection test.

**5. Outer convergence and stopping**

[I] The proposed contraction estimate is a heuristic, not an operator bound. It omits the pressure inverse and the velocity correction map. With porosity frozen, this is a stationary linear iteration. Schematically its residual propagation is
\[
r^{k+1}
=\left[I-FH\,C^\dagger\Pi E^{-1}\right]r^k,
\]
where \(\Pi\) denotes compatibility projection.

[P] Fixed counts are sufficient for the first ladder if every run records the final **unsinked four-row residual**, demonstrates convergence with increasing counts, and meets a specified threshold. One corrector is neither necessary nor sufficient mathematically to reach form C.

A tolerance key can follow the first convergence study. Before general production use, add absolute/relative tolerances, a maximum count, and a reported nonconvergence outcome. Normalize linear residuals by a cell length before combining them with constant rows; their units differ.

**6. MPI safety**

[V/P] The safe separation is:

- Nodal reconstruction and its `E013Sum2`: all workers, including workers without particles.
- Flux/residual reductions: all workers.
- Pressure solves and any loop-termination decision: all solver ranks, including master.

If adaptive stopping is introduced, communicate the worker decision to master before any rank exits the corrector loop.

[V] `E013Sum2` dispatches only communication modes 3 and 4 at [QuadSc_mpi.f90:2926](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_mpi.f90:2926); mode 3 is the declared default. Validate that supported mode explicitly.

[P] Pin `ILEV` and `/TRIAD/` at predictor, corrector, and final diagnostic evaluations. Do not sum the already normalized porosity again. Test face-, edge-, and corner-sharing, plus periodic identifications; ordinary partition sharing alone does not establish periodic correctness.

**7. Structure, keys, and maintainability**

[P] The module split is appropriate, with three clear responsibilities:

- reconstruct the transport porosity;
- evaluate conservative flux moments;
- construct and project the residual through the existing solver.

`EL_VA_FLUX` and `ELVAContinuityForm` are reasonable names. I would use `epsilon_transport_q2`, `flux_moments`, `storage_moments`, and `conservative_residual` to distinguish reconstruction, units, and purpose.

Prefer one form dispatch outside the element loop, with the legacy path left intact. The conservative path is now sufficiently different—four rows, geometry, residual semantics—that treating it only as another `adv` formula obscures the design.

[V/I] A routine calling legacy `NDFGL`/element machinery and mutating COMMON state cannot be a Fortran `PURE` procedure. For actual purity, pass connectivity, coordinates, and precomputed reference tables into a stateless numerical kernel; keep COMMON interaction in a wrapper.

[P] Other maintenance requirements:

- Separate optional reconstruction allocation from `EL_ENSURE_FIELDS`, whose size-based early return can bypass allocation.
- Invalidate caches after mesh changes, repartition, restart, or changed porosity.
- Add the source to both relevant lists in [ProjectFiles.cmake:295](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/cmake/modules/ProjectFiles.cmake:295).
- Have time integration return element rates independently of the spatial form.
- Define a later full operator as \(F H\), unless momentum correction is also deliberately changed; simply naming it `B_eps` leaves that choice ambiguous.

**8. Performance**

[I] The straightforward evaluator performs \(81\times27=2187\) basis contributions per element, plus geometry and moment work. It is much more expensive than the present eight-corner source, but no defensible deposition-cost ratio follows from that count: deposition depends on particles, support overlap, and repeated candidate-element searches.

[P] Precompute:

- once per mesh: reference tables, element maps, geometry, face area vectors, pressure-test values, storage moments, nodal averaging denominators;
- once per porosity update: reconstructed porosity and its quadrature values;
- optionally once per step: local coefficients mapping velocity DOFs directly to flux moments and weighted velocity integrals.

Only nine Q2 basis traces are nonzero on a face; exploit that after establishing correctness. Time reconstruction, flux evaluation, reductions, and pressure solves separately before deciding whether coefficient caches justify their memory cost.

**9. First cheap falsification tests**

[P] I would run these before the cluster ladder:

1. **One-cell geometry test:** integrate all six face areas and verify outward flux for \(u=(x,y,z)\). This catches the template’s sign and factor-of-four errors.
2. **Two-cell shared-face test:** use different valid local face orientations, unequal porosity, and arbitrary continuous Q2 velocity. Check pairwise cancellation before global sums.
3. **Four-row identity test:** verify
   \[
   \sum_e(z_{c,e}F_{e1}+F_{ez})+I_z-\Phi_z=0.
   \]
   Include warped trilinear and curved Q2 geometry; compare quadrature orders separately from cancellation.
4. **Reconstruction test:** include a sharp porosity jump and test values inside cells, not just at nodes.
5. **Small assembled projection test:** verify the sign of \(D\Delta u\), then manufacture a reachable conservative target and measure iteration convergence for several porosity contrasts.
6. **Compatibility negative test:** prescribe nonzero total storage rate in a closed box. The audit must expose the incompatibility even if the sink-corrected solve converges.
7. **Tiny real MPI test:** compare one worker with several workers, including a partition without particles and a periodic seam.

Use scale-aware tolerances rather than a universal absolute `1e-14`. Test all three coordinate-moment identities.

**10. Other risks: volume, clipping, boundaries, restart**

[P] Keep the initial geometry/storage policy explicit. If retaining deposition’s `dvol`, one coherent convention is
\[
Q_{em}=V_e^{storage}
\frac{\int_e q_m\,dV}{V_e^{geometry}},
\]
so all storage moments use the same discrete measure. This is a design choice, not identical to integrating the original \(d_e\) against geometric volume. A later change to double-precision geometry volumes must update deposition and storage consistently.

[I] The `1e-3` floor used in \(E\) can remain a positive preconditioning coefficient if it is used consistently throughout the residual construction. It need not define the transported porosity. Avoid conflating that floor with `el_eps_f_min`.

[V/I] Relaxation makes the rate describe relaxed porosity, so its first moment need not equal the instantaneous particle-volume flux. Clipping can produce genuine global incompatibility. Neither effect is repaired by conservative spatial fluxes.

[P] For open boundaries, derive compatibility from the actual prescribed/unknown fluxes and pressure-row treatment; do not automatically subtract the closed-box mean. For periodic E1, retain periodic mass conservation but include the seam contribution in coordinate-moment diagnostics.

[V] The mesoscale filter runs **after** the fluid core at [QuadSc_main.f90:435](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_main.f90:435). [I] A velocity filter can break the newly established conservative constraint. Audit after it, make it constraint-preserving, or reject that combination initially.

[V/P] Restart version 3 already stores porosity and rate history. Reconstructing transport porosity rather than storing it is appropriate. Rebuild it before use without advancing history again, and test continuous versus restarted runs.

[I] Finally, `int_g` is not literally the product-rule defect \(R\) quoted in the problem: it contains the divided storage and advective terms. Report \(R\), storage incompatibility, and `int_g` as distinct quantities when testing the causal explanation.

My verdict is to implement the face operator first as an independently testable four-row evaluator, using bounded Q1 porosity reconstruction and explicit geometry/storage conventions; then connect its conservative residual directly to the existing projection with fixed corrector counts. Certify a closed, unclipped, unfiltered box with `GAMMA=0`, measuring both old and new weighted-flow diagnostics and the unsinked residual. After the small geometry, algebra, and real-MPI tests pass, run the default-bitwise gates and s04 ladder. Add adaptive stopping and broader boundary/geometry support after the iteration’s convergence and the observable change have been demonstrated.