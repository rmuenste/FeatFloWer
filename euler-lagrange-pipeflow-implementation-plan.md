# Implementation Plan: Unresolved Euler–Lagrange (CFD-DEM) Module — Milestone 1: Laminar Pipe Flow, φ ≤ 20%

**Audience:** Codex/Claude Code, working inside the existing FEM-FBM DNS codebase (Q2/P1disc finite elements, fractional-step-θ time integration, multigrid solvers, coupled `pe` rigid-body physics engine, MPI domain decomposition shared between CFD and DEM).

**Goal:** Add an *unresolved* Euler–Lagrange simulation mode in which particles are smaller than the mesh resolution, tracked as Lagrangian point/parcel objects by `pe`, and coupled to the fluid via interpolated closure forces and conservative force spreading — instead of the resolved FBM rigid-body constraint and stress-integral force evaluation.

**First target case:** steady laminar pipe flow of a liquid–solid suspension, solid volume fraction up to 20%, neutrally buoyant and non-neutrally-buoyant particles.

---

## 0. Scope and non-goals

**In scope (this milestone):**
- New runtime mode `EULER_LAGRANGE` selectable per simulation, coexisting with the existing `DNS_FBM` mode. No regression to DNS behavior.
- Simplified Model-A formulation: fluid keeps `∇·u_f = 0`; void fraction ε_f enters only the drag closure and diagnostics (see §1 for exact equations and the upgrade hooks that must be built now).
- Forces on particles: hindered drag, generalized buoyancy −V_i ∇p, gravity, Saffman–Mei shear lift with wall correction, optional Magnus lift.
- Two-way coupling via conservative kernel spreading; semi-implicit drag treatment in the fluid momentum operator.
- Four-way coupling: reuse existing `pe` contact (PGS) and pairwise lubrication models, acting on particles only, with DEM sub-cycling.
- Validation suite: settling particle, Richardson–Zaki hindered settling, Segré–Silberberg migration, pressure drop vs φ against in-house DNS and μ_eff-based Hagen–Poiseuille.

**Out of scope (later milestones):**
- ε_f-weighted continuity/momentum (`∇·(ε_f u_f) = −∂ε_f/∂t`) — but the ε_f and ∂ε_f/∂t field assembly IS in scope so the upgrade is a switch, not a rewrite.
- Torque feedback to the fluid; added-mass and history forces; turbulence; heat transfer; dense regime φ > 0.25.

---

## 1. Model specification (authoritative — implement exactly this)

### 1.1 Fluid phase

Incompressible Navier–Stokes as currently discretized (Q2/P1disc, fractional-step-θ, streamline-diffusion stabilization), plus a particle momentum source:

```
ρ_f (∂u_f/∂t + u_f·∇u_f) = −∇p + ∇·τ_f + ρ_f g + f_{p→f}
∇·u_f = 0
```

with `τ_f = μ_f (∇u_f + ∇u_f^T)`. The fluid viscosity stays μ_f — do NOT use the effective-viscosity tables here (those are for Euler–Euler; using them together with hindered drag double-counts dissipation).

### 1.2 Particle phase (per particle i, integrated by `pe`)

```
m_i dU_i/dt = F_i^drag + F_i^∇p + m_i g + F_i^lift + F_i^col + F_i^lub
dX_i/dt     = U_i
I_i dω_i/dt + ω_i×(I_i ω_i) = T_i^col + T_i^lub   (+ optional hydrodynamic spin torque, see §5.4)
```

**Drag (Di Felice form):**
```
F_i^drag = (1/2) ρ_f C_D(Re_p) A_i |u_rel| u_rel · ε_f^{−χ}
u_rel    = u_f(X_i) − U_i
Re_p     = ρ_f ε_f d_i |u_rel| / μ_f
C_D      = (0.63 + 4.8/√Re_p)²
χ        = 3.7 − 0.65 exp(−(1.5 − log10 Re_p)²/2)
```
ε_f and u_f are kernel-interpolated at X_i (§3). Guard Re_p ≥ 1e−12. In the Stokes limit this reduces to hindered Stokes drag; add a unit test asserting `|F_drag − 3πμ_f d (u_f−U)| / |·| < 1%` for Re_p < 1e−3, ε_f = 1.

**Generalized buoyancy:** `F_i^∇p = −V_i ∇p(X_i)` with ∇p kernel-interpolated. Combined with `m_i g` this reproduces Archimedes buoyancy in hydrostatics AND the axial driving-gradient force in pipe flow. Do NOT additionally add `−ρ_f V_i g` (that would double count).

**Saffman–Mei lift with near-wall correction:**
```
F_i^lift = C_L ρ_f (π/8) d_i³ (u_rel × ω_fluid),   ω_fluid = ∇×u_f(X_i)
```
with C_L per Mei (1992) extension of Saffman, and a wall-distance correction factor of the Zeng et al. (2009) family. Implement the lift coefficient in an isolated, unit-tested function with the literature formulas cited in comments; expose `lift_model = {none, saffman_mei, saffman_mei_wall}` in the config. Lift is the physics that produces Segré–Silberberg migration — it must be switchable for the validation matrix.

**Magnus lift (optional, default off):** Rubinow–Keller form using particle spin relative to fluid rotation; particles acquire spin from `pe` contacts/lubrication.

**Contacts and lubrication:** unchanged `pe` PGS hard-contact solver and existing lubrication force/torque formulas (normal, sliding, torque terms with slip regularization), applied particle–particle and particle–wall. These act on particles ONLY — never spread to the fluid in this milestone.

### 1.3 Two-way coupling source

```
f_{p→f}(x) = − Σ_i (F_i^drag + F_i^lift) W_δ(x − X_i)
```

Note: `F_i^∇p` is NOT fed back (the fluid already carries −∇p — Model-A bookkeeping). Contacts/lubrication are NOT fed back. Document this convention prominently in the module header; mixing conventions is the classic CFD-DEM bug.

### 1.4 Void fraction

```
α_p(x) = Σ_i V_i W_δ(x − X_i),    ε_f = max(1 − α_p, ε_min),   ε_min = 0.4 (config)
```
At φ ≤ 0.2 the clip should never activate; emit a warning counter if it does.

---

## 2. Phase 0 — Reconnaissance (no code changes)

Codex/Claude Code should first produce a short written map (commit as `docs/el_codebase_map.md`) identifying, with file/routine names:

1. Where the fractional-step-θ substeps assemble the momentum RHS and system matrix (Q2 velocity space), and where Dirichlet/FBM constraints are imposed — the E-L mode must bypass the FBM nodal constraint path entirely.
2. Where FBM hydrodynamic forces `F_i = −∫ σ·∇α_i dx` are computed — the E-L force models replace this call site behind a common interface.
3. The CFD↔pe data exchange layer: how particle positions/velocities are sent to the CFD side, how forces are returned, how ghost/shadow particles are synchronized across MPI subdomains.
4. The element-search / point-location facility (the FBM must already classify nodes inside particles — find what exists for locating an arbitrary point X_i in the hexahedral mesh and evaluating Q2 basis functions there).
5. The time-step driver: where Δt is chosen, where pe is stepped, whether sub-cycling hooks exist.
6. How boundary conditions for the pipe geometry are configured (inflow/outflow vs periodic+body force), since the validation cases need a periodic pipe driven by a constant axial pressure gradient.

Acceptance: the map exists, names concrete routines, and lists any mismatch with assumptions in this plan. **If reality contradicts this plan, the map should say so and propose the adaptation — do not silently improvise.**

---

## 3. Phase 1 — Particle–mesh transfer infrastructure

New module (suggested name `el_transfer`), independent of force models.

### 3.1 Kernel

- Strictly positive, compactly supported kernel `W_δ(r)`, e.g. the Deen/Peskin-type polynomial or a truncated Gaussian. Width δ is a config parameter tied to particle size: default `δ = 2.5 d_p`, NOT tied to mesh h.
- **Wall renormalization:** for particles within δ of the pipe wall, renormalize so the discrete weights satisfy `Σ_K w_iK = 1` over the fluid domain. Implementation: compute raw weights over overlapped elements, divide by their sum. This is the main case in a pipe, not a corner case.
- Do NOT use Q2 shape functions `N_a(X_i)` as spreading/accumulation weights — they take negative values and will produce negative void fractions and oscillatory sources. Q2 basis evaluation is used only where genuinely interpolating FE fields (and even there, prefer kernel-weighted averages, §3.3).

### 3.2 Gather/scatter operations

Implement and unit-test, in parallel:
- `accumulate_alpha_p()` → cellwise/nodal α_p, ε_f fields. Conservation test: `Σ_K |K| α_{p,K} = Σ_i V_i` to machine precision (after wall renormalization), on 1, 2, and 8 MPI ranks, including particles whose kernel support straddles subdomain boundaries (requires contribution from ghost particles — reuse pe's ghost mechanism).
- `spread_force(F_i) → f_{p→f}` nodal RHS contributions in the Q2 velocity space via kernel-weighted projection: `b_a^p = −Σ_i F_i ∫ N_a W_δ(x−X_i) dx` (quadrature on overlapped elements). Conservation test: `Σ_a b_a^p = −Σ_i F_i` componentwise.
- `interpolate_at_particle()` for u_f, ∇u_f, p, ∇p, ε_f at X_i, kernel-weighted. Test: exact reproduction of constant fields; O(δ²) error on linear fields away from walls.
- Temporal smoothing hook for ε_f (exponential moving average, config `eps_f_relax`), default off at this φ but present.

### 3.3 Field storage

ε_f, ∂ε_f/∂t (finite-difference in time, stored even though unused by the constraint in this milestone), and f_{p→f} as registered solver fields with checkpoint/restart and VTK output support.

---

## 4. Phase 2 — One-way coupling (fluid → particles)

1. New force-model module (`el_forces`) implementing §1.2 closures as pure functions of interpolated fluid state + particle state. Each force term individually unit-tested against literature values.
2. Wire into the pe step: in E-L mode, the hydrodynamic force handed to pe per particle is `F^drag + F^∇p + m g + F^lift` instead of the FBM stress integral. Torque from hydrodynamics: zero in this milestone (contacts/lubrication still supply torque).
3. **Semi-implicit drag in the particle update.** Write drag as `F^drag = B_i (u_f(X_i) − U_i)` with `B_i ≥ 0` evaluated at the current state; update U_i with U_i treated implicitly (linear in U_i ⇒ closed-form per particle, or exponential integrator). This must be robust for `τ_p = m_i/B_i ≪ Δt_DEM`.
4. DEM sub-cycling: `n_sub` pe steps per fluid step with hydrodynamic forces frozen (re-evaluated each substep from frozen fluid fields but current particle state). Config: `n_sub` or automatic from contact stiffness criterion.

**Acceptance (Validation V1):** single sphere settling in quiescent fluid, ρ_p/ρ_f ∈ {1.2, 2.5}, Re_p ∈ {0.05, 1, 10}: terminal velocity within 2% of the implied drag-law value; velocity relaxation matches analytic exponential for the Stokes case; result independent of Δt over a 10× range (semi-implicit drag working); result independent of mesh h at fixed δ.

---

## 5. Phase 3 — Two-way coupling (particles → fluid)

1. Assemble `b^p` from §3.2 into the momentum RHS of every fractional-step-θ substep, with consistent θ-weighting (treat f_{p→f} like the body-force term is currently treated; document the choice).
2. **Semi-implicit drag on the fluid side.** Split the drag feedback: the part proportional to u_f goes into the system matrix as a nonnegative diagonal(ish) reaction term `+Σ_i B_i W-projected`, the part proportional to U_i stays on the RHS. This adds a Helmholtz-like term the multigrid solver should handle well; verify smoother robustness. Config switch `drag_coupling = {explicit, semi_implicit}`, default semi-implicit.
3. Global momentum audit (diagnostic, every N steps): fluid momentum change vs Σ forces exchanged; report drift.

**Acceptance (Validation V2):** Richardson–Zaki hindered settling. Periodic box, homogeneous suspension, φ ∈ {0.05, 0.10, 0.15, 0.20}: mean settling velocity follows `U(φ)/U_0 = ε_f^n` with n ≈ 4.65 (low Re_p) within ~10%. Additionally: zero net momentum drift in a force-free periodic suspension; no negative ε_f anywhere; results stable for at least 10⁵ steps.

---

## 6. Phase 4 — Four-way coupling

1. Enable pe contacts + lubrication in E-L mode (they should already work; verify lubrication gap h is computed from actual particle radii, unrelated to mesh).
2. Verify the combined sub-cycled loop: interpolate → forces → pe substeps (contacts/lubrication inside) → scatter → fluid solve.
3. Stress test: φ = 0.20 sheared periodic box; no particle overlap beyond contact-solver tolerance; energy bookkeeping plausible (no spontaneous heating).

---

## 7. Phase 5 — Pipe flow validation campaign

Geometry: straight circular pipe, axially periodic, driven by constant body force equivalent to a prescribed −dp/dz; wall = no-slip + pe wall for contacts/lubrication. Pipe radius R ≥ 10 d_p (kernel-width constraint; assert at startup).

- **V3 — Segré–Silberberg:** neutrally buoyant particles, dilute (φ ≈ 1–2%), channel Reynolds ~ O(10–100), lift ON: radial equilibrium position ≈ 0.6 R (±0.05 R). With lift OFF: no migration (negative control). This isolates the lift closure.
- **V4 — Pressure drop vs φ:** φ ∈ {0.05, 0.10, 0.15, 0.20}, fixed flow rate; compare apparent viscosity `μ_app = −(dp/dz) R²/(8 ⟨u⟩)`-style evaluation against (a) in-house FEM-FBM DNS of the same configuration at matched parameters, and (b) Hagen–Poiseuille with μ_eff(φ) from the published RVE viscosity tables. Expected: agreement trend-wise; E-L may undershoot DNS by single-digit percent at φ = 0.15–0.20 — record the gap per φ, it is the dense-regime readiness metric.
- **V5 — Radial concentration profiles** at φ = 0.10–0.20 vs DNS: qualitative agreement (wall depletion layer ~ d_p/2 by geometry; migration trends).

Deliverable: `docs/el_validation_report.md` with all plots, parameter files committed under `cases/el_validation/`.

---

## 8. Configuration surface (new keys, with defaults)

```
simulation_mode        = DNS_FBM | EULER_LAGRANGE
el.kernel              = gaussian | deen_poly        (default deen_poly)
el.kernel_width_factor = 2.5                          (δ = factor · d_p)
el.eps_f_min           = 0.4
el.eps_f_relax         = 0.0                          (temporal smoothing, off)
el.drag_model          = difelice | wenyu | stokes
el.lift_model          = none | saffman_mei | saffman_mei_wall
el.magnus              = false
el.pressure_force      = true                         (−V_i ∇p)
el.drag_coupling       = semi_implicit | explicit
el.n_sub               = auto | <int>                 (DEM substeps per fluid step)
el.momentum_audit_freq = 100
```

---

## 9. Known pitfalls (encode as assertions/tests where possible)

1. **Convention mixing:** drag correlation assumes Model-A decomposition (fluid carries full −∇p; particle gets −V_i∇p; only drag+lift fed back). Assert via the V1/V2 tests; document in module header.
2. **Q2 negativity:** never use raw Q2 shape values as transfer weights (assert kernel weights ≥ 0).
3. **Kernel/wall truncation:** unrenormalized kernels near the pipe wall silently destroy mass/momentum conservation — covered by the Σw=1 unit test with a wall-adjacent particle.
4. **Parallel spreading:** a particle near a subdomain boundary must deposit force/volume into neighbor-rank elements (ghost contributions) — covered by the multi-rank conservation tests.
5. **Stiff drag:** small particles in liquid ⇒ τ_p possibly ≪ Δt; both particle and fluid updates must be semi-implicit (V1 Δt-independence test).
6. **Double counting viscosity:** μ_f stays μ_f; μ_eff tables only for comparison in V4, never inside the E-L solver.
7. **∂ε_f/∂t noise:** field is assembled and stored but unused by the constraint here; if later enabled, requires temporal smoothing (hook exists).
8. **Restart correctness:** ε_f history needed for ∂ε_f/∂t must be checkpointed.

---

## 10. Suggested execution order for Codex/Claude Code

1. Phase 0 map → review checkpoint with the user.
2. Phase 1 transfer module + unit tests (no physics yet) → CI green on 1/2/8 ranks.
3. Phase 2 forces + V1 → review checkpoint (plots).
4. Phase 3 feedback + V2 → review checkpoint.
5. Phase 4 four-way + stress test.
6. Phase 5 pipe campaign + validation report.

Keep each phase a separate PR/branch; never modify FBM code paths except to factor out shared interfaces; every new physics formula carries a literature citation in a comment and a unit test.
