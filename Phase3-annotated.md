# Euler-Lagrange Phase 3: Explicit Two-Way Fluid Feedback — Annotated Plan

> **About this document.** Phase 3 plan with review annotations folded in.
> Annotations are marked `> **Review note:**` and are grounded in the current
> code (`source/src_quadLS/QuadSc_main.f90`, `QuadSc_corrections.f90`,
> `source/src_el/*`). Where an annotation *changes* the original intent, the
> original bullet is struck through and the corrected version follows.

> **Headline (the load-bearing question).** *Are the right force contributions
> spread to the fluid?* **Yes.** Spreading `feedback_force = drag + lift` and
> excluding pressure, gravity, and buoyancy is the correct dilute point-particle
> two-way reaction, and the buoyancy exclusion is exactly consistent with the
> Phase 2 gravity-free-fluid decision. See the physics note under "Key Changes".
> The remaining risks are mechanical (shared `fluid_core`, explicit-vs-implicit
> drag, stability) and are flagged inline.

## Summary

Implement opt-in dilute two-way E-L coupling for `q2p1_el_pipeflow`: use the
pre-advance particle feedback force as the Newton-pair source for the fluid
momentum RHS, while preserving Phase 2 one-way behavior by default. Keep coupling
explicit and lagged; do not add porosity-modified equations or semi-implicit
fluid-side drag.

> **Review note.** "Newton-pair" overstates it — see issue (D) below. The fluid
> receives the *explicit, old-slip* feedback, while PE advanced the particle with
> the *implicit* drag, so the per-step exchange is conserved only to O(Δt). Treat
> it as a lagged first-order approximation, not an exact pair.

## Key Changes

- Add `ELApplyFluidFeedback = Yes|No`, default `No`, parsed and printed with the
  existing E-L config; set sample pipeflow data to `No` unless a two-way testcase
  explicitly enables it.
  > **Review note.** Mirror the existing boolean parse pattern
  > (`ELApplyParticleForces`/`ELPressureForce` in `app_init.f90`: read a token,
  > compare to `"Yes"/"YES"`), add to `EL_VALIDATE_CONFIG`/`EL_PRINT_CONFIG`.

- Preserve the pre-advance deposited `feedback_force = drag + lift` in a new
  persistent E-L source buffer, separate from diagnostic `force_rhs`, because the
  post-advance diagnostic pass overwrites `force_rhs`.
  > **Review note (sign — make this explicit).** The buffer must store the
  > **deposited** quantity, i.e. the sign-flipped spread for which
  > `Σ_dof = −Σ_i feedback_i` (the EL deposit already produces `Σ_a b_a = −ΣF_i`).
  > That is the *reaction* the fluid feels, and it is what gets added with a `+`
  > sign alongside gravity (below). Storing `+feedback` would invert the coupling.
  > Capture the buffer in the **pre-advance pass only** (`advance_history`), which
  > is the pass that also arms PE — that snapshot is the force consistent with the
  > step that is about to run.

  > **Physics note — WHY drag+lift and nothing else (the answer to the review
  > question).** For dilute point-particle two-way coupling the fluid feels the
  > negative of the particle *disturbance* forces only:
  > - **drag + lift** → spread. Genuine fluid→particle interaction; reaction
  >   `−(drag+lift)` goes to the fluid.
  > - **dynamic pressure `−V∇p`** → excluded. The undisturbed-flow pressure
  >   gradient is already carried by the fluid's own `∇p`; reacting it back
  >   double-counts fluid momentum.
  > - **buoyancy** → excluded. Buoyancy is `−V∇p_hydrostatic`, but the fluid
  >   momentum equation carries **no gravity** (Phase 2 convention), so no
  >   hydrostatic `∇p` exists in the fluid — buoyancy is a driver-side closure,
  >   not a fluid-mediated force. Spreading it back would inject spurious momentum
  >   into a fluid with no compensating hydrostatic field.
  > - **gravity** → excluded. External body force, never transmitted to the fluid.
  > This is exactly `feedback_force` as defined in Phase 2, so no new force math
  > is needed — only plumbing.

- Add a small source-insertion helper that adds `tstep * feedback_source(:,dof)`
  to `QuadSc%defU/defV/defW`; do not divide by density and do not include
  pressure, gravity, or buoyancy.
  > **Review note — density/`tstep` convention VERIFIED correct.** `AddGravForce`
  > calls `Grav_QuadSc(defU,…, mgDensity, Gravity, …, tstep, E013)`, i.e. it
  > assembles `tstep · ∫N_a (ρg) dV` — a force-density load scaled by `tstep`,
  > **not** divided by density. `feedback_source` is already a *force* (drag/lift
  > in Newtons, spread to a consistent load), so `tstep * feedback_source` with no
  > density factor matches the gravity convention exactly. Note also: do **not**
  > multiply by density either (gravity needs `ρ` only to turn `g` into a force;
  > the feedback is already a force).

- In `Transport_q2p1_UxyzP_fluid_core`, apply the E-L source after
  `AddGravForce()`/`AddConstantForce()` and before `Boundary_QuadScalar_Def()`
  and before saving `QuadSc%rhs*`, so nonlinear iterations restore the same
  particle source.
  > **Review note — insertion point VERIFIED.** `QuadSc_main.f90:441-450`:
  > `AddGravForce` (441) → `AddConstantForce` (442) → [insert here] →
  > `Boundary_QuadScalar_Def` (444) → `QuadSc%rhsU = QuadSc%defU` (448-450). Exact
  > match; the source lands in `rhs*` and is restored across inner iterations.

  > **Review note — ISSUE (C): `fluid_core` is SHARED; guard the insertion.**
  > `Transport_q2p1_UxyzP_fluid_core` is called by the EL path
  > (`...UxyzP_el`, line 369, `enable_fbm=.FALSE.`) **and** by the FBM path
  > (`...UxyzP_fc_ext`, line 380, `enable_fbm=.TRUE.`); the header comment says
  > "Shared low-level Q2/P1 solve." The EL source buffer lives in modules compiled
  > only under `BUILD_Q2P1_EL_PIPEFLOW`, so the insertion must be:
  >   1. `#ifdef BUILD_Q2P1_EL_PIPEFLOW` (so non-EL builds even compile), and
  >   2. runtime-gated on `el_apply_fluid_feedback` (default `No`).
  > The FBM/non-EL path must be a guaranteed no-op (add a regression assertion —
  > see tests).

- Keep `ELApplyParticleForces` independent: it controls PE particle forcing only;
  `ELApplyFluidFeedback` controls fluid reaction only.
  > **Review note.** Good. Document the diagnostic combo
  > `ApplyParticleForces=No` + `ApplyFluidFeedback=Yes` as *intentionally*
  > breaking the Newton pair (fluid reacts to forces the particles never
  > received) so nobody mistakes it for a physical mode.

- Carry one Phase 2 cleanup: set the `EL_CONFIG` default drag model to `difelice`
  so module defaults match the pipeflow Phase 2 default.
  > **Review note.** Fine and harmless; the `.dat` already overrides to
  > `difelice`, this just aligns the module default.

## Test Plan

- Extend E-L transfer tests to assert the pre-advance source buffer equals the
  spread feedback (sign-flipped reaction), excludes pressure/gravity/buoyancy,
  and remains available after the post-advance diagnostic refresh.
- Add source-insertion tests for exact `tstep` scaling, no-op behavior when
  `ELApplyFluidFeedback=No`, and preservation through `QuadSc%rhs*` restore
  semantics.
- Add conservation tests: summed fluid source equals `-sum(feedback_force)` in
  serial, MPI-2, and MPI-8 registrations.
  > **Review note — scope of this test.** This checks **P2G self-consistency**
  > (the spreading conserves the deposited force), **not** fluid↔particle momentum
  > balance. Because the spread feedback is explicit/old-slip while PE applied the
  > implicit drag (issue D), do **not** advertise this as exact Newton-pair
  > conservation.
- Add runtime checks:
    - feedback off reproduces Phase 2 behavior within existing tolerances;
    - feedback on gives nonzero momentum RHS and stable short dilute run;
      > **Review note — issue (E): bracket the stability limit.** Explicit lagged
      > two-way coupling has the same time-Nyquist instability that motivated the
      > `α=0.35` under-relaxation on the *particle* side. Push mass loading up in
      > this test until it breaks, so the dilute bound is measured, not assumed.
    - `ELApplyParticleForces=No` with `ELApplyFluidFeedback=Yes` computes
      diagnostics but does not advance PE forces unless explicitly intended.
- **(added)** Assert the FBM/non-EL path through `fluid_core` is unchanged with
  feedback off (guards issue C — the shared-routine hazard).
- Re-run `test_el_kernel_forces`, `test_el_transfer`, the PE semi-implicit drag
  test, and a lightweight `q2p1_el_pipeflow` smoke run.

## Issues raised by review (summary)

- **(C) [Med] Shared `fluid_core`.** Guard the insertion with
  `#ifdef BUILD_Q2P1_EL_PIPEFLOW` + the runtime flag; verify the FBM path is a
  no-op. *(highest-impact correctness item)*
- **(D) [Med] Explicit vs implicit drag.** The spread reaction (old slip) is not
  the force PE applied (implicit). O(Δt), acceptable for dilute/lagged, but reword
  the "Newton-pair" framing and don't claim exact conservation. Optionally spread
  the implicit reaction for strict balance.
- **(E) [Med] Explicit-coupling stability.** No relaxation on the fluid source;
  dilute-only. Document the loading bound or add optional under-relaxation
  mirroring the particle side.
- **(B/sign) [Low] Buffer holds the deposited reaction** (`Σ = −Σ feedback`),
  captured in the pre-advance pass, added with `+tstep`.

## Assumptions (annotated)

- Phase 3 means explicit two-way dilute coupling, not only diagnostic P2G
  hardening.
- ~~The fluid source uses the same pre-advance feedback force that drove the PE
  update for per-step Newton-pair accounting.~~ → **Corrected:** the fluid source
  uses the pre-advance feedback (explicit, old slip). PE advanced the particle
  with the *implicit* drag, so this is a lagged first-order approximation of the
  Newton pair, not an exact per-step balance (issue D).
- `force_rhs` remains exported as diagnostics; the new pre-advance source buffer
  is the only field used for momentum insertion.
- No fluid-side semi-implicit drag matrix, porosity momentum correction, hindered
  settling, lubrication closure, or pressure/gravity/buoyancy feedback is added in
  this phase.
  > **Review note.** Consistent: excluding porosity is *why* spreading only
  > drag+lift (rather than a volumetric momentum-exchange term) is the correct
  > closure here.

## Suggested sequencing

1. **Config + buffer:** add `ELApplyFluidFeedback`, the persistent pre-advance
   source buffer, and the transfer-test assertions (buffer == sign-flipped spread,
   excludes pressure/gravity/buoyancy, survives the diagnostic refresh).
2. **Guarded insertion** into `fluid_core` (`#ifdef` + flag), with the FBM-no-op
   regression test (issue C) and the `tstep`/no-density scaling test.
3. **Conservation tests** (P2G self-consistency, serial/2/8 ranks).
4. **Runtime:** feedback-off ≡ Phase 2; feedback-on dilute run, pushed to find the
   stability bound (issue E).

