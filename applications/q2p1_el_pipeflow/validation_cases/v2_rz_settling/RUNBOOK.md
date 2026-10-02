# Runbook: v2_rz_settling — Richardson–Zaki hindered settling (V2, Stage 1)

Goal: the drag arbiter for the dense regime. Batch settling in a fully
periodic unit box at φ = 0.05…0.20; the lab-frame mean particle velocity
must follow U/U₀ = ε^n with n ≈ 4.65 (Richardson–Zaki, Re_p ≪ 1). Di Felice
predicts n = χ+1 → 4.7 at low Re (derivation in the campaign plan): RZ is
the empirical arbiter for the voidage exponent, complementing the ten Cate
case (single-particle end, V1b).

## Configuration

- Mesh: QBOX9 unit box [0,1]³, reused from tier2 momentum_conservation
  (mesh symlink). Production MaxMeshLevel = 3: 36³ cells, h = 1/36 →
  δ = 2.5·d_p = 0.125, δ/h = 4.5 (campaign standard).
- CFD periodicity: `SimPar@Periodic = Yes` (all axes, match length 1.0 —
  the tier2 momentum path). PE periodicity: `periodicX_/Y_/Z_: true`
  (decomposePeriodic3D, first production consumer of the per-axis
  dispatch), PE grid 3×3×3, np = 28.
- Particles: d_p = 0.05 (benchRadius_ 0.025), ρ_p = 2, ρ_f = 1, μ = 0.02
  (ρ=1 ⇒ Prop@Viscosity kinematic==dynamic). Stokes settling
  v_t = Δρ·g·d²/(18μ) = 6.8e-3 at g = 0.981 → Re_p ≈ 0.017 (RZ n = 4.65
  regime). τ_p = 1.39e-2, dt = 0.002 (dt/τ_p = 0.14), PE substeps_ = 2 for
  contact resolution (macro stepsize_ = dt = 0.002; dry contacts,
  restitution 0 defaults).
- Zero-net-flux counter-force on the fluid: `ConstantForcing` z =
  +φ·Δρ·g = 0.2·1·0.981 = 0.1962 (staged φ = 0.20 value — RESCALE per φ!).
- Seeding: `seedMode_ = random`, `volumeFraction_` per run, seed_ 12345,
  default seedMinGap (0.1·d_p), box inset per the W0.3 contract.
- N per φ (unit box, V_p = 6.545e-5): 0.05 → 764, 0.10 → 1528,
  0.15 → 2292, 0.20 → 3056. U₀ reference run: seedMode_ file with a single
  particle (particles.xyz placeholder in this case dir) and
  UseConstantForcing = No.

## Staged default = MANDATORY φ = 0.20 smoke (user-required gate)

200 steps, np = 28, before ANY long Slurm run. Stage into an OWN run dir:

```bash
tools/el_stage_rundir.sh v2_rz_settling <rundir>
# env: module load gcc/latest-v13 openmpi/…/4.1.6; OMPI_MCA_rmaps_base_oversubscribe=1
mpirun --oversubscribe --wdir <rundir> -np 28 "$BIN" > <rundir>/run_v2_smoke.log 2>&1
```

Smoke gates:
- Seeder: N = 3056 created, achieved φ within 2% (asserted), min-gap
  assertion silent, `EL_SEED_POS` list plausible.
- Startup banner `fully periodic PE decomposition` (dispatch path).
- No fatal/OOB/clipping; EL_VOLUME_CONSERVATION at machine zero;
  ⟨ε_f⟩ ≈ 0.80 (1 − φ) within 1% on EL_MEAN_SLIP.
- EL_CONTACT_STATS: max_overlap ≤ 1% d_p, no runaway Tgran.
- Per-rank balance sane (EL_HALO_RECORDS), memory headroom at 36³,
  measured wall-time/step recorded here → sizes the Slurm allocations.

## Production runs (Slurm, after the smoke passes)

Five runs: φ ∈ {0.05, 0.10, 0.15, 0.20} + U₀ (N=1). Per run adjust:
`volumeFraction_` (json), `ConstantForcing` z = φ·Δρ·g (param),
MaxNumStep ~10–20k (settling statistics; TAvg window 0.5).
- RZ metric (CORRECTED from the original plan): the periodic box with
  counter-force conserves mixture MOMENTUM (ρ_f·uf_super + φ·ρ_p·⟨u_p⟩ ≈ 0
  — that is the 1% cross-check), NOT mixture volume flux; with ρ_p ≠ ρ_f
  the two differ. RZ batch settling is defined in the zero-VOLUME-FLUX
  frame, so shift: U_RZ = ⟨u_p,z⟩ − (uf_super,z + φ·⟨u_p,z⟩), with
  ⟨u_p⟩/uf columns from `EL_MEAN_SLIP_TAVG`. Second cross-check:
  slip_intr ≈ ⟨u_p⟩ − uf_intr (backflow bookkeeping).
- Acceptance: per-φ |⟨u_p⟩/U₀ − ε^4.65|/ε^4.65 ≤ 10%; offline fit
  n = 4.65 ± 0.5; overlap ≤ 1% d_p; ⟨ε⟩ = 1−φ within 1%.
- Debug hint: n ≈ 3.7 instead of ≈ 4.7 ⇒ slip/superficial bookkeeping
  error in the drag ε usage, not the closure.
- Known bias (from V1b/W4.1): kernel self-voidage ε_eff < 1 is PART of the
  intended dense-suspension physics here — do not "correct" it; at φ ≥ 0.05
  the real neighbors dominate the deposit anyway.

## Results — φ = 0.20 smoke (2026-07-05, np=28, level 3, 200 steps)

PASS on all gates:
- Seeder: 3056 spheres, achieved φ = 0.2000147, min-gap silent;
  `fully periodic PE decomposition` banner (first runtime exercise of
  decomposePeriodic3D through the per-axis dispatch).
- ⟨ε⟩ = 0.7999853 (1−φ to 0.002%); no fatal/OOB/clipping;
  EL_CONTACT_STATS: ncontacts ≈ 5.6k (rank-summed), max_overlap = 0.0,
  Tgran 4.8e-6.
- Physics preview at t = 0.4: ⟨u_p,z⟩ = −2.271e-3; mixture momentum zero
  to 0.4%; frame-shifted U_RZ/U₀ ≈ 0.385 vs ε^4.65 = 0.354 (+8.8%,
  already inside the ±10% production gate).
- Wall-time ≈ 15 s/step at np=28 — but measured while SHARING the 6-core
  node with a 28-rank frozen run; budget ~8–10 s/step exclusive →
  production 15k steps ≈ 35–40 h ⇒ Slurm mandatory.

The FIRST smoke attempt (with `Prop@Gravity` acting on the fluid) sent the
whole mixture into free fall at −0.2: a fully periodic box has no
hydrostatic pressure to absorb ρ_f·g. This is why `ELFluidGravity = No`
exists (fluid gravity absorbed into the modified pressure; particles keep
grav_buoy analytically) — REQUIRED for all periodic-suspension cases,
harmless default Yes everywhere else.

## Results — production sweep, FIRST ATTEMPT (2026-07-05/06): ANOMALY

- U₀ reference (job 131966): slip = 8.24e-3 vs algebraic Di Felice
  7.10e-3 → +16% two-way co-flow bias (V1b effect, stronger without
  walls). Usable as the two-way-consistent U₀.
- φ = 0.05 (job 131962): settles correctly at the hindered scale
  (⟨u_p,z⟩ ≈ −9e-3) until t ≈ 25, then the WHOLE system (both phases)
  drifts upward, linear in t (+2e-3/t.u.) — a sustained internal
  momentum leak ≈ 5% of the counterforce (external forces balance to
  5e-6; ε-weighted total momentum audit grows). Onset coincides with
  clustering (ncontacts 250 → 520, Tgran 4e-5 → 1e-3). U_RZ from this
  run is NOT quotable.
- Leak investigation (task ongoing): ruled out — PE substep-contact
  drag interaction (control pair D0/D1, jobs 132110/132111: substeps
  2 vs 1 leak identically, also reproduces at level 2); post-advance
  deposit overwrite (capture is advance_history-gated); deps_f_dt
  (diagnostic-only). Under test via the EL_NEWTON_PAIR audit:
  PE-impulse-vs-mirror mismatch vs fluid-side discrete momentum error
  (non-conservative convection under developing pseudo-turbulence is
  the remaining fluid-side candidate).
- φ = 0.10/0.15/0.20 (131963-965): same substeps-2 physics — will show
  the same leak; useful as φ-scaling evidence, not as RZ points.

The RZ fit is BLOCKED on resolving the leak; the sweep will be
re-submitted afterwards.

## Leak resolution, part 1 — PE Newton-pair fix (pe commit 8b037ae)

The EL_NEWTON_PAIR audit split the leak in two. The PARTICLE-side part
(~4e-6/step under contact, machine zero without) was root-caused to four
mechanisms in HardContactEulerLagrange: (1) elRelax 1.0-vs-0.9 fold
asymmetry on cross-rank contacts, (2) pe::equal epsilon zero-test dropping
small remote corrections, (3) migration leaving stale armed hydro state
(1.9× forcing), (4) substep drag computed from contact-perturbed velocity.
Fixed in pe 8b037ae (fold-once-at-caching, exact zero test, migration
payload handover, elRefVelocity_ free-flight reference); permanent unit
reproducers: unit_straddling_contact, unit_migration_hydro,
unit_periodic_contact.

## Results — production sweep, SECOND ATTEMPT with fixed PE (2026-07-08/09, jobs 132467–132471)

All five runs (φ = 0.05/0.10/0.15/0.20 to t = 100, U₀ to t = 10)
completed on the fixed binary.

- EL_NEWTON_PAIR at machine zero at scale for the entire sweep: worst
  |mismatch| = 1.6e-18 (φ=0.05), 5.0e-18 (0.10), 1.5e-17 (0.15),
  1.9e-17 (0.20), 5.6e-22 (U₀). The particle-side leak is CLOSED.
- U₀ reference (job 132471): slip_intr,z = −8.249e-3 — reproduces the
  first-attempt value (8.24e-3, +16% two-way co-flow bias over algebraic
  Di Felice 7.10e-3). Usable as the two-way-consistent U₀.
- BUT the fluid-side leak persists and dominates by t = 100: uf_super,z
  drifts to +4.84e-2 (φ=0.05), +6.31e-2 (0.10), +1.27e-2 (0.15),
  +1.55e-2 (0.20) — non-monotone in φ, large enough that the
  frame-corrected U_RZ is unphysical (e.g. |U_RZ| ≈ 5×|U₀| at φ=0.05).
  The RZ fit remains NOT quotable from this sweep.

## Leak resolution, part 2 — fluid-side localization (EL_FLUID_PAIR, job 132698)

Dedicated φ = 0.05 audit run (20k steps to t = 40, audit every 25 steps,
src_el commits bfadfe9b + 2c6ee44d): EL_FLUID_PAIR compares the fluid
momentum change across `Transport_q2p1_UxyzP_fluid_core` against
dt·(Σ feedback source + external force·V) each audited step.

- The mismatch is REAL, z-dominant, and GROWS with the developing
  pseudo-turbulence: ~1e-11/step at t≈0 → ~3e-7/step at t = 40;
  cumulative over the 800 audited steps +7.4e-4 in z (order-consistent
  with the observed uf_super drift given the 1-in-25 audit sampling).
- EL_NEWTON_PAIR simultaneously at machine zero (worst 8e-19) in the same
  run — the two audits together prove the residual leak is INTERNAL to
  the fluid momentum solve (transfer layer and PE both conserve exactly).
- Prime suspect: non-conservative discrete convection of the growing
  velocity fluctuations (the mismatch tracks fluctuation amplitude, not
  contact activity). Next step: bisect inside the fluid core
  (convection / diffusion / pressure sub-stages).

The RZ fit stays BLOCKED on the fluid-side fix; the sweep re-runs once
EL_FLUID_PAIR is at machine zero.

## Leak resolution, part 3 — root cause found and FIXED (2026-07-23, commits 2022e126 + 879933b1)

Stage-split audit (EL_FLUID_STAGE, jobs 135770/135771/135772, φ=0.05
level 2):

1. The mismatch splits as mis_mom (momentum/Burgers stage) + mis_prj
   (projection). Measured: mis_prj at machine zero (worst 1.5e-17 over
   5000 steps); mis_mom carries the ENTIRE mismatch to every digit.
2. Tightening the momentum solver (defCrit 1e-6 → 1e-10, NLmax 8)
   leaves the mismatch unchanged to ~0.1% — not solver truncation.
3. Mass, diffusion and projection are pointwise momentum-neutral
   (partition of unity + row-sum lumping), sources audited ⇒ the leak
   is the Galerkin CONVECTIVE FORM in CONVQ2: it injects
   −ρ∫u_i(div u)dV per step, nonzero under the weakly-enforced Q2/P1
   divergence, growing with the pseudo-turbulent fluctuations — exactly
   the observed behaviour (onset with clustering, z-dominant).

Fix: `SimPar@ELConvectionForm = divergence` (commit 879933b1) switches
CONVQ2 to the integrated-by-parts form K_JI = −ρ∫(u φ_I)·∇φ_J, whose
column sums vanish POINTWISE — momentum-exact in floating point on
periodic domains regardless of quadrature and div u. Verification (job
135772, 5000 steps vs legacy job 135770): mismatch x/y ~1e-20; z flat
at 1.94e-13 the whole run (constant audit-bookkeeping offset between
the dvol volume in the expected term and the lumped-mass body-force
application — NOT dynamics; cumulative 3.9e-11 vs legacy 7.4e-4, seven
orders). Physics unchanged: slip −7.839e-3 vs −7.851e-3 (<1%, the
expected discretization-form difference), identical wall time and
solver convergence. Default stays `convective` — tier2 suite
bit-compatible, no other application affected.

REQUIRED for all periodic-suspension production runs (like
ELFluidGravity = No): `SimPar@ELConvectionForm = divergence`. The RZ
sweep re-runs with this key; the fit is UNBLOCKED.

## Leak resolution, part 4 — divergence form not energy-robust; measured-leak compensator (2026-07-24, commit 22da0d5c)

The divergence-form RZ sweep (jobs 135773–135776) FAILED: NaN in the
momentum solve at t = 41.5 (φ=0.05), ~31 (0.10), ~32 (0.15), 48.6
(0.20). Momentum conservation held to the last audit (~1e-12), but
Tgran grew ~2.5× faster than the legacy runs at equal times before
blowup — the divergence form is momentum-exact but NOT energy-robust
under developed pseudo-turbulence (the densest, most drag-damped case
survived longest; the legacy convective form ran identical configs to
t=100). Keep `divergence` for short audit/verification runs only.

Production answer: `SimPar@ELMomentumFix = Yes` with the LEGACY
convective form. The EL_FLUID_PAIR mismatch is computed every step and
its negative applied one step lagged as a uniform body force
(Grav_QuadSc pathway, deadbeat recursion C_{n+1} = C_n − mismatch_n):
cumulative fluid-momentum drift bounded by one step's leak (~1e-7)
instead of the linear ~5.7e-4/t.u. ramp. Verified (job 135865, level 2,
5000 steps): residual mismatch sign-alternating ~1e-10 (worst 6.0e-10,
cumulative 4.1e-9, no growth), physics within 0.3% of both forms,
stability unchanged.

REVISED prescription for periodic-suspension production runs:
`ELFluidGravity = No` + `ELMomentumFix = Yes` (legacy convection).
EMAC-form convection (momentum+energy conserving) is the future proper
discretization fix. The RZ sweep re-runs with the compensator.

## Sweep 3 (compensated) — infrastructure clean, RZ gate NOT met: mesoscale settling instability (2026-07-25, jobs 135866–135870)

Configuration: level-3 production sweep, `ELFluidGravity = No` +
`ELMomentumFix = Yes` (legacy convection), 50000 steps to t=100,
np=28. U₀ reference v2mfx_u0 (135870, 5000 steps to t=10). All five
runs COMPLETED (10.6h / 16.3h / 20.6h / 24.2h / short).

### Infrastructure verdict — everything the compensator sweep was gated on: PASS

| check | φ=0.05 | φ=0.10 | φ=0.15 | φ=0.20 |
|---|---|---|---|---|
| NaN/Inf | none | none | none | none |
| worst per-step mismatch_z | 2.1e-8 | 3.0e-8 | 4.7e-8 | 9.9e-8 |
| worst cumulative drift_z | 9.6e-8 | 5.2e-7 | 2.9e-7 | 3.7e-7 |
| worst Newton-pair_z | 1.3e-18 | 7.3e-18 | 1.5e-17 | 1.7e-17 |
| max_overlap (every sample) | 0.0 | 0.0 | 0.0 | 0.0 |
| ⟨ε_f⟩ vs 1−φ | 5 digits | 5 digits | 5 digits | 5 digits |

The cumulative fluid-momentum drift is bounded and sign-alternating
(vs the previous linear 5.7e-4/t.u. ramp → ~5.7e-2 by t=100: five
orders better). The energy instability of the divergence form does not
appear: legacy convection + compensator is stable to t=100 at all φ.
The momentum ledger of these runs is trustworthy; what follows is a
physics result, not an accounting artifact.

### RZ verdict — FAIL: settling ENHANCEMENT, not hindrance, at late time

Frame-corrected metric U_RZ = ⟨u_p,z⟩ − j, j = uf_super,z + φ⟨u_p,z⟩
(mixture volume flux; the balanced-body-force box settles into the
zero-total-momentum state, uf_super ≈ −(ρ_p/ρ_f)φ⟨u_p,z⟩ with
ρ_p/ρ_f = 2, so the zero-net-volume-flux RZ frame differs from the lab
frame — the correction accounts for it). U₀ = −8.2391e-3
(frame-corrected u0 TAVG; consistent with all previous references to
4–5 digits).

TAVG window t=[50,100]:

| φ | U_RZ | U/U₀ | RZ ε^4.65 | verdict |
|---|---|---|---|---|
| 0.05 | −2.179e-2 | 2.645 | 0.788 | FAIL (+236%) |
| 0.10 | −2.759e-2 | 3.349 | 0.613 | FAIL (+447%) |
| 0.15 | −2.728e-2 | 3.311 | 0.470 | FAIL (+605%) |
| 0.20 | −2.486e-2 | 3.017 | 0.354 | FAIL (+752%) |

Mean settling is 2.6–3.3× FASTER than the single particle, and
non-monotonic in φ. No RZ exponent is quotable from this state.

### Two-regime structure — the homogeneous state is unstable

Early window t=[5,15], before the instability develops:

| φ | U/U₀ | RZ ε^4.65 | n_eff |
|---|---|---|---|
| 0.05 | 0.922 | 0.788 | 1.58 |
| 0.10 | 0.812 | 0.613 | 1.98 |
| 0.15 | 0.741 | 0.470 | 1.85 |
| 0.20 | 0.511 | 0.354 | 3.01 |

Early time IS hindered (U/U₀ < 1, ordered in φ) but weaker than RZ
(fit n ≈ 2.5–3.3 vs 4.65) — and the window is itself still transient:
the box-scale backflow structure needs ~L/U ~ 100 t.u. to develop,
while the instability overtakes it from t ≈ 20.

Instability evidence (all φ, strongest at 0.20):
- Tgran grows two orders of magnitude exactly at the enhancement
  onset (φ=0.20: 5e-5 at t=15 → 9.7e-3 at t=40, decaying to 2.9e-3 by
  t=100 — not yet statistically stationary at end of run).
- Proximity-pair count rises ~50% from its seeded value (φ=0.05:
  ~240 → ~550 peak) — clustering.
- ⟨u_p,z⟩ shows slow box-scale oscillations (period ~40 t.u.,
  amplitude 2–3×) — plume/convection-cell dynamics.
- Ga = √(g'd³)/ν ≈ 0.55: far below any particle-inertia clustering
  threshold — this is the two-way-coupled concentration (plume)
  instability of a periodic settling suspension, not turbulence.

Interpretation: the RZ correlation describes wall-bounded batch
sedimentation, where container walls suppress box-scale convection.
A triply periodic box with uniform counter-forcing has no such
mechanism: concentration fluctuations self-amplify (dense regions drag
fluid down, reducing local slip and collecting more particles) and the
saturated state is cluster-induced enhanced settling. Candidate
contributors to the growth rate that are OUR modelling choices rather
than physics: no undisturbed-field correction (each particle feels its
own kernel-smoothed wake — deferred item), kernel δ ≈ mean
interparticle spacing at φ=0.05, and the box being only L/d_p = 20.

### Options forward (decision needed before Stage 5 gating)

1. Suppress the mesoscale mode: horizontal-plane-wise zero-net-flux
   forcing (subtract the xy-mean fluid velocity per plane each step),
   the standard trick for periodic hindered-settling measurements.
2. Shrink the box (L/d_p ~ 8–10): the instability's fastest-growing
   wavelength is box-limited; a smaller box may hold the homogeneous
   state long enough to measure n. Trades against ⟨ε_f⟩ statistics.
3. Early-window fit only: quote n from t=[5,15] with the caveat that
   the backflow is not fully developed (current n ≈ 2.5–3.3).
4. Reframe V2: validate against cluster-induced-settling literature
   instead of RZ (enhancement factors 2–4× are reported for two-way
   coupled point-particle sedimentation) and gate Stage 5 on the
   infrastructure checks (all PASS) plus V4's own acceptance.

Stage-5 note: the pipe cases are wall-bounded — the mesoscale mode
that breaks the periodic box is suppressed there by construction, so
Stage 5 is NOT invalidated by this result.

## Mesoscale lane filter (option 1) — probe verdict: INSUFFICIENT (2026-07-27, jobs 135956/135957)

Filtered phi=0.20 probe (ELMesoFilter=Yes, 16 bins, t=100): fluctuation
suppression works — Tgran peak 1.6e-3 vs 9.7e-3 unfiltered (6x), onset
delayed ~2x, compensator drift bounded at 1.2e-5 — but the mean settling
enhancement persists: up_z grows all run (not stationary at t=100), TAVG
U_RZ = -1.67e-2, U/U0 = 2.14 vs RZ 0.354. The diagonal/checkerboard
concentration modes, untouched by the x/y lane filter, carry the
instability. Filtered U0 = -7.830e-3 (job 135957; the filter removes
part of the self-co-flow bias). Remaining filtered phi points NOT
submitted. Proceeding to option 2: shrink L/d_p from 20 to 10 by
doubling d_p to 0.1 in the same box/mesh (delta/h = 8, Re_p ~ 0.17,
N(0.20) = 382) — the box then excludes the unstable long wavelengths.
Probe-first: phi=0.20 + U0 at d_p=0.1 before any sweep.

## Small-box probe (option 2) — verdict: FAIL, instability persists at L/d_p=10 (2026-07-27, jobs 136121/136122)

d_p=0.1 twin pair: U0 run (job 136121) gave U0 = -3.068e-2
(frame-corrected from up = -3.449e-2, uf_super = -3.810e-3). phi=0.20
run (job 136122) clean to t=100 (no NaN, worst NP 7.2e-18):

- TAVG (tail window): up_z = -2.554e-2, uf_super_z = +1.259e-2
  -> U_RZ = up_z - (uf_super,z + phi*up_z) = -3.302e-2.
- RZ target: U0 * 0.8^4.65 = -1.087e-2. Ratio 3.04x -- SAME
  enhancement band as the L/d_p=20 box (2.6-3.3x).
- |U_RZ|/|U0| = 1.08: the suspension settles FASTER than an isolated
  particle. Not residual hindrance -- net collective enhancement.
- Time series: up_z grows non-monotonically all run (plume formation /
  breakup cycles: -1.5e-2 at t=10, -1.3e-2 at t=50, -3.1e-2 at t=90);
  Tgran ramps to 3.0e-3 by t=100 with no plateau.

Halving the box (L/d_p 20 -> 10) did NOT stabilize the homogeneous
state: the fastest-growing mode still fits. Combined with the lane
filter's failure (diagonal modes) this exhausts the cheap suppression
routes. DECISION: option 4 -- V2 is reframed against the
cluster-induced-settling literature (2-4x enhancement is the documented
behavior of two-way coupled point-particle sedimentation in fully
periodic domains); RZ hindered settling is a wall-bounded/homogenized
result our periodic box does not represent. The V2 infrastructure gates
(momentum conservation, Newton pairs, compensator drift) all PASS and
stand; the RZ ratio rows stay FAIL-as-measured with this physics
explanation. Remaining small-box phi points NOT submitted.

## Post-fix re-run (sweep 4) — periodic Q2 coupling restored (staged 2026-10-02, NOT yet submitted)

### Why this re-run exists

On 2026-08-03 the DNS campaign (branch `feature/dns-validation`, D1.1b)
proved that the Q2/E013 communicator held ZERO periodic DOF pairs:
`myVERTLINK` was built with a non-periodic octree search, so every
"fully periodic" CFD run — including ALL of sweeps 1–3 above and the
lane-filter / small-box probes — was actually solved with traction-free
(free-breathing) faces. Fixed on that branch (commits `6e26474d`,
`59dc4758`, `9cc4769f`; DNS datasheet rows `d11_periodicity_absent`,
`d11b_periodic_comm`, `d11_rh_collapse`). The fixed binary was already
re-verified on the tier-2 `momentum_conservation` case (DNS row
`d11_el_momentum_recheck`, job 137542: Newton-pair 5.08e-21, ELEMINT
drift_rel 5.32e-8, 702 PERIODIC_COMM pairs on QBOX9).

Consequence for THIS case: the sweep-3 FAIL (enhancement 2.6–3.3×) and
the option-4 reframe ("a periodic box has no suppression mechanism")
were measured in a box that was not periodic. Traction-free faces admit
free in/outflow of fluid momentum and are a candidate driver of the
plume instability. The reframe is therefore UNPROVEN until the sweep is
repeated on a genuinely periodic box.

Hypothesis branches (decide from the data, never re-tune):
- (a) enhancement persists at the same band → the CIS reframe stands,
  now on solid ground; annotate the sweep-3 rows accordingly.
- (b) hindrance appears, or the enhancement changes materially →
  re-attribute the FAIL to the traction-free faces; the RZ fit becomes
  quotable again (or a new physics row is written).

What the DNS campaign found in the meantime (context, NOT a gate):
D3.2 settled N=20..120 clouds in a WALLED 6d column — stable hindered
settling, confined exponent n = 4.57 (row `d32_phi_ladder`); in a 12d
column the same clouds settle ENHANCED, growing with N (row
`d32_wide_attribution`: the 6d exponent is a confinement law driven by
collective return flow). The DNS campaign never ran the periodic RZ
box: moving particles crossing periodic faces are unsupported in the
FBM classifier (row `d31_periodic_indicator`). So DNS supports the
"walls stabilize" reading but does not test the periodic claim; this
re-run is the only direct test. (The d31 limitation is FBM-specific:
the EL kernel/halo layer wraps particles — rows `pe_z_periodic_wrap`,
`unit_periodic_contact`, Newton pairs through periodic images.)

### Configuration — identical to sweep 3, binary is the ONLY change

- Worktree `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-el`, branch
  `feature/el-rz-postfix` (forked from `feature/dns-validation`
  `d8751c7e`, which contains the EL code, the periodic fix, the
  partition guard `CheckFaceClaimDecode`, the MUMPS-on-periodic refusal
  and the pe/CFD stepsize guard; pe pin `ce3bb5f`, whose
  `HardContactEulerLagrange.h` / `el_optional_api.h` /
  `setup_el_terminal_velocity.h` are unchanged since the EL merge pin
  `e14c6cd`).
- Build `build-el-rz-postfix/` in the worktree: same options as the
  certified `build-el-phase2-pe-gcc14` (Ninja, Release, USE_PE=ON,
  PE_USE_MPI=ON, serial mode OFF, CGAL/HYPRE/MUMPS/BOOST/OPENMESH OFF,
  ENABLE_FBM_ACCELERATION=ON, SED_BENCH=OFF, gcc 14.3.0 + OpenMPI 4.1.6
  wrappers `/sfw/openmpi/gcc14.3.x/4.1.6/ucx-threaded-noverbs/bin`,
  default pe_CONSTRAINT_SOLVER = HardContactEulerLagrange, i.e. NO
  solver define — the DNS trees pass HardContactAndFluid, EL must not).
  Binary: `build-el-rz-postfix/applications/q2p1_el_pipeflow/q2p1_el_pipeflow`.
- Rundirs (staged by `stage_postfix_sweep.sh` in this directory into
  `build-el-rz-postfix/applications/`): `q2p1_el_rundir_v2pf_smoke`
  (φ=0.20, 200 steps), `..._v2pf_phi05/10/15/20` (50000 steps to
  t=100), `..._v2pf_u0` (single particle, 5000 steps to t=10).
  VERIFIED 2026-10-02: every `_data/q2p1_param.dat` is byte-identical to
  the corresponding sweep-3 rundir
  (`build-el-phase2-pe-gcc14/applications/q2p1_el_rundir_v2mfx_*`,
  jobs 135866–135870) and every `example.json` is semantically identical
  (same keys and values; only whitespace differs).
- Prescription unchanged: `ELFluidGravity = No`, `ELMomentumFix = Yes`,
  legacy convective form, `Periodic = Yes`, QBOX9 3×3×3 (axis-uniform,
  3 ranks per periodic axis — satisfies the owner's periodic-partition
  rule), np = 28. Pres coarse solver type 4 / Velo type 1 (not MUMPS).

### Gates, in order (probe-first)

- **G0 — build twin.** Stage tier-2 `momentum_conservation` from the
  new build (`ninja q2p1_el_pipeflow_tier2_momentum_conservation_stage`)
  and run it (np=28, 60 steps, ~20 s). Compare against job 137542's
  `run_mom_periodicfix.log` (Aug-3 fixed binary, same tree family):
  702 PERIODIC_COMM lines, EL_FEEDBACK_CONSERVATION residual ≤ 1e-10,
  EL_MOMENTUM_ELEMINT drift_rel < 1e-5. Expect the EL diagnostic lines
  to agree to all printed digits (the pe pin moved e14c6cd → ce3bb5f,
  but the EL solver path is untouched); any drift is investigated
  BEFORE the sweep.
- **G1 — smoke** (`v2pf_smoke`): the sweep-3 smoke gates (N=3056,
  achieved φ, `fully periodic PE decomposition` banner, ⟨ε⟩=0.80, no
  OOB/clipping, max_overlap ≤ 1% d_p) PLUS 702 PERIODIC_COMM lines in
  stdout. Record wall-time/step (sweep 3: ~1.7 s/step at np=28 on nx).
- **G2 — production** (`u0` + four φ): the sweep-3 infrastructure table
  (NaN, per-step mismatch_z, cumulative drift_z, Newton-pair_z,
  max_overlap, ⟨ε_f⟩ vs 1−φ) must hold as before; then the physics.

### Metric and verdict bookkeeping

- Same metric as sweep 3: U_RZ = ⟨u_p,z⟩ − (uf_super,z + φ⟨u_p,z⟩) from
  `EL_MEAN_SLIP_TAVG`, U₀ from the `u0` run, late window t=[50,100] and
  early window t=[5,15]; plus the instability markers (Tgran history,
  proximity-pair count, ⟨u_p,z⟩ oscillation).
- Report side by side with the sweep-3 tables above (same φ, same
  windows). Also report U₀ itself: the +16% co-flow bias over algebraic
  Di Felice was measured traction-free and may move.
- New datasheet rows carry the suffix `_postfix` and cross-reference
  `d11_periodicity_absent`; the sweep-3 rows stay, annotated. Never
  rewrite history.

### Submission (owner decision — nothing submitted as of staging)

```bash
EL=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-el
export EL_BIN=$EL/build-el-rz-postfix/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
R=$EL/build-el-rz-postfix/applications
$EL/tools/el_slurm_submit.sh $R/q2p1_el_rundir_v2pf_smoke v2pf_smoke 01:00:00 short
# after G1 passes:
$EL/tools/el_slurm_submit.sh $R/q2p1_el_rundir_v2pf_u0    v2pf_u0    02:00:00 short
for p in 05 10 15 20; do
  $EL/tools/el_slurm_submit.sh $R/q2p1_el_rundir_v2pf_phi$p v2pf_phi$p 1-23:00:00 long
done
```
Sweep-3 wall times (np=28, nx): φ=0.05 10.6 h, 0.10 16.3 h, 0.15 20.6 h,
0.20 24.2 h, U₀ 37 min. Each rundir stays ≈30 MB (OutputFreq 10 t.u.,
BackUpFreq 1000). `/data/warehouse17` was at 95% (249 GB free) when staged.

### G0 result — build twin PASS (2026-10-02, job 147442, nx-05, 14 s)

Tier-2 `momentum_conservation` on the worktree binary vs job 137542
(`run_mom_periodicfix.log`, Aug-3 fixed binary): 702 PERIODIC_COMM
lines on both; EL_MOMENTUM_ELEMINT (60 lines), EL_FEEDBACK_CONSERVATION
(60), EL_NEWTON_PAIR (60), EL_VOLUME_CONSERVATION (60) and EL_MEAN_SLIP
(120) are BITWISE IDENTICAL to the reference. Final step: drift_rel
5.319670E-08, feedback residual 5.082198E-21, Newton-pair mismatch
1.3e-29/-7.9e-30/0. The pe pin move e14c6cd -> ce3bb5f is a no-op on
the EL solver path, as expected from the unchanged pe files. Log:
`build-el-rz-postfix/applications/q2p1_el_pipeflow/run_g0_twin.log`.

### G1 — smoke SUBMITTED (job 147443, short partition, phi=0.20, 200 steps)

Binary `build-el-rz-postfix/applications/q2p1_el_pipeflow/q2p1_el_pipeflow`
via `tools/el_slurm_submit.sh` with EL_BIN set. Gates as listed above.
Production (U0 + four phi) NOT submitted — owner decision.

### G1 result on the HEAD binary — gates PASS, but a pe-pin confound surfaced (2026-10-02, job 147443, nx-05, 5:05 wall)

All listed gates hold: 702 PERIODIC_COMM pairs, `fully periodic PE
decomposition` banner, 3056 spheres at achieved φ = 0.2000147,
⟨ε_f⟩ = 0.7999853 (identical to sweep 3 to all digits), volume
conservation 1.1e-15, Newton-pair ≤ 6e-19, EL_FLUID_PAIR mismatch
1.1e-12, no NaN/OOB/clipping, ≈1.3 s/step incl. setup. At t = 0.4:
⟨u_p,z⟩ = −2.1649e-3 vs pre-fix −2.3542e-3 (8% slower settling),
uf_super,z = +8.391e-4 vs +9.519e-4, Tgran 3.17e-6 vs 4.84e-6 — the
first periodic-vs-traction-free physics signal, NOT a verdict (200 steps).

Anomaly: `EL_CONTACT_STATS ncontacts` = 0 at every sample, where sweep 3
had 5571–5574 from step 1 (seeded configuration, lubrication OFF).
Root cause traced in pe, not in the periodic fix:
- At the EL-era pins (e14c6cd … 1fef6b8) `SphereBase::calcBoundingBox`
  padded every sphere AABB by `lubricationThresh` = 1e-2 regardless of
  the lubrication switch, and `MaxContacts::collideSphereSphere` (with
  `PE_LUBRICATION_CONTACTS` defined through `pe/core/CollisionSystem.h`
  → `HardContactLubricated.h`) emitted a *lubrication-type* contact for
  every pair with gap < 1e-2. HardContactEulerLagrange's contact mask
  filters by ownership only, so those proximity pairs entered the PGS
  with positive gap: inert at these velocities (gdot_n + dist/dt > 0 ⇒
  separating ⇒ no reaction), but counted in `getNumberOfContacts()`.
  Hard contacts were and are generated only below `contactThreshold`
  = 1e-8 (pe/core/Thresholds.h); max_overlap = 0 in every old run.
- At the current pin ce3bb5f (pe PR #26, commits 4273f98/245b53b,
  2026-08-20) the padding is `lubrication::aabbPadding(r)` = 0 when
  lubrication is off, and the fine detector takes a fast path
  (`!contactGenerationEnabled()` → `dist < contactThreshold` only). So
  `ncontacts` now counts hard contacts only. Residual physics
  difference: a pair that closes its gap within one step is stopped
  pre-emptively by the old proximity constraint but caught after a
  micro-penetration (≤ one step's displacement, then erp-corrected) by
  the new pin. Sweep 3 used `ncontacts` as its clustering marker
  (proximity pairs ~240 → ~550 at φ = 0.05); that marker is gone on the
  HEAD binary.
- Frozen-box twin, lubrication OFF (`kroupa_shear` smoke deck of job
  135996, φ = 0.20, 500 steps, prescribed shear ⇒ no fluid solve ⇒
  periodic fix irrelevant), HEAD binary, job 147445: EL_MEAN_SLIP (40),
  EL_FORCE_BUDGET (500), EL_SUSP_STRESS (20) and the dp/expected fields
  of EL_NEWTON_PAIR BITWISE identical to the certified log; Tgran
  identical; differences only in `ncontacts` (0 vs 5526–5722), the
  Newton audit roundoff (|mismatch| ≤ 7.4e-18 vs 6.3e-18) and
  EL_VOLUME_CONSERVATION `rel_error` digits (both ≤ 1.8e-15, deposited
  and expected identical). The proximity contacts carried no impulse.

Decision: run sweep 4 on a binary with the EL-era pe so that the ONLY
change vs sweep 3 is the periodic coupling, and the clustering marker
keeps its meaning. That binary is the state the DNS campaign certified
on 2026-08-03 (row `d11_el_momentum_recheck`): FF `23144a01`
(= `9cc4769f` + a RUNBOOK edit) with pe `1fef6b8` (contact-relevant
files identical to e14c6cd). Rebuilt reproducibly in a second worktree:
`/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-el-aug3` (detached at
23144a01, pe 1fef6b8), build `build-aug3/` with the same cmake recipe;
binary `build-aug3/applications/q2p1_el_pipeflow/q2p1_el_pipeflow`.
G0 on that binary (job 147448, 14 s): 702 PERIODIC_COMM, all six EL
streams incl. EL_CONTACT_STATS BITWISE vs job 137542. The HEAD-pin
binary (`FeatFloWer-el/build-el-rz-postfix`) stays as the cross-check;
if the sweep-4 verdict is later re-run on it, expect `ncontacts` = 0
until genuine collisions occur and max_overlap no longer identically 0.

Submission for sweep 4 therefore uses
`EL_BIN=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-el-aug3/build-aug3/applications/q2p1_el_pipeflow/q2p1_el_pipeflow`
with the commands listed above (same rundirs; the HEAD-binary smoke
lives in `q2p1_el_rundir_v2pf_smoke`, the Aug-3-binary smoke in
`q2p1_el_rundir_v2pf_smoke_aug3`).
