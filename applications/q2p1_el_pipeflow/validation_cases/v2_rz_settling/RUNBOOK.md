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
