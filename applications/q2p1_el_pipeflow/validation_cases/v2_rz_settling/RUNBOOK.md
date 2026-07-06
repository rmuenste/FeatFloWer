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
