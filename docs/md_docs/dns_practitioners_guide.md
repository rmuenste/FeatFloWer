# DNS (FBM) practitioner's guide — v2.2

Status: **v2.2**, 2026-08-22. v2.2 rewrites §6: the "frictional stall"
of the DKT benchmark was an artifact of the `hcaf_angvel_reset` solver
defect (angular velocity zeroed every step; fixed at libs/pe ≥ de855b6)
— on the repaired binary the frictional run tumbles
(`d23_omegafix_rerun`). v2.1, 2026-08-18. v2 (2026-08-11) added D3.1 (random-array
drag law, closed) and the production-resolution probes to the v1 scope
(D1 metrology, D1.1 Hasimoto, D1.2 noise floor, D2.1 crossover rule,
D2.3 DKT). v2.1 applies the SI-hardening audit: the a_eff **sign
erratum** (a_eff = a − 0.14h, `d11_aeff_sign_erratum`), the
**transient-bias revision** of the dilute finite-Re array numbers
(`d31_transient_bias`, `d31_t12_measured`), the multiplicative
matched-Re basis (`d31_matched_re_operation`), and the completed
nine-check refinement set (`d31_l5_re9_p010`). Every number below
traces to a row in `dns_validation_datasheet.md`. Scope: spheres,
Re ≈ 0.003–32, serial-PE mode. Multi-particle collective behavior (D3),
non-spherical particles (D4) and the shared lubrication closure (D2.2)
will extend this.

## 1. Build & configuration prerequisites (non-negotiable)

- `-DUSE_PE=ON` (+ `-DUSE_PE_SERIAL_MODE=ON` for few/large particles),
  and `pe_CONSTRAINT_SOLVER=pe::response::HardContactAndFluid` via
  `CMAKE_CXX_FLAGS`. The pe-master default (`HardContactEulerLagrange`)
  silently discards FBM forces — particles stay motionless with **no
  error message** (datasheet `d0_solver`).
- Deck and PE json must agree on fluid properties (`Prop@Density`,
  `Prop@Viscosity` ↔ `fluidDensity_`, `fluidViscosity_`) — nothing
  checks this (`d0_visc`) — and on the timestep: deck `SimPar@TimeStep`
  == json `stepsize_` is **fatal-guarded** at step 1 (the serial PE
  integrates with the json value; a mismatch warps every transient —
  `pe_stepsize_mismatch`). Numeric `Prop@` lines must not carry inline
  comments (list-directed read).
- Stage cases with `tools/stage_tencate_case.py` (template-clone +
  count-checked substitutions) rather than hand-editing decks.
- **Periodic runs** (added after D1.1b): require (a) a binary containing
  the periodic Q2 coupling fix (grep the log for `PERIODIC_COMM` pairs),
  (b) an **axis-uniform Cartesian partition with ≥2 ranks per periodic
  axis** — METIS partitions are invalid and the face decode aborts,
  (c) no MUMPS coarse solver (refused at init). Before the fix every
  "periodic" face in the code base was silently traction-free.

## 2. Spatial resolution (cells per diameter, D/h)

D/h counts **Q2 elements** (h = dvol^(1/3)); nodal spacing is h/2 —
halve literature values quoted for LBM/IBM node counts when comparing.

Measured peak-velocity errors vs ten Cate Table II printed experimental
ratios, dt = 1.0 ms:

| D/h | Re=1.5 | Re=4.1 | Re=11.6 | Re=31.9 |
|---|---|---|---|---|
| 11.4 (L2) | — | +3.6% | +3.4% | +4.4% |
| 23.9 (L3) | −3.2% | +0.6% | +0.3% | +0.8% |
| 49.1 (L4) | −5.2% | −1.7% | −2.1% | −1.6% |

Cross-case structure (`d13_matrix`): the E2–E4 ladders are nearly
identical (within ~0.5 pp) despite an 8× range in Re, and the L3→L4
shift is uniform (−2.1…−2.4 pp) across all four cases — the spatial
error is **geometry-dominated** (interface representation), not
flow-regime-dependent.

The periodic-array benchmark (D1.1, closed) makes the mechanism
quantitative: after correcting each level for its measured
indicator-volume error, the drag ladder is **first-order in h** with an
effective-radius picture a_eff ≈ a − 0.14h, converging to Hasimoto's
analytic value within −0.4…−0.6% (`d11_rh_collapse`,
`d11_aeff_sign_erratum`). The FBM interface behaves like a sphere
**narrowed** by ~0.14 cells — the discrete constraint under-enforces
no-slip between velocity nodes, so flow effectively penetrates ~0.14h
into the nominal solid. (v2 recorded the sign as +0.14h; the drag
sensitivity dlnK/dln a ≈ +1.8 > 0 combined with the measured deficits
forces the minus sign. All magnitudes stand.)

Rules of thumb (this flow class):

- **D/h ≈ 24 is the workhorse**: 1%-class peaks for Re ≈ 4–32 (partly
  by S/T error cancellation at dt = 1 ms — see §3).
- **D/h ≈ 49 is spatially converged at Re ≈ 32** (residual ≈ −0.1 pp
  after removing the temporal term).
- **Low Re needs the fine mesh for the right reason** (long-range wall
  retardation); the residual ~−5% vs the E1 *experiment* is present in
  ten Cate's own LBM (−5.6%) — do not chase it with resolution.
- D/h ≈ 11 is smoke-test grade (+4%).
- D/h ≲ 1.5 spheres are **hydrodynamically transparent** (`d11_ctrl`) —
  sub-grid particles need a different method (EL), not a finer gate.

## 3. Timestep: accuracy-bounded only

All dt-dependence from SYNCED runs (deck TimeStep == json stepsize_).

- **There is NO stability floor** (`dt_stability_refuted`): synced
  dt = 0.25 ms is fully stable at ρ_p/ρ_f = 1.17 and the DKT-ratio
  probes (1.10, 1.02) are clean. The earlier "floor" was the 4× PE/CFD
  desync. The sawtooth watchdog stays as cheap insurance.
- **Genuine temporal term, sub-linear**: E4 L3 synced ladder
  +0.81/+1.41/+1.94% at dt = 1.0/0.5/0.25 ms (apparent order ~0.5–1).
  Global fit (tools/tencate_error_decomposition.py, synced data):
  T(1 ms) = −1.5/−1.1 pp for E4 and −0.8/−0.6 pp for E1 (p=1/p=2
  readings; order not pinned), S(L4) = −0.0…−0.5 pp — converged within
  uncertainty. Smaller dt moves peaks toward the (positive) spatial
  error — the good dt = 1 ms numbers ride on partial cancellation.
  Choose dt against your error budget, not against stability.
  (A previously quoted upper bound of −2.4 pp was a pre-sync carryover
  and does not reproduce from the tool — corrected 2026-08-13.)

## 4. Force noise floor (grid-crossing) — first order in h

A moving interface re-classifies indicator DOFs as it crosses cells;
the resulting force jitter is intrinsic to FBM (`d12_prescribed_ladder`,
constant-V protocol):

| D/h | 11.4 | 23.9 | 49.1 |
|---|---|---|---|
| rms(F)/mean(F) | 3.8% | 1.8% | 0.9% |

Halves per level (**O(h)**), matches the free-fall plateau numbers.
Consequences: (a) expect ~1%-class force jitter at the D/h ≈ 24
workhorse — tolerances tighter than that measure noise, not physics;
(b) trajectory jitter of 1e-5-class velocities is normal, not a defect.

## 5. Near-contact: the gap ≲ 2h crossover rule (D2.1)

Measured against Brenner's exact wall-approach solution (constant-V
protocol, Re = 0.78; figure `dns_figures/d21_brenner_crossover.png`):

- The resolved FBM force follows the lubrication divergence to within
  ~10% down to a gap of **2–3 cells**, and degrades rapidly below
  ~1.5 cells (−20% at 1.1–1.7h).
- In cell units the departure curves of D/h = 24 and 49 nearly
  **collapse** — the near-wall error is a universal function of gap/h.
  A sub-grid correction should (a) take over at gap ≲ 2h, and (b) can
  in principle apply the inverse of the collapse curve rather than a
  hard switch (D2.2 design input).
- This replaces the ad-hoc `5·h_min` threshold (2.5× over-conservative)
  and the SRR hard-coded 0.07. Note HCAF itself has **no lubrication**
  (API no-ops) — today the sub-2h regime is governed by the contact
  model alone, which leads directly to §6.

## 6. Contact parameters are physics, not stabilizers (D2.3/DKT)

**v2.2 correction.** Earlier versions of this section reported a
resolution-independent "frictional stall": tilt frozen at ~4° under the
default friction (μ_s/μ_d = 0.1/0.05) at both D/h = 8 and 16, unlocking
only when friction was zeroed. That stall was an **artifact of the
`hcaf_angvel_reset` solver defect** (datasheet row): the certified
binaries zeroed every particle's angular velocity each step, making
rolling kinematically impossible — so *any* tangential friction froze
the contact. On the repaired binary (libs/pe ≥ de855b6) the identical
frictional run proceeds through tumble onset: same drafting phase and
kiss time (t = 18.10), then tilt growing exponentially — 16.2° by
t = 25, doubling every 3 t.u. — versus 4.2° frozen before the fix
(`d23_omegafix_rerun`; the stall's apparent resolution-independence was
the defect's, not the physics').

What survives, measured on consistent binaries: friction *modulates*
the tumble rather than suppressing it — the frictionless run reaches
22.3° where the frictional one reaches 16.2° at the same t = 25. The
physical argument stands: spheres in liquid contact ride a lubrication
film that transmits little tangential traction, so a dry frictional
contact overstates tangential coupling.

- For DKT-class problems (mobile particle-particle contacts in liquid):
  prefer `staticFriction_`/`dynamicFriction_` ≈ 0 in the json (keys are
  config-driven since libs/pe `dcc35f2`), or the lubricated-contact
  add-on once D2.2's validation ladder closes — but on a repaired
  binary this is now a quantitative modeling choice, not the difference
  between tumbling and not tumbling.
- Keep friction where it is physical (dry granular contacts, resting
  beds).
- Corollary of §5: everything below ~2 cells of gap IS the contact
  model — choose its parameters as physics, not for stability.
- If your binary predates the fix (libs/pe < de855b6), all rotation-
  coupled results are suspect: check for `w = Vec3(0,0,0)` in
  `HardContactAndFluid::integratePositions`, or simply verify that
  `DNS_PART_STATE` angular velocities are not exactly zero for all
  time under nonzero torque.

## 7. Reference discipline (ten Cate specifics)

- Peak gates use the **printed Table II ratios** (0.947/0.953/0.959/
  0.955 × Table I u_∞), not digitized curve peaks (E1/E2 digitized are
  +3.4/+3.9% too fast, `tc_ref_audit`).
- Peak timing is plateau-argmin noise at low Re; use shape RMS with the
  ~40 ms PIV time-origin caveat.
- At Re ≲ 2 the honest reference band includes the paper's own
  simulations (Table II S-series).
- Hasimoto comparisons: F = f·V_cell (eq. 2.14) against the
  **superficial** (whole-cell mean, Darcy) velocity — pinned from the
  paper (`d11_convention_pinned`); the interstitial variant differs by
  O(φ) and is the wrong comparison.
- Beetstra comparisons: the correlation is for the DRAG PART only,
  F_d = (1−φ)·F_total (their p. 490 force split — our momentum-balance
  force is F_total), normalized by superficial U; Re uses superficial U
  too (`d31_convention_pinned`). Third convention pin of the campaign:
  **pin the reference's own definitions from its own text before
  believing any discrepancy.**

## 8. Prescribed-motion probes (protocol)

For controlled-velocity measurements without new solver code: heavy
sphere (ρ_p/ρ_f = 1e6) + `gravity_ = [0,0,0]` +
`initialParticleVelocity_` (+ `benchUseConfigStart_` to place it) —
velocity stays constant to ~1e-5 relative while the FBM force is logged
(`d12_prescribed_smoke`). Twin-certified no-op when the keys are absent.

## 9. The random-array drag law (D3.1, closed)

F_d*(φ, Re) measured on N=16 random fixed-sphere arrays (RSA, 0.05d min
gap, periodic cell with explicit image spheres), 3 seeds/point, vs the
Beetstra 2007 correlation (eqs. 6/17). Values are **steady-state**
(`tools/d31_steady_extrapolation.py`; the v2 table read the dilute
finite-Re runs at their t=4 endpoint, which was up to 13% transient-
biased — `d31_transient_bias`, verified by direct t=12 reruns,
`d31_t12_measured`):

| φ | Re≈0 | Re≈9 | Re≈27 |
|---|---|---|---|
| 0.05 | 1.90±0.07 | 2.25±0.14 | 2.78±0.21 |
| 0.10 | 2.48±0.17 | 2.88±0.23 | 3.57±0.34 |
| 0.20 | 4.38±0.20 | 4.95±0.20 | 6.07±0.32 |
| 0.30 | 7.54±0.20 | 8.46±0.25 | 10.57±0.32 |

Practitioner rule from that audit: weakly forced runs relax with
τ ≈ ρU_ss/(f(1−φ)) — for dilute cells this is O(1–2) time units, so
budget ≳ 6τ of integration or gate on the endpoint drift
|U(t_e)−U(t_e−1)|/U(t_e) < 5×10⁻³ before reading off a steady value.

What the deviations taught (each attributed, none unexplained):

- **Resolution**: L4→L5 shifts small and uniform — all nine refinement
  checks now **measured** (`d31_l5_re9_p010`): Stokes
  +2.3/+2.3/+2.9/+3.6% across φ, finite-Re +1.3…+1.7% — the
  a_eff = a − 0.14h picture again; L4 (D/h 13–24) is production-grade
  for array drag.
- **Finite-N**: the +2…+8% converged Stokes offset vs the correlation
  is the N=16 box talking to its own periodic images; at Beetstra's own
  N=54 the FBM value reproduces the correlation to **+0.5%** (L5,
  steady balance-force basis, `d31_n54_verdict`). Rule: for
  quantitative closure comparisons use N≈50+, or expect a few-%
  finite-N offset at N=16.
- **Microstructure dispersion at Re ≳ 20**: single-configuration drag
  scatters ±10–30% around the ensemble mean (vs 3–8% at Stokes) —
  inertia amplifies configuration sensitivity ~3×
  (`d31_corner_discriminator`). This dispersion is physics, not error:
  mean-field closures integrate it away; a DNS-informed EL closure can
  carry it as a fluctuation model (D5 input).
- **The FBM indicator is NOT periodic** (`d31_periodic_indicator`):
  spheres straddling a periodic face lose the part outside the box.
  Fixed arrays: emit explicit image spheres (generator does this).
  Moving particles crossing periodic faces remain UNSUPPORTED — flag
  for any wrapping-particle run.

**Production resolution (D/h = 3–4)** — single-sphere and array-level
evidence agree (`d11_coarse_probes`, `d31_l2_production`,
`d31_matched_re_operation`): mean drag is **84–87% of converged**
(matched-Re basis: deficits −12.1…−13.9%, Stokes → Re≈27, mildly
deepening with Re) → apply one correction factor **1.15–1.20**
(equivalently a_eff = a − 0.14h). The matched-Re comparison is
multiplicative — scale the coarse value by the correlation ratio
F_corr(Re_fine)/F_corr(Re_coarse), not by an additive increment
(the operation matters ~2 pp at Re≈27). Per-particle
forces carry ~8–15% rms grid noise there, and everything within ~2
cells of contact (= half a diameter at D/h 3–4) is the contact model's
job — tune its parameters as physics (§6).

## 10. Measured job costs (32 ranks unless noted, nx nodes)

| Config | steps | wall time |
|---|---|---|
| L2, dt=1 ms | 1300 | ~4 min |
| L3, dt=1 ms | 1300–4300 | 16–55 min |
| L4, dt=1 ms | 1300–4300 | 2–7 h (60G) |
| Hasimoto L3 (28 rk) | 400 | ~20 min |
| Hasimoto L5 (28 rk, 2 nodes) | 400 | ~6 h (93G/node peak) |
| DKT D/h=8 (109 rk) | 5000 | ~6 h |
| DKT D/h=16 (109 rk, 2 nodes) | 5000 | ~30 h |

## 11. Open items feeding v3

- D3.2 wall-bounded hindered settling — ladder CLOSED 2026-08-15
  (`d32_phi_ladder`: confined exponent n ≈ 4.3–4.6 in the 6d walled
  column vs unbounded Rowe/RZ 2.7–3.0 at swarm Re 51–77; collapse rms
  2.3%); v3 write-up pending, wall-vs-cloud attribution open
  (wider-column discriminator identified).
- D2.2 sphere–sphere approach + shared lubrication closure (uses §5's
  collapse curve and §9's dispersion finding).
- Periodic support for MOVING particles crossing faces (indicator +
  PE) if any wrapping-particle case is needed.
- D1.2 spectral characterization; E2/E3 dt points; Tenneti 2011
  comparison of the §9 surface.
- Non-spherical (D4) and the D6 capstones.
