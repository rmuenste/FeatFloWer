# DNS (FBM) practitioner's guide — v0 draft

Status: **draft**, 2026-08-01, after the D1.3 first pass (ten Cate E1–E4,
levels L2–L4, dt 1.0/0.5/0.25 ms). Every number below traces to a row in
`dns_validation_datasheet.md`. Scope so far: a single settling sphere,
Re ≈ 1–32, ρ_p/ρ_f ≈ 1.17, serial-PE mode, benchSym quarter-box mesh.
Multi-particle, non-spherical and contact-dominated guidance arrives with
D2–D4.

## 1. Build prerequisites (non-negotiable)

- `-DUSE_PE=ON` (+ `-DUSE_PE_SERIAL_MODE=ON` for few/large particles),
  and `pe_CONSTRAINT_SOLVER=pe::response::HardContactAndFluid` via
  `CMAKE_CXX_FLAGS`. The pe-master default (`HardContactEulerLagrange`)
  silently discards FBM forces — particles stay motionless with **no
  error message** (datasheet `d0_solver`).
- Deck and PE json fluid properties must agree (`Prop@Density`,
  `Prop@Viscosity` ↔ `fluidDensity_`, `fluidViscosity_`); nothing checks
  this for you (datasheet `d0_visc`). Numeric `Prop@` lines must not
  carry inline comments (list-directed read).
- Stage cases with `tools/stage_tencate_case.py` (template-clone +
  count-checked substitutions) rather than hand-editing decks.

## 2. Spatial resolution (cells per diameter, D/h)

Measured peak-velocity errors vs ten Cate Table II printed experimental
ratios, dt = 1.0 ms:

| D/h | Re=1.5 | Re=4.1 | Re=11.6 | Re=31.9 |
|---|---|---|---|---|
| 11.4 (L2) | — | — | — | +4.4% |
| 23.9 (L3) | −3.2% | +0.6% | +0.3% | +0.8% |
| 49.1 (L4) | −5.2% | — | — | −1.6% |

Rules of thumb (this flow class):

- **D/h ≈ 24 is the workhorse**: 1%-class peaks for Re ≈ 4–32 — but see
  §3; part of that accuracy is spatial/temporal error cancellation.
- **D/h ≈ 49 is spatially converged at Re ≈ 32** (residual spatial error
  ≈ −0.1 pp after removing the temporal term).
- **Low Re needs the fine mesh for the right reason**: at Re = 1.5 the
  long-range wall retardation is under-resolved at D/h = 24 (result too
  fast by ~2 pp relative to converged). The converged answer (ratio
  0.897) matches ten Cate's own LBM (0.894, +0.3%) — the residual ~−5%
  vs the *experiment* is a sim-vs-experiment gap present in the original
  paper itself (their S1: −5.6%). Do not chase it with resolution.
- D/h ≈ 11 is smoke-test grade only (+4.4%).

## 3. Timestep: the window is TWO-SIDED

The coupling (one explicit CFD↔PE force exchange per step) has both an
accuracy ceiling and a **stability floor** in dt at fixed mesh:

- **Accuracy (upper bound)**: the temporal term at dt = 1 ms is ≈ −1 pp
  (slows the particle); halving dt recovered +0.47 pp at L3 and +0.71 pp
  at L4 (E4), +0.43 pp at L3 (E1). Two-point evidence is consistent with
  a **first-order** coupling term (nominal CN order is not observed);
  order not yet pinned (see below).
- **Stability (lower bound)**: at E4/L3, dt = 0.25 ms is **unstable** — a
  growing oscillatory coupling mode from t ≈ 0.025 s, |v| reaching
  ±4 m/s, while the run "successfully finishes". Signature and regime
  (ρ_p/ρ_f = 1.17) match the added-mass instability of loosely-coupled
  partitioned FSI (Causin–Gerbeau–Nobile 2005 [in repo]: decreasing dt
  *aggravates* it). The floor is **mesh-independent** within L2–L3
  (dt = 0.25 ms unstable at both D/h = 11.4 and 23.9; dt ≥ 0.5 ms stable
  at every level tested): an absolute threshold in (0.25, 0.5] ms for
  ρ_p/ρ_f = 1.17, density-ratio-governed rather than CFL-like.
  Density-ratio dependence unmeasured — expect the floor to RISE as the
  ratio approaches 1.
- **Recommendation**: dt = 1.0 ms at D/h ≈ 24; dt = 0.5–1.0 ms at
  D/h ≈ 49 (sub-1% peaks at Re ≈ 32 with dt = 0.5 ms). Never reduce dt
  "for safety" without checking the trajectory for step-to-step sawtooth
  — plateau jitter should be ~10⁻⁵ m/s; the unstable run shows O(1)
  swings. A runtime sawtooth warning is a planned D0 hardening.

## 4. Reference discipline (ten Cate specifics)

- Peak gates use the **printed Table II ratios** (u_max/u_∞ = 0.947 /
  0.953 / 0.959 / 0.955 × Table I u_∞), not digitized curve peaks: the
  digitized E1/E2 peaks are +3.4%/+3.9% too fast (datasheet
  `tc_ref_audit`). Digitized curves remain the shape/timing reference.
- Peak *timing* comparisons are plateau-argmin noise for the flat low-Re
  curves; use shape RMS with a time-origin caveat (~40 ms PIV offset).
- At Re ≲ 2 the honest reference band must include the paper's own
  simulation values (Table II S-series), not just the experiment.

## 5. Measured job costs (32 ranks, nx nodes, this mesh family)

| Config | steps | wall time |
|---|---|---|
| L2, dt=1 ms, 1300 steps | 1300 | ~4 min |
| L3, dt=1 ms, 1300–4300 steps | 1300/4300 | 16–55 min |
| L3, dt=0.5 ms, E1 (8600 steps) | 8600 | ~1 h 45 |
| L4, dt=1 ms, 1300–4300 steps | 1300/4300 | 2 h – 7 h (60G mem) |
| L4, dt=0.5 ms, 2600 steps | 2600 | ~4 h 25 (60G mem) |

## 6. Open items feeding v1

- Pin the temporal order with a non-dyadic ladder inside the stable band
  (e.g. dt = 0.7 ms).
- Stability-floor dependence on density ratio (level dependence measured:
  none within L2–L3). The ATC/D6 regime is ~5% solids at similar density
  ratios — the floor matters there.
- E2/E3 level ladder points (only L3 run so far).
- Grid-crossing noise characterization (D1.2) and fixed-sphere array
  drag convergence order (D1.1).
