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
| 11.4 (L2) | — | +3.6% | +3.4% | +4.4% |
| 23.9 (L3) | −3.2% | +0.6% | +0.3% | +0.8% |
| 49.1 (L4) | −5.2% | −1.7% | −2.1% | −1.6% |

Cross-case structure (datasheet `d13_matrix`): the E2–E4 ladders are
nearly identical (within ~0.5 pp at every level) despite an 8× range in
Re, and the L3→L4 shift is uniform (−2.1…−2.4 pp) across **all four**
cases including E1 — the spatial error at these resolutions is
**geometry-dominated** (interface representation of the sphere), not
flow-regime-dependent. E1's overall offset is its sim-vs-experiment
reference gap (§4), not a different discretization behavior.

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

## 3. Timestep: accuracy-bounded only (corrected 2026-08-02)

All dt-dependence below is from SYNCED runs (deck TimeStep == json
stepsize_; the pe_stepsize_mismatch contamination is resolved and a
fatal guard now enforces equality).

- **There is NO stability floor.** The earlier two-sided-window finding
  is refuted: synced dt = 0.25 ms is fully stable (0 sawtooth warnings)
  at rho_p/rho_f = 1.17, and the DKT-ratio probes (1.10, 1.02) are
  clean at dt = 0.5-1.0 ms. The instability was the 4x PE/CFD desync.
  The sawtooth watchdog stays as cheap insurance.
- **Genuine temporal term, sub-linear**: E4 L3 synced ladder
  +0.81/+1.41/+1.94% at dt = 1.0/0.5/0.25 ms - near-equal increments
  per halving favor apparent order ~0.5-1. Global fit: T(1 ms) =
  -1.5..-2.4 pp (order 1..0.5 reading), S(L3) = +1.9..+3.1 pp,
  S(L4) = -0.5..+0.9 pp (converged within uncertainty). E1: +0.42 pp
  per halving (2 points).
- **Reading**: smaller dt moves peaks toward the (positive) spatial
  error; the good raw numbers at dt = 1 ms ride on partial S/T
  cancellation. Choose dt by matching the T column against your error
  budget - not by stability.
- **Config rule (fatal-guarded)**: deck SimPar@TimeStep == json
  stepsize_. The guard aborts step 1 on mismatch.

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
