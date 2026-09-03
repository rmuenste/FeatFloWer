# D6.1 — Oberbeck spheroid drag in a periodic cell (CASE_SPEC)

Drafted 2026-09-03, while the D/h=16 rung (v24f) runs on Fritz. First
non-spherical validation family. Direct heir of D1.1 (Hasimoto) — same
fixture, same driving, same lattice-correction machinery — with the sphere
replaced by a **fixed prolate spheroid** and the gate replaced by Oberbeck's
anisotropic Stokes drag.

Owner decisions already taken:
- D6.1/D6.2 run on the **local** cluster (Fritz is for node-hungry work).
- Not gated on v24f: independent physics, independent domain, disjoint
  resources (see §9).

Prerequisite state (all in place): pe ellipsoid defects D-1/D-2/D-3 fixed at
pin 6971b13 (row `pe_ellfix_twin`); torque-path review
(`docs/md_docs/dns_torque_path_review.md`) — pipeline sound, moment-arm
caveat N-1 answered here by fixing the body.

---

## 1. What is measured, and against what

A prolate spheroid (semi-axes a > b = c, aspect ratio r_e = a/b,
eccentricity e = sqrt(1 - 1/r_e^2)) held **fixed** at the center of a
triply periodic cubic cell, flow driven by a uniform body force f (D1.1
convention). Measured: the steady FBM wrench (F, T) on the body.

Unbounded-Stokes reference (Oberbeck 1876; resistance-function form of
Kim & Karrila §3.3), with L = ln((1+e)/(1-e)):

```
F_parallel      = 6 pi mu a U * X^A,   X^A = (8/3) e^3 / [ -2e + (1+e^2) L ]
F_perpendicular = 6 pi mu a U * Y^A,   Y^A = (16/3) e^3 / [ 2e + (3e^2-1) L ]
```

(a = semi-MAJOR axis; both reduce to the sphere as e -> 0.)

Numbers for the primary aspect ratio r_e = 2
(e = 0.8660254, L = 2.6339158):

| quantity | value |
|---|---|
| X^A (axial) | 0.601970 |
| Y^A (transverse) | 0.689450 |
| **anisotropy Y^A/X^A** | **1.14532** |
| Perrin ratio f_par/f_sphere (equal volume, R_s = (a b^2)^(1/3)) | 0.95557 |
| Perrin ratio f_perp/f_sphere | 1.09443 |

(Optional second aspect ratio r_e = 3: X^A = 0.468155, Y^A = 0.575874,
ratio 1.23009 — recompute from the formulas at analysis time, do not trust
these two beyond 4 digits until re-derived in the analysis script.)

Torque identity: a body with three mutually perpendicular symmetry planes
(any ellipsoid) in uniform Stokes flow experiences **zero torque about its
center at every orientation**. With the body centered in a cubic lattice and
the axis in a lattice mirror plane this survives periodization. T ~ 0 is
therefore a free, orientation-independent diagnostic of the whole torque
path (the first DNS observable ever to exercise it) — reported as
|T| / (a * |F|).

## 2. Periodic (lattice) correction

Leading order, the simple-cubic correction depends only on the Stokeslet
strength, not the body shape: the D1.1 Hasimoto machinery carries over with
the sphere radius replaced by the orientation's hydrodynamic radius
R_h = F_iso / (6 pi mu U). Procedure (certified in rows `d11_rh_collapse`,
`d11_convention`):

- Hasimoto conventions PINNED as in D1.1: F = f * V_cell (whole-cell
  balance), U = superficial velocity v_0.
- Extract R_h per orientation by inverting
  F / (6 pi mu R_h U) = 1 / (1 - 1.7601 phi^(1/3) + phi - 1.5593 phi^2),
  phi = (4/3) pi R_h^3 / L_cell^3, solved self-consistently for R_h
  (same fixed-point the d11 analysis used).
- Gate R_h(parallel) against a * X^A and R_h(perp) against a * Y^A.

**Primary gate = the anisotropy ratio** R_h(perp)/R_h(parallel) vs
Y^A/X^A = 1.14532. The lattice correction and most of the resolution bias
cancel in the ratio; the per-orientation absolute values are the secondary
gate and carry the discretization calibration (§5).

Higher-order shape dependence (the phi and phi^2 terms see the body volume
and quadrupoles) is NOT shape-corrected here — it is controlled by the V0
sphere control and the resolution rung, and the long-axis image proximity
(2a/L_cell ~ 0.53 at r_e = 2, §4) is recorded as the known largest
systematic. If the V3/V4 residuals exceed the gate band, the first
suspect is this term, and the fallback is a larger-cell rung (QBOX9 mesh at
one level lower particle size, or the box-doubling variant d11 already
exercised).

## 3. Fixture (inherited from D1.1 verbatim)

- Mesh: QBOX9 periodic unit cube fixture (from
  `q2p1_el_pipeflow/tier2_cases/momentum_conservation/mesh/`), partition
  under `_mesh/QBOX9/` — **axis_uniform layout as required for periodic
  runs** (see the axis-uniform memory/datasheet rows; METIS is invalid
  here). Reuse the d11 partition as-is.
- `_data/MG.dat` required (copy from a bench rundir, 6 levels).
- Driving: `SimPar@ConstantForcing = 0,0,1d-2`, `Prop@Viscosity = 1d0`
  (+ matching `fluidViscosity_` in the json), Re ~ 3e-3.
- Body: **fixed** (d31 fixed-sphere mechanism), centered at the cell
  center. Fixing the body sidesteps review note N-1 (no periodic
  minimum-image on the torque moment arm) and pe defect D-2 entirely, and
  removes insertion transients.
- Lubrication: **OFF** (`lubricationEnabled_: false`) — the Kroupa model is
  sphere-pair only. Assert the `DNS_LUB` line is absent from the log.

## 4. Geometry and resolution

Volume-matched to the d11 sphere (r = 1.5) so the lattice-correction
magnitude and the K ~ 1.8 regime match D1.1 exactly:

- r_e = 2: a b^2 = 1.5^3 = 3.375 with a = 2b
  -> **b = c = 1.19055, a = 2.38110**. Thin-axis diameter 2b = 2.3811.
- Orientation: axis along +z (V1, parallel to f), along +x (V2), and in
  the x-z plane at 45 degrees (V3).

Resolution is set by the THIN axis: 2b/h. QBOX9 levels give h = 0.5 / 0.25
/ 0.125 at L2/L3/L4 (d11: sphere D/h = 6/12/24):

| level | 2b/h | role |
|---|---|---|
| L3 | 9.52 | production rung (matches certified sphere-class resolution) |
| L4 | 19.05 | resolution rung (V4) |

Long-axis clearance at r_e = 2: tip-to-image gap = L_cell - 2a = 9 - 4.76
= 4.24 = 3.56 b. Acceptable for the ratio gate; recorded as the leading
systematic for the absolute gate (§2). r_e = 3 (if run) has 2a = 6.24 —
tip gap 2.76: ratio gate only, do not gate absolutes at r_e = 3 in this
cell.

## 5. Discretization calibration

d11 pinned the sphere bias: r_h/r = 1.034 / 1.028 / 1.013 at D/h = 6/12/24,
collapsing as a_eff ~ a + 0.14 h (rows `d11_rh_collapse`,
`d11_aeff_sign_erratum`). For the spheroid the same surface-smearing
argument applies per axis. Secondary-gate procedure:

- V0 control re-runs the d11 L3 sphere ON THE CURRENT BINARY (pin 6971b13)
  in the identical rundir layout — regression anchor (must reproduce the
  d11 L3 row to 5 digits) and fresh r_h/r at this exact code state.
- Apply the a_eff correction axis-wise (a + 0.14h, b + 0.14h -> recompute
  X^A, Y^A at the effective aspect ratio) for the absolute gate; the
  uncorrected values are recorded alongside.
- V4 (L4 rung) must move the absolute residuals toward zero roughly like
  the sphere ladder did — that, not any single number, is the
  resolution-consistency check.

## 6. Gates

| gate | quantity | band | rationale |
|---|---|---|---|
| G-ratio (PRIMARY) | R_h(perp)/R_h(par) at L3 | 1.14532 +- 2% | lattice + bias cancel in ratio |
| G-abs | R_h per orientation vs a*X^A, a*Y^A (a_eff-corrected) | +- 3% | tracks d11 ladder residuals |
| G-torque | |T| / (a |F|), every run | < 1e-3 (report exact) | torque-path null test |
| G-offdiag (V3) | force direction at 45 deg: F = 6 pi mu a [X^A uu + Y^A (I-uu)] U -> angle of F from the body axis = atan(Y^A/X^A * tan 45) = 48.88 deg, i.e. F drifts 3.88 deg from U toward the perpendicular | drift 3.88 +- 0.5 deg | off-diagonal mobility, first anisotropy observable beyond magnitudes |
| G0 volume | measured alpha-field solid fraction vs (4/3) pi a b^2 / V_cell | same rel. error band as d31 imaged (+-0.5%) | in-situ check of the pe D-3 volume fix + ellipsoid pointInside path |
| V0 regression | sphere K at L3 | reproduce d11 row to 5 digits | anchor |

Steadiness criterion as in D1.1: force plateau flat to 5 digits over the
final ~2 t.u.; steady by t ~ 4-5 expected (d31 precedent).

## 7. Run matrix (all local cluster; all small)

| run | body | orientation | level | est. cost |
|---|---|---|---|---|
| G0/V0 | sphere r=1.5 (control) + ellipsoid create-only smoke | — | L3 | minutes |
| V1 | r_e=2 | axis parallel to f (z) | L3 | ~10 min, ~30 ranks |
| V2 | r_e=2 | axis perpendicular (x) | L3 | ~10 min |
| V3 | r_e=2 | 45 deg in x-z | L3 | ~10 min |
| V4a/V4b | r_e=2 | par + perp | L4 | ~1-2 h each |
| V5a/V5b (optional) | r_e=3 | par + perp (ratio gate only) | L3 | ~10 min |

Total well under one node-day. Rundir naming:
`q2p1_dns_rundir_d61_{v0,v1,v2,v3,v4a,v4b,...}`.

## 8. Implementation prerequisites (the actual work; each pe/FF change
twin-gated e4_l3 before use)

1. **pe interface: ellipsoid creation in the DNS-drag serial setup.** The
   setup path (`setupDNSDragSerial` family) creates spheres from
   `particles.xyz`. Add json-driven ellipsoid support — proposed keys:
   `particleShape_: "ellipsoid"`, `semiAxes_: [a, b, c]`,
   `particleAxis_: [x, y, z]` (unit vector for the a-axis; setup rotates
   the body from +x to this axis), reusing the existing fixed-body flag.
   Sphere path must be byte-identical when `particleShape_` is absent
   (twin proves it).
2. **FF diagnostics: radius semantics.** `getParticleRadius` /
   `DNS_RESOLUTION` (`dofs_per_particle`, `D_over_h`) and the particle-CFL
   scale must use the SMALLEST semi-axis for an ellipsoid; add a
   `get_particle_semiaxes` query or map `getParticleRadius -> min(a,b,c)`
   (decide at implementation; record which). The alpha-field path
   (`pointInsideParticles` ellipsoidType branch, object_queries.cpp:572)
   already exists.
3. **Orientation output** for verification: extend the particle log/VTK
   with the body axis (pe has the quaternion; a `DNS_PART_STATE`-style
   line suffices).
4. **Deck assertion**: refuse `lubricationEnabled_: true` when any
   non-sphere body exists (guard in the setup, loud error).
5. Analysis script `tools/d61_oberbeck_analysis.py`: R_h fixed-point
   inversion (port from the d11 tool), X^A/Y^A evaluation, a_eff
   correction, gate report.

Estimated implementation: 1-2 sessions including twins; no schedule
coupling to v24f.

## 9. Relation to the D/h=16 rung (v24f) — why parallel is safe

v24f gates the RESOLUTION CONVERGENCE of the D5.2 suspension-viscosity
instrument (torque ratio in the annulus). D6.1 gates ANISOTROPIC STOKES
DRAG in the d11 periodic cell. They share no domain, no mesh, no
observable, and no closure; D6.1 carries its own resolution rung (V4) and
its own sphere anchor (V0), so nothing v24f could report would change a
D6.1 design choice. Resources are disjoint by the owner's workload split
(local vs Fritz). The only ordering that matters is INTERNAL to D6.1:
implementation (§8) -> twin -> G0/V0 -> V1..V4.

## 10. What D6.1 deliberately does not test

Free motion (translation/rotation coupling, D-2 buoyancy branch in
anger, the periodic moment-arm), orientation dynamics, and any
concentration effect. Those are D6.2 (Jeffery, local) and D6.3 (suspension,
Fritz, decision pending 6.1/6.2). The torque null (G-torque) is the only
rotational-path observable here, by design.
