# D6.3 — Oblate spheroid orbit ladder at finite Re (CASE_SPEC)

Drafted 2026-09-19 (owner decision the same day: the third non-spherical
family is the finite-Re orbit ladder on an OBLATE body, not the settling
spheroid nor the suspension of rotating bodies). Builds on D6.2 (prolate
Jeffery orbit, CLOSED 2026-09-15): same Couette box, same setup path, same
analysis discipline; new physics = fluid inertia, new body = a = b = 2c.

References (both in `literature/`, index `literature/README.md`):
- **DA00** Ding & Aidun, JFM 423 (2000) 317–344 — the primary gate: a
  neutrally buoyant oblate spheroid a = b = 2c tumbling about the vorticity
  axis; period grows with Re and diverges at Re_c = 81 as
  GT = C (Re_c − Re)^(−1/2), C = 200 (fit Re 50–81); steady orientation
  above Re_c. Read 2026-09-19 (print copy); digest with figure read-offs:
  `LITERATURE_DIGEST.md` next to this spec (page/figure-tagged; the
  read-off uncertainties quoted there are the reading uncertainty only).
- **DG25** Di Giusto, Bergougnoux & Guazzelli, JFM 1017 (2025) A41 —
  experiments on oblate bodies at Re_p = 0.03–5; one oblate spheroid
  (ELL06, r = 0.561) rotating at every Re_p tested with T_J/T ≈ 1.0 up to
  Re_p 0.73 and ≈ 0.75–0.8 at Re_p 1.03 (fig 6); the orbit-selection
  theory (Einarsson 2015b, Dabade 2016, Rosén 2015a) it relays: for oblate
  aspect ratio > 0.137 the tumbling orbit is UNSTABLE at small inertia and
  spinning (log-rolling) is the attractor, up to Re_p ≈ 5.
- Jeffery 1922 (Stokes anchor), `tools/d62_jeffery_analysis.py`.

Reynolds-number conventions (three in play — pin them):

| symbol | definition | this case (2b = 1, γ̇ = 0.2) |
|---|---|---|
| Re (DA00) | ρ_f G (2b)² / μ, 2b = major diameter | = γ̇/ν = 0.2/ν |
| Re_p (DG25) = Re_a (guide) | γ̇ b² / ν, b = major semi-axis | = Re/4 |
| Re_c(thin) (stability bookkeeping) | γ̇ c² / ν | = Re/16 |

All rung labels below use **Re (DA00)**. DA00's Re_c = 81 is Re_p ≈ 20.

---

## 1. Physics and gates

Body: oblate spheroid, semi-axes 0.5, 0.5, 0.25 (aspect ratio
λ = c/b = 0.5, thin axis 2c = 0.5 — the same thin dimension as the D6.2
prolate, so the certified 2c/h = 10.4 kit carries over). Translation locked,
centre at the box centre, planar Couette box as D6.2 (walls z = ±H/2 at
±U x̂, periodic x, y; vorticity +γ̇ ŷ). Symmetry axis n starts in the
shear plane.

Two experiments, because the literature answers two different questions:

**(A) Constrained tumbling (DA00 comparison).** DA00's model is a
one-degree-of-freedom rotation about the vorticity axis ("The axis of
rotation is always the x-axis", their §6); the small-inertia theory says a
free oblate body of this aspect ratio would leave the tumbling orbit for
log-rolling. So the DA00 rungs run with the angular DOFs locked to the
vorticity axis (fixture item 1). Observables: rotation period from the
spacing of the φ̇ peaks (two per period), the −1/2 divergence law across
rungs, the arrest above Re_c.

**(B) Free rotation (orbit selection).** The same body with all three
angular DOFs free, started in the shear plane with a deliberate small tilt.
At Re_p ≈ 1 (Re = 4) DG25's spheroid still tumbles with a ≈ 20–30 %
longer period and no completed drift within 8 periods; whether the tilt
grows (log-rolling attraction) is RECORDED, not gated — neither paper
measures the drift rate for a spheroid, and above Re_p ≈ 5 the selection
is not in the literature we hold.

Jeffery reference for the oblate body (Stokes anchor and the Re → 0 limit of
everything): with φ = atan2(n_z, n_x) the symmetry-axis angle from the flow
axis, |dφ/dt| = γ̇ (cos²φ + λ² sin²φ)/(λ² + 1): fast 0.8 γ̇ with n along the
flow, slow 0.2 γ̇ with n along the gradient (the disk face in the flow
plane); period T γ̇ = 2π(λ + 1/λ) = 15.70796 — numerically identical to the
prolate r_e = 2 period, with the slow and fast phases swapped in φ.

| gate | quantity | band | rationale |
|---|---|---|---|
| G0 | body creation with semiAxes_ = [0.25, 0.5, 0.5] (a = symmetry axis on particleAxis_), indicator DOFs vs analytic volume, DNS_PART_AXIS = symmetry axis, angular mask active (ω_x = ω_z = 0 to round-off) | DOF count ±2 % | the D-4-class check; the a < b ordering has never been run |
| G-Stokes (R0) | period vs 15.70796 and waveform vs Jeffery with λ = 0.5, rotation sense | ±1 % period, 3 % rms waveform | the oblate anchor; D6.2 reached +0.30 % |
| G-DA-period (R50, R60, R70) | GT vs DA00 fig 22: 36.7, 44.7, 59.5 (±1 read-off) | ±8 % | their body is 8–16 lattice nodes across and their confinement is inferred, ours differs (§2); an 8 % band is the honest resolution/confinement allowance, the trend is the test |
| G-DA-law | fit GT = C (Re_c − Re)^(−1/2) to our three rungs with the exponent fixed; then free exponent | Re_c within 10 % of 81; free exponent −0.5 ± 0.1 | the universal part of DA00 (saddle-node), independent of their numerics |
| G-arrest (R90) | φ̇ → 0 and stays (no π-crossing after the settle), settle time | rotation stops within Gt ≈ 20; no crossing to Gt = 40 | DA00 fig 21(d): arrest by Gt ≈ 10 with a small undershoot |
| G-free-period (F4) | T_J/T at Re = 4 (Re_p = 1.03) vs DG25 ELL06 0.75–0.8 (±0.1) | 0.65–0.95 | the only spheroid datum in DG25 |
| G-free-drift (F4, F50) | tilt of n out of the shear plane vs time; direction and rate | RECORDED | theory: toward log-rolling below Re_p ≈ 5; unknown above |
| G-wall (C50) | period at H = 4 vs H = 8, Re = 50 | measured, quoted | the D6.2 discipline: wall effect is a first-class systematic |

## 2. Geometry, units, numbers

Box units ρ = 1; ν set per rung (μ = ν). Shear γ̇ = 0.2 as D6.2 (U = 0.8 at
H = 8, 0.4 at H = 4) so that every deck differs from the D6.2 deck in
viscosity, density and the new keys only. Then Re = 0.2/ν:

| rung | Re | ν | ρ_r | angular DOFs | dt | γ̇dt | CFL U dt/h | g (thin axis) | period T (t.u.) | run length (t.u.) | steps | L4 segments |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R0 | 0.2 | 1 | 10 | y only | 0.01 | 0.002 | 0.17 | as D6.2 (0.36 sphere-equivalent) | 78.5 | 105 | 10 500 | 3 (as D6.2 V1) |
| R50 | 50 | 0.004 | 1 | y only | 0.05 | 0.01 | 0.83 | 0.10 | 184 (DA00) | 420 | 8 400 | 2 |
| R60 | 60 | 0.00333 | 1 | y only | 0.05 | 0.01 | 0.83 | 0.08 | 224 | 500 | 10 000 | 2–3 |
| R70 | 70 | 0.00286 | 1 | y only | 0.05 | 0.01 | 0.83 | 0.07 | 298 | 650 | 13 000 | 3 |
| R90 | 90 | 0.00222 | 1 | y only | 0.05 | 0.01 | 0.83 | 0.05 | arrest | 200 (Gt = 40) | 4 000 | 1 |
| C50 | 50 | 0.004 | 1 | y only | 0.05 | 0.01 | 0.83 | 0.10 | 184 | 420 | 8 400 | 1 (half the cells) |
| F4 | 4 | 0.05 | 1 | free | 0.02 | 0.004 | 0.33 | 0.48 | ≈ 98 (DG25) | 3 periods ≈ 300 | 15 000 | 4 |
| F50 | 50 | 0.004 | 1 | free, 1° tilt | 0.05 | 0.01 | 0.83 | 0.10 | 184 | 2 periods ≈ 400 | 8 000 | 2 |

Derivations and checks behind the table:

- **Density ratio.** DA00's Re_c depends on α (their fig 16: Re_c rises
  with α), so every finite-Re rung is neutrally buoyant, ρ_r = 1 — the
  value that was UNSTABLE in D6.2 at Stokes conditions. The explicit
  rotational coupling gain (guide §13) is g = 15 u (γ̇dt)/(ρ_r Re_thin)
  with u ≈ 2; it scales as 1/Re, so at Re ≥ 50 it is ≤ 0.1 even with the
  step raised 5× (column g). At Re = 4 it is 0.48 at dt = 0.02 — inside
  the g ≲ 0.5 rule; at dt = 0.05 it would be 1.2 (unsafe). Only the Stokes
  anchor keeps the D6.2 device ρ_r = 10 (τ_rot γ̇ = 0.008, zero-inertia
  Jeffery intact).
- **Time step.** Fluid CFL = (γ̇dt)(H/2h) depends on γ̇dt only: 0.83 at
  γ̇dt = 0.01 for the H = 8 L4 kit (h = 0.0481) — acceptable for the
  implicit Q2/P1 solver, and 1500+ steps per period even at the shortest
  period. Cell Reynolds number u h/ν = Re · (H/2) h ≈ 0.19 Re → 17 at
  Re = 90; the ten Cate E4 rung ran comparable values.
- **Duration.** DA00 start at χ = 0 (their b-axis along the gradient, i.e.
  the symmetry axis along the FLOW: φ = 0, the fast phase); fig 21 shows
  the periodic state established within the first pulse. Period from peak
  spacing needs 5 peaks = 2 periods; run length = 2 T + 50 t.u. transient
  (DA00 periods from fig 22: GT = 36.7 / 44.7 / 59.5 → T = GT/γ̇).
- **Cost.** D6.2 L4: 17–21 s/step on 141 ranks (2 Fritz nodes) → 41 t.u.
  per 24 h segment at dt = 0.01; at dt = 0.05 ≈ 4 300 steps per segment.
  Required set (R0, R50, R60, R70, R90, C50) ≈ 12 segments ≈ 24 node-days;
  tier 2 (F4, F50) ≈ 6 segments. Between D6.2 and the v25f rung.
- **Confinement.** DA00's 3-D domain is 40 × 200 × 80 lattice nodes with
  the body at 8 × 8 × 4 or 16 × 16 × 8 nodes; the axis assignment is not
  stated (digest §A.2). If 80 is the gradient direction, H/b = 10 (fine
  body) or 20 (coarse body). Our H = 8 kit gives H/b = 16, between the two;
  the H = 4 kit gives 8. DA00 quantify confinement only for the circular
  cylinder, and assert the −1/2 exponent is confinement-independent. So:
  the ladder runs on the existing H = 8 kit, C50 measures the wall
  sensitivity at Re = 50 with the H = 4 kit, and a matched H = 5 kit
  (H/b = 10, 13 coarse cells) is built only if C50 shows a shift larger
  than the 8 % gate band. Vorticity extent: ours L_y = 6 = 12 b, theirs
  40 nodes = 5–10 a.
- **Resolution.** Thin axis 2c = 0.5 → 2c/h = 10.4 at L4 on the H = 8 kit
  (certified class ≥ 9.5, guide §13); major axes 2b/h = 20.8. DA00 resolve
  the thin axis with 4–8 nodes. Level-3 smoke only (2c/h ≈ 5.2).
- **Initial orientation.** R-rungs: n = +x̂ (φ = 0, fast phase, as DA00).
  F4: n = (cos 5°, sin 5°, 0)-type tilt of 5° out of the shear plane
  (toward the vorticity axis) so the drift direction is readable; F50:
  1° tilt. Recorded as `particleAxis_`.

## 3. Fixture requirements (the build work)

1. **pe angular DOF mask** (`setAngularDofMask(Vec3)`, mirror of
   `setLinearDofMask`; zero components kill ω_x/ω_z after every angular
   update, the torque path untouched). New json key
   `angularDofMask_ = [0,1,0]` read by `setupDNSDragSerial` for
   `particleMotion_ = "rotationOnly"`; default [1,1,1] keeps every existing
   deck byte-identical (e4_l3 twin gates the pin, as always). The resume
   path must re-apply the mask after a checkpoint load exactly as it
   re-applies the linear lock. Unit test: a sphere spun with ω = (1,1,1)
   under mask [0,1,0] keeps only ω_y; the D6.2 sphere control (V0) with
   the mask must still spin at γ̇/2 about y.
2. **a < b ordering.** `semiAxes_ = [0.25, 0.5, 0.5]` puts the symmetry
   axis on `particleAxis_` and makes `DNS_PART_AXIS` log n directly. Audit
   everything that assumes a is the major axis: AABB padding
   (`maxSphereRadius`, `aabbPadding`), `wallCenterOffset`, the ellipsoid
   containment/inertia unit test (extend it with an oblate case), the FBM
   indicator. G0's DOF count is the acceptance check. Already covered by
   the pe branch `feature/ellipsoid-contact` (2026-09-19, not yet merged):
   the body's own bounding box (exact for any axis ordering) and the
   support function; the remaining items above are still to be audited.
3. **Deck keys.** `Prop@Viscosity` = ν per rung with `fluidViscosity_` in
   the json matching (memory `ff-viscosity-convention`: kinematic; ρ = 1 so
   μ = ν); `particleDensity_` = 1.0 (finite Re) / 10.0 (R0);
   `SimPar@TimeStep` = json `stepsize_` = dt per rung (memory
   `dns-dt-stability-floor`: they must be equal); `SimPar@GammaDot = 0.2`
   for the linear initial field; periodic lengths and the axis-uniform
   5 × 4 × 7 partition as D6.2; `PeCheckpointOnDump = Yes` and the
   checkpoint-resume chainer (guide §12) — orientation AND angular velocity
   carry across segments (the R-rungs are 2–3 segments each).
4. **Analysis.** `tools/d62_jeffery_analysis.py`: accept λ < 1 (check the
   waveform formula and any r_e ≥ 1 assumption), add the peak-spacing
   period estimator (two φ̇ peaks per period; robust near Re_c where the
   axis lingers), an arrest detector (no π-crossing after t_settle, |φ̇|
   below 1e-3 γ̇), and the tilt angle asin|n_y| vs time for the F-rungs.
   New `tools/d63_orbit_ladder.py`: GT(Re) table across rungs, the fixed-
   and free-exponent fits, the DA00 overlay (fig 22 read-offs from the
   digest), the Jeffery line.
5. **Reference discipline.** Fig 22 read-offs carry ±1 GT; do NOT gate on
   the fig 21 absolute rate extrema (their ordinate is inconsistent with
   Jeffery by ≈ 0.1 — digest §A.2); pulse widths and periods only.
6. Lubrication OFF (refused for non-spheres anyway); PE serial mode; runs on
   Fritz (2 nodes / 141 ranks per rung, several rungs in parallel).

## 4. Run matrix

| run | Re | body/DOFs | box | level | duration | est. cost | gate |
|---|---|---|---|---|---|---|---|
| G0 | 50 | oblate, y-locked | 8×6×8 | L3 | 50 steps | minutes | creation, DOF count, mask, axis record |
| R0 | 0.2 | oblate, y-locked, ρ_r = 10 | 8×6×8 | L4 | t = 105 | 3 seg | G-Stokes |
| R50 | 50 | y-locked | 8×6×8 | L4 | t = 420 | 2 seg | G-DA-period |
| R60 | 60 | y-locked | 8×6×8 | L4 | t = 500 | 2–3 seg | G-DA-period, G-DA-law |
| R70 | 70 | y-locked | 8×6×8 | L4 | t = 650 | 3 seg | G-DA-period, G-DA-law |
| R90 | 90 | y-locked | 8×6×8 | L4 | t = 200 | 1 seg | G-arrest |
| C50 | 50 | y-locked | 8×6×4 | L4 | t = 420 | 1 seg | G-wall |
| F4 (tier 2) | 4 | free, 5° tilt | 8×6×8 | L4 | 3 periods | 4 seg | G-free-period, G-free-drift |
| F50 (tier 2) | 50 | free, 1° tilt | 8×6×8 | L4 | 2 periods | 2 seg | G-free-drift |
| (optional) R75 | 75 | y-locked | 8×6×8 | L4 | t = 900 | 4 seg | strengthens G-DA-law near Re_c (DA00: GT = 82.3) |
| (optional) H = 5 kit | 50 | y-locked | 8×6×5 | L4 | t = 420 | 1 seg + meshing | only if C50 shifts > 8 % |

Order: G0 → R0 and R50 in parallel (R0 certifies the oblate anchor while
R50 already tests inertia) → R70, R90 → R60, C50 → tier 2.

## 5. What this family deliberately defers

The orbit-selection map above Re_p ≈ 5 (Rosén et al. 2015a/b are not in
the folder), any density-ratio ladder (DA00 fig 16 is 2-D), the steady
orientation angle above Re_c (DA00 give it only for the 2-D ellipse), the
suspension of rotating bodies, and flat disks/rings (DG25's alignment
results need sharp-edged bodies; see the plan's Willmarth 1964 note for a
falling-body family). Non-spherical CONTACT is a separate family opened on
the pe `feature/ellipsoid-contact` branch (2026-09-19).

## 6. Relation to the rest of the campaign

Reuses the D6.2 kits (d62_mesh_v1/v2), partition, moving-wall BC, setup
keys, checkpoint-resume chainer and analysis tool. Independent of the D5.2
dilute rung (v26/v26f) and the closure work; can start as soon as fixture
item 1 is merged and twin-gated. Closes the D4 exit gate's "finite-Re orbit
ladder" line and gives the website its first inertial non-spherical result
(held off the site until the ladder closes, per the content policy).
