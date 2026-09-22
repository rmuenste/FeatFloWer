# D6.4 — Settling spheroid: Stokes terminal velocity and the inertial orientation attractor (CASE_SPEC)

Drafted 2026-09-22 (owner request; the plan's D4.3, deferred when D6.3 took the
third-family slot). Fourth non-spherical family. First case in which an
ellipsoid TRANSLATES freely through the fluid coupling: D6.1 fixed the body,
D6.2 locked translation and freed rotation, D6.3 (designed) locks translation
too. Everything the family needs on the pe side is in the pin 8d326b8
(ellipsoid inertia/volume/buoyancy fixes D-1..D-4, contact, checkpoint
resume); the FF side has never driven a translating ellipsoid, so G0 is a real
gate.

References (all in `literature/`):
- **Happel & Brenner 1973** (`happel_brenner_1973.pdf`), the Stokes drag of
  prolate and oblate spheroids parallel and perpendicular to the symmetry
  axis (correction factors K to Stokes' law on the equal-volume sphere) —
  the Stokes-tier gate. Also Oberbeck 1876 as the origin (optional).
- **Ardekani, Costa, Breugem & Brandt 2016** (`ardekani_2016.pdf`), §3:
  isolated spheroids of aspect ratio 1/3 and 3 settling at Ga = 80–250
  (resolved IBM, 32–48 cells per D_eq, box 15 × 15 × 125 D_eq, periodic
  sides): a spheroid released from rest broadside settles steadily and
  vertically at Ga = 80 (terminal Re 64 prolate / 53 oblate, sphere 83),
  broadside is the only stable orientation at these Ga "independent of the
  initial orientation", and the oblate body's steady vertical regime ends
  in an oscillating path between Ga = 150 and 180. §4.2: spheroid pairs
  (DKT) — the bridge to the contact family, not part of this spec.
- Willmarth, Hawk & Harvey 1964 (not in the folder, optional): regime map
  for thin disks; only if a disk-like body is ever added.
- D6.1 (`d61_v4_resolution`, `d61_v5_halfsize`) for the resolution rules of
  a spheroid in this code; ten Cate E1 (`d13_*`) for the free-settling
  protocol and the wall-approach behaviour; D2.2 for the lubrication
  refusal on non-spheres.

## 1. Physics and gates

Two tiers, three bodies of equal volume (equivalent diameter D_eq = 1):
prolate a/b = 3 (semi-axes 1.0400, 0.3467, 0.3467; an earlier draft wrote 0.7211, 0.4160, 0.4160, which is a/b = 1.73 — caught at G0 staging 2026-09-22), oblate c/a = 1/3
(0.7211, 0.7211, 0.2404), and the sphere (0.5) as the control that ties the
family to the ten Cate ladder.

**Tier S (Stokes).** A spheroid released from rest at Ga ≈ 1 (Re_t ≈ 0.05)
with its symmetry axis held either along or across gravity (fixed
orientation, translation free) reaches the Happel–Brenner terminal velocity
u_t = (ρ_p − ρ_f) g V / (3π μ D_eq K), where K(λ, orientation) is the
tabulated Stokes correction. Confinement is a first-class systematic (the
owner's rule from D6.2): two box widths, and the box-size correction is
measured, not assumed.

**Tier I (inertial).** The same three bodies free in all six DOFs at Ga = 80,
released from rest with the symmetry axis tilted 30° from broadside. Gates:
the body turns broadside-on (the attractor), the terminal Reynolds number
matches Ardekani's, and the path is vertical (no lateral drift, no
oscillation) — the "steady vertical" regime of their fig. 15. A second Ga is
optional (see §4).

| gate | quantity | band | rationale |
|---|---|---|---|
| G0 | free-translation smoke at L3: body created with the semi-axes, indicator DOFs vs analytic volume (±2 %), buoyancy seeding sign (body starts to fall), no NaN over 200 steps, checkpoint written | pass/fail | first translating ellipsoid through the coupling; D-2/D-3/D-4 class checks (guide §13) |
| G-sphere | sphere at Ga = 80: Re_t vs the Ardekani/Yin–Koch value 83 and vs the ten Cate ladder's own drag calibration | ±3 % | ties the family to a certified result |
| G-Stokes-∥, G-Stokes-⊥ (tier S) | u_t of prolate and oblate, axis along and across gravity, box-extrapolated, vs Happel–Brenner | ±3 % after the box correction; the two orientations' ratio (K_⊥/K_∥) ±2 % | the analytic anchor; the ratio cancels most of the resolution bias (guide §13: a_eff does not transfer per axis, use the volume correction) |
| G-box (tier S) | u_t shift between the two box widths | measured, monotone with width | the D6.2 discipline |
| G-attractor (tier I) | tilt angle θ(t) between the symmetry axis and its broadside orientation: decays to < 3° and stays; no reversal | pass/fail + settle time RECORDED | the mechanism gate of the plan's exit criterion |
| G-Ret (tier I) | Re_t = u_t D_eq/ν on the plateau vs 64 (prolate) / 53 (oblate) | ±5 % | their resolution (32–48 cells per D_eq) and domain (125 D_eq tall) differ from ours; 5 % is the honest allowance, the ordering prolate < oblate < sphere is the test |
| G-path (tier I) | lateral displacement over the last 20 D_eq of fall; angular-velocity residual | < 0.2 D_eq; |ω| < 0.05 u_t/D_eq | steady vertical regime |
| G-mesh (exit gate of the plan's D4) | analytic `Ellipsoid` body vs an OBJ-triangulated body of the same spheroid, one run each at Ga = 80 | u_t and θ(t) agree within the resolution band | the mesh-particle classifier the production cases rely on |

## 2. Geometry, units, numbers

Box units as D6.1/D6.2: ρ_f = 1, D_eq = 1, ν chosen per Ga. Density ratio
ρ_p/ρ_f = 1.14 (Ardekani's value at Ga = 80; also close to the ten Cate 1.155).
Gravity |g| along −z; Ga = √((ρ_p/ρ_f − 1) g D_eq³) / ν.

| tier | Ga | ρ_r | g | ν | Re_t expected | u_t (≈ Re_t ν) | fall to plateau | domain (x × y × z) |
|---|---|---|---|---|---|---|---|---|
| S | 1.0 | 1.14 | 1.0 | 0.374 | 0.05 (sphere, unbounded) | ≈ 0.02 | ≈ 3 D_eq | 8 × 8 × 24 and 16 × 16 × 24 (walls) |
| I | 80 | 1.14 | 1.0 | 0.00468 | 83 / 64 / 53 | 0.39 / 0.30 / 0.25 | ≈ 40 D_eq | 15 × 15 × 60, periodic x/y, wall at the bottom, free-slip lid |

Derivations and checks:

- **Tier S ν.** With g = 1 and Δρ/ρ = 0.14, Ga = √0.14/ν = 0.374/ν, so
  ν = 0.374 gives Ga = 1; Stokes u_t of the sphere = Δρ g D²/(18 μ) = 0.0208,
  Re_t = 0.056. Creeping enough for Happel–Brenner (their table is the
  Re → 0 limit; the O(Re) Oseen correction is 3/16 Re ≈ 1 %, below the
  gate).
- **Tier S confinement.** A sphere at the centre of a square duct of side
  8 D_eq settles ≈ 20 % slower than unbounded (Faxén-type wall correction
  ≈ 1 − 1.9 d/W for a cylinder, similar for a duct); at 16 D_eq ≈ 10 %.
  The two-box extrapolation in d/W removes most of it; what remains is
  quoted, as in D6.2. The ORIENTATION RATIO K_⊥/K_∥ is nearly
  confinement-free and is the primary Stokes gate.
- **Tier S orientation lock.** Fixed orientation with free translation is a
  new pe setup mode: `particleMotion_ = "translationOnly"` (angular DOF
  mask (0,0,0) — the mask merged in pin 8d326b8 makes this a one-line
  addition next to "rotationOnly"). At Stokes conditions the broadside and
  edgewise orientations are both neutrally stable only in the exact
  Re → 0 limit; at Re ≈ 0.05 the weak inertial torque would slowly turn
  the body, which is why the lock is needed for a clean K reading.
- **Tier I ν.** Ga = 80 → ν = 0.374/80 = 0.00468; sphere Re_t = 83 gives
  u_t = 0.39, one D_eq in 2.6 t.u. The fall to the plateau (Ardekani: ≈ 30–40
  D_eq) is ≈ 100 t.u.; box height 60 D_eq with release at z = 55 and the
  plateau read between z = 25 and 10 (bottom wall influence starts ≈ 3–5
  D_eq above it — the ten Cate approach curve).
- **Tier I domain.** Ardekani: 15 × 15 D_eq periodic sides; matched. Their
  height 125 D_eq is more than we need for one plateau read.
- **Resolution.** D6.1's rule: resolve the THIN axis, 2b/h ≥ 9.5 (certified
  class). Prolate thin axis 2b = 0.693, oblate 2c = 0.481 → h ≤ 0.0506 for
  the oblate at the certified class; Ardekani used 32 cells per D_eq
  (h = 0.031) and 48 for the flattest bodies. Plan: level 3 at h ≈ 0.10
  (oblate 2c/h ≈ 4.8) for G0 and for the sphere control only; level 4 at
  h ≈ 0.05 (oblate 9.6, prolate 13.9, sphere 20) for every gated run. A
  uniform 15 × 15 × 60 box at h = 0.05 is 300 × 300 × 1200 = 108 M cells —
  too many. Hence a graded or partitioned strategy: (a) tier S boxes are
  small (8 × 8 × 24 at L4 ≈ 6 M cells, fine); (b) tier I uses a moving
  refinement or a reduced box 10 × 10 × 40 at h = 0.05 (32 M cells, ≈ 4×
  the viscometer L4; feasible on 6 Fritz nodes) with a periodic-image
  check (10 vs 15 D_eq at level 3, sphere only) — RECORDED as a
  confinement rung, not a gate. Decide at G0 with the meshing kit in hand.
- **Time step.** Explicit-coupling stability (guide §13): translational
  analogue g_t ≈ u·dt/τ_p with τ_p = ρ_p D²/(18 μ): tier I τ_p = 1.14/(18 ·
  0.00468) = 13.5 t.u., no constraint; tier S τ_p = 0.17, dt ≤ 0.02. Fluid
  CFL: u_t dt/h = 0.39 · dt/0.05 → dt = 0.05 gives 0.4. Rotational
  coupling at ρ_r = 1.14 and Re ≈ 60: g = 15 u (γ̇ dt)/(ρ_r Re) with the
  body's own rotation scale — of order 0.1, safe (the D6.2 instability was
  the creeping-flow corner). Tier I dt = 0.05 (2000 steps to the plateau);
  tier S dt = 0.01 (creeping, τ_p short), 1500 steps.
- **Cost.** Tier S: six runs (3 bodies × 2 orientations, the sphere once)
  × two boxes, each ≈ 15 t.u. on the 8 × 8 × 24 L4 mesh, local cluster,
  hours each. Tier I: three bodies + the OBJ twin + the confinement rung,
  ≈ 100 t.u. each on a 32 M-cell mesh, ≈ 2–3 Fritz segments each on 6
  nodes. Total ≈ 15 node-days Fritz + a few local days. The most expensive
  non-spherical family so far; the tier-S half is cheap and independent.

## 3. Fixture requirements (the build work)

1. **Free-translating ellipsoid through the coupling (G0).** The DNS-drag
   setup's ellipsoid path has `particleMotion_` fixed | rotationOnly; add
   `"free"` (no DOF locks, buoyancy seeding via the ellipsoid branch fixed
   in D-2) and `"translationOnly"` (angular mask (0,0,0), for tier S).
   Verify the FBM force/torque integration for a moving ellipsoid (the
   live-centre moment arm is already reviewed; what is new is the
   indicator update for a translating, rotating ellipsoid every step —
   compare `dofs_per_particle` over time against the analytic volume) and
   the pe checkpoint carries position + orientation + both velocities
   (it does: `pe_checkpoint_roundtrip_test` covers a moving spheroid).
2. **Mesh kits.** Tier S: 8 × 8 × 24 and 16 × 16 × 24 walled boxes
   (`gen_*` in the d52/d62 writer style, uniform cells, L4 target
   h ≈ 0.05); tier I: 10 × 10 × 40 (or 15 × 15 × 60 if the cell count is
   acceptable) with periodic x/y (axis-uniform partition, `SimPar@
   Periodicity`), bottom wall, free-slip lid. Partition sizes per guide
   §12.
3. **Orientation lock mode** in pe setup (`translationOnly`), one-line
   twin-gated addition; the D6.3 mask work makes it trivial.
4. **OBJ twin** (G-mesh): the same spheroid as a triangulated body via the
   `q2p1_creep` ellipsoid.obj path (`setupCreep`); one run at Ga = 80.
5. **Analysis** `tools/d64_settling_analysis.py`: u(t), Re_t plateau
   (window by height band, not time), θ(t) from `DNS_PART_AXIS` (the tilt
   from broadside), lateral drift, ω residual; Happel–Brenner K table
   (prolate/oblate, ∥/⊥) and the two-box extrapolation; Yin–Koch Ga–Re_t
   relation for the sphere control.
6. Lubrication OFF (refused for non-spheres). Contact: the bottom wall is
   reached only if a run is left too long — the plateau read ends 10 D_eq
   above it; contact is exercised deliberately only in the contact family.

## 4. Run matrix

| run | tier | body | orientation / release | box | level | duration | machine | gate |
|---|---|---|---|---|---|---|---|---|
| G0 | — | oblate, free | tilted 30° | 8 × 8 × 24 | L3 | 200 steps | local | G0 |
| S-sph | S | sphere | — | 8 × 8 × 24, 16 × 16 × 24 | L4 | 15 t.u. | local | G-sphere (Stokes), G-box |
| S-pro-∥, S-pro-⊥ | S | prolate, translationOnly | axis ∥ g, ⊥ g | both boxes | L4 | 15 t.u. | local | G-Stokes, G-box |
| S-obl-∥, S-obl-⊥ | S | oblate, translationOnly | axis ∥ g, ⊥ g | both boxes | L4 | 15 t.u. | local | G-Stokes, G-box |
| I-sph | I | sphere, free | — | 10 × 10 × 40 | L4 | 100 t.u. | Fritz, 2 seg | G-sphere (Ga 80) |
| I-pro | I | prolate, free | tilted 30° | 10 × 10 × 40 | L4 | 100 t.u. | Fritz, 2–3 seg | G-attractor, G-Ret, G-path |
| I-obl | I | oblate, free | tilted 30° | 10 × 10 × 40 | L4 | 100 t.u. | Fritz, 2–3 seg | G-attractor, G-Ret, G-path |
| I-obl-obj | I | oblate as OBJ mesh | as I-obl | 10 × 10 × 40 | L4 | 100 t.u. | Fritz, 2–3 seg | G-mesh |
| I-conf (L3) | I | sphere | — | 10 vs 15 D_eq wide | L3 | 100 t.u. | local | confinement, RECORDED |
| (optional) I-obl-160 | I | oblate, free | tilted 30° | 10 × 10 × 40 | L4 | 150 t.u. | Fritz, 4 seg | oscillating-path regime (Ardekani Ga 150–180), RECORDED |

Order: G0 → tier S in one batch (cheap, local, the analytic anchor) → I-sph
and I-obl in parallel → I-pro → I-obl-obj → I-conf.

## 5. What this family deliberately defers

Spheroid pairs (DKT, Ardekani §4.2 — the contact family's first case),
disks and other sharp-edged bodies (Willmarth), the oscillating and
chaotic regimes above Ga ≈ 150 beyond one recorded run, a full Ga ladder,
and any comparison against the experimental settling literature for
non-spheres.

## 6. Relation to the rest of the campaign

Uses the D6.1 body definitions and resolution rules, the D6.2/D6.3 setup
keys and the angular mask, the ten Cate free-settling protocol and its
wall-approach knowledge, and the checkpoint-resume chainer. Closes the
plan's D4 exit gate lines "Oberbeck terminal velocity green" (tier S) and
"orientation-attractor mechanism green" (tier I), and the "run both once"
requirement (G-mesh). Independent of D6.3; shares its analysis-tool work
(tilt trace). Website: held off until the family closes (content policy).
