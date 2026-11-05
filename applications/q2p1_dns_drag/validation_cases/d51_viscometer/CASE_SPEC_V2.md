# D5.1 V2 — through-hole numerical viscometer (case specification)

Status: DRAFT approved-in-principle by owner 2026-08-26 ("Agreed, draft the
CASE_SPEC and start the meshing/partitioning"). Supersedes the blind-bore
cell of CASE_SPEC.md for the Einstein gate and the φ ladder; V0/V1 results
stand as recorded (rows `d51_v0_baseline` … `d51_end_dilution`).

## 1. Why V2 (recorded rationale)

The blind-bore cell generates 45.7% of its baseline torque in
particle-free fluid (velocity-jump ring at the top wall + bob end face),
which dilutes η(φ) = T(φ)/T(0) by an unmeasurable factor
(row `d51_end_dilution`). The legacy FF viscometer (`Couette_20x4x10`,
recovered 2026-08-26) avoided this entirely: the coaxial hole spans the
full height and both end faces are symmetry planes, so the flow is an
axially uniform annular Couette — quasi-2D, no rings, no end faces,
torque per height directly comparable to the analytic solution and to
the 2D Prignitz–Bänsch curves. V2 mirrors that instrument in campaign
units.

## 2. Legacy correspondence

| item | legacy Couette_20x4x10 | V2 (d = 1 units) |
|---|---|---|
| inner radius r_i | 0.2 | 5 |
| outer radius r_a | 0.4 (ratio 0.5) | 10 (ratio 0.5) |
| gap / d | 0.2/0.04 = 5 | 5 |
| height | 0.8 (H/d = 20) | 10 (axial uniformity ⇒ height is statistics only) |
| z boundaries | Symmetry001 both ends | Symmetry001 both ends |
| inner-wall BC | Inflow770 (legacy rotation driver) | Inflow14 (u = Ω×r, SimPar@RPM) |
| torque mask | WallF on the OUTER wall | WallF on the INNER wall (see §5; outer-mask twin optional) |
| particle d | 0.04, φ up to 0.20 | 1.0, φ = 0.05 first (ladder to 0.20 later) |

## 3. Geometry and boundary conditions

Annulus r ∈ [5, 10], z ∈ [0, 10]. The bob is the un-meshed hole through
the FULL height (no end face anywhere).

| component | .par type | geometry line | motion |
|---|---|---|---|
| inner wall (bob) | `Inflow14` | type 7: `0 0 0 5.0 1 1 0` | u = Ω×r |
| inner wall overlay | `WallF` | same vertex set | torque mask |
| outer wall | `Wall` | type 7: `0 0 0 10.0 1 1 0` | 0 |
| z = 0 | `Symmetry001` | type 4 plane `0 0 1 0` | free slip, w = 0 |
| z = 10 | `Symmetry001` | type 4 plane `0 0 1 -10` | free slip, w = 0 |

No plane/cylinder .prj-ordering trap this time (no shared ring vertices
between bob surfaces — there is only one bob surface).
BndrForce is a single global mask ⇒ WallF tags ONLY the inner wall.

## 4. Operating point

| quantity | value | note |
|---|---|---|
| ν, ρ_f | 0.2, 1.0 | as V0/V1 |
| Ω | 0.1 (SimPar@RPM = 0.95493) | as V0/V1 |
| shear rate γ̇ ≈ Ω r_i / gap | 0.1 | |
| Re_p = γ̇ d²/ν | 0.5 | Stokes-ish |
| Taylor Ω² r_i gap³/ν² | 156 | subcritical (≪ ~1712) |
| τ_ν = gap²/ν | 125 | run budget anchor |
| dt (V2.0 / V2.1) | 0.05 / 0.005 | contacts set the V2.1 step (row d51_v1_run) |
| V2.0 budget | 4000 steps → t = 200 = 1.6 τ_ν | single job |
| V2.1 budget | segments as needed; η(t) flatness is the stop signal | restart chain certified (row d51_v1_run) |

## 5. Measurement and gates

Physical-torque protocol throughout (row `d51_torque_metrology`):
DNA (deformation form) primary, RES + transpose correction as
concordance check. Transpose correction here:
2μΩ·V_hole = 2·0.2·0.1·π·5²·10 = **31.42**.

**V2.0 (empty baseline) gates:**
1. ANALYTIC TORQUE GATE (new — the through-hole geometry earns it):
   T_exact = 4πμΩ/(r_i⁻² − r_a⁻²)·H = 4π·0.02/(0.04−0.01)·10 = **83.78**.
   Gate: |T_DNA(∞)/T_exact − 1| ≤ 3% (discretization only; there are no
   end effects to excuse).
2. u_θ(r) profile vs exact Couette solution (VTK spot check, recorded).
3. Estimator concordance ≤ 3% after the 31.42 correction.
4. Steadiness: endpoint drift < 5e-3 (τ_ν = 125; plateau expected ~t≥150).
5. OPTIONAL V2.0b: rerun with the WallF mask swapped to the outer wall
   (one .par edit, same flow) — inner/outer torque balance closes the
   loop the legacy study measured from; no code change needed.

**V2.1 (φ = 0.05):** N = 225 spheres (φ·V_ann/V_p = 0.05·π·75·10/(π/6)),
RSA, min center distance 1.05d, SURFACE clearances 0.5d to both walls and
0.5d to both symmetry planes (FBM does not mirror particles). ρ_p/ρ_f =
1.1 with gravity 0 (rheologically inert, certified coupling regime —
rows d51_v1_run / ff-deck-staging-pitfalls). Deck from the V1 rundir
deck (NOT the template), ForceScale all-ones, stripped comments.
Gate: η = T(φ)/T(0) vs the composite-Einstein target from the measured
φ(r,z) profile — UNDILUTED in this geometry (the only particle-free
zones are the modeled radial/axial margins). Naive 1 + 2.5φ = 1.124
recorded as reference. Then the φ ladder (0.10, 0.20 — legacy range) and
the D/h = 16 resolution rung become meaningful.

## 6. Mesh order (meshing workflow)

Pure O-grid annulus, NO core, NO butterfly transition:
- coarse target ≈ d/2: radial 10 layers UNIFORM across the gap (particles
  live at all radii — do not wall-grade like the legacy mesh), azimuthal
  94 (mid-radius spacing 2π·7.5/94 ≈ 0.50), axial 20 → 18,800 coarse hexes.
- production level: MaxMeshLevel 3 → mid-gap h ≈ 0.125 = d/8 (campaign
  operating point).
- .par set per §3 (type-7 cylinders flags `1 1 0`, type-4 planes).
- partition 108 subdomains (METIS acceptable — no periodicity).
- verification: refine once, r-histograms of boundary nodes at r = 5 and
  r = 10; z of symmetry-plane nodes exact.
- deliverables mirror d51_mesh_v1: mesh dir (`d52_mesh_v1/`) + README +
  REPRODUCE + partitioned `_mesh/VISCO2_108`.

## 7. Run sequence

1. V2.0 empty baseline (owner wants to inspect first output before
   continuing).
2. V2.1 φ = 0.05 (PE_SERIAL_MODE, as the whole campaign: q2p1_viscometer
   from build-dns-pe-lubaddon at pe pin f9b7115).
3. Decision point: φ ladder / D/h = 16 rung / legacy-study line-up
   (pending the old study's MaxMeshLevel and comparison-curve family).

## 8. Sign-off state

1. Through-hole + symmetry ends + legacy mirror (gap 5d, ratio 0.5,
   H = 10d): owner agreed 2026-08-26.
2. V2.0-first sequencing: owner directive 2026-08-26.
3. PE_SERIAL_MODE for particles: confirmed 2026-08-26.
4. Operating point ν/Ω unchanged from V0/V1: proposed, not yet vetoed.
5. WallF on inner wall (outer-mask twin optional): proposed in this
   draft — flag at V2.0 staging if the legacy outer-wall convention is
   preferred as primary.
