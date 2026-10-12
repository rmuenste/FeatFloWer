# FBM classification query — Layer 2 conclusion

Date: 2026-10-12. Companion to `fbm_query_layer1_conclusion.md`. Harness:
`applications/q2p1_hashgrid_test/` (its README has the build line and the construction).
Datasheet rows: `fbm_query_layer2_harness`, `fbm_query_layer2_gate`.

## What Layer 2 asked

Layer 1 compared pe's accelerated point query against the exhaustive one on synthetic
points. Layer 2 asks the same question where it matters: every Q2 DOF of a real partitioned
FeatFloWer mesh, through the real `QuadScalar_FictKnpr` sweep, on every rank, for the full
8-byte system id as well as the solid flag — with no fluid solve, and with the verdict able
to fail the run.

Mesh `QBOX9Z3`, unit cube, 9³ hexes at level 2, three z-slab workers plus the master:
37³ = 50 653 Q2 DOFs, 53 391 rank-local classifications per sweep. Serial PE, so every
worker holds all 70 bodies and pe's per-rank 32-body threshold is met globally.

## Result

| stage | mode | compared | flag viol. | owner viol. | verdict |
|---|---|---|---|---|---|
| A1 lattice (refreshed) | must pass | 53 391 | 0 | 0 | ok |
| A2 teeth (nudged, stale) | must differ | 53 391 | **336** | 0 | ok |
| G sphere / box / capsule / ellipsoid / cylinder (stale) | must differ | 53 391 | 33 / 245 / 17 / 5 / 36 | 0 | ok |
| R sphere / box / capsule / ellipsoid / cylinder (refreshed) | must pass | 53 391 | 0 | 0 | ok |
| E overlapping pairs across a plane | must pass | 53 391 | 0 | **0** | ok |
| P probes on the z = 1/3, 2/3 subdomain cuts | must pass | 53 391 | 0 | 0 | ok |

Exit 0. Every miss is `accelerated=F exhaustive=T` at the predicted DOF planes
`x = 13/36, 19/36, 25/36`, just past `2·span, 3·span, 4·span`. On the real mesh the stale
spatial index misses **hundreds of solid DOFs per step** for a body that has crossed a cell
boundary; the one `refresh_spatial_index()` after `step_simulation()` (FF `dc82d7fb`,
pe #43) removes every one of them for all five shapes. Overlap ownership agrees on the full
id under the minimum-id rule now in both paths.

Controls: `UseHashGridAccel = No` → preflight fails, exit 64. `GateRefreshIndex = No` → the
first verdict stage fails with 2507 misses, exit 62. A run that compared fewer than 1000
classifications is vacuous and fails regardless of its violation count.

## What the harness forced us to learn

1. **The construction must be polydisperse.** pe files a body in the grid whose span is
   `√hierarchyFactor = √2` times its AABB (`HashGrids.h:175, :1602`). Equal bodies fill
   70.7 % of a cell, straddle at most two adjacent cells, and the 27-cell stencil always
   covers the filed cell — a miss is geometrically impossible. The original plan's uniform
   population would have passed every stage and been called a green gate. Spacers set the
   span and trip the threshold; probes fill the span. Creation order is load-bearing.
2. **The real-mesh nudge is one DOF pitch, not a few δ.** Layer 1 chose query points just
   inside the far surface. Here the query points are the DOF planes; a miss needs a DOF
   plane *beyond* the next cell boundary inside the body. The first run nudged `2e-5`,
   compared 53 391 DOFs, found nothing — and the `A2` control refused the run. That refusal
   is the single most valuable behaviour of the gate.
3. **The lattice plane must not alias the shape.** Probe `j` had shape `j mod 5` and plane
   `j mod 5`; every shape sat on one plane, and the cylinder's plane had its far boundary
   outside the domain. `G_cylinder` showed 0 misses while every other shape differed. The
   second time a cylinder has looked immune in this project, the second time it was the
   harness. Fixed in pe `8f823c4`.
4. **Probes must be sized by x half-extent, never bounding radius**, and the real post-creation
   AABB extent must be asserted — a cylinder inscribed in a sphere has only `√3/2` of the
   diameter along its axis and would never fill the span.
5. **Nothing in pe's C surface can move a body**; `setParticle2` sets force and torque only,
   and `setParticle`/`getParticle` in `dem_query.f90` bind to symbols that do not exist. The
   gate's prescribed motion therefore lives in pe (`fbm_gate_probe.cpp`), where the span,
   the contact threshold and each body's AABB are knowable.
6. **The collective belongs in the harness loop, not in the sweep.** `COMM_SUMMN` is either a
   master-mandatory ring or an Allreduce on `MPI_COMM_WORLD`; `QuadScalar_FictKnpr` runs
   rank-asymmetrically and has already produced a deadlock and an `MPI_ERR_TRUNCATE`.
   `MPI_Abort` is not collective, so the verdict can fail the run from any rank.

Two incidental findings: FullC0ntact's `init_fc_rigid_body` reads `start/sampleRigidBody.xml`
unconditionally and the `heat` app's copy points it at a CGAL config file (the gate app now
carries `q2p1_dns_drag`'s); and the deck-echo at `param_parser.f90:1406` segfaulted once
(non-SAVEd `nFields` indexing an allocatable) and never reproduced — recorded, not fixed.

## What Layer 2 leaves to Layer 3 and after

- Stages C (rotated shapes, with a fixed box as a negative control), B (a second spacer
  size → two grids), D (asset-free triangle meshes), F (randomised, demoted) — the spine is
  green, these extend it.
- The **multiply-owned DOF count** is not measured: stage E checks that both paths agree
  on an owner, not how many DOFs have two. Sizing the CSR exception table needs a
  multi-owner query.
- **Parallel PE is untested** by this harness (serial PE throughout); its HashGrid default
  stays off.
- DistanceMap stays deferred (needs `USE_CGAL` in the FF tree; Layer 1 D2 characterised it as
  an approximation a classification gate cannot assert on).
- The `libs/pe` pin bump to a commit containing #43 and this branch is the owner's
  twin-gated call: e4_l3 on both machines, plus DKT since #43 factored
  `refreshSpatialIndex()` out of `findContacts()`.
- Layer 3: the `UseHashGridAccel` on/off A/B on a viscometer dump, deliberately last.
