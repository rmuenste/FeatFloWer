# q2p1_hashgrid_test — the Layer-2 FBM classification gate

Does pe's accelerated point query (`pointInsideParticlesAccelerated`, HashGrid) agree with the
exhaustive one (`pointInsideParticles`) on **every Q2 DOF of a real partitioned FeatFloWer
mesh**, on the solid flag *and* the full 8-byte system id, on **every rank**? This application
asks exactly that and nothing else. There is no fluid solve. It is Layer 2 of the three-layer
verification of the FBM acceleration; Layer 1 is pe's own `tests/interface/pe_fbm_query_stress*`.

## Build

A separate tree. The gate is a build option and is forced off without acceleration:

```
cmake -DUSE_PE=ON ..                                       # step 1
cmake -DUSE_PE_SERIAL_MODE=ON ..                           # step 2
cmake -DENABLE_FBM_ACCELERATION=ON -DENABLE_FBM_QUERY_GATE=ON ..   # step 3
cmake --build . --target q2p1_hashgrid_test_stage
cd applications/q2p1_hashgrid_test && ./run_gate.sh        # mpirun -np 4
```

`ENABLE_FBM_QUERY_GATE` defines `FBM_QUERY_GATE`, which is what compiles FeatFloWer's
comparison code in `fbm_getFictKnprFC2` for the first time (it used to sit behind
`VERIFY_HASHGRID`, a symbol pe only ever defined for its own targets). A binary built without the
option aborts at preflight rather than reporting a green run it could not have earned.

## What runs

Mesh `QBOX9Z3`: the unit cube, 9³ hexes refined once (`MaxMeshLevel = 2`, 18³ cells), three z-slab
subdomains, so `mpirun -np 4`. 37³ = 50 653 Q2 DOFs, ~17 800 per worker. Serial PE, so every
worker holds the whole particle world.

The program initialises exactly as a solver would, takes **one** particle step (that is what
activates pe's grid and arms the accelerated branch), runs the preflight assertions, then the
stages from `SimPar@GateStageList`:

| stage | construction | kind | expected misses |
|---|---|---|---|
| `A1` | all probes on the lattice, index refreshed | verdict | 0 |
| `A2` | A1 nudged one DOF pitch (`1/36`) past the planes, **no** refresh | **control** | **> 0** |
| `G` | per shape: lattice, refresh, nudge, no refresh | **control** | **> 0** |
| `R` | per shape: lattice, refresh, nudge, **refresh** | verdict | 0 |
| `E` | probe pairs overlapping a third of a diameter across a plane; full-id ownership | verdict | 0 |
| `P` | probes centred on the `z = 1/3`, `z = 2/3` subdomain cuts | verdict | 0 |

A **miss** is a DOF the accelerated query classified fluid while the exhaustive one found it
solid (or vice versa). Verdict stages expect none. Control stages run with a deliberately stale
index and must *find* misses — otherwise the placement has missed the cell lattice and the
passing stages prove nothing; the gate then fails the run. Each `DNS_QUERY_GATE` line carries
`misses=`, `owner_mismatches=`, `expected_misses=` and `verdict= pass | control_confirmed | FAIL`.
A stage that compared fewer than 1000 classifications is vacuous regardless of its counts.

Exit codes: `0` green · `62` a verdict stage failed · `63` the gate compared nothing · `64` the
preflight failed (grid inactive or mis-sized on some rank, fewer than 33 bodies, or
`SimPar@UseHashGridAccel = No`).

## The construction (why the bodies look the way they do)

pe files a body in the grid whose cell span is `√hierarchyFactor = √2` times the body's AABB. A
population of equal bodies therefore fills only 70.7 % of its cells, straddles at most two
adjacent cells, and the 27-cell stencil always covers the cell a body was filed in — **no miss is
possible**. The defect needs a body whose AABB fills a whole span, which only happens when
*smaller* bodies set that span. So:

- **40 spacers**, spheres of radius 0.06, created **first** and parked outside the mesh. They set
  `span = (2r + 2·contactThreshold)·√2 = 0.16970564` and trip the 32-body activation threshold.
  Order matters: only the first 32 pending bodies are re-filed when the threshold trips, and the
  first of them creates the grid.
- **30 probes**, cycling sphere / box / capsule / ellipsoid / cylinder, each sized by its x
  **half-extent** `ex = span/2 − contactThreshold − 1e-5 = 0.08484282`, so its AABB is
  `span − 2e-5`. Not by bounding radius: a cylinder inscribed in a sphere has an x-extent of only
  `√3/2` of the diameter and would never fill the span — the artefact that made Layer 1's cylinder
  look immune. pe asserts every probe's real AABB x-extent at creation.
- Probes are placed by their **real** AABB minimum, `1e-5` below the planes `k·span`, on a 5×2×3
  lattice with one z-row per subdomain. No probe is ever fixed: both queries skip fixed cylinders
  and meshes.
- The nudge is **one Q2 DOF pitch**, `1/36`, not Layer 1's few-`δ` shove. Layer 1 chose its query
  points just inside the far surface; here the query points are the mesh's DOF planes `x = m/36`,
  and a miss needs a DOF plane *beyond* the next cell boundary to lie inside the body. A smaller
  nudge leaves the far surface inside the same cell and the stencil still covers the stale one —
  the first run of this harness compared 53 391 DOFs after a `2e-5` nudge and, correctly, found
  nothing, which is what the `A2` control exists to catch.

The preflight reads the live span back from pe and aborts if it is not the computed one, so a
change to pe's `hierarchyFactor`, AABB padding or re-filing order fails the harness instead of
silently defusing it.

## Multiple ownership — the CSR exception table

Overlapping bodies are real (contact), and a DOF inside two of them has two owners that the
velocity BC and the force attribution cannot both honour; today they silently take one — the
minimum system id, since pe #43. `SimPar@FbmOwnerTable = Yes` makes that ambiguity visible and
countable without touching the hot path (module `fbm_ownership`, `source/src_fbm/`):

- the **canonical owner stays dense and unchanged** — `FictKNPR` (flag) and `FictKNPR_uint64`
  (minimum id), one per DOF;
- only DOFs inside **two or more** bodies get a row: `FictOwn_row(i)` (one integer per DOF, 0 =
  fluid or single owner, else the row), then CSR `FictOwn_dof(r)`, `FictOwn_ptr(r)`,
  `FictOwn_id(k)` with ids ascending — so a run with no overlap allocates no payload;
- invariant the gate asserts: the first id of every row equals the canonical owner, so a
  consumer that ignores the table sees exactly today's behaviour.

pe fills it through two new owner-list queries (`pointOwners`, `pointOwnersAccelerated`), the
owner-list forms of the flag queries, built on the same containment helper. With the gate on,
the harness compares the two **owner sets** per DOF — the same minimum id can hide a second
owner that only one path sees — and reports `owner_set_mismatches=` on the `DNS_QUERY_GATE`
line (a verdict stage fails on any). After each stage a `DNS_QUERY_OWNERS` line gives
`multi_owned_dofs=`, `nnz=`, `max_owners_this_rank=`, `overflow=` (more than 16 owners on one
DOF), `canonical_violations=` and a verdict; stage `E` additionally expects
`multi_owned_dofs > 0`, or its overlap construction reached no DOF.

## Deck keys

`SimPar@GateStageList = "A1,A2,G,R,E,P"` (quoted — the parser splits unquoted commas),
`SimPar@GateRefreshIndex = Yes` (`No` turns every verdict stage into a control; a debugging knob),
`SimPar@GateProbeReport = Yes` (pe prints `FBM_GATE_CONSTRUCTION` and one `FBM_GATE_PROBE` line
per probe), `SimPar@FbmOwnerTable = Yes` (fill and check the multiple-ownership table). `SimPar@UseHashGridAccel` must be `Yes`, and `stepsize_` in `pe_user_config.json`
must equal `SimPar@TimeStep` or the first particle step aborts.

## Where things live

Geometry and the C entry points: `libs/pe/pe/interface/fbm_gate_probe.h`,
`libs/pe/src/interface/fbm_gate_probe.cpp`, the `GateProbe` branch of `setupHashGridTest`.
Verdict, preflight, counters: `source/src_fbm/fbm_main.f90`. Fortran bindings:
`source/src_particles/dem_query.f90`. Stage order: `hashgrid_gate_driver.f90` here.
