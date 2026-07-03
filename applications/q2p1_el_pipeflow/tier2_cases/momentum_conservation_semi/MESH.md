# Mesh for the semi-implicit momentum-conservation Tier-2 case

This case reuses the explicit momentum-conservation QBOX9 mesh. The CMake
`q2p1_el_pipeflow_tier2_momentum_conservation_semi_stage` target copies the mesh from
`tier2_cases/momentum_conservation/mesh`; this directory intentionally carries only
the semi-implicit configuration files.

## Required mesh

- **Coarse geometry:** unit cube `[0,1]³`, structured `9³` hex mesh
  (`box.tri` + `file.prj` + the six `*.par` boundaries), refined to
  `MaxMeshLevel = 2` ⇒ `18³` at run time. Source name: `quiescent_box_9`.
- **Decomposition:** `QBOX9` — a regular `3 × 3 × 3` partition,
  **27 subdomains** ⇒ run with **`mpirun -np 28`** (1 master + 27 workers).
  This matches `example.json` (`processesX_ = processesY_ = processesZ_ = 3`)
  and `q2p1_param.dat` (`MeshFolder = "QBOX9"`, `SubMeshNumber = 27`).
  A regular 3×3×3 grid is required here so the PE decomposition aligns with the
  CFD partition and so the periodic comm (`SimPar@Periodic = Yes`) is stable —
  the z-only `QBOX9Z3` split crashes `E013_CreateComm` under periodicity.

## Boundary conditions — CRITICAL for momentum conservation

The six box faces **must NOT be tagged `Wall`** in the `.par` files. A `Wall`
tag sets a no-slip velocity Dirichlet (`bWall → knpr=1`, `u=0`) on the face,
which pins the fluid and makes the box a closed container: the net two-way
feedback momentum is absorbed at the walls and total momentum is **not**
conserved (observed drift ≈ 100 %).

Instead tag the faces with any type the parametrizer does **not** recognise as a
velocity BC — `Periodic` works (it sets neither `bWall` nor `knpr`, and does not
collide with the `Inflow`/`Temperature`/`Phase`/`Bubble` prefix parsing in
`source/Parametrization.f90`). The faces then stay free, and with
`SimPar@Periodic = Yes` the `dPeriodicity` comm couples opposite faces into a
single periodic unknown, so the fluid can carry a mean flow. No solver code
change is needed — this is purely a `.par` header retag (`Wall` → `Periodic`)
on every `.par` of the coarse mesh and of each subdomain.

`SimPar@NoOutflow = Yes` is kept to anchor the periodic pressure nullspace; once
the velocity is genuinely periodic (no walls), pinning the pressure constant is
momentum-neutral.

## Expected shared mesh layout

```
tier2_cases/momentum_conservation/mesh/
  quiescent_box_9/        file.prj, box.tri, {x,y,z}{min,max}.par   (coarse, faces tagged Periodic)
  partitions/
    QBOX9/
      sub0001/ … sub0027/  GRID.tri, GRID0001.tri, *.par           (per-subdomain, faces tagged Periodic)
```

The semi stage target copies `mesh/quiescent_box_9 → <rundir>/quiescent_box_9` and
`mesh/partitions → <rundir>/_mesh`, so the runtime reads
`_mesh/QBOX9/sub<rank>/GRID0001.tri` directly (no `./partitioner` step).

## Status / expected result

With the no-wall (`Periodic`-tagged) mesh the box should conserve momentum to leading
order with `SimPar@ELDragCoupling = semi_implicit`. The regression metric is
`EL_MOMENTUM_ELEMINT drift_rel`, which uses an element-integrated fluid momentum
and avoids shared-DOF double counting. The Tier-2 test asserts
`EL_MOMENTUM_ELEMINT drift_rel < 3.1e-3` and
`EL_FEEDBACK_CONSERVATION residual < 1.0e-10`.

The case sets `SimPar@ELMomentumAuditFreq = 1` so the permanent low-frequency
momentum audit is emitted every step. Production runs can use the default
cadence (`100`) to keep a conservation tripwire without verbose logs.
