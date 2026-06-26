# Mesh for the momentum-conservation Tier-2 case

The mesh is **not committed** (large; fetched from the cloud). Place it under
`mesh/` in this directory and the
`q2p1_el_pipeflow_tier2_momentum_conservation_stage` target copies it into the
run dir automatically.

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

## Expected local layout (for staging)

```
mesh/
  quiescent_box_9/        file.prj, box.tri, {x,y,z}{min,max}.par   (coarse, faces tagged Periodic)
  partitions/
    QBOX9/
      sub0001/ … sub0027/  GRID.tri, GRID0001.tri, *.par           (per-subdomain, faces tagged Periodic)
```

The stage target copies `mesh/quiescent_box_9 → <rundir>/quiescent_box_9` and
`mesh/partitions → <rundir>/_mesh`, so the runtime reads
`_mesh/QBOX9/sub<rank>/GRID0001.tri` directly (no `./partitioner` step).

## Status / expected result

With the no-wall (`Periodic`-tagged) mesh the box conserves momentum to leading
order: as the particle decelerates, the fluid mean flow gains the momentum
(verified over 60 steps, e.g. particle `p_z` 7.6e-6 → 0.2e-6 while fluid `p_z`
1.0e-6 → 10.1e-6). A residual **~15–20 % drift remains** — this is the
explicit/implicit drag inconsistency ("issue D"): the fluid receives the
explicitly-evaluated drag while the particle is advanced with the smaller
semi-implicit drag, so the exchange does not balance exactly. The drift is
bounded and plateaus as the drag → 0.

The test therefore asserts `drift` below a tolerance set to this issue-D band
(it quantifies the inconsistency); it cannot assert `drift → 0` until the future
"spread the **implicit** drag reaction to the fluid" fix lands.
