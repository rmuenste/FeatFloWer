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

## Expected local layout (for staging)

```
mesh/
  quiescent_box_9/        file.prj, box.tri, {x,y,z}{min,max}.par   (coarse project)
  partitions/
    QBOX9/
      sub0001/ … sub0027/  GRID.tri, GRID0001.tri, *.par           (per-subdomain)
```

The stage target copies `mesh/quiescent_box_9 → <rundir>/quiescent_box_9` and
`mesh/partitions → <rundir>/_mesh`, so the runtime reads
`_mesh/QBOX9/sub<rank>/GRID0001.tri` directly (no `./partitioner` step).

> Note: this case currently exercises the periodic two-way coupling but does not
> yet conserve momentum — the periodic box's pressure nullspace is handled with
> `NoOutflow=Yes`, which behaves enclosed and absorbs the net feedback momentum.
> A nullspace-aware (zero-mean-pressure) solve is needed for true conservation.
