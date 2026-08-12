# Mesh for the terminal-velocity Tier-2 case

The mesh is **not committed** (large; fetched from the cloud). Place it under
`mesh/` in this directory and the `q2p1_el_pipeflow_tier2_terminal_velocity_stage`
target copies it into the run dir automatically.

## Required mesh

- **Coarse geometry:** unit cube `[0,1]³`, structured `9³` hex mesh
  (`box.tri` + `file.prj` + the six `*.par` wall boundaries), refined to
  `MaxMeshLevel = 2` ⇒ `18³` at run time. Source name: `quiescent_box_9`.
- **Decomposition:** `QBOX9Z3` — a `1 × 1 × 3` (z-split) partition,
  **3 subdomains** ⇒ run with **`mpirun -np 4`** (1 master + 3 workers).
  This matches `example.json` (`processesZ_ = 3`) and `q2p1_param.dat`
  (`MeshFolder = "QBOX9Z3"`, `SubMeshNumber = 3`).

## Expected local layout (for staging)

```
mesh/
  quiescent_box_9/        file.prj, box.tri, {x,y,z}{min,max}.par   (coarse project)
  partitions/
    QBOX9Z3/
      sub0001/ … sub0003/  GRID.tri, GRID0001.tri, *.par           (per-subdomain)
```

The stage target copies `mesh/quiescent_box_9 → <rundir>/quiescent_box_9` and
`mesh/partitions → <rundir>/_mesh`, so the runtime reads
`_mesh/QBOX9Z3/sub<rank>/GRID0001.tri` directly (no `./partitioner` step).
