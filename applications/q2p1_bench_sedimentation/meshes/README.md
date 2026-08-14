# Mesh fixtures for the ten Cate sedimentation benchmark (DNS campaign)

Committed as part of the DNS validation campaign, stage D0.4
(`dns-validation-campaign-plan.md`). `.gitignore` excludes `*.tri`/`*.par`
globally, so these were added with `git add -f`; use `-f` again when
modifying them.

## `benchSym/` — quarter-box symmetry mesh

Coarse mesh of the ten Cate (2002) box exploiting two symmetry planes:
`x ∈ [-0.05, 0]`, `y ∈ [-0.05, 0]`, `z ∈ [0, 0.16]` (m); 876 elements /
1239 vertices, 36 uniform z-layers.

- `mesh.tri` + `bench.prj` + boundary `.par` files
  (`top/bot/xwall/ywall/x/y.par`). Boundary conditions live only in the
  `.par` headers.
- The symmetry planes carry no-penetration BCs; the sphere sits on the
  corner axis, so only a quarter of it is simulated. The hydrodynamic
  force must therefore be scaled ×4 in z and zeroed in x/y — runtime key
  `Prop@ForceScale = 0,0,4d0,...` in `_data/q2p1_param.dat` (historically
  also the compile-time `-DSED_BENCH` hack, retired by campaign item
  D0.3).
- **Not usable for EL runs** (the EL kernel needs the full box; see the
  EL v1b RUNBOOK).

## `benchSym/mesh12/NEWFAC/` — 1×1×12 Cartesian partition

Pre-partitioned layout for **parallel-PE** mode (np = 13: 1 CFD master +
12 workers), matching PE's `MPI_Cart_create` decomposition
`processesX_/Y_/Z_ = 1/1/12` in `example.json`. `GRID.tri`/`GRID.prj` +
`sub0001…sub0012/GRID0001.tri` with per-subdomain boundary files.

Regeneration recipe (legacy METIS-4 Cartesian path — METIS 5 graph
partitioning cannot produce ordered axis-aligned slabs):

```bash
python3 tools/partpy/PyPartitioner.py 12 -4 x1-y1-z12 NEWFAC _adc/benchSym/bench.prj
```

`pz` must divide the z element-layer count (36); see
`docs/md_docs/guide_07_q2p1_bench_sedimentation_pe_parallel_from_scratch.md`
§6 for the sizing rule and cover check.

## Full-box mesh

The full (non-symmetry) ten Cate box + its 3×3×3 27-subdomain partition
is committed under
`applications/q2p1_el_pipeflow/validation_cases/v1b_tencate_settling/mesh/ten_cate_box`
(owned by the EL campaign). The quarter box and the full box are
different discretizations of the same experiment — every DNS/EL
cross-comparison must state which one it used (campaign risk #5).

## Physical setup (ten Cate 2002, case E4 unless stated)

Box 100×100×160 mm, sphere d = 15 mm, ρ_p = 1120 kg/m³, release height
z = 0.1275 m; E1–E4 fluids (ρ_f, μ): (970, 0.373), (965, 0.212),
(962, 0.113), (960, 0.058) → Re_∞ = 1.5 / 4.1 / 11.6 / 31.9. Reference:
`ten_cate_piv.pdf` (repo root).
