# DKT mesh set (drafting–kissing–tumbling)

Box **6d × 6d × 24d**, closed container, built for the DKT resolution ladder
**D/h ∈ {8, 16, 32}**.

Generated 2026-08-02 with the `meshing_component/gmsh-learning2` toolkit
(`unit_cube_vtk.py` → `tri2vtk_converter.py` → `gen_par_from_tri.py` →
`tag_box_boundaries.py` → `PyPartitioner.py`). No bespoke `make_box_tri.py` was
needed — `unit_cube_vtk.py` already is that script.

> **Units: d = 1.** The meshes are nondimensional in the sphere diameter, so the
> bounds are x, y ∈ [−3, 3] and z ∈ [0, 24]. **Set the physical `d` before use**
> and rescale all six bounds by it (regenerate; do not scale the `.tri` by hand,
> because the `.par` plane descriptors carry the bounds too).

## Contents

```
coarseB/                12 x 12 x 48 =  6,912 el   coarse cell = d/2   D/h_coarse = 2
├── dkt_B.tri, file.prj, {x,y,z}{min,max}.par, main.vtu
└── _mesh/
    ├── DKT_B_108/      x3-y3-z12   108 subdomains,  64 coarse el each
    └── DKT_B_864/      x6-y6-z24   864 subdomains,   8 coarse el each
coarseC/                24 x 24 x 96 = 55,296 el   coarse cell = d/4   D/h_coarse = 4
├── dkt_C.tri, file.prj, {x,y,z}{min,max}.par, main.vtu
└── _mesh/
    └── DKT_C_864/      x6-y6-z24   864 subdomains,  64 coarse el each
```

All partitions verified: uniform element counts, complete cover, zero uneven
subdomains. Written with `-f v2`, i.e. `sub%04d/GRID%04d.tri`, which is what
`source/include/PartitionReader.f90` builds.

Boundaries: all six faces `Wall` with quoted type-4 plane descriptors (closed
container; `NoOutflow` handles the enclosed pressure). The quoting is required —
the CFD reads that line with a Fortran list-directed READ, which stops at the
first blank in an unquoted string.

### Header differences from benchSym (checked, benign)

These meshes differ from the DeViSoR-exported `benchSym/mesh.tri` in two header
fields. Both were traced and neither matters:

| field | benchSym | here | why it is fine |
|---|---|---|---|
| `NBCT` | 1 | 0 | `%nbct` has no consumer. Read at `mesh_refine.f90:245`, echoed back at `:37`, hard-set to 1 in the JSON path at `:1007`. Never enters a computation. |
| `KNPR` | 652 nonzero | all 0 | `mesh_refine.f90:1492` does `mesh1%knpr = 0` on every refined level, so coarse values are discarded at the first refinement. With `NLmax >= 2` they cannot matter. |

Boundary conditions come from the `.par` files via `Parametrization`, not from
`KNPR`. FeatFloWer itself writes `NBCT = 0` in one of its own export paths
(`mesh_refine.f90:101`).

Layout **does** matter and was verified: `readTriCoarse` reads `KNPR` with one
`read` per vertex (`do i=1,nvt; read(cunit,*) knpr(i)`), so the values must be
one per record. All files here — coarse and every partitioned subdomain — are
one per line.

## Resolution ladder

`D/h` counts **Q2 elements** per diameter, with `h = dvol^(1/3)` — the same
convention as `DNS_RESOLUTION` and `ComputeCFL`, and as the D1.1 runs. Velocity
nodal spacing is `h/2`, so nodal numbers are double these. LBM/IBM literature
values are nodal: **halve them before comparing.**

Because {8, 16, 32} is a doubling sequence and regular refinement halves `h`,
this is **one coarse mesh at three consecutive levels**, not three meshes.

| D/h | elements | level on B | level on C | vel. dofs |
|---|---|---|---|---|
| 8 | 442,368 | L3 | L2 | ~10.6 M |
| 16 | 3,538,944 | L4 | L3 | ~84.9 M |
| 32 | 28,311,552 | L5 | L4 | ~679.5 M |

## Partitioning plan

| rung | split | ranks | el/rank |
|---|---|---|---|
| D/h = 8 | `x3-y3-z12` | 108 | 4,096 |
| D/h = 16 | `x3-y3-z12` | 108 | 32,768 |
| D/h = 32 | `x6-y6-z24` | 864 | 32,768 |

Subdomains are 32×32×32 elements at the top rung — exactly cubic, AR 1.00, and
34% under the 50k-per-rank target rather than scraping it at 49,152.

The nicer property: **per-rank load is identical at D/h = 16 and D/h = 32**
(32,768 both). The decomposition doubles per axis exactly as the resolution
doubles, so the jump is pure weak scaling. If the old time-per-step guideline
still holds at all, the D/h = 32 step time should land close to the D/h = 16 step
time — which makes the guideline testable rather than assumed.

Coarse cells per subdomain drops to **8 (2×2×2)** on coarse B at 864 ranks.
That is thinner than benchSym's 73 or ten Cate's 12, but not degenerate — it just
means each rank's local hierarchy bottoms out at 2×2×2, with the global coarse
solve still on all 6,912. If that turns out too thin in practice, **coarse C
gives 64 per rank at the same 864 ranks**, at the cost of a 55,296-element coarse
solve and a shallower L2/L3/L4 ladder. That is the only reason `coarseC/` exists
here; it is a fallback/test, not the primary path.

Budget: **864 subdomains + 1 master = 865 MPI ranks.** Not feasible on
nx-01..06 (6 nodes × 64 cores = 384 cores total) — the D/h = 32 rung has to move
to a different cluster.

### Why not `x5-y5-z20`

It looks like a valid AR = 1.00 split of a 6:6:24 box (6/5 = 24/20 = 1.2), but it
is unrealisable: the fine element counts are `6·D/h = 2^(m+1)·3`, which never
carries a factor 5. At D/h = 32 the mesh is 192 × 192 × 768 and 192/5 = 38.4.
The AR = 1.00 family for this box is `ranks = 4k³` with `k | 6·D/h`, i.e.
108, 256, 864, 2048 — nothing between 256 and 864.

## One free check worth taking

The resolution progression stays consistent — the ladder is the mesh, not the
decomposition. But the decomposition changes at exactly the rung where the
physics gets most expensive, so anything decomposition-dependent (PE domain
assignment, halo exchange, FBM force accumulation across subdomain boundaries)
shifts simultaneously with resolution. On the r_h curve's top point, that is a
confound you would rather not have to argue about.

It is cheap to eliminate: **`x6-y6-z24` also divides the D/h = 16 mesh**
(96/6 = 16, 384/24 = 16). Same partitioned `_mesh` directory, only `NLmax`
differs. So run D/h = 16 on both 108 and 864 ranks and confirm they agree, which
pins decomposition-invariance at 3.5 M elements before spending anything at
28 M. If they agree, the D/h = 32 point inherits the credibility.

## Regenerating

```bash
module load python/3.13.5
GL2=<repo>/meshing_component/gmsh-learning2

$GL2/.venv/bin/python $GL2/pe_partpy/unit_cube_vtk.py --out /tmp/dkt_B.vtk \
  --nx 12 --ny 12 --nz 48 --x0 -3 --x1 3 --y0 -3 --y1 3 --z0 0 --z1 24
$GL2/.venv/bin/python $GL2/pe_partpy/tri2vtk_converter.py /tmp/dkt_B.vtk /tmp/dkt_B.tri
$GL2/.venv/bin/python $GL2/pe_partpy/gen_par_from_tri.py /tmp/dkt_B.tri --outdir coarseB
$GL2/.venv/bin/python $GL2/workflow/tag_box_boundaries.py --case-dir coarseB
cd coarseB
$GL2/.venv/bin/python $GL2/pe_partpy/PyPartitioner.py 108 axis_uniform x3-y3-z12 DKT_B_108 file.prj -f v2
$GL2/.venv/bin/python $GL2/pe_partpy/PyPartitioner.py 864 axis_uniform x6-y6-z24 DKT_B_864 file.prj -f v2
```

Coarse C is the same with `--nx 24 --ny 24 --nz 96` and only the 864 split.

Notes on the toolkit:
- The `x?-y?-z?` spec must equal the `(px, py, pz)` passed to PE's
  `decomposeDomain`. Subdomain numbering loops x outermost and z innermost, so
  **z varies fastest**: `subId = ((ix-1)·partY + (iy-1))·partZ + iz`.
- Numeric strategy codes (`-4`) are no longer accepted — use `axis_uniform`.
- `workflow/run_partition_to_vtu.py` always writes to `gmsh-learning2/_mesh/`
  and must be run from that repo root; `PyPartitioner.py` writes `_mesh/`
  relative to the current directory, which is why it is used directly above.
- Requires pe_partpy at or after the `axis_uniform` boundary-precision fix
  (rmuenste/pe_partpy#1). Earlier revisions silently mis-partition meshes whose
  `.tri` was written with few decimals; these meshes are full precision and are
  not affected, but a regenerated or hand-edited mesh could be.

## Wall times (this machine, for planning)

| step | wall |
|---|---|
| coarse B partition, 108 | 4.4 s |
| coarse B partition, 864 | 32.3 s |
| coarse C partition, 864 | 39.4 s |

Disk: coarseB 62 MB, coarseC 82 MB, 25,746 files total.
