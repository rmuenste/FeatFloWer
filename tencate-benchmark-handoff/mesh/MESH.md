# MESH — the two coarse meshes, refinement and partitioning

Two discretisations of the same experiment exist in the FeatFloWer repository
(campaign plan §11 item 5: "different discretizations; every cross-comparison
states which is used"). **All certified results in README §6 were computed on
the quarter box (`quarterbox_benchSym`)**; the full box kit was built for the
campaign but no certified row uses it (README §9 item 1).

Mesh files are FeatFloWer/FEAT `.tri` coarse-grid files: header line 3 gives
`NEL NVT NBCT NVE NEE NAE` (elements, vertices, boundary components, 8 vertices
per hexahedron, 12 edges, 6 faces); `DCORVG` = vertex coordinates [m]; `KVERT`
= element connectivity; `KNPR` = vertex boundary flags. The `.par` files list,
per boundary patch, the number of vertices, a type tag (`Wall`, `Outflow`,
`Symmetry100`, `Symmetry010`) and a parametrisation line
`'4 a b c d'` = plane a·x + b·y + c·z + d = 0, followed by the vertex indices.
The `.prj` file lists the mesh and its `.par` files. Each `.tri` has a `.vtk`
/ `.vtu` sibling for inspection in ParaView.

Refinement rule (both meshes): regular subdivision — every hexahedron splits
into 8 per level; level L1 = the coarse mesh, h(L) = h(L1)/2^(L−1). The campaign
reports resolution as D/h with h = (finest-level minimum element volume)^(1/3)
(`DNS_RESOLUTION … h_min … D_over_h` log lines; guide §2). Velocity is Q2
(nodal spacing h/2), pressure P1 discontinuous.

## quarterbox_benchSym/ — the certified mesh

- Files: `mesh.tri` (876 hexahedra, 1239 vertices), `bench.prj`, `top.par`,
  `bot.par`, `xwall.par`, `ywall.par`, `x.par`, `y.par`, `grid.vtu`. Origin:
  `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/` (the `_adc`
  symlink target of every certified rundir; byte-identical copies live at the
  repo root `benchSym/` and in each rundir as `_mesh/NEWFAC/GRID.tri`).
  Header: "Coarse mesh exported by DeViSoR TRI3D exporter"; file date
  2018-03-02. The generator input/script is not in the repository.
- Domain: x ∈ [−0.05, 0], y ∈ [−0.05, 0], z ∈ [0, 0.16] m — one quarter of the
  100 × 100 × 160 mm box, cut by the two vertical mid-planes.
- Boundary patches (`*.par`): `x.par` 261 vertices `Symmetry100` on x = 0;
  `y.par` 261 vertices `Symmetry010` on y = 0; `xwall.par` 78 vertices `Wall`
  on x = −0.05; `ywall.par` 78 vertices `Wall` on y = −0.05; `bot.par` 51
  vertices `Wall` on z = 0; `top.par` 51 vertices `Outflow` on z = 0.16 (the
  deck sets `SimPar@NoOutflow = Yes`; the experiment has a free surface there).
- Grading: z-planes every 4.444 mm (36 layers) in the core region around the
  axis, and every 13.33 mm (every third plane carries the full 51-vertex
  cross-section) in the outer region; lateral spacing ≈ 1.9–3.5 mm next to the
  axis (x, y ∈ [−0.02, 0]) growing to 10–14 mm at the far walls (vertex
  coordinate inventory of `mesh.tri`; inspect `grid.vtu`). At L3 the finest
  cell is h_min = 6.27e-4 m and the coarsest in-bounds cell ≈ 1.87e-3 m
  (`DNS_MESH_PROBE` log lines), i.e. the mesh is ~3× finer around the sphere
  path than in the far field.
- Sphere relative to the mesh: centre on the corner edge x = y = 0 (the
  symmetry axis), z = 0.1275 m at release; the sphere path z ∈ [0.0075, 0.1275]
  lies in the fine core; the bottom wall z = 0 is a mesh plane.
- Element counts per level (876 · 8^(L−1)) and measured h_min:

  | level | elements (quarter box) | h_min [m] | D/h | note |
  |---|---|---|---|---|
  | L1 | 876 | — | — | coarse grid |
  | L2 | 7 008 | 1.31425429e-3 | 11.413 | smoke grade (+3.4…+4.4 %) |
  | L3 | 56 064 | 6.27151451e-4 | 23.918 | workhorse |
  | L4 | 448 512 | 3.05512350e-4 | 49.098 | 60 GB on 32 ranks |

- Partition used by the certified runs: the 876-element coarse mesh split
  into **31 subdomains** (`_mesh/NEWFAC/sub0001/GRID0001…GRID0031.tri`,
  ≈ 28 coarse elements each) by the FeatFloWer partitioner ("Coarse mesh
  exported by Partitioner", log: `RecursivePartitioning = YES`,
  `PartitionFormat = legacy`), run with 32 MPI processes (1 master + 31). Deck:
  `SimPar@SubMeshNumber = 1`, `SimPar@MeshFolder = "NEWFAC"`. The partition
  kit (233 files, 1.2 MB per rundir) is FeatFloWer-specific and **not copied**;
  available on request. A second, Cartesian 1 × 1 × 12 partition of this mesh
  (`benchSym/mesh12/NEWFAC`) is mentioned in plan §D0.4 and was not inspected.

## fullbox_ten_cate_mesh_v1/ — the full-box kit (not used for certified rows)

- Files: `box.tri` (324 hexahedra, 490 vertices), `box.vtk`, `file.prj`,
  `xmin/xmax/ymin/ymax/zmin/zmax.par`. Origin: repo root `ten_cate_mesh_v1/`
  (dated 2026-07-04). Header: "Coarse mesh exported by hex_ex.py".
- Domain: x ∈ [0, 0.1], y ∈ [0, 0.1], z ∈ [0, 0.16] m — the full
  100 × 100 × 160 mm box. Structured 6 × 6 × 9 grid: cells 16.667 × 16.667 ×
  17.778 mm (x/y spacing 0.1/6, z spacing 0.16/9 — slightly anisotropic).
- Boundary patches: all six faces type `Wall` with plane parametrisation
  (`xmin.par`: `'4 1d0 0d0 0d0 0d0'` 70 vertices … `zmax.par`:
  `'4 0d0 0d0 1d0 -0.16d0'` 49 vertices). The top is a wall here (closed box),
  unlike the quarter-box `Outflow` tag.
- Sphere relative to the mesh: the experiment's sphere centre would be
  (0.05, 0.05, 0.1275): x = y = 0.05 are coarse grid lines (3 × 16.667 mm);
  z = 0.1275 lies between the coarse planes 0.1244 and 0.1422.
- Resolution per level: h = 16.7–17.8 mm / 2^(L−1) → L3 ≈ 4.2 mm (D/h ≈ 3.6),
  L4 ≈ 2.1 mm (D/h ≈ 7), L5 ≈ 1.04–1.11 mm (D/h ≈ 13.5–14.4), L6 ≈ 0.52–0.56 mm
  (D/h ≈ 27–29). Matching the certified quarter-box L3 (D/h 23.9) therefore
  needs about L6 on this kit (324 · 8^5 ≈ 10.6 million elements for the full
  box vs 56 064 for the quarter box at L3) — the kit was sized for a coarse
  Cartesian partition, not for the certified ladder. Uniform hexahedra, no
  grading.
- Partition kit in the repo: `_mesh/ten_cate_box/` with 27 subdomains
  (`sub0001…sub0027`, 12 coarse elements each = axis-uniform 3 × 3 × 3 blocks
  of 2 × 2 × 3 cells; `GRID.tri`, `GRID.prj`, per-subdomain `.par` files,
  `subdomains.pvd`/`main.vtu`), 2.1 MB; run with 28 processes (1 + 27). **Not
  copied** (FeatFloWer-specific); available on request.
- Generator: `generator/hex_ex.py` (+ its `mesh/` package `mesh.py`,
  `mesh_io.py`; `simple_ex.py`; `tri2vtk_converter.py`) from the repo's
  `tools/quadextrude/`. It reads a 2D quad mesh in gmsh `.msh` format and
  extrudes it in z in levels/layers (`-f mesh-file -l levels -e extrusion-layers
  -d distance-levels -i ids-level -o output-dir`), writing `mesh.tri`, the
  `.par` files and VTK previews. The exact invocation and the 6 × 6 quad input
  that produced `box.tri` are not recorded in the repository (gap). Python 2/3
  era script; `from mesh import *` requires running from the generator
  directory.

## For a new code

The physical box is 100 × 100 × 160 mm with the sphere on the vertical axis.
A full-box run avoids the symmetry assumption (PITFALLS P20); a quarter box
with two symmetry planes is what the FeatFloWer numbers correspond to. State
D/h in your own convention (cells vs nodes across the diameter) and the
smallest cell on the sphere path; the FeatFloWer ladder is at 11.4 / 23.9 /
49.1 Q2 elements per diameter.
