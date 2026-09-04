# chimera_meshgen — coarse shell meshes for the Chimera component

Generates coarse body-fitted "atmosphere" meshes (.tri, readTriCoarse
format) for the Chimera overlapping-mesh component
(`chimera-integration-design.md`, section 7). The refinement ladder comes
from the in-code hierarchy (`ChimeraSubmeshLev`), not from multiple
files: generate ONE coarse mesh, refine at load time (with per-level
radial boundary projection, handled by `chi_legacy_mesh_adapter`).

Boundary vertices are placed EXACTLY on the analytic surfaces — the
level-1 geometric classification in `chi_submesh.f90` relies on that.

## annulus_tri.py — z-aligned annulus (cylinder atmosphere)

O-grid ring, extruded in z, periodically closed in theta; FEAT hex
vertex ordering.

```bash
python3 annulus_tri.py --ri 1.0 --ro 2.0 --lz 0.5 --nr 2 --nt 12 --nz 1 \
        --out annulus_coarse.tri
```

The committed test fixture
`source/src_chimera/tests/fixtures/annulus_coarse.tri` is exactly this
call (r_i=1, r_o=2, L_z=0.5, nr=2, nt=12, nz=1; 24 hexes) and is the
mesh of the `chi-submesh-couette` Phase-2 gate test. Regenerate it only
together with the constants in `tests/test_chi_submesh.f90`.

## sphere_shell_tri.py — cubed-sphere shell (sphere atmosphere, Phase 5)

The six faces of the cube are mapped equiangularly onto the sphere
(near-uniform surface cells) and stacked in `nr` radial layers
(`--grading` = geometric growth factor per layer); closed, conforming
shell with `6 n^2 nr` hexes, FEAT vertex ordering, positive Jacobians.
One file serves a whole array: at load time each body's copy is fitted to
`[R, R + H_k]` and translated to the body centre.

```bash
python3 sphere_shell_tri.py --ri 1.0 --ro 2.0 --n 4 --nr 2 --out sphere_shell_coarse.tri
```

This call produced `source/src_chimera/tests/fixtures/sphere_shell_coarse.tri`
(192 hexes; 1536 at submesh level 2, 13842 Q2 dofs) used by the Hasimoto
case (`applications/q2p1_chimera/_data/q2p1_param_hasimoto*.dat`).

## seed_array.py — sphere arrays in a periodic box (Phase 5)

Emits the `SimPar@ChimeraParticleFile` body table for a simple-cubic
lattice (`--mode sc --n 1` = Hasimoto; `--offset 0,0,0` puts the sphere
on the box corner, the periodic-donor check) or a random sequential
addition (`--mode random --count N --phi --seed --mingap`), with the
non-overlap rule `H_k = min(H_max, 0.5 * nearest surface gap)` evaluated
with minimum-image distances (own periodic images count). Unlike
`tools/gen_random_array.py` (FBM) no explicit image spheres are written:
the Chimera geometry is minimum-image periodic itself.

```bash
python3 seed_array.py --mode sc --n 1 --phi 0.0193925 --hmax 1.0 --out hasimoto_center.dat
python3 seed_array.py --mode random --count 8 --phi 0.05 --seed 1 --mingap 0.5 --out random8.dat
```

## flatten_axis_partition.py — single-host layout of an axis partition

Periodic runs need an axis-aligned partition; `PyPartitioner.py 1 -123 8`
produces it as eight SUBGRIDS (`sub000k/GRID.tri`), while a single-host
run reads `sub0001/GRID000k.tri`. The script rewrites the former into
the latter (see its docstring and `docs/md_docs/chimera_usage.md`).

Phase-3 note: the FAC cylinder-annulus instantiation uses the annulus
fixture above.

## channel_tri.py — structured box background mesh (+ .par/.prj)

Uniform hex box for Chimera cases (no body in the mesh), with the DFG
2D-FAC boundary types (`Inflow2` at x=0, `Symmetry011` at x=Lx, `Wall`
at y=0/Ly, `Symmetry001` at z=0/Lz) and a project file for
`tools/PyPartitioner.py`; with `--periodic` all six faces are tagged
`Periodic` (`xmin..zmax.par`) for a periodic box
(`applications/q2p1_chimera/_data/CHIMERA_BOX6`: 6^3 unit cube).

```bash
python3 channel_tri.py --lx 2.2 --ly 0.41 --lz 0.05 --nx 44 --ny 8 --nz 1 \
        --name chimera_fac --outdir ../../applications/q2p1_chimera/_data/CHIMERA_FAC
```

This call produced the vendored milestone-1 case
`applications/q2p1_chimera/_data/CHIMERA_FAC` (352 hexes, h = 0.05 at the
coarse level; 0.025 at `MaxMeshLevel = 2`, 0.0125 at level 3).  The
atmosphere for that case is the Phase-2 annulus fixture, fitted at load
time to r_i = 0.05, r_o = 0.1 (`CHI_SUBMESH_FIT_COARSE`).
