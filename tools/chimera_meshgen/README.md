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

Planned (design roadmap): sphere shells (Phase 5 array closures; the
`applications/mesh_ref/Particle/*.tri` 24-hex shells serve as interim
fixtures), per-particle atmosphere width tables
`H_k = min(H_max, 0.5 * nearest-surface gap)` for the seeding rule, and
the FAC cylinder-annulus instantiation for the Phase-3 milestone.

## channel_tri.py — structured box background mesh (+ .par/.prj)

Uniform hex box for Chimera cases (no body in the mesh), with the DFG
2D-FAC boundary types (`Inflow2` at x=0, `Symmetry011` at x=Lx, `Wall`
at y=0/Ly, `Symmetry001` at z=0/Lz) and a project file for
`tools/PyPartitioner.py`.

```bash
python3 channel_tri.py --lx 2.2 --ly 0.41 --lz 0.05 --nx 44 --ny 8 --nz 1 \
        --name chimera_fac --outdir ../../applications/q2p1_chimera/_data/CHIMERA_FAC
```

This call produced the vendored milestone-1 case
`applications/q2p1_chimera/_data/CHIMERA_FAC` (352 hexes, h = 0.05 at the
coarse level; 0.025 at `MaxMeshLevel = 2`, 0.0125 at level 3).  The
atmosphere for that case is the Phase-2 annulus fixture, fitted at load
time to r_i = 0.05, r_o = 0.1 (`CHI_SUBMESH_FIT_COARSE`).
