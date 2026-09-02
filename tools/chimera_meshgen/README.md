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
