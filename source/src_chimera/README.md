# src_chimera — Chimera overlapping-mesh component

Implementation of the Chimera domain-decomposition method of
arXiv:2506.22831 (`literature/2506.22831.pdf`): background Q2/P1 mesh +
per-particle body-fitted "atmosphere" submeshes with Dirichlet–Robin
coupling. Authoritative design document: `chimera-integration-design.md`
(repo root, v3); runtime keys: `docs/md_docs/parameter_reference.md`
("Chimera Overlapping-Mesh Component").

## Ground rules (from the design, enforced in review)

- No `#ifdef`s anywhere; the component is always compiled and enabled only
  via `SimPar@ChimeraEnable` (default `No`). A disabled run must be
  bit-identical to a build without the component (off-regression gate:
  `q2p1_fc_ext_cylinder` at tolerance 0, protocol byte-identity included).
- Existing solver files interact with the component through the single
  facade `chimera_api.f90` (`MODULE CHIMERA_API`) only; `param_parser`'s
  direct import of `CHIMERA_CONFIG` is the sole exception.
- Layer discipline:
  - **Layer L** (`chimera_config.f90`, built into `ff_util`): runtime
    configuration + validation; no other dependencies.
  - **Layer M** (`chi_*.f90`, built into `ff_chimera`): COMMON-free and
    `var_QuadScalar`-free **transitively**; instance-based state only, no
    module singletons; never calls `SETLEV` or legacy F77 assembly
    kernels.
  - **Layer H** (`chimera_api.f90`, later `chi_coupling.f90`/
    `chi_penalty.f90`/`chi_legacy_mesh_adapter.f90`, built into
    `ff_quadLS_app`): may USE `var_QuadScalar`/`def_FEAT`, owns every
    call into legacy global-state code.
- No print-and-continue stubs: an enabled-but-unimplemented or
  enabled-but-uninitialized path aborts with a clear message (`STOP 1`).

## Current state (Phase 0–1)

| File | Phase | Purpose |
|---|---|---|
| `chimera_config.f90` | 0 | Runtime keys + `CHIMERA_VALIDATE_CONFIG` |
| `chimera_api.f90` | 0 | Facade; `Chimera_IsEnabled`, lifecycle contract, fatal stubs for later phases |
| `chi_geometry.f90` | 1 | Q1 map (pinned vs `EL_Q1_MAP`), Newton inverse map, 3×3×3 Gauss, 3×3 inverse |
| `chi_fem_eval.f90` | 1 | Q2 basis in FeatFloWer local ordering, DOF map, physical gradients, P1 eval |
| `chi_locator.f90` | 1 | Instance-based element-bbox bucket-grid point locator |
| `chi_sparse_direct.f90` | 1 | Instance-based UMFPACK wrapper (per-instance handles, owned CSR copies) |
| `tests/test_chi_*.f90` | 1 | Serial ctests (`chi-*-serial`) |

Phase 2+ (per the design roadmap): `chi_kernels.f90` (new reentrant
assembly kernels — the legacy F77 kernels read `/ELEM/ /CUB/ /COAUX1/
/TRIAD/` and must not be called for submeshes), `chi_submesh/solver/
forces`, `chi_legacy_mesh_adapter`, `chi_exchange`, `chi_coupling`,
`chi_penalty`, application `q2p1_chimera`, tool `tools/chimera_meshgen`.

## Hooks in existing code (complete list as of Phase 0)

- `source/src_quadLS/QuadSc_main.f90` — hook H1 in
  `Transport_q2p1_UxyzP_fluid_core` (guarded `Chimera_BeginStep`).
- `source/src_util/param_parser.f90` — `CASE ("Chimera*")` branches,
  validation call, enabled-only echo.
- `cmake/modules/ProjectFiles.cmake`, `GenerateLinkerFlags.cmake` —
  library + test wiring.
