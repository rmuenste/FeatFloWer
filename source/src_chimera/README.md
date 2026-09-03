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

## Current state (Phases 0–3)

| File | Phase | Purpose |
|---|---|---|
| `chimera_config.f90` | 0 | Runtime keys + `CHIMERA_VALIDATE_CONFIG` |
| `chimera_api.f90` | 0 | Facade; `Chimera_IsEnabled`, lifecycle contract, fatal stubs for later phases |
| `chi_geometry.f90` | 1 | Q1 map (pinned vs `EL_Q1_MAP`), Newton inverse map, 3×3×3 Gauss, 3×3 inverse |
| `chi_fem_eval.f90` | 1 | Q2 basis in FeatFloWer local ordering, DOF map, physical gradients, P1 eval |
| `chi_locator.f90` | 1 | Instance-based element-bbox bucket-grid point locator |
| `chi_sparse_direct.f90` | 1 | Instance-based UMFPACK wrapper (per-instance handles, owned CSR copies) |
| `chi_kernels.f90` | 2 | Reentrant Q2/P1 saddle-point assembly (deformation-form NS, B/Bᵀ, Robin, Dirichlet rows, face quadrature); the legacy F77 kernels (`/ELEM/ /CUB/ /TRIAD/`) are never used for submeshes |
| `chi_submesh.f90` | 2 | `tChimeraSubmesh`: own `tMultiMesh`, bitmask boundary classification with hierarchical propagation, radial projection, Q2 coords, boundary face lists |
| `chi_solver.f90` | 2 | Per-submesh monolithic Picard/direct solve (Dirichlet or Robin outer BC, pressure gauge) |
| `chi_forces.f90` | 2 | Surface-stress force/torque with the design §5 normal/sign conventions |
| `chi_legacy_mesh_adapter.f90` | 2 (Layer H) | Sole bridge to `mesh_structures` (readTriCoarse/refineMeshLevel/genMeshStructures), with coarse-shell fitting and per-level classify+project |
| `chi_markers.f90` | 3 | Hole/fringe classification kernel (paper §6 definition; two-array `kind`/`pid` scheme with the MAX merge rule) |
| `chi_exchange.f90` | 3 | MPI service on a passed communicator: replicated query lists, MIN-rank ownership, SUM reduction; broadcast |
| `chi_output.f90` | 3 | Legacy-VTK dump of a submesh solution |
| `chi_penalty.f90` | 4 | Weak-variant penalty operator: paper's damping function beta (parametrised ramp), 27-point tabulation + consistent Q2 penalty matrix D and vector g on the level pattern, nodal (Lobatto) lumped penalty D_L (positive, constant-exact), Jacobi-PCG for paper eq. (12) with sum/filter/allsum callbacks |
| `chi_coupling.f90` | 3-4 (Layer H) | The SAVEd coupling state and hook bodies: body table, replicated submeshes + solvers, background locator, markers synchronised with `E013Max_SUPER`, donor caches, per-step Robin exchange → owner solve → broadcast → fringe values (strong) or nodal Dirichlet data g (weak) → forces; weak variant: per-level D_L (E013Sum'd copies), H8/H10/H11/H15 bodies, coupling under-relaxation; restart write/read of the replicated submesh states |
| `tests/` | 1–4 | ctests; `chi-submesh-couette` is the Phase-2 analytic gate; `chi-markers-cutcell` and `chi-exchange-np{1,2,3}` are the Phase-3 tests; `chi-algebra-serial` (symmetry, D·1 = g(const), penalised volume, PSD, nodal lumping positivity/volume, PCG residual + filter rows, D = 0 bitwise) is the Phase-4 test |

Implementation notes recorded in Phase 2: the FEAT 1:8 refinement yields
child elements of mixed orientation, so boundary normals are oriented
geometrically (`CHI_FACE_GEOM` centroid test); with Q1 geometry the
curved-boundary L2 convergence limit is 2nd order (production shares
this); direct surface-traction torque converges sub-quadratically.

Implementation notes recorded in Phase 3: the submesh problem is the
time-discrete one (backward Euler with the background Δt from its own
previous level) — a steady submesh solve coupled to the impulsively
started background blows up on the initial pressure transient; the
matrix-row filter follows the `FictKNPR` practice of indexing the
finest-level marker array on every multigrid level (coarse levels only
affect the preconditioner); cylinder z-faces of an atmosphere carry the
slab symmetry condition (w = 0 only); fringe nodes in the chord gap of
the polygonal inner surface are located by a relaxed nearest-element
search (`CHI_LOCATE_NEAREST`, extrapolation reported at init).

Implementation notes recorded in Phase 4 (weak variant, Chimera-W):
the velocity multigrid's coarse solver is a damped Jacobi iteration on
the stored global diagonal, which diverges for the consistent
(mass-like) Q2 penalty matrix once `dt*gamma` is large — the production
penalty is therefore the nodal Lobatto-lumped `D_L` (row-sum lumping of
a beta-weighted Q2 mass matrix is not positive); the correction of paper
eq. (12) and the pressure Poisson operator use the same (optionally
capped) penalised lumped mass, default cap 0 = plain projection
(momentum-only penalty; a penalised Poisson operator steepens the
interior pressure and diverges at the ramp edge on coarse backgrounds;
the paper's damped correction with the plain operator drifts); the
step-iterated Robin ↔ penalty loop is stable only if the free band
between the last penalised node and the atmosphere boundary is at least
about one background cell (keys `ChimeraBetaFull/Zero`; the paper's
0.5/0.75 assume h ≪ H). The consistent matrix + CG path is kept as the
diagnostic `ChimeraPenaltyLumped = No`.

## Hooks in existing code (complete list as of Phase 4)

- `source/src_quadLS/QuadSc_main.f90` — hook H1 in
  `Transport_q2p1_UxyzP_fluid_core` (guarded `Chimera_BeginStep`).
- `source/src_quadLS/QuadSc_boundary.f90` — hooks H2/H3/H4 in
  `Boundary_QuadScalar_Def/_Val/_Mat/_Mat_9` (guarded facade calls,
  siblings of the `FictKNPR` branches).
- `source/src_quadLS/QuadSc_main.f90` — hook H10 beside the body-force
  terms (`Chimera_AddMomentumRHS`, weak variant only).
- `source/src_quadLS/QuadSc_def.f90` — hook H8 at the end of the per-level
  matrix assembly in `Matdef_general_QuadScalar` (`Chimera_AddMomentumMatrix`)
  and its defect sibling in the part-wise (non-Newtonian, S = 0) defect
  branch (`Chimera_AddMomentumDefect`); hook H15 in `Create_ParCMat`
  (`Chimera_AddPressureMass`).
- `source/src_quadLS/QuadSc_assembly.f90` — hook H15 in `Create_CMat`
  (penalised lumped mass in `B^T M^-1 B`).
- `source/src_quadLS/QuadSc_corrections.f90` — hook H11 in
  `Velocity_Correction` (`Chimera_CorrectVelocity`; the existing diagonal
  loop runs verbatim when the call reports `applied = .FALSE.`) plus the
  helper `Chimera_DefectFilter3` (the module's Dirichlet defect filter on
  caller-owned vectors, used by the CG of the consistent-penalty mode).
- `applications/q2p1_chimera/` — hooks H5/H6 (app-local
  `Chimera_Initialize`/`Chimera_Finalize`), H12 (restart write beside every
  flow dump, read after initialisation), and the weak-variant rebuild of
  the projection operators (`Create_CMat`/`Create_ParCMat`) after the
  penalty exists; drives the shared fluid core with `enable_fbm = .FALSE.`.
- `source/src_util/param_parser.f90` — `CASE ("Chimera*")` branches,
  validation call, enabled-only echo.
- `cmake/modules/ProjectFiles.cmake`, `GenerateLinkerFlags.cmake` —
  library + test wiring.
