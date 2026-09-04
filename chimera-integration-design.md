# Chimera component — software design & integration plan (v3)

Drafted 2026-09-01 (v1); revised 2026-09-02 (v2, findings 1–8 of
`chimera-plan-review.md`); **revised 2026-09-02 to v3 in response to
`chimera-plan-review-v2.md`** (both remaining blockers and all smaller
items adopted; disposition in §12/§13). Design document for integrating
the Chimera overlapping-mesh method (arXiv:2506.22831,
`literature/2506.22831.pdf`) into FeatFloWer. Per the v2 verdict,
**Phases 0–1 are cleared for implementation; Phase 2 is gated on this
revision.**

---

## 1. Context and goals

The paper introduces a multimesh method for particulate-flow DNS: the fixed
background Q2/P1 mesh solves Navier–Stokes over the whole fictitious
domain, and each particle carries a small body-fitted annular/shell
submesh ("atmosphere") solving its own Navier–Stokes problem with a
rigid-body Dirichlet BC on the inner boundary and a Robin BC (data
interpolated from the background solution) on the outer boundary.
Forces/torques come from surface-stress integration on the submesh.

Two coupling variants (paper Section 6):

- **Chimera-S (strong)** — hole and fringe nodes of the background mesh
  get strongly imposed Dirichlet values; the paper's per-time-step
  algorithm is background solve → submesh solves → background re-solve,
  **at least two outer iterations within one time step**.
- **Chimera-W (weak)** — a distributed interior-penalty matrix `D` plus
  vector `g` enter the background momentum equation and the
  velocity-correction step; the paper's algorithm **may terminate after
  one outer iteration** (submesh solves with old background data →
  background solve with the penalty).

**Coupling-fidelity ground rule:** *static particle geometry does not
imply steady fluid physics.* A single coupling pass per time step
realizes the paper's algorithms only for (a) genuinely steady /
pseudo-time-converged problems (successive time steps act as fixed-point
coupling iterations) and (b) time-dependent problems under **Chimera-W**
(its one-iteration variant is the paper's own algorithm). Time-accurate
**Chimera-S** requires the in-step outer driver (H13). Unsteady cases
(DFG 2D-2 vortex shedding) are excluded from milestone 1.

Strategic context: this component is the next-generation force/torque
engine for the DNS-ground-truth role in the FF-EL "DNS-informed EL"
closure program. **Milestone 1 is static, steady Chimera-S**: the steady
3D flow-around-cylinder (FAC) configuration already regression-pinned in
this repo (`q2p1_fc_ext_cylinder`, C_D 5.5795 / C_L 0.010619), then
Hasimoto arrays and random-array drag closures.

### Hard design requirements (project lead)

1. The component plugs in easily.
2. The standard operational mode is NOT disturbed — machine-checkable.
3. Encapsulation wherever possible.
4. Chimera functionality lives in separate files/compilation units.

---

## 2. Governing architectural decisions

**All Chimera code lives in `source/src_chimera/` as proper Fortran
modules** (`chi_`/`Chimera_` namespace, `IMPLICIT NONE`, `PRIVATE` with
explicit `PUBLIC` lists). No textual `include`s into `MODULE
Transport_Q2P1`.

**No compile-time option, no `#ifdef`s** *(project-lead decision)*: the
component is always compiled (UMFPACK is already mandatory) and enabled
purely by the runtime key `SimPar@ChimeraEnable` (default `No`). The
shipped default binary must reproduce existing baselines bit-identically
(§9).

**Single facade, internally rank-safe (v2-review smaller items 1):**
existing solver files interact with Chimera through exactly one module,
**`CHIMERA_API`**; all state stays `PRIVATE` behind it. **Every facade
operation is rank-safe by itself** — operations that must not run on the
master (or on any rank) return immediately after an internal check;
callers never need scattered `myid` guards. The one necessary exception
is documented at H1: the *call site* sits inside an existing
`myid.ne.master` region because argument association
(`LinSc%valP(NLMAX)%x`) is itself only safe on workers. `Chimera_BeginStep`
on an enabled-but-uninitialized state **aborts with a clear message**
("ChimeraEnable=Yes requires an application that initializes Chimera,
e.g. q2p1_chimera") — a hard `STOP 1`, never print-and-continue.

```fortran
Chimera_IsEnabled() / Chimera_VariantIsWeak()
Chimera_Initialize(mfile)                             ! app-local (H5)
Chimera_BeginStep(valU, valV, valW, valP)             ! H1
Chimera_ApplyBoundaryDef(defU, defV, defW, ndof)      ! H2
Chimera_ApplyBoundaryValues(valU, valV, valW, ndof)   ! H3
Chimera_FilterMatrixRows(A11, A22, A33, KLD, ndof)    ! H4 (+ 9-block form)
Chimera_AddMomentumMatrix(A11, A22, A33, qMat, ilev)  ! H8 (W)
Chimera_AddMomentumRHS(defU, defV, defW)              ! H10 (W)
Chimera_CorrectVelocity(applied)                      ! H11 (W; see §5 fast path)
Chimera_WriteRestart(unit) / Chimera_ReadRestart(unit)! H12
Chimera_Finalize()                                    ! H6, app-local
```

The sole facade bypass is `param_parser`, which imports `CHIMERA_CONFIG`
directly (the established `bUseHashGridAccel` pattern).

**Lifecycle pairing (v2-review smaller items 2):** initialization and
finalization are BOTH application-local — `q2p1_chimera` calls
`Chimera_Initialize` after `init_q2p1_app` and `Chimera_Finalize` after
the time loop. No hook in shared `Init_QuadScalar_Stuctures` at all
(one fewer shared-file touch point than v2). Any other application that
sets `ChimeraEnable=Yes` without initializing hits the H1 fatal check.
`Chimera_Finalize` is idempotent and safe on partial initialization.

**Zero new fields in `var_QuadScalar`, zero COMMON usage in new code.**

### Precedents imitated / avoided

| Imitate | Where |
|---|---|
| Config module + validation (`STOP 1`) | `source/src_el/el_config.f90` |
| SAVEd state + lifecycle + versioned restart | `source/src_el/el_fields.f90` |
| Runtime flag / `LOGICAL` argument, never `#ifdef` | `fluid_core(..., enable_fbm)` |
| Additive physics term behind a runtime flag | EL drag block, `QuadSc_def.f90:2417/2516` |
| Central parameter parsing into low-layer config | `param_parser.f90` |
| Unit tests beside the subsystem, serial + MPI | `source/src_el/tests/` |

| Avoid | Where it went wrong |
|---|---|
| `#ifdef` leakage into core files / app_inits | PE integration |
| App-local parser copies | `q2p1_el_pipeflow/app_init.f90` |
| Undefined preprocessor symbols | `QuadSc_main.f90:1079` |
| `include 'X.f90'` pseudo-modularity | `QuadSc_main.f90:71-79` |
| Print-and-continue stubs | `fbm_main.f90:187` |
| State in `var_QuadScalar` / COMMON | `QuadSc_var.f90` |
| Module-global solver handles | `source/UmfpackSolver.f90:1-73` |
| **Legacy F77 assembly kernels for submesh work** | **they read `/ELEM/ /CUB/ /COAUX1/ /TRIAD/` and `var_QuadScalar` (v2-review blocker 1)** |

---

## 3. Module decomposition (three dependency layers)

### Layer L — appended to the `src_util` source list

**`chimera_config.f90` — `MODULE CHIMERA_CONFIG`** — runtime keys with
defaults and `CHIMERA_VALIDATE_CONFIG()` (`STOP 1` on violations; the
FBM-exclusivity cross-check runs in `Chimera_Initialize`, where runtime
context exists):

```fortran
LOGICAL           :: chimera_enable = .FALSE.
CHARACTER(LEN=8)  :: chimera_variant = 'strong'   ! 'strong' | 'weak'
CHARACTER(LEN=16) :: chimera_outer_bc = 'robin'   ! 'robin' | 'dirichlet' (diagnostic)
CHARACTER(LEN=256):: chimera_particle_file = '', chimera_submesh_file = ''
INTEGER :: chimera_submesh_nlmax = 3
REAL*8  :: chimera_robin_alpha = 1d0, chimera_gamma_max = 0d0
INTEGER :: chimera_outer_iters = 1
INTEGER :: chimera_sub_nl = 3
LOGICAL :: chimera_write_vtk = .FALSE.
LOGICAL :: bChimeraS = .FALSE., bChimeraW = .FALSE.
```

### Layer M — new library `ff_chimera`
(`add_library(ff_chimera ...)` linking `ff_util ${FF_DEFAULT_LIBS}`;
appended to `FF_APPLICATION_LIBS`)

COMMON-free and `var_QuadScalar`-free **verified transitively** — the two
v2-review counterexamples are resolved by construction: no legacy
assembly kernels are called at all (below), and legacy mesh I/O lives in
a Layer-H adapter.

**`chi_geometry.f90`** — dependency-free geometric kernel: Q1 trilinear
forward map + Jacobian (functional twin of `EL_Q1_MAP` — pinned against
it by unit test; `el_quadrature.f90` itself is NOT compiled into
ff_chimera, it USEs `EL_HALO`/`EL_FIELDS`/`PP3D_MPI` and already belongs
to `ff_quadLS_app`), Newton inverse map (damped, `GetPointFromElement`
idiom, `part_step.f90:26`), 3×3×3 Gauss cubature, 3×3 inverse.

**`chi_fem_eval.f90`** — Q2 basis (FeatFloWer 27-node local ordering,
verified against the hard-coded basis in `RETURN_Velo`,
`part_step.f90:1711`): values + reference derivatives at arbitrary xi,
physical gradients via the Q1 Jacobian, local→global Q2 DOF map
(vertices/edges/faces/center), P1 centroid-linear pressure evaluation.

**`chi_locator.f90`** — instance-based element-bbox bucket grid,
`BUILD / LOCATE / RELEASE`. (The `OctTreeSearch.f90` singleton is not
touched.)

**`chi_sparse_direct.f90`** — instance-based direct solver
`TYPE tSparseDirectSolver` (v2-review blocker-5 lineage): owns its
0-based CSR copy (caller arrays never mutated), per-instance
`symbolic`/`numeric` handles — safe because the umf4 F77 wrapper stores
pointers in a handle table
(`extern/libraries/umfpack4/src/umf4_f77wrapper_port.c`, `StorePointer`),
so concurrent factorizations are supported by design. `INIT / FACTORIZE
/ SOLVE / FREE`; `FREE` idempotent. The existing `UMFPackSolver` module
keeps serving the background coarse solver untouched.

**`chi_kernels.f90`** *(Phase 2; replaces v2's `chi_assembly` F77-reuse
plan — v2-review blocker 1)* — **new reentrant Chimera assembly
kernels**, written on `chi_geometry`/`chi_fem_eval`: element loops with
`nel`, connectivity, coordinates, cubature, basis data, physical
parameters (ρ, μ, Δt, θ), mesh velocity, and target CSR all passed
explicitly. Implementation finding (Phase 2): the FEAT 1:8 refinement
produces child elements of MIXED orientation (negative detJ — harmless
in production, which uses |detJ| and parametrized normals), so boundary
normals are oriented geometrically (centroid→face test in
`CHI_FACE_GEOM`), never from reference sign tables alone. Terms: mass, diffusion, ALE convection (`u − u_mesh` with
constant per-submesh mesh velocity — one subtraction; nothing is lost by
dropping `CONVQ2`), B/Bᵀ, Robin surface terms (§5), inner Dirichlet
rows, pressure gauge (one pinned pressure DOF) for the all-Dirichlet
diagnostic mode. The legacy kernels (`QuadSc_laplace.f`, `QuadSc_conv.f`,
`QuadSc_massrho.f`, `QuadSc_BMatrix.f`) are **never called for
submeshes**: they read `/ELEM/ /CUB/ /COAUX1/ /TRIAD/` and
`var_QuadScalar` configuration, so "passing mesh arrays" is insufficient
and COMMON-switching would break reentrancy. Verified by the Phase-2
annular-Couette analytic gate.

**`chi_submesh.f90`** *(Phase 2)* — `TYPE tChimeraSubmesh`: own
`tMultiMesh`, placement, rigid-body state, inner/outer boundary DOF+face
lists **with stored outward normals** (§5 conventions), donor caches,
a `tSparseDirectSolver`. Everything downstream of the raw coarse read
(which arrives via the Layer-H adapter).

**`chi_solver.f90`** *(Phase 2)* — per-submesh monolithic saddle-point
solve (`3·ndofQ2 + 4·nel`): symbolic once, Picard loop with numeric
refactorization.

**`chi_forces.f90`** *(Phase 2)* — force/torque per §5;
`ChimeraForce:` protocol lines (only when enabled).

**`chi_exchange.f90`** *(Phase 3)* — MPI service on **`MPI_COMM_SUBS`**
(workers only; the master returns immediately): identical replicated
query lists on all workers; owner selection by `MPI_Allreduce(MPI_MIN)`
over proposed ranks; value reduction by `MPI_Allreduce(MPI_SUM)` over
owner-masked contributions; equal counts by construction.

### Layer H — appended to the `ff_quadLS_app` source lists

**`chi_legacy_mesh_adapter.f90`** *(Phase 2)* — owns every call into
legacy mesh I/O (`readTriCoarse`/`refineMesh` `USE var_QuadScalar`,
`mesh_refine.f90:210-218,1205-1213`); a temporary bridge by contract,
never an assembly path.

**`chi_coupling.f90`** *(Phase 3)* — the private SAVEd state and hook
implementations. **Marker representation (v2-review blocker 2):** two
arrays, never a signed encoding —

```fortran
INTEGER*1, ALLOCATABLE :: marker_kind(:)   ! 0=free, 1=fringe, 2=hole
INTEGER,   ALLOCATABLE :: marker_pid(:)    ! particle id where kind>0, else 0
```

`E013Max_SUPER`-style numeric MAX on `marker_kind` implements exactly the
required precedence (hole > fringe > free); a second pass synchronizes
`marker_pid` among ranks that agree on the winning kind (MAX over pid
where the local kind equals the reduced kind, 0 elsewhere — unambiguous
under M1's non-overlapping atmospheres). The v2 signed scheme (`−k`
fringe) would have been destroyed by MAX(−k, 0)=0 on partition
interfaces. **Normative fringe/hole definition (paper Section 6):** for
a background cell crossed by `∂B_k`, a node `x_i` is a *hole node* if
`x_i ∈ B_k` and a *fringe node* if `x_i ∈ Ω̂_k`; in the paper's
Chimera-S sequence, the first outer iteration imposes hole nodes only,
subsequent updates hole+fringe. M1's step-iterated variant imposes
hole+fringe with previous-step submesh values, whose fixed point is the
converged paper iteration. A dedicated MPI test covers a cut cell
spanning two partitions (§9).

**Phase-3 implementation notes** (deviations/refinements of the text
above, all reviewer-visible):

- *Time-discrete submesh problem.* The atmosphere step is the same
  backward-Euler step as the background (`ρ/Δt·M` on the matrix, the
  atmosphere's own previous level on the right-hand side,
  `CHI_ASM_MASS_RHS`). A steady submesh solve coupled to the impulsively
  started background blows up on the initial pressure transient
  (observed: O(10²) atmosphere velocities, NaN in the background after
  six steps); the time-discrete form is stable and its steady limit is
  the intended coupled steady state. `chi_solver` keeps the steady mode
  (`dtinv = 0`) for the analytic tests.
- *Marker arrays* are default `INTEGER` (not `INTEGER*1`) so that the
  production `E013Max_SUPER` synchronises them without copies; the
  two-array semantics are unchanged. The second (id) pass is exactly
  the conditional MAX described above.
- *Matrix-row filter on coarse levels* follows the `FictKNPR` practice
  (finest-level marker array indexed with coarse-level dof numbers);
  coarse levels only shape the preconditioner, the converged solution
  is unaffected. A per-level marker hierarchy is a possible later
  refinement.
- *Atmosphere z-faces* (cylinder in the slab) carry the slab symmetry
  condition `w = 0` only; u, v are free (`CHI_DIR_*` component masks in
  the tabulated solve). Cylinders are treated as infinite along z.
- *Coarse shell fitting.* One generic shell `.tri` serves all bodies:
  `CHI_SUBMESH_FIT_COARSE` remaps the shell's own radial (and axial)
  range affinely onto `[r, r+H]` (and `[zlo, zhi]`) before
  classification, so the Phase-2 fixture is the milestone-1 atmosphere.
- *Fringe nodes in the chord gap* between the analytic inner surface
  and its polygonal approximation are located by a relaxed nearest-
  element search (`CHI_LOCATE_NEAREST`), i.e. a slight Q2 extrapolation;
  the maximum reference-coordinate excess is reported at init (0 on the
  milestone-1 case).
- *Body table* format and the application layout are documented in
  `docs/md_docs/parameter_reference.md` / `chimera_usage.md`.
- *Milestone-1 steady FAC results* (DFG 2D-1, uniform channel
  background h = 0.05/2^(L-1), annulus atmosphere r_o = 2 r, Robin
  α = 1, backward Euler; forces from the atmosphere surface stress):

  | Case | C_D | C_L | Note |
  |---|---|---|---|
  | Chimera L2 / atmosphere L2 (h = D/4), dt 0.01, t = 4.56 | 5.515 | −0.040 | lift sign wrong: hole/fringe staircase asymmetric by up to h = D/4 top vs bottom |
  | Chimera L3 / atmosphere L3 (h = D/8), dt 0.05, t = 4.0 | 5.565 | +0.013 (decaying) | drag −0.3 % vs body-fitted L3; lift error shrinks 3–4× per refinement |
  | body-fitted `q2p1_fc_ext` L2 / L3 (t = 10) | 5.6013 / 5.5808 | 0.00995 / 0.01067 | measured baselines |
  | DFG reference | 5.5795 | 0.01062 | Schäfer–Turek band C_D 5.57–5.59, C_L 0.0104–0.0110 |

  Conclusion: the coupled Chimera-S solution converges to the DFG
  values under joint refinement; on the coarse L2 background (four
  cells per diameter) the lift is dominated by the staircase
  representation of the body in the background (as expected for the
  strong variant — paper §6 uses h ≪ D). The pinned regression value
  (`q2p1_chimera_cylinder`, L2, dt 0.05, t = 6) is therefore a
  *regression anchor*, not an accuracy claim; the accuracy statement is
  the L3 row. Atmosphere-only refinement (background fixed) does not
  fix the lift — the limiting error is in the background near the
  fringe, which is the motivation for the weak variant / finer
  backgrounds in Phase 4/5.
- *Milestone-1 restrictions enforced at init:* `strong` only; no FBM
  particles (`myFBM%nParticles > 0` aborts — the `calculateFBM()` switch
  is on even in the plain FAC run and is therefore not a usable
  criterion).

**`chi_penalty.f90`** *(Phase 4, DONE 2026-09-03)* — Layer-M penalty
kernels: the paper's damping function β with a parametrised ramp
(`ChimeraBetaFull/Zero`), 27-point tabulation + consistent Q2 penalty
matrix `D` and vector `g` on the level pattern, the **nodal
(Lobatto-quadrature) lumped penalty `D_L`** (positive; constant-exact
with `g_i = D_L(i)·û(x_i)`) and a Jacobi-PCG for eq. (12) with
sum/filter/allsum callbacks. Layer H (`chi_coupling`) keeps per-level
`D_L` (rank-partial and E013Sum'd copies; static: built once) and the
finest-level nodal donor list; `g` is rebuilt every coupling update.

*Phase-4 implementation notes (what the FAC case taught):*
- The velocity multigrid's coarse solver (`Velo@MGCrsSolverType = 1`) is
  a damped Jacobi iteration on the stored global diagonal
  (`MGE013%UE11`); the consistent penalty matrix (mass-like, not
  diagonally dominant) makes it diverge once `dt·γ ≳ 100`. The
  production operator is therefore `D_L` (row-sum lumping of a
  β-weighted Q2 mass matrix is *not* positive; nodal quadrature is).
  The consistent matrix + CG correction remains as the diagnostic
  `ChimeraPenaltyLumped = No`.
- Projection: the paper's damped correction (12) with the plain Poisson
  operator (11) leaves the divergence in the penalised zone uncorrected
  and was observed to drift; using the penalised lumped mass in both
  (11) and (12) — the variable-density mechanism, `ChimeraProjCap > 0` —
  steepens the interior pressure gradient by `1+κ` and blows up at the
  ramp edge on `h = D/4`. Default: **plain projection** (`ProjCap = 0`,
  momentum-only penalty, standard correction loop = fast path); the
  corrected velocity is exactly discretely divergence-free and the fixed
  point is the penalised steady state with an O(1/γ) interior leak.
- The step-iterated Robin ↔ penalty loop has a slowly growing mode on
  the L2 background (growth ≈ ×1.2 per step, onset t ≈ 4, for γ = 1e3,
  1e4, 1e5 alike) when the paper's β support (up to `R + 0.75H`) leaves
  less than a background cell free before the atmosphere boundary Γ
  where the Robin data are sampled; with the support shrunk to
  `R + 0.5H` (`ChimeraBetaFull/Zero = 0.25/0.5`, one free cell) the run
  converges monotonically. `ChimeraCouplingRelax` (under-relaxation of
  the Dirichlet data, both variants) is available as an additional
  damping knob. The strong variant is stable on the same case.
- The weak variant freezes `dt` in the projection operator (constant
  `TimeStep` enforced at the correction).

**`chi_periodic.f90`** *(Phase 5, Layer M)* — periodic-box geometry as a
value type `tChiPeriodic` (per-axis flag, period, box origin):
minimum-image displacement/distance and box wrap, exact identities when
inactive. The box is derived in `chi_coupling` from the base code's
`dPeriodicity` (now settable from the shared parser with
`SimPar@PeriodicLength = Lx,Ly,Lz`; the E013 communicator pairs the
opposite-face dofs, there is no periodic BC in the discretisation) plus
the global background bounding box (allreduce; the background must span
exactly one period per periodic axis, checked). It is passed as an
OPTIONAL argument into the Layer-M kernels (`CHI_POINT_IN_BODY`,
`CHI_CLASSIFY_MARKERS`, `CHI_PENALTY_BETA/TABULATE/NODAL`), so every
"x − X_k" is the minimum image; donor searches map background points into
the body frame (`image_point`), Robin/Dirichlet sample points of an
atmosphere straddling a periodic face are wrapped into the box before
`CHI_EXCHANGE_BG_EVAL`. `check_atmospheres` enforces the paper's
non-overlap assumption (atmosphere ∩ other body = ∅ for both variants,
pairwise-disjoint atmospheres for Chimera-S, own periodic images
included). Element prefilters (`CHI_BODY_NEAR_BOX`) make the marker
classification O(nel) per body for arrays.

**Phase-5 implementation notes (DONE 2026-09-04, milestone: Hasimoto /
random arrays):**

- Tooling (`tools/chimera_meshgen/`): `sphere_shell_tri.py` (equiangular
  cubed-sphere shell, one file fitted per body), `seed_array.py` (SC
  lattice / RSA in a periodic box, emits the body table with
  `H_k = min(H_max, ½ nearest gap)`, minimum-image), `channel_tri.py
  --periodic`, `partition_periodic_box.py` + `flatten_axis_partition.py`
  (axis-aligned 2×2×2 partition in the single-host `sub0001/GRID000k`
  layout that the host-name rank mapping of `get_pid.f90` expects).
- `ChimeraSubStokes`: linear (Stokes) submesh operator factorised once,
  rhs-only coupling updates (`chi_solver::solve_frozen`). On one core the
  Picard path costs ≈ 2 min/step (two UMFPACK factorisations of 47670
  unknowns); the frozen path ≈ 1 s/step, identical forces to 8 digits at
  Re 3e-3.
- The submesh momentum equation receives the background's uniform body
  forces (`ConstantForcing`, fluid gravity) as `rho·M·f`; without it the
  steady force of a periodic array misses exactly the atmosphere volume
  fraction (observed −14 % before the fix).
- Composite bulk velocity (`ChimeraBulk:` lines, `bulk_diagnostic`): cell
  average with 0 in the bodies, the submesh field in the atmospheres,
  the background elsewhere — the superficial velocity of the drag
  closures; plus the summed surface forces for the momentum balance.
- Hasimoto convention (RUNBOOK `m2_hasimoto`): the surface traction of a
  fluid-only body force is `f·V_fluid` at steady state; Hasimoto's
  mean-gradient force is `F_H = F_s + f·V_solid` (the FBM constraint force
  contains this term because the hole fluid is forced).
- Results at D/h = 4 (6³ box, L2, atmosphere level 2): strong
  K_meas 1.8020 (−1.6 % vs 1.8317), K_bal 1.8165 (−0.8 %), balance
  residual −0.82 %; weak (γ 1e5, Phase-4 settings) 1.8008 / 1.8157;
  corner-placed sphere (atmosphere across all six periodic faces) equal
  to the centred one to 5e-8. FBM (D1.1) needs D/h = 24 for −2.3 %.
  D/h = 8 (L3 background, level-2 atmosphere): strong 1.8034 / 1.8151,
  weak 1.8024 / 1.8140, balance residual −0.66 % (halved) — the
  remaining ≈ 1 % is the atmosphere discretisation; the level-3 shell
  (362598 unknowns) exceeds the direct solver (UMFPACK failure; Vanka/MG
  or 64-bit UMFPACK is the upgrade path).
- Random array, 8 spheres, φ = 0.05, L3 background: strong mean K
  2.356 / K_bal 2.364, weak 2.354 / 2.362, per-sphere forces equal to
  4 digits (spread ±8 % = the FF-EL dispersion datum). The value sits
  between the SC lattice (2.51) and the random-array closures
  (Beetstra 1.76, Hill–Koch–Ladd 1.95, Tenneti 1.53): the strong
  variant's fringe rule forces `gap ≳ 4h` (0.73 d here), so the seeded
  configuration is a perturbed lattice, not hard-sphere microstructure —
  the Beetstra–Tenneti band is a resolution question (h ≲ 0.012 d for
  0.05 d gaps) for Chimera-S. The Chimera-W probe on a configuration
  with 0.22 d gaps and atmospheres down to 0.6 background cells runs
  stably (balance +0.2 %) and gives mean K 2.003 — 2.5 % above
  Hill–Koch–Ladd, 14 % above Beetstra, i.e. the upper edge of the band
  for one realisation: random arrays with closure-fit microstructure
  are a Chimera-W job on coarse backgrounds.
- Off-regression: FAC strong/weak anchors bit-identical after the
  Phase-5 changes; `q2p1_fc_ext` protocol unchanged (parser key inert
  when absent).
- Deferred: the halo (neighbourhood) exchange upgrade — the replicated
  collective is adequate for the single-host runs; the interface
  (`CHI_EXCHANGE_BG_EVAL`) is unchanged.

**`chimera_api.f90`** *(Phase 0 onward)* — the facade; thin delegation.

---

## 4. Parallel strategy

Submeshes worker-replicated (all workers, M1); solve on the owner
worker; solution broadcast over `MPI_COMM_SUBS`; deterministic
replicated evaluation ⇒ bit-identical interface values without
arbitration. No second set of E013/MGE013 comm structures, ever. The
master holds no Chimera runtime state; every facade operation no-ops
there internally. Scaling upgrades (el_halo-style exchange, intersecting
-worker replication) slot in behind unchanged interfaces (Phase 5).

---

## 5. Mathematical conventions (normative)

**Discrete-time treatment of the penalty** (fully implicit `D`, paper
eqs. (9)–(12); FeatFloWer rows carry an overall factor Δt relative to
the paper's `M/Δt` scaling):

- H8: `A11/A22/A33 += tstep·ChiDMat` — `tstep`, not `thstep` — on every
  level inside the `DO ILEV` loop; the defect then contains the implicit
  term automatically; **no explicit old-velocity `D` term exists**.
- H10: `rhs += tstep·g` once per coupling update, before the RHS freeze.
- H11, paper eq. (12): solve `[M_L + tstep·ChiD] δu = def` where
  `def = Δt·B(Δp)` (already `E013Sum3`'d and Dirichlet-filtered), then
  `u = ũ − δu`.

**Bit-identity mechanism for the correction (v2-review smaller items 3):**
`Chimera_CorrectVelocity` takes a **direct fast path** — when the weak
variant is inactive or `ChiD` is absent/zero it executes the *existing
diagonal loop verbatim* and reports `applied=.FALSE.` so the hook's ELSE
path semantics are preserved; bit identity is guaranteed by this fast
path, **not** by CG convergence (a CG solve is never bit-identical).
When the CG runs, its **scalar-product contract** is: partition-of-unity
weights `w_i = 1/c_i`, with `c_i` the share count obtained once at init
by `E013Sum` of a ones-vector; local weighted dots reduced with
`COMM_Summ`; Jacobi preconditioner `diag(M_L) + tstep·diag(ChiD)`;
mat-vec followed by `E013Sum3`.

**Normals** (stored with the boundary lists at classification):
`n_B` outward from the particle on `∂B_k`; `n_Γ` outward from the
atmosphere on `Γ_k`.

**Robin coupling, paper eq. (5d)** (all terms with `n = n_Γ`):
`σ̂·n − α(û·n)û = σ·n − α(u·n)u`. Weak form of the submesh momentum
equation (test function `v̂`, `v̂ = 0` on `∂B_k`): integration by parts
gives `… + ∫ σ̂:D(v̂) − ∮_Γ (σ̂·n)·v̂ ds = 0`; substituting (5d),

```
LHS  += − α ∮_Γ (û_k·n_Γ) (û_{k+1}·v̂) ds     (Picard: transport velocity û_k frozen)
RHS  += + ∮_Γ h·v̂ ds,   h = σ(u_bg,p_bg)·n_Γ − α (u_bg·n_Γ) u_bg
```

— i.e. the Picard-linearized α-matrix contribution enters with a
**minus** sign; for inflow across `Γ` (`û·n_Γ < 0`) it is a coercive
(backflow-stabilizing) addition, which is the intent of the
Dirichlet–Robin coupling. This complete residual replaces v2's
ambiguous "enters the matrix" (v2-review smaller items 4).

**Force/torque, paper eqs. (4a)/(4b):**
`F_k = −∮_{∂B_k} σ̂·n_B ds`, `T_k = −∮_{∂B_k} (x−X_k)×(σ̂·n_B) ds`.
Sign/orientation unit tests: hydrostatic linear pressure ⇒ discrete
buoyancy `F = −∇p·V_B`; cylindrical Couette ⇒ known traction and torque
sign; constant pressure over a closed surface ⇒ 0.

---

## 6. The hook contract

Milestones: **M1** static steady Chimera-S, **M2** static Chimera-W
(+ unsteady W validation), **M3** moving-W, **M4** time-accurate S.
All calls go through `CHIMERA_API`.

| # | File : site | Change | Ms |
|---|---|---|---|
| H1 | `QuadSc_main.f90` : `fluid_core`, after the `IF (enable_fbm)` block (~:597); call site inside a worker-only region because `LinSc%valP(NLMAX)%x` argument association is worker-only | `IF (Chimera_IsEnabled()) CALL Chimera_BeginStep(...)`; fatal if enabled-but-uninitialized | M1 (stub P0) |
| H2 | `QuadSc_boundary.f90` : `Boundary_QuadScalar_Def` (~:557) | `Chimera_ApplyBoundaryDef(...)` — sibling of the `FictKNPR` branch | M1 |
| H3 | same : `Boundary_QuadScalar_Val` (~:610) | `Chimera_ApplyBoundaryValues(...)` | M1 |
| H4 | same : `Boundary_QuadScalar_Mat` (~:700) + `_Mat_9` | `Chimera_FilterMatrixRows(...)` | M1 |
| H5 | `applications/q2p1_chimera/q2p1_chimera.f90` : after `init_q2p1_app` | `CALL Chimera_Initialize(ufile)` — **app-local** (paired with H6; no shared-init hook) | M1 |
| H6 | same : after the time loop | `CALL Chimera_Finalize()` — idempotent, partial-init-safe | M1 |
| P1 | `param_parser.f90` : `GDATNEW` SELECT CASE (~:1100) + validate after file close + echo (~:1345) | `CASE ("Chimera*")` → `CHIMERA_CONFIG`; `CHIMERA_VALIDATE_CONFIG()` when enabled; echo only when enabled | P0 |
| C1 | `ProjectFiles.cmake`, `GenerateLinkerFlags.cmake` | source lists, `ff_chimera`, ctests, layer-graph comment | P0 |
| H8 | `QuadSc_def.f90` : `Matdef_General_QuadScalar`, `idef==-1`, inside the level loop | `Chimera_AddMomentumMatrix(...)` (`tstep·ChiD`, §5) | M2 |
| H10 | `QuadSc_main.f90` : beside `AddGravForce()` (:630) | `Chimera_AddMomentumRHS(...)` (`tstep·g`) | M2 |
| H11 | `QuadSc_corrections.f90` : `Velocity_Correction` (~:30) | `Chimera_CorrectVelocity(applied)`; existing loop when `.NOT.applied` (fast path, §5); helper `Chimera_DefectFilter3` | M2 |
| H12 | driving app restart path | `Chimera_Write/ReadRestart` — versioned (`CHIMERA_RESTART_V1`); donor caches rebuilt, not restored | M2 |
| H15 | `QuadSc_assembly.f90` : `Create_CMat`; `QuadSc_def.f90` : `Create_ParCMat` | `Chimera_AddPressureMass(...)` — penalised lumped mass (capped) in `B^T M^-1 B`; the app rebuilds both after `Chimera_Initialize` | M2 |
| H13 | `QuadSc_main.f90` : in-step outer driver (paper §6 S-sequence) with full time-level state restoration (`valU^n`, `valP_old`/:785-786, `thstep`, stats) | default 1 ⇒ identical control flow | M4 |
| H14 | `QuadSc_handlers.f90` : `Init_Chimera_Handlers()` | forces → PE stepper | M3 |

**ZERO changes:** `Solve_General_QuadScalar`, the pressure path,
`FAC_GetForces`, `updateFBMGeometry`, `Init_QuadScalar_Stuctures`,
`OctTreeSearch.f90`, `UMFPackSolver`, `var_QuadScalar`, all COMMON
blocks, all existing `applications/*/app_init.f90`, all legacy F77
assembly kernels.

---

## 7. Submesh subsystem and data flow

Per coupling update (`Chimera_BeginStep`, workers only): cached outer
quadrature points → `CHI_EVAL_BG_AT_POINTS` → Robin data per §5 →
inner Dirichlet rows → owner solves → broadcast → fringe values →
forces (owner, reported once). `ChimeraOuterBC=dirichlet` (with pressure
gauge) is a diagnostic mode only; Robin is the acceptance path.

Submesh meshes: coarse shell `.tri` per particle shape from
`tools/chimera_meshgen/` *(project-lead decision)* — annulus (O-grid,
extruded; FAC) and sphere shells; later per-particle
`H_k = min(H_max, ½·nearest gap)` tables.
`applications/mesh_ref/Particle/*.tri` are unit-test fixtures.

---

## 8. Build, configuration, application

CMake: `chimera_config.f90` → `src_util` list; `add_library(ff_chimera)`
(Layer M); `chimera_api.f90` (+ later Layer-H files) → the two
`ff_quadLS_app` lists; `ff_chimera` → `FF_APPLICATION_LIBS`; tests under
`if(BUILD_TESTING)`. Runtime keys as in §3-L, documented in
`docs/md_docs/parameter_reference.md` + `chimera_usage.md`. Application:
`applications/q2p1_chimera/` (Phase 3) with app-local H5/H6.

---

## 9. Testing and verification

1. **Off-regression gate (every phase):** `q2p1_fc_ext_cylinder` passes
   unchanged (tolerance 0) on the default build; disabled-run protocol
   output byte-identical.
2. **Build-layering check (Phase 0):** `ff_chimera` builds standalone
   with no `ff_quadLS_app` module dependencies.
3. **Unit ctests** (`source/src_chimera/tests/`): `test_chi_geometry`
   (pinned vs `EL_Q1_MAP`; Newton on distorted hexes), `test_chi_eval`
   (Q2 delta property/partition of unity/quadratic exactness incl.
   gradients; P1), `test_chi_locator` (vs brute force),
   `test_chi_sparse_direct` (two interleaved instances — the
   independent-handles regression), *(Phase 2, green)* `test_chi_kernels`
   + annular Couette + §5 sign tests — convergence gate 2nd-order L2
   velocity (the Q1-geometry curved-boundary limit, shared with the
   production discretization; 3rd order would need isoparametric Q2
   geometry) and direct surface-traction torque monotone to < 1.5 %
   (boundary-flux evaluation is sub-quadratic; a variationally
   consistent force evaluation is a Phase-3+ accuracy option),
   *(Phase 3, green)* `test_chi_exchange`
   (1/2/3 ranks: exactness for quadratic fields, MIN-rank ownership,
   bit identity across ranks, missing-point report) + the
   **two-partition cut-cell marker test** `test_chi_markers` (serial
   emulation of the MAX merge on the shared interface dofs of a cut
   cell, for a body inside one partition and for a body crossing the
   interface), *(Phase 4)*
   `test_chi_algebra` (manufactured comparison vs paper eqs. (10)/(12);
   fast-path bit-identity of the correction).
4. **Physics case:** `q2p1_chimera_cylinder.yaml` pinning `ChimeraForce:`
   C_D/C_L for steady FAC vs the DFG band + FBM baseline cross-check.

---

## 10. Milestone roadmap

| Phase | Content | Exit criteria |
|---|---|---|
| 0 — Scaffold | `chimera_config`, `chimera_api` (fatal-stub `BeginStep`), P1, C1, H1 | builds; off-regression tol 0; layering check; disabled deck byte-identical |
| 1 — Services | `chi_geometry`, `chi_fem_eval`, `chi_locator`, `chi_sparse_direct` + tests | Phase-1 ctests green |
| 2 — Submesh subsystem *(reviewer-cleared; DONE 2026-09-02)* | `chi_kernels` (new reentrant kernels), `chi_legacy_mesh_adapter`, `chi_submesh/solver/forces`, meshgen | annular Couette 2nd-order L2 (Q1-geometry limit; measured 2.08), torque → −8π/3 (1.0 % at L3), Robin consistency, sign tests — all green (`chi-submesh-couette`) |
| 3 — **M1: static steady Chimera-S** *(DONE 2026-09-02; see §3 Phase-3 notes)* | `chi_markers`, `chi_exchange`, `chi_output`, `chi_coupling` (two-array markers); H1–H6 live; `q2p1_chimera` + vendored channel/annulus case; steady FAC | `chi-markers-cutcell`, `chi-exchange-np{1,2,3}` green; steady FAC `q2p1_chimera_cylinder` pinned (values in the baseline yaml, compared against the body-fitted `q2p1_fc_ext_cylinder` and the DFG band); worker-count invariance; off-regression exact |
| 4 — Static Chimera-W + unsteady validation *(DONE 2026-09-03; see §3 Phase-4 notes)* | `chi_penalty`; H8/H10/H11/H12/H15 (+ defect sibling); `test_chi_algebra` | `chi-algebra-serial` green; W vs S on steady FAC (`q2p1_chimera_cylinder_weak` pinned); worker-count invariance; fast path = the disabled/strong path runs the original loop (off-regression exact, strong anchor unchanged); restart round trip bit-identical; **unsteady Re 100 one-pass W NOT achieved at L2 — the strong variant diverges identically (RUNBOOK); carried to Phase 5/6 with the L3 background** |
| 5 — Arrays + periodicity *(DONE 2026-09-04; see §3 Phase-5 notes and `applications/q2p1_chimera/validation_cases/m2_hasimoto/RUNBOOK.md`)* | `chi_periodic` (minimum image / wrap), `SimPar@PeriodicLength`, H_k seeding + sphere-shell + periodic-partition tooling, `ChimeraSubStokes`, submesh body force, `ChimeraBulk` diagnostic, `test_chi_periodic`; halo upgrade deferred (interface unchanged) | `chi-periodic-serial` green; Hasimoto at D/h = 4: K_bal −0.8 % / K_meas −1.6 % (both variants; FBM needs D/h 24–48 for that); corner-placed sphere = centred to 5e-8; FAC anchors bit-identical; 8-sphere array S = W to 0.1 % with ±8 % per-sphere dispersion; Beetstra–Tenneti band: the quasi-ordered configuration the strong atmosphere rule allows at this resolution gives 2.36 (between SC 2.51 and the closures), the Chimera-W thin-atmosphere probe on 0.22 d gaps gives 2.00 (upper edge of the band, one realisation); ensembles carried to the campaign |
| 6 — Moving | submesh ALE, per-step reclassification, H14 (M3); H13 (M4) | ten Cate / FBM cross-checks; force continuity; `outer_iters=1` flow byte-identical |

---

## 11. Resolved decisions

- Always-compile + runtime key; meshgen tool; steady FAC as M1 target
  (project lead, 2026-09-01/02).
- Robin = acceptance path; Dirichlet = diagnostic (v2 reviewer).
- Penalty cadence via dirty flags (v2 reviewer).
- Hole-region background pressure DOFs left free (FBM precedent).
- "Undisturbed" = protocol-file byte identity when disabled.
- M1: Chimera and FBM particle mode mutually exclusive (fatal check in
  `Chimera_Initialize`).

## 12. Review disposition — v2 review (`chimera-plan-review-v2.md`, 2026-09-02)

| Item | Disposition |
|---|---|
| Blocker 1 — legacy kernels not COMMON-free | Adopted: new reentrant `chi_kernels` on `chi_geometry`/`chi_fem_eval`; legacy kernels never called for submeshes; `chi_legacy_mesh_adapter` restricted to mesh I/O as a temporary Layer-H bridge (§3). |
| Blocker 2 — signed markers vs `E013Max_SUPER` | Adopted: `marker_kind`/`marker_pid` two-array scheme with MAX-precedence + conditional identity sync; paper-exact cut-cell hole/fringe definition made normative; two-partition cut-cell test added (§3, §9). |
| Master guard / rank safety | Adopted: every facade op internally rank-safe; H5/H6 both app-local, no shared-init hook; H1 call-site placement justified by argument association (§2, §6). |
| CG contract / bit identity | Adopted: direct fast path guarantees identity; weighted (partition-of-unity) scalar products + `COMM_Summ` + Jacobi preconditioner specified (§5). |
| Robin residual written out | Adopted: full weak-form residual with the minus-signed Picard α-term (§5). |

## 13. Review disposition — v1 review (`chimera-plan-review.md`, summary)

All eight findings adopted in v2 and carried into v3: M1 coupling-scope
restriction; fully implicit penalty algebra (former H9 deleted);
eq.-(5d) sign and `−∮σ̂n_B` force convention; layer extraction
(`chi_geometry`) and mesh-I/O adapter; instance-based direct solver;
worker-only MPI contract; single facade; lifecycle/restart contract.
