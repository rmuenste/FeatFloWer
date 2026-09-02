# Review of `chimera-integration-design.md`

The plan has a strong direction—separate modules, runtime isolation,
instance-owned state, narrow hooks, and staged validation—but it is not
implementation-ready. Several design assumptions conflict with the paper or
the current codebase.

## Blocking findings

### 1. Milestone 1 does not implement Chimera-S's coupling algorithm

The plan calls `Chimera_BeginStep` before the background solve, performs only
one background solve, and treats later time steps as outer coupling iterations
(`chimera-integration-design.md:314-334`). The paper's Section 6 requires,
within the same time step:

1. a background solve;
2. a submesh solve;
3. another background solve with fringe constraints.

At least two outer iterations are required.

Time steps cannot substitute for these iterations in a time-accurate problem.
Moreover, DFG 2D-2 is explicitly unsteady because of vortex shedding, despite
being called "steady" in the plan. The current solver advances pressure and
history state during each call (`source/src_quadLS/QuadSc_main.f90:785-786`),
so simply calling it repeatedly would also be incorrect without restoring
time-level state.

Recommendation: either move the real outer-coupling driver into M1 or remove
DFG 2D-2 from M1 and restrict that milestone to genuinely steady or
pseudo-time-converged Stokes cases. Static particle geometry does not imply
steady fluid physics.

### 2. The proposed Chimera-W time discretization differs from the paper

The plan adds `thstep*D` to the implicit matrix and an explicit
`(1-theta)*D*u^n` defect term (`chimera-integration-design.md:321-325`). In the
paper's equations (9)-(12) and Appendix A, `D^(n+1)` is fully implicit:

```text
A(u^(n+1)) + D^(n+1)
M_L + delta_t D^(n+1)
```

It is not split with the convection and diffusion theta terms. In this solver,
`thstep` is first `delta_t*(1-theta)` and later `delta_t*theta`
(`source/src_quadLS/QuadSc_main.f90:598,651`). H8 and H9 would therefore
implement a different method.

Recommendation: specify the exact scaled algebra in the plan:

- Momentum matrix: `A += tstep*D`.
- Right-hand side: `rhs += tstep*g`.
- Do not add an explicit old-velocity `D` term.
- Correction system: `M_L + tstep*D`.

Add a manufactured algebraic test comparing the assembled equations with the
paper's equations (10) and (12).

### 3. The Robin and force conventions are incorrect or underspecified

The plan defines

```text
h = sigma n + alpha (u.n) u
```

at `chimera-integration-design.md:353`, whereas paper equation (5d) uses

```text
sigma_hat n - alpha (u_hat.n) u_hat
  = sigma n - alpha (u.n) u.
```

The sign is wrong. The atmospheric alpha term is also nonlinear; saying it
"enters the matrix" is insufficient. The Picard or Newton linearization must
be defined.

Likewise, the plan's `F = integral sigma_hat*n` at
`chimera-integration-design.md:190-193` is only correct if `n` is specifically
the atmosphere-domain outward normal on the inner surface. The paper defines
particle force as `-integral_(boundary B) sigma*n_B ds`.

Recommendation: define one canonical normal convention for every boundary
list and provide explicit formulas and sign tests using constant pressure and
Couette flow.

### 4. Layer M cannot currently have the claimed dependency properties

The plan proposes compiling all of `el_quadrature.f90` into `ff_chimera`
merely to reuse `EL_Q1_MAP` (`chimera-integration-design.md:144-152`). That
module imports `EL_HALO`, `EL_KERNEL_FUNCTIONS`, `EL_FIELDS`, MPI, and
configuration state, and contains several COMMON blocks
(`source/src_el/el_quadrature.f90:1-45`). It is already compiled into
`ff_quadLS_app` (`cmake/modules/ProjectFiles.cmake:294-330`). Adding it to
`ff_chimera` creates module-order or circular-dependency problems and duplicate
objects.

The proposed submesh loader also uses `readTriCoarse` and `refineMesh`, both of
which currently `USE var_QuadScalar`
(`source/src_mesh/mesh_refine.f90:216-218,1210-1213`). Thus the statement that
Layer M contains everything that "does not need `var_QuadScalar`" is not true
transitively.

Recommendation:

- Extract the Q1 map into a small dependency-free geometry module.
- Either refactor the mesh loader/refiner to receive their required
  configuration explicitly or put a `chi_legacy_mesh_adapter` in Layer H.
- Add a CMake target-dependency diagram and a build-only Phase 0 test before
  designing further modules around the proposed layering.

### 5. The existing UMFPACK wrapper cannot provide per-submesh solver handles

The plan assumes each `tChimeraSubmesh` owns an independent symbolic and
numeric factorization (`chimera-integration-design.md:181-188`). The current
wrapper has module-global `symbolic` and `numeric` handles, imports
`var_QuadScalar`, always performs symbolic and numeric factorization together,
and mutates CSR indices in place from one-based to zero-based
(`source/UmfpackSolver.f90:1-73`).

Multiple submeshes would overwrite each other's factorization. The submesh
solver can also conflict with the background coarse solver.

Recommendation: introduce a new instance-based `tSparseDirectSolver` wrapper
with owned CSR storage and independent symbolic and numeric handles. Do not use
`UMFPackSolver` directly as the submesh backend. Also define a pressure gauge
for the all-Dirichlet fallback, whose saddle-point matrix otherwise has a
pressure nullspace.

### 6. The M1 MPI protocol is internally inconsistent

The interface says every rank has a replicated query list and then performs an
`MPI_Allreduce`, but also says the master participates with zero points
(`chimera-integration-design.md:195-202`). All participants in an
`MPI_Allreduce` must use the same count. In addition, `MPI_COMM_SUBS` explicitly
excludes rank 0 (`source/src_mpi/pp3d_mpi.f90:122-134`), so the plan must choose
between worker-only and world collectives.

A "lowest-rank tie-break" also cannot be achieved by a simple sum reduction
without first selecting an owner.

Recommendation for M1:

- Use `MPI_COMM_SUBS`.
- Exclude the master from the service with an immediate no-op return.
- Replicate identical point ordering on all worker ranks.
- First reduce the candidate owner with `MPI_MIN`; then sum values masked by
  the selected owner.
- Replicate all small static submeshes on every worker for M1. Defer
  intersecting-rank communicators until scaling requires them.

## Important modularity findings

### 7. The facade boundary is contradicted by the hook design

The plan says `chi_coupling` is the only Chimera module existing code uses
(`chimera-integration-design.md:211-225`), but the dependency graph exposes
`chi_coupling`, `chi_penalty`, and `chimera_config`, while H2-H4 reach directly
into `chi_state%KNPR` (`chimera-integration-design.md:272-275,315-317`). This
leaks representation into legacy solver files.

Use a single public facade such as `CHIMERA_API` exposing operations like:

- `Chimera_IsEnabled`
- `Chimera_Initialize`
- `Chimera_ApplyBoundaryDef`
- `Chimera_ApplyBoundaryValues`
- `Chimera_FilterMatrixRows`
- `Chimera_AddMomentumTerms`
- `Chimera_CorrectVelocity`
- `Chimera_Finalize`

Keep `tChimeraState`, markers, penalty storage, and donor caches private. The
parser can be the sole exception that imports `CHIMERA_CONFIG`.

### 8. Lifecycle, restart, and ownership contracts are incomplete

`Chimera_Release` is listed but has no hook in the complete hook table.
Restart and checkpoint behavior is absent, despite the submesh solution being
part of the coupled time-dependent state. The existing EL subsystem has
versioned state serialization (`source/src_el/el_fields.f90:120-169`).

The plan should define:

- where finalization is called;
- idempotent partial-initialization cleanup;
- restart format and versioning;
- whether donor caches are restored or rebuilt;
- which fields must be restored for time-accurate continuation;
- resource ownership and finalization for mesh and direct-solver handles.

## Answers to the open review questions

### 1. Is the Allreduce exchange acceptable for M1?

Yes, after correcting the communicator, equal-count, and owner-selection
contract described above.

### 2. Are nodal basis derivatives sufficient for the Robin data?

Q2 basis derivatives and element-local physical P1 pressure evaluation are the
correct evaluation of the discrete background solution. Robin should remain
the M1 acceptance path; Dirichlet is useful as a diagnostic comparison, not as
a substitute for validating the paper's method.

### 3. Should the penalty cadence follow `MatrixRenewal`?

No. Tie it to explicit state invalidation:

- Static geometry: assemble `D` once.
- Moving geometry: mark `D` dirty after particle motion or reclassification.
- Updated submesh solution: rebuild `g` every coupling update.
- Multigrid hierarchy: rebuild coarse penalty operators whenever the finest
  geometry changes.

## Overall verdict

Retain the module-decomposition concept, runtime-off regression gate, and
phased tests. Revise the coupling driver, mathematical scaling, direct-solver
abstraction, MPI contract, and facade boundary before starting Phase 0.
