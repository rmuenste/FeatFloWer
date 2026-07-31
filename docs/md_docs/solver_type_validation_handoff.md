# Handoff: Early Validation of Solver Types vs Compiled Capabilities

Goal: when `q2p1_param.dat` selects a coarse/fine solver type whose library is
not compiled in, the run must fail **immediately after parameter read** with a
clear message naming the required CMake option — instead of the current
behavior, which is a late and sometimes messy failure.

## Current behavior (verified 2026-07-26)

Solver type dispatch for the pressure Q2/P1 system lives in
`source/src_quadLS/QuadSc_mg.f90` (`mgCoarseGridSolver_P`, and a velocity
counterpart `mgCoarseGridSolver_U`), the type-10 fine-level dispatch in
`source/src_quadLS/QuadSc_def.f90` (`Solve_General_LinScalar`, ~line 2044),
and the direct-solver factorization in
`source/src_quadLS/QuadSc_solver_coarse.f90` (`Setup_UMFPACK_CoarseSolver`).

Per type, in a build without the corresponding library:

- Types 1–4: always available (vendored UMFPACK/AMD built unconditionally).
- Type 5 (MUMPS, `#ifdef MUMPS_AVAIL`): `'MUMPS is not available!'` + `STOP`
  on all ranks at first coarse solve (`QuadSc_mg.f90` ~2175).
- Types 7/8 (Hypre coarse, `#ifdef HYPRE_AVAIL`): message + `STOP` all ranks
  (`QuadSc_mg.f90` ~2276).
- Type 10 (Hypre fine, `#ifdef HYPRE_AVAIL`): message + `STOP` all ranks
  (`QuadSc_def.f90` ~2049).
- Type 9 (MKL PARDISO, `#ifdef MKL_PARDISO_AVAIL`): `STOP` inside
  **master-only** branches (`QuadSc_solver_coarse.f90` ~113 factorize,
  `QuadSc_mg.f90` ~2089 solve). Workers stay blocked in MPI; Open MPI's
  runtime eventually kills the job, but it is not a clean abort and can hang
  on other MPI stacks.
- **Unknown type (e.g. 6, 11, typos): silent no-op.** The dispatch is a chain
  of independent IFs; no branch matches, the coarse solve does nothing, and
  MG runs without a coarse correction — degraded/diverging solves with no
  error message. This is the worst failure mode and must be caught.

All failures happen at the first solve (or factorization), i.e. after
startup, partitioning, and matrix assembly — minutes into a large job.

## Task

1. **Add a validation routine**, e.g. `ValidateSolverTypes(...)` in a
   sensible existing utility/module location under `source/` (or a small new
   file added to `cmake/modules/ProjectFiles.cmake` — but see the
   coordination section: if adding a source file requires editing CMake
   modules, prefer placing the routine in an existing compiled source file
   such as the module that owns the parameter structures, to avoid touching
   the build system at all).
   - Called on **all ranks** as soon as both `Pres@MGCrsSolverType` and
     `Velo@MGCrsSolverType` (fields `...%prm%MGprmIn%CrsSolverType` of the
     pressure/velocity solver structures) are populated from
     `q2p1_param.dat` — find the parameter-read call site used by the q2p1
     applications (all apps share the same init path; hook it in shared
     init code, not per-app).
   - Validates against the known type list per component. Determine the
     velocity dispatcher's actually-supported set from
     `mgCoarseGridSolver_U` (do not assume it equals the pressure set).
     Type 10 is valid for pressure only.
   - Availability per `#ifdef`: 5 → `MUMPS_AVAIL` (CMake `USE_MUMPS` /
     external MUMPS), 7/8/10 → `HYPRE_AVAIL` (`USE_HYPRE`), 9 →
     `MKL_PARDISO_AVAIL` (`USE_MKL_PARDISO`).
   - On violation: master rank writes a self-explanatory message (which type
     was requested, for which component, what it needs — e.g. "rebuild with
     -DUSE_HYPRE=ON" — or, for unknown types, the list of valid values),
     then the job aborts via `MPI_Abort(MPI_COMM_WORLD, ...)` with nonzero
     code. Use a barrier-free pattern: every rank evaluates the same check
     on the same values, master prints, all call MPI_Abort (MPI_Abort from
     any rank kills the job, so ordering is not critical — just make sure
     the message is flushed before the abort on the printing rank).
   - Validation must run **before** the type-7/8/10 hierarchy reconfiguration
     in `QuadSc_initialization.f90` (lines ~24/709/959: master NLMAX=NLMIN,
     bMasterTurnedON=.FALSE.) so a doomed run never reshapes anything.
2. **Convert the existing 'not available' `STOP`s to `MPI_Abort`** (the ones
   listed above). They become defense-in-depth behind the new validation;
   the master-only type-9 ones are the important conversions. Keep the
   messages. Do NOT touch other STOPs (e.g. the Phase-2 velocity-MUMPS
   numbering guard in `mgCoarseGridSolver_U` stays a STOP or may also become
   MPI_Abort — your choice, it is a collective code path).
3. **Doc**: write `docs/md_docs/solver_type_validation.md` (what is
   validated, when, example messages). Do NOT edit `docs/md_docs/README.md`
   (see coordination); note in your report that the index entry is pending.

## Coordination — another agent is working in this checkout RIGHT NOW

A second agent is concurrently implementing a FetchContent build-system
change. Hard rules to avoid collisions:

- Do NOT edit: root `CMakeLists.txt`, anything under `cmake/modules/`,
  `docs/md_docs/README.md`, `docs/md_docs/fetchcontent*.md`. If your
  preferred design needs a new source file registered in
  `ProjectFiles.cmake`, change the design: put the code in an
  already-registered file.
- Do NOT reconfigure or build in `stage1-work/build-phase2`,
  `build-gpu-phase3`, or any `build-fetchdeps*` tree.
- Create your own build tree(s), e.g. `stage1-work/build-typevalid`.
  If `cmake` configure fails with errors that look like half-edited CMake
  files (the other agent's work in flight), wait and retry a couple of
  times before reporting it as a blocker.

## General constraints

- Work directly in `/data/warehouse17/rmuenste/code/FF-GPU/FeatFloWer`.
  No git worktree, no branch, no commits, no checkout/restore/reset/stash —
  the tree holds uncommitted work that must survive.
- No Slurm submissions (no sbatch/srun). Verification is compile-based on
  the login node.
- Fortran style: match the surrounding fixed-habit F90 style of
  `src_quadLS` (uppercase keywords, existing USE patterns). The MPI module
  usage pattern in these files is `mpif.h`/`USE var_QuadScalar` etc. — copy
  whatever the file you edit already does.
- `#ifdef FF_SOLVER_TIMING` instrumentation must remain untouched.

## Environment (no modules needed)

```
export PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/bin:/sfw/openmpi/gcc13.2.x/4.1.6/ucx-threaded-noverbs/bin:$PATH
export LD_LIBRARY_PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/lib64:$LD_LIBRARY_PATH
```

Known-good configure lines: read `stage1-work/build-phase2/CMakeCache.txt`
(full-featured: SuiteSparse/PARDISO/MUMPS/Hypre/timing) and
`docs/md_docs/solver_phase1_direct_solvers.md` (CGAL_DIR/GMP/MPFR flags are
mandatory). For the minimal build, drop the optional solver options but keep
the CGAL flags if configure demands them.

## Verification (definition of done)

1. `stage1-work/build-typevalid-full`: configure with all solver options ON
   (mirror build-phase2's cache values), build `q2p1_fc_ext` and
   `q2p1_fac3d`. Both must compile and link; exit codes checked directly,
   never through a `| tail` pipe.
2. `stage1-work/build-typevalid-min`: configure with NO optional solver
   libraries (no USE_HYPRE/USE_MUMPS/USE_MKL_PARDISO/USE_EXTERNAL_*), build
   `q2p1_fc_ext`. This proves the `#else` paths of your validation compile.
3. If cheaply possible without Slurm (small mpirun on the login node using
   an existing app `_data` setup), demonstrate one negative case: minimal
   build + parameter file requesting type 9 or 10 → immediate abort with the
   new message. If a local run is not practical, say so; compile-level
   verification is the requirement, the runtime demo is a bonus.

## Deliverables

- Validation routine + call site wiring (shared init path, all q2p1 apps).
- STOP → MPI_Abort conversions listed above.
- `docs/md_docs/solver_type_validation.md` (README index entry deferred).
- Final report: files touched, exact configure lines, build exit codes,
  the exact wording of the new error messages, any deviations and why.
