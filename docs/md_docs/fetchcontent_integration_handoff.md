# Handoff: FetchContent Integration for Solver Libraries

Goal: make the solver dependencies introduced during the Phase 0–3 evaluation
reproducible via CMake FetchContent instead of hand-pointed site installs,
without breaking any existing `USE_EXTERNAL_*` path.

Decided scope (user-approved 2026-07-26):

1. **Hypre 2.33.0 via FetchContent** — primary deliverable.
2. **SuiteSparse 7.x (UMFPACK) via FetchContent** — secondary deliverable.
3. **MKL PARDISO → `find_package(MKL CONFIG)`** — never fetched (proprietary
   binary); replace the hand-rolled static library list with the `MKL::MKL`
   target from `MKLConfig.cmake`.
4. **MUMPS: unchanged.** Keep the existing `USE_EXTERNAL_MUMPS` find-path
   exactly as is. Do not add a fetch path (Phase 2 verdict: not competitive).

## Background you must read first

- `docs/md_docs/solver_phase1_direct_solvers.md` — how `USE_EXTERNAL_SUITESPARSE`
  and `USE_MKL_PARDISO` are currently wired, and the **exact working configure
  command lines** (CGAL_DIR / GMP / MPFR flags are mandatory on this machine).
- `docs/md_docs/solver_phase2_velocity_coarse.md` — the external MUMPS wiring
  (context only; do not modify).
- `docs/md_docs/solver_phase3_fine_level_pressure.md` — why Hypre matters and
  the upcoming `HYPRE_ENABLE_GPU_AWARE_MPI` experiment this work must enable.
- Current CMake integration points: root `CMakeLists.txt` (option blocks for
  `USE_EXTERNAL_SUITESPARSE`, `USE_MKL_PARDISO`, `USE_EXTERNAL_MUMPS`, and the
  Hypre section), `cmake/modules/GenerateLinkerFlags.cmake` (string-based
  library lists: `FF_UMFPACK_LINK_LIBS`, `FF_MKL_PARDISO_LIBS`,
  `MUMPS_LIBRARY_LIST` branches), `cmake/modules/GenerateIncludeFlags.cmake`,
  `cmake/modules/ProjectFiles.cmake` (le_solvers sources incl.
  `ff_timing.f90`, `PardisoSolver.f90`, `umf4_f77wrapper_port.c`).

## Design

### Provider pattern

Create `cmake/modules/FFDependencies.cmake`, included from the root
`CMakeLists.txt`. For each fetchable dependency use the CMake ≥3.24 form:

```cmake
FetchContent_Declare(HYPRE
  URL https://github.com/hypre-space/hypre/archive/refs/tags/v2.33.0.tar.gz
  URL_HASH SHA256=<compute and pin — download once and sha256sum it>
  SOURCE_SUBDIR src
  FIND_PACKAGE_ARGS NAMES HYPRE)
```

Semantics: a system/site install found by `find_package` wins; otherwise the
pinned tarball is fetched and built as a subproject. Preserve the existing
explicit-dir options (`EXTERNAL_SUITESPARSE_DIR`, hypre dir/prefix variables)
as `CMAKE_PREFIX_PATH`/hint inputs so current workflows keep working
unchanged. Add an opt-out (`FF_FETCH_DEPENDENCIES=OFF` or per-dep
`FF_HYPRE_PROVIDER=system|fetch`) so a build can *require* the system install
and fail loudly instead of silently fetching.

`cmake_minimum_required` is currently 3.18 (login node has 3.31.8). Either
bump the project minimum to 3.24 or gate the FetchContent module on
`CMAKE_VERSION VERSION_GREATER_EQUAL 3.24` with a clear fatal error when a
fetch is requested on older CMake. Bumping is acceptable; note it in the doc.

### Hypre specifics

- Hypre 2.33.0's CMakeLists lives in `src/` → `SOURCE_SUBDIR src`.
- The option set that is known-good on this machine is captured in the cache
  of the manual build: `stage1-work/hypre-2.33.0-build/CMakeCache.txt`
  (`HYPRE_ENABLE_MPI=ON`, `HYPRE_ENABLE_HYPRE_BLAS/LAPACK=ON`, CUDA family
  toggles). Drive the same options from FeatFloWer-level settings:
  - `USE_HYPRE_CUDA` → `HYPRE_ENABLE_CUDA`, `HYPRE_ENABLE_CUBLAS`,
    `HYPRE_ENABLE_CURAND`, `HYPRE_ENABLE_CUSPARSE`, `HYPRE_ENABLE_CUSOLVER`,
    `HYPRE_ENABLE_CUDA_STREAMS` (all ON as in the reference cache), plus the
    existing sm_80 arch setting.
  - New option `FF_HYPRE_GPU_AWARE_MPI` (default OFF) → `HYPRE_ENABLE_GPU_AWARE_MPI`.
    This is the enabler for the next GPU experiment; it only needs to plumb
    through, no further logic.
- A local source tree already exists at `stage1-work/hypre-2.33.0-src`; use it
  via `FETCHCONTENT_SOURCE_DIR_HYPRE` to develop/verify without re-downloading,
  but the committed declaration must carry the URL + SHA256 so it works on a
  clean clone.
- Link the fetched target (`HYPRE`) rather than extending the string list in
  `GenerateLinkerFlags.cmake`; for the found-package branch keep whatever the
  current mechanism is. A thin interface target (e.g. `ff::hypre`) that both
  providers populate, consumed where the hypre libs are currently linked, is
  the preferred shape — but keep the refactor minimal: do not convert
  unrelated libraries to targets.

### SuiteSparse specifics

- Fetch SuiteSparse v7.x (pick the latest stable 7.x tag, pin URL + SHA256)
  with `SUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;camd;colamd;ccolamd;cholmod;umfpack"`.
- The fetched build provides `SuiteSparse::UMFPACK` etc.; wire these into
  `FF_UMFPACK_LINK_LIBS` / a target equivalent. The vendored
  `umf4_f77wrapper_port.c` compiles against the fetched headers — reuse the
  existing `USE_EXTERNAL_SUITESPARSE` plumbing in `ProjectFiles.cmake`
  (include dir + wrapper source) with the fetched include dir.
- Migration note for the doc: this moves from site UMFPACK 5.7.8 to 7.x; the
  `umf4*` C API is stable, but flag it so a benchmark re-bless can be run
  later. Needs a BLAS; on this machine OpenBLAS is available (see the MUMPS
  block in `CMakeLists.txt` for how it is found) — prefer letting SuiteSparse's
  own find logic run, and only intervene if configure fails.

### MKL specifics

- Replace the hand-assembled `--start-group`/`--end-group` static list with
  `find_package(MKL CONFIG)` using `MKL_DIR=$MKL_PARDISO_DIR/lib/cmake/mkl`
  (site root: `/sfw/intel/2024.0.1/mkl/latest`), `MKL_LINK=static`,
  `MKL_INTERFACE=lp64`, `MKL_THREADING=sequential`, then link `MKL::MKL`.
  Keep `USE_MKL_PARDISO` as the user-facing option and `MKL_PARDISO_AVAIL`
  as the preprocessor define. If `MKLConfig.cmake` proves incompatible with
  this project's Fortran link line, it is acceptable to keep the existing
  list and document why — this sub-task is opportunistic cleanup, not a
  blocker.

### Offline/cluster behavior (must document)

In the new doc `docs/md_docs/fetchcontent_dependencies.md` (create it, index
it in `docs/md_docs/README.md`):

- Downloads happen at configure time on the login node (network exists there;
  compute nodes have none). Configure must therefore never run inside a Slurm
  job unless sources are pre-staged.
- Pre-staging: `FETCHCONTENT_SOURCE_DIR_<UCNAME>` overrides, and a shared
  `FETCHCONTENT_BASE_DIR` (e.g. `stage1-work/_deps-cache`) so multiple build
  trees don't re-download/re-build.
- The provider/opt-out options and the pinned versions + hashes.

## Constraints

- **Work directly in `/data/warehouse17/rmuenste/code/FF-GPU/FeatFloWer` on
  the current checkout. Do NOT create a worktree or branch and do NOT
  commit** — the tree carries uncommitted Phase 0–3 changes that are part of
  this work's context; a worktree would silently lose them.
- Do not submit Slurm jobs. Verification is configure + build on the login
  node only.
- Do not modify the benchmark harness
  (`tools/cluster_scripts/solver_baseline_phase0.slurm`), the solver Fortran
  sources, or anything under `source/` except where a link/include change
  strictly requires it.
- Every pre-existing configure path must keep working:
  `USE_EXTERNAL_SUITESPARSE + EXTERNAL_SUITESPARSE_DIR`, `USE_MKL_PARDISO`,
  `USE_EXTERNAL_MUMPS`, plain vendored builds with all options OFF.

## Environment (no modules needed)

```
export PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/bin:/sfw/openmpi/gcc13.2.x/4.1.6/ucx-threaded-noverbs/bin:$PATH
export LD_LIBRARY_PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/lib64:$LD_LIBRARY_PATH
# CUDA builds additionally: /sfw/cuda/12.4/bin and its lib64
```

CGAL configure flags (mandatory, exact values in
`solver_phase1_direct_solvers.md`): `-DCGAL_DIR=/sfw/cgal/gcc13.2.0/lib64/cmake/CGAL`
plus explicit `-DGMP_*`/`-DMPFR_*`/`-DGMPXX_*` pointing into the gcc prefix.

Reference build trees (do not modify them; use as configure-line references):
`stage1-work/build-phase2` (current CPU evaluation build),
`stage1-work/build-gpu-phase3` (CUDA build).

## Verification (definition of done)

1. New build tree `stage1-work/build-fetchdeps-cpu`: configure with fetched
   Hypre + fetched SuiteSparse + `USE_MKL_PARDISO` + `ENABLE_SOLVER_TIMING=ON`
   (no `EXTERNAL_SUITESPARSE_DIR`, no pre-set hypre paths). Build targets
   `q2p1_fc_ext` and `q2p1_fac3d` to completion. Check exit codes directly —
   never through a `| tail` pipe.
2. Regression: reconfigure-from-scratch a copy of the `build-phase2`
   configure line (external/site everything) in a scratch build dir and
   confirm it still configures and builds one app. Delete the scratch dir
   afterwards.
3. `ldd` (or `readelf -d`) the fetched-build `q2p1_fc_ext` binary to confirm
   it links the fetched libraries, not stray site ones.
4. Optional if time permits: `stage1-work/build-fetchdeps-gpu` configure-only
   check with `USE_HYPRE_CUDA=ON -DFF_HYPRE_GPU_AWARE_MPI=ON` to prove the
   toggle plumbs through (full CUDA build not required).

## Deliverables

- `cmake/modules/FFDependencies.cmake` + root `CMakeLists.txt` integration.
- Linker/include module updates as needed (minimal).
- `docs/md_docs/fetchcontent_dependencies.md` + README index entry.
- A short final report: what was built, exact configure lines used, any
  deviations from this plan and why, and anything left open.
