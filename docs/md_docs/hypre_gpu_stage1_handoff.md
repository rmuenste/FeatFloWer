# Handoff Plan — Stage 1: GPU-Backed Hypre Coarse-Grid Solver

**Target system:** remote node with one NVIDIA A100 (compute capability `sm_80`), full CUDA
stack available via the `module` system, GCC/GFortran + MPI toolchain.

**Goal:** the existing Hypre coarse-grid pressure solve (`Pres@MGCrsSolverType = 7` or `8`)
runs on the A100 instead of the CPU. No solver-logic changes beyond GPU-compatible
BoomerAMG parameters.

**Scope guard:** MUMPS/UMFPACK paths and everything else stay untouched. The CPU build
(`-DUSE_HYPRE_CUDA=OFF`) must remain identical in behavior.

---

## Background (read first)

- Hypre is vendored at `extern/libraries/hypre/src` (v2.25.0, `HYPRE_NUMBER 22500` in its
  `CMakeLists.txt:16`). FeatFloWer builds it via `add_subdirectory` when `-DUSE_HYPRE=ON`
  (root `CMakeLists.txt:431-438`), defines `-DHYPRE_AVAIL`, and compiles the wrapper
  `source/HypreSolver.f90`.
- **Version strategy:** v2.25.0 (mid-2022) belongs to hypre's early GPU era; many GPU
  fixes (IJ host/device pointer handling, device SpGemm/setup, PMIS on device,
  `HYPRE_SetGpuAwareMPI`) landed in later releases. Stage 1 therefore uses a **current
  hypre release (>= 2.31, ideally latest 2.3x) built externally on the remote system as
  the primary path** (Step 2a/2b). Building the vendored 2.25 with CUDA (Step 2c) is the
  fallback only. The wrapper's API surface (`HYPRE_Init`, IJ interface,
  BoomerAMG/GMRES/PCG, `HYPRE_SetMemoryLocation/ExecutionPolicy/SpGemmUseVendor`,
  `HYPREf.h`) is stable across 2.25 -> 2.3x, so `source/HypreSolver.f90` should compile
  unchanged; report any signature drift you do hit.
- For reference, the vendored hypre's own CMake exposes the GPU knobs: `HYPRE_WITH_CUDA`,
  `HYPRE_ENABLE_UNIFIED_MEMORY`, `HYPRE_CUDA_SM` (default `70`; the A100 needs `80`).
- The wrapper **already requests device execution**:
  `HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE)`, `HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE)`,
  `HYPRE_SetSpGemmUseVendor(0)` at `source/HypreSolver.f90:24-29` (repeated at 182-187 and
  337-342 for the GMRES/PCG variants). In the current CPU-only build these are no-ops; once
  hypre is compiled with CUDA they become live.
- Hypre runs on a **sub-communicator excluding rank 0** (`MPI_COMM_split` at
  `source/HypreSolver.f90:33-39`); rank 0 is the FeatFloWer master and does no solve work.
  All worker ranks will share the single A100.
- Matrix/RHS marshalling happens in `source/src_quadLS/QuadSc_solver_hypre.f90`
  (`Setup_HYPRE_CoarseLevel_Full` = type 7, `Setup_HYPRE_CoarseLevel_Geometric` = type 8)
  into plain Fortran arrays `myHYPRE%rows/ncols/cols/values/rhs/sol` (`tHYPRE` defined in
  `source/src_util/types.f90:433`). These **host** arrays are passed to
  `HYPRE_IJMatrixSetValues` / `HYPRE_IJVectorSetValues` — this is the one real risk point,
  see Step 4.
- Smoke-test app exists: `applications/hypre_test/`. Runtime solver selection:
  `Pres@MGCrsSolverType` in `_data/q2p1_param.dat` (7 = GMRES+AMG on the full coarse
  matrix, 8 = AMG on the /16-coarsened auxiliary system).

## Step 1 — Environment

```bash
module load cuda gcc openmpi   # adjust to site names; CUDA >= 11.x
nvidia-smi                     # confirm A100 visible, note driver/CUDA version
which nvcc && nvcc --version
```

CUDA-aware MPI is *not* required (hypre stages communication buffers through the host by
default), but note whether the MPI module is CUDA-aware for later tuning.

## Step 2a — Obtain a current GPU-enabled hypre (primary path)

Check the module system first (`module avail hypre` — only use a module that was built
with CUDA *and* unified memory; verify with `grep HYPRE_USING_UNIFIED_MEMORY` on its
installed `HYPRE_config.h`). Otherwise build the latest 2.3x release from source — this
is quick (~minutes) and gives full control:

```bash
git clone --depth 1 --branch v2.33.0 https://github.com/hypre-space/hypre.git hypre-gpu
cd hypre-gpu/src
./configure --prefix=$HOME/opt/hypre-2.33-cuda \
            --with-cuda --with-gpu-arch=80 \
            --enable-unified-memory \
            CC=mpicc CXX=mpicxx FC=mpif90
make -j16 install
```

(Adjust the tag to the newest release available. `--enable-unified-memory` matters: the
FeatFloWer wrapper passes host-side Fortran arrays, see Step 4.)

## Step 2b — CMake wiring in FeatFloWer for an external hypre

Extend the hypre block in the root `CMakeLists.txt` (lines 431-438) with a
`USE_EXTERNAL_HYPRE` path so the vendored copy is bypassed:

```cmake
option(USE_EXTERNAL_HYPRE "Link a pre-installed hypre instead of the vendored copy" OFF)
set(EXTERNAL_HYPRE_DIR "" CACHE PATH "Install prefix of the external hypre")
if(USE_HYPRE)
  add_definitions(-DHYPRE_AVAIL)
  if(USE_EXTERNAL_HYPRE)
    find_library(HYPRE_LIB HYPRE PATHS ${EXTERNAL_HYPRE_DIR}/lib NO_DEFAULT_PATH)
    set(HYPRE_LIBRARIES ${HYPRE_LIB})
    include_directories(${EXTERNAL_HYPRE_DIR}/include)
  else()
    add_subdirectory(extern/libraries/hypre/src)
    set(HYPRE_LIBRARIES HYPRE)
  endif()
  set(src_q2p1 ${src_q2p1} ${CMAKE_SOURCE_DIR}/source/HypreSolver.f90)
endif()
```

Two auxiliary CMake files reference the vendored paths and need the same switch:
`cmake/modules/GenerateIncludeFlags.cmake:28-39` (hypre include dirs — must point at
`${EXTERNAL_HYPRE_DIR}/include`, where `HYPREf.h` also lives) and
`cmake/modules/GenerateLinkerFlags.cmake:125-128` (link line — use `${HYPRE_LIBRARIES}`
and append the CUDA runtime libs: `cudart cusparse curand cublas stdc++`, from the
CUDA toolkit's `lib64`).

Configure a fresh build directory:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON \
      -DUSE_HYPRE=ON -DUSE_EXTERNAL_HYPRE=ON \
      -DEXTERNAL_HYPRE_DIR=$HOME/opt/hypre-2.33-cuda ..
make -j16 2>&1 | tee build.log
```

C++ device code linked into Fortran executables typically needs `-lstdc++`; add it to the
link line if you see missing `__cxa_*`/`std::` symbols.

## Step 2c — Fallback: build the vendored 2.25 with CUDA

Only if the external build hits a blocker (e.g. wrapper/API incompatibility). Add a
`USE_HYPRE_CUDA` option that forwards to the vendored hypre's own options **before**
`add_subdirectory`:

```cmake
option(USE_HYPRE_CUDA "Build vendored hypre with CUDA backend" OFF)
...
if(USE_HYPRE_CUDA)
  set(HYPRE_WITH_CUDA ON CACHE BOOL "" FORCE)
  set(HYPRE_ENABLE_UNIFIED_MEMORY ON CACHE BOOL "" FORCE)
  set(HYPRE_CUDA_SM "80" CACHE STRING "" FORCE)   # A100
endif()
add_subdirectory(extern/libraries/hypre/src)
```

then `cmake -DUSE_HYPRE=ON -DUSE_HYPRE_CUDA=ON ..`. Expect more GPU rough edges at
runtime (this is a 2022 release) and the same link-line fixes as in Step 2b if hypre's
target does not propagate CUDA libs transitively.

## Step 3 — Smoke test with the standalone app

Build and run `applications/hypre_test/` with 2-4 ranks. Verify with `nvidia-smi dmon`
(or `nsys profile`) that kernels actually execute on the GPU during the solve. If the app
does not exercise the same code path as the wrapper, proceed to Step 6's application test
regardless — this step is only a fast sanity check of the library build.

## Step 4 — The host-pointer risk (the one place code changes may be needed)

With a device-memory build, hypre's `IJMatrixSetValues`/`IJVectorSetValues` may assume the
input arrays live in device-accessible memory. The wrapper passes ordinary Fortran heap
arrays. Three outcomes, in order of likelihood:

1. **It just works** — the unified-memory build plus hypre's internal pointer-location
   checks handle host input. Current 2.3x releases are much better at this than 2.25.
   Test first; do nothing if clean.
2. **Segfault / illegal access inside SetValues** — fix by making the marshalling arrays
   managed. Cleanest approach: a small C shim compiled into the project:

   ```c
   void* ff_malloc_managed(size_t bytes) { void* p; cudaMallocManaged(&p, bytes); return p; }
   void  ff_free_managed(void* p)        { cudaFree(p); }
   ```

   and in `source/src_quadLS/QuadSc_solver_hypre.f90` replace the `ALLOCATE` of
   `myHYPRE%rows/ncols/cols/values/rhs/sol` with `c_f_pointer`-wrapped managed allocations
   (guard with `#ifdef HYPRE_CUDA` so the CPU build is untouched). The arrays are written
   on the host once per solve and read by hypre — managed memory handles both sides.
3. **Silent wrong results** — compare iteration counts and final residual norms (both
   printed by the wrapper) against a CPU-only reference build on identical input; treat
   divergence beyond AMG-variant noise as outcome 2 in disguise.

## Step 5 — GPU-compatible BoomerAMG parameters

Hypre executes only certain AMG components on device; the current settings in
`source/HypreSolver.f90` were chosen for CPU and will silently fall back to host paths.
Change, in all three solver routines (`myHypre_Solve`, `myHypreGMRES_Solve`,
`myHyprePCG_Solve`):

| Setting | Current | GPU-capable value |
|---|---|---|
| Coarsening | HMIS (`HYPRE_BoomerAMGSetCoarsenType` 10) | **PMIS (8)** |
| Interpolation | type 13 | **type 6** (extended+i) or 3, with `SetTruncFactor` ~0.25 |
| Relaxation | type 8 (hybrid sym-G-S) | **type 18 (l1-Jacobi)** or 7 (Jacobi) |
| Aggressive coarsening | (as set) | keep, with `SetAggInterpType` 5 or 7 if levels > 0 |

Keep `HYPRE_SetSpGemmUseVendor(0)` initially (hypre's own SpGemm has historically beaten
cuSPARSE for AMG setup); on a current release it is worth benchmarking both settings.
These parameter changes alter iteration counts slightly; that
is expected and acceptable as long as the solve converges to the same tolerance.

## Step 6 — Validate in the real application

1. Build a q2p1 application (e.g. `q2p1_devel`), run a case with
   `Pres@MGCrsSolverType = 7`, then `8`, on e.g. 4 and 8 ranks (`mpirun -np N`; remember
   rank 0 idles during the hypre solve, so N ranks = N-1 GPU clients).
2. Since multiple ranks share one A100, start **CUDA MPS**
   (`nvidia-cuda-mps-control -d`) for the multi-rank runs — without it, context switching
   between ranks serializes badly and will mask any speedup.
3. Correctness gates: identical (or near-identical) nonlinear convergence history vs. the
   CPU-hypre reference build; pressure field visually/numerically consistent; no growth in
   outer MG iterations.
4. Performance measurement: wall-clock of the coarse solve (the wrapper reports
   iterations; wrap the call site `source/src_quadLS/QuadSc_mg.f90:2097-2196` with timers
   if none print), plus `nsys profile` on one representative time step.

## Step 7 — Report back

Deliver:

- (a) the diff (CMake wiring, any managed-memory shim, AMG parameter changes),
- (b) build log confirmations and the exact hypre version/configure line used,
- (c) a small table — CPU-hypre vs GPU-hypre coarse-solve time and iteration counts at
  4/8 ranks for types 7 and 8,
- (d) any version-specific bugs or wrapper/API incompatibilities hit.

## Known limitations to state up front (not to fix in Stage 1)

- **Modest expected speedup.** The coarse-grid system is small; GPU AMG on it is
  latency-bound. Stage 1's purpose is to establish a working, validated GPU-hypre build
  and the managed-memory pattern — the payoff comes in Stage 2, when the *fine-level*
  pressure system is handed to the same machinery.
- **The vendored hypre stays at 2.25 for now.** Stage 1 links an external current hypre
  on the remote system; replacing the vendored `extern/libraries/hypre` tree with a
  current release (and revalidating all existing CPU builds) is a separate follow-up
  task, informed by which version Stage 1 validates.
