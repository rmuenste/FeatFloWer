# Reproducible Solver Dependencies via CMake FetchContent

The solver-library evaluation (Phases 0–3) introduced several dependencies
that were originally satisfied by hand-pointed site installs
(`/sfw/hypre/...`, `/sfw/intel/.../mkl`, the EL9 system UMFPACK). This page
documents the provider mechanism that makes those dependencies reproducible
from a clean clone while keeping every pre-existing configure path working
byte-for-byte.

Implementation: `cmake/modules/FFDependencies.cmake`, included from the root
`CMakeLists.txt` right after the option block.

## Pinned versions

| Dependency | Version | URL | SHA256 |
| --- | --- | --- | --- |
| hypre | 2.33.0 | `https://github.com/hypre-space/hypre/archive/refs/tags/v2.33.0.tar.gz` | `0f9103c34bce7a5dcbdb79a502720fc8aab4db9fd0146e0791cde7ec878f27da` |
| SuiteSparse | 7.12.2 | `https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.12.2.tar.gz` | `679412daa5f69af96d6976595c1ac64f252287a56e98cc4a8155d09cc7fd69e8` |

Both hashes were computed from an actual download of the tarball. MKL is
never fetched (proprietary binary distribution); MUMPS is deliberately left
on its existing `USE_EXTERNAL_MUMPS` find-path (Phase 2 verdict: not
competitive, so it was not worth a fetch path).

## Providers

Each fetchable dependency is resolved by exactly one *provider*:

| Provider | Meaning |
| --- | --- |
| `vendored` | The copy in `extern/libraries` (hypre submodule, UMFPACK 4 + AMD). |
| `external` | A site install named by the legacy `EXTERNAL_*_DIR` variables. |
| `system` | Whatever `find_package` locates on the default search paths; fails loudly if absent. |
| `fetch` | The pinned tarball above, downloaded and built as a sub-project. |
| `auto` | *(default)* `external` if the legacy `USE_EXTERNAL_*` option is ON, else `vendored` if the vendored copy exists, else `fetch`. |

Options:

| Option | Default | Purpose |
| --- | --- | --- |
| `FF_HYPRE_PROVIDER` | `auto` | `auto\|vendored\|external\|system\|fetch` |
| `FF_SUITESPARSE_PROVIDER` | `auto` | `auto\|vendored\|external\|system\|fetch` |
| `FF_FETCH_DEPENDENCIES` | `ON` | Set to `OFF` to make any configuration that *would* download fail with a fatal error instead. Use this on machines that must build strictly against site installs. |
| `FF_HYPRE_GPU_AWARE_MPI` | `OFF` | Drives `HYPRE_ENABLE_GPU_AWARE_MPI` in the fetched hypre; warns if a pre-built hypre does not report it. |
| `FF_MKL_USE_CMAKE_CONFIG` | `OFF` | Link MKL via `MKLConfig.cmake` (`MKL::MKL`) instead of the explicit static link group. Off by default — see [Why MKLConfig.cmake is not the default](#why-mklconfigcmake-is-not-the-default). |

Because `auto` prefers the legacy resolution, existing command lines keep
their exact previous behaviour:

- `-DUSE_HYPRE=ON -DUSE_EXTERNAL_HYPRE=ON -DEXTERNAL_HYPRE_DIR=...` → `external`
- `-DUSE_HYPRE=ON` alone → `vendored` (`extern/libraries/hypre/src`)
- `-DUSE_EXTERNAL_SUITESPARSE=ON [-DEXTERNAL_SUITESPARSE_DIR=...]` → `external`
- no SuiteSparse options → `vendored` (`extern/libraries/umfpack4` + `amd`)

Setting a provider explicitly is what opts a build into FetchContent, e.g.
`-DFF_HYPRE_PROVIDER=fetch -DFF_SUITESPARSE_PROVIDER=fetch`.

The `FetchContent_Declare` calls carry `FIND_PACKAGE_ARGS`, so with
`FF_*_PROVIDER=auto` on a machine without the vendored copies a system
install found through `CMAKE_PREFIX_PATH` still wins over a download.
`FF_*_PROVIDER=fetch` forces the tarball
(`FETCHCONTENT_TRY_FIND_PACKAGE_MODE` is not consulted because the fetch
branch calls `FetchContent_MakeAvailable` after a hard provider decision).

### Interface targets

Every provider populates one thin `INTERFACE` target:

- `ff::hypre` — also exported as `HYPRE_LIBRARIES` (consumed by
  `GenerateLinkerFlags.cmake`), with `HYPRE_STRUCTURES_INCLUDE_PATH` set for
  `GenerateIncludeFlags.cmake` so `include 'HYPREf.h'` resolves.
- `ff::umfpack` — also exported as `FF_UMFPACK_LINK_LIBS`.

`GenerateLinkerFlags.cmake` and `GenerateIncludeFlags.cmake` are therefore
unchanged. `ProjectFiles.cmake` now keys the `umf4_f77wrapper_port.c`
compilation on `FF_UMFPACK_WRAPPER_PORT` (true for the `external`, `system`
and `fetch` providers) instead of `USE_EXTERNAL_SUITESPARSE`.

## hypre build options

The fetched hypre is configured to match the manual reference build in
`stage1-work/hypre-2.33.0-build`:

```
HYPRE_ENABLE_MPI=ON  HYPRE_ENABLE_HYPRE_BLAS=ON  HYPRE_ENABLE_HYPRE_LAPACK=ON
HYPRE_ENABLE_BIGINT=OFF  HYPRE_ENABLE_MIXEDINT=OFF
HYPRE_BUILD_EXAMPLES=OFF  HYPRE_BUILD_TESTS=OFF
HYPRE_ENABLE_GPU_AWARE_MPI=${FF_HYPRE_GPU_AWARE_MPI}
```

With `-DUSE_HYPRE_CUDA=ON` the CUDA family is switched on as in the
reference CUDA build:

```
HYPRE_ENABLE_CUDA=ON  HYPRE_ENABLE_UNIFIED_MEMORY=ON
HYPRE_ENABLE_CUBLAS=ON  HYPRE_ENABLE_CURAND=ON
HYPRE_ENABLE_CUSPARSE=ON  HYPRE_ENABLE_CUSOLVER=ON
HYPRE_ENABLE_CUDA_STREAMS=ON
CMAKE_CUDA_ARCHITECTURES=80   (only if not already set)
```

hypre 2.33.0 keeps its CMake project in `src/`, hence `SOURCE_SUBDIR src`
in the declaration, and it provides a `HYPRE::HYPRE` alias for
FetchContent/`add_subdirectory` consumers.

`FF_HYPRE_GPU_AWARE_MPI` is the enabler for the GPU experiment described in
`solver_phase3_fine_level_pressure.md`; it only plumbs the hypre option
through, there is no FeatFloWer-side logic behind it.

## SuiteSparse build options

```
SUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;camd;colamd;ccolamd;cholmod;umfpack"
SUITESPARSE_USE_CUDA=OFF  SUITESPARSE_USE_PYTHON=OFF  SUITESPARSE_DEMOS=OFF
BUILD_SHARED_LIBS=OFF  BUILD_STATIC_LIBS=ON
```

Static libraries are deliberate: the applications then need no
`LD_LIBRARY_PATH` entry for the fetched UMFPACK on compute nodes.

Two guards are applied around the sub-configure, both of which are needed
because FeatFloWer would otherwise leak state into the sub-project:

1. `CMAKE_MODULE_PATH` is emptied. FeatFloWer ships its own
   `FindBLAS.cmake` / `FindLAPACK.cmake` / `FindMKL.cmake`, which
   SuiteSparse (and hypre) must not pick up instead of CMake's modules.
2. `BLAS_LIBRARIES` / `LAPACK_LIBRARIES` are unset. FeatFloWer sets them to
   its vendored netlib `blas`/`lapack` targets, and
   `SuiteSparseBLAS.cmake` uses a pre-set `BLAS_LIBRARIES` *as-is* — which
   would have silently linked UMFPACK against reference BLAS. With the
   variables cleared, SuiteSparse's own search runs and finds the system
   OpenBLAS (`/usr/lib64/libopenblas.so`, 32-bit integer interface).

Both are restored immediately afterwards.

### Migration note: UMFPACK 5.7.8 → 7.12.2

`USE_EXTERNAL_SUITESPARSE` previously pointed at the EL9 system UMFPACK
5.7.8. The fetch provider moves to SuiteSparse 7.12.2 (UMFPACK 6.x). The
`umf4*` Fortran entry points used by `source/UmfpackSolver.f90` come from
the vendored shim `extern/libraries/umfpack4/src/umf4_f77wrapper_port.c`,
which only calls the `umfpack_di_*` API — stable since UMFPACK 4 — so no
solver code changes are required. The Phase 0/1 coarse-solver benchmarks
(`crs-p-t2`, `crs-p-t4`) were blessed against 5.7.8; they have now been
re-run against the fetched 7.12.2 build — **the `crs-p-t2` result does not
carry over**, see
[Benchmark re-bless: UMFPACK 7.12.2 vs 5.7.8](#benchmark-re-bless-umfpack-7122-vs-578)
below.

## MKL PARDISO

MKL is never fetched, and `USE_MKL_PARDISO` / `MKL_PARDISO_DIR` / the
`MKL_PARDISO_AVAIL` define are unchanged. **The hand-assembled static link
group remains the default**:

```
-Wl,--start-group libmkl_gf_lp64.a libmkl_sequential.a libmkl_core.a -Wl,--end-group pthread m dl
```

The `MKLConfig.cmake` route (`MKL::MKL`, `MKL_LINK=static`,
`MKL_INTERFACE=lp64`, `MKL_THREADING=sequential`) is implemented but sits
behind the opt-in option `FF_MKL_USE_CMAKE_CONFIG` (default `OFF`).

### Why MKLConfig.cmake is not the default

`FF_MKL_PARDISO_LIBS` sits in the middle of `FF_DEFAULT_LIBS` /
`FF_APPLICATION_LIBS`, directly after the vendored BLAS/LAPACK archives. As
a list of raw flags and `.a` paths it is emitted at exactly that position in
the link line. As an imported target, CMake resolves it through the target
graph and emits it near the **end** of the link line instead.

On this toolchain that reordering is not cosmetic. The build mixes the
static `libgfortran.a` from the `/sfw` GCC 13.2.0 prefix with the system
`libgfortran.so.5`, so which `libgfortran.a` members the linker pulls in
depends on the undefined-symbol set at the moment `-lgfortran` is scanned.
With MKL moved to the end, the linked application ended up with
`_gfortran_st_open` / `_gfortran_st_read` bound to statically linked members
while `_gfortran_st_close` bound to the shared library — two Fortran unit
tables in one process. `CLOSE` then did not close the unit `OPEN`/`READ`
were using, the next sequential `READ` past EOF returned `iostat=5001`
("READ after EOF") instead of `-1`, and the parameter-parser loops that exit
only on `-1` spun forever: every rank hung deterministically at startup in
`GetPresParameters`.

This was proven by relinking the *identical* object files with only the MKL
group moved back to its original position — the hang disappears and
`_gfortran_st_close` becomes locally defined again.

Quick check on any new build:

```bash
nm applications/q2p1_fc_ext/q2p1_fc_ext | grep -E ' _gfortran_st_(open|read|close)$'
```

All three must be `T` (locally defined). A `U` on any of them means the
Fortran runtime is split and startup will hang.

If you enable `FF_MKL_USE_CMAKE_CONFIG=ON`, the resolved interface is still
checked for the GNU Fortran layer (`mkl_gf_*`) and a warning about the above
hazard is printed.

## Offline and cluster behaviour

**Downloads happen at configure time.** The login nodes have network access;
compute nodes do not. Never run `cmake` inside a Slurm job for a build that
still has to fetch — configure on the login node and only build/run in the
job, or pre-stage the sources.

### Pre-staging sources

Point FetchContent at an already-unpacked source tree; no network access,
no hash check:

```bash
cmake -S . -B <build> \
  -DFF_HYPRE_PROVIDER=fetch \
  -DFETCHCONTENT_SOURCE_DIR_HYPRE=$PWD/stage1-work/hypre-2.33.0-src \
  ...
```

The variable name is `FETCHCONTENT_SOURCE_DIR_<UCNAME>`, i.e.
`FETCHCONTENT_SOURCE_DIR_HYPRE` and `FETCHCONTENT_SOURCE_DIR_SUITESPARSE`.

### Sharing downloads between build trees

```bash
-DFETCHCONTENT_BASE_DIR=$PWD/stage1-work/_deps-cache
```

`FETCHCONTENT_BASE_DIR` holds `<dep>-src`, `<dep>-build` and
`<dep>-subbuild`. Sharing it avoids re-downloading *and* re-building, but
because the sub-project **build** trees live there too, only share it
between build trees that configure the dependency identically (same
compilers, same `USE_HYPRE_CUDA`, same `FF_HYPRE_GPU_AWARE_MPI`). For a CPU
tree and a CUDA tree, either use separate base dirs or share only the
sources via `FETCHCONTENT_SOURCE_DIR_*`.

### Forbidding downloads

```bash
-DFF_FETCH_DEPENDENCIES=OFF
```

Any provider resolution that would need a download then aborts with a fatal
error naming the dependency, instead of silently going to the network.

## CMake version

`FetchContent`'s `FIND_PACKAGE_ARGS` requires CMake ≥ 3.24. The project
minimum stays at **3.18** so that no existing non-fetch workflow breaks;
instead the fetch branch aborts with an explicit message on older CMake:

```
Fetching HYPRE needs CMake >= 3.24 (FetchContent FIND_PACKAGE_ARGS);
this is CMake <version>. Use a site install and set FF_HYPRE_PROVIDER=external/system.
```

`CMP0135` (extraction timestamps of downloaded archives) is set to `NEW`
inside `FFDependencies.cmake` so that re-running CMake does not needlessly
re-configure the sub-projects.

## Reference configure lines (warehouse/tardis)

```bash
export PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/bin:/sfw/openmpi/gcc13.2.x/4.1.6/ucx-threaded-noverbs/bin:$PATH
export LD_LIBRARY_PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/lib64:$LD_LIBRARY_PATH
gccpfx=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl
```

Fetched hypre + fetched SuiteSparse + MKL PARDISO (CPU):

```bash
cmake -S . -B stage1-work/build-fetchdeps-cpu -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON \
  -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort \
  -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=ON -DCGAL_DIR=/sfw/cgal/gcc13.2.0/lib64/cmake/CGAL \
  -DMPFR_INCLUDE_DIR=$gccpfx/include -DMPFR_LIBRARIES=$gccpfx/lib64/libmpfr.a \
  -DGMP_INCLUDE_DIR=$gccpfx/include -DGMP_LIBRARIES=$gccpfx/lib64/libgmp.a \
  -DGMPXX_INCLUDE_DIR=/usr/include -DGMPXX_LIBRARIES=/usr/lib64/libgmpxx.so \
  -DUSE_HYPRE=ON -DFF_HYPRE_PROVIDER=fetch \
  -DUSE_HYPRE_GPU_AMG=ON -DENABLE_HYPRE_TIMING=ON -DENABLE_SOLVER_TIMING=ON \
  -DFF_SUITESPARSE_PROVIDER=fetch \
  -DUSE_MKL_PARDISO=ON -DMKL_PARDISO_DIR=/sfw/intel/2024.0.1/mkl/latest
cmake --build stage1-work/build-fetchdeps-cpu --target q2p1_fc_ext q2p1_fac3d --parallel 12
```

Fetched hypre with CUDA and GPU-aware MPI. hypre's
`HYPRE_SetupCUDAToolkit.cmake` aborts with *"CUDA_PATH or CUDA_HOME not
set"* unless one of those **environment** variables points at the toolkit,
so export them next to `PATH`/`LD_LIBRARY_PATH`:

```bash
export PATH=/sfw/cuda/12.4/bin:$PATH
export LD_LIBRARY_PATH=/sfw/cuda/12.4/lib64:$LD_LIBRARY_PATH
export CUDA_HOME=/sfw/cuda/12.4
export CUDA_PATH=/sfw/cuda/12.4

cmake -S . -B stage1-work/build-fetchdeps-gpu -G Ninja \
  ... same CGAL/GMP/MPFR flags ... \
  -DUSE_HYPRE=ON -DFF_HYPRE_PROVIDER=fetch \
  -DUSE_HYPRE_CUDA=ON -DFF_HYPRE_GPU_AWARE_MPI=ON \
  -DFETCHCONTENT_SOURCE_DIR_HYPRE=$PWD/stage1-work/hypre-2.33.0-src \
  -DENABLE_HYPRE_TIMING=ON -DENABLE_SOLVER_TIMING=ON
```

The unchanged site/external line (Phase 2 evaluation build) still works
exactly as before:

```bash
cmake -S . -B <build> -G Ninja \
  ... same CGAL/GMP/MPFR flags ... \
  -DUSE_HYPRE=ON -DUSE_EXTERNAL_HYPRE=ON \
  -DEXTERNAL_HYPRE_DIR=/sfw/hypre/v2.33.0/gcc13.2.0-32bit \
  -DUSE_HYPRE_GPU_AMG=ON -DENABLE_HYPRE_TIMING=ON -DENABLE_SOLVER_TIMING=ON \
  -DUSE_EXTERNAL_SUITESPARSE=ON \
  -DUSE_MKL_PARDISO=ON -DMKL_PARDISO_DIR=/sfw/intel/2024.0.1/mkl/latest \
  -DUSE_MUMPS=ON -DUSE_EXTERNAL_MUMPS=ON -DEXTERNAL_MUMPS_DIR=/sfw/mumps/gcc13.2.0 \
  -DEXTERNAL_SCALAPACK_DIR=/sfw/scalapack/gcc13.2.0 \
  -DEXTERNAL_SCOTCH_DIR=/sfw/scotch/gcc13.2.0 \
  -DEXTERNAL_PARMETIS_DIR=/sfw/parmetis/gcc13.2.0-32bit
```

## Benchmark re-bless: UMFPACK 7.12.2 vs 5.7.8

The Phase 1 numbers for the gathered direct pressure coarse solvers
(`Pres@MGCrsSolverType` 2 and 4) were measured against the EL9 system
**UMFPACK 5.7.8**. The fetch provider builds **SuiteSparse 7.12.2
(UMFPACK 6.3.7)** instead, so the Phase 0/1 harness was re-run against
`stage1-work/build-fetchdeps-cpu` on the same node (tardis) with the same
rank counts, levels and step counts.

| Run | Result tree | Job | Build |
| --- | --- | --- | --- |
| 7.12.2 (fetch) | `stage1-work/results/rebless-umfpack712` | 136211 | `build-fetchdeps-cpu` |
| 5.7.8 same-day control | `stage1-work/results/rebless-control-suitesparse578` | 136213 | `build-phase1` |
| 5.7.8 original | `stage1-work/results/phase1-direct` | 135970 | `build-phase1` |
| vendored UMFPACK 4 | `stage1-work/results/phase0-baseline-v2` | 135968 | `build-matched-cpu` |

```bash
sbatch --export=ALL,FF_P0_RUN_ID=rebless-umfpack712,\
FF_P0_BUILD=$PWD/stage1-work/build-fetchdeps-cpu,FF_P0_SOLVER_TYPES="2 4 9" \
  tools/cluster_scripts/solver_baseline_phase0.slurm
```

Type 9 (MKL PARDISO) is carried along as a **cross-build control**: it uses
the same gather (`Setup_UMFPACK_CoarseSolver`) but never calls UMFPACK, so
its timings must be unchanged for the comparison to be apples-to-apples.
The 5.7.8 control run repeats the Phase 1 measurement on the same day to
separate build effects from run-to-run drift.

### Median time per pressure coarse call (ms)

Median over repetitions of the per-case median of `crs-p-t<N>`
(`total_median_s`), warmup 2 calls, all rows `status=ok`.

| Bench | np | Type | vendored UMF4 | 5.7.8 (Phase 1) | 5.7.8 (control) | 7.12.2 (fetch) | 7.12.2 / 5.7.8 | 7.12.2 / vendored |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cylinder | 4 | 2 | 0.3487 | 0.1849 | 0.1814 | 0.3474 | **1.92x** | 1.00x |
| cylinder | 4 | 4 | 0.5664 | 0.5607 | 0.5622 | 0.5616 | 1.00x | 0.99x |
| cylinder | 4 | 9 | – | 0.1166 | 0.1114 | 0.1139 | 1.02x | – |
| cylinder | 8 | 2 | 0.3529 | 0.1507 | 0.1506 | 0.3408 | **2.26x** | 0.97x |
| cylinder | 8 | 4 | 0.5787 | 0.5746 | 0.5753 | 0.5735 | 1.00x | 0.99x |
| cylinder | 8 | 9 | – | 0.1219 | 0.1235 | 0.1220 | 0.99x | – |
| FAC3D | 4 | 2 | 6.8654 | 4.6331 | 4.5070 | 5.1150 | **1.13x** | 0.75x |
| FAC3D | 4 | 4 | 21.4026 | 20.9646 | 20.9719 | 21.3832 | 1.02x | 1.00x |
| FAC3D | 4 | 9 | – | 1.5278 | 1.6308 | 1.8846 | 1.16x | – |

Ratios are against the same-day 5.7.8 control; > 1 means 7.12.2 is slower.
Per-repetition spread within a run is below 2 % everywhere, so the deltas
below are systematic, not noise. The 5.7.8 control reproduces the Phase 1
numbers to 0.95–1.07x, which confirms `phase1-direct` is still a valid
reference.

Stage-level records are unaffected: `mg-velocity`, `mg-pressure` and
`crs-u-t1` agree with Phase 1 within 0–4 % for every type.

### Forces

`BenchForce` trajectories were compared line by line against the Phase 1
runs (not just the last value):

- Types 4 and 9: **bitwise identical** on both benchmarks, all time steps.
- Type 2 on FAC3D: bitwise identical.
- Type 2 on cylinder: identical except the lift at one intermediate step
  (`7.6501786E-03` → `7.6501785E-03` at np4, `7.4825736E-03` →
  `7.4825735E-03` at np8) — one unit in the last printed digit, max
  relative deviation **1.3e-8**. The final-step forces are identical. The
  same one-digit difference appears against the vendored baseline, i.e.
  UMFPACK 6.3.7 rounds one solve marginally differently; 5.7.8 and
  UMFPACK 4 happen to agree there.

### Verdict

**Not a clean re-bless: the Phase 1 `crs-p-t2` speedup does not survive the
move to SuiteSparse 7.12.2.**

- `crs-p-t2` on cylinder is **1.9–2.3x slower** with UMFPACK 6.3.7 than
  with 5.7.8, landing exactly back at the vendored UMFPACK 4 level
  (1.00x / 0.97x of the vendored baseline). The Phase 1 headline
  "external UMFPACK t2 is 1.9–2.3x faster than vendored on cylinder" is
  therefore **specific to 5.7.8** and must not be quoted for a fetch-provider
  build.
- `crs-p-t2` on FAC3D keeps part of the win: 5.115 ms vs 6.865 ms vendored
  (**1.34x faster**), against 1.52x for 5.7.8. Note that the t9 control on
  FAC3D is itself 1.16x between the two builds — at the top of the accepted
  10–20 % control band — so the FAC3D t2 delta of 1.13x is at the level of
  the build-to-build drift on that case and should be read as "roughly
  unchanged, slightly worse", not as a precise 13 % regression. The cylinder
  control (0.99–1.02x) leaves no such ambiguity, and the cylinder regression
  is far outside it.
- `crs-p-t4` is unchanged everywhere (1.00–1.02x), as expected: it is
  dominated by the /16 extraction and Q1→P1 reconstruction, not by the
  factorization.
- Type 9 (PARDISO) remains the fastest gathered coarse solver on both
  benchmarks and is unaffected by the SuiteSparse version. With a fetched
  SuiteSparse it is now **3.0x faster than t2** on cylinder (0.114 vs
  0.347 ms) instead of 1.6x, which strengthens the Phase 1 recommendation
  to prefer type 9 where MKL is available.
- Correctness is not affected (see Forces above).

Cause: not the BLAS threading layer. A diagnostic re-run of the cylinder
np4/t2 case with `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` gives
0.348–0.352 ms, i.e. no change (job 136212). The fetched UMFPACK is built
`-O3 -DNDEBUG` and links the same system OpenBLAS as the site 5.7.8, and
the gather path is shared with types 4 and 9 which are unchanged — so the
extra ~0.17 ms per call sits inside `umfpack_di_solve` itself on this very
small coarse system. Cheapest mitigations, in order: use type 9 where MKL
is available, or keep `USE_EXTERNAL_SUITESPARSE=ON` pointing at 5.7.8 for
t2-heavy production runs.
