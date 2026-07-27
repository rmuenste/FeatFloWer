# Phase 1: Modern Gathered Direct Solvers for the Pressure Coarse System

Phase 1 of the solver-library evaluation (see `solver_baseline_phase0.md` for
the harness and baseline) upgrades the gathered direct solvers behind the
pressure coarse dispatch:

1. **External SuiteSparse UMFPACK** replaces the vendored 2005-era UMFPACK 4
   behind the existing solver types 2 and 4, with no solver-code changes.
2. **MKL PARDISO** becomes the new pressure coarse solver type 9
   (`Pres@MGCrsSolverType = 9`), solving the same gathered full coarse matrix
   as type 2.

## External SuiteSparse (types 2/4)

CMake options:

| Option | Purpose |
| --- | --- |
| `USE_EXTERNAL_SUITESPARSE` | Link an external SuiteSparse UMFPACK instead of building the vendored `extern/libraries/umfpack4` + `amd`. |
| `EXTERNAL_SUITESPARSE_DIR` | Install prefix; empty uses system paths (`/usr/include/suitesparse`, `libumfpack`). |

The Fortran code keeps calling the `umf4*` interface. The handle-table
wrapper `extern/libraries/umfpack4/src/umf4_f77wrapper_port.c` is compiled
against the external `umfpack.h` (it only calls the standard `umfpack_di_*`
API, which is unchanged since UMFPACK 4), so `INTEGER*4` handles remain safe.
On the EL9 nodes the system package provides UMFPACK 5.7.8; pointing
`EXTERNAL_SUITESPARSE_DIR` at a self-built SuiteSparse 7.x works the same
way.

## MKL PARDISO (type 9)

CMake options:

| Option | Purpose |
| --- | --- |
| `USE_MKL_PARDISO` | Compile `source/PardisoSolver.f90`, define `MKL_PARDISO_AVAIL`, link static sequential MKL. |
| `MKL_PARDISO_DIR` | MKL install root (defaults to `$MKLROOT`; on the warehouse nodes `/sfw/intel/2024.0.1/mkl/latest`). |

Type 9 mirrors type 2: the master gathers the full coarse pressure matrix
into `UMF_CMat`/`UMF_lMat` (`Setup_UMFPACK_CoarseSolver`), PARDISO factorizes
it once (phase 12, METIS ordering, weighted matching, scaling), and each
coarse solve is a phase-33 forward/backward substitution with iterative
refinement. The factorization is reused across time steps exactly like the
UMFPACK one and is released under the same `myMatrixRenewal%C >= 2`
condition. Sequential MKL is linked deliberately: the master rank is bound
to one core during the coarse solve.

Restrictions: type 9 requires `Pres@MGMedLev == Pres@MGMinLev` (the nested
coarse-MG path used by types 1–4 when the levels differ is not wired for
it); the run stops with a clear message otherwise. Selecting type 9 in a
build without `USE_MKL_PARDISO` stops with "MKL PARDISO is not available!".

## Build (warehouse/tardis nodes)

```bash
export PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/bin:/sfw/openmpi/gcc13.2.x/4.1.6/ucx-threaded-noverbs/bin:$PATH
gccpfx=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl
cmake -S . -B stage1-work/build-phase1 -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON \
  -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort \
  -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=ON -DCGAL_DIR=/sfw/cgal/gcc13.2.0/lib64/cmake/CGAL \
  -DMPFR_INCLUDE_DIR=$gccpfx/include -DMPFR_LIBRARIES=$gccpfx/lib64/libmpfr.a \
  -DGMP_INCLUDE_DIR=$gccpfx/include -DGMP_LIBRARIES=$gccpfx/lib64/libgmp.a \
  -DGMPXX_INCLUDE_DIR=/usr/include -DGMPXX_LIBRARIES=/usr/lib64/libgmpxx.so \
  -DUSE_HYPRE=ON -DUSE_EXTERNAL_HYPRE=ON \
  -DEXTERNAL_HYPRE_DIR=/sfw/hypre/v2.33.0/gcc13.2.0-32bit \
  -DUSE_HYPRE_GPU_AMG=ON -DENABLE_HYPRE_TIMING=ON -DENABLE_SOLVER_TIMING=ON \
  -DUSE_EXTERNAL_SUITESPARSE=ON \
  -DUSE_MKL_PARDISO=ON -DMKL_PARDISO_DIR=/sfw/intel/2024.0.1/mkl/latest
cmake --build stage1-work/build-phase1 --target q2p1_fc_ext q2p1_fac3d --parallel 12
```

(The GMP/MPFR/CGAL variables replicate what the `gcc` and `cgal` environment
modules provide on the compute nodes.)

## Benchmark and acceptance gate

Run the Phase 0 harness against the Phase 1 build, adding type 9:

```bash
sbatch --export=ALL,FF_P0_RUN_ID=phase1-direct,FF_P0_BUILD=$PWD/stage1-work/build-phase1,FF_P0_SOLVER_TYPES="1 2 3 4 7 8 9" \
  tools/cluster_scripts/solver_baseline_phase0.slurm
```

Gate (against `stage1-work/results/phase0-baseline-v2/summary.tsv`):

- All statuses `ok`; `last_bench_force` identical to the printed precision
  for matched cases.
- `crs-p-t2`/`crs-p-t4` medians with external SuiteSparse at or below the
  vendored-UMFPACK baseline (cylinder np4: 0.35 ms / 0.57 ms; FAC3D: 6.9 ms /
  21.4 ms).
- `crs-p-t9` compared against `crs-p-t2` decides whether PARDISO becomes a
  recommended option.
