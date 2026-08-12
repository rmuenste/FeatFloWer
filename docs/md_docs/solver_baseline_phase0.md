# Solver Baseline Harness (Phase 0)

This note describes the Phase 0 measurement harness for the solver-library
evaluation. Its purpose is to make every later solver comparison trivial:
all coarse-solver backends and both MG stage solves report per-call
setup/solve/total timings in one machine-readable format, collected by one
summarizer over a fixed benchmark set.

The benchmark set is the `q2p1_fc_ext` flow-around-a-cylinder case and the
undeformed FAC3D case. SSE is deliberately excluded from this harness.
The harness explicitly sets `SimPar@ApplyFAC3DMeshDeformation = No` so FAC3D
mesh behavior cannot depend on the staged application's parameter-file default.

## Instrumentation

The CMake option `ENABLE_SOLVER_TIMING` (independent of `USE_HYPRE`) defines
`FF_SOLVER_TIMING` and activates `FF_TIMING` records. The module
`source/ff_timing.f90` prints one record per call in the same field layout
as the existing `HYPRE_TIMING` records:

```
FF_TIMING solver=<label> call=<n> setup_s=<s> solve_s=<s> total_s=<s> iterations=<n> rel_res=<r>
```

Times are the **minimum** across `MPI_COMM_WORLD` and are printed by the
master rank; iterations and residual are the maximum. The minimum is used
because ranks enter an instrumented routine at different moments — the idle
master arrives at the coarse solver first and waits there while the workers
finish smoothing — so the smallest per-rank duration measures the
straggler's actual pass through the solver, excluding synchronization wait.
`ff_report_timing` is collective; every instrumented site is reached by all
ranks.

Instrumented sites and labels:

| Label | Site | Meaning |
| --- | --- | --- |
| `mg-velocity` | `Solve_General_QuadScalar` around `MG_Solver` | One full momentum MG solve (per nonlinear iteration). `iterations` = MG cycles, `rel_res` = DefFinal/DefInitial. |
| `mg-pressure` | `Solve_General_LinScalar` around `MG_Solver` | One full pressure MG solve. Same fields. |
| `crs-u-t<N>` | `mgCoarseGridSolver_U` (whole routine) | One velocity coarse solve with `Velo@MGCrsSolverType = N`. For MUMPS (N=5), `setup_s` covers init+matrix setup and `solve_s` the factorize/solve/gather. |
| `crs-p-t<N>` | `mgCoarseGridSolver_P` (whole routine) | One pressure coarse solve with `Pres@MGCrsSolverType = N`, including gather/scatter and, for types 3/4/8, the /16 extraction and Q1-to-P1 reconstruction. |

For coarse types 7/8 the inner `HYPRE_TIMING` records continue to appear;
the `crs-p-t7`/`crs-p-t8` wrapper records additionally include the
marshalling and reconstruction cost around the Hypre call, which the inner
records do not capture. `rel_res` is 0 where no residual is available
(direct solvers, wrapper records).

Build example (matched CPU build):

```bash
cmake -S . -B stage1-work/build-matched-cpu -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON \
  -DUSE_HYPRE=ON -DUSE_EXTERNAL_HYPRE=ON \
  -DEXTERNAL_HYPRE_DIR=/sfw/hypre/v2.33.0/gcc13.2.0-32bit \
  -DUSE_HYPRE_GPU_AMG=ON -DENABLE_HYPRE_TIMING=ON \
  -DENABLE_SOLVER_TIMING=ON
cmake --build stage1-work/build-matched-cpu --target q2p1_fc_ext q2p1_fac3d --parallel 8
```

`ENABLE_SOLVER_TIMING=OFF` (the default) compiles all instrumentation away;
the default build is unchanged.

## Job harness

`tools/cluster_scripts/solver_baseline_phase0.slurm` stages and runs the CPU
baseline for both benchmarks. It follows the Stage 1 script conventions
(independent case directories, `case.env`, per-case `run.log`), but is
CPU-only: no CUDA module, no MPS, no `nvidia-smi` requirements. A failing
case records its exit code in `case.env` and does not abort the job.

```bash
sbatch --export=ALL,FF_P0_RUN_ID=phase0-baseline \
  tools/cluster_scripts/solver_baseline_phase0.slurm
```

Defaults and overrides (`FF_P0_*` environment variables):

| Variable | Default | Meaning |
| --- | --- | --- |
| `FF_P0_BUILD` | `stage1-work/build-matched-cpu` | Build directory with instrumented executables. |
| `FF_P0_APPS` | `cylinder fac3d` | Benchmark selection. |
| `FF_P0_SOLVER_TYPES` | `1 2 3 4 7 8` | `Pres@MGCrsSolverType` values to sweep. Type 5 requires a MUMPS build. |
| `FF_P0_REPEATS` | `3` | Repetitions per matrix entry. |
| `FF_P0_WARMUP_CALLS` | `2` | Per-solver calls discarded from timing medians (correctness checks use all calls). |
| `FF_P0_STEPS_CYLINDER` / `FF_P0_STEPS_FAC3D` | `10` / `4` | Time steps per run. FAC3D needs more steps than `FF_P0_WARMUP_CALLS`, since `mg-pressure`/`mg-velocity` fire only a few times per step. |
| `FF_P0_RANKS_CYLINDER` / `FF_P0_RANKS_FAC3D` | `4 8` / `4` | MPI rank counts. |
| `FF_P0_FAC3D_LEVEL` | `4` | FAC3D `MaxMeshLevel`. |
| `FF_P0_KEEP_OUTPUT` | `0` | Keep `_vtk`/`_dump` files per case. |

Case directories are named `<app>-cpu-np<ranks>-type<type>-rep<rep>` under
`stage1-work/results/<run-id>/runs/`.

## Summarizer

`tools/cluster_scripts/summarize_hypre_gpu_stage1.py` now parses both
`HYPRE_TIMING` and `FF_TIMING` records and emits one TSV row per
(case, solver label) with median setup/solve/total times and iteration
counts. New columns: `app` and `solver`. Old Stage 1 result trees (case
names without the app prefix, Hypre records only) still parse.

Correctness gates: application exit code and residual finiteness apply to
every row; the residual thresholds (`1e-5` for `gmres-amg`, `1e-7` for
`amg`/`pcg-amg`) apply only to Hypre records, whose solver tolerances are
known. Rows are labeled `ok`, `residual_failed`, `nonfinite`, `app_failed`,
or `no_samples`.

## What "baseline" means for later phases

Phase 1+ candidates (external SuiteSparse, PARDISO, upgraded MUMPS, new
Hypre configurations) are accepted or rejected against the medians and
correctness gates produced by this harness on the same node class, rank
counts, and step counts. The relevant comparisons are:

- `crs-p-t2`/`crs-p-t4` medians: the current UMFPACK 4 coarse solve that a
  modern SuiteSparse or PARDISO must beat (or match) in Phase 1.
- `crs-p-t7`/`crs-p-t8` vs the inner `gmres-amg`/`amg` records: the
  marshalling overhead a fine-level Hypre path (Phase 3) must amortize.
- `mg-velocity` and `mg-pressure`: the stage-level costs that bound what any
  coarse-solver improvement can be worth end-to-end, and the reference for a
  Phase 4 fine-level momentum solve.
