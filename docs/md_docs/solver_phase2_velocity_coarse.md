# Phase 2: Velocity Coarse Solver Options and External MUMPS

Phase 2 of the solver-library evaluation targets the velocity (Burgers)
coarse solve, which the Phase 0 baseline identified as ~95% of the momentum
stage on the cylinder benchmark (type-1 Jacobi, 50 iterations per call,
32.6 ms of the 34 ms `mg-velocity` stage at np4).

## External MUMPS (removes the Intel-compiler requirement)

The vendored MUMPS 5.0.1 required the Intel compiler and MKL. The new CMake
path links a pre-installed GCC-compatible MUMPS instead, making velocity and
pressure coarse solver type 5 available in GCC builds:

| Option | Purpose |
| --- | --- |
| `USE_EXTERNAL_MUMPS` | (requires `USE_MUMPS=ON`) Link an external MUMPS instead of building the vendored 5.0.1 with MKL. |
| `EXTERNAL_MUMPS_DIR` | MUMPS install prefix (site: `/sfw/mumps/gcc13.2.0`, version 5.7.3). |
| `EXTERNAL_SCALAPACK_DIR` | Matching ScaLAPACK prefix (site: `/sfw/scalapack/gcc13.2.0`). |
| `EXTERNAL_SCOTCH_DIR` | Matching SCOTCH/esmumps prefix (site: `/sfw/scotch/gcc13.2.0`). |

The link list is `dmumps mumps_common pord` + ScaLAPACK + SCOTCH variants +
the vendored `metis` target + OpenBLAS. `MumpsSolver.f90` compiles unchanged
against the 5.7.3 `dmumps_struc.h` (all used fields verified present).
`MUMPS_STRUCTURES_INCLUDE_PATH` automatically points at the external
headers.

Build: add to the Phase 1 configuration (see
`solver_phase1_direct_solvers.md`):

```bash
  -DUSE_MUMPS=ON -DUSE_EXTERNAL_MUMPS=ON \
  -DEXTERNAL_MUMPS_DIR=/sfw/mumps/gcc13.2.0 \
  -DEXTERNAL_SCALAPACK_DIR=/sfw/scalapack/gcc13.2.0 \
  -DEXTERNAL_SCOTCH_DIR=/sfw/scotch/gcc13.2.0
```

The evaluation build tree is `stage1-work/build-phase2`.

## Velocity coarse sweep in the harness

`tools/cluster_scripts/solver_baseline_phase0.slurm` gained
`FF_P0_VELO_TYPES` (default `1`), which sets `Velo@MGCrsSolverType` per
case. Case directories are now named
`<app>-cpu-np<ranks>-type<pres>-vt<velo>-rep<rep>`; the summarizer emits a
`velo_type` column (old case names parse as `velo_type = 1`).

Velocity coarse types: 1 = Jacobi/block-SSOR (baseline), 2 = BiCGStab,
5 = MUMPS (distributed, full analyze+factorize+solve per call because the
momentum matrix changes every nonlinear iteration).

Phase 2 benchmark (pressure fixed at the Phase 1 winner, type 9):

```bash
sbatch --export=ALL,FF_P0_RUN_ID=phase2-velocity,FF_P0_BUILD=$PWD/stage1-work/build-phase2,FF_P0_SOLVER_TYPES=9,FF_P0_VELO_TYPES="1 2 5" \
  tools/cluster_scripts/solver_baseline_phase0.slurm
```

## Gate

Against `crs-u-t1` and `mg-velocity` in the Phase 0 baseline
(`phase0-baseline-v2`): a velocity coarse option is accepted if all
correctness gates pass (statuses `ok`, forces unchanged) and it reduces the
`mg-velocity` stage median. The Hypre velocity marshalling (types 7/8 for
the velocity dispatch) is deferred until these existing options are
measured; if MUMPS or BiCGStab already removes the coarse bottleneck, the
marshalling work may not be justified at current problem sizes.

## Results (jobs 135971, 135972, 136026, 136027)

**Level constraint discovered and guarded.** The global Q2 numbering used
by the velocity MUMPS data path (`Create_GlobalNumbering`) is built once at
`Pres@MGMedLev`, while the velocity coarse solve runs at `Velo@MGMinLev`.
With the standard configs (velocity 2, pressure 1) the mismatch killed
every velocity-MUMPS run in an opaque `MPI_ERR_TRUNCATE`;
`mgCoarseGridSolver_U` now detects it collectively and stops with an
actionable message. Even `applications/mumps_test` ships with the
mismatched levels, so this path cannot have worked as configured.

Per-call medians with the velocity coarse level aligned to 1
(`FF_P0_VELO_MINLEV=1`), pressure type 9:

| Case | mg-velocity, L2+Jacobi | mg-velocity, L1+Jacobi | mg-velocity, L1+MUMPS |
| --- | ---: | ---: | ---: |
| cylinder np4 | 33.9 ms | **13.5 ms** | 35.8 ms |
| cylinder np8 | 15.6 ms | **6.5 ms** | 28.2 ms |
| FAC3D np4 | 5.71 s | 5.87 s | 6.58 s |

Conclusions:

1. **The dominant win is the level, not the library**: aligning
   `Velo@MGMinLev/MGMedLev = 1` cuts the cylinder momentum stage ~2.5×
   with the existing Jacobi coarse solver (coarse solve 32.6 → 4.5 ms).
   On FAC3D, whose momentum stage is fine-smoothing dominated, the effect
   is neutral (~3% slower).
2. **Velocity MUMPS (type 5) is fixed and correct** — forces bitwise
   identical to the Jacobi path at the same level — but not competitive:
   the momentum matrix changes every nonlinear iteration, so MUMPS
   re-runs analysis+factorization per call. Keep it as the robust direct
   fallback (e.g. stiff non-Newtonian coupled systems).
3. **BiCGStab (type 2) is ~6× slower** than Jacobi as configured — not a
   contender.
4. Force gate: at a fixed level, all coarse solvers agree bitwise;
   changing the level shifts the lift coefficient by ~5e-8, comparable to
   the accepted Stage 1 CPU/GPU differences (hierarchy change, different
   iteration path).
5. **Hypre-for-velocity marshalling is not justified** at these problem
   sizes: after the level fix the velocity coarse solve is 4.5 ms of a
   13.5 ms stage; the remaining cost is fine-level smoothing (a Phase 4
   PETSc topic, not a coarse-solver one).
