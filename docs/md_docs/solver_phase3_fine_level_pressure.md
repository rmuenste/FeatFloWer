# Phase 3: Fine-Level Pressure Solve via Hypre GMRES+BoomerAMG

Phase 3 of the solver-library evaluation adds a pressure solver path that
hands the **finest-level** (NLMAX) pressure system to Hypre, bypassing the
geometric multigrid entirely. It targets the FAC3D-class workload where the
Phase 0 baseline showed the pressure MG stage dominating the time step
(25.4 s per solve vs 5.7 s for the whole momentum stage).

## Usage

```
Pres@MGCrsSolverType = 10
```

Type 10 is not a coarse-grid solver: `Solve_General_LinScalar` dispatches to
`Solve_LinScalar_HYPRE_Fine` instead of `MG_Solver`. Requires a Hypre build
(`USE_HYPRE=ON`); selecting it otherwise stops with a message. The
harness sweeps it like any other type (`FF_P0_SOLVER_TYPES="9 10"`).

## Implementation

- `Setup_HYPRE_CoarseLevel_Full` (`QuadSc_solver_hypre.f90`) gained an
  optional `iSetupLev` argument; the validated CSR+halo marshalling now runs
  at any level. The global-numbering handshake
  (`GetMyHYPRENumberingLimits`) ignores the master's own element count, so
  the master may remain at its reduced level.
- The fine matrix is marshalled **lazily at the first pressure solve**
  (after the boundary operators are applied), then reused — the same
  matrix-reuse semantics as coarse types 7/8, valid while
  `myMatrixRenewal%C` keeps the pressure matrix constant.
- The solve itself is the unchanged Stage-1-validated `myHypreGMRES_Solve`
  (tol 1e-5, maxit 80, one BoomerAMG V-cycle as preconditioner, PMIS
  coarsening / ext+i interpolation / symmetric hybrid relaxation in
  `HYPRE_GPU_AMG` builds).
- Like types 7/8, type 10 turns the master's fine levels off
  (`NLMAX=NLMIN`, `bMasterTurnedON=.FALSE.` in `QuadSc_initialization.f90`).
- Limitation: the singular no-outflow fix is coarse-level specific and is
  not applied at the fine level; type 10 is for configurations with an
  outflow (both harness benchmarks qualify).

## Results (jobs 136107 smoke, 136108 full; results `phase3-finelevel`)

`mg-pressure` stage medians per solve, type 9 (geometric MG + PARDISO
coarse) vs type 10 (fine-level Hypre):

| Case | Geometric MG | Hypre fine level | GMRES iters |
| --- | ---: | ---: | ---: |
| cylinder np4 | 10.7 ms | 20.9 ms (0.5×) | 19 |
| cylinder np8 | 5.2 ms | 12.9 ms (0.4×) | 22 |
| **FAC3D np4** | **25.46 s** | **14.49 s (1.76× faster)** | 17 |

All 72 correctness rows `ok`. FAC3D forces identical between paths;
cylinder forces agree to ~1e-6 relative in drag (tolerance-level difference
from the 1e-5 GMRES stopping criterion — the geometric MG and GMRES stop on
different criteria).

Conclusions:

1. **Fine-level AMG is the right tool exactly where the plan predicted**:
   large fine-level pressure systems (FAC3D-class). For small systems the
   geometric MG remains superior — keep type 9 (or 2) there.
2. This is the true Stage 2 GPU workload: ~17 GMRES+AMG iterations over
   the distributed fine system (~1.2M rows at FAC3D level 4) per pressure
   solve, instead of the 204-row coarse solves that made Stage 1
   latency-bound. Re-running the `USE_HYPRE_CUDA` backend against type 10
   is the natural next experiment.
3. Refinement candidates, in order: PCG instead of GMRES (the system is
   near-SPD; `myHyprePCG_Solve` exists), a configurable tolerance (the
   1e-5 default is hard-coded in `HypreSolver.f90`), re-marshalling values
   when `myMatrixRenewal%C` rebuilds the matrix (needed before using
   type 10 in applications with time-varying pressure operators).

## GPU rerun (job 136117, build `stage1-work/build-gpu-phase3`)

The CUDA Hypre 2.33 build was rebuilt with the current sources
(`USE_HYPRE_CUDA=ON` + `ENABLE_SOLVER_TIMING=ON`) and run with type 10 on
FAC3D (np4, 3 workers sharing the A100 through MPS; the harness gained
`FF_P0_BACKEND=gpu` for CUDA module + MPS + `nvidia-smi dmon` handling).

| Backend | mg-pressure per solve | GMRES iters | per iteration |
| --- | ---: | ---: | ---: |
| CPU type 9 (geometric MG) | 25.46 s | — | — |
| CPU type 10 | 14.49 s | 17 | 0.85 s |
| GPU type 10 | 29.40 s | 17.5 | 1.66 s |

Numerically the GPU path is flawless (identical iteration counts and
residuals to CPU). Performance-wise it is 2× slower than CPU — but for a
different reason than Stage 1: average SM utilization is now 21% with
bursts to 100% (Stage 1: active in 22 of 1,933 samples), so the device
finally has real work. The bottleneck has moved to communication: the
Open MPI build is not CUDA-aware, so every halo exchange in every
BoomerAMG level of every iteration round-trips device→host→MPI→host→device
across three ranks sharing one GPU.

Next GPU steps, in order of expected impact:

1. **CUDA-aware Open MPI + `HYPRE_ENABLE_GPU_AWARE_MPI=ON`** — the site
   module tree exposes an `openmpi/options/cuda/...` switch; Hypre must be
   rebuilt with GPU-aware MPI enabled. This attacks the identified
   bottleneck directly.
2. Fewer ranks per GPU. The clean 1-worker isolation experiment is blocked
   by pre-existing single-worker bugs: two invalid `0(...)` runtime formats
   in `get_pid.f90` (fixed) and a remaining segfault in the parallel
   structure setup at `subnodes=1` (not pursued).
3. Larger mesh levels (level 5) — more work per transfer.
