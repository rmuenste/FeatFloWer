# GPU pressure-Poisson test on one A100 80 GB (tardis) — report

Executed 2026-10-02 per `gpu-poisson-a100/README.md` §9 (stages 0–4), route 1 (CuPy + PyAMG).
Harness `gpu_amg_cg.py`; Slurm scripts `jobs/stage{0,1,2,2b,3,3c}.sbatch`; raw results `results/gpu_*.json`,
`results/stage*_*.json`, logs `results/stage*_<jobid>.log`, device-memory traces `results/*_nvsmi_*.csv`.

## Headline

**One pressure solve at the D6.4 size (160 × 160 × 480, 49.15 M unknowns, 5.21 G non-zeros) on one A100 80 GB:
2.4 s with AMG-CG (11 iterations, V-cycle 0.19 s), 24.4 s with Jacobi-CG (520 iterations); peak device
memory 54.9 GB in symmetric storage (31.9 GB fine matrix + 17.3 GB AMG hierarchy + vectors). The full size
fits.** Relative error vs x_true 8.4e-9 (AMG-CG) / 2.9e-7 (Jacobi-CG) at tolerance 1e-8 on the relative residual.

Caveats (§ "Caveats" below): synthetic matrix values of the same construction (Q2/P1-disc, C = Bᵀ Ml⁻¹ B, uniform
grid, no FBM, box walls, constant mode pinned), AMG setup done once on the host CPU (PyAMG, single thread,
248 s at full size, not included in the per-solve time), fine operator and CG in f64, hierarchy in f32.

## Per-stage status

| stage | what | status | key numbers |
|---|---|---|---|
| 0 | environment on tardis | done | A100 80GB PCIe (79.2 GiB usable), driver 580.159.04, toolkit cuda/12.6, python 3.13.5, cupy-cuda12x 14.2.0, pyamg 5.3.0, numpy 2.3.1, scipy 1.16.0; node-local `/scratch` 1.7 TB (1.5 TB free), 502 GB host RAM, 16 cores/32 threads — `results/environment.txt` |
| 1 | correctness + baseline, shipped instances | done | 40x40x120: GPU Jacobi-CG 141 it on stored b (err 7.2e-8), 274 it on the cpu_reference rough RHS = CPU count exactly (20x20x60: 139 = 139); symmetric matvec = full matvec to 1.8e-16; AMG-CG 11–30 it by variant — `results/gpu_p1_20x20x60.json`, `results/gpu_p1_40x40x120.json`, `results/stage1_*.json` |
| 2 | half size 128x128x384, full storage | done, fits | peak device 44.1 GB (31.5 GB matrix); Jacobi-CG 423 it / 11.6 s; AMG-CG 11 it / 1.41 s; host setup 130 s / 140 GB — `results/gpu_p1_128x128x384.json` (+ `stage2_p1_128x128x384_{sym,full_kernel}.json`) |
| 3 | full size 160x160x480, symmetric storage | done, fits (no ¾ fallback needed) | peak device 54.9 GB; Jacobi-CG 520 it / 24.4 s; AMG-CG 11 it / 2.41 s; host setup 248 s / 291 GB — `results/gpu_p1_160x160x480.json`; full-storage attempts in `stage3_p1_160x160x480_full_kernel*.json` |
| 4 | report | done | this file |

## Stage 0 — environment and submission notes

`module load python/3.13.5 cuda/12.6`; `python3 -m pip install --user cupy-cuda12x pyamg` (the bare `pip` on
PATH is the system Python-3.9 pip and installs cp39 wheels; use `python3 -m pip`). Slurm: `--nodelist=tardis`
is rejected by the site's lua job-submit plugin ("Requested node configuration is not available", also for
the other GPU nodes); `--constraint=ampere --gres=gpu:1 --nodes=1 --partition=med` selects tardis (the only
ampere node) and was used for all jobs; whole node (`-c 32 --mem=0`) for stages 2–3. All GPU work ran inside
Slurm jobs on tardis (147429, 147430, 147432, 147433, 147436, 147441); the queue was empty at every submission.

## Stage 1 — correctness and baseline on the shipped instances

| instance | nu | nnz | Jacobi-CG, stored b (GPU) | Jacobi-CG, rough RHS GPU / CPU reference | error vs x_true | CPU ref time | GPU time |
|---|---|---|---|---|---|---|---|
| 20x20x60 | 96 000 | 9.49 M | 72 it | 139 / 139 | 4.2e-8 | 1.10 s | 0.06 s |
| 40x40x120 | 768 000 | 79.0 M | 141 it | 274 / 274 | 7.2e-8 | 19.7 s | 0.28 s |

The rough RHS is exactly the one of `cpu_reference.py` (rng(1) normal, ones-mean removed); the stored b is
smooth, hence the lower count. 4x4 block-Jacobi-CG gives the same counts (72 / 140), as in the CPU reference.
The symmetric-storage matvec (upper triangle U, y = U x + (U − D)ᵀ x, warp-per-row CUDA kernel with f64
atomics) agrees with the full cuSPARSE matvec to 1.8e-16 relative (brief note 3 asked for 1e-12).

AMG variants on 40x40x120 (PyAMG setup on the host; device V-cycle: damped Jacobi ω = 4/3·1/ρ(D⁻¹A), 2 pre /
2 post, coarse levels and transfers in f32, dense coarse inverse on the device; flexible PCG in f64, tol 1e-8):

| variant | description | levels / sizes | op. complexity | AMG-CG it | solve s | error vs x_true |
|---|---|---|---|---|---|---|
| sa_csr_spec | README spec: SA on dof-level CSR, B = ones, strength symmetric, aggregate standard, smooth jacobi, max_coarse 1000 | 3: 768k / 7840 / 133 | 1.005 | 30 | 0.128 | 9.3e-8 |
| sa_csr_mean | as spec with B = e (the true constant-pressure mode: 1 on cell means, 0 on slopes) | 3 | 1.003 | 28 | 0.120 | 3.7e-8 |
| sa_bsr_mean | SA on the 4x4 cell blocks (BSR), B = e | 3 | 1.002 | 28 | 0.119 | 3.7e-8 |
| **sa_bsr_lin** | SA on the 4x4 cell blocks, B = [1, x, y, z] in the P1-disc basis | 3: 768k / 31 360 / 1336 | 1.040 | **11** | **0.046** | 9.7e-9 |
| geo_lin | predefined 2x2x2 cell aggregation, B = linears, unsmoothed P (nested P1-disc coarse spaces) | 5 | 1.136 | 49 | 0.248 | 1.7e-7 |
| geo_lin_jac | geo_lin with Jacobi-smoothed P | 5 | 1.673 | 6 | 0.031 | 2.0e-9 |

Decisions taken from this: the operator's near-nullspace is the cell-mean indicator, not `ones`, and the
piecewise-linear candidate set on the cell blocks (sa_bsr_lin: 11 iterations, complexity 1.04) is the best
preconditioner per byte — it became the primary variant for stages 2–3. geo_lin_jac converges fastest but its
complexity 1.67 (≈ 2.5 × the fine matrix in f32 at full size) cannot fit; geo_lin's count grows with the
number of levels (27 → 49 → 105 → 169). pyamg's own CG solve on the spec hierarchy needs 11–13 iterations with
its block-Gauss-Seidel smoothers vs our 23–30 with Jacobi. ω = 1/ρ instead of 4/3/ρ costs 2–8 iterations
(`stage1_p1_40x40x120_omega1.json`). Matvec on 40x40x120: cuSPARSE 0.70 ms, own CSR warp kernel 0.73 ms,
symmetric kernel 0.57 ms.

## Stage 2 — half size 128x128x384, full storage

Instance built slab-wise on tardis by the harness (`--build`; identical construction to `synth_p1_pressure.py`,
verified bitwise on both shipped instances with `--verify-build`): 25 165 824 unknowns, 2 610 809 224 nnz
(103.7/row), 92 s with 8 workers / 62 GB host RAM, written to node-local
`/scratch/rmuenste/gpu-poisson/p1_128x128x384.npz` (31.9 GB, 41 s). Nothing beyond the two shipped instances
was written to the warehouse volume.

| storage / matvec | fine operator on device | matvec | peak device memory | Jacobi-CG | AMG-CG sa_bsr_lin | V-cycle |
|---|---|---|---|---|---|---|
| full, cuSPARSE (3 row slabs ≤ 2³⁰ nnz) | 31.5 GB | 23.6 ms (1334 GB/s) | **44.1 GB** (nvidia-smi 42.1 GB) | 423 it, 11.6 s, err 2.5e-7 | **11 it, 1.41 s** (128 ms/it), err 1.1e-8 | 110 ms |
| full, own CSR warp kernel | 31.5 GB | 24.9 ms | 43.5 GB | – | 11 it, 1.47 s | 115 ms |
| symmetric, atomics kernel | 16.0 GB | 20.7 ms | 27.9 GB (nvidia-smi 30.6 GB) | 423 it, 10.4 s | 11 it, 1.27 s | 99 ms |

sa_bsr_lin hierarchy: 4 levels (25.2 M / 946 688 / 38 700 / 1436), operator complexity 1.038, 8.8 GB on the
device, PyAMG host setup 129 s, host RSS 140–172 GB. geo_lin: 6 levels, complexity 1.141, 105 it / 14.7 s.
Block-Jacobi-CG 422 it / 15.5 s. Memory agrees with the README's estimate (~45 GB).

Two facts recorded rather than hidden:
1. cuSPARSE `SpMV` returned `CUSPARSE_STATUS_INTERNAL_ERROR` on a row slab of 2.1 G nnz with int32 indices
   (job 147432); slabs were capped at 2³⁰ nnz and the run repeated (job 147433). Row pointers are int64
   throughout and column indices int32, so the matrix itself has no 2³¹ limit; CuPy's own CSR would force
   int64 column indices above 2³¹ nnz (+20–40 GB), which the slab/kernel paths avoid.
2. The README-spec variant (SA on the dof-level CSR) cannot run above 2³¹ nnz: scipy upcasts the indices to
   int64 and pyamg 5.3.0's compiled core accepts int32 only (`TypeError: symmetric_strength_of_connection():
   incompatible function arguments`, in `amg_variants.sa_csr_spec` of the half- and full-size JSONs). The
   cell-block (BSR 4x4) variants have 16 × fewer index entries, stay in int32 and are how PyAMG handles this
   operator at scale; at the shipped sizes they converge like the spec variant (28 vs 30 it) or better.

## Stage 3 — full size 160x160x480 (D6.4 fine level)

Instance: 49 152 000 unknowns, 5 208 269 628 nnz (106.0/row; the README's 103/row estimate was slightly low),
built in 151 s with 8 workers / 124 GB host RAM, saved to `/scratch/rmuenste/gpu-poisson/p1_160x160x480.npz`
(63.7 GB, 83 s); reload 55 s.

| storage / matvec | fine operator on device | matvec | peak device memory | Jacobi-CG | AMG-CG sa_bsr_lin | AMG-CG geo_lin |
|---|---|---|---|---|---|---|
| **symmetric** (upper triangle 2.63 G nnz, atomics kernel) — deliverable | 31.9 GB | 39.4 ms (811 GB/s on stored bytes) | **54.9 GB** driver view (nvidia-smi 80.8 GB incl. transients during the hierarchy upload; mempool 33.1 GB) | **520 it, 24.4 s** (46.9 ms/it), relres 9.8e-9, err 2.9e-7 | **11 it, 2.41 s** (219 ms/it, V-cycle 189 ms), relres 5.8e-9, err 8.4e-9 | 169 it, 41.2 s, err 5.6e-7 |
| full, own CSR warp kernel, with sa_bsr_lin | 62.9 GB | 49.5 ms (1272 GB/s) | **OOM**: `Out of memory allocating 393,216,000 bytes (allocated so far: 84,043,158,016 bytes)` while uploading the 17.3 GB hierarchy | 520 it, 29.5 s | does not fit | – |
| full, own CSR warp kernel, geo_lin only (job 147441) | 62.9 GB | 49.3 ms | 78.9 GB (nvidia-smi 75.3 GB) — fits, 1 GB to spare | – | – | 169 it, 49.1 s (V-cycle 234 ms) |

sa_bsr_lin hierarchy at full size: 4 levels (49.15 M / 1 866 240 / 69 172 / 2592), nnz 5.26 G / 196 M / 7.4 M /
0.27 M, operator complexity 1.039, 17.3 GB on the device in f32 (P, R, coarse operators), PyAMG host setup
248 s single-threaded, device setup 16 s, host RSS 291 GB (341 GB when the full CSR is also resident). The
symmetric path therefore fits the full size with ~25 GB headroom, as the README's table predicted (46 GB
estimated vs 54.9 GB measured peak; the fine matrix is 31.9 GB, the hierarchy 17.3 GB instead of the assumed
×0.4, vectors ~3 GB). The ¾-grid fallback was not needed. In full (non-symmetric) storage the full size fits
only with the lean geo_lin hierarchy (complexity 1.14) and then converges 20 × slower; the sa_bsr_lin
hierarchy does not fit beside the 63 GB matrix (one earlier geo_lin OOM in job 147436 was an artefact of a
not-yet-released failed hierarchy in the same process; the harness was fixed and the clean rerun is the row
above).

Jacobi-CG counts scale as expected with h⁻¹: 141 (h = 1/40) → 423 (1/128) → 520 (1/160) on the smooth RHS,
far below the 20 000 cap; AMG-CG stays at 11 iterations from 96 k to 49 M unknowns (scalable preconditioner).

## Headline interpretation

Seconds per pressure solve at the D6.4 size on one A100: **2.4 s (AMG-CG) / 24 s (Jacobi-CG)**, of which the
fine-level matvec (6 per AMG-CG iteration: 1 in CG, 5 in the V-cycle) is 39 ms × 66 ≈ 2.6 s minus overlap
— the solve is bandwidth-bound on the 32 GB symmetric matrix (0.8 TB/s effective; the A100's HBM2e is
~2 TB/s, the atomics path and the 4-dof row structure leave room for a ~2 × faster kernel). For comparison the
campaign spends ~113 s per time step on 121 CPU ranks for the whole step; the pressure-solve share of that was
not measured here (README §7), so the ratio is not given.

Caveats (from the brief): matrix values are those of the synthetic construction (uniform grid, unit aspect,
no FBM, pure Neumann box walls, constant mode pinned by a 1e-6 shift); sparsity, block layout, symmetry and
conditioning are those of the real operator. The AMG setup ran on the host (PyAMG, single thread, 248 s at
full size) and is excluded from the per-solve time; on a fixed mesh without FBM the matrix does not change
between steps, so the setup is a one-off. CG vectors and the fine matvec are f64; hierarchy and coarse
smoothers f32 (flexible CG absorbs the slight non-symmetry; the true residual equals the recursive one to
1e-16). The symmetric-storage matvec uses f64 atomics, so runs are not bitwise reproducible (histories agree
to ~1e-15). Timings are from a quiet, exclusively used node.

## What a GPU BoomerAMG test would add, and what it needs (not executed, per README §9)

BoomerAMG is the classical-AMG counterpart of what was measured here and the integration-relevant one
(FeatFloWer already carries `USE_HYPRE`). It would answer three things the CuPy/PyAMG route cannot: (i) setup
on the device — seconds instead of the 248 s single-threaded host setup — and no int32 index restriction
(`--enable-bigint` or mixed-int builds); (ii) whether classical coarsening (PMIS + ext+i interpolation,
l1-Jacobi relaxation, aggressive first-level coarsening) matches the 11 iterations of smoothed aggregation on
the 4-dofs-per-cell P1-disc stencil — with the block structure passed as `HYPRE_BoomerAMGSetNumFunctions(4)`
(unknown-based coarsening, the analogue of the sa_bsr variants) or the nodal-coarsening option, and the
near-nullspace handled through the interpolation rather than explicit candidates; (iii) the memory of a
classical hierarchy (operator complexity expected 1.3–1.6, i.e. 40–50 GB on top of the matrix): hypre holds
the fine matrix in full (non-symmetric) storage, 63 GB at full size, so on one 80 GB A100 the realistic
BoomerAMG test is the ½ grid (31 GB matrix, as in stage 2) or the ¾ grid, or the full size only with a
unified-memory build. Needs: a CUDA build of hypre 3.1.0 in the user's space (`--with-cuda
--with-gpu-arch=80 --enable-unified-memory` [or device memory], cuda/12.6 + gcc 13.2, optionally
`--enable-bigint`) — every installed `/sfw/hypre/v3.1.0/*` is CPU-only; a ~150-line C driver that reads the
instance (the `.npz` arrays, or the harness's slab files) into a `HYPRE_IJMatrix` with
`HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE)` and runs `HYPRE_PCG` + BoomerAMG with the GPU-recommended
settings; one tardis job of the shape of stage 2/3 (instances are already on `/scratch/rmuenste/gpu-poisson/`,
~1 h including build). The stage-2/3 JSONs give the Jacobi-CG and SA-AMG-CG numbers to compare against.

## Files

```
results/environment.txt                                stage 0 (job 147429)
results/gpu_p1_20x20x60.json, gpu_p1_40x40x120.json   stage 1 (job 147430): 6 AMG variants, Jacobi/block-Jacobi, rough RHS, sym validation
results/stage1_p1_40x40x120_{sym,full_kernel,omega1}.json
results/gpu_p1_128x128x384.json                       stage 2 deliverable (job 147433, full storage, cuSPARSE)
results/stage2_p1_128x128x384_{sym,full_kernel}.json  (jobs 147432, 147433)
results/gpu_p1_160x160x480.json                       stage 3 deliverable (job 147436, symmetric storage)
results/stage3_p1_160x160x480_full_kernel.json        full storage, sa_bsr_lin OOM (job 147436)
results/stage3_p1_160x160x480_full_kernel_geo_lin.json  full storage, geo_lin fits (job 147441)
results/stage{0,1,2,2b,3,3c}_<jobid>.log, results/*_nvsmi_<jobid>.csv
instances (NOT in the repo, node-local on tardis): /scratch/rmuenste/gpu-poisson/p1_128x128x384.npz (31.9 GB), p1_160x160x480.npz (63.7 GB)
```
