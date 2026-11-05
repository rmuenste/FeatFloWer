# GPU pressure-Poisson test on the A100 (tardis) — hand-off

Owner request (2026-10-02): solve a synthetic system of the SIZE and STRUCTURE of the
FeatFloWer pressure-Poisson operator of the largest campaign run on ONE A100 80 GB, as a
test of GPU solvers. No integration into FeatFloWer, no production use. This folder holds
the generator, two reduced instances, a CPU reference, and this brief. The campaign
context is in `docs/md_docs/dns_practitioners_guide.md` and the datasheet
(`docs/md_docs/dns_validation_datasheet.csv`); nothing here touches them.

## 1. What is being modelled

FeatFloWer solves incompressible Navier–Stokes with Q2/P1 finite elements and a discrete
projection scheme: per time step a nonlinear solve for the intermediate velocity (Q2, three
components), then a **pressure-Poisson solve** on the P1-discontinuous pressure space,
then the velocity correction. The pressure operator is the Schur-complement
C = Bᵀ M_l⁻¹ B (B = discrete divergence, M_l = lumped velocity mass), assembled as the
"C matrix" in `source/src_quadLS/QuadSc_struct.f90` (`Create_LinMatStruct` /
`Get_CMatStruct`):

- 4 unknowns per hexahedral cell (cell mean + three slopes, P1-discontinuous);
- row of cell K couples to all 4 unknowns of every cell sharing a VERTEX with K
  (27-cell neighbourhood → 108 entries per interior row; measured average 103 with
  boundaries, from the D6.4 G0 log: 256 000 unknowns, 26.3 M non-zeros);
- symmetric positive semi-definite, one constant mode (pure Neumann pressure);
- conditioning O(h⁻²) like any Laplacian; solved today by the code's own geometric
  multigrid with a UMFPACK coarse grid.

**The largest campaign system** (D6.4 tier-S settling-spheroid runs, `q2p1_dns_rundir_d64_s_*`,
8 × 8 × 24 box at mesh level 4, h = 0.05): 12.288 M cells = **160 × 160 × 480**, i.e.
**49.15 M pressure unknowns, ≈ 5.0 G non-zeros** (see `tools/gpu_poisson` sizing in §4).
Velocity unknowns are ~300 M on top, not part of this test.

## 2. The synthetic system (`synth_p1_pressure.py`)

Built the same way as the real operator, C = Bᵀ M_l⁻¹ B, on a uniform hexahedral grid of
the unit-aspect box, with:
- P1-discontinuous pressure (identical layout: unknown 4K+0 = mean of cell K, 4K+1..3 = slopes),
- a triquadratic Q2 velocity on the 27-node lattice, as in the code (the Q2/P1-disc pair is
  inf-sup stable; a Q1 velocity was tried first and rejected: it leaves spurious pressure
  modes at the shift level and thousands of CG iterations),
- B assembled exactly with a 2 × 2 × 2 Gauss rule, M_l the lumped Q1 mass,
- the constant mode pinned by a diagonal shift (`--shift`, default 1e-6) → strictly SPD.

Values are those of this construction (uniform grid, no FBM, no walls other than the box);
sparsity pattern, block layout, symmetry and conditioning are what a solver experiences.
Row length 98.9 per row on the 20x20x60 grid (walls lower the average; the code logs 103
at level 3). Measured on 20x20x60: lambda_min 2.8e-3, lambda_max 1.26, condition 456 at
h^-2 = 400, Jacobi-CG 139 iterations to 1e-8 - a Poisson-class operator.

Output (`.npz`): `indptr` int64 (nu+1), `indices` int32, `data` float64 — CSR as a 32-bit-index
GPU solver consumes it — plus `b` (= C·x_true for a smooth x_true) and `x_true`, so a solve is
verifiable, and `nx, ny, nz, h`.

```
module load python/3.13.5
python3 synth_p1_pressure.py NX NY NZ [--out F.npz] [--shift 1e-6] [--dry-run]
```
`--dry-run` prints the memory budget only. Building with scipy: measured 13 s and 4.2 GB peak host RAM for 768 000 unknowns, i.e.
~5.4 GB per million unknowns -> the half-size instance (25 M) needs ~140 GB and ~10 min,
the full size (49 M) ~270 GB and ~20 min: both fit tardis's 480 GB host RAM. If memory
becomes the limit, add slab-wise assembly (the construction is cell-local).

Instances already built (`instances/`):

| file | cells | unknowns | nnz | CSR size |
|---|---|---|---|---|
| `p1_20x20x60.npz` | 24 000 | 96 000 | 9.49 M (98.9/row) | 0.11 GB |
| `p1_40x40x120.npz` | 192 000 | 768 000 | 79.0 M (102.9/row) | 0.95 GB |

`figures/p1_structure.png`: sparsity of two x-planes, one interior cell's 4 × 108 block row,
and the row-length histogram of the 40×40×120 instance.

## 3. CPU reference (`cpu_reference.py`, results in `results/`)

Per instance: λ_min (inverse iteration), λ_max (Lanczos), condition number, Jacobi-CG and
4×4-block-Jacobi-CG iteration counts to 1e-8 on a rough mean-free random RHS, and the
smooth-RHS solve checked against x_true. The GPU solver must reproduce the smooth-RHS
solution to the same tolerance and should beat the Jacobi-CG count by an order of
magnitude with AMG. See `results/cpu_reference.log` / `results/*.cpu.json`.

## 4. Sizing on one A100 80 GB (fine level only)

Full-size grid 160 × 160 × 480, 49.15 M unknowns, ~5.0 G nnz at 103/row:

| storage | matrix | + 8 vectors | + AMG hierarchy (×1.4) |
|---|---|---|---|
| full CSR, f64 values + i32 indices | 60 GB | 63 GB | 87 GB — does not fit |
| upper triangle only (symmetric) | 30 GB | 34 GB | 46 GB — fits |
| upper triangle, f32 hierarchy | 30 GB | 34 GB | 42 GB — fits |

Reduced grids for full (non-symmetric) storage:

| grid | unknowns | full CSR | AMG-CG total |
|---|---|---|---|
| 150 × 150 × 450 (¾) | 40.5 M | 49 GB | 72 GB — fits, tight |
| 128 × 128 × 384 (½) | 25.2 M | 31 GB | 45 GB — fits |
| 100 × 100 × 300 (¼) | 12.0 M | 15 GB | 21 GB — fits |

Vectors: 0.39 GB each at full size. Indices MUST stay 32-bit (5 G nnz × 8 B would be 40 GB
alone); 49 M unknowns are well inside int32.

## 5. Solver recommendation and library options

**Method:** conjugate gradient (the operator is SPD) preconditioned by one V-cycle of an
algebraic multigrid, mixed precision (Krylov vectors and residual in f64, the AMG hierarchy
and smoothers in f32), with a smoother that respects the 4 × 4 cell blocks (block-Jacobi or
Chebyshev/Jacobi; damped Jacobi is the safe GPU default). Plain Jacobi-CG is the baseline:
its iteration count grows as h⁻¹ and will be in the thousands at full size.

**Libraries, in the order to try on tardis** (CUDA 12.6 is installed: `module load cuda/12.6`;
`libcusparse`/`libcusolver` under `/sfw/cuda/12.6/lib64`; no AmgX, Ginkgo, PETSc or CuPy
preinstalled; PyPI is reachable and `pip download` of `cupy-cuda12x` and `pyamg` wheels for
python/3.13.5 works):

1. **CuPy + PyAMG (fastest to a first number, Python only).** `pip install --user cupy-cuda12x
   pyamg` under `module load python/3.13.5 cuda/12.6`. Load the `.npz`, build the smoothed-
   aggregation hierarchy with PyAMG on the CPU (host RAM is plentiful), copy each level's
   CSR to the GPU as CuPy sparse matrices (f32 for the hierarchy), run CG in CuPy with the
   V-cycle as the preconditioner (CuPy's `cupyx.scipy.sparse.linalg.cg` accepts a
   `LinearOperator` M). Memory-wise identical to the table above. Also gives the plain
   Jacobi-CG GPU baseline in a few lines. Limitation: PyAMG's setup is CPU-only and
   single-threaded; at full size expect tens of minutes of setup, acceptable for a test.
2. **hypre BoomerAMG on the GPU.** BoomerAMG is the right algorithm class (classical AMG,
   proven on Poisson, GPU-capable since hypre 2.20). But EVERY installed hypre
   (`/sfw/hypre/v3.1.0/*`, 14 compiler variants, `module avail hypre`) is a CPU-only build —
   no `HYPRE_USING_CUDA` in any `HYPRE_config.h`. Using it on the A100 means building hypre
   3.1.0 from source with `--with-cuda --enable-unified-memory` (or `--with-cuda` + device
   memory) against cuda/12.6 and gcc 13.2, in the user's space; then a ~150-line C driver
   that reads the CSR, builds a `HYPRE_IJMatrix`, and runs `HYPRE_PCG` + `BoomerAMG` with the
   GPU-recommended settings (coarsening PMIS, interpolation ext+i, Jacobi/l1-Jacobi
   relaxation, aggressive coarsening on the first level, `HYPRE_SetMemoryLocation(DEVICE)`).
   This is the route that answers "would BoomerAMG do it" for a later integration, since
   FeatFloWer already has a `USE_HYPRE` option. Memory: BoomerAMG's operator complexity on
   this stencil is the unknown; with aggressive coarsening expect 1.3–1.6.
3. **NVIDIA AmgX** (classical or aggregation AMG + PCG, CUDA-native, best raw speed, takes CSR
   with 32-bit indices directly through its C API; needs a source build, Apache licence) —
   equally valid, slightly more build effort than hypre, no path into FeatFloWer.
4. **Ginkgo** (C++, PGM-AMG + CG, mixed precision built in) — a good third option if 1–3 fail.

Recommendation: do 1 first (hours, gives the baseline and the AMG-CG number at ½ and full
size), then 2 (BoomerAMG is the integration-relevant one), report both.

## 6. The tardis node

Slurm node `tardis`: 1 × NVIDIA Ampere (A100 80 GB expected — verify with `nvidia-smi`),
32 cores (16 physical), 480 GB RAM, in partitions short (2 h), med (8 h), long (2 d),
ultralong. Request it with `--nodelist=tardis --gres=gpu:ampere:1 --partition=long` (or med);
check `squeue -w tardis` first, it is shared. The warehouse volume (`/data/warehouse17`) is
99 % full — write instances to the scratch area or the node's local disk, not here, and keep
the full-size `.npz` (60 GB) out of the repo.

## 7. Deliverables for the agent

- `results/gpu_<instance>.json` per instance run: library, solver settings, memory used on
  the device (peak), setup time, solve time, iteration count, final residual, relative error
  vs `x_true` on the smooth RHS, and the Jacobi-CG baseline count on the same instance.
- Instances to run: ½ size (128×128×384, build it) in full storage with route 1, then full
  size (160×160×480) in symmetric storage if the library supports it, else ¾ size.
- A short `results/REPORT.md`: did it fit, how many iterations, wall time per solve, and the
  one number the owner wants — seconds per pressure solve at the D6.4 size on one A100
  versus the code's CPU multigrid (the campaign runs spend ~113 s per TIME STEP on 121 ranks
  for the whole step at this size; the pressure solve is a fraction of that, to be measured
  separately if needed).

## 8. Files

```
gpu-poisson-a100/
  README.md                 this brief
  synth_p1_pressure.py      generator (same as tools/gpu_poisson/synth_p1_pressure.py)
  cpu_reference.py          CPU spectrum + CG reference
  instances/                p1_20x20x60.npz, p1_40x40x120.npz (built here; bigger ones: build on tardis)
  figures/p1_structure.png  structure of the 40x40x120 instance
  results/                  cpu_reference.log, *.cpu.json; the agent adds gpu_*.json and REPORT.md
```

## 9. Testing plan (owner decisions, 2026-10-02)

Walltime on tardis: **8 h per job** (`--partition=med --time=08:00:00 --nodelist=tardis
--gres=gpu:ampere:1`), one job per stage below; check `squeue -w tardis` before submitting.

**Stage 0 — environment (login node, no GPU needed).**
`module load python/3.13.5 cuda/12.6`; `pip install --user cupy-cuda12x pyamg` (wheels fetch
from PyPI; both verified downloadable for cp313). On tardis confirm `nvidia-smi` shows the
A100 and its memory, and `python3 -c "import cupy; cupy.cuda.Device(0).compute_capability"`.
Record versions in `results/environment.txt`.

**Stage 1 — correctness and baseline on the shipped instances (minutes).**
On `instances/p1_40x40x120.npz`: (a) GPU Jacobi-CG in CuPy to 1e-8 on the stored `b`,
compare with `x_true` and with `results/p1_40x40x120.cpu.json` (same iteration count ±1,
same error); (b) PyAMG smoothed-aggregation hierarchy on the host (`strength='symmetric'`,
`aggregate='standard'`, `smooth='jacobi'`, `max_coarse=1000`), levels copied to the GPU as
CuPy CSR in f32, one V-cycle (Jacobi smoother, 2 pre/2 post, direct coarse solve on the
host or dense on the device) as the CG preconditioner, CG in f64; report iterations, setup
and solve time, device memory, operator complexity. This is the harness every later stage
reuses; put it in `gpu_amg_cg.py` in this folder.

**Stage 2 — half size, full storage (the first real number).**
Build `p1_128x128x384.npz` on tardis (25.2 M unknowns, ~140 GB host RAM, ~10 min; write it
to local scratch, not the warehouse volume). Run Jacobi-CG (expect thousands of iterations
— cap at 20 000 and record) and AMG-CG as in stage 1. Expected device footprint ~45 GB
(31 GB matrix in f64/i32 + hierarchy + vectors). Deliverable `results/gpu_p1_128x128x384.json`.

**Stage 3 — full size, symmetric storage.**
Build `p1_160x160x480.npz` on tardis (49.15 M unknowns, ~270 GB host RAM, ~20 min). Keep
ONLY the upper triangle on the device (CuPy: `cupyx.scipy.sparse.triu` then a symmetric
matvec y = U x + Uᵀ x − diag·x; verify it against the full matvec on the 40x40x120
instance first). Jacobi-CG and AMG-CG as before; the AMG hierarchy in f32 built on the host
from the full matrix, only its coarse levels need to live on the device in full storage
(they are small). Expected footprint ~42–46 GB. If the symmetric path proves awkward in
CuPy, fall back to the ¾-size grid 150x150x450 in full storage (49 GB matrix, ~72 GB total)
and say so. Deliverable `results/gpu_p1_160x160x480.json`.

**Stage 4 — report.** `results/REPORT.md`: per stage, did it fit (peak device memory), iteration
counts (Jacobi-CG vs AMG-CG), setup and solve wall time, error vs `x_true`; the headline =
**seconds per pressure solve at the D6.4 size on one A100**, with the caveats (synthetic
values, uniform grid, no FBM, setup done on the host). Then the one paragraph the owner
asked for: what a GPU BoomerAMG test would add and what it needs (a CUDA build of hypre
3.1.0, see §5 route 2) — not executed in this hand-off.

**Not in scope:** any change to FeatFloWer, any run on the campaign queues other than tardis,
the hypre build (deferred until the CuPy/PyAMG numbers are in).

Reporting: every stage appends to `results/`; no git commits (the lead commits the folder);
the instance files and anything over 100 MB stay out of the repo.
