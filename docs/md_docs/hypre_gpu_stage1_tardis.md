# Hypre GPU Stage 1 on Tardis

This note records the implemented Stage 1 path for running the existing Hypre
coarse-grid pressure solvers on the NVIDIA A100 in `tardis`. The validation
case is the `q2p1_fc_ext` flow-around-a-cylinder benchmark from
`guide_01_q2p1_fc_ext_cylinder_benchmark_from_scratch.md`.

## Result

Both supported coarse-solver selections run to completion on the GPU:

- `Pres@MGCrsSolverType = 7`: GMRES with BoomerAMG preconditioning on the full
  coarse matrix.
- `Pres@MGCrsSolverType = 8`: BoomerAMG on the geometrically reduced auxiliary
  system.

The validated matrix covers 4 and 8 MPI ranks, where rank 0 remains the
FeatFloWer master and 3 or 7 worker ranks share one A100 through CUDA MPS.
CPU and GPU runs produce matching final cylinder force values to the printed
precision. The coarse systems in this benchmark are small and latency-bound,
so GPU execution is slower than CPU execution; this stage establishes the
working data path for a later fine-level GPU solve.

## Implementation

The default CPU build remains on the vendored Hypre and keeps its existing
solver parameters. The new CMake controls are:

| Option | Purpose |
| --- | --- |
| `USE_EXTERNAL_HYPRE` | Use a Hypre install of version 2.31 or newer. |
| `EXTERNAL_HYPRE_DIR` | External Hypre install prefix. |
| `USE_HYPRE_CUDA` | Validate CUDA and unified-memory support and enable the CUDA data path. |
| `USE_HYPRE_GPU_AMG` | Select GPU-supported BoomerAMG parameters. Enabled automatically by `USE_HYPRE_CUDA`. |
| `ENABLE_HYPRE_TIMING` | Print one machine-readable `HYPRE_TIMING` record per coarse solve. |

The GPU BoomerAMG settings are PMIS coarsening (8), extended+i interpolation
(6) with truncation factor 0.25, and symmetric hybrid
Gauss-Seidel/Jacobi relaxation (6). The latter is important for this case:
l1-Jacobi (18) and Jacobi (7) both made the type-8, 8-rank case diverge, while
the GPU-supported symmetric hybrid smoother converged reliably. Hypre's nodal
system expansion is disabled in the GPU-only GMRES/PCG parameter blocks because
Hypre 2.33 performs that scalar C/F expansion on the host.

Ordinary Fortran allocations passed to the GPU IJ interface caused invalid
device copies and illegal-address failures. With `HYPRE_CUDA` defined, the six
marshalling arrays (`rows`, `ncols`, `cols`, `values`, `rhs`, and `sol`) are
therefore allocated by a small `cudaMallocManaged` C shim and exposed to
Fortran with `c_f_pointer`. The original allocatable arrays remain unchanged in
non-CUDA builds.

## Environment and builds

The validated node reported an NVIDIA A100 80 GB PCIe, driver 580.159.04, and
CUDA driver capability 13.0. The build used CUDA Toolkit 12.4, GCC 13.2.0, and
Open MPI 4.1.6. This MPI build is not CUDA-aware; Hypre was configured with
GPU-aware MPI disabled.

```bash
module purge
module load gcc/13.2.0
module load openmpi/options/interface/ethernet openmpi/options/cuda/no openmpi/4.1.6
module load cuda/12.4
```

The external library is Hypre v2.33.0, built as a static library for `sm_80`:

```bash
cmake -S stage1-work/hypre-2.33.0-src/src \
  -B stage1-work/hypre-2.33.0-cuda-build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/stage1-work/hypre-2.33.0-cuda" \
  -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_CUDA_COMPILER=/sfw/cuda/12.4/bin/nvcc \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DHYPRE_ENABLE_CUDA=ON \
  -DHYPRE_ENABLE_UNIFIED_MEMORY=ON \
  -DHYPRE_ENABLE_GPU_AWARE_MPI=OFF \
  -DHYPRE_ENABLE_BIGINT=OFF -DHYPRE_ENABLE_MIXEDINT=OFF \
  -DHYPRE_BUILD_TESTS=OFF -DHYPRE_BUILD_EXAMPLES=OFF
cmake --build stage1-work/hypre-2.33.0-cuda-build --parallel 16
cmake --install stage1-work/hypre-2.33.0-cuda-build
```

Configure the GPU application build with:

```bash
cmake -S . -B stage1-work/build-gpu -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON \
  -DUSE_HYPRE=ON -DUSE_EXTERNAL_HYPRE=ON \
  -DEXTERNAL_HYPRE_DIR="$PWD/stage1-work/hypre-2.33.0-cuda" \
  -DUSE_HYPRE_CUDA=ON -DENABLE_HYPRE_TIMING=ON
cmake --build stage1-work/build-gpu --target q2p1_fc_ext hypre_test --parallel 8
```

The matched CPU reference uses the site Hypre 2.33.0/32-bit module and
`USE_HYPRE_GPU_AMG=ON`, so the performance comparison changes only the backend:

```bash
module load hypre/2.33.0/32bit
cmake -S . -B stage1-work/build-matched-cpu -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON \
  -DUSE_HYPRE=ON -DUSE_EXTERNAL_HYPRE=ON \
  -DEXTERNAL_HYPRE_DIR=/sfw/hypre/v2.33.0/gcc13.2.0-32bit \
  -DUSE_HYPRE_GPU_AMG=ON -DENABLE_HYPRE_TIMING=ON
cmake --build stage1-work/build-matched-cpu \
  --target q2p1_fc_ext hypre_test --parallel 8
```

A separate vendored-Hypre build with all new options disabled successfully
built `ff_le_solvers`, guarding the original CPU path.

## Run and summarize

The job harness partitions the cylinder mesh for both rank counts, stages an
independent application directory for every run, starts CUDA MPS, records
`nvidia-smi dmon`, and generates a TSV summary. Submit it with the required
Tardis resources using:

```bash
sbatch --export=ALL,FF_STAGE1_RUN_ID=measured-final \
  tools/cluster_scripts/hypre_gpu_stage1_cylinder.slurm
```

Useful overrides are `FF_STAGE1_STEPS`, `FF_STAGE1_REPEATS`,
`FF_STAGE1_WARMUP_CALLS`, `FF_STAGE1_BACKENDS`, `FF_STAGE1_RANKS`,
`FF_STAGE1_SOLVER_TYPES`, and `FF_STAGE1_USE_MPS`. Results are written below
`stage1-work/results/<run-id>/`; this local working tree is ignored by Git.

Each timing is the maximum across the worker communicator. Setup, solve, total,
iteration count, and final relative residual are emitted separately. The
summarizer reports the median solve call after discarding the configured warmup
calls.

## Measured cylinder results

The final job (`135340`) ran 10 time steps and three repetitions per matrix
entry, discarding the first five Hypre calls in each repetition. Times below
are the median of the three per-run medians, in seconds per Hypre coarse-solver
call; iteration counts are medians. All 24 application runs completed in 129
seconds of allocated-node wall time.

| Ranks | Solver | CPU setup | CPU solve | GPU setup | GPU solve | CPU iters | GPU iters |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | type 7, GMRES+AMG | 0.000353 | 0.006253 | 0.012470 | 0.045192 | 48 | 16 |
| 8 | type 7, GMRES+AMG | 0.000364 | 0.005919 | 0.012850 | 0.048854 | 47 | 16 |
| 4 | type 8, AMG | 0.000325 | 0.000314 | 0.034762 | 0.043282 | 12 | 10 |
| 8 | type 8, AMG | 0.000462 | 0.000787 | 0.035981 | 0.066211 | 12 | 11 |

The maximum reported relative residual is `9.90e-6` for type 7 (whose GMRES
tolerance is `1e-5`) and `9.84e-8` for CPU / `9.69e-8` for GPU type 8 (whose
AMG tolerance is `1e-7`). At the final simulated time `0.1`, CPU and GPU drag
coefficients agree to the printed precision for every matched case; the largest
printed lift-coefficient difference is `4.43e-8`. One-second `nvidia-smi dmon`
sampling observed up to 70% SM activity, confirming device execution (memory
utilization peaked at 1%, consistent with the tiny coarse systems).

## Failure history and limitations

- Hypre accepted ordinary Fortran host arrays during IJ setup but later issued
  invalid device-to-device copies and an illegal address. Managed marshalling
  allocations fixed this interface boundary.
- Enabling `SetNumFunctions(4)` and `SetNodal(3)` entered a host-only Hypre 2.33
  scalar C/F path. These calls remain untouched for CPU and are omitted for GPU.
- Relaxation types 18 and 7 were numerically unstable for direct AMG at 8 ranks
  in this benchmark. The supported symmetric hybrid option 6 passes the full
  validation matrix.
- Multiple MPI worker processes share one A100. MPS is required to avoid severe
  context-serialization overhead, but cannot remove the launch latency of these
  very small coarse systems.
- Stage 1 does not update the vendored Hypre 2.25 tree. CUDA runs use the external
  Hypre 2.33 installation; a dependency upgrade is separate work.

## FAC3D follow-up

The larger `q2p1_fac3d` cylinder setup was staged from
`_adc/3D_FAC_FBM/file.prj` with the same `_data`, `start`, partitioned `_mesh`,
and output-directory layout as the established regression build. The reusable
job is `tools/cluster_scripts/hypre_gpu_stage1_fac3d.slurm`; its default
`MaxMeshLevel=5` matches the full reference hierarchy and it reserves enough
Tardis memory for that case.

A level-4 diagnostic already reaches 864,049 finest-level DOFs per worker at
4 ranks. The level-1 pressure storage has 816 local entries, but the wrapper
submits only 204 rows per worker to Hypre (`lPMat%nu/4`). This is the key
scaling detail: `SetUp_HYPRE_Solver` constructs the matrix at
`Pres@MGMinLev`, which is 1 in the FAC3D configuration. Increasing
`MaxMeshLevel` to 5 makes the surrounding FeatFloWer calculation much larger
but leaves the Hypre problem unchanged. It therefore cannot improve the GPU
comparison unless the pressure level handed to Hypre is changed in a later
stage.

With the original attraction-weight mesh deformation enabled, neither Hypre
coarse-solver option supplied a valid FAC3D baseline:

| Configuration | Type 7 result | Type 8 result |
| --- | --- | --- |
| Original CPU AMG parameters | Hit 80 iterations; relative residual `2.4e2` to `3.7e2` | `NaN` residual and cylinder force |
| Matched GPU-compatible parameters on CPU | Hit 80 iterations; residual up to `2.26e2` | `NaN` residual and cylinder force |
| CUDA Hypre 2.33 | Hit 80 iterations; residual up to `2.21e2` | `NaN` residual and cylinder force |

The application reaches its normal finalization message in these one-step
runs, so process exit status alone is not a correctness gate. The result
summarizer now labels rows `ok`, `residual_failed`, `nonfinite`, `app_failed`,
or `no_samples`. No performance result was reported from those runs because
there was no converged CPU result to compare against. These failures were
subsequently traced to the application-specific mesh deformation described
below and should not be used to assess the undeformed FAC3D configuration.

### FAC3D mesh-deformation follow-up

`q2p1_fac3d` forced `bFAC3D_CylUmbrellaWeight` on during
`General_init_ext`, applied eight `CylinderAttraction` passes, and followed them
with two weighted umbrella-smoothing passes. A separate initialization
umbrella loop could also modify the coordinates. All three deformation loops
are now disabled and the weighting flag is forced off. Standard refinement and
boundary parametrization remain enabled.

Native CPU job 135763 ran the undeformed level-4 mesh for one step at four MPI
ranks. Type 7 now passes with a maximum relative residual of `6.85e-6`, where
the attraction-enabled run previously diverged by orders of magnitude. The
initial visualization is stored under result ID
`fac3d-no-deformation-initial` as `_vtk/initial_mesh.pvtu`.

For a representative worker piece, comparing the new coordinates with the
previous attraction-enabled VTK shows 516,857 of 864,049 points moving by more
than `1e-7`; the maximum displacement caused by the old path was `5.76e-2` and
the RMS displacement was `1.69e-2`. The deformation was therefore material.

After visual confirmation of the undeformed mesh, correctness job 135764 ran
both solver types with the matched CPU and CUDA builds. All four cases passed.
Timing job 135765 then ran two time steps, three repetitions, and four MPI ranks
(three workers) for both backends and solver types. Each repetition produced
12 Hypre calls; the first two were treated as warm-up. The table reports the
median of the three per-repetition medians:

| Solver | Backend | Setup | Solve | Total | Iterations | Maximum residual |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Type 7 | CPU | 0.001735 s | 0.014305 s | 0.016055 s | 19 | `9.87e-6` |
| Type 7 | GPU | 0.013708 s | 0.095774 s | 0.111152 s | 20 | `9.45e-6` |
| Type 8 | CPU | 0.000997 s | 0.000663 s | 0.001685 s | 11 | `6.84e-8` |
| Type 8 | GPU | 0.036849 s | 0.060725 s | 0.099035 s | 11.5 | `7.63e-8` |

Both GPU configurations are numerically valid, but type 7 is 6.92 times slower
and type 8 is 58.78 times slower than the matched CPU backend by total Hypre
time. The A100 reached only 35% sampled SM utilization and was active in 22 of
1,933 one-second samples. This confirms the expected launch-latency behavior
for the 204-row local Hypre problem: FAC3D is large outside the coarse solver,
but it is not a larger Stage 1 GPU workload.

All runs printed identical zero cylinder forces at the short final time
`0.002`; this is consistent across backends but is not a discriminating force
validation. The residual and iteration checks are the useful correctness gates
for this short run. Job 135765 completed with exit code 0 in 32:32. Detailed
results are in `stage1-work/results/fac3d-nodeform-timed/summary.tsv`. Set
`FF_FAC3D_KEEP_OUTPUT=0` for repeated timing jobs to discard their redundant
initial VTK and restart files after each case while retaining logs and timing
summaries.

## SSE/TSE follow-up

The `q2p1_sse` TSE workflow from Guide 03 was also built and staged in the
matched CPU and CUDA build trees. Both use GCC 13.2 and the site CGAL 6.0.1
package; the only Hypre difference remains the configured execution policy.
The reusable job is
`tools/cluster_scripts/hypre_gpu_stage1_sse.slurm`. It preserves the documented
`e3d_start.py` setup, meshing, JSON partitioning, 32-rank launch, and ROMIO
environment while allowing the number of steps, mesh quality, pressure coarse
level, and solver type to be controlled through `FF_SSE_*` variables.

The documented coarse short test is correct with type 8 on both backends, but
the representative worker owns only 30 Hypre rows. Increasing the TSE
`MeshQuality` grows the level-1 matrix without changing the physical case:

| Mesh quality | Local Hypre rows | Type 8 CPU | Type 8 GPU |
| --- | ---: | --- | --- |
| `coarse` | 30 | Pass | Pass |
| `medium` | 73 | Pass | Pass |
| `fine` | 136 | Pass | Intermittent residual failure |

Moving `Pres@MGMinLev` from 1 to 2 is not a usable scaling shortcut for this
application. With either `MaxMeshLevel=2` or 3, both CPU and GPU stop in
`SETLEV` before constructing Hypre. The harness retains this setting as an
explicit experimental control, but defaults to the valid level 1.

Job 135353 exercised the fine mesh for 12 steps, three repeats, 32 ranks, and
both solver types. Each run produced 24 Hypre calls; timing medians below omit
the first two calls, while correctness checks include every call:

| Solver | Backend | Median setup | Median solve | Maximum residual | Status |
| --- | --- | ---: | ---: | ---: | --- |
| Type 7 | CPU | 0.00919 s | 0.11752 s | `1.20e-4` | Failed |
| Type 7 | GPU | 0.08853 s | 4.13503 s | `3.68e-2` | Failed |
| Type 8 | CPU | 0.00450 s | 0.00629 s | `9.46e-8` | Pass |
| Type 8 | GPU | 0.11299 s | 0.43154 s | `2.22e-1` | Failed |

The type 8 GPU failure is repeatable rather than a one-time initialization
artifact: several calls hit the 50-iteration cap in every repeat, although
other calls converge. As with FAC3D, `q2p1_sse` still prints its successful
finalization message. Therefore no SSE speed comparison is reported. The
summarizer now treats application exit and all residuals as correctness gates;
`--warmup-calls` affects timing samples only.

After the CGAL-enabled SSE builds relinked the shared application libraries,
job 135354 reran the 2D cylinder type 8 smoke test. CPU and GPU both passed,
and the printed drag and lift values remained identical, confirming that the
validated Stage 1 cylinder baseline was preserved.
