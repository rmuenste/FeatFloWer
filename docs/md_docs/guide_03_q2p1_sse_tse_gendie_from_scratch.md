# FeatFloWer Guide 03: Complete From-Scratch Run of `q2p1_sse` for SSE/TSE (and gendie context)

This guide is the third runbook in the from-scratch series.
It documents the full workflow to configure, build, stage, and submit a **short TSE test run** with the
historically named `q2p1_sse` application.

Important scope clarification:

- application name: `q2p1_sse`
- physical use cases: **SSE**, **TSE**, and (with dedicated setup) **gendie** workflows
- this guide run command: `python3 ./e3d_start.py -n 32 -f _ianus/TSE/Conv -a 0 --short-test`

This command already triggers mesh preparation/partitioning logic inside the workflow.

---

## 1) Why `q2p1_sse` Needs Special Configure Flags

For this application family, these options are essential:

- `-DUSE_CGAL=ON` (**mandatory** for `q2p1_sse` and related targets to be added)
- `-DUSE_HYPRE=ON`
- `-DUSE_PE=OFF` (current SSE/TSE setup in this workflow)
- `-DENABLE_FBM_ACCELERATION=OFF`

Important dependency rule:

- `ENABLE_FBM_ACCELERATION` depends on PE/HashGrid integration.
- In current build logic, if `USE_PE=OFF`, FBM acceleration is forced/stays `OFF`.

---

## 2) Environment Setup (RHEL module workflow)

### For human users (interactive shell)

```bash
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13
module load openmpi/options/interface/ethernet
module load openmpi/4.1.6
module list
```

Expected loaded modules:

- `gcc/latest-v13`
- `openmpi/options/interface/ethernet`
- `openmpi/4.1.6`

### For AI agents / non-interactive commands

Each command should include module loading in the same shell invocation:

```bash
bash -lc '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13
  module load openmpi/options/interface/ethernet
  module load openmpi/4.1.6

  # actual command here
  cmake --version
'
```

---

## 3) Configure Build Directory

From repository root:

```bash
cmake -S . -B build-sse-ninja-mod -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DUSE_CGAL=ON \
  -DUSE_HYPRE=ON \
  -DUSE_PE=OFF \
  -DENABLE_FBM_ACCELERATION=OFF
```

Recommended validation:

```bash
rg -n "ENABLE_FBM_ACCELERATION:BOOL|USE_CGAL:BOOL|USE_HYPRE:BOOL|USE_PE:BOOL" \
  build-sse-ninja-mod/CMakeCache.txt
```

Expected values:

- `ENABLE_FBM_ACCELERATION:BOOL=OFF`
- `USE_CGAL:BOOL=ON`
- `USE_HYPRE:BOOL=ON`
- `USE_PE:BOOL=OFF`

---

## 4) Build Required Targets

Build the solver and METIS:

```bash
cmake --build build-sse-ninja-mod --target q2p1_sse metis -- -j8
```

Key expected artifacts:

- `build-sse-ninja-mod/applications/q2p1_sse/q2p1_sse`
- `build-sse-ninja-mod/extern/libraries/metis-4.0.3/Lib/libmetis.so`

---

## 5) Use the Dedicated Staging Target

A dedicated staging target is available (analogous to Guide 02 style):

```bash
cmake --build build-sse-ninja-mod --target q2p1_sse_stage -- -j8
```

What `q2p1_sse_stage` prepares in `build-sse-ninja-mod/applications/q2p1_sse`:

- runtime directories (`_data`, `_data_BU`, `_ianus`, etc.)
- `partitioner/` scripts
- `e3d_start.py`, `q2p1_sse_start.py`, `conv_check.sh`
- `_data/MG.dat`
- `libmetis.so`
- `s3d_mesher` executable (critical dependency)

Quick checks:

```bash
cd build-sse-ninja-mod/applications/q2p1_sse

ls -lh q2p1_sse s3d_mesher libmetis.so _data/MG.dat
ls -ld _ianus _data_BU partitioner
```

---

## 6) Create SLURM Submission Script (2h, `nx`, 32 tasks)

Create `run_sse_short.sh` in `build-sse-ninja-mod/applications/q2p1_sse`:

```bash
#!/bin/bash
#SBATCH --partition=short
#SBATCH --constraint=nx
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=02:00:00
#SBATCH --job-name=q2p1_sse_short

source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13
module load openmpi/options/interface/ethernet
module load openmpi/4.1.6

# Critical MPI-IO tuning for this workflow
export OMPI_MCA_io=romio321
export ROMIO_CB_BUFFER_SIZE=16777216
export ROMIO_DS_WRITE=enable

cd "$SLURM_SUBMIT_DIR"

python3 ./e3d_start.py -n 32 -f _ianus/TSE/Conv -a 0 --short-test
```

Make executable:

```bash
chmod +x run_sse_short.sh
```

---

## 7) Submit and Monitor

Submit from the runtime folder:

```bash
cd build-sse-ninja-mod/applications/q2p1_sse
sbatch run_sse_short.sh
```

Monitor:

```bash
squeue -u $USER
tail -f slurm-<jobid>.out
```

---

## 8) Critical Pitfalls and Fast Fixes

### Pitfall A: Missing `s3d_mesher`

Symptom in SLURM log:

```text
/bin/sh: line 1: ./s3d_mesher: No such file or directory
```

Cause:

- Runtime directory not fully staged, or `s3d_mesher` not copied.

Fix:

1. Run staging target again:
```bash
cmake --build build-sse-ninja-mod --target q2p1_sse_stage -- -j8
```
2. Verify:
```bash
ls -lh build-sse-ninja-mod/applications/q2p1_sse/s3d_mesher
```
3. Resubmit job.

### Pitfall B: Poor write performance due missing ROMIO environment

Symptom:

- Run works but file I/O throughput is unexpectedly poor.

Cause:

- Missing exports in submission script.

Required exports:

```bash
export OMPI_MCA_io=romio321
export ROMIO_CB_BUFFER_SIZE=16777216
export ROMIO_DS_WRITE=enable
```

### Pitfall C: `q2p1_sse` target missing during configure

Cause:

- `USE_CGAL=OFF`.

Fix:

- Reconfigure with `-DUSE_CGAL=ON`.

---

## 9) Minimal Reproducible Command Block

```bash
# 1) Module environment
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13
module load openmpi/options/interface/ethernet
module load openmpi/4.1.6

# 2) Configure
cmake -S . -B build-sse-ninja-mod -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DUSE_CGAL=ON \
  -DUSE_HYPRE=ON \
  -DUSE_PE=OFF \
  -DENABLE_FBM_ACCELERATION=OFF

# 3) Build solver + metis
cmake --build build-sse-ninja-mod --target q2p1_sse metis -- -j8

# 4) Stage runtime assets
cmake --build build-sse-ninja-mod --target q2p1_sse_stage -- -j8

# 5) Submit
cd build-sse-ninja-mod/applications/q2p1_sse
sbatch run_sse_short.sh
```

---

## 10) Notes for gendie Workflows

`q2p1_sse` also underpins gendie-style pipelines (historical naming).
This guide focuses on TSE short-test execution, but the same build baseline (`USE_CGAL=ON`, staged runtime,
module-loaded batch scripts) is the foundation for gendie runs.

When extending to full gendie production runs, ensure gendie-specific tools/scripts are staged consistently
with your run directory requirements.

---

## 11) Guide Series Context

- Guide 01: baseline `q2p1_fc_ext` benchmark from scratch
- Guide 02: `q2p1_bench_sedimentation` + PE serial mode from scratch
- Guide 03 (this guide): `q2p1_sse` workflow for SSE/TSE execution with robust staging and SLURM submission

