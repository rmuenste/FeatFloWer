# GEMINI.md - FeatFloWer Repository Guide

This file documents the knowledge and workflows for the FeatFloWer repository, specifically tailored for the Gemini agent. It consolidates information from `AGENTS.md`, `CLAUDE.md`, and practical execution experience.

## 1. System Environment (RHEL 9.7)

This repository relies on specific environment modules. You **MUST** load these in every shell session or script execution.

```bash
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
```

## 2. Build System

- **Generator:** Ninja is preferred (`-G Ninja`).
- **Compilers:** Use MPI wrappers (`mpicc`, `mpicxx`, `mpifort`).

### Standard Configuration (PE Serial Mode)

This is the standard configuration for running benchmarks like `q2p1_bench_sedimentation` with rigid body physics.

```bash
cmake -S . -B build-ninja-release -G Ninja 
  -DCMAKE_BUILD_TYPE=Release 
  -DBUILD_APPLICATIONS=ON 
  -DUSE_PE=ON 
  -DUSE_PE_SERIAL_MODE=ON 
  -DUSE_JSON=ON 
  -DCMAKE_C_COMPILER=mpicc 
  -DCMAKE_CXX_COMPILER=mpicxx 
  -DCMAKE_Fortran_COMPILER=mpifort
```

**Critical Flags:**
- `-DUSE_PE=ON`: Enables the Physics Engine.
- `-DUSE_PE_SERIAL_MODE=ON`: Enables Serial PE mode (essential for large particles/benchmarks).
- `-DUSE_JSON=ON`: **Mandatory** for Serial Mode to load `example.json` configuration.

## 3. Workflow: Building & Partitioning

### 1. Build Application & Metis
```bash
cmake --build build-ninja-release --target q2p1_bench_sedimentation -- -j8
cmake --build build-ninja-release --target metis -- -j8
```

### 2. Mesh Partitioning (The "Metis Dance")
The partitioning process requires `libmetis.so` to be in the *application* directory, not just the build directory.

```bash
# 1. Copy library
cp build-ninja-release/extern/libraries/metis-4.0.3/Lib/libmetis.so 
   build-ninja-release/applications/q2p1_bench_sedimentation/

# 2. Run Python Partitioner (from app dir)
cd build-ninja-release/applications/q2p1_bench_sedimentation
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" 
python3 ../../../tools/PyPartitioner.py <N_PARTS> 1 1 NEWFAC <PROJECT_FILE>
```
*Note: For 32-core nodes, usually use 31 partitions (1 master + 31 workers).*

## 4. Simulation Runtime Setup

### Required Files in Run Directory
1.  **Executable**: `q2p1_bench_sedimentation`
2.  **CFD Config**: `_data/q2p1_param.dat` (Controls `MaxNumStep` for CFD)
3.  **PE Config**: `example.json` (Required if `-DUSE_JSON=ON`. Controls particle physics parameters).
    *   *Source:* `libs/pe/pe/interface/example.json`
4.  **Multigrid**: `_data/MG.dat` (Often staged automatically, but verify).
5.  **Metis Lib**: `libmetis.so` (If re-partitioning is needed).

### SLURM Submission
Always use explicit module loads in the SBATCH script.

```bash
#!/bin/bash
#SBATCH ...
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
mpirun -np <N_PROCS> ./application
```

## 5. Key Learnings & Pitfalls

- **Simulation Duration**: Controlled by `SimPar@MaxNumStep` in `_data/q2p1_param.dat` (CFD steps), but `example.json` also has `timesteps_`. Usually, CFD drives the loop.
- **Serial PE Mode**: In this mode, each CFD domain runs a local PE instance. They are synchronized via FBM/CFD communication. It is critical for benchmarks where particles span multiple domains.
- **Param File**: `q2p1_param.dat` is often binary or text. Use `cat` to inspect if `read_file` fails on "binary" detection.
- **Logs**: Always redirect build and run output to log files to avoid flooding the context.
