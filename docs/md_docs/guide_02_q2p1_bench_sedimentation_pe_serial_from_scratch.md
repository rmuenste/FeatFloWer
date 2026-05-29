# FeatFloWer Guide 02 (Draft, Part 1): Configure and Build `q2p1_bench_sedimentation` with PE Serial Mode

This is the first part of Guide 02.
It extends Guide 01 by enabling rigid-body coupling via the PE (physics engine) library for a sphere
sedimentation benchmark.

Scope of this part:

- load required RHEL 9.7 module environment
- configure FeatFloWer with PE enabled in **serial mode**
- build `q2p1_bench_sedimentation`
- build METIS and partition the sedimentation mesh
- stage PE runtime JSON config in the benchmark run directory
- keep logs in files (avoid terminal/context flooding)

## 1) Environment Setup on RHEL 9.7 (module-based)

### For Human Users (Interactive Shell)

Load modules once at the start of your session:

```bash
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
module list
```

The loaded modules persist for your entire shell session. You can now run all subsequent commands (cmake, ninja, etc.) normally.

### For Agents/Automation (Non-Interactive)

Each command runs in a separate shell session, so module loading must be included in each command:

```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  # Your actual command here
  cmake --build build-ninja-release --target q2p1_bench_sedimentation
'
```

This guide shows **human-friendly commands** by default. Agent/automation examples are noted where relevant.

### Notes

- MPI wrappers (`mpicc`, `mpicxx`, `mpifort`) come from the loaded OpenMPI module.
- `mpifort` is preferred (modern wrapper name); `mpif90` may also exist and often maps to the same compiler.

### Agent Sandbox Requirements (Codex CLI / Similar)

Some sandboxed agents run with limited network access and cannot reach external Git hosts or the Slurm controller.
To execute the full workflow autonomously (FetchContent + build + Slurm submission), the agent needs:

- **Outbound network access** to fetch dependencies (e.g., `https://github.com/nlohmann/json.git`,
  `https://gitlab.com/libeigen/eigen.git`).
- **Cluster controller reachability** so `sbatch`, `squeue`, and `scontrol` can contact `slurmctld`.

If the agent is limited to a restricted mode (e.g., “write-workspace” only), expect failures like:
- `Could not resolve host: github.com` (FetchContent)
- `slurm_set_addr: Unable to resolve <controller>` (Slurm submission)

Practical workaround if network is blocked:
- Use `-DEIGEN=OFF` and point JSON to a pre-downloaded source via
  `-DFETCHCONTENT_SOURCE_DIR_JSON=/path/to/json-src`.

This guide assumes a mode with **both filesystem write access and outbound network access**. For Codex CLI,
the **recommended** setting for fully autonomous execution of this workflow is:

```
codex --sandbox workspace-write -c 'sandbox_workspace_write.network_access=true'
```

In earlier tests, `danger-full-access` was required, but with `workspace-write` plus explicit network access,
the full workflow (FetchContent + build + Slurm submission) works end-to-end.

A middle-ground permission mode (between “write-workspace” and “danger-full-access”) would still be ideal
for safe autonomous runs, but the setting above is currently sufficient.

## 2) Configure (CMake + Ninja) for PE Serial Mode

From repository root (assumes modules are loaded):

```bash
cmake -S . -B build-ninja-release -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DUSE_PE=ON \
  -DUSE_PE_SERIAL_MODE=ON \
  -DUSE_JSON=ON \
  -DSED_BENCH=ON \
  -DVERIFY_HASHGRID=OFF \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort
```

**Important:** `-DUSE_JSON=ON` is **required** for PE Serial Mode. `-DSED_BENCH=ON` enables the sedimentation benchmark output (e.g. `SED_BENCH_VEL` lines). `-DVERIFY_HASHGRID=OFF` disables hashgrid verification overhead. The setup functions (`setupParticleBenchSerial`, etc.) load runtime configuration from `example.json` via `SimulationConfig::loadFromFile()`. Without JSON support, this function is a no-op and all parameters remain at default values.

The configure banner should show:

- `Adding pe library to build`
- `PE SERIAL MODE ENABLED`
- `SED_BENCH is ON`
- `nlohmann/json will be downloaded to: ...`
- `Found MPI_C`, `Found MPI_CXX`, `Found MPI_Fortran`

**Agent/automation example:**

```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  cd /path/to/FeatFloWer
  cmake -S . -B build-ninja-release -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_APPLICATIONS=ON \
    -DUSE_PE=ON \
    -DUSE_PE_SERIAL_MODE=ON \
    -DUSE_JSON=ON \
    -DSED_BENCH=ON \
    -DVERIFY_HASHGRID=OFF \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_Fortran_COMPILER=mpifort
' 2>&1 | tee /tmp/configure_guide02.log
```

## 3) Build the Sedimentation Application

**Human workflow** (assumes modules loaded in current shell):

```bash
cmake --build build-ninja-release --target q2p1_bench_sedimentation -- -j8
```

To capture output to a log file:

```bash
cmake --build build-ninja-release --target q2p1_bench_sedimentation -- -j8 \
  > build_q2p1_bench_sedimentation.log 2>&1
```

Quick failure scan:

```bash
rg -n "error:|undefined reference|FAILED:|ninja: build stopped" \
  build_q2p1_bench_sedimentation.log
```

**Agent/automation workflow:**

```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  cd /path/to/FeatFloWer
  cmake --build build-ninja-release --target q2p1_bench_sedimentation -- -j8
' > /tmp/build_q2p1_bench_sedimentation.log 2>&1 && \
  echo "✅ Build completed successfully" || \
  echo "❌ Build failed - check log"
```

Verify the executable was created:

```bash
ls -lh build-ninja-release/applications/q2p1_bench_sedimentation/q2p1_bench_sedimentation
```

Expected output: executable file around 40-45 MB.

## 4) Build METIS and partition the benchmark mesh

### Build METIS

**Human workflow:**

```bash
cmake --build build-ninja-release --target metis -- -j8
```

Verify METIS library was built:

```bash
ls -lh build-ninja-release/extern/libraries/metis-4.0.3/Lib/libmetis.so
```

Expected: `libmetis.so` file around 412K.

**Agent workflow:**

```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  cd /path/to/FeatFloWer
  cmake --build build-ninja-release --target metis -- -j8
  ls -lh build-ninja-release/extern/libraries/metis-4.0.3/Lib/libmetis.so
'
```

### Partition to 31 Subdomains

This benchmark uses **31 partitions** for a 32-core node:
- 31 partitions = CFD worker subdomains
- 1 remaining MPI rank = coordinator/master rank (control + coarse problem)

**Human workflow:**

```bash
# Copy METIS library to application directory
cp build-ninja-release/extern/libraries/metis-4.0.3/Lib/libmetis.so \
   build-ninja-release/applications/q2p1_bench_sedimentation/

# Change to application directory
cd build-ninja-release/applications/q2p1_bench_sedimentation

# Run partitioner (uses relative paths from application directory)
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
python3 ../../../tools/PyPartitioner.py 31 1 1 NEWFAC _adc/benchSym/bench.prj
```

Expected output should include:
- `Calling Metis...`
- Edge cut count
- `The partitioning was successful!`

### One-Shot Staging + Partitioning Target (Recommended)

You can now do build + METIS + staging + partitioning in one command:

```bash
cmake --build build-ninja-release --target q2p1_bench_sedimentation_stage
```

To change the partition count:

```bash
cmake -S . -B build-ninja-release -DQ2P1_BENCH_SEDIMENTATION_PARTS=31
cmake --build build-ninja-release --target q2p1_bench_sedimentation_stage
```

**Agent workflow:**

```bash
bash -c '
  cd /path/to/FeatFloWer
  cp build-ninja-release/extern/libraries/metis-4.0.3/Lib/libmetis.so \
     build-ninja-release/applications/q2p1_bench_sedimentation/
  cd build-ninja-release/applications/q2p1_bench_sedimentation
  LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
  python3 ../../../tools/PyPartitioner.py 31 1 1 NEWFAC _adc/benchSym/bench.prj
' 2>&1 | tee /tmp/partition_31.log
```

### Verify Partition Files

Count partition files (should be exactly 31):

```bash
ls _mesh/NEWFAC/sub0001/ | grep -c "GRID00"
```

List all partition files:

```bash
ls _mesh/NEWFAC/sub0001/GRID*.tri
```

Expected files:
- `GRID0001.tri` through `GRID0031.tri` (31 partitions)
- `GRID.tri` (base/reference grid file, not counted as a partition)

## 5) Stage PE JSON config in runtime directory

**Critical:** PE Serial Mode setup functions (`setupParticleBenchSerial`, etc.) load runtime configuration from `example.json` via `SimulationConfig::loadFromFile()`. This file must be present in the application run directory.

### Copy example.json

The template `example.json` is located in the PE library source:

```bash
cp libs/pe/pe/interface/example.json \
   build-ninja-release/applications/q2p1_bench_sedimentation/
```

### Verify Configuration

Check the file was copied:

```bash
ls -lh build-ninja-release/applications/q2p1_bench_sedimentation/example.json
```

Display contents to verify parameters:

```bash
cat build-ninja-release/applications/q2p1_bench_sedimentation/example.json
```

Expected contents include (example values for the benchmark):

```json
{
    "timesteps_": 20000,
    "stepsize_": 0.002,
    "benchRadius_": 0.0075,
    "particleDensity_": 1120.0,
    "fluidViscosity_": 0.058,
    "fluidDensity_": 960.0,
    "gravity_": [0.0, 0.0, -9.81],
    ...
}
```

**Note:** `libs/pe/pe/interface/example.json` now contains **all JSON-supported fields with default values**
(from `SimulationConfig::SimulationConfig()`). Edit this file for the sedimentation benchmark.

**Physical parameters for Nylon sphere in silicone oil:**
- `benchRadius_`: 0.0075 m (sphere radius = 7.5 mm, diameter = 15 mm)
- `particleDensity_`: 1120.0 kg/m³ (Nylon)
- `fluidViscosity_`: 0.058 Pa·s = 58 mPa·s (silicone oil, ~58× water viscosity)
- `fluidDensity_`: 960.0 kg/m³ (silicone oil)
- `gravity_`: [0.0, 0.0, -9.81] m/s² (standard gravity, negative z-direction)

**Benchmark characteristics (from publication):**
- Density ratio: ρ_particle/ρ_fluid = 1120/960 ≈ 1.167
- Terminal velocity: **v_terminal ≈ 0.128 m/s**
- Reynolds number: **Re ≈ 31.9** (intermediate regime, not Stokes flow)
- Flow regime: Transition between Stokes and inertial regimes

## 6) Ensure `_data/MG.dat` exists in runtime directory

`MG.dat` is now staged automatically by CMake for this benchmark. Verify it exists:

```bash
ls build-ninja-release/applications/q2p1_bench_sedimentation/_data/MG.dat
```

## 7) Understanding JSON and Eigen Dependencies

### JSON is Required

**Critical:** `-DUSE_JSON=ON` is **mandatory** for PE Serial Mode. The setup functions load runtime configuration from `example.json`, and without JSON support, all parameters remain at default values.

CMake's `FetchContent` automatically downloads nlohmann/json from GitHub:
- Source: `https://github.com/nlohmann/json.git`
- Downloaded to: `build-ninja-release/_deps/json-src/`
- This is header-only and works reliably in most network environments

### Eigen (Optional)

Eigen support can be disabled with `-DEIGEN=OFF` if you encounter network issues or don't need it for your simulation:
- Source (when enabled): `https://gitlab.com/libeigen/eigen.git`
- May fail in restricted network environments
- Not required for basic sphere sedimentation benchmarks

### Troubleshooting FetchContent Network Issues

If you encounter `Could not resolve host` errors during configure:

1. **Check network/proxy settings** - CMake FetchContent needs internet access
2. **Use pre-downloaded dependencies** - Advanced: vendor the dependencies locally
3. **Disable optional components** - Use `-DEIGEN=OFF` if Eigen fetch fails

In this guide's tested environment (RHEL 9.7 with standard network access), JSON FetchContent succeeded without issues.

## 8) Build Summary and Verification Checklist

This section provides a quick verification checklist for the completed build.

### Verification Commands

From repository root:

```bash
# 1. Verify executable built
ls -lh build-ninja-release/applications/q2p1_bench_sedimentation/q2p1_bench_sedimentation
# Expected: ~41 MB executable

# 2. Verify JSON support enabled in PE
grep "HAVE_JSON" build-ninja-release/libs/pe/config.h
# Expected: #define HAVE_JSON 1

# 3. Verify 31 partition files created
ls build-ninja-release/applications/q2p1_bench_sedimentation/_mesh/NEWFAC/sub0001/ | grep -c "GRID00"
# Expected: 31

# 4. Verify example.json staged (then edit for benchmark values)
ls -lh build-ninja-release/applications/q2p1_bench_sedimentation/example.json
grep -E "particleDensity_|fluidViscosity_|fluidDensity_" \
  build-ninja-release/applications/q2p1_bench_sedimentation/example.json
# Expected (after editing example.json for the benchmark):
#   "particleDensity_": 1120.0  (Nylon, kg/m³)
#   "fluidViscosity_": 0.058    (Silicone oil, Pa·s = 58 mPa·s)
#   "fluidDensity_": 960.0      (Silicone oil, kg/m³)

# 5. Verify MG.dat staged
ls -lh build-ninja-release/applications/q2p1_bench_sedimentation/_data/MG.dat

# 6. Verify parameter file references correct mesh
cd build-ninja-release/applications/q2p1_bench_sedimentation
grep "MeshFolder\|ProjectFile" _data/q2p1_param.dat
# Expected:
#   SimPar@MeshFolder = "NEWFAC"
#   SimPar@ProjectFile = "_adc/benchSym/bench.prj"
```

### Build Artifacts Summary

```
build-ninja-release/
├── applications/
│   └── q2p1_bench_sedimentation/
│       ├── q2p1_bench_sedimentation     # 41 MB executable
│       ├── example.json                  # PE configuration
│       ├── libmetis.so                   # 412K (for partitioner)
│       ├── _data/
│       │   ├── q2p1_param.dat            # CFD parameters
│       │   └── MG.dat                    # Multigrid config (auto-staged)
│       ├── _mesh/NEWFAC/sub0001/
│       │   ├── GRID0001.tri ... GRID0031.tri  # 31 partitions
│       │   └── GRID.tri                  # Base grid
│       └── _adc/benchSym/
│           └── bench.prj                 # Project file
├── libs/pe/
│   ├── config.h                          # HAVE_JSON 1 ✓
│   └── lib/libpe.a                       # PE library
└── extern/libraries/metis-4.0.3/Lib/
    └── libmetis.so                       # 412K
```

### Troubleshooting Common Issues

**Issue:** `Could not load the Metis library!`
- **Cause:** `libmetis.so` not in application directory or `LD_LIBRARY_PATH`
- **Fix:** Copy METIS library to application directory (see Section 4)

**Issue:** All parameters at default values despite `example.json`
- **Cause:** PE built without JSON support (`HAVE_JSON` not defined)
- **Fix:** Reconfigure with `-DUSE_JSON=ON` and rebuild (see Section 2)

**Issue:** Partition count mismatch errors at runtime
- **Cause:** Mesh partitioned to different count than parameter file expects
- **Fix:** Re-partition with correct count (31 for this guide)

**Issue:** `MG.dat` missing errors
- **Cause:** Staging step not run or build directory incomplete
- **Fix:** Run `cmake --build build-ninja-release --target q2p1_bench_sedimentation_stage`
  or manually copy from repository (see Section 6)

## 9) SLURM Job Submission and Monitoring

### Create SLURM Submission Script

Create `run_sedimentation.sh` in the application directory:

```bash
#!/bin/bash
#SBATCH --partition=short
#SBATCH --constraint=epyc2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=q2p1_bench_sed

# Load modules explicitly (instead of aliases)
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13
module load openmpi/options/interface/ethernet
module load openmpi/4.1.6

# Change to job submission directory
cd $SLURM_SUBMIT_DIR

# Print environment info for debugging
echo "=== Job started at $(date) ==="
echo "Working directory: $(pwd)"
echo "Hostname: $(hostname)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_NTASKS: $SLURM_NTASKS"
echo ""

# Verify critical files exist
echo "=== Verifying files ==="
ls -lh q2p1_bench_sedimentation example.json _data/MG.dat
echo ""

# Run the simulation
echo "=== Starting simulation at $(date) ==="
mpirun -np 32 ./q2p1_bench_sedimentation

# Report completion
echo ""
echo "=== Job completed at $(date) ==="
```

**Key improvements over alias-based scripts:**
- Explicit `module load` commands (reliable in batch sessions)
- Debug output for troubleshooting
- File verification before execution
- No dependency on ~/.bashrc aliases

### Submit the Job

**Human workflow:**

```bash
cd build-ninja-release/applications/q2p1_bench_sedimentation
sbatch run_sedimentation.sh
```

Expected output: `Submitted batch job 100380`

**Agent workflow:**

```bash
bash -c '
  cd /path/to/build-ninja-release/applications/q2p1_bench_sedimentation
  JOB_ID=$(sbatch --parsable run_sedimentation.sh)
  echo "Submitted job: $JOB_ID"
'
```

### Monitor Job Status

**Check queue status:**

```bash
squeue -u $USER
```

Output example:
```
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST
100380     short q2p1_ben rmuenste  R       1:10      1 worldgames
```

**Job states:**
- `PD` = Pending (waiting for resources)
- `R` = Running
- `CG` = Completing
- `CD` = Completed

**Detailed job information:**

```bash
scontrol show job 100380
```

Key fields to check:
- `JobState`: Current state
- `NumCPUs`: Should be 32
- `NumNodes`: Should be 1
- `Partition`: Should be `short`
- `Features`: Should include `epyc2`
- `WorkDir`: Application directory path

### Monitor Simulation Output

**Real-time output (while job runs):**

```bash
tail -f slurm-100380.out
```

**Check for completion:**

```bash
grep "Job completed" slurm-100380.out
```

**Extract key physics results:**

```bash
# Terminal velocity and forces
grep "Force:\|BenchForce:\|Velocity:" slurm-100380.out | tail -50

# Particle position over time
grep "Position:" slurm-100380.out

# Check for PE Serial Mode
grep "SERIAL PE mode" slurm-100380.out
```

### Expected Simulation Output

**PE Serial Mode confirmation:**
```
Force calculation: SERIAL PE mode
Force with SED BENCH settings!
```

**Typical physics output:**
```
Position:  0.000  0.000  0.126     (m)
Velocity:  0.000  0.000 -0.035     (m/s, negative = settling)
Forces:    0.000  0.000  0.002     (N)
```

**Terminal velocity (benchmark publication):**
- Expected: v_terminal ≈ 0.128 m/s
- Reynolds number: Re ≈ 31.9
- Settling time: ~1.2 seconds (for 0.15 m height)

### Agent Monitoring Loop

For automated workflows, agents can monitor until completion:

```bash
bash -c '
  JOB_ID=100380
  echo "Monitoring job $JOB_ID..."

  # Wait for job to complete
  while squeue -j $JOB_ID 2>/dev/null | grep -q $JOB_ID; do
    echo "Job still running... ($(date +%H:%M:%S))"
    sleep 30
  done

  echo "Job completed at $(date)"

  # Extract results
  echo "=== Terminal velocity ==="
  grep "Velocity:" slurm-${JOB_ID}.out | tail -5

  echo "=== Final forces ==="
  grep "BenchForce:" slurm-${JOB_ID}.out | tail -5
'
```

### Troubleshooting SLURM Jobs

**Job fails immediately:**
- Check: `cat slurm-<jobid>.out` for error messages
- Common: Module loading failed → verify module names
- Common: Files missing → check `ls` output in SLURM log

**Job stuck in PD (Pending) state:**
- Check: `squeue -j <jobid> -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"`
- Common reasons: No nodes available, partition constraints not met

**Simulation crashes:**
- Check for MPI errors: `grep -i "mpi\|error\|fail" slurm-<jobid>.out`
- Check memory usage: Job may need more than 3G/CPU
- Check partition logs: `grep "GRID00" slurm-<jobid>.out` (verify mesh loaded)

### Job Completion Checklist

After job completes successfully:

```bash
# 1. Verify job ran to completion
grep "Job completed" slurm-100380.out

# 2. Check final particle state
tail -100 slurm-100380.out | grep -E "Position:|Velocity:"

# 3. Verify PE Serial Mode was active
grep "SERIAL PE mode" slurm-100380.out

# 4. Extract benchmark forces
grep "BenchForce:" slurm-100380.out > forces_timeseries.dat

# 5. Check for errors
grep -i "error\|fail\|abort" slurm-100380.out
```

---

**Guide 02 Part 1 Complete!**

This guide covered:
- ✅ Environment setup (human vs agent workflows)
- ✅ Configuration with JSON support (required!)
- ✅ Building application and METIS
- ✅ Mesh partitioning (31 subdomains)
- ✅ JSON configuration staging
- ✅ SLURM job submission and monitoring

**Next:** Part 2 will cover result analysis, validation against benchmark publication, and post-processing visualization.
