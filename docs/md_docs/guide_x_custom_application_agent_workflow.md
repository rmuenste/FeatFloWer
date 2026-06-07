# FeatFloWer Guide X (WIP): Custom Application Build and Deployment - Agent Workflow

**Status:** Work In Progress
**Audience:** AI agents and automated workflows
**Example Application:** `q2p1_hashgrid_test`

This guide documents the workflow for configuring, building, staging, and deploying custom FeatFloWer applications with PE Serial Mode using an AI agent-based approach. It captures lessons learned from successful autonomous builds and serves as a template for similar applications.

---

## Overview

This guide demonstrates building a custom CFD application (`q2p1_hashgrid_test`) with:
- **PE Serial Mode** for rigid body physics (< 20 large particles)
- **Dual build configurations** (baseline vs accelerated) for performance comparison
- **Agent-friendly workflow** with module loading in each command
- **63-domain mesh partitioning** for 64-core nodes (63 workers + 1 master)

### When to Use This Guide

Use this workflow when:
- Building custom Q2/P1 applications with PE coupling
- Particles are large relative to domain size (span multiple domains)
- Testing new algorithms with baseline/accelerated comparisons
- Automating builds in non-interactive shell environments
- Running on HPC systems with module-based environments

---

## Prerequisites

### System Requirements
- **OS:** RHEL 9.7 or compatible Linux distribution
- **Build Tools:** CMake ≥ 3.18, Ninja build system
- **Compilers:** GCC 13+ with MPI wrappers (mpicc, mpicxx, mpifort)
- **Python:** Python 3 for mesh partitioner
- **Network:** Outbound access for FetchContent (JSON, Eigen dependencies)

### Environment Modules (RHEL 9.7)
```bash
module load gcc/latest-v13
module load openmpi/options/interface/ethernet
module load openmpi/4.1.6
```

### Application Resources Checklist

Before building, ensure you have:
- [ ] **Mesh files** (`.tri` format and boundary `.par` files)
- [ ] **Project file** (`.prj` defining mesh metadata)
- [ ] **Parameter file** (`q2p1_param.dat` for CFD configuration)
- [ ] **FullC0ntact start folder** (required if using FullC0ntact library)
  - Location: `FullC0ntact/start/` directory
  - Initialized during `init_fc_rigid_body()` call
  - Contains initialization data for rigid body interface

---

## Part 1: Dual Build Strategy (Baseline + Accelerated)

For algorithm validation and performance testing, we build **two versions**:

### Build Configuration Matrix

| Build Type | Directory | HashGrid | Verification | Purpose |
|------------|-----------|----------|--------------|---------|
| **Baseline** | `build-ninja-baseline` | OFF | OFF | Reference implementation (O(N) linear search) |
| **Accelerated** | `build-ninja-accelerated` | ON | ON | Optimized implementation with timing comparison |

This dual-build approach allows:
- ✅ Correctness validation (accelerated vs baseline)
- ✅ Performance measurement (speedup factor)
- ✅ Regression testing after code changes

---

## Part 2: Configure Baseline Build

### Step 1: CMake Configuration (Baseline)

**Command Pattern:**
```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  cd /path/to/repository
  cmake -S . -B build-ninja-baseline -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_APPLICATIONS=ON \
    -DUSE_PE=ON \
    -DUSE_PE_SERIAL_MODE=ON \
    -DUSE_JSON=ON \
    -DUSE_ACCELERATED_POINT_QUERY=OFF \
    -DVERIFY_HASHGRID=OFF \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_Fortran_COMPILER=mpifort
' 2>&1 | tee /tmp/configure_baseline.log | tail -50
```

**Expected Output Indicators:**
```
-- PE SERIAL MODE ENABLED
-- Accelerated point query is DISABLED (using linear O(N) search)
-- HashGrid verification is DISABLED
-- Found MPI_C: .../mpicc (found version "3.1")
-- nlohmann/json will be downloaded to: ...
-- Configuring done
-- Generating done
```

**Configuration Flags Explained:**

| Flag | Value | Purpose |
|------|-------|---------|
| `USE_PE` | ON | Enable physics engine library |
| `USE_PE_SERIAL_MODE` | ON | Serial PE per domain (large particles) |
| `USE_JSON` | ON | **MANDATORY** - Load PE config from `example.json` |
| `USE_ACCELERATED_POINT_QUERY` | OFF | Use baseline O(N) algorithm |
| `VERIFY_HASHGRID` | OFF | No verification overhead |

**Critical:** Without `-DUSE_JSON=ON`, PE configuration remains at default values even if `example.json` exists.

---

## Part 3: Build Baseline Executable

### Step 2: Compile Application

**Command Pattern:**
```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  cd /path/to/repository
  cmake --build build-ninja-baseline --target YOUR_APP_NAME -- -j8
' > /tmp/build_baseline.log 2>&1 && \
  echo "✅ Baseline build completed" || \
  echo "❌ Baseline build failed"
```

**Replace:** `YOUR_APP_NAME` with your application target (e.g., `q2p1_hashgrid_test`)

**Verification:**
```bash
# Check executable was created
ls -lh build-ninja-baseline/applications/YOUR_APP_NAME/YOUR_APP_NAME

# Expected: ~40-50 MB executable
# Example: -rwxr-xr-x 1 user group 46M Feb 6 11:44 q2p1_hashgrid_test

# Scan build log for errors (should return no critical errors)
grep -i "error\|undefined reference\|failed" /tmp/build_baseline.log | \
  grep -v "error.f-pp.f\|error_category.hpp"
```

---

## Part 4: Stage Resources and Partition Mesh (Baseline)

### Step 3: Build METIS

```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  cd /path/to/repository
  cmake --build build-ninja-baseline --target metis -- -j8
'

# Verify METIS library
ls -lh build-ninja-baseline/extern/libraries/metis-4.0.3/Lib/libmetis.so
# Expected: ~336-412K shared library
```

### Step 4: Stage Application Resources

**Directory Structure to Create:**
```
build-ninja-baseline/applications/YOUR_APP_NAME/
├── YOUR_APP_NAME              # Executable
├── _data/
│   └── q2p1_param.dat         # CFD parameters
├── _mesh/NEWFAC/sub0001/
│   └── GRID.tri               # Base mesh (before partitioning)
├── your_mesh_case/            # Mesh source directory
│   ├── mesh_file.tri
│   ├── file.prj
│   └── *.par                  # Boundary parameter files
├── libmetis.so                # For partitioner
└── FullC0ntact/start/         # If using FullC0ntact library (optional)
```

**Staging Commands:**
```bash
APP_DIR=build-ninja-baseline/applications/YOUR_APP_NAME

# Create directory structure
mkdir -p ${APP_DIR}/_data
mkdir -p ${APP_DIR}/_mesh/NEWFAC/sub0001

# Copy parameter file
cp /path/to/resources/q2p1_param.dat ${APP_DIR}/_data/

# Copy mesh source directory
cp -r /path/to/resources/mesh_case ${APP_DIR}/

# Copy METIS library
cp build-ninja-baseline/extern/libraries/metis-4.0.3/Lib/libmetis.so ${APP_DIR}/

# Copy FullC0ntact start folder (if needed)
# Required if app_init.f90 calls init_fc_rigid_body()
if [ -d "FullC0ntact/start" ]; then
  cp -r FullC0ntact/start ${APP_DIR}/
fi

# Copy base mesh to partition directory
cp ${APP_DIR}/mesh_case/mesh_file.tri ${APP_DIR}/_mesh/NEWFAC/sub0001/GRID.tri
```

### Step 5: Partition Mesh

**For 64-core node:** Use 63 partitions (63 workers + 1 master/coordinator)

```bash
cd build-ninja-baseline/applications/YOUR_APP_NAME

# Run PyPartitioner
# Syntax: PyPartitioner.py <num_partitions> <min_level> <max_level> <folder_name> <project_file>
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
python3 ../../../tools/PyPartitioner.py 63 1 1 NEWFAC mesh_case/file.prj

# Verify partition count
ls _mesh/NEWFAC/sub0001/ | grep -c "GRID00"
# Expected output: 63
```

**Partition Files Created:**
- `GRID0001.tri` through `GRID0063.tri` (63 partitions)
- `GRID.tri` (base mesh, not counted as partition)

**Common Partitioning Configurations:**

| Node Type | Total Cores | Partitions | Master Ranks |
|-----------|-------------|------------|--------------|
| Small | 32 | 31 | 1 |
| Standard | 64 | 63 | 1 |
| Large | 128 | 127 | 1 |

**PyPartitioner Arguments:**
- `63` - Number of partitions
- `1` - Minimum mesh refinement level
- `1` - Maximum mesh refinement level
- `NEWFAC` - Mesh folder name (matches `MeshFolder` in param.dat)
- `mesh_case/file.prj` - Project file path

---

## Part 5: Configure and Build Accelerated Version

### Step 6: CMake Configuration (Accelerated)

**Key Differences from Baseline:**
- `USE_ACCELERATED_POINT_QUERY=ON` (enable HashGrid)
- `VERIFY_HASHGRID=ON` (enable verification and timing)

```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  cd /path/to/repository
  cmake -S . -B build-ninja-accelerated -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_APPLICATIONS=ON \
    -DUSE_PE=ON \
    -DUSE_PE_SERIAL_MODE=ON \
    -DUSE_JSON=ON \
    -DUSE_ACCELERATED_POINT_QUERY=ON \
    -DVERIFY_HASHGRID=ON \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_Fortran_COMPILER=mpifort
' 2>&1 | tee /tmp/configure_accelerated.log | tail -50
```

**Expected Output Indicators:**
```
-- PE SERIAL MODE ENABLED
-- Accelerated point query (HashGrid spatial hashing) is ENABLED
-- HashGrid verification is ENABLED - will compare accelerated vs baseline
   Note: This significantly impacts performance and should only be used for testing
```

### Step 7: Build Accelerated Executable

```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
  cd /path/to/repository
  cmake --build build-ninja-accelerated --target YOUR_APP_NAME -- -j8
' > /tmp/build_accelerated.log 2>&1 && \
  echo "✅ Accelerated build completed" || \
  echo "❌ Accelerated build failed"
```

### Step 8: Stage Resources (Accelerated)

**Additional file needed:** `example.json` for PE configuration

```bash
APP_DIR=build-ninja-accelerated/applications/YOUR_APP_NAME

# Create directory structure
mkdir -p ${APP_DIR}/_data
mkdir -p ${APP_DIR}/_mesh/NEWFAC/sub0001

# Copy all resources (same as baseline)
cp /path/to/resources/q2p1_param.dat ${APP_DIR}/_data/
cp -r /path/to/resources/mesh_case ${APP_DIR}/
cp build-ninja-accelerated/extern/libraries/metis-4.0.3/Lib/libmetis.so ${APP_DIR}/

# Copy PE JSON configuration
cp libs/pe/pe/interface/example.json ${APP_DIR}/

# Copy FullC0ntact start folder if needed
if [ -d "FullC0ntact/start" ]; then
  cp -r FullC0ntact/start ${APP_DIR}/
fi

# Copy base mesh
cp ${APP_DIR}/mesh_case/mesh_file.tri ${APP_DIR}/_mesh/NEWFAC/sub0001/GRID.tri

# Partition (same as baseline)
cd ${APP_DIR}
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
python3 ../../../tools/PyPartitioner.py 63 1 1 NEWFAC mesh_case/file.prj
```

---

## Part 6: SLURM Job Submission

### Step 9: Create Baseline SLURM Script

**File:** `build-ninja-baseline/applications/YOUR_APP_NAME/run_baseline.sh`

```bash
#!/bin/bash
#SBATCH --partition=short
#SBATCH --constraint=epyc2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=YOUR_APP_baseline

# Load modules explicitly (critical for batch jobs)
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13
module load openmpi/options/interface/ethernet
module load openmpi/4.1.6

# Change to job submission directory
cd $SLURM_SUBMIT_DIR

# Print environment info
echo "=== Job started at $(date) ==="
echo "Working directory: $(pwd)"
echo "Hostname: $(hostname)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_NTASKS: $SLURM_NTASKS"
echo ""

# Verify critical files
echo "=== Verifying files ==="
ls -lh YOUR_APP_NAME _data/q2p1_param.dat
ls _mesh/NEWFAC/sub0001/ | grep -c "GRID00"
echo ""

# Run the simulation (baseline - linear search)
echo "=== Starting BASELINE simulation at $(date) ==="
mpirun -np 64 ./YOUR_APP_NAME

# Report completion
echo ""
echo "=== Job completed at $(date) ==="
```

### Step 10: Create Accelerated SLURM Script

**File:** `build-ninja-accelerated/applications/YOUR_APP_NAME/run_accelerated.sh`

```bash
#!/bin/bash
#SBATCH --partition=short
#SBATCH --constraint=epyc2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=YOUR_APP_accel

# Load modules explicitly
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13
module load openmpi/options/interface/ethernet
module load openmpi/4.1.6

# Change to job submission directory
cd $SLURM_SUBMIT_DIR

# Print environment info
echo "=== Job started at $(date) ==="
echo "Working directory: $(pwd)"
echo "Hostname: $(hostname)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_NTASKS: $SLURM_NTASKS"
echo ""

# Verify critical files
echo "=== Verifying files ==="
ls -lh YOUR_APP_NAME example.json _data/q2p1_param.dat
ls _mesh/NEWFAC/sub0001/ | grep -c "GRID00"
echo ""

# Run the simulation (accelerated with verification)
echo "=== Starting ACCELERATED simulation at $(date) ==="
echo "Verification enabled - will compare HashGrid vs baseline"
mpirun -np 64 ./YOUR_APP_NAME

# Report completion
echo ""
echo "=== Job completed at $(date) ==="
```

### Step 11: Submit Jobs

```bash
# Submit baseline
cd build-ninja-baseline/applications/YOUR_APP_NAME
sbatch run_baseline.sh
# Output: Submitted batch job 100456

# Submit accelerated
cd build-ninja-accelerated/applications/YOUR_APP_NAME
sbatch run_accelerated.sh
# Output: Submitted batch job 100457

# Monitor jobs
squeue -u $USER
```

---

## Part 7: Verification and Output Analysis

### Expected Output (Baseline)

```
==================================================
FBM Alpha Field Computation Performance
==================================================
Wall-clock time:    2.456789 s
Total point queries:123456
Performance:        50234.1 queries/sec

> Total dofs inside: 45678
> Dofs per Particle: 234
```

### Expected Output (Accelerated with Verification)

```
==================================================
FBM Alpha Field Computation Performance
==================================================
Wall-clock time:    0.234567 s
Total point queries:123456
Performance:        526234.7 queries/sec

==================================================
HashGrid Verification & Performance Summary
==================================================

--- Correctness ---
Total queries:         123456
Result:                ✓ Perfect match (0 mismatches)

--- Performance Comparison ---
Accelerated (HashGrid):      0.123456 s
  Per query:                 1.000 μs

Baseline (Linear):           2.345678 s
  Per query:                19.000 μs

Speedup:                    19.00x
==================================================

> Total dofs inside: 45678
> Dofs per Particle: 234
```

### Key Metrics to Extract

**From SLURM output:**
```bash
# Baseline performance
grep "Wall-clock time" slurm-100456.out
grep "Performance:" slurm-100456.out

# Accelerated verification results
grep "Perfect match\|mismatches" slurm-100457.out
grep "Speedup:" slurm-100457.out

# Compare dofs inside (should be identical)
grep "Total dofs inside" slurm-100456.out
grep "Total dofs inside" slurm-100457.out
```

---

## Part 8: Complete Verification Checklist

### Build Verification

Run these checks after building:

```bash
# 1. Baseline executable
ls -lh build-ninja-baseline/applications/YOUR_APP_NAME/YOUR_APP_NAME
# Expected: ~40-50 MB

# 2. Accelerated executable
ls -lh build-ninja-accelerated/applications/YOUR_APP_NAME/YOUR_APP_NAME
# Expected: ~40-50 MB (similar size)

# 3. Baseline partitions
ls build-ninja-baseline/applications/YOUR_APP_NAME/_mesh/NEWFAC/sub0001/ | grep -c "GRID00"
# Expected: 63

# 4. Accelerated partitions
ls build-ninja-accelerated/applications/YOUR_APP_NAME/_mesh/NEWFAC/sub0001/ | grep -c "GRID00"
# Expected: 63

# 5. JSON configuration present (accelerated only)
ls -lh build-ninja-accelerated/applications/YOUR_APP_NAME/example.json
# Expected: ~700-800 bytes

# 6. Parameter file in both builds
ls -lh build-ninja-*/applications/YOUR_APP_NAME/_data/q2p1_param.dat
# Expected: exists in both

# 7. Verify CMake flags in config headers
grep "USE_ACCELERATED_POINT_QUERY\|VERIFY_HASHGRID" \
  build-ninja-baseline/libs/pe/config.h
grep "USE_ACCELERATED_POINT_QUERY\|VERIFY_HASHGRID" \
  build-ninja-accelerated/libs/pe/config.h
```

### Expected CMake Config Headers

**Baseline** (`build-ninja-baseline/libs/pe/config.h`):
```cpp
/* #undef USE_ACCELERATED_POINT_QUERY */
/* #undef VERIFY_HASHGRID */
#define HAVE_JSON 1
```

**Accelerated** (`build-ninja-accelerated/libs/pe/config.h`):
```cpp
#define USE_ACCELERATED_POINT_QUERY 1
#define VERIFY_HASHGRID 1
#define HAVE_JSON 1
```

---

## Part 9: Troubleshooting Common Issues

### Configuration Issues

**Issue:** `HAVE_JSON` not defined in `config.h`

**Symptoms:**
```
All parameters remain at default values despite example.json
```

**Solution:**
```bash
# Reconfigure with JSON explicitly enabled
cmake -S . -B build-ninja-accelerated \
  -DUSE_JSON=ON \
  [other flags...]
```

---

**Issue:** PE Serial Mode not enabled

**Symptoms:**
```
-- PE standard mode enabled
(Missing "PE SERIAL MODE ENABLED" message)
```

**Solution:**
```bash
# Two-step configuration for serial mode
cmake -S . -B build -DUSE_PE=ON ...
cmake -S . -B build -DUSE_PE_SERIAL_MODE=ON
```

---

**Issue:** MPI compilers not found

**Symptoms:**
```
Could NOT find MPI_C
Could NOT find MPI_Fortran
```

**Solution:**
```bash
# Verify modules are loaded
module list | grep -i "gcc\|openmpi"

# Verify MPI wrappers are in PATH
which mpicc mpicxx mpifort

# Explicitly specify compilers
cmake ... \
  -DCMAKE_C_COMPILER=$(which mpicc) \
  -DCMAKE_CXX_COMPILER=$(which mpicxx) \
  -DCMAKE_Fortran_COMPILER=$(which mpifort)
```

---

### Build Issues

**Issue:** METIS library not found during partitioning

**Symptoms:**
```
Could not load the Metis library!
```

**Solution:**
```bash
# Verify METIS was built
ls -lh build-ninja-*/extern/libraries/metis-4.0.3/Lib/libmetis.so

# Copy to application directory
cp build-ninja-*/extern/libraries/metis-4.0.3/Lib/libmetis.so \
   build-ninja-*/applications/YOUR_APP_NAME/

# Use LD_LIBRARY_PATH when running partitioner
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" python3 ../../../tools/PyPartitioner.py ...
```

---

**Issue:** FetchContent fails (network/proxy)

**Symptoms:**
```
Could not resolve host: github.com
Failed to download nlohmann/json
```

**Solution:**
```bash
# Option 1: Check network connectivity
curl -I https://github.com/nlohmann/json.git

# Option 2: Use pre-downloaded dependencies
cmake ... \
  -DFETCHCONTENT_SOURCE_DIR_JSON=/path/to/json-src

# Option 3: Disable optional components
cmake ... -DEIGEN=OFF
```

---

### Runtime Issues

**Issue:** Partition count mismatch

**Symptoms:**
```
Expected 63 partitions but found 31
Mesh loading failed
```

**Solution:**
```bash
# Re-partition with correct count
cd build-ninja-*/applications/YOUR_APP_NAME
rm _mesh/NEWFAC/sub0001/GRID00*.tri
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
python3 ../../../tools/PyPartitioner.py 63 1 1 NEWFAC mesh_case/file.prj

# Verify count
ls _mesh/NEWFAC/sub0001/ | grep -c "GRID00"
```

---

**Issue:** FullC0ntact initialization fails

**Symptoms:**
```
Error in init_fc_rigid_body()
Cannot find start directory
```

**Solution:**
```bash
# Verify FullC0ntact/start exists
ls -la FullC0ntact/start/

# Copy to application directory
cp -r FullC0ntact/start build-ninja-*/applications/YOUR_APP_NAME/
```

---

**Issue:** SLURM job fails immediately

**Symptoms:**
```
Job completes in < 1 second
Exit code: 1 or 127
```

**Debugging:**
```bash
# Check SLURM output
cat slurm-<jobid>.out

# Common causes:
# 1. Module loading failed
grep -i "module: command not found" slurm-<jobid>.out

# 2. Missing files
grep -i "no such file" slurm-<jobid>.out

# 3. MPI errors
grep -i "mpi\|error\|fail" slurm-<jobid>.out
```

---

## Part 10: Adaptation Template for New Applications

### Creating a New Custom Application

**Step-by-step adaptation:**

1. **Replace application name throughout:**
   ```bash
   # In all commands, replace:
   YOUR_APP_NAME → your_new_app_name
   ```

2. **Update resource paths:**
   ```bash
   # Update these paths to your resources:
   /path/to/resources/q2p1_param.dat → /path/to/your/params.dat
   /path/to/resources/mesh_case → /path/to/your/mesh_dir
   ```

3. **Adjust partition count:**
   ```bash
   # For different node sizes:
   # 32-core: 31 partitions
   # 64-core: 63 partitions
   # 128-core: 127 partitions
   python3 PyPartitioner.py <partitions> 1 1 NEWFAC mesh_case/file.prj
   ```

4. **Modify SLURM resources:**
   ```bash
   # In SLURM script, adjust:
   #SBATCH --ntasks-per-node=64  # Match partition count + 1
   #SBATCH --time=02:00:00       # Adjust for problem size
   #SBATCH --mem-per-cpu=3G      # Adjust for memory needs
   ```

5. **Add application-specific CMake flags:**
   ```bash
   # Example: Add custom flags to cmake command
   cmake -S . -B build-ninja-baseline \
     -DUSE_YOUR_FEATURE=ON \
     [standard flags...]
   ```

### Checklist for New Application

- [ ] Application target exists in `applications/CMakeLists.txt`
- [ ] Application has `app_init.f90` with proper initialization
- [ ] Mesh files and project file are prepared
- [ ] Parameter file (`q2p1_param.dat`) configured for mesh
- [ ] FullC0ntact start folder available (if using rigid bodies)
- [ ] Partition count matches available cores
- [ ] SLURM script resources match computational needs

---

## Part 11: Performance Baseline Data (Example: q2p1_hashgrid_test)

### Test Configuration
- **Mesh:** Unit cube (25 refinement)
- **Partitions:** 63 domains
- **Particles:** 1000 spheres (10×10×10 grid)
- **Queries:** ~123,456 point queries per alpha field computation

### Expected Performance (Reference)

| Metric | Baseline (Linear) | Accelerated (HashGrid) | Speedup |
|--------|-------------------|------------------------|---------|
| **Wall-clock time** | ~2.5 s | ~0.13 s | ~19x |
| **Queries/second** | ~50,000 | ~950,000 | ~19x |
| **Per-query time** | ~20 μs | ~1 μs | ~20x |

**Note:** Actual performance varies with:
- Particle count (O(N) vs O(1) scaling)
- Mesh resolution (number of queries)
- Hardware (CPU, memory bandwidth)
- MPI configuration (domain decomposition)

---

## Part 12: Summary and Best Practices

### Agent Workflow Best Practices

1. **Always use module loading pattern:**
   ```bash
   bash -c '
     source /etc/profile.d/modules.sh
     module purge
     module load [required modules]
     [your command]
   '
   ```

2. **Log all output to files:**
   ```bash
   command > /tmp/output.log 2>&1
   ```

3. **Verify at each step:**
   - After configure: Check for expected CMake output
   - After build: Verify executable exists
   - After staging: Confirm all files present
   - After partition: Count partition files

4. **Use tee for important commands:**
   ```bash
   command | tee /tmp/output.log | tail -50
   ```

5. **Check errors explicitly:**
   ```bash
   command && echo "✅ Success" || echo "❌ Failed"
   ```

### File Organization

```
project_root/
├── build-ninja-baseline/
│   └── applications/YOUR_APP_NAME/
│       ├── YOUR_APP_NAME
│       ├── run_baseline.sh
│       └── [resources]
├── build-ninja-accelerated/
│   └── applications/YOUR_APP_NAME/
│       ├── YOUR_APP_NAME
│       ├── run_accelerated.sh
│       ├── example.json
│       └── [resources]
└── docs/
    └── md_docs/
        └── this_guide.md
```

### Automation Scripts

For fully autonomous builds, create a master script:

```bash
#!/bin/bash
# auto_build_and_deploy.sh

set -e  # Exit on error

APP_NAME="your_app_name"
PARTITIONS=63

# Configure and build baseline
bash -c "source /etc/profile.d/modules.sh && module load [modules] && \
  cmake -S . -B build-ninja-baseline [baseline flags] && \
  cmake --build build-ninja-baseline --target ${APP_NAME}"

# Stage and partition baseline
bash -c "[staging commands for baseline]"

# Configure and build accelerated
bash -c "source /etc/profile.d/modules.sh && module load [modules] && \
  cmake -S . -B build-ninja-accelerated [accelerated flags] && \
  cmake --build build-ninja-accelerated --target ${APP_NAME}"

# Stage and partition accelerated
bash -c "[staging commands for accelerated]"

echo "✅ Build and deployment complete"
```

---

## Appendix A: Quick Reference Commands

### Build Commands (One-liners)

**Configure Baseline:**
```bash
bash -c 'source /etc/profile.d/modules.sh && module purge && module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6 && cmake -S . -B build-ninja-baseline -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON -DUSE_PE=ON -DUSE_PE_SERIAL_MODE=ON -DUSE_JSON=ON -DUSE_ACCELERATED_POINT_QUERY=OFF -DVERIFY_HASHGRID=OFF -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort'
```

**Build Application:**
```bash
bash -c 'source /etc/profile.d/modules.sh && module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6 && cmake --build build-ninja-baseline --target YOUR_APP_NAME -- -j8'
```

**Partition Mesh:**
```bash
cd build-ninja-baseline/applications/YOUR_APP_NAME && LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" python3 ../../../tools/PyPartitioner.py 63 1 1 NEWFAC mesh_case/file.prj
```

### Verification One-liners

```bash
# Verify all builds complete
ls -lh build-ninja-*/applications/YOUR_APP_NAME/YOUR_APP_NAME

# Verify partition counts
for dir in build-ninja-*/applications/YOUR_APP_NAME; do echo "$dir: $(ls $dir/_mesh/NEWFAC/sub0001/ | grep -c GRID00) partitions"; done

# Compare results
diff <(grep "Total dofs inside" slurm-<baseline-id>.out) <(grep "Total dofs inside" slurm-<accel-id>.out)
```

---

## Appendix B: Example Resource Files

### Example q2p1_param.dat (minimal)

```fortran
----------- Simulation Parameters --------------
SimPar@MeshFolder = "NEWFAC"
SimPar@SubMeshNumber = 1
SimPar@ProjectFile = 'mesh_case/file.prj'
SimPar@MinMeshLevel = 1
SimPar@MaxMeshLevel = 2
SimPar@TimeStep = 0.001d0
SimPar@MaxSimTime = 1.0d0
```

### Example file.prj

```
'mesh_file.tri'
```

---

## Appendix C: Changelog

**2026-02-06:** Initial version based on successful q2p1_hashgrid_test build
- Documented dual-build workflow (baseline + accelerated)
- Added agent-friendly command patterns
- Included verification procedures
- Added troubleshooting section

---

**Guide Status:** Work In Progress (WIP)
**Last Updated:** 2026-02-06
**Maintainer:** AI Agent Workflow Documentation
**Related Guides:** Guide 02 (Sedimentation Benchmark)
