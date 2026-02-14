# FeatFloWer Guide 01: Complete From-Scratch Run of the `q2p1_fc_ext` Cylinder Benchmark

This guide is the first in a practical runbook series.
It shows how to configure, compile, partition, run, and validate a **basic FeatFloWer benchmark**:

- application: `q2p1_fc_ext`
- case: 2D flow around a cylinder (`_adc/2D_FAC/2Dbench.prj`)
- quantity of interest: hydrodynamic drag/lift (written to protocol output)

## 1) What You Need Before You Start

- FeatFloWer source tree checked out.
- Submodules available.
- CMake + Ninja + GNU toolchain.
- MPI (OpenMPI or equivalent) and wrappers (`mpicc`, `mpicxx`, `mpifort`).
- Python 3 for partitioning tools.
- Mesh repository available and pointed to by `Q2P1_MESH_DIR`.

### Why `Q2P1_MESH_DIR` matters

During configure, FeatFloWer creates `_adc` links for application run folders (see `cmake/modules/CreateDataDirectories.cmake`).
If `Q2P1_MESH_DIR` is set, `_adc` is symlinked to that mesh repo.
The benchmark project file used in this guide is:

- `_adc/2D_FAC/2Dbench.prj`

## 2) One-Time Environment Setup

From the FeatFloWer repository root:

```bash
export Q2P1_MESH_DIR=/absolute/path/to/your/mesh-repository
```

Optional sanity check:

```bash
echo "$Q2P1_MESH_DIR"
```

## 3) Configure a Release Build (Ninja + HYPRE)

```bash
cd /path/to/FeatFloWer

git submodule update --init --recursive

cmake -S . -B build-release \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_HYPRE=ON \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort
```

Notes:

- `-DUSE_HYPRE=ON` requires working MPI detection.
- Using MPI wrapper compilers avoids common `FindMPI` failures.

## 4) Build the Benchmark Application and METIS

Build the benchmark executable:

```bash
cmake --build build-release --target q2p1_fc_ext -- -j3
```

Build METIS shared library (needed by Python partitioner):

```bash
cmake --build build-release --target metis -- -j3
```

After this, `libmetis.so` should be available at least in:

- `build-release/applications/q2p1_fc_ext/libmetis.so`
- `build-release/extern/libraries/metis-4.0.3/Lib/libmetis.so`

## 5) Verify the Benchmark Run Directory

Go to the benchmark runtime folder produced by CMake:

```bash
cd build-release/applications/q2p1_fc_ext
```

Check key files:

```bash
ls q2p1_fc_ext _data/q2p1_param.dat _adc/2D_FAC/2Dbench.prj
```

If `_adc/...` is missing, re-check `Q2P1_MESH_DIR` and re-configure.

## 6) Partition the Mesh with `PyPartitioner`

For this example, create **3 partitions** into mesh folder `NEWFAC`:

```bash
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
python3 /path/to/FeatFloWer/tools/PyPartitioner.py \
  3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj
```

Parameter meaning:

- `3` = number of partitions
- `1` = METIS method (default/recommended)
- `1` = number of submeshes (no recursive coarse split)
- `NEWFAC` = output folder under `_mesh/`
- project file = `_adc/2D_FAC/2Dbench.prj`

Expected output includes:

- `Calling Metis...`
- `The partitioning was successful!`

Check generated files:

```bash
ls _mesh/NEWFAC/sub0001/GRID000*.tri
```

You should see:

- `GRID0001.tri`
- `GRID0002.tri`
- `GRID0003.tri`

## 7) Ensure Runtime Parameters Match the Partitioned Mesh

Default parameter file is:

- `build-release/applications/q2p1_fc_ext/_data/q2p1_param.dat`

Relevant entries should be:

- `SimPar@MeshFolder = "NEWFAC"`
- `SimPar@SubMeshNumber = 1`
- `SimPar@ProjectFile = "_adc/2D_FAC/2Dbench.prj"`

If needed, edit `_data/q2p1_param.dat` accordingly.

## 8) Run the Benchmark

Run with MPI from the same folder (`build-release/applications/q2p1_fc_ext`).

Important convention used by this app launcher flow:

- solver mesh partitions are typically `NumProcessor - 1`
- example: 3 mesh partitions -> `mpirun -np 4`

Run command (recommended: write full stdout/stderr to a logfile):

```bash
mpirun -np 4 ./q2p1_fc_ext > run_q2p1_fc_ext_np4.log 2>&1
```

Quick check of important runtime lines:

```bash
rg -n "Force:|BenchForce:|error|failed|stopped" run_q2p1_fc_ext_np4.log | tail -n 20
```

Note to self (assistant/automation): avoid terminal/context flooding by logging full output to file and only printing filtered summaries.

## 9) Extract Drag/Lift Information

Protocol file location (from default parameters):

- `_data/prot.txt`

Quick extraction of latest force line:

```bash
grep 'Force acting' _data/prot.txt | tail -1
```

A compact drag/lift extraction pattern used in existing helper scripts:

```bash
grep 'Force acting' _data/prot.txt | tail -1 | awk '{print $7, $8}'
```

## 10) How to Read a Typical Timestep Log Block

Example markers you will see:

- `time: ... | itns: ... | dt: ...` -> current simulation time, step index, and timestep size.
- `FBM computation step` / `FBM time` -> immersed-boundary/FBM coupling work and wall time for this phase.
- `INL ...` block -> nonlinear velocity solver outer iterations with residual/defect information.
- pressure block (`nMGcycPres`, `DefInitPres`, `DefFinalPres`) -> pressure solve convergence.
- `Force:` and `BenchForce:` -> benchmark force coefficients (`C_D`, `C_L`) plus split contributions.

Important number-format note:

- FeatFloWer logs use Fortran exponents:
  - `D` and `E` both denote powers of ten in these outputs.
  - Example: `0.7350D-09 = 7.35e-10`, `0.1940D-09 = 1.94e-10`.

### Interpreting the inner velocity-solver lines

Inside the `INL` section you may see lines like:

```text
   0  0.1281E-08   1.000
   1  0.1255E-20  0.9797E-12 |    50 |    2
```

These are inner multigrid/coarse-solver diagnostics for the current nonlinear stage:

- first column (`0`, `1`) = inner iteration/cycle index
- second value = current inner defect/residual estimate
- third value = normalized ratio vs initial inner defect (`1.000` at start)
- trailing `| 50 | 2` = coarse-solver work counters / metadata (code-specific print format)

In the example above, the inner solve reduces the defect by ~12 orders of magnitude (`~1e-12` ratio), indicating a very strong linear solve in that stage.

### Force line consistency check

For lines of the form:

```text
Force: Time C_D C_L ForceVx ForceVy ForcePx ForcePy
```

you can verify:

- `C_D = ForceVx + ForcePx`
- `C_L = ForceVy + ForcePy`

If both equalities hold (within roundoff), force assembly is internally consistent for that step.

## 11) Runtime Pitfalls, Symptoms, and Fast Checks

This section is intentionally user-facing and developer-facing: these failures are common and should be reduced with better checks and defaults.

### 1) Partitioning failed, but the user did not notice

Partitioning currently assumes basic CLI and Python workflow familiarity. If partitioning silently failed or was skipped, runtime failures can be confusing.

Typical symptoms:

- solver starts but fails early when loading partitioned mesh files
- missing expected files under `_mesh/<MeshFolder>/sub0001/`

Fast checks:

```bash
test -f _mesh/NEWFAC/sub0001/GRID0001.tri && echo "partition files exist"
rg -n "The partitioning was successful|Could not load the Metis library|Traceback|error" partition_*.log -i
```

#### 1a) `libmetis.so` not found at runtime by Python partitioner

Even if METIS is built, dynamic loader search paths may not include it.

Typical symptom:

- `Could not load the Metis library!`

Fast fix:

```bash
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
python3 /path/to/FeatFloWer/tools/PyPartitioner.py 3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj
```

Alternative explicit path:

```bash
LD_LIBRARY_PATH="/path/to/build-release/extern/libraries/metis-4.0.3/Lib:$LD_LIBRARY_PATH"
```

#### 1b) Partition-count mismatch (very common and hard to debug)

Example failure mode:

- mesh was partitioned to 64
- later run is started with settings/process layout corresponding to 32 without re-partitioning

This can produce hard-to-interpret runtime errors because the coarse problem and partitioned subproblems no longer match consistently.

Fast checks:

```bash
ls _mesh/NEWFAC/sub0001/GRID*.tri | wc -l
```

Ensure that partition count in `_mesh/<MeshFolder>/sub0001/GRIDxxxx.tri` matches your intended run setup.

### 2) `_data/MG.dat` missing in runtime folder

Some application/runtime setups may miss `_data/MG.dat`, causing failures during solver initialization.

Fast check:

```bash
test -f _data/MG.dat || echo "MISSING: _data/MG.dat"
```

Fast fix:

```bash
cp /path/to/FeatFloWer/_data/MG.dat _data/
```

Developer note: verify and unify CMake copy behavior for `MG.dat` across all applications/build paths.

### 3) Mesh path under `_adc` referenced, but mesh repository was never configured

If `Q2P1_MESH_DIR` was not set during configure, `_adc` linkage can be missing or wrong.

Typical symptom:

- `_adc/2D_FAC/2Dbench.prj` not found

Fast checks:

```bash
echo "$Q2P1_MESH_DIR"
ls -la _adc
test -f _adc/2D_FAC/2Dbench.prj || echo "MISSING: project file"
```

Fast fix:

- set `Q2P1_MESH_DIR` to mesh-repo root
- re-run CMake configure so symlinks are regenerated

### 4) Build-time MPI detection issue (configure-stage pitfall)

If `-DUSE_HYPRE=ON` is enabled and MPI is not detected correctly, configure may fail.

Use MPI wrappers explicitly:

```bash
-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort
```

### Developer hardening checklist (recommended)

1. Add explicit pre-run validation in launcher scripts for:
   - existence of `_data/MG.dat`
   - readability of `SimPar@ProjectFile`
   - existence of partition files for expected partition count
2. Add a partition-manifest file (e.g. `partition_meta.json`) written by partitioner and checked by launcher.
3. Improve `PyPartitioner` error messages for missing `libmetis.so` with exact path hints.
4. Consider failing fast in application init when partition/coarse-grid consistency checks fail.

## 12) Minimal Reproducible Command Block

```bash
# from repo root
export Q2P1_MESH_DIR=/absolute/path/to/mesh-repo

git submodule update --init --recursive

cmake -S . -B build-release -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DUSE_HYPRE=ON \
  -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort

cmake --build build-release --target q2p1_fc_ext -- -j3
cmake --build build-release --target metis -- -j3

cd build-release/applications/q2p1_fc_ext

LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
python3 /path/to/FeatFloWer/tools/PyPartitioner.py 3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj

mpirun -np 4 ./q2p1_fc_ext > run_q2p1_fc_ext_np4.log 2>&1

rg -n "Force:|BenchForce:|error|failed|stopped" run_q2p1_fc_ext_np4.log | tail -n 20

grep 'Force acting' _data/prot.txt | tail -1
```

## 13) Where to Go Next in the Series

Recommended follow-up guides:

1. Guide 02: parameter sweeps (Reynolds number, viscosity, time step).
2. Guide 03: scaling runs (partition/submesh strategy, MPI placement).
3. Guide 04: postprocessing automation and regression checks.

Additional documentation recommendation:

- Yes, create a dedicated reference file (for example `docs/md_docs/log_output_reference.md`) that explains all major output sections and columns once, then keep this guide focused on the benchmark workflow.
