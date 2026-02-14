# FeatFloWer Test System: Usage Guide

This guide covers installation, configuration, and practical use of the
`featflower-test` CLI tool introduced in Phase 1a of the automated test system.

- tool location: `tools/featflower_test/`
- first supported test case: Guide 01 (`q2p1_fc_ext` cylinder benchmark)
- supported runners: local (mpirun) and SLURM (sbatch)

## 1) Prerequisites

- Python 3.8+
- PyYAML >= 5.0 (installed automatically)
- FeatFloWer source tree with submodules
- `Q2P1_MESH_DIR` environment variable pointing to the mesh repository
- For local runs: MPI (OpenMPI or equivalent)
- For SLURM runs: access to a SLURM cluster with `sbatch`/`sacct`

## 2) Installation

From the repository root:

```bash
cd tools/featflower_test
pip install -e ".[dev]"
```

This installs:
- `featflower-test` CLI command
- `pytest` for running the unit test suite

Verify the installation:

```bash
featflower-test --version
```

If your `$PATH` does not include `~/.local/bin`, you can also invoke via:

```bash
python -m featflower_test --version
```

## 3) Quick Start: Validate a Test Definition

The fastest way to check that a YAML test definition is well-formed:

```bash
featflower-test validate \
  tools/featflower_test/testcases/definitions/q2p1_fc_ext_cylinder.yaml
```

Expected output:

```
OK: q2p1_fc_ext_cylinder (Q2P1 FC_EXT Cylinder Benchmark)
  schema_version: 0.2
  suite: smoke
  build steps: 2
  metrics: 1
  baseline: testcases/baselines/q2p1_fc_ext_cylinder.yaml
```

Exit code is `0` on success, `1` on validation failure.

### What validation checks

- `schema_version` is `"0.2"` (currently the only supported version)
- All required top-level keys are present: `schema_version`, `id`, `name`,
  `setup`, `build`, `run`, `metrics`, `references`
- Every build step `kind` is one of `cmake_configure` or `cmake_build`
- A warning is printed if the referenced baseline file does not exist

## 4) Dry Run: Walk Through the Pipeline Without Executing

A dry run exercises the full pipeline (validate, setup, build, partition, run)
but does not execute any compilation, partitioning, or simulation commands.
Use this to verify that the YAML definition, environment, and directory
structure are consistent before committing to a real run.

```bash
featflower-test run \
  tools/featflower_test/testcases/definitions/q2p1_fc_ext_cylinder.yaml \
  --dry-run
```

Expected output:

```
Run ID: 20260210-091500-a1b2c3d4
Test:   q2p1_fc_ext_cylinder (Q2P1 FC_EXT Cylinder Benchmark)

Stage: setup ... passed
Stage: build ... passed
Stage: partition ... passed
Stage: run ... passed
Stage: metrics ... failed (1 metrics)
Stage: compare ... skipped (no baseline)
```

Notes:

- **Setup stage** still runs for real: it checks `Q2P1_MESH_DIR` and syncs
  submodules.  This ensures the environment is validated even in dry-run mode.
- **Build/partition/run** stages log the commands they *would* execute but skip
  actual execution.
- **Metrics extraction** will report an error because no simulation output
  (`_data/prot.txt`) exists yet. This is expected.

### Result directory

Every run (including dry runs) creates a result directory:

```
results/runs/<run-id>/tests/q2p1_fc_ext_cylinder/
  stages/
    setup.log
    build.log
    partition.log
    run-level-2.log
  metadata.json
  metrics.json
```

You can inspect the build log to see exactly what cmake commands would be issued:

```bash
cat results/runs/<run-id>/tests/q2p1_fc_ext_cylinder/stages/build.log
```

## 5) Full Local Run (Guide 01 Cylinder Benchmark)

This section walks through a complete local run of the Guide 01 test case.

### 5.1) Environment setup

```bash
export Q2P1_MESH_DIR=/absolute/path/to/your/mesh-repository
```

### 5.2) Execute the test

```bash
featflower-test run \
  tools/featflower_test/testcases/definitions/q2p1_fc_ext_cylinder.yaml
```

This executes the following pipeline stages in order:

1. **Setup** -- Validates `Q2P1_MESH_DIR`, runs `git submodule sync` and
   `git submodule update --init --recursive`, records submodule SHAs.
2. **Build** -- Runs `cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
   -DUSE_HYPRE=ON ...`, then builds targets `q2p1_fc_ext` and `metis`.
3. **Partition** -- Invokes `PyPartitioner.py 3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj`
   from the application's build directory, verifies 3 GRID*.tri files are created.
4. **Run** -- Patches `SimPar@MaxMeshLevel = 2` into the parameter file, then
   runs `mpirun -np 4 ./q2p1_fc_ext`.
5. **Metrics** -- Extracts drag and lift from the last `BenchForce:` line in
   `_data/prot.txt` (columns 2 and 3, 0-indexed after `split()`).
6. **Compare** -- Compares extracted values against the baseline
   (`testcases/baselines/q2p1_fc_ext_cylinder.yaml`) with tolerance `1.0e-4`
   per metric.

### 5.3) Expected output on success

```
Run ID: 20260210-093000-e5f6a7b8
Test:   q2p1_fc_ext_cylinder (Q2P1 FC_EXT Cylinder Benchmark)

Stage: setup ... passed
Stage: build ... passed
Stage: partition ... passed
Stage: run ... passed
Stage: metrics ... passed (1 metrics)
Stage: compare ... pass

Run:    20260210-093000-e5f6a7b8
Test:   q2p1_fc_ext_cylinder
Commit: 2aa93c5f...
Branch: master

Stages:
  setup: passed (1.2s)
  build: passed (245.3s)
  partition: passed (3.1s)
  run-level-2: passed (180.6s)

Metrics:
  drag_lift: drag=5.579535, lift=0.010618

Comparison: pass
  [PASS] drag_lift.drag: measured=5.579535 ref=5.579535 tol=0.0001 diff=0.0
  [PASS] drag_lift.lift: measured=0.010618 ref=0.010618 tol=0.0001 diff=0.0
```

Exit code is `0` on pass, `1` on comparison failure or pipeline error.

### 5.4) Skip the build stage

If you have already built the application and just want to re-run from
partitioning onwards:

```bash
featflower-test run \
  tools/featflower_test/testcases/definitions/q2p1_fc_ext_cylinder.yaml \
  --skip-build
```

## 6) SLURM Execution

To submit the simulation step via SLURM instead of running locally with
`mpirun`:

```bash
featflower-test run \
  tools/featflower_test/testcases/definitions/q2p1_fc_ext_cylinder.yaml \
  --slurm
```

The SLURM configuration is taken from the `run.slurm` section in the YAML:

```yaml
slurm:
  partition: med
  nodes: 1
  ntasks: 4
  cpus_per_task: 1
  time: "00:30:00"
```

### What happens under the hood

1. A job script is generated at
   `build-release/applications/q2p1_fc_ext/job_level_2.sh` containing the
   `#SBATCH` directives and `mpirun` command.
2. The script is submitted via `sbatch`.
3. The tool polls `sacct --format=State --jobs=<id> --noheader --parsable2`
   every 10 seconds until the job reaches a terminal state.
4. Terminal states are mapped to error categories:
   - `COMPLETED` --> success
   - `FAILED`, `NODE_FAIL`, `OUT_OF_MEMORY` --> `RUN_ERROR`
   - `TIMEOUT` --> `TIMEOUT`
   - `CANCELLED`, `PREEMPTED` --> `SLURM_ERROR`

### Customizing SLURM parameters

Edit the `slurm:` block in the YAML definition. All standard SLURM directives
are supported: `partition`, `nodes`, `ntasks`, `cpus_per_task`, `time`,
`constraint`, `mem_per_cpu`.

## 7) Inspecting Results After a Run

### 7.1) Check run status

```bash
featflower-test status <run-id> --results-dir results
```

Example output:

```
Run:  20260210-093000-e5f6a7b8
Test: q2p1_fc_ext_cylinder
  setup: passed
  build: passed
  partition: passed
  run-level-2: passed
```

### 7.2) Generate a text report

```bash
featflower-test report <run-id> --format text
```

### 7.3) Generate a JSON report

```bash
featflower-test report <run-id> --format json
```

The JSON report includes full metadata (commit SHA, branch, submodule SHAs,
cmake options), extracted metric values, and per-metric comparison details.
This is suitable for machine consumption, dashboards, or CI artifact storage.

### 7.4) Result directory layout

```
results/runs/<run-id>/tests/<test-id>/
  stages/
    setup.log         # submodule sync output, SHA listing
    build.log         # cmake configure + build stdout
    partition.log     # partitioner output, file verification
    run-level-2.log   # mpirun/SLURM simulation output
  artifacts/          # collected output files (future use)
  metadata.json       # full run metadata
  metrics.json        # extracted metric values
  compare.json        # comparison results
  errors.json         # structured errors (if any)
```

### 7.5) Reading individual stage logs

Stage logs are plain text and contain the exact commands that were executed:

```bash
cat results/runs/<run-id>/tests/q2p1_fc_ext_cylinder/stages/build.log
```

Example content:

```
=== cmake configure ===
$ cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON -DUSE_HYPRE=ON ... -S /path/to/FeatFloWer -B /path/to/build-release

...
--- exit code: 0 ---

=== cmake build: q2p1_fc_ext ===
$ cmake --build /path/to/build-release --target q2p1_fc_ext -- -j8

...
--- exit code: 0 ---
```

## 8) Understanding the YAML Test Definition

The test definition is the central configuration file. Here is an annotated
walkthrough of `testcases/definitions/q2p1_fc_ext_cylinder.yaml`:

```yaml
# Schema version (required, must be "0.2")
schema_version: "0.2"

# Unique test identifier and human-readable name
id: q2p1_fc_ext_cylinder
name: Q2P1 FC_EXT Cylinder Benchmark

# Classification (informational)
suite: smoke
priority: normal
enabled: true

# --- Setup stage ---
setup:
  # Environment modules to load (SLURM job scripts)
  modules:
    - gcc/latest-v13
    - openmpi/options/interface/ethernet
    - openmpi/4.1.6
  # Submodule handling
  git:
    submodules:
      mode: strict           # fail on submodule mismatch
      update_recursive: true # --recursive flag

# --- Build stage ---
build:
  workspace: build-release
  cache_policy: reuse_if_compatible
  steps:
    # Step 1: cmake configure
    - kind: cmake_configure
      build_dir: build-release
      source_dir: .
      options:
        - CMAKE_BUILD_TYPE=Release
        - BUILD_APPLICATIONS=ON
        - USE_HYPRE=ON
        - CMAKE_C_COMPILER=mpicc
        - CMAKE_CXX_COMPILER=mpicxx
        - CMAKE_Fortran_COMPILER=mpifort
    # Step 2: cmake build (both targets)
    - kind: cmake_build
      build_dir: build-release
      targets: [q2p1_fc_ext, metis]
      jobs: 8

# --- Run stage ---
run:
  # Mesh level iteration
  levels:
    start: 2        # starting mesh level
    count: 1        # number of levels to run
    parameter_patch:
      file_in: _adc/2D_FAC/q2p1_param_2D.dat   # source param file
      file_out: _data/q2p1_param.dat             # destination (in workdir)
      key: SimPar@MaxMeshLevel                   # key to patch

  # Mesh partitioning
  partition:
    tool: python
    # {repo_root} is resolved at runtime to the repository root
    command: "python3 {repo_root}/tools/PyPartitioner.py 3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj"
    expected_partition_count: 3  # verified after partitioning

  # Simulation launch
  launch:
    workdir: build-release/applications/q2p1_fc_ext
    mpi_ranks: 4
    command: "./q2p1_fc_ext"

  # SLURM parameters (used when --slurm flag is given)
  slurm:
    partition: med
    nodes: 1
    ntasks: 4
    cpus_per_task: 1
    time: "00:30:00"

# --- Metric extraction ---
metrics:
  - id: drag_lift
    parser: keyword_columns      # extraction method
    file: _data/prot.txt         # relative to workdir
    keyword: "BenchForce:"       # line selection keyword
    columns: {drag: 2, lift: 3}  # 0-indexed columns after split()
    occurrence: last             # use last matching line
    numeric_format: fortran_d_or_e
    compare:
      type: tolerance
      tolerance: {drag: 1.0e-4, lift: 1.0e-4}

# --- Baseline reference ---
references:
  baseline: testcases/baselines/q2p1_fc_ext_cylinder.yaml
```

### Understanding `keyword` + `columns`

The `keyword_columns` parser works as follows:

1. Read the file and find all lines containing the `keyword` string.
2. Select a line based on `occurrence` (`"last"` = last matching line,
   `"first"` = first matching line).
3. Split the selected line by whitespace: `line.split()`.
4. Extract values at the specified 0-indexed column positions.
5. Parse each value through the Fortran float normalizer (handles `D`/`d`
   exponents like `0.7350D-09`).

Example: for a line like

```
BenchForce: 0.1000E+01 5.579535 0.010618 0.2100E+01 0.0062E-01
```

After `split()`, the tokens are:

```
[0] BenchForce:
[1] 0.1000E+01
[2] 5.579535      <-- drag (column 2)
[3] 0.010618      <-- lift (column 3)
[4] 0.2100E+01
[5] 0.0062E-01
```

### Understanding the baseline file

The baseline file (`testcases/baselines/q2p1_fc_ext_cylinder.yaml`) stores
reference values keyed by metric ID and value name:

```yaml
metrics:
  drag_lift:
    drag: 5.579535
    lift: 0.010618
```

Comparison checks: `abs(measured - reference) <= tolerance` per value.

## 9) Writing a New Test Definition

To add a new test case:

### 9.1) Create the YAML definition

```bash
cp tools/featflower_test/testcases/definitions/q2p1_fc_ext_cylinder.yaml \
   tools/featflower_test/testcases/definitions/my_new_test.yaml
```

Edit the file, changing:
- `id` and `name`
- Build options and targets
- Partition command and expected count
- Launch command, MPI ranks, and working directory
- Metric keyword, columns, and file path

### 9.2) Validate

```bash
featflower-test validate \
  tools/featflower_test/testcases/definitions/my_new_test.yaml
```

### 9.3) Dry run

```bash
featflower-test run \
  tools/featflower_test/testcases/definitions/my_new_test.yaml \
  --dry-run
```

### 9.4) First real run (to establish baseline)

```bash
featflower-test run \
  tools/featflower_test/testcases/definitions/my_new_test.yaml
```

Without a baseline file, the comparison stage is skipped. Inspect the metrics
output in the report to verify correctness, then create the baseline file:

```bash
featflower-test report <run-id> --format json | \
  python3 -c "
import json, sys, yaml
data = json.load(sys.stdin)
baseline = {'metrics': {}}
for m in data['metrics']:
    mid = m['metric_id']
    if m.get('error'):
        continue
    baseline['metrics'][mid] = {v['name']: v['value'] for v in m['values']}
yaml.dump(baseline, sys.stdout, default_flow_style=False)
" > tools/featflower_test/testcases/baselines/my_new_test.yaml
```

### 9.5) Regression run

Now that the baseline exists, subsequent runs will compare against it:

```bash
featflower-test run \
  tools/featflower_test/testcases/definitions/my_new_test.yaml
```

## 10) Error Handling and Troubleshooting

### Error categories

Every failure is categorized. The category tells you which stage to
investigate:

| Category | Meaning | Where to look |
|----------|---------|---------------|
| `CONFIGURATION_ERROR` | Missing env var, bad YAML | `stages/setup.log`, CLI output |
| `SUBMODULE_ERROR` | Submodule sync/update failed | `stages/setup.log` |
| `BUILD_ERROR` | CMake configure or build failed | `stages/build.log` |
| `PARTITION_ERROR` | Partitioner failed or wrong file count | `stages/partition.log` |
| `RUN_ERROR` | Simulation returned non-zero exit | `stages/run-level-N.log` |
| `METRICS_ERROR` | Metric extraction failed | `metrics.json` |
| `SLURM_ERROR` | sbatch/sacct failure or cancelled job | `stages/run-level-N.log` |
| `TIMEOUT` | Simulation exceeded time limit | `stages/run-level-N.log` |

### Common issues

**Q2P1_MESH_DIR not set:**

```
CONFIGURATION_ERROR: Q2P1_MESH_DIR environment variable is not set
```

Fix: `export Q2P1_MESH_DIR=/path/to/mesh-repo` and re-run.

**Submodule update fails:**

```
SUBMODULE_ERROR: git submodule update failed (rc=1)
```

Check network connectivity and SSH keys. Inspect `stages/setup.log` for
the full git error output.

**Build fails:**

```
BUILD_ERROR: cmake build target 'q2p1_fc_ext' failed (rc=2)
```

Check `stages/build.log` for compiler errors. Common causes: missing MPI
wrappers, unsupported compiler flags, missing dependencies.

**Partition count mismatch:**

```
PARTITION_ERROR: Expected 3 partition files (GRID*.tri), found 0
```

Check that `_adc` symlinks are in place (requires `Q2P1_MESH_DIR` set during
cmake configure) and that `libmetis.so` is discoverable. The tool
automatically sets `LD_LIBRARY_PATH` to include the working directory.

**Keyword not found in output:**

```
METRICS_ERROR: Keyword 'BenchForce:' not found in .../prot.txt
```

The simulation may have failed before producing the expected output.
Check `stages/run-level-N.log` for errors.

## 11) Running the Unit Test Suite

The package includes a comprehensive unit test suite:

```bash
cd tools/featflower_test
pip install -e ".[dev]"
pytest tests/ -v
```

Expected: 41 tests passing across 6 test modules:

| Module | Tests | What it covers |
|--------|-------|----------------|
| `test_fortran_numbers.py` | 11 | Fortran D/E exponent parsing |
| `test_keyword_columns.py` | 6 | Keyword + column metric extraction |
| `test_comparison.py` | 5 | Tolerance comparison engine |
| `test_parameter_patch.py` | 4 | SimPar@MaxMeshLevel patching |
| `test_config.py` | 9 | YAML loading and validation |
| `test_cli.py` | 5 | CLI entry point and commands |

These tests run without MPI, SLURM, or a build environment. They use fixture
files in `tests/fixtures/` to verify parsing, extraction, and comparison
logic in isolation.

## 12) CLI Command Reference

```
featflower-test validate <yaml> [--repo-root DIR]
```
Validate a YAML test definition. Exit 0 on success, 1 on failure.

```
featflower-test run <yaml> [--slurm] [--dry-run] [--skip-build]
                           [--repo-root DIR] [--results-dir DIR]
```
Execute the full pipeline. Options:
- `--slurm`: submit simulation via SLURM instead of local mpirun
- `--dry-run`: walk through stages without executing commands
- `--skip-build`: skip cmake configure and build (reuse existing build)
- `--repo-root`: override auto-detected repository root
- `--results-dir`: override default `results/` directory

```
featflower-test status <run-id> [--test-id ID] [--results-dir DIR]
```
Print stage statuses for a completed or in-progress run.

```
featflower-test compare <run-id> [--test-id ID] [--results-dir DIR]
```
Re-run comparison from stored metrics (placeholder in Phase 1a).

```
featflower-test report <run-id> [--format text|json] [--test-id ID]
                                [--results-dir DIR]
```
Generate a human-readable or machine-readable report.

## 13) Architecture Overview

```
                    q2p1_fc_ext_cylinder.yaml
                              |
                              v
                    +-------------------+
                    |  config.py        |  YAML loading + validation
                    |  (TestDefinition) |
                    +--------+----------+
                             |
                    +--------v----------+
                    |     cli.py        |  Pipeline orchestration
                    |   cmd_run()       |
                    +--------+----------+
                             |
              +--------------+--------------+
              |              |              |
     +--------v---+  +------v-----+  +-----v------+
     | stages/    |  | stages/    |  | stages/    |
     | setup.py   |  | build.py   |  | partition  |
     +------------+  +------------+  +------------+
                             |
                    +--------v----------+
                    | stages/run.py     |  parameter patch + delegation
                    +--------+----------+
                             |
              +--------------+--------------+
              |                             |
     +--------v---------+       +-----------v--------+
     | runners/local.py |       | runners/slurm.py   |
     | (mpirun)         |       | (sbatch + sacct)   |
     +------------------+       +--------------------+
                             |
                    +--------v----------+
                    | parsers/          |  metric extraction
                    | keyword_columns   |
                    +--------+----------+
                             |
                    +--------v----------+
                    | comparison.py     |  tolerance check
                    +--------+----------+
                             |
                    +--------v----------+
                    | results.py        |  JSON persistence
                    | report.py         |  text/JSON output
                    +-------------------+
```

## 14) Relationship to Legacy System

The `featflower-test` tool replaces the following legacy components:

| Legacy | Replacement |
|--------|-------------|
| `tools/dashboard/q2p1_ctest_start.py` | `cli.py` + stages |
| `applications/*/tests/test-*.json` | `testcases/definitions/*.yaml` |
| `getLogEntry()` / `get_col_data()` | `parsers/keyword_columns.py` |
| `moveAndSetLevel()` | `stages/run.py:patch_parameter_file()` |
| `submitAndObserveSync()` | `runners/slurm.py` |
| `generateSlurmScript()` | `runners/slurm.py:_generate_job_script()` |
| CTest/CDash submission | `results.py` + `report.py` |

The extraction logic (keyword search, column indexing, Fortran float parsing)
is intentionally behaviour-compatible with the legacy code to ensure
continuity of metric values during migration.
