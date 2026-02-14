# FeatFloWer Automated Test System Plan V2 (Backend-First, Migration-Aware)

This is the next planning iteration for FeatFloWer's automated test system. It
integrates the review feedback in `docs/md_docs/automated_test_system_review.md`
and keeps the implementation grounded in existing repository realities.

The first implementation target remains:

- `docs/md_docs/guide_01_q2p1_fc_ext_cylinder_benchmark_from_scratch.md`
- `docs/md_docs/guide_02_q2p1_bench_sedimentation_pe_serial_from_scratch.md`
- `docs/md_docs/guide_03_q2p1_sse_tse_gendie_from_scratch.md`
- `build_guide.md`

## Scope and Principles

**Goals (Phase 1)**
- Automate benchmark build/run/compare on RHEL 9.7 through SLURM.
- Keep definitions declarative and versioned in-repo.
- Store all outputs in a dashboard-ready, machine-readable format.
- Replace legacy dashboard-coupled test execution with a maintainable backend.
- Explicitly replace CTest/CDash-based test reporting with the new system.

**Non-goals (Phase 1)**
- No mandatory web service layer yet.
- No full schema freeze yet; schema stays intentionally minimal while core
  feature coverage is validated.

## Existing Infrastructure and Migration Strategy

The repository already contains working but outdated automation:

- JSON test definitions in `applications/*/tests/test-*.json` (17 files).
- SLURM orchestrator in `tools/dashboard/q2p1_ctest_start.py`.
- CTest/CDash integration in CMake (legacy path to retire).
- App-specific launcher scripts.

These are treated as **migration input**, not as target architecture.

### Migration Plan

1. Keep legacy jobs runnable during transition.
2. Introduce new YAML-based runner in parallel for selected tests.
3. Add a converter for legacy JSON test files to new YAML with manual review.
4. Migrate guide 01 + guide 02 first, then remaining tests by suite.
5. Remove legacy dashboard-specific execution path only after feature parity on:
   SLURM submission, level loops, extraction, baselines, and result reporting.
6. Decommission CTest/CDash submission/reporting path after migration cutover.

## Architecture (Phase 1)

```
+-------------------------+      +--------------------------+
| Test Definitions (YAML) | ---> | Orchestrator CLI         |
| + Baselines/References  |      | (validation + planning)  |
+-------------------------+      +------------+-------------+
                                              |
                                              v
                                   +----------+-----------+
                                   | SLURM Runner         |
                                   | (sbatch/sacct/squeue)|
                                   +----------+-----------+
                                              |
                                              v
                                   +----------+-----------+
                                   | Collectors/Parsers   |
                                   | logs + solution data |
                                   +----------+-----------+
                                              |
                                              v
                                   +----------+-----------+
                                   | Results Store        |
                                   | filesystem + JSON    |
                                   +----------------------+
```

## YAML Model: Minimal but Concrete

The schema is intentionally not "final" yet. Still, the runner must validate a
minimal required set and reject invalid tests deterministically.

### Required fields (Phase 1)

- `schema_version`
- `id`
- `name`
- `setup`
- `build.steps`
- `run.levels`
- `run.launch`
- `metrics`
- `references`

### Example (Phase 1 target shape)

```yaml
schema_version: 0.2
id: q2p1_fc_ext_cylinder
name: Q2P1 FC_EXT Cylinder
suite: smoke
priority: normal
enabled: true

setup:
  modules:
    - gcc/latest-v13
    - openmpi/options/interface/ethernet
    - openmpi/4.1.6
  git:
    submodules:
      mode: strict
      update_recursive: true

build:
  workspace: build-release
  cache_policy: reuse_if_compatible
  steps:
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
    - kind: cmake_build
      build_dir: build-release
      targets: [q2p1_fc_ext, metis]
      jobs: 8

run:
  levels:
    start: 2
    count: 1
    parameter_patch:
      file_in: _adc/2D_FAC/q2p1_param_2D.dat
      file_out: _data/q2p1_param.dat
      key: SimPar@MaxMeshLevel
  partition:
    tool: python
    command: "python3 tools/PyPartitioner.py 3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj"
    expected_partition_count: 3
  launch:
    workdir: build-release/applications/q2p1_fc_ext
    mpi_ranks: 4
    command: "./q2p1_fc_ext"
  slurm:
    partition: med
    nodes: 1
    ntasks: 4
    cpus_per_task: 1
    time: "00:30:00"

outputs:
  logs: [run.log, build.log]
  artifacts:
    - _data/prot.txt
    - _data/Statistics.txt

metrics:
  - id: drag_lift
    parser: keyword_columns
    file: _data/prot.txt
    keyword: "Force acting"
    columns: {drag: 7, lift: 8}
    occurrence: last
    numeric_format: fortran_d_or_e
    compare:
      type: tolerance
      tolerance: {drag: 1.0e-4, lift: 1.0e-4}

references:
  baseline: testcases/baselines/q2p1_fc_ext_cylinder.yaml
```

### Parser types in Phase 1

- `keyword_columns` (primary; mirrors proven legacy behavior)
- `regex_scalar` (fallback for irregular lines)
- `paraview_python` (structured post-processing from solution outputs)
- `image_compare` (reference image comparison from ParaView-produced outputs)

Free-form shell parsing is allowed only as explicit fallback with opt-in, not
as default extraction mode.

## Setup and Build Specification

### Environment and data-root setup (mandatory)

Before build/run, runner must validate and record:

1. `Q2P1_MESH_DIR` is set (for guide-01/guide-02 style cases).
2. `_adc` links resolve and required project file exists.
3. Required module environment is loaded in every batch shell.

This matches guide behavior where `_adc` links depend on configure-time mesh
directory setup.

### Submodule handling (mandatory)

Submodules must be handled explicitly and reproducibly:

1. `git submodule sync --recursive`
2. `git submodule update --init --recursive`
3. Record each submodule SHA in run metadata.
4. Fail fast on mismatched/unavailable submodule state (`CONFIGURATION_ERROR`).

### Ordered CMake steps (mandatory)

Build is modeled as ordered `steps` (not a single flat option list), so cases
like PE serial mode can be expressed safely (single-step or two-step configure
is allowed, but final option state is strictly validated):

1. Configure with `USE_PE=ON`.
2. Re-configure enabling `USE_PE_SERIAL_MODE=ON`.
3. Build target set.

For `q2p1_bench_sedimentation`, orchestration should support the staged target
path from guide 02:

- `cmake --build <build-dir> --target q2p1_bench_sedimentation_stage`

This target bundles build + runtime staging + partition preparation and reduces
manual staging failures.

### Build caching policy

Phase 1 default: `reuse_if_compatible`:

- Reuse build dir only if compiler/MPI/CMake options match.
- Otherwise create fresh build dir.
- Record cache decision and cache key in metadata.

## Metrics, References, and Validation

### Metric extraction sources

- Logs (`prot.txt`, run logs, `Statistics.txt`)
- Solution-derived data (via ParaView Python pipeline)
- Optional reference images produced from solution fields

Guide-aligned extraction defaults:

- `q2p1_fc_ext`: primary force extraction from `_data/prot.txt` line
  containing `Force acting`.
- `q2p1_bench_sedimentation`: extract from SLURM/application logs
  (`BenchForce:`, `Velocity:`, `Position:`) and optional benchmark series files.

### Reference data model

References are versioned in-repo under `testcases/baselines/` and can include:

- Scalar baselines (drag/lift/force/velocity/time)
- Curve/table baselines (time series)
- Image baselines for view-based checks

Each reference declares comparison method and tolerance/threshold.

### Example reference fragment

```yaml
metrics:
  drag_lift:
    drag: 5.579535
    lift: 0.010618
    tolerance:
      drag_abs: 1.0e-4
      lift_abs: 1.0e-4
images:
  pressure_slice:
    file: refs/pressure_slice.png
    compare: ssim
    min_score: 0.995
```

Baseline updates are explicit via CLI promotion command (no implicit promotion
on passing runs).

## Error Model, Categorization, and Traceability

All failures are persisted with category, stage, and trace pointers.

### Status categories

- `SUCCESS`
- `CONFIGURATION_ERROR`
- `SUBMODULE_ERROR`
- `BUILD_ERROR`
- `RUN_ERROR`
- `METRICS_ERROR`
- `REFERENCE_MISMATCH`
- `SLURM_ERROR`
- `TIMEOUT`
- `INFRA_ERROR`

Categories map to user-facing groups:

- configuration: `CONFIGURATION_ERROR`, `SUBMODULE_ERROR`
- compilation/build: `BUILD_ERROR`
- runtime/execution: `RUN_ERROR`, `SLURM_ERROR`, `TIMEOUT`
- postprocessing/validation: `METRICS_ERROR`, `REFERENCE_MISMATCH`
- infrastructure: `INFRA_ERROR`

### Stored error metadata (minimum)

- `category`
- `stage` (`setup|build|run|collect|compare`)
- `exit_code`
- `slurm_job_id` and SLURM state/reason when available
- `log_path`
- `first_error_line` (if detectable)
- `timestamp`

This supports dashboard trace-back to exact logs and failing stage.

Pre-run hard failures from guide workflows are explicitly captured as
`CONFIGURATION_ERROR`:

- missing `_data/MG.dat` in runtime directory
- missing/invalid `_adc/.../*.prj` reference
- partition count mismatch between generated `GRID*.tri` and launch layout
- missing `example.json` for PE serial mode runs
- PE serial mode requested without `USE_JSON=ON`

## SLURM Concurrency Model

Parallel execution is a first-class requirement:

- Submit independent tests in parallel through SLURM.
- Support job-array mode for homogeneous level runs.
- Respect suite-level `max_in_flight` limits to avoid oversubscription.
- Use dependency chains only when tests explicitly declare dependencies.

Scheduler behavior is delegated to SLURM; orchestrator handles submission,
tracking, and aggregation.

## External Submission Adapter: `mini-sr`

The repository `mini-sr` is a strong candidate as a remote submission adapter
for resource-constrained environments where test execution cannot run locally.

### What is directly useful

- SSH-based remote staging and submission workflow (`auto_sim_submit.py`):
  local case prep, remote case directory creation, script generation, and
  scheduler submission.
- Scheduler abstraction for SLURM/PBS with per-case scheduler parameters.
- Practical config features already proven useful in production workflows:
  `additional_files`, `template_dir` override, `start_commands`, `env_setup`,
  and optional local metadata export with remote case/job mapping.
- Real and mock integration tests, including optional real-cluster smoke tests
  and a Docker-based SLURM test environment.

### Integration boundaries (recommended)

- Keep FeatFloWer test orchestration, metric extraction, baseline comparison,
  and result storage in this new framework.
- Use `mini-sr` patterns/components only for the **remote submission backend**
  (staging + `sbatch`/`qsub` submission) when backend mode is `remote`.
- Treat existing JSON master-config format as an interchange format initially;
  YAML frontends can be added without architecture changes.

### Gaps to close before production integration

- Add job lifecycle tracking (queued/running/completed/failed) and SLURM reason
  capture; current flow is submission-focused.
- Enforce strict exit propagation in generated scripts (`set -euo pipefail` and
  robust pipeline handling) so failed start commands are not masked.
- Add artifact retrieval/pullback from remote case directories into the
  framework result store (`prot.txt`, `Statistics.txt`, ParaView outputs, logs).
- Harden command/path quoting for shell-safe execution from config data.
- Extend job-id parsing and metadata behavior to non-SLURM schedulers where
  required.
- Add preflight validation parity with FeatFloWer needs (partition count,
  `_data/MG.dat`, `example.json`, mesh/project file checks).

### Integration plan for `mini-sr`

1. Add a `remote_submitter=minisr` backend mode in the new runner.
2. Generate a per-run master submission file from FeatFloWer test YAML.
3. Invoke adapter and capture returned mapping: test-id -> remote-case-dir/job-id.
4. Track job completion, fetch artifacts, then run normal compare/report stages.
5. Incrementally upstream/port hardening improvements listed above.

## Result Storage Layout (Phase 1)

```
results/
  runs/
    <run-id>/
      metadata.json
      plan.json
      summary.json
      tests/
        <test-id>/
          test_metadata.json
          stages/
            setup.log
            build.log
            partition.log
            run-level-2.log
          artifacts/
            prot.txt
            Statistics.txt
            slurm/
              slurm-<jobid>.out
              slurm-<jobid>.err
            partition_meta.json
            paraview/
              metrics.json
              images/
          metrics.json
          compare.json
          errors.json
```

`metadata.json` includes commit SHA, branch, dirty flag, submodule SHAs, loaded
modules, compiler/MPI versions, and all SLURM job IDs.

For reproducibility, store the resolved CMake option set, including PE/JSON
flags and any FetchContent source overrides used in restricted environments.

## Guide 01 / Guide 02 / Guide 03 Test Blueprint

### Guide 01 (`q2p1_fc_ext`)

- Build: `q2p1_fc_ext` and `metis` in release mode with MPI wrappers and
  `USE_HYPRE=ON`.
- Partition: 3 partitions with `PyPartitioner`, validate `GRID0001..0003`.
- Run: `mpirun -np 4 ./q2p1_fc_ext`.
- Compare: drag/lift from latest `Force acting` line in `_data/prot.txt`.

### Guide 02 (`q2p1_bench_sedimentation`, PE serial)

- Configure/build with `USE_PE=ON`, `USE_PE_SERIAL_MODE=ON`, `USE_JSON=ON`.
- Prefer stage target: `q2p1_bench_sedimentation_stage`.
- Validate runtime prerequisites: `example.json`, `_data/MG.dat`, partition set
  (`31` worker partitions for `mpirun -np 32` default guide case).
- Run under SLURM with explicit module loads inside job script.
- Compare using extracted force/velocity/position metrics and optional
  ParaView-derived outputs.

### Guide 03 (`q2p1_sse` for SSE/TSE, gendie context)

- Configure/build with `USE_CGAL=ON`, `USE_HYPRE=ON`, `USE_PE=OFF`,
  `ENABLE_FBM_ACCELERATION=OFF`.
- Build and stage with `q2p1_sse_stage` so runtime dependencies are present
  (`s3d_mesher`, `libmetis.so`, `_ianus`, `_data/MG.dat`, scripts).
- Run under SLURM using `python3 ./e3d_start.py -n 32 -f _ianus/TSE/Conv -a 0 --short-test`.
- Require explicit ROMIO/OpenMPI I/O environment in submit script:
  `OMPI_MCA_io=romio321`, `ROMIO_CB_BUFFER_SIZE=16777216`,
  `ROMIO_DS_WRITE=enable`.
- Compare via log-state checks and selected convergence/output metrics.

## CTest/CDash Replacement

CTest/CDash is explicitly out of scope as an execution/reporting backend for
the new system and will be replaced.

Cutover policy:

- During migration, CTest/CDash may remain temporarily for legacy jobs only.
- New YAML-runner results are the authoritative source for pass/fail and trends.
- After defined migration completion, CDash submission is disabled and the old
  dashboard/CTest glue code is retired.

## CLI Contract (Phase 1)

```bash
featflower-test list
featflower-test validate <test-or-suite>
featflower-test run <test-or-suite> [--slurm] [--max-in-flight N]
featflower-test status <run-id>
featflower-test compare <run-id>
featflower-test report <run-id> --format json
featflower-test promote-baseline <run-id> --test <id>
featflower-test migrate-json <legacy-json> --out <yaml>
```

## Security and Robustness Rules

- Prefer structured parsers over ad-hoc shell pipelines.
- Restrict executable commands to declared stage fields.
- Store exact command lines and environment in metadata.
- Treat missing artifacts and parser mismatch as typed errors, not silent skips.

## Phased Implementation Plan

1. Implement runner skeleton (`validate`, `run`, `status`) and storage layout.
2. Implement setup/build/run stages including strict submodule checks.
3. Implement `keyword_columns` parser and tolerance comparison.
4. Implement guide 01 and guide 02 YAML definitions and baselines.
5. Implement typed error model + traceable failure reporting.
6. Add parallel SLURM submission (`max_in_flight`, optional job arrays).
7. Add ParaView extraction path (`paraview_python`) and optional image compare.
8. Add legacy JSON migration helper and migrate remaining test definitions.

## Acceptance Criteria for Phase 1

- Guide 01, guide 02, and guide 03 run fully through new system on SLURM.
- Submodule state is verified and logged for every run.
- Failures are categorized and traceable to logs/stages.
- Scalar comparisons and at least one ParaView-derived comparison are supported.
- Multiple tests can be submitted in parallel through SLURM.
- Legacy JSON-to-YAML migration path exists and is documented.
- CTest/CDash is no longer required for benchmark test reporting.
