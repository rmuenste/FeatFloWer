# FeatFloWer Automated Test System

This document describes the current architecture of the `featflower-test`
runner. For installation and command examples, see
[featflower_test_usage_guide.md](featflower_test_usage_guide.md).

## Scope

The test system provides a Python CLI for reproducible benchmark setup, build,
partition, execution, metric extraction, baseline comparison, and result
reporting. It lives in `tools/featflower_test/` and complements the repository's
CTest coverage; it does not replace CTest for all targets.

The implemented commands are:

```text
featflower-test validate <definition.yaml>
featflower-test run <definition.yaml> [--slurm] [--dry-run]
featflower-test status <run-id>
featflower-test compare <run-id>
featflower-test report <run-id> --format text|json
```

## Repository Layout

```text
tools/featflower_test/
├── src/featflower_test/
│   ├── cli.py
│   ├── config.py
│   ├── models.py
│   ├── results.py
│   ├── comparison.py
│   ├── report.py
│   ├── parsers/
│   ├── runners/
│   └── stages/
├── testcases/
│   ├── definitions/
│   └── baselines/
└── tests/
```

The included benchmark definitions are:

- `q2p1_fc_ext_cylinder.yaml`
- `q2p1_bench_sedimentation_pe_serial.yaml`

## Execution Pipeline

A `run` executes these stages:

1. Load and validate the YAML definition.
2. Validate setup requirements and submodule configuration.
3. Run ordered CMake configure and build steps.
4. Partition the mesh when a partition stage is configured.
5. Launch locally with MPI or submit through SLURM.
6. Extract configured metrics from produced files.
7. Compare metrics with the referenced baseline.
8. Persist stage logs, metadata, metrics, and reports.

`--dry-run` resolves the pipeline and prints commands without compiling,
partitioning, or launching the solver.

## YAML Model

The supported schema version is `0.2`. A definition contains:

- `id`, `name`, `suite`, `priority`, and `enabled`
- `setup` for modules and submodule handling
- ordered `build.steps`
- `run.levels`, optional `run.partition`, `run.launch`, and `run.slurm`
- `outputs` and `metrics`
- `references.baseline`

The source of truth for concrete syntax is the checked-in definitions under
`tools/featflower_test/testcases/definitions/`.

## Runners

### Local

The local runner executes the configured MPI command directly and records its
combined output and exit status.

### SLURM

The SLURM runner generates a batch script, submits it with `sbatch`, and tracks
the resulting job through scheduler commands. Module setup and launch
parameters come from the test definition.

## Metrics and Baselines

The implemented parser is `keyword_columns`. It:

1. selects a line containing a configured keyword,
2. splits the line into columns,
3. converts Fortran `D` or `E` notation,
4. stores named values, and
5. compares them with absolute tolerances from the test definition.

Baselines are versioned under:

```text
tools/featflower_test/testcases/baselines/
```

Missing baselines generate a validation warning and cause comparison to be
skipped rather than silently creating a new reference.

## Result Storage

Runs are stored below the configured results directory. The result store keeps
run metadata, per-stage status and logs, extracted metrics, comparison results,
and generated text or JSON reports. The `status`, `compare`, and `report`
commands operate on this stored state.

## Tests

The Python unit tests cover:

- configuration validation,
- CLI behavior,
- Fortran number parsing,
- keyword-column extraction,
- metric comparison, and
- parameter-file patching.

Run them from `tools/featflower_test/` with:

```bash
pytest
```

## Remaining Work

The following capabilities are not yet part of the implemented CLI:

- suite discovery and a `list` command,
- concurrent multi-test submission,
- job-array support,
- legacy JSON migration,
- ParaView-derived metrics and image comparison,
- baseline promotion,
- remote submission adapters,
- broader benchmark migration, and
- a formal stable schema beyond version `0.2`.

CTest and the YAML runner should remain in parallel until migrated benchmark
coverage and operational reliability justify retiring any legacy path.
