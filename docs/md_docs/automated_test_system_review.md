# Review: FeatFloWer Automated Test System Design

**Document under review:** `docs/md_docs/automated_test_system.md`
**Commit:** `fa65bfd71122c7079a8330554e3f490828fe3bf3`

---

## Executive Summary

The design document proposes a backend-first automated test system using YAML
definitions, a Python CLI orchestrator, SLURM job submission, and structured
result storage. The overall direction is sound: declarative test definitions,
separation of concerns between definition and execution, and a storage model
designed for future dashboard consumption are all good architectural choices.

However, the document does not sufficiently account for the **existing test
infrastructure** already present in the repository, and several areas need
deeper specification before implementation can begin. This review identifies
what works well, what is missing, and concrete recommendations.

---

## Strengths

### 1. Declarative YAML test definitions

Replacing the current deeply-nested JSON format (e.g.,
`applications/q2p1_fc_ext/tests/test-fac2d.json`) with flat, stage-oriented
YAML is a clear improvement. The existing JSON conflates test configuration
with dashboard visualization concerns (Google API chart schemas, column
definitions, row templates). The proposed YAML cleanly separates build, run,
and metrics stages.

### 2. Backend-first philosophy

Starting with filesystem storage and optional SQLite, then scaling to
PostgreSQL, is pragmatic. The `results/runs/<run-id>/` layout with
`metadata.json` and `metrics.json` is well-structured and easy to query
programmatically.

### 3. Clear SLURM integration plan

The five-stage job script structure (module setup, build, run, capture, exit
status) maps directly to the existing workflow in `q2p1_ctest_start.py` and
the guide documents. Parameterizing SLURM directives per test is the right
approach.

### 4. Metric extraction with tolerance-based comparison

Defining `tolerance` bands for drag/lift regression detection is appropriate
for CFD benchmarks where exact reproducibility is not guaranteed across
compiler versions or MPI configurations.

---

## Issues and Gaps

### Issue 1: No acknowledgment of existing test infrastructure

**Severity: High**

The document reads as a greenfield proposal, but the repository already has
substantial test infrastructure:

- **17 JSON test definition files** across `applications/*/tests/test-*.json`
  covering multiple applications (q2p1_fc_ext, q2p1_bench_sedimentation,
  q2p1_creep, q2p1_drill, q2p1_fac_visco, q2p1_fc2, q2p1_fsi_bench,
  q2p1_rotation, q2p1_span, q2p1_xParticles).
- **`tools/dashboard/q2p1_ctest_start.py`** (~506 lines): a working
  orchestrator that reads JSON configs, generates SLURM scripts, submits via
  `sbatch`, polls with `sacct`, extracts metrics from `prot.txt` and
  `Statistics.txt`, and POSTs results to a dashboard API.
- **CTest integration**: `enable_testing()` in the root CMakeLists.txt,
  `CTestConfig.cmake` pointing to a CDash server, and per-application
  `add_test()` calls.
- **Per-application Python launchers** (e.g.,
  `applications/q2p1_fc_ext/q2p1_fc_ext_start.py`) handling partitioning,
  parameter file modification, and MPI launch.

**Recommendation:** The design document should include a migration strategy:
which parts of the existing system are being replaced, which are being wrapped,
and how the 15+ existing test definitions will be converted. Without this, there
is a risk of building a parallel system that duplicates existing functionality.

### Issue 2: The YAML schema is underspecified

**Severity: Medium**

The example YAML is labeled "conceptual" but several critical details are
missing:

- **No schema definition or validation strategy.** What happens when a
  required field is missing? How are defaults handled?
- **The `metrics.extract.command` field uses a shell pipeline** (`grep ... |
  tail -1`), but there is no specification for how the output is parsed into
  structured fields. The existing system uses keyword matching
  (`benchKeyword`) plus column index arrays (`DataValues.value`) to extract
  specific numbers from output lines. The YAML example loses this precision.
- **No support for multi-level testing.** The existing system iterates over
  mesh refinement levels (`testLevels` in the JSON config, implemented as the
  `range(2, 2 + testLevels)` loop in `q2p1_ctest_start.py`). This is a core
  workflow pattern for CFD validation (convergence studies) and is absent from
  the YAML schema.
- **Missing fields for parameter file modification.** The existing workflow
  copies a template parameter file and modifies `SimPar@MaxMeshLevel` per
  level (see `moveAndSetLevel()`). The YAML schema has no equivalent.
- **No specification for test dependencies.** Some tests may require specific
  build configurations (e.g., `USE_PE=ON` + `USE_PE_SERIAL_MODE=ON` requires
  a two-step CMake configuration as documented in CLAUDE.md). The YAML
  `build.cmake.options` list does not capture this ordering requirement.

**Recommendation:** Define a formal schema (JSON Schema or a Python dataclass
model) with required/optional field annotations. Include multi-level test
support and parameter file templating as first-class features.

### Issue 3: Build stage design does not match actual build complexity

**Severity: Medium**

The YAML example shows a simple `cmake` + `targets` structure, but the actual
build process involves:

1. **Submodule initialization** (`git submodule update --init --recursive`) —
   not mentioned.
2. **Two-step CMake for PE serial mode** — the `USE_PE_SERIAL_MODE` option
   depends on `USE_PE=ON` being set first due to `cmake_dependent_option`.
   A flat list of options will fail for this case.
3. **Module loading** — the SLURM section mentions `module load` but the
   build stage does not account for compiler/MPI environment setup.
4. **Build caching** — no mention of whether builds are cached between test
   runs or rebuilt from scratch each time. For a CI-like system, this has
   major implications for turnaround time.

**Recommendation:** Add a `setup` stage before `build` for environment
preparation (modules, submodules). Support ordered CMake configuration steps.
Address build caching strategy explicitly.

### Issue 4: Metric extraction is too loosely defined

**Severity: Medium**

The existing `q2p1_ctest_start.py` has a precise extraction mechanism:
- Search for lines matching a keyword (e.g., `"BenchForce:"`)
- Split the line and extract values at specific column indices
- This handles the structured output format of FeatFloWer applications

The proposed YAML uses a freeform `command` field for extraction, which:
- Loses the structured keyword + index approach
- Makes comparison harder (how are `drag` and `lift` fields mapped to grep
  output?)
- Introduces shell injection risk if commands are not sanitized

**Recommendation:** Keep the keyword + column-index extraction model as a
built-in parser type. Reserve the freeform `command` approach as a fallback
for non-standard outputs. Example:

```yaml
metrics:
  drag_lift:
    parser: keyword_columns
    file: _data/prot.txt
    keyword: "BenchForce:"
    columns:
      drag: 2
      lift: 3
    last_occurrence: true
    compare:
      type: tolerance
      tolerance:
        drag: 1.0e-4
        lift: 1.0e-4
```

### Issue 5: No error handling or retry strategy

**Severity: Medium**

The document mentions that exit status is "recorded and propagated" but does
not address:

- **SLURM job failures**: timeout, OOM, node failure. The existing
  `submitAndObserveSync()` only checks for `COMPLETED`, `FAILED`, or
  `CANCELLED` but does not distinguish failure modes.
- **Build failures**: how are they captured and reported?
- **Partial failures**: what if the build succeeds but the run fails? What if
  metrics extraction fails on an otherwise successful run?
- **Retry policy**: should transient failures (e.g., SLURM scheduler
  congestion) trigger automatic retries?

**Recommendation:** Define a status model with at least: `SUCCESS`,
`BUILD_FAILURE`, `RUN_FAILURE`, `METRICS_FAILURE`, `TIMEOUT`, `INFRA_ERROR`.
Store per-stage exit codes in `metadata.json`.

### Issue 6: No baseline management strategy

**Severity: Medium**

The "Comparison & Alerting" section mentions comparing against "a baseline
with tolerance bands" but does not specify:

- How baselines are established (first successful run? manually curated?)
- Where baselines are stored (in-repo under `testcases/`? in the results
  store?)
- How baselines are updated (manually? automated promotion of passing runs?)
- Whether baselines are per-commit, per-branch, or per-configuration

**Recommendation:** Store baselines as versioned YAML files alongside test
definitions (e.g., `testcases/q2p1_fc_ext_cylinder.baseline.yaml`). Include a
CLI command to promote a successful run's metrics to become the new baseline.

### Issue 7: Missing concurrency and resource management

**Severity: Low**

The document does not address:

- Running multiple tests concurrently via SLURM (job arrays or independent
  submissions)
- Resource contention on shared clusters
- Prioritization of tests (fast smoke tests vs. long-running convergence
  studies)

**Recommendation:** Support a `priority` field in test definitions and
implement a simple scheduling strategy (e.g., submit all tests, let SLURM
handle scheduling; or define test suites with dependency ordering).

### Issue 8: Dashboard API is premature

**Severity: Low**

The architecture diagram includes a "Dashboard API (read-only, JSON)" layer,
and the tech stack recommends FastAPI. For the initial backend-first phase,
this is unnecessary complexity. The `metrics.json` and `metadata.json` files
are already machine-readable.

**Recommendation:** Defer the API layer entirely until there is a concrete
frontend consumer. Document the JSON schema so a future API can be built
against it, but do not implement it in phase 1.

---

## Additional Recommendations

### Integrate with existing CTest/CDash infrastructure

The project already has CTest configured with a CDash server at
`arthurdayne.mathematik.tu-dortmund.de/CDash/public`. Rather than building a
completely separate results store, consider whether CTest can serve as the
execution driver with the new YAML definitions feeding into it. This would
preserve the existing `add_test()` integration and CDash reporting while
improving test definition ergonomics.

### Add a `testcases/` migration path

The existing `testcases/` directory contains 6 test case directories (fac,
fac_nnewt, fac_visco, fac_visco_elastic, fallingparticle, gendie) with
parameter files and initial data. The design should clarify the relationship
between these existing test cases and the new YAML-defined tests.

### Define the CLI interface early

The document mentions `typer` or `click` but does not sketch the CLI
interface. Defining the command structure early helps validate the design:

```
featflower-test run <test-id> [--slurm | --local]
featflower-test list
featflower-test status <run-id>
featflower-test compare <run-id> --baseline <baseline-id>
featflower-test report <run-id> --format [json|markdown]
```

### Consider environment reproducibility

The document mentions RHEL 9.7 and module loading but does not address
reproducibility across environments. Consider capturing the full environment
(module list, compiler versions, MPI version) in `metadata.json` to enable
meaningful cross-environment comparisons.

---

## Summary of Recommendations

| # | Issue | Severity | Action |
|---|-------|----------|--------|
| 1 | No existing infrastructure acknowledgment | High | Add migration strategy section |
| 2 | Underspecified YAML schema | Medium | Define formal schema with multi-level support |
| 3 | Build stage oversimplified | Medium | Add setup stage, ordered CMake, build caching |
| 4 | Loose metric extraction | Medium | Keep keyword+column parser as primary method |
| 5 | No error handling strategy | Medium | Define status model with per-stage exit codes |
| 6 | No baseline management | Medium | Version baselines alongside test definitions |
| 7 | Missing concurrency model | Low | Add priority field and scheduling strategy |
| 8 | Premature API layer | Low | Defer to phase 2; document JSON schema only |

---

## Conclusion

The design document establishes a reasonable architectural direction for
FeatFloWer's automated test system. The key gap is the disconnect between the
proposal and the existing, functional test infrastructure. Addressing this gap
— through an explicit migration plan and by incorporating proven patterns from
the current `q2p1_ctest_start.py` (keyword-based extraction, multi-level
iteration, SLURM job management) — would significantly strengthen the design
and reduce implementation risk.

The recommended next step is to revise the design document to incorporate the
migration strategy and formalize the YAML schema before beginning
implementation.
