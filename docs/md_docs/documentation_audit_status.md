# Documentation Audit Status

**Date:** 2026-06-07

## Scope

This audit removed obsolete implementation-update notes, corrected confirmed
factual and path errors, and refreshed the automated-test, QuadSc, and
provenance documentation. The Euler-Lagrange documents remain pending manual
review.

## Removed Documents

The following files were removed because they were historical implementation
records, duplicated maintained documentation, or described behavior that no
longer matches the repository:

- `automated_test_system_review.md`
- `guide_x_custom_application_agent_workflow.md`
- `hashgrid_initialization_fix.md`
- `hashgrid_test_implementation.md`
- `kvel_implementation_summary.md`
- `repartition_q2_load_overwrite_bug.md`
- `testing_guide.md`

The retained replacements are:

- Automated tests: `automated_test_system.md` and
  `featflower_test_usage_guide.md`
- FBM acceleration: `fbm_acceleration_usage.md`,
  `hashgrid_verification.md`, and `kvel_force_acceleration.md`
- Application creation: `guide_05_agent_new_application_workflow.md`
- Repartition analysis: `repartition_dump_workflow_notes.md` and
  `provenance_dump_format.md`

## Corrected Information

- Documented that `ENABLE_FBM_ACCELERATION` defaults to `OFF` and must be
  enabled explicitly.
- Reconciled Guide 05 with the repository's active CTest support.
- Listed both included `featflower-test` definitions.
- Corrected the Guide 01 baseline path in its documentation and YAML
  definition.
- Replaced the removed `tools/combine_fields.py` reference with the maintained
  `tools/featflower_combine` package and CLI.
- Replaced workstation-specific absolute links with repository-relative links.
- Removed a reference to the nonexistent `el-frozen/README.md`.
- Replaced references to nonexistent or incorrectly cased QuadSc/UMFPACK
  documentation and source files.
- Removed stale `Last Updated` fields that did not reflect later repository
  changes.

## Validation

- No references to the seven removed documents remain.
- No `/data/warehouse17` or `FF-ATC-NEW` paths remain in `docs/md_docs`.
- The corrected Guide 01 test definition validates through `featflower-test`.
- Markdown and YAML changes pass `git diff --check`.

## Follow-up Work Completed

- Recast `automated_test_system.md` as current architecture plus remaining work.
- Consolidated the QuadSc planning notes into
  `quadsc_refactoring_status.md`.
- Refreshed `quadsc_current_implementation.md`.
- Recast the provenance status note as `provenance_dump_format.md`.
- Added `README.md` as the documentation index and linked it from `AGENTS.md`.
- Moved `docs/CONTAINER.md` to `containerization.md` and updated it to match
  the current root `Dockerfile`.

The Euler-Lagrange documents were intentionally left for manual review.
