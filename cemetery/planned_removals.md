# Planned Removals Registry

This directory is a holding area for code, data structures, modules, source files,
and configuration fragments that are suspected to be obsolete or superseded.

The intent is to keep a visible, version-controlled backlog of ballast in the
active codebase without removing it prematurely.

## Rules

1. Entries in this file are candidates only.
2. No item listed here should be removed autonomously.
3. Every removal requires an explicit author review and confirmation.
4. When a suspiciously unused or superseded item is found, add it here together
   with the reason, current status, and blocking dependencies.
5. Prefer listing exact subroutines, modules, structures, files, or parameter
   fragments instead of vague areas.

## Entry Format

Use this structure for new items:

```text
### Item Name
- Type: subroutine | module | structure | source file | parameter block | other
- Location: path:line or path
- Status: candidate | blocked | ready-for-review | removed
- Reason:
- Known dependencies:
- Notes:
```

## Current Candidates

### FilterColdElements
- Type: subroutine
- Location: `source/src_quadLS/QuadSc_material_properties.f90`
- Status: candidate
- Reason: currently appears unused in the active code path; the known call site in
  `QuadSc_main.f90` is commented out.
- Known dependencies: none identified so far.
- Notes: verify there are no out-of-tree users or compile-time variants that still
  enable it before removal.

### HEATALPHA entry points for multi-material alpha transport
- Type: subroutine group
- Location: `source/src_LinSc/GenLinSc_transport_extensions.f90`,
  `source/src_LinSc/GenLinSc_def_extension.f90`,
  `source/src_LinSc/GenLinSc_user.f90`
- Status: blocked
- Reason: `q1_scalar_multimat` and `q2p1_sse` have been moved toward alpha-only
  workflows, so the old `temp + alpha...` path is becoming legacy ballast.
- Known dependencies: likely still used by `q1_scalar_melt`, `q2p1_sse_melt`,
  and `q2p1_sse_mesh`; these users must be audited first.
- Notes: candidate routines include `Init_GenLinSc_HEATALPHA_Q1`,
  `Matdef_HEATALPHA_GenLinSc_Q1`, `Knpr_GenLinSc_HEATALPHA_Q1`,
  `Boundary_GenLinSc_HEATALPHA_Q1_Val`, `CheckAlphaConvergence`,
  `Correct_GenLinSc_Q1_ALPHA`, and `Add_DissipativeEnergy_HEATALPHA_Q1`.

### Temperature-first field layout assumptions in legacy multimat transport
- Type: structure / indexing convention
- Location: legacy `GenLinScalar` users in `source/src_LinSc` and `source/src_quadLS`
- Status: blocked
- Reason: several historical paths assume `GenLinScalar%Fld(1)` is temperature and
  alpha fields start at index 2. This convention is already being phased out in
  alpha-only workflows and should eventually disappear.
- Known dependencies: any remaining apps that still rely on the old HEATALPHA
  transport layout.
- Notes: keep this as a tracking item until all remaining legacy users are either
  migrated or explicitly retained.

### Stale temperature-oriented field naming in q1_scalar_multimat templates
- Type: parameter block
- Location: `applications/q1_scalar_multimat/_data/q2p1_param.dat`
- Status: candidate
- Reason: the template still reflects the old temperature-led field naming even
  though the solver path is being moved to alpha-only semantics.
- Known dependencies: restart and postprocessing conventions should be checked
  before changing or removing the old naming fragment.
- Notes: this is not a code removal item by itself, but it is part of the same
  legacy track and should stay visible.

### q2p1_sse_mesh
- Type: application / source subtree
- Location: `applications/q2p1_sse_mesh`
- Status: candidate
- Reason: suspected legacy development track that should be reviewed for removal
  from the active codebase.
- Known dependencies: application-specific build wiring, startup path, and any
  remaining users or benchmark cases need to be checked before deletion.
- Notes: keep this item visible until its runtime purpose and maintenance status
  are explicitly reviewed by the author.

### run_q2p1_fc_ext helper scripts
- Type: source file (shell script), 5 byte-identical copies
- Location: `applications/q2p1_creep/run_q2p1_fc_ext`,
  `applications/q2p1_dns_drag/run_q2p1_fc_ext`,
  `applications/q2p1_fc2/run_q2p1_fc_ext`,
  `applications/q2p1_fc_ext/run_q2p1_fc_ext`,
  `applications/q2p1_xParticles/run_q2p1_fc_ext`
- Status: removed
- Reason: superseded copy-paste template cruft; all five copies share one md5
  (`b599620e50ff88b7fa61428b46e5979d`) and none can run as written:
  1. partitioning invokes `python /home/rafa/bin/PyPartitioner.py` — a hardcoded
     path into a foreign user's home directory, via the Python 2 `python` binary;
  2. the run step is `mpirun -np 5 ./q2p1_fc_ext`, but `add_executable(q2p1_fc_ext ...)`
     is commented out in four of the five hosting apps, which build `q2p1_creep`,
     `q2p1_dns_drag`, `q2p1_fc2` and `q2p1_xParticles` instead — so that binary
     exists in only one of the five directories;
  3. results extraction greps `_data/prot.txt` for `Force acting`, a string no
     source file emits. grep matches nothing, so the awk pipeline appends an
     empty line to `results.txt` rather than failing — a silent-wrong-answer
     path. (This same dead grep was propagated into guide 01 and has been fixed
     there; the scripts are its last remaining home.)
- Known dependencies: none found. No reference from any `CMakeLists.txt`,
  `*.cmake`, CI config, test harness, or documentation, and CMake never copies
  the scripts into a build tree — verified against three configured builds.
  The live replacement is `q2p1_fc_ext_start.py` together with
  `tools/dashboard/q2p1_ctest_start.py`: both are staged via `file(COPY ...)`
  and driven by `add_test(q2p1-fac-newt python ./q2p1_ctest_start.py)`.
- Notes: each copy's only commit is the one that created its application
  (2017-07-24 for `q2p1_fc_ext`, through 2026-03-10 for `q2p1_dns_drag`); none
  has ever been edited afterwards. Marked ready-for-review rather than candidate
  because the supersession is unambiguous and no dependency was found, but per
  rule 3 removal still needs explicit author confirmation.

### test-config.xml placeholder templates
- Type: configuration fragment, 9 copies
- Location: `applications/*/test-config.xml` (9 copies)
- Status: removed
- Reason: adjacent finding from the `run_q2p1_fc_ext` review. All nine copies
  still contain the unfilled template values (`meshFolderPath="/path/to/folder"`,
  `meshProjectFile="/path/to/file"`, `testDataFile="/path/to/file"`), so none
  describes a real test configuration.
- Known dependencies: none found — no reference from any `CMakeLists.txt`,
  `*.cmake`, Python tooling, or documentation. The CTest path in use reads
  `tests/*.json` instead (e.g. `applications/q2p1_fc_ext/tests/test-fac2d.json`).
- Notes: listed separately from the shell scripts because the two may have been
  part of the same abandoned harness; confirm whether an XML-driven test runner
  was ever completed before removing.
