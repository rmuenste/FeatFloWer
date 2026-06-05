# Provenance Dump Status

## Scope

This note summarizes the current status of the provenance-aware dump writer and
reader introduced for deterministic frozen-field restart experiments.

It focuses on:

- the new `_dump_prov/<idx>/` format
- how to enable it in applications
- what is currently working
- what is still provisional

## Current Format

The provenance dump path writes text-based files into:

- `_dump_prov/<idx>/manifest.txt`
- `_dump_prov/<idx>/q2_ownership.csv`
- `_dump_prov/<idx>/q2_ownership_audit.csv`
- `_dump_prov/<idx>/velocity.csv`
- `_dump_prov/<idx>/coordinates.csv`
- `_dump_prov/<idx>/pressure.csv`
- `_dump_prov/<idx>/time.txt`
- optional Q2 scalar field CSVs such as `MaterialDistribution.csv`

The current implementation is in:

- [source/postprocessing/solution_io_provenance.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io_provenance.f90)

## Runtime Selection

Runtime selection is controlled by:

- `SimPar@UseProvDump = Yes|No`

The shared flag is held in:

- [source/src_util/prov_dump_config.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/src_util/prov_dump_config.f90)

The parameter is parsed centrally in:

- [source/src_util/param_parser.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/src_util/param_parser.f90)

`q2p1_el_frozen_trace` also has a local parser path, so its local
`myGDATNEW(...)` was extended to recognize `UseProvDump` as well:

- [applications/q2p1_el_frozen_trace/app_init.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/applications/q2p1_el_frozen_trace/app_init.f90)

Applications currently wired to the provenance path:

- `q2p1_fc_ext`
- `q2p1_fac3d`
- `q2p1_el_frozen_trace`

## Restart Semantics

For the provenance path, the current frozen-field application behavior is:

- `StartingProc = 1` uses `SolFromFileProv(...)`
- `StartingProc = 2` uses `SolFromFileProv(...)` on the lower-level structures
  and then prolongates
- `StartingProc = 3` uses `SolFromFileRepartProv(...)`

Important note:

- `read_sol_from_file_prov(...)`
- `read_sol_from_file_repart_prov(...)`

currently dispatch to the same common reader routine.

So at the moment:

- `StartingProc = 1` and `StartingProc = 3` are conceptually different
- but operationally they are effectively aliases in the provenance reader path

This is acceptable for the current prototype because the new format is intended
to be partition-independent by construction.

## What Is Working

### Writer

The provenance writer now completes successfully for `q2p1_fc_ext`.

Important fixes that were required:

- rank 0 must not assume its own local `myDump` shape
- root-side gather logic must adopt Q2/P1 slot counts from worker ranks
- provenance metadata must use the same dump-vertex identity space that the
  legacy dump machinery uses

The writer now produces complete `_dump_prov/<idx>/` directories.

### Validator

The validator script:

- [tools/postprocessing_scripts/provenance_dump_validate.py](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/tools/postprocessing_scripts/provenance_dump_validate.py)

successfully validates the generated dump structure and ownership tables.

This is a structural validation pass, not a full physics validation.

### Reader

`q2p1_el_frozen_trace` now actually enters the provenance reader path when:

- `SimPar@UseProvDump = Yes`

and the local app parser recognizes it.

### Frozen-Field Result Quality

The current observed result is:

- frozen-field serial with provenance looks plausible
- PE-parallel with provenance also works
- PE-parallel still differs slightly from PE-serial
- but the result is qualitatively similar, unlike the earlier legacy
  repartition path failures

This is the strongest indication so far that the provenance path improves
restart consistency materially.

## Known Caveats

### 1. Pressure Manifest Bug in Older Dumps

An early version of the provenance writer emitted:

- `p1_slots_per_coarse = 1`

in `manifest.txt`, while `pressure.csv` actually contained one row per slot of
the coarse dump structure.

This caused reader-side buffer overruns when loading older provenance dumps.

Current status:

- the reader now scans `pressure.csv` and sizes itself from the actual file
  content
- the manifest writer has been corrected for newly written dumps

So:

- old provenance dumps may print a pressure-slot mismatch warning
- but they are still loadable with the current reader

### 2. Ownership Meaning

The current Q2 ownership grouping uses the same `myDump%Vertices` identity that
the legacy dump writer uses for coarse-row field collection.

This has two consequences:

- large duplicate counts in `q2_ownership.csv` are not automatically invalid
- the ownership key is an output-grid dump-point identity, not a separate new
  FE-global dof numbering

This is internally consistent with the legacy dump chain, but the format should
still be treated as experimental until more restart comparisons are collected.

### 3. Diagnostics Output

`q2p1_el_frozen_trace` originally wrote too many per-step CSV diagnostics by
default.

Current status:

- `ELWriteDiagnostics` now defaults to off
- sampling debug files remain available behind:
  - `SimPar@ELWriteSamplingDebug = Yes`

So the debug writers are preserved, but normal runs are no longer flooded by
rank- or step-based CSV output.

## Recommended Usage

For a provenance-based frozen-field restart:

```text
SimPar@StartFile = 1
SimPar@StartingProc = 1
SimPar@UseProvDump = Yes
```

or, conceptually for repartition:

```text
SimPar@StartFile = 1
SimPar@StartingProc = 3
SimPar@UseProvDump = Yes
```

The current reader path treats both as the same provenance load operation.

## Recommended Next Checks

1. Generate short serial and PE-parallel frozen-field runs from the same
   `_dump_prov/1`.
2. Compare the first particle pass again.
3. Re-run the provenance validator on newly written dumps after any further
   writer changes.
4. If the format is adopted more broadly, add a dedicated regression test for:
   - write provenance dump
   - read provenance dump on same partitioning
   - read provenance dump on different partitioning
   - compare resulting frozen-field samples

## Current Assessment

The provenance dump path is no longer just a prototype that compiles.

It is now:

- writable
- structurally valid
- loadable in `q2p1_el_frozen_trace`
- and it produces the most plausible frozen-field restart behavior observed so
  far in this investigation

That is sufficient to preserve as a checkpointed experimental path, but not yet
enough to declare it fully validated for all applications and restart cases.
