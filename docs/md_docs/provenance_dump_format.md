# Experimental Provenance Dump Format

The provenance dump format is an experimental, text-based restart format
designed to make Q2 field ownership explicit enough for deterministic
cross-partition frozen-field experiments.

It is implemented and usable, but it is not yet a fully validated replacement
for the established dump formats.

## Enabling the Format

Set:

```text
SimPar@UseProvDump = Yes
```

The shared option is parsed in `source/src_util/param_parser.f90` and stored by
`source/src_util/prov_dump_config.f90`.

Applications currently wired to this path include:

- `q2p1_fc_ext`
- `q2p1_fac3d`
- `q2p1_el_frozen_trace`

The implementation is in:

- [solution_io_provenance.f90](../../source/postprocessing/solution_io_provenance.f90)

## Directory Layout

One output index is written below `_dump_prov/<idx>/`:

```text
_dump_prov/<idx>/
├── manifest.txt
├── q2_ownership.csv
├── q2_ownership_audit.csv
├── velocity.csv
├── coordinates.csv
├── pressure.csv
├── time.txt
└── <optional-q2-field>.csv
```

Optional Q2 scalar files can include fields such as
`MaterialDistribution.csv`.

## File Roles

### `manifest.txt`

Records format dimensions and metadata needed by the reader, including coarse
element and Q2/P1 slot counts.

### `q2_ownership.csv`

Records the selected source for each dump-point identity used by Q2 fields.

### `q2_ownership_audit.csv`

Records ownership candidates and duplicate information for diagnostics. Large
duplicate counts are not automatically invalid because the identity space is
derived from the established dump-grid mapping.

### Field CSV Files

`velocity.csv`, `coordinates.csv`, and optional scalar CSV files store values
using the Q2 ownership map. `pressure.csv` stores pressure data using its
coarse-grid slot layout.

### `time.txt`

Stores the simulation time associated with the output index.

## Writer Behavior

The writer gathers metadata and field rows from worker ranks. Rank 0 adopts the
Q2 and P1 slot dimensions supplied by the workers rather than assuming that its
own local dump structures describe the distributed field.

The ownership identity currently uses the same `myDump%Vertices` space as the
legacy coarse-row dump machinery. It is an output-grid dump-point identity, not
a new globally unique finite-element DOF numbering.

## Reader Behavior

The public entry points are:

- `SolFromFileProv(...)`
- `SolFromFileRepartProv(...)`

Both currently dispatch to the same common provenance reader. Therefore,
same-partition and repartition modes are conceptually distinct at the
application layer but operationally share the same format reader.

In `q2p1_el_frozen_trace`:

- `StartingProc = 1` loads through `SolFromFileProv`.
- `StartingProc = 2` loads lower-level structures and then prolongates.
- `StartingProc = 3` loads through `SolFromFileRepartProv`.

## Structural Validation

The validator is:

- [provenance_dump_validate.py](../../tools/postprocessing_scripts/provenance_dump_validate.py)

It checks the manifest, ownership tables, field structure, and expected row
relationships. This is structural validation only; it does not prove that a
restart reproduces the same physical solution.

## Compatibility

Older experimental dumps may contain:

```text
p1_slots_per_coarse = 1
```

even when `pressure.csv` contains more pressure slots. The current reader scans
the pressure file, sizes its buffers from the file contents, and emits a
warning when the manifest disagrees. Newly generated dumps use the corrected
manifest value.

## Recommended Usage

For a same-partition frozen-field load:

```text
SimPar@StartFile = 1
SimPar@StartingProc = 1
SimPar@UseProvDump = Yes
```

For a cross-partition experiment:

```text
SimPar@StartFile = 1
SimPar@StartingProc = 3
SimPar@UseProvDump = Yes
```

Always run the validator after producing a new dump.

## Validation Limits

The format has demonstrated:

- successful writing from `q2p1_fc_ext`,
- structurally valid ownership and field files,
- successful loading in `q2p1_el_frozen_trace`, and
- materially better frozen-field restart behavior than the ambiguous legacy
  repartition reconstruction in the investigated cases.

It has not yet demonstrated:

- exact serial versus PE-parallel equivalence,
- broad coverage across applications and mesh layouts,
- same-partition and cross-partition round-trip regression tests, or
- long-term format stability.

Treat the format as experimental until those checks are automated.

## Required Next Checks

1. Add a write/read round-trip test on the same partitioning.
2. Add a write/read test with a different partitioning.
3. Compare reconstructed Q2 fields before running particle dynamics.
4. Compare first-step frozen-field samples between serial and parallel runs.
5. Decide whether the dump-point identity is sufficient as a permanent
   ownership key or whether a dedicated global FE identity is required.
