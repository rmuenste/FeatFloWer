# MPI-PRF checkpoints

The shared checkpoint policy lets participating applications select their dump
backend in `q2p1_param.dat`:

```text
SimPar@DumpFormat = MPI_PRF
```

Accepted values are `MPI_PRF`, `DMP`, `LST`, and `PROVENANCE`; `MPI` is an
alias for `MPI_PRF`, matching the particle-side spelling. Applications register
the formats they implement before the common parameter parser runs. An
unsupported value or an unknown token is fatal. If neither the new key nor the
legacy key is present, the selection is `MPI_PRF`. Maintained `q2p1_fc_ext`
decks explicitly select `DMP`, while maintained SSE decks explicitly select
`MPI_PRF`.

`SimPar@UseProvDump` remains a compatibility input. `Yes` maps to
`PROVENANCE`, and `No` maps to `DMP`. If both keys occur, equivalent values are
accepted and conflicting values are fatal. New decks should use only
`DumpFormat`.

## Application support

| Application | Supported formats |
| --- | --- |
| `q2p1_sse` | `MPI_PRF`, `LST` |
| `q2p1_sse_temp` | `MPI_PRF`, `LST` |
| `q2p1_fc_ext` | `MPI_PRF`, `DMP`, `PROVENANCE` |

Registration, parsing, aliases, and compatibility policy live in
`source/src_util/checkpoint_config.f90`. The backend-neutral field/context
types are in `source/src_util/checkpoint_types.f90`; MPI-PRF metadata and
transaction handling are in `source/src_util/checkpoint_mpi_prf.f90`, and
`source/src_util/checkpoint_service.f90` is the application-facing read/write
facade. This keeps application integration to registration plus typed field
selection and read/write dispatch. A future `q2p1_ATC` adopter should use the
same boundary rather than adding another application-local dump switch.

## On-disk contract

MPI-PRF uses `_dump/<slot>/` and preserves the established key and component
payload names. A successfully committed generation also contains
`checkpoint_meta.prf`. `time.prf` is deliberately not claimed because the
legacy SmartDump writer already owns that name.

The versioned metadata records:

- simulation time and completed step;
- generation, writer rank count, native integer/real sizes, and byte order;
- coarse-element count and source mesh level;
- every field's layout, component count, chunks per component, global entry
  count, and payload bytes per component.

Readers use the recorded chunk layout and validate representation, manifest,
and aggregate payload sizes before reading. Same-level and repartition
restarts restore the clock; lower-level MPI-PRF restarts restore the clock and
then prolongate. `StartFile` may be a numeric slot or a path whose final
component is numeric, such as `_dump/1`.

Writes are transactional at slot level: create an incomplete sentinel, remove
the old commit marker, clear recognized old MPI-PRF payloads for a replacement,
write/truncate all selected payloads, write metadata to a temporary file, then
rename the metadata into place and remove the sentinel. This ordering prevents
a partially overwritten slot from looking valid and removes orphan chunks when
a new generation is smaller. SSE-temperature uses field-merge mode so its
temperature update retains the momentum fields and manifest entries.

Marker-less MPI-PRF data remains readable with one warning and without clock
restoration. A marker-less legacy SmartDump directory whose `time.prf` starts
with `timens` is identified separately rather than treated as malformed MPI-PRF
metadata. A sentinel without a commit marker, malformed metadata, incompatible
native representation, or missing/truncated payload is fatal.

## OpenMPI/ROMIO performance

The E3D launchers share `tools/e3d_scripts/mpi_io_environment.py`. It queries
`ompi_info`, selects the installed versioned ROMIO component when
`OMPI_MCA_io` is unset, and supplies these application defaults:

```text
ROMIO_CB_BUFFER_SIZE=16777216
ROMIO_DS_WRITE=enable
```

Explicit environment values are preserved. An explicit `OMPI_MCA_io` that is
not provided by the active OpenMPI installation fails early, avoiding a stale
`romio321` setting after an MPI upgrade. Selecting a non-ROMIO component leaves
ROMIO-specific variables untouched.

The matrix in `testcases/checkpoint_mpi_prf/` compares untouched OpenMPI,
ROMIO's internal defaults, several collective-buffer sizes, and data-sieving
settings. Screen at 8 ranks and rerun the best two profiles at 16 ranks. The
correctness repartition test is 4 ranks writing and 8 ranks reading; large
72-to-288 runs are deferred until queue conditions make them economical.

## Current portability boundary

Payloads remain native-endian MPI-IO data. Metadata detects incompatible kind
sizes or byte order and fails rather than converting. Cross-architecture
conversion, SSE parity for classic `.dmp`, and particle `REPART` vocabulary are
outside this component's current scope.
