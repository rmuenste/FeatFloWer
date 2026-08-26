# MPI-PRF checkpoint implementation plan

## Goal

Provide MPI-IO `.prf` checkpoints as a user-selectable dump and restart
backend, configured through `SimPar@DumpFormat` in `q2p1_param.dat`. The
implementation must expose a shared application-facing component instead of
adding application-specific calls to the low-level MPI-IO routines.

The first adopters are `q2p1_sse`, `q2p1_sse_temp`, and `q2p1_fc_ext`.
`q2p1_ATC` is a deferred adopter that should be able to register through the
same interface.

## Configuration policy

1. Add `SimPar@DumpFormat` with the values:

   - `MPI_PRF`, with `MPI` accepted as an alias;
   - `DMP`;
   - `LST`;
   - `PROVENANCE`.

2. Let applications register their supported formats before the common
   parameter parser runs:

   - `q2p1_sse`: `MPI_PRF`, `LST`;
   - `q2p1_sse_temp`: `MPI_PRF`, `LST`;
   - `q2p1_fc_ext`: `MPI_PRF`, `DMP`, `PROVENANCE`.

3. Use `MPI_PRF` when neither `DumpFormat` nor the compatibility parameter is
   present. Maintained `q2p1_fc_ext` input files explicitly select `DMP` to
   preserve their intended behavior, while maintained SSE inputs explicitly
   select `MPI_PRF`.

4. Preserve `SimPar@UseProvDump` as a compatibility input:

   - `Yes` maps to `PROVENANCE`;
   - `No` maps to `DMP`;
   - equivalent old and new selections are accepted;
   - conflicting selections are fatal.

5. Report the resolved format and whether it came from the default, the new
   parameter, or the legacy parameter in the protocol output. Reject unknown
   or application-unsupported formats during initialization.

## Shared component design

Implement the feature in four layers:

1. `checkpoint_config.f90`

   Owns format constants, parsing, aliases, compatibility mapping,
   application registration, validation, and conversion of `StartFile` to a
   numeric MPI-PRF slot.

2. `checkpoint_types.f90`

   Defines the checkpoint context, write/restart modes, layout identifiers,
   and typed field selections. Application code selects semantic fields such
   as velocity, pressure, coordinates, temperature, or generic scalars rather
   than constructing character-code lists itself.

3. `checkpoint_service.f90`

   Provides the application-facing read/write facade and translates typed
   requests to the MPI-PRF backend. It handles replacement writes,
   field-merge writes, same-level reads, lower-level prolongation, and
   repartition reads.

4. `checkpoint_mpi_prf.f90` and the existing MPI dump routines

   Own the on-disk metadata contract, transaction lifecycle, payload
   validation, MPI-IO operations, clock restoration, and compatibility with
   markerless legacy MPI-PRF checkpoints.

The low-level field implementation must check optional arrays before use.
Unavailable optional fields are skipped with a field-specific warning. Core
application fields remain selected through the typed facade.

## On-disk format and transaction

Use `_dump/<slot>/` and retain the established MPI-PRF payload filenames. Add
`checkpoint_meta.prf` as the versioned commit marker. Do not claim `time.prf`,
because legacy SmartDump already uses that name.

The metadata records:

- format version and generation;
- simulation time and completed time-step number;
- writer rank count;
- native integer and real sizes and byte order;
- coarse-element count and source mesh level;
- field name, layout, component count, chunk count, global entry count, and
  payload bytes per component.

Readers use the recorded chunk count instead of recomputing it from
`DataSizeThresholdMPI`. Before MPI-IO starts, they validate the metadata
version, native representation, manifest, aggregate payload sizes, coarse
mesh, and expected source level.

Replacement writes follow this order:

1. create `checkpoint_incomplete.prf`;
2. remove the previous commit marker;
3. clear recognized MPI-PRF payloads from the reused slot;
4. write and truncate all selected key and component files;
5. write metadata to `checkpoint_meta.prf.tmp`;
6. atomically rename it to `checkpoint_meta.prf`;
7. remove the incomplete sentinel.

This makes interrupted generations invalid and removes orphaned files if a
later generation contains fewer fields or chunks. SSE-temperature uses a
field-merge transaction so that updating temperature retains the momentum
fields and their manifest entries. Field cleanup must match the complete
field-name prefix so similarly named fields cannot remove one another.

Markerless MPI-PRF checkpoints remain readable with a warning and without
clock restoration. A legacy SmartDump `time.prf` beginning with `timens` is
recognized separately. Incomplete, malformed, incompatible, missing, or
truncated versioned checkpoints fail before field loading.

## Simulation clock and mesh-level handling

Store `timens` and the completed `istep_ns` value in the metadata at each
committed replacement write. Same-level and repartition restarts restore both
values. Lower-level restarts restore them before prolongating the loaded
solution.

Finalization must pass `istep_ns - 1` as the completed step without mutating
the global counter. The MPI checkpoint refactor must likewise use local level
variables instead of temporarily incrementing and restoring `NLMAX`, so an
error path cannot leave global mesh state corrupted.

For lower-level `q2p1_fc_ext` restarts, prolongate every allocated generic
scalar after loading it, matching the existing SSE treatment.

## Application integration

### q2p1_fc_ext

- Register `MPI_PRF`, `DMP`, and `PROVENANCE` before common initialization.
- Dispatch intermediate and final writes through the selected backend.
- Support MPI-PRF `StartingProc` modes 1, 2, and 3.
- Restore simulation time and completed step for MPI-PRF restarts.
- Preserve provenance and classic-DMP restart paths.

### q2p1_sse

- Register `MPI_PRF` and `LST` before parsing.
- Replace the integer `myTransientSolution%DumpFormat` decisions at the
  checkpoint integration points with the shared policy.
- Use the typed SSE field set for same-level and repartition MPI-PRF restarts.
- Use the existing lower-level prolongation sequence after the shared reader.
- Preserve the existing non-MPI restart behavior.

### q2p1_sse_temp

- Register `MPI_PRF` and `LST` before parsing.
- Load each transient angle through one typed MPI-PRF field request.
- Select segment and generic scalar fields only when configured and allocated.
- Merge the temperature field into an existing MPI-PRF slot instead of
  replacing its other fields or clock metadata.

## OpenMPI and ROMIO setup

Move launcher-side MPI-IO environment handling into a shared Python helper.
Query `ompi_info` for the active OpenMPI I/O components and select the
installed ROMIO component when `OMPI_MCA_io` is not already set.

Supply these application defaults only when the user has not provided values:

```text
ROMIO_CB_BUFFER_SIZE=16777216
ROMIO_DS_WRITE=enable
```

Preserve explicit values. Reject an explicit `OMPI_MCA_io` value that is not
available in the active OpenMPI installation, and leave ROMIO variables alone
when a non-ROMIO component is explicitly selected.

## Verification plan

### Automated and build-time checks

- Test default, explicit, alias, legacy, and conflicting configuration paths.
- Round-trip simulation time and completed step exactly.
- Recognize markerless MPI-PRF and legacy SmartDump directories.
- Reject incomplete transactions and truncated payloads.
- Verify replacement cleanup removes orphaned chunks.
- Verify field cleanup preserves prefix-related field names.
- Unit-test ROMIO detection, defaults, overrides, non-ROMIO selection, and
  stale component rejection.
- Compile and link `q2p1_fc_ext`, `q2p1_sse`, and `q2p1_sse_temp`.

### Reduced-rank runtime correctness

Use a prepared small case to write with 4 MPI ranks and restart the same slot
with 8 ranks using `StartingProc = 3`. Compare restored fields bit-for-bit on
the same architecture and verify the restored time and step against
`checkpoint_meta.prf`. Also exercise same-rank same-level restart and
lower-level prolongation.

The former 72-to-288-rank case is deferred because of queue cost.

### ROMIO performance matrix

Screen these configurations at 8 ranks with at least three repetitions:

- untouched OpenMPI defaults;
- selected ROMIO with ROMIO's internal defaults;
- 4, 16, and 64 MiB collective buffers with data sieving enabled;
- 16 MiB with data sieving disabled.

Retain solver logs and record total runtime, checkpoint-only runtime, payload
bytes, and exit status. Rerun only the two fastest configurations at 16 ranks.
Use equivalent prepared states for every measurement before drawing a default
configuration conclusion.

## Deferred work

- Register `q2p1_ATC` as the next non-SSE adopter.
- Add cross-endian or cross-kind payload conversion; current payloads remain
  native MPI-IO data and incompatible metadata is rejected.
- Add classic `.dmp` parity to SSE if it becomes a requirement.
- Reconcile particle-side `REPART` vocabulary beyond accepting `MPI` as the
  shared MPI-PRF alias.
- Run large 72-to-288-rank scalability measurements when queue conditions and
  campaign needs justify them.
