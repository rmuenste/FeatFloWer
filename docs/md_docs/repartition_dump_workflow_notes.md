# Repartition Dump Workflow Notes

## Scope

This note summarizes the current understanding of how FeatFloWer dump files are
generated and consumed in the repartition workflow that is relevant for
`q2p1_el_frozen_trace`.

It focuses on:

- partitioned dump files in `_dump/processor_*/<idx>/`
- merged repartition dumps in `_dump/<idx>/`
- the function chains involved in writing and reading them
- the currently open ambiguity around reconstructing fine Q2 fields from a
  repartition dump

## File Layout

Two dump layouts are relevant:

1. Partitioned dump layout

   - `_dump/processor_<rank>/<idx>/velocity.dmp`
   - `_dump/processor_<rank>/<idx>/pressure.dmp`
   - `_dump/processor_<rank>/<idx>/coordinates.dmp`
   - `_dump/processor_<rank>/<idx>/MaterialDistribution.dmp`
   - `_dump/processor_<rank>/<idx>/time.dmp`

2. Merged repartition dump layout

   - `_dump/<idx>/velocity.dmp`
   - `_dump/<idx>/pressure.dmp`
   - `_dump/<idx>/coordinates.dmp`
   - `_dump/<idx>/MaterialDistribution.dmp`
   - `_dump/<idx>/time.dmp`

The merged layout is produced from the partitioned layout by
`tools/combine_fields.py`.

## Writer Chains

### 1. Legacy compact dump path

This is the legacy coarse-row dump path built around `WriteSol`.

Call chain:

- `SolToFile_Compact`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:220)
- `WriteSol_Velo`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:568)
- `WriteSol_Pres`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:615)
- `WriteSol_Coor`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:651)
- `WriteSol`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:682)

Purpose:

- write compact field dumps for one output index
- encode field data per coarse element row
- store Q2/Q1-related field values in `.dmp` files

Important clarification:

- `WriteSol` itself does not write per-rank `processor_*` files.
- In [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:705),
  nonzero ranks assemble local coarse-element rows and send them to rank 0.
- Rank 0 gathers those rows through `coarse%pELEMLINK` and writes a single
  merged file `_dump/<idx>_<field>.dmp`.

So `WriteSol` is still useful to study because it shows the coarse-row encoding,
but it is not by itself the writer of the `orig_dump/processor_*/<idx>/*.dmp`
layout.

### 2. Newer merged dump path

This is the path currently used by the application-side merged dump writers.

Call chain:

- `SolToFile`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:1)
- `write_vel_sol`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:494)
- `write_pres_sol`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:303)
- `write_q2_sol`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:1014)
- `write_time_sol`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:1323)

Purpose:

- write merged `.dmp` field files directly
- rely on `coarse%myELEMLINK`, `myDump%Vertices`, and `myDump%Elements`
- serialize field data according to the dump structures built by
  `CreateDumpStructures`

## Reader Chains

### 1. Repartitioned merged dump read path

This is the path used when loading a merged repartition dump into a possibly
different partitioning.

Call chain:

- `SolFromFileRepart`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:101)
- `read_vel_sol_single`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:623)
- `read_pres_sol_single`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:749)
- `read_q2_sol_single`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:885)
- `read_time_sol_single`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:1374)

Purpose:

- read one merged dump file per field from `_dump/<idx>/`
- reconstruct the local FE field arrays on the current partitioning

### 2. Non-repartition merged dump read path

Call chain:

- `SolFromFile`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:298)
- `read_vel_sol`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:575)
- `read_pres_sol`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:431)
- `read_q2_sol`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:1123)
- `read_time_sol`
  - [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90:1347)

Purpose:

- read merged field files with the usual merged dump path
- typically used when no repartition reconstruction is needed

## Dump Topology Structures

The common topology map for merged dump writing and reading is built in:

- `CreateDumpStructures`
  - [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:1738)

Purpose:

- build `myDump%Elements`
- build `myDump%Vertices`
- define how one coarse element row is expanded into fine subelements and Q2
  point locations

Inputs:

- local FE mesh topology
- `_data/MG.dat`

Key mechanism:

- `myDump%Elements(IEL,:)` stores all fine subelements for coarse element `IEL`
- `myDump%Vertices(IEL,IVT)` uses `MG.dat` to select a subelement and local
  vertex/point index

More concretely:

- `nLengthE = 8**(NLMAX-1)` is the number of fine P1 subelements inside one
  coarse element.
- `nLengthV = (2**(NLMAX-1)+1)**3` is the number of Q2 sample points inside one
  coarse element.
- `MG.dat` provides, for each coarse-row point `IVT`, a pair
  `(ke(IVT), jv(IVT))`.
- `ke(IVT)` selects which fine subelement inside the coarse element is used.
- `jv(IVT)` selects which local point of that fine subelement is used.
- `CreateDumpStructures` then stores
  `myDump%Vertices(IEL,IVT) = mg_mesh%level(NLMAX)%kvert(jv(IVT), kel)`.

This means the coarse-row Q2 representation is purely:

- an ordered list of Q2 point values per coarse element

It does not store:

- which coarse element owns a shared fine Q2 point
- which neighboring coarse row should be authoritative if the same reconstructed
  fine Q2 dof can be reached from more than one coarse row

That missing provenance is the current reason to suspect that repartitioned
reconstruction is underdetermined.

## What `WriteSol` Actually Writes

For Q2-like fields, `WriteSol` uses `CollectVertField`.

Mechanism:

- for each local coarse element `IEL`
- for each coarse-row point `IVT = 1..nLengthV`
- it looks up `JVT = myDump%Vertices(IEL,IVT)`
- and writes `xField(JVT)` into the coarse-row buffer slot `Field(IVT,IEL)`

So a velocity/coodinate/material row is written as an ordered vector of length:

- `DofsInElement = (2**(NLMAX+iiLev)+1)**3`

For `NLMAX+iiLev = 2`, that is:

- `DofsInElement = 5^3 = 125`

The row therefore stores:

- values at the 125 structured Q2 sample points of one coarse element

It does not store:

- the source fine-element id `kel`
- the `MG.dat` pair `(ke,jv)` used to obtain each slot
- any ownership flag for sample points lying on coarse-element interfaces
- any provenance saying which coarse row should win during repartitioned
  reconstruction

For P1/pressure-like fields, `WriteSol` uses `CollectElemField`, which is
analogous but based on `myDump%Elements(IEL,:)`.

This structure is central to both:

- writing coarse-element rows
- reconstructing fine local fields from merged repartition dumps

## `combine_fields.py`

Script:

- [tools/combine_fields.py](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/tools/combine_fields.py)

Purpose:

- read partitioned dump files from `_dump/processor_*/<idx>/`
- merge them into `_dump/<idx>/`

Current observed behavior:

- reads every `(coarse element id, component payload)` tuple from all
  `processor_*` folders
- sorts the collected rows by coarse element id
- writes them into the merged dump folder

Important finding from `orig_dump/processor_*/1/velocity.dmp`:

- `31` processor folders were checked
- `612` unique coarse element ids were found
- `0` duplicate coarse element ids were found across processors
- recomputing the merged `orig_dump/1/velocity.dmp` with `combine_fields.py`
  reproduced the existing file byte-for-byte

So for `velocity.dmp` in `orig_dump/1`, `combine_fields.py` does not appear to
be corrupting the merge through duplicate coarse-element rows.

## Current Reconstruction Findings

### Established

- The serial and PE-parallel frozen-field discrepancy appears already at step 1.
- Particle positions and particle velocities match exactly at step 1.
- The discrepancy is in the sampled carrier field, not in PE stepping.
- The discrepancy is already present in the local sampled Q2 field values before
  interpolation.
- The earlier blind overwrite behavior in the repartition reader was a real bug.

### Patched already

In:

- [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90)

we tried two reader-side reconstruction policies for overlapping Q2 dofs:

1. first write wins
2. smallest `elemmap(i_local)` wins

Both improved the situation significantly compared to the original overwrite
behavior, but neither reproduced the serial reference exactly.

### Current interpretation

The remaining discrepancy suggests that the repartition dump representation is
still underdetermined for unique fine-level Q2 reconstruction.

That means:

- the merged repartition dump is likely a faithful union of coarse-element rows
- but the dump format does not encode enough provenance or ownership information
  for the reader to know which coarse-row source is authoritative for a shared
  reconstructed fine Q2 dof

The direct `WriteSol` audit supports this:

- the 125-slot coarse-element row stores values only
- the mapping from slot to fine Q2 dof is implicit in local topology and
  `MG.dat`
- no explicit ownership/provenance is written for shared interface points

This is currently the leading hypothesis.

## Open Questions

1. Which exact writer routine produces the partitioned `processor_*/<idx>/*.dmp`
   files in the workflow that created `orig_dump`?
2. Does legacy `WriteSol` already imply a hidden ownership convention that the
   merged repartition reader does not reproduce?
3. Is additional ownership/provenance metadata needed in the repartition dump
   format to reconstruct Q2 fields uniquely?
4. Does the non-repartition dump path avoid the issue entirely, as expected?

## Recommended Next Steps

1. Audit legacy `WriteSol` in [source/OutputProfiles.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/OutputProfiles.f90:682)
   as the likely partitioned `.dmp` writer.
2. Compare repartition versus non-repartition dump loading for the same frozen
   field.
3. If the issue is repartition-specific, extend the dump format or writer logic
   so the authoritative source for shared reconstructed Q2 dofs is explicit.
