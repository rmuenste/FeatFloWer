# Repartition Q2 Load Overwrite Bug

## Summary

There is a clear bug in the repartitioned Q2 dump reader path used by FeatFloWer applications that load a velocity or other Q2 field from `_dump/<StartFile>` with a different partitioning.

The affected logic is in the repartitioned single-file readers:

- [source/postprocessing/solution_io.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/source/postprocessing/solution_io.f90)
  - `read_vel_sol_single`
  - `read_q2_sol_single`

The core issue is that the reader blindly scatters coarse-element buffer values into local Q2 dofs:

```fortran
jvt = edofs(i_local,ivt)
u(jvt) = buf(ivt)
```

and similarly for `v`, `w`, and generic Q2 fields.

If the same local Q2 dof appears in multiple local coarse rows, it is written multiple times. The last writer wins.

There is currently no ownership rule, no first-write rule, and no deduplication.

## Why This Is A Bug

For repartitioned Q2 loading, the same local Q2 dof can be referenced by multiple coarse-element rows in `myDump%Vertices`.

That is not inherently wrong by itself. The bug is that the reader treats every occurrence as an authoritative assignment and overwrites the same dof repeatedly.

This makes the final loaded field depend on:

- local coarse-row ordering
- partition layout
- overlap pattern between coarse rows

So the final loaded Q2 field is not uniquely defined by the dump file alone.

## Observed Symptom

In `q2p1_el_frozen_trace`, the repartition-loaded frozen velocity field showed:

- plausible values on some sampled Q2 dofs
- exact zeros or wrong values on neighboring sampled Q2 dofs
- a discrepancy already at the first particle pass between serial and PE-parallel runs

This was first visible through sampled carrier-velocity differences, and then localized with dedicated diagnostics.

## Debugging Result

The following points were verified:

1. Serial and parallel sample the same physical element.
2. The local reference coordinates `xi1,xi2,xi3` match.
3. The sampled hexahedron vertex coordinates match.
4. The repartition map `myDump%Vertices` contains the relevant dofs.
5. The loaded `QuadSc%ValU/V/W` values do not match the corresponding dump-buffer values for the identified map slots.
6. The same local Q2 dofs are written from multiple coarse rows during load.

This confirms that the corruption comes from repeated writes during repartitioned field loading.

## Concrete Evidence

Examples from the frozen-field investigation:

- A single local Q2 dof such as `135`, `184`, `492`, `5188`, `6592`, or `11920` can have multiple sources.
- The same dof may be written by 2 to 8 different coarse elements during one load.
- The final field value on that dof matches neither a unique source nor a stable ownership rule.

This explains why one run configuration can appear acceptable while another shows obvious corruption: the bug is order-dependent.

## Why Serial Can Still Look Fine

The serial-reference case can still appear correct even though it uses the same repartitioned reader because:

- the local overlap pattern differs
- the local coarse-row order differs
- therefore the final last-writer overwrite differs

So this is not evidence that the reader is correct. It only means the bug is partially masked in that configuration.

## Affected Scope

This is not specific to `q2p1_el_frozen_trace`.

Any application that uses the repartitioned Q2 dump reader path may be affected, especially when:

- loading from `_dump/<StartFile>` with changed partitioning
- reading velocity through `SolFromFileRepart`
- reading custom Q2 fields through `read_q2_sol_single`

Potentially affected applications include any app that uses:

- `SolFromFileRepart(...)`
- `read_vel_sol_single(...)`
- `read_q2_sol_single(...)`

## Likely Fix

The reader must not overwrite the same local Q2 dof multiple times without a canonical rule.

The simplest robust fix is:

1. Allocate a `written(ndof)` mask in `read_vel_sol_single` and `read_q2_sol_single`.
2. On the first write to a local dof, store the value and mark it.
3. On later occurrences of the same dof, skip the write.

Pseudo-pattern:

```fortran
if (.not. written(jvt)) then
  u(jvt) = buf(ivt)
  written(jvt) = .true.
end if
```

The same logic applies to `v`, `w`, and generic Q2 field components.

## Important Note

The right long-term solution is a clearly defined ownership rule for repartitioned dump loading.

A first-write mask is the simplest immediate stabilization, but the chosen policy should be made explicit and shared consistently between:

- dump writing
- repartition mapping
- repartition reading

## Recommended Follow-Up

1. Patch `read_vel_sol_single`.
2. Patch `read_q2_sol_single`.
3. Add a regression test for repartitioned Q2 loading.
4. Recheck applications that use `SolFromFileRepart`.
5. Verify that loaded Q2 dofs are not multi-written without a defined rule.
