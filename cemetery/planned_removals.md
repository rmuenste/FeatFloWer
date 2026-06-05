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
