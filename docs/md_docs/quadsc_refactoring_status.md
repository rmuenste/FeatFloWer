# QuadSc Refactoring Status

This document consolidates the former QuadSc refactoring, modernization, and
phase-roadmap notes. It records what is complete and what remains; it is not a
step-by-step historical implementation diary.

## Completed Work

The structural refactoring was completed in phases during December 2025:

1. HYPRE setup moved to `QuadSc_solver_hypre.f90`.
2. UMFPACK coarse setup moved to `QuadSc_solver_coarse.f90`.
3. Repeated mass, diffusion, and parallel assembly loops gained generic helpers.
4. CSR structure creation moved to `QuadSc_struct.f90`.
5. matrix assembly moved to `QuadSc_assembly.f90`.
6. `QuadSc_def.f90` became a compatibility facade that re-exports the existing
   public entry points.
7. The `/16` coarse-matrix allocation rule was diagnosed and documented as the
   result of the supported four-by-four row/column subsampling.

The current architecture is documented in
[quadsc_current_implementation.md](quadsc_current_implementation.md).

## Current Dependency Direction

The intended dependency direction is:

```text
var_QuadScalar
  ├── QuadSc_struct
  ├── QuadSc_solver_hypre
  └── QuadSc_solver_coarse

QuadSc_struct + solver modules
  └── QuadSc_assembly

focused modules
  └── def_QuadScalar facade
```

In practice, the focused modules still depend on shared state and MPI globals.
Avoid adding reverse dependencies from lower-level modules back into the facade.

## Remaining Modernization Work

### State and Interfaces

- Reduce direct access to mutable globals in `var_QuadScalar`.
- Pass mesh, matrix, and solver objects explicitly where this produces a clear
  testing or ownership benefit.
- Keep compatibility wrappers in `def_QuadScalar` while call sites migrate.

### Allocation Safety

- Centralize checked allocation and reallocation patterns.
- Propagate allocation failures with context instead of relying on unchecked
  `ALLOCATE` and `DEALLOCATE`.
- Validate dimensions before reusing level-dependent arrays.

### Sparse Structure Construction

- Replace broad fixed preallocation estimates where count-first construction is
  practical.
- Document every layout-specific constant near the code that depends on it.
- Preserve the `/16` rule only where the four-by-four geometric mapping is
  explicitly guaranteed.

### Error Handling

- Introduce a consistent solver error-reporting convention.
- Attach rank, level, matrix family, and allocation size to fatal diagnostics.
- Make HYPRE and UMFPACK failures visible at the owning abstraction boundary.

### Naming

- Correct `OperatorRegenaration` without breaking callers.
- Rename "Matrix Renewal scheme" in user-facing documentation to clarify that
  it assigns update groups rather than frequencies.

### Testing

Focused tests are still needed for:

- CSR row-pointer and column-index invariants,
- matrix dimensions across levels,
- generic assembly equivalence with established paths,
- UMFPACK coarse extraction,
- HYPRE numbering conversion, and
- sequential versus MPI matrix consistency.

The repository's benchmark and CTest workflows remain the practical regression
barrier until narrower tests exist.

## Change Rules

When changing QuadSc structure, assembly, or coarse-solver code:

1. Keep the facade API stable unless the migration is part of the same change.
2. Build with bounds checking during development when the compiler supports it.
3. Verify at least one sequential and one MPI benchmark.
4. Compare solver convergence and representative output metrics.
5. Record required rank counts and CMake options in the change description.
6. Measure performance when changing level loops, communication, or sparse
   allocation.

## Deferred Opportunities

These are valid research or modernization topics, but are not committed plans:

- OpenMP parallelization of element assembly,
- stronger derived-type ownership for complete solver state,
- replacement of legacy COMMON/global state,
- richer performance instrumentation, and
- alternative sparse construction strategies.

Implement these only with benchmark coverage and an explicit ownership model.
