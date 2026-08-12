# QuadSc Current Implementation

This document describes the current organization and key conventions of the
QuadSc matrix and solver implementation after the 2025 structural refactoring.

## Module Organization

The public compatibility module remains:

- `source/src_quadLS/QuadSc_def.f90` (`def_QuadScalar`)

It acts as a facade and re-exports routines from focused implementation modules:

| Module | Responsibility |
|---|---|
| `QuadSc_struct.f90` | CSR sparsity structures for Q2, P1, and parallel matrices |
| `QuadSc_assembly.f90` | Mass, diffusion, convection, stabilization, coupling, and pressure-matrix assembly |
| `QuadSc_solver_hypre.f90` | Full and geometrically coarsened HYPRE coarse-level setup |
| `QuadSc_solver_coarse.f90` | UMFPACK coarse-grid pressure solver setup |
| `QuadSc_var.f90` | Shared QuadSc types, matrices, pointers, and runtime state |
| `QuadSc_mg.f90` | Multigrid solver machinery |

Existing callers can continue to `USE def_QuadScalar`; new internal work should
be placed in the focused module that owns the behavior.

## CSR Matrix Representation

`TMatrix` in `source/src_util/types.f90` stores the CSR structure:

```fortran
TYPE TMatrix
  INTEGER :: nu, na
  INTEGER, DIMENSION(:), ALLOCATABLE :: ColA, LdA
END TYPE
```

- `nu` is the row count.
- `na` is the number of stored entries.
- `LdA(1:nu+1)` contains one-based row offsets.
- `ColA(1:na)` contains one-based column indices.

Matrix values are stored separately in level-dependent arrays such as
`mg_MMat`, `mg_DMat`, `mg_KMat`, `mg_BXMat`, and related structures.

For row `i`, entries occupy:

```fortran
LdA(i):LdA(i+1)-1
```

## Multigrid Conventions

- `NLMIN` is the coarsest active level.
- `NLMAX` is the finest active level.
- `ILEV` identifies the level currently selected by `SETLEV(2)`.
- `mg_mesh%level(ILEV)` owns the mesh topology for that level.
- `mg_*` matrix arrays store one matrix or value structure per level.

The Q2 velocity-space degree count is:

```text
NVT + NET + NAT + NEL
```

corresponding to corner vertices, edge midpoints, face midpoints, and element
centers of the 27-node hexahedral Q2 element.

Pressure uses a discontinuous P1 representation with four coefficients per
element in the main Q2/P1 solver path.

## Structure Creation

`QuadSc_struct` exports:

- `Create_QuadMatStruct`
- `Create_QuadLinMatStruct`
- `Create_LinMatStruct`
- `Create_ParLinMatStruct`

These routines allocate and populate CSR patterns. They still rely on shared
mesh and solver state from `var_QuadScalar` and `PP3D_MPI`.

## Matrix Assembly

`QuadSc_assembly` exports generic helpers:

- `Assemble_Mass_Generic`
- `Assemble_Diffusion_Alpha_Generic`
- `Assemble_ParallelMatrix_Generic`

and the established high-level entry points:

- `Create_MRhoMat`
- `Create_MMat`
- `Create_hDiffMat`
- `Create_ConstDiffMat`
- `Create_DiffMat`
- `Create_BMat`
- `Create_CMat`
- `Create_SMat`
- `Create_KMat`

The generic helpers remove repeated level loops while retaining the existing
element kernels and shared-state conventions.

## Coarse Solvers

### HYPRE

`QuadSc_solver_hypre` provides:

- `Setup_HYPRE_CoarseLevel_Full`
- `Setup_HYPRE_CoarseLevel_Geometric`

The interface converts FeatFloWer's one-based CSR structures into the numbering
and distribution required by HYPRE.

### UMFPACK

`QuadSc_solver_coarse` provides:

- `Setup_UMFPACK_CoarseSolver`

The routine supports the existing coarse solver modes and delegates numerical
factorization to `source/UmfpackSolver.f90`.

The geometric coarse representation uses every fourth pressure row and column.
Consequently, its unknown count is reduced by four and its structured nonzero
count by sixteen for the supported layout. This is an implementation invariant,
not a general sparse-matrix coarsening rule.

## Matrix Update Groups

The parameter commonly called the "Matrix Renewal scheme" assigns matrices to
update groups; the numbers are not timestep frequencies.

- Group `1`: generated during initialization.
- Group `2`: optional secondary update hook.
- Group `3`: generated during the timestep update path.
- Group `0`: excluded from regeneration.

`OperatorRegenaration` in `QuadSc_corrections.f90` compares each matrix's group
assignment with the requested group and invokes the relevant `Create_*`
routine.

A common constant-viscosity configuration assigns mass, diffusion, and
coupling matrices to group 1 and convection to group 3.

## Compatibility Boundary

`def_QuadScalar` intentionally re-exports the extracted routines. This keeps
application and solver call sites stable while allowing implementation modules
to evolve independently.

The facade still contains compatibility logic and some high-level routines, so
the refactoring is structurally complete but not a full removal of global
state.

## Known Technical Debt

- Heavy dependence on mutable module state in `var_QuadScalar`.
- Manual allocation and reallocation without a uniform error policy.
- Fixed preallocation estimates and layout-specific constants.
- Sparse direct and HYPRE setup tied closely to global numbering conventions.
- Limited focused tests for structure creation and assembly.
- Legacy naming such as `OperatorRegenaration` and "Matrix Renewal scheme".
- Remaining implementation in the facade that could move to narrower modules.

## Related Documentation

- [matrix_structures_guide.md](matrix_structures_guide.md)
- [matrix_assembly_guide.md](matrix_assembly_guide.md)
- [hypre_csr_analysis.md](hypre_csr_analysis.md)
- [umfpack_hypre_analysis.md](umfpack_hypre_analysis.md)
- [quadsc_refactoring_status.md](quadsc_refactoring_status.md)
