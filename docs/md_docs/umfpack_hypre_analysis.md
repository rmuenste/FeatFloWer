# UMFPACK vs HYPRE solver interface analysis

## Context from the refactoring plan

The refactoring plan splits solver interfaces into dedicated modules: `QuadSc_solver_hypre.f90`
for HYPRE and `QuadSc_solver_coarse.f90` for UMFPACK/coarse logic. Solver modules are
restricted to `var_QuadScalar`, `pp3d_mpi`, and solver libraries to avoid circular dependencies,
and UMFPACK logic is extracted from `Create_CMat` into a dedicated coarse solver module.

## UMFPACK solver interface analysis

### Dedicated module and entry point

UMFPACK coarse-grid setup is implemented in `QuadSc_solver_coarse.f90` with the public entry
point `Setup_UMFPACK_CoarseSolver`, documented as the extracted UMFPACK coarse solver
interface.

### Execution model and process ownership

UMFPACK setup runs only on the coarse-grid process (`myid == 0`); other ranks return
immediately, making the coarse solver a root-rank operation.

### Solver types handled

`Setup_UMFPACK_CoarseSolver` supports solver types 2, 3, and 4:
- Type 2: full coarse grid matrix, no geometric coarsening, with LU factorization.
- Type 3: geometric coarsening (/16) without factorization.
- Type 4: geometric coarsening (/16) with LU factorization.

### bNoOutflow (singular matrix) handling

The routine accepts `bNoOutFlow` and, when set, zeroes the row slice between
`lMat_work%LdA(1)+1` and `lMat_work%LdA(2)-1` before factorization or coarsening, applying the
singular matrix fix for no-outflow configurations on the UMFPACK path.

### Data path and coarse structure construction

- The routine activates the coarse level via `ILEV = coarse_lev` and `SETLEV`.
- Type 2 copies the full CSR structure into `UMF_lMat` and factors `UMF_CMat`.
- Types 3/4 apply stride-4 access to build the /16-reduced CSR (`crsSTR%A`), including row
  pointer construction and column/value extraction from the fine matrix.
- Boundary constraints (`knprP`) are mapped onto the coarse diagonal with a small value
  (`1d-16`).
- LU factorization is performed only for type 4.

### Call site

`Setup_UMFPACK_CoarseSolver` is invoked from `Create_CMat` in `QuadSc_assembly.f90` and
explicitly receives `bNoOutflow`.

## Differences vs the HYPRE solver interface

### Parallel ownership model

- UMFPACK runs only on `myid == 0`.
- HYPRE runs on worker ranks (`myid != 0`) and requires collective numbering operations.

### Solver type mapping

- UMFPACK uses types 2/3/4 (with /16 geometric coarsening for 3/4).
- HYPRE uses types 7/8 with separate full and geometric coarse setup routines.

### Data structures and conversion work

- UMFPACK copies or coarsens CSR data into `UMF_lMat` or `crsSTR%A` and factors with
  `myUmfPack_Factorize`.
- HYPRE constructs global numbering, off-partition indices, and HYPRE-specific CSR arrays,
  then handles zero-based conversion and column sorting.

### bNoOutflow (singular fix) differences

- UMFPACK applies the fix on the coarse rank by zeroing a row slice before factorization.
- HYPRE applies the fix on worker ranks; in full coarse it uses the element link mask, and in
  geometric coarse it limits the fix to `myid == 1` and the first row slice.

### Integration point

- UMFPACK is called from `Create_CMat` in `QuadSc_assembly.f90`.
- HYPRE is called from `SetUp_HYPRE_Solver` in `QuadSc_def.f90`.
