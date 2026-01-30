# HYPRE CSR Interface Analysis (C Matrix)

This note documents how the internal CSR data structures for the pressure Schur complement (C matrix) are built and how they are converted into HYPRE/BoomerAMG structures on the coarse level.

## 1) Internal CSR data model

The pressure system uses the generic CSR type `TMatrix` from `source/src_util/types.f90`:

- `nu`: number of rows (pressure DOFs)
- `na`: number of nonzeros
- `LdA(1:nu+1)`: row pointer (1-based)
- `ColA(1:na)`: column indices (1-based)

The CSR *values* are stored separately:

- Local values: `mg_CMat(ILEV)%a`
- Parallel/off-partition values: `mg_CPMat(ILEV)%a`

The CSR *structure* (row pointers and column indices) is stored in:

- `mg_lMat(ILEV)` for the local C matrix
- `mg_lPMat(ILEV)` for the parallel/off-partition C matrix

Pointers are held in `source/src_quadLS/QuadSc_var.f90` and are bound per level by `SETLEV(2)`.

## 2) Building the C-matrix CSR structure

### Local structure (`mg_lMat`)

Built in `Create_LinMatStruct` in `source/src_quadLS/QuadSc_struct.f90`:

1. `Get_CMatLen(..., iP=1)` computes `na` (total nonzeros).
2. `mg_lMat(ILEV)%nu = 4 * nel` (P1: 4 DOFs per hexahedron).
3. `Get_CMatStruct(..., iP=1)` fills `LdA` and `ColA`.

`Get_CMatStruct` (in `source/src_quadLS/QuadSc_proj.f`) uses vertex/element adjacency to assemble a 4x4 block for each element and its neighbor elements, sorted in ascending order.

### Parallel/off-partition structure (`mg_lPMat`)

Built in `Create_ParLinMatStruct` in `source/src_quadLS/QuadSc_struct.f90`:

1. `Get_CMatLen(..., iP=2)` for off-partition adjacency only.
2. `mg_lPMat(ILEV)%nu = 4 * NEL`.
3. `Get_CMatStruct(..., iP=2)` fills `LdA` and `ColA` for off-process contributions.

## 3) Assembling the C-matrix values

Values are computed via `Get_CMat` in `source/src_quadLS/QuadSc_proj.f`:

- Formula: `C = B^T * M^{-1} * B`
- Uses the P1/Q2 coupling matrices `B` and `B^T`.
- Uses `MlRhoPmat` as the lumped mass matrix for inversion.
- Dirichlet velocity DOFs are excluded by setting inverse mass to 0 if `KNPRU/V/W == 1`.

Assembly flow:

- `Create_CMat` in `source/src_quadLS/QuadSc_assembly.f90` allocates `mg_CMat(ILEV)%a` and calls `Get_CMat` for each level.
- Parallel values are assembled into `mg_CPMat(ILEV)%a` in `source/src_quadLS/QuadSc_def.f90` using the parallel B matrices and `lPMat` structure.

## 4) Coarse-level selection for HYPRE

HYPRE setup always uses the coarsest multigrid level:

- `ILEV = lScalar%prm%MGprmIn%MinLev`
- `CALL SETLEV(2)` updates global level state across all ranks

Pointers are then bound:

- `lMat => mg_lMat(ILEV)` and `CMat => mg_CMat(ILEV)%a` (local)
- `lPMat => mg_lPMat(ILEV)` and `CPMat => mg_CPMat(ILEV)%a` (parallel)

## 5) Conversion to HYPRE structures (Full matrix, type 7)

`Setup_HYPRE_CoarseLevel_Full` in `source/src_quadLS/QuadSc_solver_hypre.f90` performs the conversion.

### Global numbering

- `GetMyHYPRENumberingLimits` computes `ilower/iupper` per rank using MPI
- `Numbering(1:nrows)` becomes contiguous global IDs for local DOFs

### Off-partition numbering

- Scan `lPMat` to find maximum remote column ID
- Allocate `OffPartitionNumbering` and fill via `GetHYPREParPressureIndices`

### HYPRE row-wise arrays

Allocates and fills:

- `ncols(i)` = nnz in row i (local + parallel)
- `rows(i)` = global row index (`Numbering(i)`)
- `cols(1:nonzeros)` and `values(1:nonzeros)` from local + parallel data

Columns are mapped to global indices using `Numbering` (local) and `OffPartitionNumbering` (remote). Each row is sorted by column index via `SORT_DOFs`.

### Index base

If `myHYPRE%ZeroBased` is `.true.`, `ilower/iupper`, `rows`, and `cols` are shifted to 0-based indexing for HYPRE.

## 6) Geometric coarsening path (type 8)

`Setup_HYPRE_CoarseLevel_Geometric` applies explicit coarsening:

- DOFs reduced by /4: `nrows = lPMat%nu / 4`
- Nonzeros reduced by /16: `nonzeros = lPMat%na/16 + lMat%na/16`

It then extracts a coarse stencil using stride-4 access:

- Row mapping: `JEQ = 4*(IEQ-1)+1`
- Column mapping: `(global_index + 3) / 4`

Local and parallel contributions are merged, sorted, and optionally converted to 0-based indexing as in the full path. A diagnostic block reports actual nnz vs allocated nnz for the /16 assumption.

## 7) No-outflow (singular) matrix handling

The coarse C matrix can be singular in no-outflow configurations. Both HYPRE setup routines apply a localized fix before conversion:

- **Full path**: zeroes off-diagonal entries on selected rows associated with `coarse%myELEMLINK(iel) == 1`.
- **Geometric path**: zeroes off-diagonal entries for the first row on rank 1.

This imposes a Dirichlet pressure anchor required for convergence with BoomerAMG.

---

### Key files and entry points

- CSR structure: `source/src_quadLS/QuadSc_struct.f90`
- CSR values: `source/src_quadLS/QuadSc_assembly.f90`, `source/src_quadLS/QuadSc_def.f90`
- Element kernels: `source/src_quadLS/QuadSc_proj.f`
- HYPRE conversion: `source/src_quadLS/QuadSc_solver_hypre.f90`
- HYPRE numbering: `source/src_quadLS/QuadSc_mpi.f90`
- Global type definitions: `source/src_util/types.f90`

## 8) HYPRE/BoomerAMG call sites (coarse solve)

The coarse-grid solve dispatch for HYPRE happens in `mgCoarseGridSolver_P` inside
`source/src_quadLS/QuadSc_mg.f90` under `#ifdef HYPRE_AVAIL` and when
`MyMG%CrsSolverType` is 7 or 8:

- **Type 7 (full matrix, no geometric coarsening)**:
  - On worker ranks (`myid != 0`), `myHypre%rhs` and `myHypre%sol` are set from
    `myMG%B(mgLev)%x` and `myMG%X(mgLev)%x`.
  - The solve is executed via `CALL myHypreGMRES_Solve(CoarseIter)` (the older
    `myHypre_Solve` and `myHyprePCG_Solve` calls are present but commented out).
  - The solution is copied back into `myMG%X(mgLev)%x` on workers.

- **Type 8 (geometric coarsening)**:
  - Builds a coarse RHS/initial guess by sampling every 4th P1 DOF:
    `j = 4*(i-1)+1` and copies into `myHypre%rhs/sol`.
  - Calls `CALL myHypre_Solve(CoarseIter)` (not GMRES here).
  - Post-processes the solution (scaled by 3) and interpolates back to P1 via
    `INTPVB` and `IntQ1toP1`, then applies an SOR smooth.

This is the runtime point where the CSR data prepared in
`Setup_HYPRE_CoarseLevel_*` is consumed by the HYPRE solver wrappers.
