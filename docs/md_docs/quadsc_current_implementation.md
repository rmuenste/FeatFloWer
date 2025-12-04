# QuadSc Current Implementation Documentation

**Purpose**: Document the current state of QuadSc_def.f90 before refactoring
**Created**: Phase 0 Baseline (December 2025)
**Status**: Living document - update as understanding improves

---

## Overview

The `QuadSc_def.f90` module (~4500 lines) implements matrix structure creation, assembly, and solving for quadratic scalar transport equations on hierarchical multigrid Q2/P1 finite element meshes. It handles:

- Q2 (quadratic, 27-point stencil) scalar matrices
- Q2-P1 coupling matrices (B and B^T)
- Multigrid hierarchy management (NLMIN to NLMAX)
- HYPRE algebraic multigrid solver interface
- UMFPACK direct coarse grid solver
- Parallel MPI matrix assembly

---

## Matrix Storage Format: CSR (Compressed Sparse Row)

### TMatrix Structure

Defined in `source/src_util/types.f90:358-361`:

```fortran
TYPE TMatrix
  INTEGER :: nu, na
  INTEGER, DIMENSION(:), ALLOCATABLE :: ColA, LdA
END TYPE
```

**Fields**:
- `nu`: Number of unknowns (rows)
- `na`: Number of allocated non-zero entries
- `ColA(1:na)`: Column indices of non-zero entries (1-based Fortran indexing)
- `LdA(1:nu+1)`: Row pointer array (CSR format)

**CSR Format Explanation**:

The CSR (Compressed Sparse Row) format stores a sparse matrix using three arrays:
1. **LdA (Lead Array)**: Size `nu+1`, where `LdA(i)` points to the first entry of row `i` in ColA
2. **ColA (Column Array)**: Size `na`, stores column indices of non-zero entries
3. **Actual values**: Stored separately (e.g., `DMat`, `KMat`, `CMat`)

**Example**: 3×3 matrix with 5 non-zeros
```
Matrix:          CSR Storage:
[1.0  0   2.0]   LdA  = [1, 3, 5, 6]
[3.0  4.0  0 ]   ColA = [1, 3, 1, 2, 1]
[5.0  0    0 ]   Values=[1.0, 2.0, 3.0, 4.0, 5.0]
```

Row `i` spans entries from `LdA(i)` to `LdA(i+1)-1` in ColA and Values arrays.

**Accessing row `i`**:
```fortran
DO j = LdA(i), LdA(i+1)-1
  column_index = ColA(j)
  value = SomeMatrixArray(j)  ! e.g., DMat(j), CMat(j)
END DO
```

**Important**: The actual matrix values are NOT stored in TMatrix. They are stored in separate arrays like:
- `DMat(:)` - Diffusion matrix values
- `KMat(:)` - Stiffness matrix values
- `CMat(:)` - Combined/assembled matrix values
- `BXMat(:), BYMat(:), BZMat(:)` - Coupling matrices

---

## Multigrid Hierarchy Conventions

### Level Indexing

**Global Parameters** (in `var_QuadScalar`):
- `NLMIN`: Coarsest multigrid level (typically 1 or 2)
- `NLMAX`: Finest multigrid level (typically 3-7 depending on refinement)
- `ILEV`: Current active level (used in loops and SETLEV calls)

**Convention**:
- Level 1 (NLMIN): Coarsest mesh (~100-1000 elements)
- Level NLMAX: Finest mesh (~100K-10M elements)
- Refinement: Each level has ~8× more elements than previous (3D octree refinement)

**Multigrid Arrays**:
```fortran
TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_qMat   ! (NLMIN:NLMAX)
TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_lMat   ! (NLMIN:NLMAX)
TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_qlMat  ! (NLMIN:NLMAX)
TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_lqMat  ! (NLMIN:NLMAX)
```

**Active Level Pointer**:
```fortran
qMat => mg_qMat(NLMAX)  ! Active Q2 matrix pointer (set after structure creation)
lMat => mg_lMat(NLMAX)  ! Active P1 matrix pointer
```

### SETLEV(2) Call

The `CALL SETLEV(2)` routine (from FEAT library) sets the active mesh level:
- Sets global mesh pointers to `mg_mesh%level(ILEV)`
- Updates element/vertex/edge counts for current level
- Prepares connectivity arrays (kvert, kedge, karea)

**Pattern** (seen throughout QuadSc_def.f90:18-62):
```fortran
DO ILEV = NLMIN, NLMAX
  CALL SETLEV(2)         ! Activate level ILEV
  ! ... work on level ILEV ...
  CALL Create_Structure_For_Level_ILEV()
END DO
```

---

## Degree of Freedom Calculation

### Q2 Scalar System (27-node hexahedral element)

**DOF Count** (`QuadSc_def.f90:29-32`):
```fortran
ndof = mg_mesh%level(ilev)%nvt +  ! vertices (8 per element)
       mg_mesh%level(ilev)%net +  ! edges (12 per element)
       mg_mesh%level(ilev)%nat +  ! faces (6 per element)
       mg_mesh%level(ilev)%nel    ! element centers (1 per element)
```

Total: **27 nodes per Q2 hexahedral element** (8 corners + 12 edges + 6 faces + 1 center)

### P1 Scalar System (4-node tetrahedral or 8-node hexahedral)

**DOF Count** (for pressure, `QuadSc_def.f90:1141`):
```fortran
ndof_P1 = mg_mesh%level(ilev)%nel    ! For P1 discontinuous (cell-centered)
! OR
ndof_P1 = mg_mesh%level(ilev)%nvt    ! For P1 continuous (vertex-based)
```

In FeatFloWer, pressure uses **P1 discontinuous** (one DOF per element).

---

## Matrix Structure Allocation

### Q2 Matrix Structure (`Create_QuadMatStruct`, line 18)

**Initial Estimate** (`QuadSc_def.f90:34`):
```fortran
MatSize = 300 * NDOF
nERow = 300              ! Estimated entries per row
```

**Why 300?**
- Q2 element has 27 nodes
- Each node couples to ~10-15 neighboring elements in 3D
- Worst case: ~27 × 10 = 270 entries per row
- 300 provides safety margin

**Actual Structure** (created by `AP7` routine):
- `AP7` traverses all elements and builds sparsity pattern
- Only actually connected entries are stored
- Final `na` is typically much less than initial estimate

### Q2-P1 Coupling Structure (`Create_QuadLinMatStruct`, line 66)

**B Matrix (Q2 to P1 divergence)** (`QuadSc_def.f90:84`):
```fortran
MatSize = 16 * 27 * mg_mesh%level(ilev)%nel
nERow = 16
```

**Why 16×27?**
- 27 Q2 nodes per element
- Up to 16 entries per Q2 node when coupling to P1 (conservative estimate)
- After structure creation, na is multiplied by 4 (line 99) for x/y/z components + padding

**B^T Matrix (P1 to Q2 gradient)** (`QuadSc_def.f90:111`):
```fortran
MatSize = 4 * 27 * mg_mesh%level(ilev)%nel
nERow = 27
```

---

## The /16 Mystery

### Location

**File**: `QuadSc_def.f90:1142`
**Context**: UMFPACK coarse grid solver setup

```fortran
IF (coarse_solver.EQ.4.OR.coarse_solver.EQ.3) THEN
  crsSTR%A%nu = lMat%nu/4
  crsSTR%A%na = lMat%na/16   !!!! /16????????????????
```

### Current Understanding (Hypothesis)

**Geometric Coarsening from P1 to Coarser P1**:

When setting up the UMFPACK coarse grid solver, the code creates a coarsened version of the P1 pressure matrix `lMat`:

1. **Unknowns coarsening** (`nu/4`):
   - P1 discontinuous: 1 DOF per element
   - Geometric coarsening in 3D: Each coarse element combines 2×2×2 = 8 fine elements
   - But FeatFloWer appears to use a 4:1 coarsening ratio (possibly 2D or structured)
   - Thus: `coarse_nu = fine_nu / 4`

2. **Non-zeros coarsening** (`na/16`):
   - If unknowns reduce by 4×, and sparsity pattern also becomes sparser
   - **Row reduction**: 4× fewer rows
   - **Column reduction**: 4× fewer columns per row (assuming similar connectivity)
   - **Total**: 4 × 4 = 16× fewer non-zero entries
   - Thus: `coarse_na ≈ fine_na / 16`

**Evidence from code** (`QuadSc_def.f90:1151-1164`):
```fortran
DO i = 1, crsSTR%A%nu
  j = 4*i - 3                                    ! Take every 4th fine row
  crsSTR%A%LdA(i+1) = crsSTR%A%LdA(i) + (lMat%LdA(j+1)-lMat%LdA(j))/4
END DO

DO i = 1, crsSTR%A%nu
  j = 4*(i-1) + 1
  DO jCol = lMat%LdA(j), lMat%LdA(j+1)-1, 4     ! Take every 4th column
    iEntry = iEntry + 1
    crsSTR%A%ColA(iEntry) = (lMat%ColA(jCol)-1)/4 + 1
    crsSTR%A_MAT(iEntry) = CMat(jCol)
  END DO
END DO
```

The code explicitly:
- Takes every 4th row (stride 4 in row loop)
- Takes every 4th column entry (stride 4 in jCol loop)
- Result: 4 × 4 = 16× reduction in entries

**Why /16 is not always exact**:
- Boundary effects (boundary rows may have different connectivity)
- Irregular meshes may not have exactly 4:1 coarsening
- The `/16` is a **conservative allocation estimate**, actual count is determined by the loop

### Diagnostic Strategy (Phase 0.4)

Add logging to verify this hypothesis:
```fortran
IF (myid == 0 .and. ILEV == NLMIN) THEN
  WRITE(*,'(A)') '========== UMFPACK Coarse Grid Diagnostic =========='
  WRITE(*,'(A,I0)') 'Fine grid unknowns (lMat%nu):     ', lMat%nu
  WRITE(*,'(A,I0)') 'Fine grid non-zeros (lMat%na):    ', lMat%na
  WRITE(*,'(A,I0)') 'Coarse unknowns (lMat%nu/4):      ', lMat%nu/4
  WRITE(*,'(A,I0)') 'Coarse non-zeros (lMat%na/16):    ', lMat%na/16
  WRITE(*,'(A,I0)') 'Actual coarse non-zeros (iEntry): ', iEntry
  WRITE(*,'(A,F10.4)') 'nu ratio (fine/coarse):           ', REAL(lMat%nu)/REAL(lMat%nu/4)
  WRITE(*,'(A,F10.4)') 'na ratio (fine/coarse):           ', REAL(lMat%na)/REAL(iEntry)
  WRITE(*,'(A)') '===================================================='
END IF
```

---

## Matrix Assembly Overview

### Element Assembly Pattern

**Q2 Stiffness/Mass Matrix** (typical pattern):
```fortran
DO iel = 1, nel
  ! 1. Get element connectivity
  CALL E013(..., kve, ...)  ! Returns 27 local node indices for element iel

  ! 2. Compute local element matrix (27×27)
  CALL Local_Element_Routine(iel, A_local)

  ! 3. Assemble into global CSR matrix
  DO i = 1, 27
    iglob = kve(i)                          ! Global node index
    DO j = LdA(iglob), LdA(iglob+1)-1
      jglob = ColA(j)                       ! Global column index
      ! Find local index k such that kve(k) == jglob
      DMat(j) = DMat(j) + A_local(i, k)    ! Accumulate
    END DO
  END DO
END DO
```

**MPI Parallelization**:
- Each process owns a subset of elements
- Ghost elements on processor boundaries require MPI communication
- `COMM_SUMMN` and `COMM_Maximum` synchronize matrix values across processes

---

## HYPRE Solver Interface

### Matrix Transfer to HYPRE (`QuadSc_def.f90:684-689`)

```fortran
myHYPRE%nrows = lPMat%nu/4           ! Coarsen P1 matrix
myHYPRE%ilower = (myHYPRE%ilower+3)/4
myHYPRE%iupper = myHYPRE%iupper/4
myHYPRE%nonzeros = lPMat%na/16 + lMat%na/16  ! Both matrices contribute
```

**HYPRE Setup**:
1. Coarsen the P1 matrix by 4× (geometric coarsening)
2. Estimate non-zeros using /16 rule
3. Transfer matrix to HYPRE's IJMatrix format
4. HYPRE performs algebraic multigrid (BoomerAMG)

**Note**: HYPRE uses 0-based C indexing, FeatFloWer uses 1-based Fortran indexing. Conversion happens in HYPRE interface routines.

---

## Subroutine Organization (Current State)

### Matrix Structure Creation
- `Create_QuadMatStruct()` (line 18): Q2 scalar matrix structure
- `Create_QuadLinMatStruct()` (line 66): Q2-P1 coupling structures (B, B^T)

### Matrix Assembly
- `Create_M_Matrix()`: Mass matrix assembly
- `Create_D_Matrix()`: Diffusion matrix assembly
- `Create_K_Matrix()`: Convection matrix assembly
- `Create_C_Matrix()`: Combined system matrix assembly
- `Create_A_MKDC()`: General assembly routine

### Solver Setup
- `Create_CoarseSolverStructure()`: UMFPACK/HYPRE coarse grid setup
- `myUmfPack_Factorize()`: Direct LU factorization

### MPI Communication
- Numerous `COMM_SUMMN` calls for matrix synchronization
- `E011Sum`, `E011DMat`: Element-level MPI communication

---

## Known Issues (To Be Addressed in Refactoring)

1. **Unsafe Memory Reallocation** (lines throughout):
   ```fortran
   IF (ALLOCATED(arr)) DEALLOCATE(arr)  ! No STAT checking
   ALLOCATE(arr(n))                     ! No STAT checking
   ```

2. **Magic Numbers**:
   - `300` - Q2 matrix entries per row estimate (line 34)
   - `16` - Q2-P1 coupling estimate (line 91)
   - `/16` - Coarse grid sparsity estimate (line 1142)

3. **Code Duplication**:
   - Similar assembly patterns repeated for M/D/K/A matrices
   - Nearly identical structure creation for different element types

4. **Mixed Responsibilities**:
   - Matrix structure + assembly + solver setup all in one 4500-line file
   - Difficult to unit test individual components

5. **Lack of Error Handling**:
   - No propagation of ALLOCATE errors
   - No validation of matrix dimensions before operations

---

## References

### Related Files
- `source/src_util/types.f90`: TMatrix definition
- `source/src_quadLS/QuadSc_var.f90`: Global variables and parameters
- `source/src_quadLS/QuadSc_mg.f90`: Multigrid solver
- `source/UMFPackSolver.f90`: UMFPACK interface

### Element Routines
- `E013`: Q2 hexahedral element (27 nodes)
- `E010`: P1 element (8 nodes for continuous, 1 for discontinuous)
- `E011`: Communication element routine

### FEAT Library Calls
- `AP7`: Assemble Q2 matrix structure (square)
- `AP9`: Assemble rectangular coupling structure (Q2-P1)
- `SETLEV`: Set active multigrid level

---

**Last Updated**: December 4, 2025 (Phase 0.3 - Initial documentation)
