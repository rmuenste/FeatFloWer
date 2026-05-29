# Refactoring Plan: QuadSc_def.f90

**Date:** November 28, 2025
**Status:** Proposal
**Last Updated:** November 28, 2025 (Enhanced with dependency graph and testing framework)

The module `def_QuadScalar` in `source/src_quadLS/QuadSc_def.f90` has grown into a monolithic file (~4500 lines) mixing matrix allocation, assembly, solver interfacing, and coarse grid logic. This document outlines a strategy to split it into cohesive compilation units and reduce code duplication.

## 1. Proposed Module Split

We recommend breaking `QuadSc_def.f90` into 5 smaller, focused modules. This improves compile times, readability, and testability.

| New Module Name | Content / Responsibility | Key Subroutines |
| :--- | :--- | :--- |
| **`QuadSc_struct.f90`** | **Sparsity Patterns**. Allocation of CSR structures (`ColA`, `LdA`). | `Create_QuadMatStruct`<br>`Create_QuadLinMatStruct`<br>`Create_LinMatStruct`<br>`Create_QuadLinParMatStruct` |
| **`QuadSc_assembly.f90`** | **Matrix Assembly**. Filling matrix values (`%a`). | `Create_MMat`<br>`Create_DiffMat`<br>`Create_KMat`<br>`Create_SMat`<br>`Create_BMat` |
| **`QuadSc_solver_hypre.f90`** | **HYPRE Interface**. Setup and data conversion for HYPRE. | `SetUp_HYPRE_Solver` |
| **`QuadSc_solver_coarse.f90`** | **Coarse Grid**. UMFPACK factorization and coarse matrix logic. | *Extract from* `Create_CMat`<br>`Create_CoarseSolver_UMFPACK` (new) |
| **`QuadSc_system.f90`** | **System Definition**. High-level defect and matrix definition routines. | `Matdef_general_QuadScalar`<br>`Matdef_General_LinScalar` |

### 1.1 Module Dependency Graph (NEW)

To avoid circular dependencies, modules must follow a strict hierarchical dependency structure:

```
                    ┌─────────────────────┐
                    │  var_QuadScalar     │  (Existing global state)
                    │  (types & globals)  │
                    └──────────┬──────────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                ▼
    ┌──────────────┐  ┌──────────────┐  ┌──────────────┐
    │ QuadSc_      │  │ QuadSc_      │  │ QuadSc_      │
    │ struct.f90   │  │ solver_      │  │ solver_      │
    │              │  │ hypre.f90    │  │ coarse.f90   │
    └──────┬───────┘  └──────────────┘  └──────────────┘
           │
           │ (depends on structures)
           ▼
    ┌──────────────┐
    │ QuadSc_      │
    │ assembly.f90 │
    └──────┬───────┘
           │
           │ (high-level coordination)
           ▼
    ┌──────────────┐
    │ QuadSc_      │
    │ system.f90   │
    └──────┬───────┘
           │
           │ (facade re-exports)
           ▼
    ┌──────────────┐
    │ def_         │
    │ QuadScalar   │  (Facade module)
    │ (facade)     │
    └──────────────┘
           │
           │ (consumers unaffected)
           ▼
    ┌──────────────┐
    │ QuadSc_      │
    │ main.f90     │
    │ & others     │
    └──────────────┘
```

### 1.2 Dependency Rules (NEW)

**CRITICAL: These rules MUST be enforced to prevent circular dependencies:**

1. **Lower modules CANNOT depend on higher modules** in the hierarchy above.

2. **Allowed dependencies**:
   - `QuadSc_struct` may USE: `var_QuadScalar`, `pp3d_mpi`
   - `QuadSc_assembly` may USE: `QuadSc_struct`, `var_QuadScalar`, `pp3d_mpi`, element routines (`E013`, `AP7`)
   - `QuadSc_solver_hypre` may USE: `var_QuadScalar`, `pp3d_mpi`, `hypre` library
   - `QuadSc_solver_coarse` may USE: `var_QuadScalar`, `pp3d_mpi`, `UMFPACK` library
   - `QuadSc_system` may USE: All of the above
   - `def_QuadScalar` (facade) may USE: All of the above, re-exports PUBLIC interfaces

3. **FORBIDDEN dependencies** (will cause circular dependency):
   - `QuadSc_struct` CANNOT USE `QuadSc_assembly` or `QuadSc_system`
   - `QuadSc_assembly` CANNOT USE `QuadSc_system`
   - Solver modules (`hypre`, `coarse`) CANNOT USE `QuadSc_assembly` or `QuadSc_system`

4. **Verification**: After implementing split, run:
   ```bash
   # Check for circular dependencies
   grep -n "USE QuadSc" source/src_quadLS/QuadSc_*.f90 | grep -v "!.*USE"
   ```
   Manually verify the hierarchy is preserved.

## 2. Code Duplication Candidates

The following routines exhibit high duplication and should be refactored into generic wrappers.

### 2.1 Mass Matrix Assembly
*   **Current:** `Create_MRhoMat` and `Create_MMat` are nearly identical.
*   **Refactoring:** Create a generic `Assemble_Mass_Generic(MatrixArray, ElementBuilderFunc)` subroutine.
    *   Pass `BuildMRhoMat` or `BuildMMat` as an argument.
    *   Pass the target matrix array (`mg_Mmat` or `mg_MlRhomat`).

    **Implementation Location**: `QuadSc_assembly.f90`

    **Example Interface**:
    ```fortran
    ! Generic mass matrix assembly
    SUBROUTINE Assemble_Mass_Generic(mg_Mat, NLEV, ElementKernel)
      TYPE(T_Matrix), INTENT(INOUT) :: mg_Mat(:)
      INTEGER, INTENT(IN) :: NLEV
      INTERFACE
        SUBROUTINE ElementKernel(iel, local_matrix)
          INTEGER, INTENT(IN) :: iel
          REAL(KIND=dp), INTENT(OUT) :: local_matrix(:,:)
        END SUBROUTINE
      END INTERFACE

      ! Common loop structure
      DO ILEV = NLEV, NLMIN, -1
        ! Allocate if needed
        CALL SafeReallocate(mg_Mat(ILEV)%a, ...)

        ! Assembly loop
        DO iel = 1, nel
          CALL ElementKernel(iel, local_mat)
          ! Insert into global matrix
        END DO
      END DO
    END SUBROUTINE
    ```

### 2.2 Diffusion Matrix Assembly
*   **Current:** `Create_hDiffMat`, `Create_ConstDiffMat`, and `Create_DiffMat` share the same loop over levels and allocation checks. They only differ in the low-level kernel called (`DIFFQ2_alpha` vs `DIFFQ2_NEWT`).
*   **Refactoring:** Create `Assemble_Diff_Generic(MatrixArray, DiffKernelFunc, AlphaValue)`.

    **Implementation Location**: `QuadSc_assembly.f90`

### 2.3 Parallel Matrix Setup
*   **Current:** `Create_QuadLinParMatStruct` and `Fill_QuadLinParMat` share the exact same loop and call to `Create_ParB_COLMAT`.
*   **Refactoring:** Merge into a single routine `Manage_QuadLinParMat(bAllocate)` where a logical flag controls whether allocation occurs.

    **Implementation Location**: `QuadSc_struct.f90` (since it deals with structure setup)

## 3. Logic Extraction

### 3.1 UMFPACK Logic
In `Create_CMat`, lines ~1370-1430 contain explicit UMFPACK setup logic inside the matrix creation loop.
*   **Issue:** Mixing matrix assembly with solver-specific setup.
*   **Action:** Extract this block to `QuadSc_solver_coarse.f90`. This allows `Create_CMat` to focus purely on the Schur complement / Pressure matrix assembly.

    **New routine signature**:
    ```fortran
    ! In QuadSc_solver_coarse.f90
    SUBROUTINE Setup_UMFPACK_Coarse(ILEV, lMat, crsSTR, UMF_CMat)
      INTEGER, INTENT(IN) :: ILEV
      TYPE(T_Matrix), INTENT(IN) :: lMat
      TYPE(T_CSR_Structure), INTENT(INOUT) :: crsSTR
      TYPE(T_UMFPACK_Matrix), INTENT(INOUT) :: UMF_CMat

      ! Extract UMFPACK-specific setup code here
      ! Handle the /16 mystery (see def_modernization.md Section 1.3.1)
    END SUBROUTINE
    ```

### 3.2 HYPRE Setup
`SetUp_HYPRE_Solver` is a 270-line subroutine inside the matrix definition module. It depends on `myHYPRE` types which are conceptually distinct from the core FEM matrices.
*   **Action:** Move entire subroutine to `QuadSc_solver_hypre.f90`.
*   **Benefit:** Isolates HYPRE-specific code, making it easier to update or replace the solver interface.

## 4. State Management Refactoring

**Current Issue:** The module relies on global pointers defined in `var_QuadScalar`:
```fortran
qMat => mg_qMat(ILEV)  ! Implicit global state dependency
```

**Target State:** Subroutines should accept derived types as arguments:
```fortran
SUBROUTINE Create_MMat(Mesh, MMat_Array)
   TYPE(T_Mesh), INTENT(IN) :: Mesh
   TYPE(T_Matrix), INTENT(INOUT) :: MMat_Array(:)
   ...
END SUBROUTINE
```
*   **Benefit:** Allows unit testing of matrix assembly without initializing the entire global simulation state.
*   **Step 1:** Keep global pointers but move them to the *top* of the new modules to maintain current behavior.
*   **Step 2:** Gradually change interfaces to pass objects explicitly (Phase 4 in roadmap).

## 5. Unit Testing Framework (NEW)

Each new module should have corresponding unit tests.

### 5.1 Test Structure

Create `source/src_quadLS/tests/` directory with test programs:

```
source/src_quadLS/tests/
├── test_quadsc_struct.f90     # Test sparsity pattern allocation
├── test_quadsc_assembly.f90   # Test matrix assembly
├── test_quadsc_solver_hypre.f90  # Test HYPRE interface (optional)
└── CMakeLists.txt
```

### 5.2 Example: Test Matrix Structure Allocation

**File**: `source/src_quadLS/tests/test_quadsc_struct.f90`

```fortran
PROGRAM test_quadsc_struct
!===============================================================================
! DESCRIPTION:
!   Unit test for QuadSc_struct module.
!   Tests CSR structure allocation and resizing behavior.
!===============================================================================
  USE QuadSc_struct, ONLY: Create_QuadMatStruct
  USE var_QuadScalar, ONLY: mg_qMat, NLMIN, NLMAX
  IMPLICIT NONE
  INTEGER :: ILEV, expected_size, ierr
  LOGICAL :: test_passed

  WRITE(*,'(A)') '========== Testing QuadSc_struct =========='

  ! Initialize test data
  NLMIN = 1
  NLMAX = 3
  ALLOCATE(mg_qMat(NLMIN:NLMAX))

  ! Test 1: Initial allocation
  WRITE(*,'(A)') 'Test 1: Initial CSR structure allocation'
  CALL Create_QuadMatStruct(NLMAX, 100, ierr)  ! 100 DOFs
  IF (ierr /= 0) THEN
    WRITE(*,*) 'FAIL: Structure allocation failed'
    STOP 1
  END IF

  ! Verify structure was created
  DO ILEV = NLMIN, NLMAX
    IF (.NOT. ALLOCATED(mg_qMat(ILEV)%ColA)) THEN
      WRITE(*,*) 'FAIL: ColA not allocated at level', ILEV
      STOP 1
    END IF
    IF (.NOT. ALLOCATED(mg_qMat(ILEV)%LdA)) THEN
      WRITE(*,*) 'FAIL: LdA not allocated at level', ILEV
      STOP 1
    END IF
  END DO
  WRITE(*,'(A)') '  PASS'

  ! Test 2: Verify CSR structure size
  WRITE(*,'(A)') 'Test 2: Verify CSR sizes'
  ! Check that structure has reasonable size (heuristic: na > 0)
  DO ILEV = NLMIN, NLMAX
    IF (mg_qMat(ILEV)%na <= 0) THEN
      WRITE(*,*) 'FAIL: Invalid na at level', ILEV
      STOP 1
    END IF
  END DO
  WRITE(*,'(A)') '  PASS'

  ! Test 3: Reallocation with different size
  WRITE(*,'(A)') 'Test 3: Reallocation with different DOF count'
  CALL Create_QuadMatStruct(NLMAX, 200, ierr)  ! 200 DOFs
  IF (ierr /= 0) THEN
    WRITE(*,*) 'FAIL: Reallocation failed'
    STOP 1
  END IF
  WRITE(*,'(A)') '  PASS'

  ! Cleanup
  DO ILEV = NLMIN, NLMAX
    IF (ALLOCATED(mg_qMat(ILEV)%ColA)) DEALLOCATE(mg_qMat(ILEV)%ColA)
    IF (ALLOCATED(mg_qMat(ILEV)%LdA)) DEALLOCATE(mg_qMat(ILEV)%LdA)
    IF (ALLOCATED(mg_qMat(ILEV)%a)) DEALLOCATE(mg_qMat(ILEV)%a)
  END DO
  DEALLOCATE(mg_qMat)

  WRITE(*,*)
  WRITE(*,'(A)') 'All tests PASSED!'
  WRITE(*,'(A)') '=========================================='

END PROGRAM test_quadsc_struct
```

### 5.3 Integration with CMake

**File**: `source/src_quadLS/tests/CMakeLists.txt`

```cmake
# Unit tests for QuadSc modules

# Test 1: Structure allocation
add_executable(test_quadsc_struct test_quadsc_struct.f90)
target_link_libraries(test_quadsc_struct
                      src_quadLS  # Link against QuadSc modules
                      src_util)   # For memory utilities
add_test(NAME QuadSc_StructureAllocation
         COMMAND test_quadsc_struct)

# Test 2: Assembly routines
add_executable(test_quadsc_assembly test_quadsc_assembly.f90)
target_link_libraries(test_quadsc_assembly
                      src_quadLS
                      src_util
                      src_pp3d)  # For element routines
add_test(NAME QuadSc_MatrixAssembly
         COMMAND test_quadsc_assembly)

# Optional: HYPRE interface test (only if HYPRE enabled)
if(USE_HYPRE)
  add_executable(test_quadsc_solver_hypre test_quadsc_solver_hypre.f90)
  target_link_libraries(test_quadsc_solver_hypre
                        src_quadLS
                        hypre)
  add_test(NAME QuadSc_HYPREInterface
           COMMAND test_quadsc_solver_hypre)
endif()
```

**Parent CMakeLists.txt**: Add to `source/CMakeLists.txt`:
```cmake
# Enable testing
if(BUILD_TESTING)
  add_subdirectory(src_quadLS/tests)
endif()
```

### 5.4 Running Tests

```bash
# Configure with testing enabled
cmake -DBUILD_TESTING=ON ..

# Build tests
make test_quadsc_struct test_quadsc_assembly

# Run individual tests
./source/src_quadLS/tests/test_quadsc_struct

# Run all QuadSc tests
ctest -R QuadSc
```

## 6. Documentation Requirements (NEW)

**Before refactoring**, document the current behavior to preserve knowledge.

### 6.1 Current State Documentation

Create `docs/md_docs/quadsc_current_implementation.md`:

**Contents**:
1. **Matrix Storage Format**:
   - CSR (Compressed Sparse Row) structure layout
   - How `ColA`, `LdA`, and `a` arrays relate
   - Indexing conventions (0-based vs 1-based)

2. **Multigrid Level Indexing**:
   - Convention for `ILEV`, `NLMIN`, `NLMAX`
   - Direction of coarsening (fine→coarse or coarse→fine)

3. **MPI Decomposition**:
   - How matrices are distributed across processes
   - Halo/ghost region handling
   - Parallel matrix assembly synchronization

4. **Element Numbering**:
   - Local vs global element indexing
   - Mapping between elements and DOFs

### 6.2 Generate Call Graph

Use automated tools to document current dependencies:

```bash
# Using Ford (Fortran documentation generator)
ford quadsc_doc_project.md

# Or using doxygen (if configured for Fortran)
doxygen Doxyfile

# Or manual grep-based analysis
grep -r "CALL Create_" source/src_quadLS/QuadSc_def.f90 > call_graph.txt
```

### 6.3 Inline Documentation

Add documentation headers to each major routine BEFORE moving:

```fortran
!===============================================================================
SUBROUTINE Create_MMat(...)
!-------------------------------------------------------------------------------
! DESCRIPTION:
!   Assembles the mass matrix for Q2 velocity elements across all multigrid
!   levels. Uses E013 element integration routine.
!
! ALGORITHM:
!   1. Loop over multigrid levels (fine to coarse)
!   2. For each level:
!      a. Allocate CSR structure if needed
!      b. Loop over elements, assemble local mass matrix
!      c. Insert into global CSR matrix
!
! INPUTS:
!   - Uses global mesh data from mg_mesh(ILEV)
!   - Uses global Q2 structures from mg_qMat(ILEV)
!
! OUTPUTS:
!   - Fills mg_Mmat(ILEV)%a for each level
!
! DEPENDENCIES:
!   - E013 (element mass matrix routine)
!   - var_QuadScalar (global state)
!
! PERFORMANCE:
!   - Typically 10-20% of total setup time
!   - Parallel via MPI domain decomposition
!
! NOTES:
!   - Matrix is symmetric
!   - Uses CSR storage format
!-------------------------------------------------------------------------------
  ...
END SUBROUTINE
```

### 6.4 Documentation Checklist

Before splitting modules, ensure:
- [ ] All major subroutines have documentation headers
- [ ] Global variables are documented with purpose and lifecycle
- [ ] CSR storage format is documented
- [ ] Multigrid indexing conventions are documented
- [ ] MPI parallelization strategy is documented
- [ ] Call graph has been generated and reviewed

## 7. Verification After Refactoring

After each module split:

1. **Compilation**: Ensure clean build with no warnings
   ```bash
   make clean && make -j8 2>&1 | tee build.log
   grep -i "warning\|error" build.log
   ```

2. **Unit Tests**: All new unit tests must pass
   ```bash
   ctest -R QuadSc --output-on-failure
   ```

3. **Integration Tests**: Full regression suite
   ```bash
   ctest --output-on-failure
   ```

4. **Benchmark Comparison**: See `refactoring_roadmap.md` Phase 0 for baseline comparison strategy

## 8. Notes

*   **Module Granularity**: The proposed 5-module split balances maintainability with complexity. If a module grows beyond ~1000 lines, consider further splitting.

*   **Backward Compatibility**: The facade pattern (Section 1.1) ensures that existing code using `USE def_QuadScalar` continues to work without modification.

*   **Performance**: Module boundaries should not impact runtime performance when compiler optimization is enabled (`-O3`). Modern Fortran compilers can inline across module boundaries.

*   **Testing Strategy**: Unit tests are essential for catching regressions early. Integration tests verify overall system behavior.
