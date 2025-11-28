# Refactoring Plan: QuadSc_def.f90

**Date:** November 28, 2025
**Status:** Proposal

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

## 2. Code Duplication Candidates

The following routines exhibit high duplication and should be refactored into generic wrappers.

### 2.1 Mass Matrix Assembly
*   **Current:** `Create_MRhoMat` and `Create_MMat` are nearly identical.
*   **Refactoring:** Create a generic `Assemble_Mass_Generic(MatrixArray, ElementBuilderFunc)` subroutine.
    *   Pass `BuildMRhoMat` or `BuildMMat` as an argument.
    *   Pass the target matrix array (`mg_Mmat` or `mg_MlRhomat`).

### 2.2 Diffusion Matrix Assembly
*   **Current:** `Create_hDiffMat`, `Create_ConstDiffMat`, and `Create_DiffMat` share the same loop over levels and allocation checks. They only differ in the low-level kernel called (`DIFFQ2_alpha` vs `DIFFQ2_NEWT`).
*   **Refactoring:** Create `Assemble_Diff_Generic(MatrixArray, DiffKernelFunc, AlphaValue)`.

### 2.3 Parallel Matrix Setup
*   **Current:** `Create_QuadLinParMatStruct` and `Fill_QuadLinParMat` share the exact same loop and call to `Create_ParB_COLMAT`.
*   **Refactoring:** Merge into a single routine `Manage_QuadLinParMat(bAllocate)` where a logical flag controls whether allocation occurs.

## 3. Logic Extraction

### 3.1 UMFPACK Logic
In `Create_CMat`, lines ~1370-1430 contain explicit UMFPACK setup logic inside the matrix creation loop.
*   **Issue:** Mixing matrix assembly with solver-specific setup.
*   **Action:** Extract this block to `QuadSc_solver_coarse.f90`. This allows `Create_CMat` to focus purely on the Schur complement / Pressure matrix assembly.

### 3.2 HYPRE Setup
`SetUp_HYPRE_Solver` is a 270-line subroutine inside the matrix definition module. It depends on `myHYPRE` types which are conceptually distinct from the core FEM matrices.
*   **Action:** Move entire subroutine to `QuadSc_solver_hypre.f90`.

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
*   **Step 2:** Gradually change interfaces to pass objects explicitly.
