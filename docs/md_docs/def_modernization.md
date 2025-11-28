# Modernization Plan: QuadSc_def.f90

**Date:** November 28, 2025  
**Status:** Planning / Analysis Phase

This document tracks the analysis, pitfalls, and modernization strategy for the `QuadSc_def.f90` module (specifically `MODULE def_QuadScalar`). This module is critical as it handles the allocation and assembly of the system matrices (Mass, Stiffness, Gradient, Divergence) for the FEM solver.

## 1. Robustness & Implementation Pitfalls (Current State)

Analysis reveals several high-risk patterns that compromise stability, particularly for dynamic mesh scenarios or varying problem sizes.

### 1.1 Dangerous Memory Re-allocation Logic (CRITICAL)
The code protects against double-allocation but **fails to handle resizing**.
*   **Pattern:**
    ```fortran
    IF (.not.ALLOCATED(mg_hDmat(ILEV)%a)) THEN
       ALLOCATE(mg_hDmat(ILEV)%a(qMat%na))
    END IF
    ```
*   **Risk:** If `qMat%na` (number of non-zeros) changes—e.g., during mesh refinement (AMR) or moving mesh simulations—the code will **silently use the old, wrong-sized array** or crash if bounds are exceeded.
*   **Remediation:** Implement a safe reallocation utility.
    ```fortran
    IF (ALLOCATED(arr)) THEN
       IF (SIZE(arr) /= new_size) DEALLOCATE(arr)
    END IF
    IF (.NOT. ALLOCATED(arr)) ALLOCATE(arr(new_size), STAT=ierr)
    ```

### 1.2 "Magic Number" Pre-allocations
*   **Pattern:** `MatSize = 300*NDOF` (in `Create_QuadMatStruct`) and `MatSize = 16*27*...` (in `Create_QuadLinMatStruct`).
*   **Risk:** Relies on heuristic multipliers to estimate sparse matrix storage. If the mesh topology is irregular or element order increases, this will cause buffer overflows.
*   **Remediation:**
    1.  Short term: Ensure these constants are named parameters (e.g., `MAX_NNZ_PER_ROW`).
    2.  Long term: Implement a "count-first" pass or dynamic structure (e.g., linked list) to calculate exact NNZ before final CSR allocation.

### 1.3 The "Mystery Divider" (UMFPACK Interface)
*   **Pattern:** `crsSTR%A%na = lMat%na/16 !!!! /16????????????????` in `Create_CMat`.
*   **Risk:** Indicates "programming by coincidence." This likely hardcodes a relationship between fine grid and coarse grid sparsity or Q2 vs P1 connectivity, which is fragile.
*   **Remediation:** Reverse engineer the exact mapping between `lMat` (fine/parallel matrix) and the UMFPACK serial structure to replace the magic number with a calculated dimension.

### 1.4 Fragile MPI Logic
*   **Pattern:** Hardcoded rank checks like `IF (myid.eq.1)` for debug prints or specific logic.
*   **Risk:** Fails if the MPI topology changes or if rank 1 is not special.
*   **Remediation:** Use named logicals (e.g., `bIsMaster`, `bIsDebugRank`) or standard root rank checks (`IF (myid == root_rank)`).

### 1.5 Unchecked Allocations
*   **Pattern:** `ALLOCATE(...)` without `STAT=ierr`.
*   **Risk:** Immediate segfault/crash on Out-Of-Memory (OOM) instead of a graceful error message.

## 2. Modernization Roadmap

This effort will be split into phases to avoid destabilizing the solver.

### Phase 1: Safety & Types (Low Risk)
- [ ] **Type Modernization**: Replace `REAL*8` with `REAL(KIND=dp)` and `CHARACTER*10` with `CHARACTER(LEN=10)`.
- [ ] **Implicit None**: Ensure strict typing (already present, but verify scope).
- [ ] **Allocation Safety**: wrapper routine for allocation that checks `STAT` and handles resizing logic correctly. Apply this to all matrix allocations.

### Phase 2: Architecture & State (Medium Risk)
- [ ] **Remove Global Pointer Dependency**: The module currently modifies global pointers (e.g., `qMat => ...`) defined in `var_QuadScalar`.
    - *Goal*: Refactor subroutines to accept matrices as `INTENT(INOUT)` arguments.
    - *Benefit*: Makes subroutines thread-safe and testable.

### Phase 3: Logic Repair (High Risk)
- [ ] **Fix UMFPACK Setup**: Investigate and document the `/16` magic number. Replace with deterministic logic.
- [ ] **Dynamic Sparsity**: Replace `300*NDOF` heuristics with a pre-computation loop to determine exact matrix sizes.

## 3. Notes for Developers

*   **Testing:** Any changes to this file require running the full regression suite (`ctest`) because these matrices underpin the entire flow solver.
*   **Dependencies:** This module depends heavily on `pp3d_mpi` and low-level FEM routines (`AP7`, `AP9`, `E013`). Changes here may require interface updates in those assembly routines.
