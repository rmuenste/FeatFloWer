# Refactoring Roadmap: QuadSc_def.f90

**Strategy:** Phased Refactoring (Lowest Risk $\to$ Highest Risk)

This document outlines the execution plan for refactoring `source/src_quadLS/QuadSc_def.f90`. The guiding principle is to maintain a working build at every step, minimizing "big bang" changes.

## Phase 0: Preparation & Baseline

Before touching code, establish a safety net.

1.  **Baseline Generation**: Run a standard benchmark (e.g., `q2p1_devel` or `q2p1_bench_sedimentation`) and save the `pe1.log` and any `.gmv`/`.vtu` output. This is our "Golden Output."
2.  **Test Check**: Ensure `ctest` passes 100%.

## Phase 1: Logic Extraction (Low Risk)

**Goal**: Move isolated, heavy logic into separate modules *without* changing the dependency graph of the main application. `QuadSc_def.f90` will essentially become thinner by offloading work to new modules, but it will still be the primary entry point for consumers.

*   **Step 1.1: Extract HYPRE Interface**
    *   **Action**: Move `SetUp_HYPRE_Solver` (and its helper `SORT_DOFs`) to a new module `QuadSc_solver_hypre.f90`.
    *   **Integration**: `def_QuadScalar` will `USE mod_QuadSc_solver_hypre`.
    *   **Risk**: Low. This routine is self-contained.

*   **Step 1.2: Extract Coarse Grid Solver**
    *   **Action**: Identify the UMFPACK setup code block inside `Create_CMat`. Move it to `QuadSc_solver_coarse.f90`.
    *   **Integration**: `Create_CMat` calls the new external routine.
    *   **Risk**: Low/Medium. Requires careful handling of local variables (`crsSTR`, `UMF_CMat`) passed to the new routine.

## Phase 2: Internal Deduplication (Medium Risk)

**Goal**: Reduce code size by identifying common patterns. We will implement "Generic Assembly Routines" **inside** the current file first (or a helper module) to prove logic correctness before splitting files.

*   **Step 2.1: Mass Matrix Unification**
    *   **Action**: Create `Assemble_Mass_Generic`. Refactor `Create_MMat` and `Create_MRhoMat` to call this generic routine.
    *   **Risk**: Medium. Must ensure `E013` element integration behaves identically when passed as an argument vs hardcoded.

*   **Step 2.2: Diffusion Matrix Unification**
    *   **Action**: Create `Assemble_Diff_Generic`. Refactor `Create_hDiffMat` and `Create_ConstDiffMat`.
    *   **Risk**: Low. These routines are structurally identical.

*   **Step 2.3: Parallel Matrix Setup**
    *   **Action**: Merge logic from `Create_QuadLinParMatStruct` and `Fill_QuadLinParMat`.
    *   **Risk**: Medium. Logic flow for allocation vs filling needs to be handled by a clean flag or argument.

## Phase 3: Structural Split (High Risk)

**Goal**: Physically break the file into multiple compilation units. This changes `CMakeLists.txt` and potentially `USE` statements in other files.

*   **Step 3.1: Create `QuadSc_struct.f90`**
    *   **Action**: Move `Create_QuadMatStruct`, `Create_QuadLinMatStruct`, `Create_LinMatStruct` to a new module.
    *   **Dependency**: This module will likely depend on `pp3d_mpi` and `var_QuadScalar`.

*   **Step 3.2: Create `QuadSc_assembly.f90`**
    *   **Action**: Move the newly deduplicated assembly routines (`Create_MMat`, `Create_KMat`, etc.) to this module.
    *   **Dependency**: Will depend on `QuadSc_struct` (conceptually) or just share `var_QuadScalar`.

*   **Step 3.3: Facade Pattern (The Switch)**
    *   **Action**: Turn `def_QuadScalar` into a "Facade Module" that simply `USE`s and `PUBLIC`s the new modules.
    *   **Benefit**: This prevents breaking the hundreds of calls in `QuadSc_main.f90` that expect `USE def_QuadScalar`.
    *   **Code Example**:
        ```fortran
        MODULE def_QuadScalar
          USE mod_QuadSc_struct, ONLY: Create_QuadMatStruct, ...
          USE mod_QuadSc_assembly, ONLY: Create_MMat, ...
          USE mod_QuadSc_solver_hypre
          IMPLICIT NONE
          PUBLIC  ! Re-export everything
        END MODULE
        ```

## Phase 4: API Modernization (Highest Risk)

**Goal**: Remove global pointer dependency (`qMat => ...`). This is technically "Modernization" but fits the refactoring arc.

*   **Step 4.1**: Update interfaces in `QuadSc_assembly.f90` to accept `INTENT(INOUT)` matrix types.
*   **Step 4.2**: Update the "Facade" module to pass the globals into these clean interfaces.

---

**Verification Strategy:**
After *every* step in Phase 1 and 2, and after the major split in Phase 3:
1.  Compile (`make -j`).
2.  Run `q2p1_devel` (benchmark case).
3.  Diff the output logs against the "Golden Output" from Phase 0.
4.  Run `ctest`.
