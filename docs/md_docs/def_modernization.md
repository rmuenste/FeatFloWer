# Modernization Plan: QuadSc_def.f90

**Date:** November 28, 2025
**Status:** Planning / Analysis Phase
**Last Updated:** November 28, 2025 (Enhanced with safety utilities and performance monitoring)

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

*   **Remediation:** Implement a robust safe reallocation utility with full error checking.

    **Implementation:** Create `source/src_util/QuadSc_memory_utils.f90`:
    ```fortran
    MODULE QuadSc_memory_utils
    !===============================================================================
    ! MODULE: QuadSc_memory_utils
    !
    ! DESCRIPTION:
    !   Provides safe memory allocation/deallocation utilities with comprehensive
    !   error checking and automatic resizing capabilities.
    !===============================================================================
      IMPLICIT NONE
      PRIVATE

      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)

      ! Generic interface for different array types
      PUBLIC :: SafeReallocate
      INTERFACE SafeReallocate
        MODULE PROCEDURE SafeReallocate_Real
        MODULE PROCEDURE SafeReallocate_Integer
        MODULE PROCEDURE SafeReallocate_Logical
      END INTERFACE

    CONTAINS

    !===============================================================================
    SUBROUTINE SafeReallocate_Real(arr, new_size, routine_name, ierr)
    !-------------------------------------------------------------------------------
    ! DESCRIPTION:
    !   Safely allocate or reallocate a REAL array with full error checking.
    !   - Checks if reallocation is actually needed (size changed)
    !   - Handles deallocation errors
    !   - Handles allocation errors with informative messages
    !   - Initializes new arrays to zero
    !
    ! ARGUMENTS:
    !   arr          - Array to allocate/reallocate (INOUT)
    !   new_size     - Desired size (IN)
    !   routine_name - Name of calling routine for error messages (IN)
    !   ierr         - Error code: 0=success, >0=error (OUT)
    !-------------------------------------------------------------------------------
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT) :: arr(:)
      INTEGER, INTENT(IN) :: new_size
      CHARACTER(len=*), INTENT(IN) :: routine_name
      INTEGER, INTENT(OUT) :: ierr

      ierr = 0

      ! Check if array needs resizing
      IF (ALLOCATED(arr)) THEN
        IF (SIZE(arr) /= new_size) THEN
          DEALLOCATE(arr, STAT=ierr)
          IF (ierr /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'ERROR: DEALLOCATE failed in ', &
                                   TRIM(routine_name), ', STAT=', ierr
            RETURN
          END IF
        ELSE
          ! Already correct size, nothing to do
          RETURN
        END IF
      END IF

      ! Allocate if not already allocated
      IF (.NOT. ALLOCATED(arr)) THEN
        ALLOCATE(arr(new_size), STAT=ierr)
        IF (ierr /= 0) THEN
          WRITE(*,'(A,I0,A,A,A,I0)') 'ERROR: Cannot allocate ', new_size, &
                                      ' REAL elements in ', TRIM(routine_name), &
                                      ', STAT=', ierr
          RETURN
        END IF
        arr = 0.0_dp  ! Initialize to zero
      END IF

    END SUBROUTINE SafeReallocate_Real


    !===============================================================================
    SUBROUTINE SafeReallocate_Integer(arr, new_size, routine_name, ierr)
    !-------------------------------------------------------------------------------
    ! DESCRIPTION: Safe reallocation for INTEGER arrays
    !-------------------------------------------------------------------------------
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: arr(:)
      INTEGER, INTENT(IN) :: new_size
      CHARACTER(len=*), INTENT(IN) :: routine_name
      INTEGER, INTENT(OUT) :: ierr

      ierr = 0

      IF (ALLOCATED(arr)) THEN
        IF (SIZE(arr) /= new_size) THEN
          DEALLOCATE(arr, STAT=ierr)
          IF (ierr /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'ERROR: DEALLOCATE failed in ', &
                                   TRIM(routine_name), ', STAT=', ierr
            RETURN
          END IF
        ELSE
          RETURN
        END IF
      END IF

      IF (.NOT. ALLOCATED(arr)) THEN
        ALLOCATE(arr(new_size), STAT=ierr)
        IF (ierr /= 0) THEN
          WRITE(*,'(A,I0,A,A,A,I0)') 'ERROR: Cannot allocate ', new_size, &
                                      ' INTEGER elements in ', TRIM(routine_name), &
                                      ', STAT=', ierr
          RETURN
        END IF
        arr = 0  ! Initialize to zero
      END IF

    END SUBROUTINE SafeReallocate_Integer


    !===============================================================================
    SUBROUTINE SafeReallocate_Logical(arr, new_size, routine_name, ierr)
    !-------------------------------------------------------------------------------
    ! DESCRIPTION: Safe reallocation for LOGICAL arrays
    !-------------------------------------------------------------------------------
      LOGICAL, ALLOCATABLE, INTENT(INOUT) :: arr(:)
      INTEGER, INTENT(IN) :: new_size
      CHARACTER(len=*), INTENT(IN) :: routine_name
      INTEGER, INTENT(OUT) :: ierr

      ierr = 0

      IF (ALLOCATED(arr)) THEN
        IF (SIZE(arr) /= new_size) THEN
          DEALLOCATE(arr, STAT=ierr)
          IF (ierr /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'ERROR: DEALLOCATE failed in ', &
                                   TRIM(routine_name), ', STAT=', ierr
            RETURN
          END IF
        ELSE
          RETURN
        END IF
      END IF

      IF (.NOT. ALLOCATED(arr)) THEN
        ALLOCATE(arr(new_size), STAT=ierr)
        IF (ierr /= 0) THEN
          WRITE(*,'(A,I0,A,A,A,I0)') 'ERROR: Cannot allocate ', new_size, &
                                      ' LOGICAL elements in ', TRIM(routine_name), &
                                      ', STAT=', ierr
          RETURN
        END IF
        arr = .FALSE.  ! Initialize to false
      END IF

    END SUBROUTINE SafeReallocate_Logical

    END MODULE QuadSc_memory_utils
    ```

    **Usage Example:**
    ```fortran
    USE QuadSc_memory_utils, ONLY: SafeReallocate

    ! Old dangerous pattern:
    ! IF (.not.ALLOCATED(mg_hDmat(ILEV)%a)) ALLOCATE(mg_hDmat(ILEV)%a(qMat%na))

    ! New safe pattern:
    CALL SafeReallocate(mg_hDmat(ILEV)%a, qMat%na, 'Create_hDiffMat', ierr)
    IF (ierr /= 0) THEN
      ! Handle error appropriately
      RETURN
    END IF
    ```

### 1.2 "Magic Number" Pre-allocations
*   **Pattern:** `MatSize = 300*NDOF` (in `Create_QuadMatStruct`) and `MatSize = 16*27*...` (in `Create_QuadLinMatStruct`).
*   **Risk:** Relies on heuristic multipliers to estimate sparse matrix storage. If the mesh topology is irregular or element order increases, this will cause buffer overflows.
*   **Remediation:**
    1.  **Short term**: Ensure these constants are named parameters (e.g., `MAX_NNZ_PER_ROW`).
        ```fortran
        INTEGER, PARAMETER :: MAX_NNZ_PER_DOF_Q2 = 300  ! Heuristic for Q2 elements
        INTEGER, PARAMETER :: MAX_NNZ_PER_DOF_Q2P1 = 16 * 27  ! Q2-P1 coupling

        MatSize = MAX_NNZ_PER_DOF_Q2 * NDOF
        ```
    2.  **Long term**: Implement a "count-first" pass or dynamic structure (e.g., linked list) to calculate exact NNZ before final CSR allocation.
        - See **Performance Monitoring** section below for impact assessment.

### 1.3 The "Mystery Divider" (UMFPACK Interface)
*   **Pattern:** `crsSTR%A%na = lMat%na/16 !!!! /16????????????????` in `Create_CMat`.
*   **Risk:** Indicates "programming by coincidence." This likely hardcodes a relationship between fine grid and coarse grid sparsity or Q2 vs P1 connectivity, which is fragile.

#### 1.3.1 Reverse Engineering Strategy for /16 Mystery

**Context:** The `/16` divider appears in two critical locations:
1. **HYPRE setup** (line 689): `myHYPRE%nonzeros = lPMat%na/16 + lMat%na/16`
2. **UMFPACK setup** (line 1142): `crsSTR%A%na = lMat%na/16`

Also note `/4` for unknowns: `crsSTR%A%nu = lMat%nu/4`

**Hypothesis**: This is related to Q2 vs P1 element structure:
- Q2 (quadratic) elements: 27 nodes in 3D, 9 nodes in 2D
- P1 (linear) elements: 8 nodes in 3D, 4 nodes in 2D
- DOF ratio in 3D: 27/8 ≈ 3.375; velocity vs pressure: 3×27 / 8 = 10.125
- Sparsity pattern ratio may involve (3×3)/(something) = 9/something or 16/something

**Investigation Steps:**

**Step 1**: Add diagnostic logging BEFORE refactoring. Create a diagnostic version:

```fortran
! In Create_CMat, around line 1140
IF (myid == 0 .and. ILEV == NLMIN) THEN
  WRITE(*,*)
  WRITE(*,'(A)') '========== UMFPACK Matrix Size Diagnostic =========='
  WRITE(*,'(A,I0)') 'Multigrid level (NLMIN): ', NLMIN
  WRITE(*,'(A,I0)') 'Fine grid unknowns (lMat%nu):     ', lMat%nu
  WRITE(*,'(A,I0)') 'Fine grid non-zeros (lMat%na):    ', lMat%na
  WRITE(*,'(A,I0)') 'Coarse unknowns (lMat%nu/4):      ', lMat%nu/4
  WRITE(*,'(A,I0)') 'Coarse non-zeros (lMat%na/16):    ', lMat%na/16
  WRITE(*,'(A,F10.4)') 'Actual nu ratio:                  ', REAL(lMat%nu)/REAL(lMat%nu/4)
  WRITE(*,'(A,F10.4)') 'Actual na ratio:                  ', REAL(lMat%na)/REAL(lMat%na/16)
  WRITE(*,'(A,I0)') 'Average NNZ per row (fine):       ', lMat%na / lMat%nu
  WRITE(*,'(A,I0)') 'Average NNZ per row (coarse):     ', (lMat%na/16) / (lMat%nu/4)
  WRITE(*,'(A)') '===================================================='
  WRITE(*,*)
END IF
```

**Step 2**: Run diagnostic on multiple test cases:
- 2D uniform mesh
- 3D uniform mesh
- 2D refined mesh
- 3D refined mesh
- Irregular mesh

**Step 3**: Document findings in table:

| Mesh Type | Dimension | lMat%nu | lMat%na | nu/4 | na/16 | Actual Ratio |
|-----------|-----------|---------|---------|------|-------|--------------|
| Uniform   | 2D        | ???     | ???     | ???  | ???   | ???          |
| Uniform   | 3D        | ???     | ???     | ???  | ???   | ???          |

**Step 4**: Based on findings, replace magic number with calculated value:
```fortran
! Instead of:
! crsSTR%A%na = lMat%na/16

! Use documented formula (example):
! For Q2-P1 3D: Coarse grid (P1 pressure) has 1/4 the unknowns
! and approximately 1/16 the sparsity due to reduced connectivity
INTEGER, PARAMETER :: Q2_TO_P1_NU_RATIO = 4  ! Velocity vs Pressure DOFs
INTEGER, PARAMETER :: Q2_TO_P1_NA_RATIO = 16 ! Sparsity pattern ratio

crsSTR%A%nu = lMat%nu / Q2_TO_P1_NU_RATIO
crsSTR%A%na = lMat%na / Q2_TO_P1_NA_RATIO
```

**Action Item**: Run diagnostics in Phase 0 before any refactoring.

### 1.4 Fragile MPI Logic
*   **Pattern:** Hardcoded rank checks like `IF (myid.eq.1)` for debug prints or specific logic.
*   **Risk:** Fails if the MPI topology changes or if rank 1 is not special.
*   **Remediation:** Use named logicals (e.g., `bIsMaster`, `bIsDebugRank`) or standard root rank checks (`IF (myid == root_rank)`).

    **Implementation:**
    ```fortran
    ! In var_QuadScalar or PP3D_MPI module
    INTEGER, PARAMETER :: ROOT_RANK = 0
    INTEGER, PARAMETER :: DEBUG_RANK = 0  ! Can be changed for debugging
    LOGICAL :: bIsMaster, bIsDebugRank

    ! Initialize once
    bIsMaster = (myid == ROOT_RANK)
    bIsDebugRank = (myid == DEBUG_RANK)

    ! Use throughout code
    IF (bIsMaster) WRITE(*,*) "Master process: ..."
    IF (bIsDebugRank) WRITE(*,*) "Debug info: ..."
    ```

### 1.5 Unchecked Allocations
*   **Pattern:** `ALLOCATE(...)` without `STAT=ierr`.
*   **Risk:** Immediate segfault/crash on Out-Of-Memory (OOM) instead of a graceful error message.
*   **Remediation:** Use `SafeReallocate` utility (Section 1.1) or at minimum add STAT checking to all allocations.

### 1.6 Const-correctness (NEW)
*   **Pattern:** Many subroutines don't specify `INTENT` for parameters.
*   **Risk:** Accidental modification of read-only data, unclear interfaces.
*   **Remediation:** Add `INTENT(IN)`, `INTENT(OUT)`, `INTENT(INOUT)` to all subroutine parameters.
    ```fortran
    ! Before:
    SUBROUTINE Create_MMat(NLEV, NDOF)

    ! After:
    SUBROUTINE Create_MMat(NLEV, NDOF)
      INTEGER, INTENT(IN) :: NLEV, NDOF
    ```

## 2. Modernization Roadmap

This effort will be split into phases to avoid destabilizing the solver.

### Phase 1: Safety & Types (Low Risk)

#### Phase 1.1: Type Modernization
- [ ] Replace `REAL*8` with `REAL(KIND=dp)` where `dp = SELECTED_REAL_KIND(15, 307)`
- [ ] Replace `CHARACTER*10` with `CHARACTER(LEN=10)`
- [ ] Replace old-style operators: `.EQ.` → `==`, `.NE.` → `/=`, `.LT.` → `<`, `.GT.` → `>`

#### Phase 1.2: Implicit None Verification
- [ ] Ensure `IMPLICIT NONE` is present in module and all subroutines
- [ ] Verify no implicit variables exist (compile with `-fimplicit-none` or `-u`)

#### Phase 1.3: Allocation Safety
- [ ] **Create `QuadSc_memory_utils.f90` module** with `SafeReallocate` generic interface
- [ ] **Add unit tests** for memory utilities (see example in Phase 1.5)
- [ ] Apply `SafeReallocate` to all matrix allocations:
  - `mg_hDmat(ILEV)%a`
  - `mg_Mmat(ILEV)%a`
  - `mg_qMat(ILEV)%a`
  - `mg_KMat(ILEV)%a`
  - All CSR structure arrays (`ColA`, `LdA`)

#### Phase 1.4: Add INTENT Declarations
- [ ] Add `INTENT(IN)`, `INTENT(OUT)`, `INTENT(INOUT)` to all subroutine parameters
- [ ] Document which parameters are modified vs read-only

#### Phase 1.5: Unit Test Framework (NEW)

Create `source/src_quadLS/tests/test_memory_utils.f90`:

```fortran
PROGRAM test_safe_reallocation
  USE QuadSc_memory_utils, ONLY: SafeReallocate
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
  REAL(KIND=dp), ALLOCATABLE :: arr_real(:)
  INTEGER, ALLOCATABLE :: arr_int(:)
  INTEGER :: ierr, i

  WRITE(*,'(A)') '========== Testing SafeReallocate =========='

  ! Test 1: Initial allocation
  WRITE(*,'(A)') 'Test 1: Initial allocation'
  CALL SafeReallocate(arr_real, 100, 'Test1', ierr)
  IF (ierr /= 0 .OR. .NOT. ALLOCATED(arr_real)) THEN
    WRITE(*,*) 'FAIL: Initial allocation'
    STOP 1
  END IF
  WRITE(*,'(A)') '  PASS'

  ! Test 2: Same size (should not reallocate)
  WRITE(*,'(A)') 'Test 2: Same size - no reallocation'
  arr_real = 42.0_dp
  CALL SafeReallocate(arr_real, 100, 'Test2', ierr)
  IF (arr_real(1) /= 42.0_dp) THEN
    WRITE(*,*) 'FAIL: Unnecessary reallocation occurred'
    STOP 1
  END IF
  WRITE(*,'(A)') '  PASS'

  ! Test 3: Resize larger
  WRITE(*,'(A)') 'Test 3: Resize larger'
  CALL SafeReallocate(arr_real, 200, 'Test3', ierr)
  IF (SIZE(arr_real) /= 200) THEN
    WRITE(*,*) 'FAIL: Resize larger'
    STOP 1
  END IF
  WRITE(*,'(A)') '  PASS'

  ! Test 4: Resize smaller
  WRITE(*,'(A)') 'Test 4: Resize smaller'
  CALL SafeReallocate(arr_real, 50, 'Test4', ierr)
  IF (SIZE(arr_real) /= 50) THEN
    WRITE(*,*) 'FAIL: Resize smaller'
    STOP 1
  END IF
  WRITE(*,'(A)') '  PASS'

  ! Test 5: Integer array
  WRITE(*,'(A)') 'Test 5: Integer array allocation'
  CALL SafeReallocate(arr_int, 100, 'Test5', ierr)
  IF (ierr /= 0 .OR. SIZE(arr_int) /= 100) THEN
    WRITE(*,*) 'FAIL: Integer array allocation'
    STOP 1
  END IF
  WRITE(*,'(A)') '  PASS'

  WRITE(*,*)
  WRITE(*,'(A)') 'All tests PASSED!'
  WRITE(*,'(A)') '=========================================='

  DEALLOCATE(arr_real, arr_int)
END PROGRAM test_safe_reallocation
```

**Integration**: Add to `CMakeLists.txt`:
```cmake
# Unit tests for QuadSc utilities
add_executable(test_memory_utils source/src_quadLS/tests/test_memory_utils.f90)
target_link_libraries(test_memory_utils src_util)
add_test(NAME QuadSc_MemoryUtils COMMAND test_memory_utils)
```

### Phase 2: Architecture & State (Medium Risk)
- [ ] **Remove Global Pointer Dependency**: The module currently modifies global pointers (e.g., `qMat => ...`) defined in `var_QuadScalar`.
    - *Goal*: Refactor subroutines to accept matrices as `INTENT(INOUT)` arguments.
    - *Benefit*: Makes subroutines thread-safe and testable.
    - *Implementation*: See `def_refactoring.md` Section 4 for details.

### Phase 3: Logic Repair (High Risk)
- [ ] **Fix UMFPACK Setup**: Investigate and document the `/16` magic number (see Section 1.3.1). Replace with deterministic logic.
- [ ] **Dynamic Sparsity**: Replace `300*NDOF` heuristics with a pre-computation loop to determine exact matrix sizes.
    - **Caution**: Monitor performance impact (see Section 4 below).

## 3. Performance Monitoring

Any changes to allocation logic or assembly routines must be monitored for performance impact.

### 3.1 Baseline Timing

Add to Phase 0 baseline generation:

```fortran
! In Create_QuadMatStruct
REAL(KIND=dp) :: t_start, t_end

CALL CPU_TIME(t_start)
! ... matrix structure allocation code ...
CALL CPU_TIME(t_end)

IF (myid == 0) THEN
  WRITE(*,'(A,F8.4,A)') 'Create_QuadMatStruct time: ', t_end - t_start, ' sec'
END IF
```

Repeat for all major allocation/assembly routines:
- `Create_MMat`
- `Create_DiffMat`
- `Create_KMat`
- `Create_SMat`
- `Create_BMat`

### 3.2 Generic Assembly Performance Test

When implementing generic assembly routines (Phase 2, see `def_refactoring.md`):

1. **Keep both versions temporarily**:
   ```fortran
   ! Old hardcoded version
   SUBROUTINE Create_MMat_Old(...)

   ! New generic version
   SUBROUTINE Create_MMat_Generic(...)
   ```

2. **Add timing comparison**:
   ```fortran
   CALL CPU_TIME(t1)
   CALL Create_MMat_Old(...)
   CALL CPU_TIME(t2)

   CALL CPU_TIME(t3)
   CALL Create_MMat_Generic(...)
   CALL CPU_TIME(t4)

   IF (myid == 0) THEN
     WRITE(*,'(A,F8.4,A)') 'Old version: ', t2-t1, ' sec'
     WRITE(*,'(A,F8.4,A)') 'New version: ', t4-t3, ' sec'
     WRITE(*,'(A,F8.2,A)') 'Overhead:    ', 100.0*(t4-t3)/(t2-t1) - 100.0, '%'
   END IF
   ```

3. **Performance criteria**:
   - If new version is <5% slower → acceptable
   - If 5-10% slower → investigate compiler optimizations (`-O3`, `-ffast-math`)
   - If >10% slower → consider alternative approaches:
     - ELEMENTAL functions
     - Preprocessor macros instead of procedure pointers
     - OpenMP parallelization (see Section 5)

### 3.3 Count-First Dynamic Sparsity

**Before implementing** (Phase 3):
- Profile current allocation time vs total assembly time
- If allocation is <1% of total time, exact counting may not be worth the added complexity
- If allocation is >5% of total time OR if buffer overflows occur, implement count-first

## 4. Additional Modernization Opportunities

### 4.1 OpenMP Parallelization (NEW)

Matrix assembly loops are embarrassingly parallel. After Phase 2 (removal of global pointers), consider:

```fortran
! In Create_MMat or Create_KMat
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel, kve, ia, ib, ...)
DO iel = 1, nel
  ! Element assembly
  CALL E013(...)  ! Mass matrix element

  ! Assembly into global matrix (atomic add)
  !$OMP CRITICAL
  DO ia = 1, nvt
    DO ib = 1, nvt
      qMat%a(index) = qMat%a(index) + local_matrix(ia, ib)
    END DO
  END DO
  !$OMP END CRITICAL
END DO
!$OMP END PARALLEL DO
```

**Prerequisites**:
- Ensure element routines (`E013`, `AP7`, `AP9`) are thread-safe
- Use `!$OMP CRITICAL` or `!$OMP ATOMIC` for global matrix updates
- Test with `OMP_NUM_THREADS=1,2,4,8` for scaling

### 4.2 Centralized Error Handling Module (NEW)

Create `source/src_util/QuadSc_errors.f90`:

```fortran
MODULE QuadSc_errors
!===============================================================================
! MODULE: QuadSc_errors
!
! DESCRIPTION:
!   Centralized error handling for QuadSc modules.
!   Provides consistent error logging and optional program termination.
!===============================================================================
  USE PP3D_MPI, ONLY: myid, master
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: HandleError, ErrorCode

  ! Error codes
  TYPE :: ErrorCode
    INTEGER :: ALLOC_FAILED = 1
    INTEGER :: SIZE_MISMATCH = 2
    INTEGER :: INVALID_LEVEL = 3
    INTEGER :: MPI_ERROR = 4
    INTEGER :: SOLVER_FAILED = 5
  END TYPE

  TYPE(ErrorCode), PARAMETER :: ERR = ErrorCode()

CONTAINS

  SUBROUTINE HandleError(code, routine, message, should_abort)
    INTEGER, INTENT(IN) :: code
    CHARACTER(len=*), INTENT(IN) :: routine, message
    LOGICAL, INTENT(IN), OPTIONAL :: should_abort
    LOGICAL :: do_abort

    do_abort = .TRUE.
    IF (PRESENT(should_abort)) do_abort = should_abort

    ! Log error to stderr and log file
    IF (myid == master) THEN
      WRITE(*,'(A)') '========== ERROR =========='
      WRITE(*,'(A,I0)') 'Code:    ', code
      WRITE(*,'(A,A)') 'Routine: ', TRIM(routine)
      WRITE(*,'(A,A)') 'Message: ', TRIM(message)
      WRITE(*,'(A)') '==========================='

      ! Also log to file if available
      OPEN(UNIT=999, FILE='error.log', POSITION='APPEND')
      WRITE(999,'(A,I0,A,A,A,A)') 'ERROR ', code, ' in ', TRIM(routine), ': ', TRIM(message)
      CLOSE(999)
    END IF

    IF (do_abort) THEN
      CALL MPI_Abort(MPI_COMM_WORLD, code, ierr)
      STOP
    END IF

  END SUBROUTINE HandleError

END MODULE QuadSc_errors
```

**Usage**:
```fortran
USE QuadSc_errors, ONLY: HandleError, ERR

IF (ierr /= 0) THEN
  CALL HandleError(ERR%ALLOC_FAILED, 'Create_MMat', &
                   'Failed to allocate mass matrix', should_abort=.TRUE.)
END IF
```

### 4.3 Array Bounds Checking During Development

Add to CMake configuration:

```cmake
# Debug build flags
set(CMAKE_Fortran_FLAGS_DEBUG
    "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=bounds -fcheck=array-temps")

# Release build (bounds checking OFF for performance)
set(CMAKE_Fortran_FLAGS_RELEASE
    "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -march=native")
```

## 5. Notes for Developers

*   **Testing:** Any changes to this file require running the full regression suite (`ctest`) because these matrices underpin the entire flow solver.

*   **Dependencies:** This module depends heavily on `pp3d_mpi` and low-level FEM routines (`AP7`, `AP9`, `E013`). Changes here may require interface updates in those assembly routines.

*   **Build Configuration:** Add new utility modules to `cmake/modules/ProjectFiles.cmake`:
    ```cmake
    set(src_util
      ...
      ${CMAKE_SOURCE_DIR}/source/src_util/QuadSc_memory_utils.f90
      ${CMAKE_SOURCE_DIR}/source/src_util/QuadSc_errors.f90
      ...
    )
    ```

*   **Performance Baseline:** Always run baseline timing (Section 3.1) before and after major changes.

*   **Verification Strategy:** See `refactoring_roadmap.md` for detailed verification procedures after each phase.
