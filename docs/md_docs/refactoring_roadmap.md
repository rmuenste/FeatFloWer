# Refactoring Roadmap: QuadSc_def.f90

**Strategy:** Phased Refactoring (Lowest Risk → Highest Risk)
**Last Updated:** November 28, 2025 (Enhanced with verification tools and rollback strategy)

This document outlines the execution plan for refactoring `source/src_quadLS/QuadSc_def.f90`. The guiding principle is to maintain a working build at every step, minimizing "big bang" changes.

## Git Branch Strategy

**CRITICAL**: All refactoring work MUST happen on feature branches, never directly on `master`.

```bash
# Create main refactoring branch from master
git checkout master
git checkout -b feature/quadsc-refactoring

# Create phase-specific sub-branches
git checkout -b feature/quadsc-phase0-baseline
# ... work on Phase 0 ...
git checkout feature/quadsc-refactoring
git merge feature/quadsc-phase0-baseline

git checkout -b feature/quadsc-phase1-extraction
# ... work on Phase 1 ...
git checkout feature/quadsc-refactoring
git merge feature/quadsc-phase1-extraction

# And so on for each phase...
```

**Branch Naming Convention**:
- `feature/quadsc-phase0-baseline` - Baseline generation and diagnostics
- `feature/quadsc-phase1-extraction` - Logic extraction (HYPRE, UMFPACK)
- `feature/quadsc-phase2-deduplication` - Generic assembly routines
- `feature/quadsc-phase3-module-split` - Physical file splitting
- `feature/quadsc-phase4-api-modern` - API modernization

**Rollback Strategy**:
```bash
# If a phase fails verification:
git checkout feature/quadsc-refactoring
git branch -D feature/quadsc-phase2-deduplication  # Delete failed branch
git checkout -b feature/quadsc-phase2-deduplication-v2  # Start over

# Emergency rollback to stable state:
git checkout master  # Always remains stable
```

## Phase 0: Preparation & Baseline

Before touching code, establish a safety net and diagnostic baseline.

### Phase 0.1: Baseline Generation

Run a standard benchmark and save golden outputs:

```bash
# Navigate to build directory
cd build

# Run benchmark (choose appropriate application)
mpirun -np 4 applications/q2p1_devel/q2p1_devel > pe1_golden.log 2>&1

# Save outputs
cp pe1.log ../golden_outputs/pe1_golden.log
cp _data/*.gmv ../golden_outputs/ 2>/dev/null || true
cp _data/*.vtu ../golden_outputs/ 2>/dev/null || true
cp _data/sol_*.dat ../golden_outputs/ 2>/dev/null || true

# Record system state
git log -1 --oneline > ../golden_outputs/git_commit.txt
git diff > ../golden_outputs/git_diff.txt
make --version > ../golden_outputs/compiler_versions.txt
mpirun --version >> ../golden_outputs/compiler_versions.txt
```

### Phase 0.2: Test Check

Ensure all tests pass before starting:

```bash
ctest --output-on-failure > ../golden_outputs/ctest_golden.log 2>&1

# Verify 100% pass rate
if [ $? -ne 0 ]; then
  echo "ERROR: Baseline tests failing! Fix before refactoring."
  exit 1
fi
```

### Phase 0.3: Document Current Behavior (NEW)

**Action Items**:
1. **Document Matrix Storage Format**:
   - Create `docs/md_docs/quadsc_current_implementation.md`
   - Document CSR (Compressed Sparse Row) structure
   - Explain `ColA`, `LdA`, and `a` array relationships
   - Document indexing conventions (0-based vs 1-based)

2. **Document Multigrid Conventions**:
   - `ILEV`, `NLMIN`, `NLMAX` meanings
   - Coarsening direction (fine→coarse)
   - Level numbering scheme

3. **Generate Call Graph**:
   ```bash
   # Extract all CALL statements
   grep -n "CALL Create_" source/src_quadLS/QuadSc_def.f90 > docs/call_graph_quadsc.txt

   # Generate dependency graph (requires ford or doxygen)
   ford quadsc_doc.md  # If using Ford

   # Or manual analysis
   grep -n "^SUBROUTINE\|^FUNCTION" source/src_quadLS/QuadSc_def.f90 | \
       awk '{print $2}' > docs/quadsc_routines.txt
   ```

4. **Add Inline Documentation Headers**:
   - See `def_refactoring.md` Section 6.3 for template
   - Document every major subroutine BEFORE moving it

### Phase 0.4: Add Diagnostic Logging for /16 Mystery (NEW)

Before any refactoring, add diagnostic code to understand the `/16` divider:

**File**: `source/src_quadLS/QuadSc_def.f90` (temporary diagnostic version)

```fortran
! In Create_CMat, around line 1140, ADD THIS BEFORE the /16 division:
IF (myid == 0 .and. ILEV == NLMIN) THEN
  OPEN(UNIT=9999, FILE='_data/umfpack_diagnostic.log', STATUS='REPLACE')
  WRITE(9999,*)
  WRITE(9999,'(A)') '========== UMFPACK Matrix Size Diagnostic =========='
  WRITE(9999,'(A,I0)') 'Multigrid level (NLMIN): ', NLMIN
  WRITE(9999,'(A,I0)') 'Fine grid unknowns (lMat%nu):     ', lMat%nu
  WRITE(9999,'(A,I0)') 'Fine grid non-zeros (lMat%na):    ', lMat%na
  WRITE(9999,'(A,I0)') 'Coarse unknowns (lMat%nu/4):      ', lMat%nu/4
  WRITE(9999,'(A,I0)') 'Coarse non-zeros (lMat%na/16):    ', lMat%na/16
  WRITE(9999,'(A,F10.4)') 'Actual nu ratio:                  ', REAL(lMat%nu)/REAL(lMat%nu/4)
  WRITE(9999,'(A,F10.4)') 'Actual na ratio:                  ', REAL(lMat%na)/REAL(lMat%na/16)
  WRITE(9999,'(A,I0)') 'Average NNZ per row (fine):       ', lMat%na / lMat%nu
  WRITE(9999,'(A,I0)') 'Average NNZ per row (coarse):     ', (lMat%na/16) / (lMat%nu/4)
  WRITE(9999,'(A)') '===================================================='
  WRITE(9999,*)
  CLOSE(9999)
END IF
```

**Run diagnostic on multiple test cases**:
```bash
# Test case 1: 2D uniform mesh
mpirun -np 1 q2p1_devel  # config for 2D
cat _data/umfpack_diagnostic.log

# Test case 2: 3D uniform mesh
mpirun -np 1 q2p1_devel  # config for 3D
cat _data/umfpack_diagnostic.log

# Test case 3: Refined mesh
# ... etc
```

Document findings in `docs/md_docs/def_modernization.md` Section 1.3.1.

### Phase 0.5: Performance Baseline (NEW)

Add timing instrumentation to measure current performance:

```fortran
! Add to beginning of major routines:
REAL(KIND=8) :: t_start_routine, t_end_routine

CALL CPU_TIME(t_start_routine)
! ... routine code ...
CALL CPU_TIME(t_end_routine)

IF (myid == 0) THEN
  OPEN(UNIT=9998, FILE='_data/performance_baseline.log', POSITION='APPEND')
  WRITE(9998,'(A,F10.4)') 'Create_MMat time (sec): ', t_end_routine - t_start_routine
  CLOSE(9998)
END IF
```

**Target routines for timing**:
- `Create_QuadMatStruct`
- `Create_MMat`
- `Create_DiffMat`
- `Create_KMat`
- `Create_SMat`
- `Create_BMat`
- `SetUp_HYPRE_Solver`

**Save baseline**:
```bash
cp _data/performance_baseline.log ../golden_outputs/performance_baseline.log
```

## Enhanced Verification Strategy

### Comparison Tools

Create `tools/compare_outputs.sh`:

```bash
#!/bin/bash
#===============================================================================
# compare_outputs.sh
#
# DESCRIPTION:
#   Compares test outputs against golden baseline with tolerance for
#   floating-point differences.
#
# USAGE:
#   ./tools/compare_outputs.sh <golden_log> <test_log>
#===============================================================================

GOLDEN=$1
TEST=$2
REL_TOL=1e-12
ABS_TOL=1e-14

python3 << EOF
import re
import sys

def compare_logs(golden_file, test_file, rel_tol=$REL_TOL, abs_tol=$ABS_TOL):
    """
    Compare log files allowing for floating-point differences.
    Returns True if files match within tolerance.
    """
    with open(golden_file, 'r') as f1, open(test_file, 'r') as f2:
        line_num = 0
        for line1, line2 in zip(f1, f2):
            line_num += 1

            # Skip lines that are purely textual (no numbers)
            nums1 = [float(x) for x in re.findall(r'[-+]?\d*\.?\d+[eE]?[-+]?\d*', line1)]
            nums2 = [float(x) for x in re.findall(r'[-+]?\d*\.?\d+[eE]?[-+]?\d*', line2)]

            if len(nums1) != len(nums2):
                print(f"MISMATCH at line {line_num}: Different number of values")
                print(f"  Golden: {line1.strip()}")
                print(f"  Test:   {line2.strip()}")
                return False

            # Compare each number with tolerance
            for n1, n2 in zip(nums1, nums2):
                # Relative tolerance: |n1 - n2| <= abs_tol + rel_tol * max(|n1|, |n2|)
                max_val = max(abs(n1), abs(n2))
                tolerance = abs_tol + rel_tol * max_val

                if abs(n1 - n2) > tolerance:
                    print(f"NUMERICAL DIFF at line {line_num}:")
                    print(f"  Expected: {n1}")
                    print(f"  Got:      {n2}")
                    print(f"  Diff:     {abs(n1-n2)}")
                    print(f"  Tolerance: {tolerance}")
                    print(f"  Line (golden): {line1.strip()}")
                    print(f"  Line (test):   {line2.strip()}")
                    return False

    return True

# Run comparison
if compare_logs('$GOLDEN', '$TEST'):
    print("✓ Output matches golden baseline within tolerances")
    sys.exit(0)
else:
    print("✗ Output differs from golden baseline")
    sys.exit(1)
EOF
```

Make executable:
```bash
chmod +x tools/compare_outputs.sh
```

### Verification Checklist (After Each Phase)

Create `tools/verify_phase.sh`:

```bash
#!/bin/bash
#===============================================================================
# verify_phase.sh <phase_number>
#
# DESCRIPTION:
#   Runs comprehensive verification after completing a refactoring phase.
#===============================================================================

PHASE=$1
GOLDEN_DIR="golden_outputs"
TEST_DIR="_data"

echo "=========================================="
echo "Verifying Phase $PHASE"
echo "=========================================="

# Step 1: Compilation
echo "Step 1: Clean rebuild..."
make clean
if ! make -j8 2>&1 | tee build_phase${PHASE}.log; then
  echo "✗ FAIL: Compilation failed"
  exit 1
fi

# Check for warnings
WARNINGS=$(grep -i "warning" build_phase${PHASE}.log | wc -l)
echo "  Compiler warnings: $WARNINGS"
if [ $WARNINGS -gt 0 ]; then
  echo "  WARNING: New compiler warnings detected"
  grep -i "warning" build_phase${PHASE}.log | head -20
fi
echo "✓ PASS: Compilation successful"

# Step 2: Unit tests
echo "Step 2: Running unit tests..."
if ! ctest -R QuadSc --output-on-failure; then
  echo "✗ FAIL: Unit tests failed"
  exit 1
fi
echo "✓ PASS: Unit tests passed"

# Step 3: Full regression suite
echo "Step 3: Running full test suite..."
if ! ctest --output-on-failure > ctest_phase${PHASE}.log 2>&1; then
  echo "✗ FAIL: Regression tests failed"
  cat ctest_phase${PHASE}.log
  exit 1
fi
echo "✓ PASS: All tests passed"

# Step 4: Run benchmark
echo "Step 4: Running benchmark..."
if ! mpirun -np 4 applications/q2p1_devel/q2p1_devel > pe1_phase${PHASE}.log 2>&1; then
  echo "✗ FAIL: Benchmark execution failed"
  exit 1
fi
echo "✓ PASS: Benchmark completed"

# Step 5: Compare outputs
echo "Step 5: Comparing outputs against golden baseline..."
if ! ./tools/compare_outputs.sh $GOLDEN_DIR/pe1_golden.log pe1_phase${PHASE}.log; then
  echo "✗ FAIL: Output differs from baseline"
  exit 1
fi
echo "✓ PASS: Outputs match within tolerance"

# Step 6: Performance check
echo "Step 6: Performance comparison..."
if [ -f "$GOLDEN_DIR/performance_baseline.log" ] && [ -f "$TEST_DIR/performance_baseline.log" ]; then
  python3 << 'EOF'
import re

golden = {}
test = {}

with open('golden_outputs/performance_baseline.log') as f:
    for line in f:
        match = re.match(r'(.+) time \(sec\):\s+([\d.]+)', line)
        if match:
            golden[match.group(1)] = float(match.group(2))

with open('_data/performance_baseline.log') as f:
    for line in f:
        match = re.match(r'(.+) time \(sec\):\s+([\d.]+)', line)
        if match:
            test[match.group(1)] = float(match.group(2))

print("\nPerformance Comparison:")
print(f"{'Routine':<30} {'Golden (s)':<12} {'Test (s)':<12} {'Change (%)':<12}")
print("-" * 70)

regression_detected = False
for routine in golden:
    if routine in test:
        golden_time = golden[routine]
        test_time = test[routine]
        change_pct = 100.0 * (test_time - golden_time) / golden_time if golden_time > 0 else 0.0

        print(f"{routine:<30} {golden_time:<12.4f} {test_time:<12.4f} {change_pct:+11.2f}%")

        if change_pct > 10.0:  # More than 10% slower
            print(f"  WARNING: Performance regression detected!")
            regression_detected = True

if regression_detected:
    print("\n⚠ Performance regression detected (>10% slower)")
    exit(1)
else:
    print("\n✓ Performance within acceptable range (±10%)")
    exit(0)
EOF

  if [ $? -ne 0 ]; then
    echo "✗ FAIL: Performance regression detected"
    exit 1
  fi
else
  echo "  (Performance comparison skipped - missing baseline)"
fi

# Step 7: Memory usage check
echo "Step 7: Memory usage check..."
# Simple check: executable size shouldn't grow significantly
GOLDEN_SIZE=$(stat -c%s "$GOLDEN_DIR/../build/applications/q2p1_devel/q2p1_devel" 2>/dev/null || echo "0")
TEST_SIZE=$(stat -c%s "applications/q2p1_devel/q2p1_devel" 2>/dev/null || echo "0")

if [ $GOLDEN_SIZE -gt 0 ] && [ $TEST_SIZE -gt 0 ]; then
  SIZE_CHANGE=$(python3 -c "print(100.0 * ($TEST_SIZE - $GOLDEN_SIZE) / $GOLDEN_SIZE)")
  echo "  Executable size change: ${SIZE_CHANGE}%"
  if (( $(echo "$SIZE_CHANGE > 20" | bc -l) )); then
    echo "  WARNING: Executable size increased >20%"
  fi
else
  echo "  (Memory check skipped - missing baseline)"
fi

echo ""
echo "=========================================="
echo "✓ Phase $PHASE verification PASSED"
echo "=========================================="
```

Make executable:
```bash
chmod +x tools/verify_phase.sh
```

## Phase 1: Logic Extraction (Low Risk)

**Goal**: Move isolated, heavy logic into separate modules *without* changing the dependency graph of the main application. `QuadSc_def.f90` will essentially become thinner by offloading work to new modules, but it will still be the primary entry point for consumers.

### Step 1.1: Extract HYPRE Interface

**Branch**: `feature/quadsc-phase1-extraction`

**Action**: Move `SetUp_HYPRE_Solver` (and its helper `SORT_DOFs`) to a new module `QuadSc_solver_hypre.f90`.

```fortran
! Create: source/src_quadLS/QuadSc_solver_hypre.f90
MODULE QuadSc_solver_hypre
  USE var_QuadScalar
  USE PP3D_MPI
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: SetUp_HYPRE_Solver

CONTAINS

  SUBROUTINE SetUp_HYPRE_Solver(...)
    ! Move entire routine here
  END SUBROUTINE

  SUBROUTINE SORT_DOFs(...)
    ! Move helper routine here
  END SUBROUTINE

END MODULE QuadSc_solver_hypre
```

**Integration**: In `QuadSc_def.f90`:
```fortran
MODULE def_QuadScalar
  USE QuadSc_solver_hypre, ONLY: SetUp_HYPRE_Solver
  ! ... rest of module ...
END MODULE
```

**CMake**: Add to `cmake/modules/ProjectFiles.cmake`:
```cmake
set(src_quadLS
  ...
  ${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_solver_hypre.f90
  ...
)
```

**Risk**: Low. This routine is self-contained.

**Verification**:
```bash
./tools/verify_phase.sh 1.1
```

### Step 1.2: Extract Coarse Grid Solver

**Action**: Identify the UMFPACK setup code block inside `Create_CMat`. Move it to `QuadSc_solver_coarse.f90`.

**New routine**:
```fortran
! Create: source/src_quadLS/QuadSc_solver_coarse.f90
MODULE QuadSc_solver_coarse
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Setup_UMFPACK_Coarse

CONTAINS

  SUBROUTINE Setup_UMFPACK_Coarse(ILEV, lMat, crsSTR, UMF_CMat)
    ! Extract UMFPACK setup from Create_CMat
    ! This handles the /16 mystery - see diagnostics from Phase 0.4
  END SUBROUTINE

END MODULE
```

**Integration**: `Create_CMat` calls the new external routine.

**Risk**: Low/Medium. Requires careful handling of local variables (`crsSTR`, `UMF_CMat`) passed to the new routine.

**Verification**:
```bash
./tools/verify_phase.sh 1.2
```

## Phase 2: Internal Deduplication (Medium Risk)

**Goal**: Reduce code size by identifying common patterns. We will implement "Generic Assembly Routines" **inside** the current file first (or a helper module) to prove logic correctness before splitting files.

### Step 2.1: Mass Matrix Unification

**Action**: Create `Assemble_Mass_Generic`. Refactor `Create_MMat` and `Create_MRhoMat` to call this generic routine.

**Risk**: Medium. Must ensure `E013` element integration behaves identically when passed as an argument vs hardcoded.

**Performance Check**: Use timing comparison (see `def_modernization.md` Section 3.2):
```fortran
! Keep both versions temporarily
CALL CPU_TIME(t1)
CALL Create_MMat_Old(...)
CALL CPU_TIME(t2)

CALL CPU_TIME(t3)
CALL Create_MMat_Generic(...)
CALL CPU_TIME(t4)

IF (myid == 0) THEN
  WRITE(*,'(A,F8.4,A)') 'Old version: ', t2-t1, ' sec'
  WRITE(*,'(A,F8.4,A)') 'New version: ', t4-t3, ' sec'
  overhead = 100.0 * (t4-t3) / (t2-t1) - 100.0
  WRITE(*,'(A,F8.2,A)') 'Overhead:    ', overhead, '%'

  IF (overhead > 10.0) THEN
    WRITE(*,*) 'WARNING: Performance regression >10%'
  END IF
END IF
```

**Verification**:
```bash
./tools/verify_phase.sh 2.1
```

### Step 2.2: Diffusion Matrix Unification

**Action**: Create `Assemble_Diff_Generic`. Refactor `Create_hDiffMat` and `Create_ConstDiffMat`.

**Risk**: Low. These routines are structurally identical.

**Verification**:
```bash
./tools/verify_phase.sh 2.2
```

### Step 2.3: Parallel Matrix Setup

**Action**: Merge logic from `Create_QuadLinParMatStruct` and `Fill_QuadLinParMat`.

**Risk**: Medium. Logic flow for allocation vs filling needs to be handled by a clean flag or argument.

**Verification**:
```bash
./tools/verify_phase.sh 2.3
```

## Phase 3: Structural Split (High Risk)

**Goal**: Physically break the file into multiple compilation units. This changes `CMakeLists.txt` and potentially `USE` statements in other files.

**CRITICAL**: Follow dependency graph from `def_refactoring.md` Section 1.1 to avoid circular dependencies.

### Step 3.1: Create `QuadSc_struct.f90`

**Action**: Move `Create_QuadMatStruct`, `Create_QuadLinMatStruct`, `Create_LinMatStruct` to a new module.

**Dependency Check**:
```bash
# After creating module, verify no circular dependencies
grep -n "USE QuadSc" source/src_quadLS/QuadSc_*.f90 | grep -v "!.*USE"

# Manually verify hierarchy is preserved
```

**Verification**:
```bash
./tools/verify_phase.sh 3.1
```

### Step 3.2: Create `QuadSc_assembly.f90`

**Action**: Move the newly deduplicated assembly routines (`Create_MMat`, `Create_KMat`, etc.) to this module.

**Dependency**: Will depend on `QuadSc_struct` or just share `var_QuadScalar`. Must NOT create circular dependency.

**Verification**:
```bash
./tools/verify_phase.sh 3.2
```

### Step 3.3: Facade Pattern (The Switch)

**Action**: Turn `def_QuadScalar` into a "Facade Module" that simply `USE`s and `PUBLIC`s the new modules.

**Benefit**: This prevents breaking the hundreds of calls in `QuadSc_main.f90` that expect `USE def_QuadScalar`.

**Code Example**:
```fortran
MODULE def_QuadScalar
!===============================================================================
! MODULE: def_QuadScalar (Facade)
!
! DESCRIPTION:
!   Facade module that re-exports all QuadSc functionality.
!   This allows existing code to continue using "USE def_QuadScalar" without
!   modification, while internally the implementation is split across multiple modules.
!===============================================================================
  USE QuadSc_struct, ONLY: Create_QuadMatStruct, Create_QuadLinMatStruct, &
                            Create_LinMatStruct, Create_QuadLinParMatStruct
  USE QuadSc_assembly, ONLY: Create_MMat, Create_DiffMat, Create_KMat, &
                              Create_SMat, Create_BMat, Assemble_Mass_Generic
  USE QuadSc_solver_hypre, ONLY: SetUp_HYPRE_Solver
  USE QuadSc_solver_coarse, ONLY: Setup_UMFPACK_Coarse
  USE QuadSc_system, ONLY: Matdef_general_QuadScalar, Matdef_General_LinScalar

  IMPLICIT NONE
  PUBLIC  ! Re-export everything

END MODULE def_QuadScalar
```

**Verification**:
```bash
# This is the most critical verification - entire system must still work
./tools/verify_phase.sh 3.3

# Additional check: ensure no external code needs modification
grep -r "USE def_QuadScalar" source/ applications/ | wc -l
# This number should be unchanged from before refactoring
```

## Phase 4: API Modernization (Highest Risk)

**Goal**: Remove global pointer dependency (`qMat => ...`). This is technically "Modernization" but fits the refactoring arc.

*   **Step 4.1**: Update interfaces in `QuadSc_assembly.f90` to accept `INTENT(INOUT)` matrix types.
*   **Step 4.2**: Update the "Facade" module to pass the globals into these clean interfaces.

**Implementation**: See `def_refactoring.md` Section 4 for details.

**Verification**:
```bash
./tools/verify_phase.sh 4
```

---

## Code Review Checkpoints (NEW)

After each phase, perform code review:

### Code Review Checklist

- [ ] **Compilation**:
  - [ ] Clean build with no errors
  - [ ] No new compiler warnings
  - [ ] All targets build successfully

- [ ] **Error Handling**:
  - [ ] All `ALLOCATE` statements have `STAT=` checking
  - [ ] All file I/O has `IOSTAT=` checking
  - [ ] All error paths properly propagate to caller

- [ ] **Code Quality**:
  - [ ] No magic numbers (constants are named `PARAMETER`s)
  - [ ] All subroutines have documentation headers
  - [ ] All module variables have explanatory comments
  - [ ] `INTENT` declared for all subroutine parameters

- [ ] **Dependency Management**:
  - [ ] No circular dependencies (verify with `grep` check)
  - [ ] Dependency graph matches design (see `def_refactoring.md` Section 1.1)
  - [ ] No global pointers modified outside designated modules (except Phase 4)

- [ ] **Testing**:
  - [ ] All unit tests pass
  - [ ] Full regression suite passes
  - [ ] Benchmark output matches golden baseline
  - [ ] Performance within ±10% of baseline

- [ ] **Documentation**:
  - [ ] Changes documented in commit message
  - [ ] Module documentation updated
  - [ ] README updated if public API changed

- [ ] **Git Hygiene**:
  - [ ] Commits are atomic (one logical change per commit)
  - [ ] Commit messages are descriptive
  - [ ] No debug code or commented-out blocks committed

## Emergency Rollback Procedure

If verification fails and cannot be quickly fixed:

```bash
# 1. Save diagnostic information
cp build_phase${PHASE}.log ../failed_phase${PHASE}/
cp pe1_phase${PHASE}.log ../failed_phase${PHASE}/
git diff > ../failed_phase${PHASE}/changes.diff

# 2. Rollback to last stable state
git checkout feature/quadsc-refactoring  # Parent branch
git branch -D feature/quadsc-phase${PHASE}  # Delete failed branch

# 3. Analyze failure
less ../failed_phase${PHASE}/build_phase${PHASE}.log
less ../failed_phase${PHASE}/pe1_phase${PHASE}.log

# 4. Document what went wrong
echo "Phase ${PHASE} failed due to: [REASON]" > ../failed_phase${PHASE}/postmortem.txt

# 5. Start over with lessons learned
git checkout -b feature/quadsc-phase${PHASE}-v2
```

## Final Merge to Master

Only after ALL phases complete successfully:

```bash
# Ensure feature branch is up to date with master
git checkout master
git pull origin master
git checkout feature/quadsc-refactoring
git rebase master

# Final comprehensive verification
./tools/verify_phase.sh FINAL

# Merge to master
git checkout master
git merge --no-ff feature/quadsc-refactoring -m "Refactor QuadSc_def.f90 into modular components

Complete refactoring of QuadSc_def.f90 (~4500 lines) into 5 focused modules:
- QuadSc_struct.f90: Sparsity pattern allocation
- QuadSc_assembly.f90: Matrix assembly with generic routines
- QuadSc_solver_hypre.f90: HYPRE interface
- QuadSc_solver_coarse.f90: Coarse grid solver (UMFPACK)
- QuadSc_system.f90: High-level system definition
- def_QuadScalar: Facade module (maintains backward compatibility)

All phases verified against golden baseline.
Performance within ±5% of original implementation.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"

# Tag the successful refactoring
git tag -a v2.0-quadsc-refactored -m "Completed QuadSc_def.f90 refactoring"

# Push to remote
git push origin master
git push origin v2.0-quadsc-refactored
```

## Notes

*   **Patience is Key**: This refactoring may take weeks. Don't rush phases.
*   **Communication**: If multiple developers are involved, coordinate to avoid conflicts.
*   **Backup**: Keep the original `QuadSc_def.f90` in a backup directory for reference.
*   **Testing is Non-negotiable**: Never skip verification steps.
