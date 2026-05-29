# Matrix Assembly Guide for FeatFloWer

**Author:** Generated during QuadSc refactoring (Phase 3.2)
**Date:** 2025
**Purpose:** Guide for understanding matrix assembly in the Q2/P1 finite element solver

---

## Table of Contents

1. [Overview](#overview)
2. [Module Architecture](#module-architecture)
3. [Matrix Types and Purpose](#matrix-types-and-purpose)
4. [Assembly Workflow](#assembly-workflow)
5. [COMMON Blocks (Legacy Global State)](#common-blocks-legacy-global-state)
6. [Generic vs Specific Assembly](#generic-vs-specific-assembly)
7. [Multigrid Hierarchy](#multigrid-hierarchy)
8. [Parallel Matrix Handling](#parallel-matrix-handling)
9. [CSR Matrix Format](#csr-matrix-format)
10. [Code Examples](#code-examples)
11. [Troubleshooting](#troubleshooting)

---

## Overview

FeatFloWer uses a **Q2/P1 finite element discretization** for incompressible Navier-Stokes equations:

- **Q2 (Quadratic):** Velocity components (u, v, w) - 27 DOFs per hexahedral element
- **P1 (Linear):** Pressure (p) - 4 DOFs per hexahedral element (vertices only)

**Matrix assembly** refers to the process of:
1. **Structure allocation** - Creating sparsity patterns (which entries are non-zero)
2. **Coefficient filling** - Computing actual matrix values from element contributions

**Storage format:** Compressed Sparse Row (CSR) with `LdA` (row pointers), `ColA` (column indices), and coefficient arrays.

---

## Module Architecture

### Refactored Structure (Post-Phase 3)

```
var_QuadScalar (data module)
    ↓
QuadSc_struct (structure allocation)
    ↓
QuadSc_assembly (coefficient assembly)
    ↓
def_QuadScalar (facade module - backward compatibility)
```

### Key Modules

#### `var_QuadScalar` (Data Structures)
```fortran
! Multigrid matrix structures (all levels NLMIN:NLMAX)
TYPE(TMatrix), ALLOCATABLE :: mg_qMat(:)     ! Q2 velocity matrices
TYPE(TMatrix), ALLOCATABLE :: mg_lMat(:)     ! P1 pressure matrices (local)
TYPE(TMatrix), ALLOCATABLE :: mg_qlMat(:)    ! Q2→P1 gradient (B matrix)
TYPE(TMatrix), ALLOCATABLE :: mg_lqMat(:)    ! P1→Q2 divergence (B^T matrix)
TYPE(TParMatrix), ALLOCATABLE :: mg_lPMat(:) ! P1 parallel matrices

! Coefficient arrays (multigrid hierarchy)
TYPE(mg_Matrix), ALLOCATABLE :: mg_Mmat(:)   ! Mass matrix
TYPE(mg_Matrix), ALLOCATABLE :: mg_Dmat(:)   ! Diffusion matrix
TYPE(mg_Matrix), ALLOCATABLE :: mg_Kmat(:)   ! Convection matrix
TYPE(mg_Matrix), ALLOCATABLE :: mg_CMat(:)   ! Pressure Schur complement
TYPE(mg_Matrix), ALLOCATABLE :: mg_BXMat(:)  ! Gradient (X-component)
TYPE(mg_Matrix), ALLOCATABLE :: mg_BYMat(:)  ! Gradient (Y-component)
TYPE(mg_Matrix), ALLOCATABLE :: mg_BZMat(:)  ! Gradient (Z-component)
```

#### `QuadSc_struct.f90` (Structure Allocation)
```fortran
PUBLIC :: Create_QuadMatStruct        ! Q2 velocity matrix structure
PUBLIC :: Create_QuadLinMatStruct     ! Q2↔P1 gradient/divergence structure
PUBLIC :: Create_LinMatStruct         ! P1 pressure matrix structure (local)
PUBLIC :: Create_ParLinMatStruct      ! P1 pressure matrix structure (parallel)
```

**Purpose:** Allocate sparsity patterns without computing coefficients.

#### `QuadSc_assembly.f90` (Coefficient Assembly)
```fortran
! Generic routines (from Phase 2 deduplication)
PUBLIC :: Assemble_Mass_Generic
PUBLIC :: Assemble_Diffusion_Alpha_Generic
PUBLIC :: Assemble_ParallelMatrix_Generic

! Specific matrix builders
PUBLIC :: Create_MMat, Create_MRhoMat           ! Mass matrices
PUBLIC :: Create_DiffMat, Create_hDiffMat, ...  ! Diffusion matrices
PUBLIC :: Create_BMat, Create_CMat              ! Velocity-pressure coupling
PUBLIC :: Create_KMat, Create_SMat              ! Convection matrices
```

**Purpose:** Fill matrix coefficients by calling element kernels.

---

## Matrix Types and Purpose

### 1. Mass Matrices (Inertia)

| Matrix | Symbol | Equation Role | Routine |
|--------|--------|---------------|---------|
| `mg_Mmat` | M | ∂u/∂t → M·u_t | `Create_MMat` |
| `mg_MlMat` | M_l | Lumped mass (diagonal) | (computed in `Create_MMat`) |
| `mg_MlRhomat` | ρM_l | Density-weighted lumped | `Create_MRhoMat` |

**Element kernel:** `BuildMMat`, `BuildMRhoMat` (call `AB07` - mass matrix integration)

**Lumping:** Row sums for time-stepping stability:
```fortran
DO I=1,qMat%nu
  DML = 0d0
  DO J=qMat%LdA(I), qMat%LdA(I+1)-1
    DML = DML + mg_Mmat(ILEV)%a(J)
  END DO
  mg_MlMat(ILEV)%a(I) = DML
END DO
```

### 2. Diffusion Matrices (Viscosity)

| Matrix | Symbol | Equation Role | Routine |
|--------|--------|---------------|---------|
| `mg_Dmat` | D | ν∇²u → D·u | `Create_DiffMat` |
| `mg_hDmat` | D_h | h-dependent (LBB stabilization) | `Create_hDiffMat` |
| `mg_ConstDMat` | D_const | Constant diffusion | `Create_ConstDiffMat` |

**Element kernel:** `DIFFQ2_alpha`, `DIFFQ2_NEWT`, `DIFFQ2_NNEWT` (call `AD07` - diffusion integration)

**Alpha parameter:**
- α = 0.0 → constant diffusion
- α = 1.0 → h-dependent (mesh size scaling)
- α = f(x,y,z) → user-supplied scalar field

### 3. Convection Matrices (Advection)

| Matrix | Symbol | Equation Role | Routine |
|--------|--------|---------------|---------|
| `mg_Kmat` | K | (u·∇)u → K·u | `Create_KMat` |
| `mg_Smat` | S | Streamline diffusion (9 components) | `Create_SMat` |

**Element kernel:** `CONVQ2` (convection), `CUBATURESTRESS` (stress tensor)

**Streamline diffusion:** 9 matrices (S11, S12, ..., S33) for anisotropic stabilization.

### 4. Velocity-Pressure Coupling

| Matrix | Symbol | Equation Role | Routine |
|--------|--------|---------------|---------|
| `mg_BXMat` | B_x | ∇p → B^T·p (X-gradient) | `Create_BMat` |
| `mg_BYMat` | B_y | ∇p → B^T·p (Y-gradient) | `Create_BMat` |
| `mg_BZMat` | B_z | ∇p → B^T·p (Z-gradient) | `Create_BMat` |
| `mg_BTXMat` | B_x^T | ∇·u → B·u (X-divergence) | `Create_BMat` |
| `mg_CMat` | C | Schur: B^T M^{-1} B | `Create_CMat` |

**Element kernel:** `Build_BMatP1`, `Build_BTMatP1`, `Get_CMat`

**Schur complement:**
```
C = B^T · M_l^{-1} · B
```
Used in pressure Poisson solve.

---

## Assembly Workflow

### Typical Initialization Sequence

```fortran
! 1. Structure allocation (done once)
CALL Create_QuadMatStruct()        ! Allocate mg_qMat sparsity
CALL Create_QuadLinMatStruct()     ! Allocate mg_qlMat, mg_lqMat
CALL Create_LinMatStruct()         ! Allocate mg_lMat
CALL Create_ParLinMatStruct(...)   ! Allocate mg_lPMat (MPI)

! 2. Coefficient assembly (done once or per time step)
CALL Create_MMat()                 ! Fill mass matrix
CALL Create_DiffMat(myScalar)      ! Fill diffusion (depends on flow)
CALL Create_BMat()                 ! Fill gradient/divergence
CALL Create_CMat(...)              ! Fill Schur complement

! 3. Time-dependent updates (every time step or Newton iteration)
CALL Create_KMat(myScalar)         ! Convection depends on velocity
```

### Matrix Renewal (Time-Stepping)

Some matrices change during simulation:
- **Mass matrix:** Constant (unless density changes)
- **Diffusion:** Changes if viscosity is non-Newtonian
- **Convection:** **Always changes** (depends on current velocity)
- **Gradient/Divergence:** Constant (unless mesh adapts)

---

## COMMON Blocks (Legacy Global State)

### ⚠️ IMPORTANT: Global Variables via COMMON Blocks

FeatFloWer uses **Fortran COMMON blocks** extensively. These are **global variables** accessible across all routines without explicit passing.

#### Key COMMON Blocks

**From `CBasicMG.inc`:**
```fortran
COMMON /MG/  NLMIN,NLMAX,ILEV
INTEGER      NLMIN     ! Coarsest multigrid level
INTEGER      NLMAX     ! Finest multigrid level (working level)
INTEGER      ILEV      ! Current active level
```

**From `CBasicElem.inc`:**
```fortran
COMMON /ELEM/ DU1,DU2,DU3,DCMASS,KDFG,IDFL,NEL,NVT,NAT,NET,NVEL,NEEL,NBCT
INTEGER       NEL      ! Number of elements
INTEGER       NVT      ! Number of vertices
INTEGER       NEEL     ! Number of edges
INTEGER       KDFG(*)  ! Degrees of freedom per element
```

#### How COMMON Blocks Affect Matrix Assembly

**Problem:** Routines don't explicitly show their dependencies.

**Example:**
```fortran
SUBROUTINE Create_MMat()
  ! Appears to take no arguments, but actually uses:
  ! - ILEV (from COMMON /MG/)
  ! - NEL, NVT (from COMMON /ELEM/)
  ! - mg_qMat (from var_QuadScalar module)
  ! - mg_mesh (from var_QuadScalar module)

  DO ILEV=NLMIN,NLMAX  ! Uses COMMON block variables
    CALL SETLEV(2)      ! Sets COMMON block ELEM variables for this level
    qMat => mg_qMat(ILEV)

    ! NEL is now set to this level's element count by SETLEV
    CALL BuildMMat(..., NEL, ...)
  END DO
END SUBROUTINE
```

#### `SETLEV(2)` - The Magic Updater

**Purpose:** Updates COMMON block variables for the current multigrid level.

```fortran
CALL SETLEV(2)  ! Argument 2 = "QuadScalar" context
! After this call, COMMON blocks contain:
! - NEL = mg_mesh%level(ILEV)%nel
! - NVT = mg_mesh%level(ILEV)%nvt
! - KVERT = mg_mesh%level(ILEV)%kvert (pointer)
! - etc.
```

**Critical pattern in all assembly routines:**
```fortran
DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)         ! MUST call before using NEL, NVT, etc.
  qMat => mg_qMat(ILEV)  ! Get matrix structure for this level

  ! Now safe to use COMMON block variables
  CALL SomeKernel(..., NEL, NVT, ...)
END DO

! Restore to finest level
ILEV=NLMAX
CALL SETLEV(2)
```

#### Why This Is Confusing

1. **Hidden dependencies:** Function signatures don't show all inputs
2. **Global state mutation:** `SETLEV` changes many variables simultaneously
3. **Order matters:** Forget `SETLEV` → wrong NEL value → crash or wrong results
4. **Hard to trace:** IDEs can't track COMMON block usage

#### Best Practices

✅ **DO:**
- Always call `SETLEV(2)` when changing `ILEV`
- Restore `ILEV=NLMAX` and `SETLEV(2)` after loops
- Document COMMON block usage in routine headers

❌ **DON'T:**
- Assume `NEL` has a certain value without checking `ILEV`
- Modify `ILEV` without calling `SETLEV`
- Pass COMMON variables as arguments (redundant, confusing)

---

## Generic vs Specific Assembly

### Phase 2 Deduplication: Generic Routines

**Problem:** Original code had massive duplication:
- `Create_MMat` and `Create_MRhoMat`: 95% identical
- `Create_hDiffMat`, `Create_ConstDiffMat`, `Create_DiffMat`: similar structure

**Solution:** Extract common logic into generic routines.

#### Example: `Assemble_Mass_Generic`

```fortran
SUBROUTINE Assemble_Mass_Generic(use_density, mg_MlMatrix, mg_MlPMatrix, &
                                  density_opt, label)
  LOGICAL, INTENT(IN) :: use_density
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_MlMatrix(NLMIN:NLMAX)
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_MlPMatrix(NLMIN:NLMAX)
  TYPE(mg_dVector), INTENT(IN), OPTIONAL :: density_opt(NLMIN:NLMAX)
  CHARACTER(LEN=*), INTENT(IN) :: label

  DO ILEV=NLMIN,NLMAX
    CALL SETLEV(2)
    qMat => mg_qMat(ILEV)

    ! Allocate and zero
    IF (.not.ALLOCATED(mg_Mmat(ILEV)%a)) ALLOCATE(mg_Mmat(ILEV)%a(qMat%na))
    mg_Mmat(ILEV)%a = 0d0

    ! Call appropriate kernel
    IF (use_density) THEN
      CALL BuildMRhoMat(density_opt(ILEV)%x, mg_Mmat(ILEV)%a, ...)
    ELSE
      CALL BuildMMat(mg_Mmat(ILEV)%a, ...)
    END IF

    ! Lump matrix (row sums) - identical for both variants
    DO I=1,qMat%nu
      DML = SUM(mg_Mmat(ILEV)%a(qMat%LdA(I):qMat%LdA(I+1)-1))
      mg_MlMatrix(ILEV)%a(I) = DML
    END DO

    ! Parallel sync
    mg_MlPMatrix(ILEV)%a = mg_MlMatrix(ILEV)%a
    CALL E013SUM(mg_MlPMatrix(ILEV)%a)
  END DO
END SUBROUTINE
```

#### Specific Wrappers

```fortran
SUBROUTINE Create_MMat()
  IF (.not.ALLOCATED(mg_MlMat))  ALLOCATE(mg_MlMat(NLMIN:NLMAX))
  IF (.not.ALLOCATED(mg_MlPMat)) ALLOCATE(mg_MlPMat(NLMIN:NLMAX))

  ! Call generic with use_density=.FALSE.
  CALL Assemble_Mass_Generic(.FALSE., mg_MlMat, mg_MlPMat, label="[M] & [Ml]")

  ! Set pointers to NLMAX
  qMat   => mg_qMat(NLMAX)
  Mmat   => mg_Mmat(NLMAX)%a
  MlMat  => mg_MlMat(NLMAX)%a
  MlPMat => mg_MlPMat(NLMAX)%a
END SUBROUTINE

SUBROUTINE Create_MRhoMat()
  IF (.not.ALLOCATED(mg_MlRhomat))  ALLOCATE(mg_MlRhomat(NLMIN:NLMAX))
  IF (.not.ALLOCATED(mg_MlRhoPmat)) ALLOCATE(mg_MlRhoPmat(NLMIN:NLMAX))

  ! Call generic with use_density=.TRUE., pass density array
  CALL Assemble_Mass_Generic(.TRUE., mg_MlRhomat, mg_MlRhoPmat, &
                             mgDensity, "[MRho] & [MlRho]")

  ! Special case: zero mass matrix for steady-state
  IF (bSteadyState) THEN
    DO ILEV_save=NLMIN,NLMAX
      mg_Mmat(ILEV_save)%a = 0d0
    END DO
  END IF

  ! Set pointers
  qMat      => mg_qMat(NLMAX)
  Mmat      => mg_Mmat(NLMAX)%a
  MlRhomat  => mg_MlRhomat(NLMAX)%a
  MlRhoPmat => mg_MlRhoPmat(NLMAX)%a
END SUBROUTINE
```

**Benefits:**
- Reduced duplication: ~60 lines → ~20 lines per variant
- Single point of maintenance for common logic
- Clear separation: wrapper handles allocation, generic handles assembly

---

## Multigrid Hierarchy

### Level Indexing

```
NLMIN = 1  (coarsest level, few DOFs, used in MG preconditioner)
   ↓
NLMIN+1
   ↓
  ...
   ↓
NLMAX = 3  (finest level, working mesh, where solution lives)
```

**Key insight:** All matrices exist on all levels for multigrid efficiency.

### Level-Dependent Data

```fortran
! Mesh data (different for each level)
mg_mesh%level(ILEV)%nel      ! Number of elements at level ILEV
mg_mesh%level(ILEV)%nvt      ! Number of vertices
mg_mesh%level(ILEV)%kvert    ! Element→vertex connectivity
mg_mesh%level(ILEV)%dcorvg   ! Vertex coordinates

! Matrix structures (different sparsity for each level)
mg_qMat(ILEV)%nu             ! Number of DOFs (unknowns)
mg_qMat(ILEV)%na             ! Number of non-zeros
mg_qMat(ILEV)%LdA            ! Row pointers
mg_qMat(ILEV)%ColA           ! Column indices

! Coefficient arrays (different values for each level)
mg_Mmat(ILEV)%a              ! Mass matrix coefficients at level ILEV
mg_Dmat(ILEV)%a              ! Diffusion coefficients
```

### Typical Multigrid Loop

```fortran
DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)  ! Update COMMON blocks for this level
  qMat => mg_qMat(ILEV)

  ! Allocate for this level
  IF (.not.ALLOCATED(mg_Mmat(ILEV)%a)) THEN
    ALLOCATE(mg_Mmat(ILEV)%a(qMat%na))
  END IF

  ! Assemble for this level (uses level-specific mesh)
  CALL BuildMMat(mg_Mmat(ILEV)%a, qMat%na, qMat%ColA, qMat%LdA, &
                 mg_mesh%level(ILEV)%kvert, &  ! Level-specific
                 mg_mesh%level(ILEV)%karea, &
                 mg_mesh%level(ILEV)%kedge, &
                 mg_mesh%level(ILEV)%dcorvg, &
                 E013)
END DO

! Restore to finest level (CRITICAL!)
ILEV=NLMAX
CALL SETLEV(2)
qMat => mg_qMat(NLMAX)
Mmat => mg_Mmat(NLMAX)%a  ! Pointers now point to finest level
```

---

## Parallel Matrix Handling

### Local vs Parallel Matrices

**Problem:** In MPI domain decomposition, some DOFs are on **different processes**.

| Matrix Type | Scope | Example |
|-------------|-------|---------|
| Local (`mg_lMat`) | Only DOFs owned by this process | On-process pressure |
| Parallel (`mg_lPMat`) | Off-process contributions | Pressure at halo cells |

### Parallel Gradient Matrix (B Matrix)

**Structure allocation:**
```fortran
CALL Create_QuadLinParMatStruct(myPLinSc)
! Internally calls:
CALL Assemble_ParallelMatrix_Generic(allocate_structure=.TRUE., myPLinSc_opt=myPLinSc)
```

**What it does:**
```fortran
DO ILEV = NLMIN, NLMAX
  CALL SETLEV(2)

  ! Initialize MPI communication pattern
  CALL ParPresComm_Init(qlMat%ColA, qlMat%LdA, qlMat%nu, &
                        mg_mesh%level(ILEV)%nel, ILEV)

  ! Allocate parallel matrix structure
  mg_qlPMat(ILEV)%nu = qlMat%nu
  ALLOCATE(mg_qlPMat(ILEV)%LdA(mg_qlPMat(ILEV)%nu+1))

  ! Get row pointers for off-process entries
  CALL Create_ParB_LD(mg_qlPMat(ILEV)%LdA, qlMat%LdA, qlMat%nu, ILEV)

  MatSize = mg_qlPMat(ILEV)%LdA(mg_qlPMat(ILEV)%nu+1)

  ! Allocate column indices and coefficients
  ALLOCATE(mg_qlPMat(ILEV)%ColA(MatSize))
  ALLOCATE(mg_BXPMat(ILEV)%a(MatSize))
  ALLOCATE(mg_BYPMat(ILEV)%a(MatSize))
  ALLOCATE(mg_BZPMat(ILEV)%a(MatSize))
END DO
```

**Coefficient filling:**
```fortran
CALL Fill_QuadLinParMat()
! Internally calls:
CALL Assemble_ParallelMatrix_Generic(allocate_structure=.FALSE.)
! Only refills mg_BXPMat, mg_BYPMat, mg_BZPMat coefficients
```

### MPI Synchronization

**Mass matrix lumping + sync:**
```fortran
! Local lumped mass
mg_MlMatrix(ILEV)%a(I) = DML  ! Row sum

! Parallel lumped mass (synchronized across processes)
mg_MlPMatrix(ILEV)%a = mg_MlMatrix(ILEV)%a
CALL E013SUM(mg_MlPMatrix(ILEV)%a)  ! MPI_Allreduce equivalent
```

---

## CSR Matrix Format

### Structure

**Compressed Sparse Row (CSR)** stores only non-zero entries:

```
Example 4×4 matrix:
  [2.0  0.0  1.0  0.0]
  [0.0  3.0  0.0  4.0]
  [5.0  0.0  6.0  0.0]
  [0.0  7.0  0.0  8.0]

CSR representation:
LdA  = [1, 3, 5, 7, 9]           (row pointers, length nu+1)
ColA = [1, 3, 2, 4, 1, 3, 2, 4]  (column indices, length na)
a    = [2.0, 1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]  (coefficients, length na)

nu = 4 (number of rows/DOFs)
na = 8 (number of non-zeros)
```

### FeatFloWer Matrix Type

```fortran
TYPE TMatrix
  INTEGER :: nu              ! Number of unknowns (rows)
  INTEGER :: na              ! Number of non-zero entries
  INTEGER, POINTER :: LdA(:) ! Row pointers (nu+1 entries)
  INTEGER, POINTER :: ColA(:)! Column indices (na entries)
END TYPE

TYPE mg_Matrix
  REAL*8, POINTER :: a(:)    ! Coefficient array (na entries)
END TYPE

! Usage:
TYPE(TMatrix) :: qMat
TYPE(mg_Matrix) :: mg_Mmat(NLMIN:NLMAX)

! Access row I:
DO J = qMat%LdA(I), qMat%LdA(I+1)-1
  column = qMat%ColA(J)
  value  = mg_Mmat(ILEV)%a(J)
  ! This is matrix entry (I, column) = value
END DO
```

### Why Separate Structure and Coefficients?

**Structure (`qMat`):** Sparsity pattern (which entries exist)
- Allocated once during initialization
- Shared by multiple matrices (M, D, K all use `mg_qMat` structure)

**Coefficients (`mg_Mmat%a`):** Actual values
- Can be recomputed without reallocating `LdA`, `ColA`
- Different physical operators share the same sparsity

**Example:**
```fortran
! Structure allocated once
CALL Create_QuadMatStruct()  ! Allocates mg_qMat(ILEV)%LdA, %ColA

! Multiple coefficient arrays use same structure
CALL Create_MMat()           ! Allocates mg_Mmat(ILEV)%a(qMat%na)
CALL Create_DiffMat(...)     ! Allocates mg_Dmat(ILEV)%a(qMat%na)
CALL Create_KMat(...)        ! Allocates mg_Kmat(ILEV)%a(qMat%na)

! All three use qMat%LdA and qMat%ColA for sparsity!
```

---

## Code Examples

### Example 1: Adding a New Matrix Type

**Goal:** Create a custom diffusion matrix with spatially-varying coefficient α(x,y,z).

```fortran
! Step 1: Declare in var_QuadScalar.f90
TYPE(mg_Matrix), ALLOCATABLE :: mg_CustomDmat(:)

! Step 2: Add assembly routine to QuadSc_assembly.f90
SUBROUTINE Create_CustomDiffMat(alpha_field)
  TYPE(mg_dVector), INTENT(IN) :: alpha_field(NLMIN:NLMAX)
  EXTERNAL E013

  IF (.not.ALLOCATED(mg_CustomDmat)) ALLOCATE(mg_CustomDmat(NLMIN:NLMAX))

  DO ILEV=NLMIN,NLMAX
    CALL SETLEV(2)
    qMat => mg_qMat(ILEV)

    ! Allocate coefficient array (shares qMat structure)
    IF (.not.ALLOCATED(mg_CustomDmat(ILEV)%a)) THEN
      ALLOCATE(mg_CustomDmat(ILEV)%a(qMat%na))
    END IF
    mg_CustomDmat(ILEV)%a = 0d0

    ! Call element kernel with custom alpha
    CALL DIFFQ2_CustomAlpha(alpha_field(ILEV)%x, &
         mg_CustomDmat(ILEV)%a, qMat%na, qMat%ColA, qMat%LdA, &
         mg_mesh%level(ILEV)%kvert, &
         mg_mesh%level(ILEV)%karea, &
         mg_mesh%level(ILEV)%kedge, &
         mg_mesh%level(ILEV)%dcorvg, &
         E013)
  END DO

  ! Restore to NLMAX
  ILEV=NLMAX
  CALL SETLEV(2)
  qMat => mg_qMat(NLMAX)
  CustomDmat => mg_CustomDmat(NLMAX)%a
END SUBROUTINE
```

### Example 2: Matrix Renewal During Time-Stepping

```fortran
! In QuadSc_main.f90 time-stepping loop:
DO iTimeStep = 1, nTimeSteps

  ! Update velocity-dependent matrices
  CALL Create_KMat(myScalar)  ! Convection depends on current velocity

  ! For non-Newtonian flows, update diffusion too
  IF (bNonNewtonian) THEN
    CALL Create_DiffMat(myScalar)  ! Viscosity depends on strain rate
  END IF

  ! Assemble system matrix: A = M/dt + K + D
  CALL AssembleSystemMatrix()

  ! Solve
  CALL MG_Solver(...)

  ! Update solution
  myScalar%valU = myScalar%valU + deltaU
END DO
```

### Example 3: Debugging Matrix Assembly

```fortran
! In QuadSc_assembly.f90, after assembly:
SUBROUTINE DebugPrintMatrix(label, mat_struct, coeffs, level)
  CHARACTER(LEN=*) :: label
  TYPE(TMatrix) :: mat_struct
  REAL*8 :: coeffs(:)
  INTEGER :: level
  INTEGER :: I, J, col

  IF (myid.ne.0) RETURN  ! Only master process

  WRITE(*,*) "=== Matrix Debug: ", TRIM(label), " at level ", level, " ==="
  WRITE(*,*) "Dimensions: nu=", mat_struct%nu, " na=", mat_struct%na

  ! Print first 5 rows
  DO I=1,MIN(5,mat_struct%nu)
    WRITE(*,'(A,I4,A)', advance='no') "Row ", I, ": "
    DO J=mat_struct%LdA(I), mat_struct%LdA(I+1)-1
      col = mat_struct%ColA(J)
      WRITE(*,'(A,I4,A,E12.4,A)', advance='no') "(", col, ",", coeffs(J), ") "
    END DO
    WRITE(*,*)
  END DO

  ! Check for NaN or Inf
  IF (ANY(ieee_is_nan(coeffs))) WRITE(*,*) "WARNING: NaN detected!"
  IF (ANY(.not. ieee_is_finite(coeffs))) WRITE(*,*) "WARNING: Inf detected!"
END SUBROUTINE

! Usage in Create_MMat:
CALL DebugPrintMatrix("Mass", mg_qMat(NLMAX), mg_Mmat(NLMAX)%a, NLMAX)
```

---

## Troubleshooting

### Common Errors

#### 1. Segmentation Fault in Assembly

**Symptom:** Crash with `Segmentation fault (core dumped)` during matrix assembly.

**Likely causes:**
```fortran
! ❌ Forgot to call SETLEV
DO ILEV=NLMIN,NLMAX
  qMat => mg_qMat(ILEV)
  CALL BuildMMat(..., NEL, ...)  ! NEL has wrong value! (from old ILEV)
END DO

! ✅ Correct version
DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)  ! Updates NEL for this level
  qMat => mg_qMat(ILEV)
  CALL BuildMMat(..., NEL, ...)
END DO
```

**Debugging:**
```fortran
WRITE(*,*) "ILEV=", ILEV, " NEL=", NEL, " qMat%na=", qMat%na
CALL FLUSH(6)
```

#### 2. Wrong Results After Refactoring

**Symptom:** Forces differ from baseline after moving routines.

**Check:**
1. **All USE statements present?**
   ```fortran
   USE QuadSc_struct, ONLY: Create_LinMatStruct  ! Forgot this?
   ```

2. **COMMON blocks still accessible?**
   - Check that `CBasicMG.inc`, `CBasicElem.inc` are included

3. **ILEV restoration:**
   ```fortran
   ! ❌ Forgot to restore
   DO ILEV=NLMIN,NLMAX
     ...
   END DO
   ! Now ILEV=NLMAX+1 → wrong!

   ! ✅ Always restore
   ILEV=NLMAX
   CALL SETLEV(2)
   ```

#### 3. Matrix Structure Not Allocated

**Symptom:** `Error: Trying to access unallocated pointer`

**Cause:** Structure allocation must happen before coefficient assembly.

**Fix:**
```fortran
! Initialization order matters!
CALL Create_QuadMatStruct()        ! First: structure
CALL Create_MMat()                 ! Then: coefficients
```

#### 4. Parallel Matrix Mismatch

**Symptom:** MPI deadlock or wrong results in parallel.

**Check:**
```fortran
! All processes must call assembly routines
IF (myid.ne.0) THEN  ! ❌ WRONG! Workers skip assembly
  CALL Create_BMat()
END IF

! ✅ Correct: everyone assembles (but master may skip kernel)
CALL Create_BMat()  ! Internal worker guards if needed
```

### Performance Tips

#### 1. Avoid Redundant Assembly

```fortran
! ❌ Reassembles every time step (slow!)
DO iTimeStep = 1, nTimeSteps
  CALL Create_MMat()  ! Mass matrix doesn't change!
  CALL Create_KMat(myScalar)
END DO

! ✅ Assemble once outside loop
CALL Create_MMat()
DO iTimeStep = 1, nTimeSteps
  CALL Create_KMat(myScalar)  ! Only update what changes
END DO
```

#### 2. Early Exit for Existing Matrices

```fortran
SUBROUTINE Create_hDiffMat()
  IF (ALLOCATED(mg_hDmat)) THEN
    IF (myid.eq.showID) WRITE(MTERM,'(A)') " [hD]: Exists |"
    RETURN  ! Already assembled, skip
  END IF
  ! ... assembly code ...
END SUBROUTINE
```

#### 3. Multigrid Efficiency

- Assemble all levels (NLMIN:NLMAX), not just NLMAX
- Multigrid solver needs coarse-level matrices for efficiency

---

## Summary Checklist

When adding/modifying matrix assembly:

- [ ] Declare matrix in `var_QuadScalar.f90`
- [ ] Add structure allocation to `QuadSc_struct.f90` (if new sparsity pattern)
- [ ] Add assembly routine to `QuadSc_assembly.f90`
- [ ] Export from module PUBLIC interface
- [ ] USE new routine in `def_QuadScalar.f90` (facade)
- [ ] Call `SETLEV(2)` before accessing COMMON variables
- [ ] Restore `ILEV=NLMAX; CALL SETLEV(2)` after multigrid loops
- [ ] Set NLMAX-level pointers for backward compatibility
- [ ] Test with baseline forces: `ForcesOnSurfInt` should match
- [ ] Document COMMON block usage in routine header

---

## Further Reading

- `def_refactoring.md` - Refactoring roadmap and dependency analysis
- `matrix_structures_guide.md` - Detailed CSR format and structure allocation
- `refactoring_roadmap.md` - Phase 1-4 refactoring plan
- FeatFloWer User Manual (if available)

**Questions?** Check `QuadSc_assembly.f90` for working examples of all matrix types.

---

**END OF GUIDE**
