# Matrix Structures and Allocation Guide

**Authors:** Claude Code (AI Assistant) + FeatFloWer Development Team
**Date:** 2025-12-05
**Version:** 1.0
**Purpose:** Guide for understanding matrix structures, allocation patterns, and memory management in the QuadScalar transport equation solver

---

## Table of Contents

1. [Overview](#overview)
2. [Fundamental Data Types](#fundamental-data-types)
3. [Matrix Naming Conventions](#matrix-naming-conventions)
4. [Mass Matrix Family](#mass-matrix-family)
5. [Diffusion Matrix Family](#diffusion-matrix-family)
6. [Two-Level Allocation Pattern](#two-level-allocation-pattern)
7. [Memory Management Best Practices](#memory-management-best-practices)
8. [Common Pitfalls and Debugging](#common-pitfalls-and-debugging)
9. [Code Examples](#code-examples)

---

## Overview

The FeatFloWer QuadScalar solver uses a sophisticated multigrid matrix structure to efficiently solve transport equations on hierarchical mesh levels. Understanding the allocation patterns and data types is crucial for:

- Implementing new matrix assembly routines
- Debugging segmentation faults and memory issues
- Extending the solver with new physics
- Refactoring and maintaining existing code

This document captures knowledge gained during Phase 1 and Phase 2 refactoring efforts (December 2025).

---

## Fundamental Data Types

The solver defines several key data types in `source/src_util/types.f90` and `source/src_quadLS/QuadSc_var.f90`:

### 1. `mg_Matrix` - Multigrid Matrix Container

**Definition:**
```fortran
TYPE mg_Matrix
  REAL*8, DIMENSION(:), ALLOCATABLE :: a
END TYPE mg_Matrix
```

**Purpose:** Stores matrix coefficients for a single multigrid level.

**Component:**
- `a`: 1D array of matrix entries in CSR (Compressed Sparse Row) format

**Usage:** All assembled matrices (mass, diffusion, stiffness, etc.)

**Example:**
```fortran
TYPE(mg_Matrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_MMat
! mg_MMat(ILEV)%a contains matrix values at level ILEV
```

---

### 2. `mg_dVector` - Multigrid Double Precision Vector

**Definition:**
```fortran
TYPE mg_dVector
  REAL*8, DIMENSION(:), ALLOCATABLE :: x
END TYPE mg_dVector
```

**Purpose:** Stores vector data (density, coefficients, solution fields) per multigrid level.

**Component:**
- `x`: 1D array of scalar values at mesh nodes

**Usage:** Density fields, diffusion coefficients, RHS vectors

**Example:**
```fortran
TYPE(mg_dVector), DIMENSION(:), ALLOCATABLE :: mgDensity
! mgDensity(ILEV)%x contains density values at level ILEV
```

---

### 3. `mg_kVector` - Multigrid Integer Vector

**Definition:**
```fortran
TYPE mg_kVector
  INTEGER, DIMENSION(:), ALLOCATABLE :: x
END TYPE mg_kVector
```

**Purpose:** Stores integer data (boundary markers, DOF constraints) per level.

**Component:**
- `x`: 1D array of integer flags/indices

**Usage:** Pressure DOF constraints (`knprP`), boundary condition flags

---

### 4. `TMatrix` - Matrix Structure Descriptor

**Key Components:**
```fortran
TYPE TMatrix
  INTEGER :: nu      ! Number of rows (unknowns)
  INTEGER :: na      ! Number of non-zero entries
  INTEGER, DIMENSION(:), ALLOCATABLE :: LdA   ! Row pointer (nu+1)
  INTEGER, DIMENSION(:), ALLOCATABLE :: ColA  ! Column indices (na)
END TYPE TMatrix
```

**Purpose:** Describes CSR matrix sparsity pattern (structure only, no values).

**Usage:** Paired with coefficient arrays (e.g., `lMat` structure + `CMat` values)

---

## Matrix Naming Conventions

Understanding the naming scheme is critical for navigating the codebase:

### Prefix Conventions

| Prefix | Meaning | Example | Description |
|--------|---------|---------|-------------|
| `mg_`  | Multigrid array | `mg_MMat` | Array indexed by level (NLMIN:NLMAX) |
| `l`    | Local/Linear | `lMat` | Local to process, linear element |
| `q`    | Quadratic | `qMat` | Quadratic finite element |
| `C`    | Coefficients | `CMat` | Actual matrix values |
| `P`    | Parallel | `lPMat` | Off-process coupling entries |

### Suffix Conventions

| Suffix | Meaning | Example | Description |
|--------|---------|---------|-------------|
| `Mat`  | Full matrix | `MMat` | Assembled matrix |
| `lMat` | Lumped matrix | `MlMat` | Row-sum diagonal |
| `PMat` | Parallel-synchronized | `MlPMat` | After MPI communication |

### Special Naming Patterns

- **`mg_Mmat`**: Full mass matrix (per level)
- **`mg_MlMat`**: Lumped mass matrix (diagonal, per level)
- **`mg_MlPMat`**: Parallel-synchronized lumped mass matrix
- **`mg_MlRhomat`**: Density-weighted lumped mass matrix
- **`Mmat` (no prefix)**: Pointer to current level's mass matrix

---

## Mass Matrix Family

The mass matrix family is central to time-dependent transport equations. Understanding their relationships is essential:

### 1. Full Mass Matrix (`mg_MMat`)

**Mathematical Form:**
```
M[i,j] = ∫_Ω φ_i(x) φ_j(x) dx
```

**Storage:** CSR format, sparse matrix

**Assembly:**
- `BuildMMat()`: Standard mass matrix (constant density ρ = 1)
- `BuildMRhoMat(DENS, ...)`: Density-weighted mass matrix

**Allocation:**
```fortran
! Two-level allocation (see below)
ALLOCATE(mg_MMat(NLMIN:NLMAX))              ! Level array
ALLOCATE(mg_MMat(ILEV)%a(qMat%na))          ! Per-level coefficients
```

**Usage:** Time discretization, velocity transport

---

### 2. Lumped Mass Matrix (`mg_MlMat`, `mg_MlRhomat`)

**Mathematical Form:**
```
M_lumped[i] = Σ_j M[i,j]  (row sum)
```

**Storage:** 1D array (diagonal matrix stored as vector)

**Computation:**
```fortran
DO I = 1, qMat%nu
  DML = 0d0
  DO J = qMat%LdA(I), qMat%LdA(I+1)-1
    DML = DML + mg_MMat(ILEV)%a(J)
  END DO
  mg_MlMat(ILEV)%a(I) = DML
END DO
```

**Purpose:** Mass-lumping for explicit time integration, efficient diagonal inversion

**Variants:**
- **`mg_MlMat`**: From standard mass matrix
- **`mg_MlRhomat`**: From density-weighted mass matrix

---

### 3. Parallel-Synchronized Lumped Mass (`mg_MlPMat`, `mg_MlRhoPmat`)

**Purpose:** Accumulates contributions from all MPI processes sharing DOFs

**Synchronization:**
```fortran
mg_MlPMat(ILEV)%a = mg_MlMat(ILEV)%a
CALL E013SUM(mg_MlPMat(ILEV)%a)  ! MPI reduction
```

**Why needed:** Domain decomposition creates duplicate DOFs at subdomain boundaries. The MPI sum ensures global consistency.

**Usage:** Global mass conservation, parallel preconditioners

---

### Mass Matrix Relationships

```
┌─────────────────┐
│  mgDensity(x)   │  (input: density field)
└────────┬────────┘
         │
         ▼
┌─────────────────────────────┐
│  BuildMRhoMat / BuildMMat   │  (assembly)
└────────┬────────────────────┘
         │
         ▼
┌─────────────────┐
│   mg_MMat(lev)  │  (full CSR matrix)
└────────┬────────┘
         │
         ├───── (row sum) ──────► mg_MlMat(lev)  ──── (E013SUM) ──►  mg_MlPMat(lev)
         │                        (local lumped)                     (parallel lumped)
         │
         └───── (row sum) ──────► mg_MlRhomat(lev) ── (E013SUM) ──►  mg_MlRhoPmat(lev)
                                   (density lumped)                   (parallel density lumped)
```

---

## Diffusion Matrix Family

The diffusion matrix family represents spatial diffusion operators (viscous/thermal transport) in the transport equations. Understanding their structure and assembly is crucial for:

- Implementing viscosity models (Newtonian, non-Newtonian)
- Temperature-dependent diffusion
- Multi-material simulations
- Phase 2.2 refactoring insights

### 1. Standard Diffusion Matrix (`mg_hDmat`)

**Mathematical Form:**
```
D[i,j] = ∫_Ω ∇φ_i(x) · ∇φ_j(x) dx
```

**Physical Meaning:** Standard Laplacian operator (heat diffusion, constant viscosity)

**Assembly:**
- Kernel: `DIFFQ2_alpha(..., alpha=1.0)`
- Alpha parameter: 1.0 (standard weighting)
- Process guard: Only `myid.ne.0` (worker processes only)

**Allocation:**
```fortran
! Two-level allocation
ALLOCATE(mg_hDmat(NLMIN:NLMAX))              ! Top-level
ALLOCATE(mg_hDmat(ILEV)%a(qMat%na))          ! Per-level coefficients
```

**Usage:**
- Thermal diffusion with constant conductivity
- Momentum diffusion with constant viscosity
- Baseline diffusion operator for multigrid smoothers

**Early-Exit Behavior:**
```fortran
IF (ALLOCATED(mg_hDmat)) THEN
  ! Matrix already exists, skip reassembly
  WRITE(MTERM,'(A)') " [hD]: Exists |"
  RETURN
END IF
```

**Why "hDmat"?** The "h" likely refers to "homogeneous" or "heat" diffusion (standard Laplacian).

---

### 2. Constant Viscous Diffusion Matrix (`mg_ConstDMat`)

**Mathematical Form:**
```
D[i,j] = ∫_Ω ν ∇φ_i(x) · ∇φ_j(x) dx  (where ν = constant)
```

**Physical Meaning:** Viscous diffusion with constant viscosity coefficient

**Assembly:**
- Kernel: `DIFFQ2_alpha(..., alpha=0.0)`
- Alpha parameter: 0.0 (constant viscosity weighting)
- Process guard: None (all processes assemble)

**Allocation:**
```fortran
! Two-level allocation
ALLOCATE(mg_ConstDMat(NLMIN:NLMAX))          ! Top-level
ALLOCATE(mg_ConstDMat(ILEV)%a(qMat%na))      ! Per-level coefficients
```

**Usage:**
- Newtonian fluid flow (constant dynamic viscosity)
- Constant thermal conductivity problems
- Reference operator for viscous problems

**Early-Exit Behavior:**
```fortran
IF (ALLOCATED(mg_ConstDMat)) THEN
  ! Matrix already exists, skip reassembly
  WRITE(MTERM,'(A)') " [VD]: Exists |"
  RETURN
END IF
```

**Why "ConstDMat"?** "Const" = constant viscosity, "D" = diffusion operator

**Label "[VD]":** "Viscous Diffusion" in progress messages

---

### 3. General Diffusion Matrix (`mg_Dmat`)

**Mathematical Form:**
```
D[i,j] = ∫_Ω ν(u, ∇u, T, ...) ∇φ_i(x) · ∇φ_j(x) dx
```

**Physical Meaning:** Variable-viscosity diffusion (velocity-dependent, temperature-dependent, or material-dependent)

**Assembly:** Complex, multiple kernels based on physics:

**Newtonian Case (`bNonNewtonian = .FALSE.`):**
```fortran
CALL DIFFQ2_NEWT(mg_Dmat(ILEV)%a, ...)
```
- Standard Newtonian viscosity
- No velocity dependence

**Non-Newtonian Case (`bNonNewtonian = .TRUE.`):**

**Single Material:**
```fortran
CALL DIFFQ2_NNEWT(myScalar%valU, myScalar%valV, myScalar%valW, &
                  Temperature, mg_Dmat(ILEV)%a, ...)
```
- Velocity-dependent viscosity (shear-thinning, shear-thickening)
- Temperature-dependent viscosity
- Generalized Newtonian models (power-law, Carreau, etc.)

**Multi-Material (`bMultiMat = .TRUE.`):**
```fortran
CALL DIFFQ2_AlphaNNEWT(myScalar%valU, myScalar%valV, myScalar%valW, &
                       MaterialDistribution(ILEV)%x, &
                       mg_Dmat(ILEV)%a, ...)
```
- Material-dependent viscosity
- Interface tracking (multi-phase flows)
- Heterogeneous material properties

**Allocation:**
```fortran
! Two-level allocation
IF (.NOT. ALLOCATED(mg_Dmat)) ALLOCATE(mg_Dmat(NLMIN:NLMAX))
ALLOCATE(mg_Dmat(ILEV)%a(qMat%na))
```

**Usage:**
- Non-Newtonian fluid flows (polymer melts, blood, etc.)
- Temperature-dependent viscosity (thermal coupling)
- Multi-material simulations (solid-fluid, fluid-fluid)
- Viscoelastic flows

**Why Not Refactored in Phase 2.2?**
- Complex branching logic (bNonNewtonian, bMultiMat)
- Three different assembly kernels
- Velocity field dependencies (`myScalar%valU/V/W`)
- Temperature/material distribution dependencies
- Would require sophisticated generic routine
- Left unchanged for clarity and stability

---

### Diffusion Matrix Comparison Table

| Property | `mg_hDmat` | `mg_ConstDMat` | `mg_Dmat` |
|----------|------------|----------------|-----------|
| **Alpha Parameter** | 1.0 | 0.0 | N/A (multiple kernels) |
| **Assembly Kernel** | `DIFFQ2_alpha` | `DIFFQ2_alpha` | `DIFFQ2_NEWT` / `DIFFQ2_NNEWT` / `DIFFQ2_AlphaNNEWT` |
| **Process Guard** | `myid.ne.0` only | All processes | All processes |
| **Early Exit** | Yes (if exists) | Yes (if exists) | No |
| **Velocity Dependent** | No | No | Yes (non-Newtonian) |
| **Material Dependent** | No | No | Yes (multi-material) |
| **Phase 2.2 Refactored** | ✓ Yes | ✓ Yes | ✗ No (too complex) |
| **Label** | `[hD]` | `[VD]` | `[D]` |
| **Typical Use** | Baseline operator | Constant viscosity | Variable viscosity |

---

### Assembly Kernel: `DIFFQ2_alpha`

**Signature:**
```fortran
SUBROUTINE DIFFQ2_alpha(DA, NA, KCOLA, KLDA, KVERT, KAREA, KEDGE, &
                        DCORVG, ELE, alpha)
  REAL*8 :: DA(NA)              ! Output matrix coefficients
  INTEGER :: NA                 ! Number of non-zeros
  INTEGER :: KCOLA(NA)          ! Column indices (CSR)
  INTEGER :: KLDA(NU+1)         ! Row pointers (CSR)
  INTEGER :: KVERT(8,NEL)       ! Element-vertex connectivity
  INTEGER :: KAREA(6,NEL)       ! Element-face connectivity
  INTEGER :: KEDGE(12,NEL)      ! Element-edge connectivity
  REAL*8 :: DCORVG(3,NVT)       ! Vertex coordinates
  EXTERNAL :: ELE               ! Finite element routine
  REAL*8 :: alpha               ! Alpha weighting parameter
END SUBROUTINE
```

**Purpose:** Assembles diffusion matrix with alpha-weighted quadrature

**Alpha Parameter Meaning:**
- `alpha = 1.0`: Standard diffusion (hDmat)
- `alpha = 0.0`: Constant viscous diffusion (ConstDMat)
- Other values: User-defined weighting (not used in current code)

**Implementation Location:** `source/assemblies/` (exact file TBD - would need search)

---

### Phase 2.2 Refactoring: `Assemble_Diffusion_Alpha_Generic`

**Purpose:** Unified assembly for `Create_hDiffMat` and `Create_ConstDiffMat`

**Signature:**
```fortran
SUBROUTINE Assemble_Diffusion_Alpha_Generic(mg_OutMatrix, alpha, label, &
                                             require_worker)
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_OutMatrix(NLMIN:NLMAX)
  REAL*8, INTENT(IN) :: alpha
  CHARACTER(LEN=*), INTENT(IN) :: label
  LOGICAL, INTENT(IN) :: require_worker
```

**Parameters:**
- `mg_OutMatrix`: Output diffusion matrix (must be allocated by caller)
- `alpha`: Parameter for `DIFFQ2_alpha` kernel
- `label`: Progress message label (`"[hD]"` or `"[VD]"`)
- `require_worker`: If `.TRUE.`, only `myid.ne.0` assembles

**Refactored Routines:**

**`Create_hDiffMat` (Before: 57 lines, After: 18 lines):**
```fortran
! Caller handles allocation and early-exit
IF (.NOT. ALLOCATED(mg_hDmat)) THEN
  ALLOCATE(mg_hDmat(NLMIN:NLMAX))
ELSE
  RETURN  ! Already exists
END IF

! Call generic with alpha=1.0, worker-only
CALL Assemble_Diffusion_Alpha_Generic(mg_hDmat, 1.0d0, "[hD]", .TRUE.)
```

**`Create_ConstDiffMat` (Before: 55 lines, After: 18 lines):**
```fortran
! Caller handles allocation and early-exit
IF (.NOT. ALLOCATED(mg_ConstDMat)) THEN
  ALLOCATE(mg_ConstDMat(NLMIN:NLMAX))
ELSE
  RETURN  ! Already exists
END IF

! Call generic with alpha=0.0, all processes
CALL Assemble_Diffusion_Alpha_Generic(mg_ConstDMat, 0.0d0, "[VD]", .FALSE.)
```

**Code Reduction:**
- `Create_hDiffMat`: 68% reduction (57 → 18 lines)
- `Create_ConstDiffMat`: 67% reduction (55 → 18 lines)
- Total duplication eliminated: ~94 lines

**Design Benefits:**
- Single source of truth for alpha-based diffusion
- Easy to add new alpha variants in future
- Consistent timing, progress messages, allocation patterns
- Clear separation: simple (alpha-based) vs complex (`mg_Dmat`)

---

### Diffusion Matrix Relationships

```
                    ┌──────────────────────────┐
                    │  Physics Requirements    │
                    └────────────┬─────────────┘
                                 │
                ┌────────────────┼────────────────┐
                │                │                │
                ▼                ▼                ▼
        ┌───────────────┐  ┌──────────────┐  ┌────────────────┐
        │  Constant     │  │  Constant    │  │  Variable      │
        │  Diffusion    │  │  Viscosity   │  │  Viscosity     │
        └───────┬───────┘  └──────┬───────┘  └────────┬───────┘
                │                 │                    │
                ▼                 ▼                    ▼
        ┌───────────────┐  ┌──────────────┐  ┌────────────────┐
        │  mg_hDmat     │  │ mg_ConstDMat │  │   mg_Dmat      │
        │  (alpha=1.0)  │  │  (alpha=0.0) │  │  (NEWT/NNEWT)  │
        └───────┬───────┘  └──────┬───────┘  └────────┬───────┘
                │                 │                    │
                └─────────┬───────┴────────────────────┘
                          ▼
                ┌──────────────────────────┐
                │  Multigrid Hierarchy     │
                │  NLMIN → NLMAX           │
                │  Each level: CSR matrix  │
                └──────────────────────────┘
```

---

### When to Use Which Diffusion Matrix?

**Use `mg_hDmat` when:**
- Implementing baseline multigrid smoothers
- Constant thermal conductivity problems
- Need standard Laplacian operator
- Debugging viscous solvers (simplest case)

**Use `mg_ConstDMat` when:**
- Newtonian fluid flow (constant dynamic viscosity)
- Constant material properties
- Benchmark problems with known viscosity
- Reference solutions for validation

**Use `mg_Dmat` when:**
- Non-Newtonian flows (shear-thinning/thickening)
- Temperature-dependent viscosity
- Multi-material simulations
- Polymer processing (extrusion, injection molding)
- Blood flow simulations
- Any velocity-dependent constitutive law

---

### Common Usage Patterns

**Pattern 1: Matrix Renewal (Newtonian Flow)**
```fortran
! Called once at initialization
CALL Create_ConstDiffMat()  ! Allocates and assembles

! During time-stepping: reuse unless mesh changes
! (ConstDMat doesn't change if viscosity is constant)
```

**Pattern 2: Matrix Renewal (Non-Newtonian Flow)**
```fortran
! Called at initialization
CALL Create_DiffMat(QuadSc)  ! Initial assembly

! During time-stepping: reassemble when needed
IF (myMatrixRenewal%D == 1) THEN
  ! Velocity has changed significantly
  CALL Create_DiffMat(QuadSc)  ! Reassemble with new velocity
END IF
```

**Pattern 3: Mixed Usage**
```fortran
! Use constant diffusion for preconditioner
CALL Create_ConstDiffMat()

! Use variable diffusion for exact operator
CALL Create_DiffMat(QuadSc)

! Multigrid: smooth with ConstDMat, correct with Dmat
```

---

### Memory Footprint

For a typical 3D problem with `NEL` elements and Q2 finite elements:

**Number of unknowns per level:**
```
NU = NVT + NET + NAT + NEL
   ≈ (2^lev + 1)^3  (vertices)
   + 3*(2^lev)*(2^lev + 1)^2  (edges)
   + 3*(2^lev)^2*(2^lev + 1)  (faces)
   + (2^lev)^3  (cells)
```

**Non-zero entries (CSR format):**
```
NA ≈ 27 * NU  (Q2 stencil: 27-point in 3D)
```

**Memory per diffusion matrix:**
```
Memory ≈ 8 bytes/entry * 27 * NU
       ≈ 216 * NU bytes
```

**All three diffusion matrices:**
```
Total ≈ 3 * 216 * NU bytes
      = 648 * NU bytes
```

**Example (NLMAX with NU ≈ 1M unknowns):**
```
Per matrix: 216 MB
All three:  648 MB
```

**Optimization tip:** If using only one type, comment out calls to others to save memory.

---

### Debugging Diffusion Matrices

**Check if assembled:**
```fortran
IF (ALLOCATED(mg_hDmat)) THEN
  WRITE(*,*) 'hDmat allocated on', LBOUND(mg_hDmat), 'to', UBOUND(mg_hDmat)
  DO ILEV = NLMIN, NLMAX
    IF (ALLOCATED(mg_hDmat(ILEV)%a)) THEN
      WRITE(*,*) '  Level', ILEV, ': size =', SIZE(mg_hDmat(ILEV)%a)
    END IF
  END DO
END IF
```

**Verify symmetry (should be symmetric for diffusion):**
```fortran
! Check if D[i,j] ≈ D[j,i]
symmetric_error = 0d0
DO IEQ = 1, qMat%nu
  DO IA = qMat%LdA(IEQ), qMat%LdA(IEQ+1)-1
    JCOL = qMat%ColA(IA)
    ! Find D[j,i] and compare with D[i,j]
    ! (implementation left as exercise)
  END DO
END DO
```

**Check diagonal dominance:**
```fortran
! Diffusion matrices should be diagonally dominant
DO IEQ = 1, qMat%nu
  diag_val = mg_hDmat(ILEV)%a(qMat%LdA(IEQ))  ! Diagonal entry
  off_diag_sum = 0d0
  DO IA = qMat%LdA(IEQ)+1, qMat%LdA(IEQ+1)-1
    off_diag_sum = off_diag_sum + ABS(mg_hDmat(ILEV)%a(IA))
  END DO
  IF (diag_val < off_diag_sum) THEN
    WRITE(*,*) 'WARNING: Row', IEQ, 'not diagonally dominant!'
  END IF
END DO
```

---

## Two-Level Allocation Pattern

**Critical concept:** Multigrid arrays require **two allocation steps**:

### Level 1: Top-Level Array Structure

**Purpose:** Allocate the array of containers (one per multigrid level)

**Responsibility:** **Caller function** (e.g., `Create_MMat`, `Create_MRhoMat`)

**Example:**
```fortran
IF (.NOT. ALLOCATED(mg_MMat)) ALLOCATE(mg_MMat(NLMIN:NLMAX))
```

**Memory Layout:**
```
mg_MMat(NLMIN:NLMAX)  →  [mg_Matrix, mg_Matrix, ..., mg_Matrix]
                          ↓          ↓                ↓
                          NLMIN      NLMIN+1          NLMAX
```

---

### Level 2: Per-Level Component Allocation

**Purpose:** Allocate the actual data arrays within each level's container

**Responsibility:** **Assembly routine** or **generic helper**

**Example:**
```fortran
DO ILEV = NLMIN, NLMAX
  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  ! Allocate coefficient array for this level
  IF (.NOT. ALLOCATED(mg_MMat(ILEV)%a)) THEN
    ALLOCATE(mg_MMat(ILEV)%a(qMat%na))
  END IF
END DO
```

**Memory Layout:**
```
mg_MMat(ILEV)%a  →  [coeff₁, coeff₂, ..., coeff_na]
                     (CSR format: qMat%na entries)
```

---

### Why Two Levels?

1. **Multigrid hierarchy:** Different levels have different mesh sizes (`qMat%na` varies)
2. **Fortran allocatable arrays:** Cannot directly allocate `mg_MMat(NLMIN:NLMAX)%a(na)` in one step
3. **Flexibility:** Per-level allocation allows dynamic sizing based on matrix structure

---

### Common Allocation Mistake

**WRONG (causes segfault):**
```fortran
SUBROUTINE Create_MMat()
  ! Missing top-level allocation!
  DO ILEV = NLMIN, NLMAX
    ALLOCATE(mg_MMat(ILEV)%a(qMat%na))  ! SEGFAULT: mg_MMat not allocated
  END DO
END SUBROUTINE
```

**CORRECT:**
```fortran
SUBROUTINE Create_MMat()
  ! Step 1: Allocate top-level array
  IF (.NOT. ALLOCATED(mg_MMat)) ALLOCATE(mg_MMat(NLMIN:NLMAX))

  ! Step 2: Allocate per-level components
  DO ILEV = NLMIN, NLMAX
    IF (.NOT. ALLOCATED(mg_MMat(ILEV)%a)) THEN
      ALLOCATE(mg_MMat(ILEV)%a(qMat%na))
    END IF
  END DO
END SUBROUTINE
```

---

## Memory Management Best Practices

### 1. Allocation Responsibility Separation

When creating generic/helper routines, clearly define allocation responsibilities:

**Example: Phase 2.1 Generic Mass Matrix Assembly**

```fortran
! Caller (Create_MMat): Allocates top-level arrays
ALLOCATE(mg_MlMat(NLMIN:NLMAX))
ALLOCATE(mg_MlPMat(NLMIN:NLMAX))

! Generic routine (Assemble_Mass_Generic): Allocates per-level components
DO ILEV = NLMIN, NLMAX
  IF (.NOT. ALLOCATED(mg_MlMat(ILEV)%a)) THEN
    ALLOCATE(mg_MlMat(ILEV)%a(qMat%nu))
  END IF
END DO
```

**Rationale:**
- Caller knows the full context (which arrays to allocate)
- Generic routine handles uniform per-level allocation
- Clear separation prevents duplicate allocation checks

---

### 2. Parameter Intent Declarations

Use correct `INTENT` for subroutine parameters:

```fortran
SUBROUTINE Assemble_Mass_Generic(use_density, mg_MlMatrix, mg_MlPMatrix, &
                                  density_opt, label)
  LOGICAL, INTENT(IN) :: use_density
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_MlMatrix(NLMIN:NLMAX)   ! Will be modified
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_MlPMatrix(NLMIN:NLMAX)  ! Will be modified
  TYPE(mg_dVector), INTENT(IN), OPTIONAL :: density_opt(NLMIN:NLMAX)  ! Read-only
  CHARACTER(LEN=*), INTENT(IN) :: label
```

**Never check `ALLOCATED()` on dummy arguments:**
```fortran
! WRONG - causes compiler error
IF (.NOT. ALLOCATED(mg_MlMatrix)) ALLOCATE(mg_MlMatrix(NLMIN:NLMAX))

! CORRECT - check per-level components
IF (.NOT. ALLOCATED(mg_MlMatrix(ILEV)%a)) ALLOCATE(mg_MlMatrix(ILEV)%a(nu))
```

---

### 3. Deallocation Before Reallocation

Avoid memory leaks when regenerating matrices:

```fortran
! Safe reallocation pattern
IF (ALLOCATED(mg_MMat(ILEV)%a)) DEALLOCATE(mg_MMat(ILEV)%a)
ALLOCATE(mg_MMat(ILEV)%a(qMat%na))
```

**When to deallocate:**
- Matrix structure changes (mesh refinement)
- Solver reconfiguration
- Memory optimization (freeing unused levels)

---

### 4. Nullifying Pointers After Use

Prevent dangling pointers:

```fortran
! After using pointers
Mmat      => mg_MMat(NLMAX)%a
MlMat     => mg_MlMat(NLMAX)%a

! At end of routine (optional but safe)
NULLIFY(Mmat, MlMat)
```

---

## Common Pitfalls and Debugging

### 1. Segmentation Fault: Unallocated Top-Level Array

**Symptom:**
```
Caught signal 11 (Segmentation fault: address not mapped to object at address (nil))
Backtrace: at QuadSc_def.f90:345
```

**Cause:** Accessing `mg_MMat(ILEV)%a` when `mg_MMat` itself is not allocated.

**Debug Strategy:**
```fortran
! Add diagnostic check
IF (.NOT. ALLOCATED(mg_MMat)) THEN
  WRITE(*,*) 'ERROR: mg_MMat not allocated before use!'
  STOP
END IF
```

**Fix:** Ensure caller allocates top-level array before passing to subroutines.

---

### 2. Type Mismatch: `mg_Matrix` vs `mg_dVector`

**Symptom:**
```
Error: 'x' at (1) is not a member of the 'mg_matrix' structure
```

**Cause:** Wrong type used in function signature.

**Correct Usage:**
- **Matrices:** Use `TYPE(mg_Matrix)` → access via `%a`
- **Vectors (density, RHS):** Use `TYPE(mg_dVector)` → access via `%x`
- **Integer flags:** Use `TYPE(mg_kVector)` → access via `%x`

**Example Fix:**
```fortran
! WRONG
TYPE(mg_Matrix) :: density(NLMIN:NLMAX)
val = density(ILEV)%x  ! ERROR: mg_Matrix has no %x

! CORRECT
TYPE(mg_dVector) :: density(NLMIN:NLMAX)
val = density(ILEV)%x  ! OK
```

---

### 3. MPI Synchronization Issues

**Symptom:** Results differ between sequential and parallel runs, or non-deterministic errors.

**Cause:** Forgetting to synchronize shared DOF contributions.

**Solution:** Always use parallel-synchronized arrays for global operations:

```fortran
! Local computation
mg_MlMat(ILEV)%a = ... (row sums)

! Synchronize for global consistency
mg_MlPMat(ILEV)%a = mg_MlMat(ILEV)%a
CALL E013SUM(mg_MlPMat(ILEV)%a)  ! MPI_Allreduce equivalent

! Use parallel version for global operations
global_mass = SUM(mg_MlPMat(ILEV)%a)  ! Correct
! global_mass = SUM(mg_MlMat(ILEV)%a)  ! WRONG (local only)
```

---

### 4. SETLEV Side Effects

**Issue:** Legacy COMMON block global variables (see Phase 1 HYPRE bug documentation).

**Affected Variables:** `ILEV`, `NEL`, `NVT`, etc.

**Rule:** Call `SETLEV(2)` on **all MPI processes** before MPI collective operations, even if only worker processes (myid ≠ 0) perform matrix work.

**Example:**
```fortran
! CORRECT: All processes call SETLEV
ILEV = lScalar%prm%MGprmIn%MinLev
CALL SETLEV(2)  ! Updates global NEL, etc.

IF (myid .NE. 0) THEN
  ! Only workers assemble matrices
  CALL BuildMatrix(...)
END IF

! MPI collective depends on NEL being set on ALL processes
CALL GetMyNumberingLimits(ilower, iupper, NEL)  ! Uses global NEL
```

---

## Code Examples

### Example 1: Creating a New Matrix Type

```fortran
! Step 1: Declare module-level variable (QuadSc_var.f90)
TYPE(mg_Matrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_MyNewMat

! Step 2: Create assembly routine (QuadSc_def.f90)
SUBROUTINE Create_MyNewMat()
  EXTERNAL E013
  INTEGER :: I, J

  ! Early exit if not needed
  IF (.NOT. bMasterTurnedOn) RETURN

  ! Step 2.1: Allocate top-level array
  IF (.NOT. ALLOCATED(mg_MyNewMat)) ALLOCATE(mg_MyNewMat(NLMIN:NLMAX))

  ! Step 2.2: Loop over multigrid levels
  DO ILEV = NLMIN, NLMAX
    CALL SETLEV(2)
    qMat => mg_qMat(ILEV)

    ! Step 2.3: Allocate per-level matrix
    IF (.NOT. ALLOCATED(mg_MyNewMat(ILEV)%a)) THEN
      ALLOCATE(mg_MyNewMat(ILEV)%a(qMat%na))
    END IF

    mg_MyNewMat(ILEV)%a = 0d0

    ! Step 2.4: Call assembly kernel
    CALL BuildMyNewMat(mg_MyNewMat(ILEV)%a, qMat%na, &
                       qMat%ColA, qMat%LdA, &
                       mg_mesh%level(ILEV)%kvert, &
                       mg_mesh%level(ILEV)%dcorvg, E013)
  END DO

  ! Step 2.5: Set pointer to finest level
  ILEV = NLMAX
  CALL SETLEV(2)
  MyNewMat => mg_MyNewMat(NLMAX)%a

END SUBROUTINE Create_MyNewMat
```

---

### Example 2: Generic Routine with Optional Parameters

```fortran
SUBROUTINE Assemble_Generic(use_feature, mg_Output, feature_data)
  LOGICAL, INTENT(IN) :: use_feature
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_Output(NLMIN:NLMAX)
  TYPE(mg_dVector), INTENT(IN), OPTIONAL :: feature_data(NLMIN:NLMAX)

  ! Validate optional parameter
  IF (use_feature) THEN
    IF (.NOT. PRESENT(feature_data)) THEN
      WRITE(*,*) 'ERROR: use_feature=.TRUE. requires feature_data!'
      STOP
    END IF
  END IF

  ! Allocate and assemble
  DO ILEV = NLMIN, NLMAX
    CALL SETLEV(2)
    qMat => mg_qMat(ILEV)

    IF (.NOT. ALLOCATED(mg_Output(ILEV)%a)) THEN
      ALLOCATE(mg_Output(ILEV)%a(qMat%na))
    END IF

    IF (use_feature) THEN
      CALL BuildWithFeature(mg_Output(ILEV)%a, feature_data(ILEV)%x, ...)
    ELSE
      CALL BuildStandard(mg_Output(ILEV)%a, ...)
    END IF
  END DO

END SUBROUTINE Assemble_Generic
```

---

### Example 3: Debugging Allocation Issues

```fortran
SUBROUTINE Debug_MatrixAllocation()
  INTEGER :: ILEV

  WRITE(*,*) '===== Matrix Allocation Diagnostics ====='

  ! Check top-level allocation
  WRITE(*,*) 'mg_MMat allocated:', ALLOCATED(mg_MMat)
  IF (ALLOCATED(mg_MMat)) THEN
    WRITE(*,*) 'mg_MMat bounds:', LBOUND(mg_MMat), UBOUND(mg_MMat)
  END IF

  ! Check per-level allocation
  IF (ALLOCATED(mg_MMat)) THEN
    DO ILEV = NLMIN, NLMAX
      WRITE(*,'(A,I2,A,L1)') ' Level ', ILEV, ' allocated: ', &
                              ALLOCATED(mg_MMat(ILEV)%a)
      IF (ALLOCATED(mg_MMat(ILEV)%a)) THEN
        WRITE(*,'(A,I10)') '   Size: ', SIZE(mg_MMat(ILEV)%a)
      END IF
    END DO
  END IF

  WRITE(*,*) '=========================================='

END SUBROUTINE Debug_MatrixAllocation
```

---

## Summary

**Key Takeaways:**

1. **Understand the data types:**
   - `mg_Matrix` (component: `%a`) for matrices
   - `mg_dVector` (component: `%x`) for scalar fields
   - `mg_kVector` (component: `%x`) for integer flags

2. **Two-level allocation is mandatory:**
   - Caller allocates: `ALLOCATE(mg_MMat(NLMIN:NLMAX))`
   - Assembly allocates: `ALLOCATE(mg_MMat(ILEV)%a(na))`

3. **Mass matrix family:**
   - Full (`mg_MMat`) → Lumped (`mg_MlMat`) → Parallel (`mg_MlPMat`)
   - Density variant: `mg_MlRhomat`, `mg_MlRhoPmat`

4. **Naming conventions encode structure:**
   - `mg_` = multigrid array
   - `l` suffix = lumped (diagonal)
   - `P` suffix = parallel-synchronized

5. **Common pitfalls:**
   - Forgetting top-level allocation → segfault
   - Type mismatch (`mg_Matrix` vs `mg_dVector`)
   - Missing MPI synchronization → wrong results
   - SETLEV scope issues (see Phase 1 HYPRE bug)

6. **Debugging checklist:**
   - Check `ALLOCATED()` status at both levels
   - Verify correct type and component access (`%a` vs `%x`)
   - Ensure SETLEV called on all processes before MPI collectives
   - Use diagnostic output to trace allocation flow

---

## References

- **Phase 1 Documentation:** `docs/md_docs/phase1_extraction_complete.md` (HYPRE/UMFPACK solver extraction, SETLEV pitfalls)
- **Phase 2 Roadmap:** `docs/md_docs/refactoring_roadmap.md` (Deduplication strategy)
- **Source Files:**
  - Type definitions: `source/src_util/types.f90`
  - Matrix variables: `source/src_quadLS/QuadSc_var.f90`
  - Mass matrix routines: `source/src_quadLS/QuadSc_def.f90` (lines 263-550)
  - Assembly kernels: `source/assemblies/QuadSc_massrho.f`

---

**Document History:**
- 2025-12-05: Initial version (Phase 2.1 mass matrix refactoring insights)
- 2025-12-05: Added Diffusion Matrix Family section (Phase 2.2 refactoring insights)

**Contributing:**
If you discover additional allocation patterns or encounter new pitfalls, please update this document and commit with a descriptive message.
