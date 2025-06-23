# FeatFloWer Mesh Structure Documentation

This document provides comprehensive documentation of the mesh organization, data structures, connectivity arrays, and access patterns used in the FeatFloWer computational fluid dynamics framework.

---

## 1. Overview

FeatFloWer uses a sophisticated hierarchical mesh structure designed for multigrid finite element methods with 3D hexahedral elements. The mesh supports:

- **Hexahedral elements**: 8-node brick elements for volume discretization
- **Hierarchical refinement**: Multi-level grids for multigrid solvers
- **Mixed finite elements**: Q2/P1 discretization (27-node velocity, 4-node pressure)
- **Complex connectivity**: Comprehensive element-vertex-edge-face relationships
- **Parallel decomposition**: Domain decomposition for MPI parallelization

---

## 2. Mesh Data Structures

### 2.1 Primary Mesh Type (`tMesh`)

**Location**: `source/src_util/types.f90`

```fortran
TYPE tMesh
  ! Mesh dimensions
  integer :: NEL = 0    ! Number of elements
  integer :: NVT = 0    ! Number of vertices
  integer :: NET = 0    ! Number of edges
  integer :: NAT = 0    ! Number of faces (areas)
  integer :: NBCT = 0   ! Number of boundary conditions
  integer :: NVBD = 0   ! Number of boundary vertices
  integer :: NEBD = 0   ! Number of boundary edges
  integer :: NABD = 0   ! Number of boundary faces

  ! Element topology constants
  integer :: NVE = 8    ! Vertices per element (hexahedral)
  integer :: NEE = 12   ! Edges per element
  integer :: NAE = 6    ! Faces per element

  ! Coordinate arrays
  real*8,  pointer, dimension(:,:) :: dcorvg => null()  ! Vertex coordinates
  real*8,  pointer, dimension(:,:) :: dcorag => null()  ! Face center coordinates

  ! Primary connectivity arrays
  integer, allocatable, dimension(:,:) :: kvert         ! Element-vertex connectivity
  integer, allocatable, dimension(:,:) :: kedge         ! Element-edge connectivity
  integer, allocatable, dimension(:,:) :: kadj          ! Element adjacency
  integer, allocatable, dimension(:,:) :: karea         ! Element-face connectivity
  
  ! Secondary connectivity arrays
  integer, allocatable, dimension(:,:) :: kvel          ! Vertex-element connectivity
  integer, allocatable, dimension(:,:) :: kved          ! Edge-vertex connectivity
  integer, allocatable, dimension(:,:) :: kvar          ! Face-vertex connectivity
  
  ! Properties and markers
  integer, allocatable, dimension(:) :: knpr            ! Vertex/face properties
  real*8,  allocatable, dimension(:) :: dvol            ! Element volumes
END TYPE
```

### 2.2 Hierarchical Mesh Structure (`tMultiMesh`)

```fortran
TYPE tMultiMesh
  integer :: nlmin         ! Minimum refinement level
  integer :: nlmax         ! Maximum refinement level
  integer :: maxlevel      ! Usually nlmax + 1
  
  type(tMesh), allocatable, dimension(:) :: level  ! Array of mesh levels
END TYPE
```

### 2.3 COMMON Block Organization

FeatFloWer uses FORTRAN COMMON blocks for efficient mesh data management:

```fortran
! Current mesh dimensions
COMMON /TRIAD/ NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
               NVAR,NEAR,NBCT,NVBD,NEBD,NABD

! Array pointers for current level
COMMON /TRIAA/ LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
               LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
               LNPR,LBCT,LVBD,LEBD,LABD

! Multi-level array pointers
COMMON /MGTRA/ KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
               KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
               KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
               KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
               KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
               KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
               KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
```

---

## 3. Connectivity Arrays

### 3.1 KVERT - Element-Vertex Connectivity

**Purpose**: Maps each element to its 8 vertices  
**Dimensions**: `KVERT(8, NEL)`  
**Usage**: `KVERT(IVE, IEL)` = global vertex number for local vertex `IVE` of element `IEL`

**Hexahedral Vertex Ordering** (FEAT convention):
```
Bottom face (z-):     Top face (z+):
4 ——————— 3           8 ——————— 7
|         |           |         |
|         |           |         |
1 ——————— 2           5 ——————— 6
```

**Example Usage**:
```fortran
DO IEL = 1, NEL
  DO IVE = 1, 8
    IVT = KVERT(IVE, IEL)          ! Global vertex number
    X = DCORVG(1, IVT)             ! X-coordinate
    Y = DCORVG(2, IVT)             ! Y-coordinate  
    Z = DCORVG(3, IVT)             ! Z-coordinate
  END DO
END DO
```

### 3.2 KADJ - Element Adjacency

**Purpose**: Neighboring elements through faces  
**Dimensions**: `KADJ(6, NEL)`  
**Usage**: `KADJ(IFACE, IEL)` = neighboring element across face `IFACE` of element `IEL`  
**Boundary Indicator**: `KADJ(IFACE, IEL) = 0` indicates boundary face

**Face Numbering Convention**:
```
Face 1: Bottom (vertices 1,2,3,4)
Face 2: Front  (vertices 1,2,6,5)  
Face 3: Right  (vertices 2,3,7,6)
Face 4: Back   (vertices 3,4,8,7)
Face 5: Left   (vertices 4,1,5,8)
Face 6: Top    (vertices 5,6,7,8)
```

**Example - Finding Boundary Faces**:
```fortran
DO IEL = 1, NEL
  DO IFACE = 1, 6
    JNEL = KADJ(IFACE, IEL)
    IF (JNEL == 0) THEN
      ! This is a boundary face
      WRITE(*,*) 'Element', IEL, 'Face', IFACE, 'is on boundary'
    ELSE
      ! Process neighbor element JNEL
    END IF
  END DO
END DO
```

### 3.3 KEDGE - Element-Edge Connectivity

**Purpose**: Maps each element to its 12 edges  
**Dimensions**: `KEDGE(12, NEL)`  
**Usage**: `KEDGE(IEE, IEL)` = global edge number for local edge `IEE` of element `IEL`

**Edge Numbering Convention**:
```
Bottom edges (1-4): 1→2, 2→3, 3→4, 4→1
Vertical edges (5-8): 1→5, 2→6, 3→7, 4→8  
Top edges (9-12): 5→6, 6→7, 7→8, 8→5
```

**Static Local Edge Definition**:
```fortran
! In source/src_mesh/mesh_refine.f90
DIMENSION KIV(2,12)
DATA KIV /1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5/
```

### 3.4 KAREA - Element-Face Connectivity

**Purpose**: Maps each element to its 6 faces  
**Dimensions**: `KAREA(6, NEL)`  
**Usage**: `KAREA(IFACE, IEL)` = global face number for local face `IFACE` of element `IEL`

**Static Local Face Definition**:
```fortran
! In source/src_mesh/mesh_refine.f90
DIMENSION KIAD(4,6)
DATA KIAD/1,2,3,4, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/
```

### 3.5 KVEL - Vertex-Element Connectivity

**Purpose**: Lists all elements containing each vertex  
**Dimensions**: `KVEL(NVEL, NVT)` where `NVEL` is max elements per vertex  
**Usage**: `KVEL(J, IVT)` = J-th element containing vertex `IVT`  
**Termination**: `KVEL(J, IVT) = 0` indicates no more elements

**Example - Elements Around a Vertex**:
```fortran
DO IVT = 1, NVT
  J = 1
  DO WHILE (KVEL(J, IVT) /= 0 .AND. J <= NVEL)
    IEL = KVEL(J, IVT)
    ! Process element IEL containing vertex IVT
    J = J + 1
  END DO
END DO
```

### 3.6 KVED - Edge-Vertex Connectivity

**Purpose**: Maps each edge to its 2 vertices  
**Dimensions**: `KVED(2, NET)`  
**Usage**: `KVED(1:2, IEDG)` = global vertex numbers of edge `IEDG`

### 3.7 KVAR - Face-Vertex Connectivity

**Purpose**: Maps each face to its 4 vertices  
**Dimensions**: `KVAR(4, NAT)`  
**Usage**: `KVAR(1:4, IFACE)` = global vertex numbers of face `IFACE`

---

## 4. Coordinate Arrays

### 4.1 DCORVG - Vertex Coordinates

**Purpose**: 3D coordinates of all mesh vertices  
**Dimensions**: `DCORVG(3, NVT)`  
**Usage**: `DCORVG(1:3, IVT)` = (x,y,z) coordinates of vertex `IVT`

### 4.2 DCORAG - Face Center Coordinates

**Purpose**: 3D coordinates of face centers  
**Dimensions**: `DCORAG(3, NAT)`  
**Usage**: `DCORAG(1:3, IFACE)` = (x,y,z) coordinates of face center `IFACE`

---

## 5. Input Mesh File Format (.tri)

### 5.1 File Structure

FeatFloWer uses a custom `.tri` format for 3D hexahedral meshes. The format is parsed by `readTriCoarse()` in `source/src_mesh/mesh_refine.f90`.

```
Coarse mesh 3D
modified by tr2to3
   NEL   NVT  NBCT  NVE  NEE  NAE NEL,NVT,NBCT,NVE,NEE,NAE
 DCORVG
  x1 y1 z1
  x2 y2 z2
  ...
  xNVT yNVT zNVT
 KVERT
  v1_1 v1_2 v1_3 v1_4 v1_5 v1_6 v1_7 v1_8
  v2_1 v2_2 v2_3 v2_4 v2_5 v2_6 v2_7 v2_8
  ...
  vNEL_1 vNEL_2 vNEL_3 vNEL_4 vNEL_5 vNEL_6 vNEL_7 vNEL_8
 KNPR
  prop1
  prop2
  ...
  propNVT
```

### 5.2 Format Details

**Header Lines**:
- Line 1: Comment/description ("Coarse mesh 3D")
- Line 2: Comment/description ("modified by tr2to3")
- Line 3: Six integers defining mesh dimensions

**Sections**:
- **DCORVG**: Vertex coordinates (NVT lines, 3 floats per line)
- **KVERT**: Element-vertex connectivity (NEL lines, 8 integers per line)
- **KNPR**: Vertex properties/boundary markers (NVT lines, 1 integer per line)

**Example**:
```
Coarse mesh 3D
modified by tr2to3
      8      12       0       8      12       6 NEL,NVT,NBCT,NVE,NEE,NAE
 DCORVG
   0.0000000E+00   0.0000000E+00   0.0000000E+00
   1.0000000E+00   0.0000000E+00   0.0000000E+00
   1.0000000E+00   1.0000000E+00   0.0000000E+00
   0.0000000E+00   1.0000000E+00   0.0000000E+00
   0.0000000E+00   0.0000000E+00   1.0000000E+00
   1.0000000E+00   0.0000000E+00   1.0000000E+00
   1.0000000E+00   1.0000000E+00   1.0000000E+00
   0.0000000E+00   1.0000000E+00   1.0000000E+00
 KVERT
      1     2     3     4     5     6     7     8
 KNPR
 0
 0
 0
 0
 0
 0
 0
 0
```

---

## 6. Mesh Generation and Refinement

### 6.1 Connectivity Generation

**Primary Function**: `genMeshStructures()` in `source/src_mesh/mesh_refine.f90`

```fortran
subroutine genMeshStructures(mesh, extended, icurr, noe)
  ! Generates all connectivity arrays from basic KVERT
  call genKVEL(mesh%level(icurr))           ! Vertex-to-element connectivity
  call genKEDGE3(mesh%level(icurr), icurr, noe)  ! Edge numbering and connectivity
  call genKADJ2(mesh%level(icurr))          ! Element adjacency through faces
  call genKAREA(mesh%level(icurr))          ! Face numbering
  call genKVAR(mesh%level(icurr))           ! Face-to-vertex connectivity
  call genDCORAG(mesh%level(icurr))         ! Face center coordinates
end subroutine
```

### 6.2 Hierarchical Refinement

**Refinement Strategy**: 1:8 hexahedral subdivision
- Each parent hexahedron is divided into 8 child hexahedra
- Vertex coordinates are computed for edge midpoints, face centers, and element centers
- Connectivity is inherited and extended from parent level

**Key Functions**:
- `refineMesh()`: Controls multi-level refinement
- `refineMeshLevel()`: Single level 1:8 refinement
- `genMeshStructures()`: Regenerates connectivity after refinement

---

## 7. Mesh Access Patterns and Examples

### 7.1 Basic Element Iteration

```fortran
! Loop over all elements
DO IEL = 1, NEL
  ! Get element vertices
  DO IVE = 1, 8
    IVT = KVERT(IVE, IEL)
    X(IVE) = DCORVG(1, IVT)
    Y(IVE) = DCORVG(2, IVT)
    Z(IVE) = DCORVG(3, IVT)
  END DO
  
  ! Compute element volume, etc.
  CALL SETARE(DVOL, IEL, KVERT, DCORVG)
END DO
```

### 7.2 Neighbor Element Processing

```fortran
! Process element neighbors
DO IEL = 1, NEL
  DO IFACE = 1, 6
    JNEL = KADJ(IFACE, IEL)
    IF (JNEL > 0) THEN
      ! Process interior face between IEL and JNEL
      IAREA = KAREA(IFACE, IEL)
      ! Get face vertices: KVAR(1:4, IAREA)
    ELSE
      ! Process boundary face
      CALL ApplyBoundaryCondition(IEL, IFACE)
    END IF
  END DO
END DO
```

### 7.3 Vertex-Based Operations

```fortran
! Loop over vertices and their elements
DO IVT = 1, NVT
  NELEM = 0
  DO J = 1, NVEL
    IEL = KVEL(J, IVT)
    IF (IEL == 0) EXIT
    NELEM = NELEM + 1
    ! Process element IEL containing vertex IVT
  END DO
  WRITE(*,*) 'Vertex', IVT, 'belongs to', NELEM, 'elements'
END DO
```

### 7.4 Finite Element Assembly Pattern

```fortran
! Typical finite element assembly loop
DO IEL = 1, NEL
  ! Get degrees of freedom for this element
  CALL NDFGL(IEL, 1, IELTYP, KVERT, KEDGE, KAREA, KDFG, KDFL)
  
  ! Loop over local DOFs
  DO I = 1, IDFL
    IG = KDFG(I)    ! Global DOF number
    LDOF = KDFL(I)  ! Local basis function index
    
    ! Access solution values
    U_LOCAL(I) = U_GLOBAL(IG)
  END DO
  
  ! Perform element computations
  CALL ElementAssembly(...)
END DO
```

### 7.5 Boundary Condition Application

```fortran
! Apply boundary conditions using face properties
DO IEL = 1, NEL
  DO IFACE = 1, 6
    IF (KADJ(IFACE, IEL) == 0) THEN  ! Boundary face
      IAREA = KAREA(IFACE, IEL)
      BC_TYPE = KNPR(NVT + IAREA)     ! Boundary condition marker
      
      SELECT CASE(BC_TYPE)
        CASE(1)  ! Dirichlet boundary
          CALL ApplyDirichletBC(IEL, IFACE)
        CASE(2)  ! Neumann boundary  
          CALL ApplyNeumannBC(IEL, IFACE)
        CASE(3)  ! Mixed boundary
          CALL ApplyMixedBC(IEL, IFACE)
      END SELECT
    END IF
  END DO
END DO
```

---

## 8. Performance Considerations

### 8.1 Memory Layout

**Optimized Access Patterns**:
- Arrays are dimensioned for cache-friendly access: `ARRAY(local_dim, global_dim)`
- Sequential element processing minimizes cache misses
- Connectivity arrays are pre-computed and stored for fast access

### 8.2 Vectorization and Parallelization

**Loop Structures**:
- Inner loops over local element entities (vertices, faces, edges)
- Outer loops over global entities (elements, vertices)
- Suitable for OpenMP parallelization and compiler vectorization

### 8.3 Memory Management

**COMMON Block Strategy**:
- Uses pointer-based workspace management (`L(NNARR)` arrays)
- Efficient memory allocation for multi-level structures
- Supports dynamic mesh adaptation and refinement

---

## 9. Integration with Finite Element Operations

### 9.1 DOF Management

The mesh connectivity seamlessly integrates with DOF numbering through `NDFGL()`:

```fortran
CALL NDFGL(IEL, 1, IELTYP, KVERT, KEDGE, KAREA, KDFG, KDFL)
! KDFG: Global DOF numbers for element IEL
! KDFL: Local basis function indices
```

### 9.2 Quadrature Point Processing

```fortran
! Element loop with quadrature
DO IEL = 1, NEL
  ! Set up element geometry
  DO IVE = 1, 8
    IVT = KVERT(IVE, IEL)
    DX(IVE) = DCORVG(1, IVT)
    DY(IVE) = DCORVG(2, IVT)
    DZ(IVE) = DCORVG(3, IVT)
  END DO
  
  ! Quadrature loop
  DO ICUBP = 1, NCUBP
    ! Evaluate basis functions and derivatives
    CALL ELE(XI1, XI2, XI3, -3)
    
    ! Interpolate solution at quadrature point
    ! using DBAS and local DOF values
  END DO
END DO
```

---

## 10. Debugging and Validation

### 10.1 Mesh Consistency Checks

```fortran
! Verify element-vertex connectivity
DO IEL = 1, NEL
  DO IVE = 1, 8
    IVT = KVERT(IVE, IEL)
    IF (IVT < 1 .OR. IVT > NVT) THEN
      WRITE(*,*) 'Invalid vertex reference in element', IEL
    END IF
  END DO
END DO

! Verify adjacency symmetry
DO IEL = 1, NEL
  DO IFACE = 1, 6
    JNEL = KADJ(IFACE, IEL)
    IF (JNEL > 0) THEN
      ! Find corresponding face in neighbor
      FOUND = .FALSE.
      DO JFACE = 1, 6
        IF (KADJ(JFACE, JNEL) == IEL) FOUND = .TRUE.
      END DO
      IF (.NOT. FOUND) THEN
        WRITE(*,*) 'Adjacency not symmetric:', IEL, JNEL
      END IF
    END IF
  END DO
END DO
```

### 10.2 Mesh Quality Assessment

```fortran
! Compute element aspect ratios and volumes
DO IEL = 1, NEL
  CALL SETARE(VOL, IEL, KVERT, DCORVG)
  IF (VOL <= 0.0) THEN
    WRITE(*,*) 'Negative volume element:', IEL
  END IF
  
  ! Check for degenerate elements
  CALL CheckElementQuality(IEL, KVERT, DCORVG, QUALITY)
  IF (QUALITY < THRESHOLD) THEN
    WRITE(*,*) 'Poor quality element:', IEL, QUALITY
  END IF
END DO
```

---

## 11. Summary

The FeatFloWer mesh structure provides a comprehensive and efficient framework for 3D hexahedral finite element computations. Key features include:

- **Complete connectivity**: Full element-vertex-edge-face relationships
- **Hierarchical support**: Multi-level grids for multigrid methods
- **Efficient access**: Cache-friendly data layouts and access patterns
- **Robust I/O**: Well-defined `.tri` file format for mesh exchange
- **Integration**: Seamless connection with finite element DOF management

This mesh organization enables FeatFloWer to handle complex 3D fluid flow simulations with sophisticated geometries while maintaining computational efficiency and numerical accuracy.