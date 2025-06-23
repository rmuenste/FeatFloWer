# Domain Decomposition and Parallelization in FeatFloWer

This document describes the domain decomposition and parallelization strategy implemented in FeatFloWer, a parallel Q2/P1 finite element framework for computational fluid dynamics. The approach combines METIS-based graph partitioning with hierarchical mesh refinement and sophisticated communication patterns optimized for multigrid solvers.

---

## 1. Overview

FeatFloWer employs a **master-worker parallelization paradigm** with geometric domain decomposition. The approach is specifically designed for:

- **Q2/P1 mixed finite elements**: Quadratic velocity (27-node) with linear pressure (4-node discontinuous)
- **Multigrid solvers**: Hierarchical communication across multiple grid levels
- **Large-scale 3D CFD**: Efficient scaling for realistic fluid flow simulations
- **Fictitious Boundary Method**: Parallel handling of moving particles and complex geometries

The parallelization workflow follows this sequence:

1. **Initial partitioning** using METIS graph partitioner
2. **Hierarchical mesh refinement** with automatic connectivity generation
3. **Communication structure setup** for parallel matrix-vector operations
4. **Global DOF numbering** for consistent linear algebra operations

---

## 2. Mesh Partitioning Strategy

### 2.1 METIS-Based Graph Partitioning

**Location**: `/partitioner/` executable  
**Dependencies**: METIS 4.0.3 (`extern/libraries/metis-4.0.3/`)

The initial coarse mesh is partitioned using METIS graph partitioning:

```bash
./partitioner <nprocs> 1 <nSubCoarseMesh> "<mesh_name>" "<project_file>" 1> /dev/null
```

**Partitioning Process:**
- Converts mesh connectivity to graph representation
- Minimizes edge cuts (inter-partition communication)
- Balances element distribution across processes
- Generates partition files: `_mesh/<mesh_name>/sub<XXX>/GRID<YYYY>.tri`

### 2.2 Partition Assignment

**Files**: `PartitionReader.f90` and `PartitionReader2.f90`

Two strategies for assigning partitions to MPI processes:

**Simple Assignment (PartitionReader.f90):**
```fortran
! Master process (rank 0): GRID.tri
! Worker process i: sub<i>/GRID0001.tri
```

**Load-Balanced Assignment (PartitionReader2.f90):**
```fortran
kSubPart = FLOOR(DBLE(subnodes)/DBLE(nSubCoarseMesh)-1d-10)+1
iSubPart = FLOOR(DBLE(myid)/DBLE(kSubPart)-1d-10)+1
iPart    = myid - (iSubPart-1)*kSubPart
```

This distributes sub-partitions more evenly when `nSubCoarseMesh` > number of processes.

---

## 3. Hierarchical Mesh Refinement

### 3.1 Multi-Level Mesh Generation

**Key Functions** (`source/src_mesh/mesh_refine.f90`):
- `readTriCoarse()`: Reads initial partitioned mesh
- `refineMesh()`: Performs hierarchical 1:8 hexahedral refinement
- `genMeshStructures()`: Generates connectivity arrays (KEDGE, KADJ, KAREA)

**Refinement Strategy:**
```fortran
call readTriCoarse(CMESH1, mg_mesh)
call refineMesh(mg_mesh, mg_Mesh%maxlevel, .true.)
```

Each refinement level:
- Subdivides hexahedra 1:8 (each hex → 8 child hexes)
- Maintains geometric consistency across partition boundaries
- Generates required connectivity for Q2/P1 elements

### 3.2 Level Distribution

**Master Process (rank 0):**
- `NLMAX = LinSc%prm%MGprmIn%MedLev`
- Handles fewer levels to balance computational load
- Coordinates global operations and coarse grid solves

**Worker Processes (rank ≠ 0):**
- `NLMAX = LinSc%prm%MGprmIn%MedLev + 1`
- Perform local computations on finest mesh level
- Handle majority of finite element assembly and local solves

---

## 4. Communication Structure Setup

### 4.1 Communication Hierarchy

The parallel communication setup follows a specific sequence in `General_init_ext()`:

```fortran
! 1. Initial geometric mapping
CALL PARENTCOMM(...)

! 2. Multi-level communication patterns  
DO II=NLMIN+1,NLMAX
  CALL CREATECOMM(II,...)
END DO

! 3. Q2 element communication (coarse level)
CALL E013_CreateComm_coarse(...)

! 4. Global DOF numbering
CALL Create_GlobalNumbering(...)

! 5. Q2 communication (refined levels)
DO ILEV=NLMIN+1,NLMAX
  CALL E013_CreateComm(...)
END DO

! 6. Communication optimization
CALL ExtraxtParallelPattern()

! 7. Final DOF communication
CALL E011_CreateComm(NDOF)
```

### 4.2 Key Communication Subroutines

#### PARENTCOMM (`source/src_mpi/parentcomm.f90`)
**Purpose**: Master-worker mesh mapping and initial communication setup

- Creates mapping between master mesh and distributed worker meshes
- Uses octree search for efficient geometric matching
- Establishes element-to-element and vertex-to-vertex correspondence
- Sets up initial communication patterns for face-based data exchange
- Handles periodic boundary conditions

**Key Data Structures:**
```fortran
coarse%pELEMLINK  ! Element mapping: local → global
coarse%pVERTLINK  ! Vertex mapping: local → global  
coarse%pFACELINK  ! Face mapping for communication
```

#### CREATECOMM (`source/src_mpi/pp3d_mpi.f90`)
**Purpose**: Multigrid level communication pattern creation

- Creates communication structures for different multigrid levels
- Refines face-based communication patterns from coarse to fine grids
- Each level multiplies communication size by 4 (for Q2 elements)
- Links local face indices to element connectivity

#### E013_CreateComm (`source/src_quadLS/QuadSc_mpi.f90`)
**Purpose**: Q2 element communication for quadratic finite elements

- Builds communication structures specific to 27-node Q2 elements
- Links local DOFs to global DOFs through geometric coordinate matching
- Creates send/receive patterns for matrix-vector operations
- Handles vertices, edge midpoints, face centers, and element centers

#### Create_GlobalNumbering (`source/src_quadLS/QuadSc_mpi.f90`)
**Purpose**: Establishes global DOF numbering for parallel linear algebra

- Creates consistent global numbering across all processes
- Handles all Q2/P1 DOF types:
  - **Vertices**: Corner nodes of elements
  - **Edges**: Midpoint nodes on element edges
  - **Faces**: Center nodes on element faces  
  - **Elements**: Center nodes within elements
- Essential for parallel matrix assembly and solver operations

---

## 5. Master-Worker Coordination

### 5.1 Process Roles

**Master Process (rank 0):**
- Coordinates mesh partitioning via external METIS partitioner
- Reads full coarse mesh and broadcasts to workers
- Manages fewer refinement levels for load balancing
- Handles global matrix assembly and coarse grid solves
- Coordinates data collection and output operations

**Worker Processes (rank ≠ 0):**
- Read individual partition mesh files
- Perform local mesh refinement to target level
- Handle local finite element assembly and computation
- Participate in parallel communication via MPI structures
- Contribute to global reductions and data exchange

### 5.2 Load Balancing

**Computational Load Distribution:**
- Master handles coordination overhead vs. direct computation
- Workers handle majority of finite element operations
- Refinement levels adjusted per process type
- Communication patterns optimized to minimize synchronization

**Memory Distribution:**
- Each process stores only local mesh partition
- Communication buffers sized for local interface elements
- Global data structures minimized to reduce memory footprint

---

## 6. Communication Patterns and Data Exchange

### 6.1 Face-Based Communication

**Primary Communication Entity**: Shared faces between partition boundaries

- Face identification during mesh connectivity generation
- Mapping of local face indices to neighboring partition faces
- Exchange of DOF values across partition interfaces
- Consistent boundary condition application

### 6.2 DOF Communication Types

**Vertex Communication (P1 elements):**
- Simpler pattern for linear pressure elements
- Direct vertex-to-vertex mapping between neighbors
- Used for pressure degrees of freedom

**Q2 Element Communication:**
- Complex pattern for 27-node quadratic elements
- Includes vertices, edge midpoints, face centers, element centers
- Communication size scales with element refinement (×4 per level)

**Communication Optimization:**
- `ExtraxtParallelPattern()` analyzes communication dependencies
- Creates optimized message passing schedules
- Minimizes synchronization overhead between processes

---

## 7. Global Numbering and Linear Algebra

### 7.1 Local-to-Global Mapping

**Local Numbering**: Within each process partition
- Elements numbered 1 to `nel_local`
- Vertices numbered 1 to `nvt_local`
- Consistent with local mesh connectivity

**Global Numbering**: Across all processes
- Established through coordinate matching via `Create_GlobalNumbering()`
- Enables consistent parallel matrix assembly
- Required for parallel iterative solvers

### 7.2 DOF Management

**Q2/P1 Mixed Elements:**
```fortran
NDOF = nvt + nat + nel + net  ! Total DOFs per level
```

- `nvt`: Number of vertices (pressure P1 + velocity Q2 corners)
- `nat`: Number of faces (velocity Q2 face centers)  
- `nel`: Number of elements (velocity Q2 element centers)
- `net`: Number of edges (velocity Q2 edge midpoints)

---

## 8. Multigrid Parallelization

### 8.1 Hierarchical Communication

**Multi-Level Strategy:**
- Each multigrid level has independent communication pattern
- Coarse grids handled by master or specialized communication  
- Fine grids use distributed communication patterns
- Smooth transitions between levels via prolongation/restriction

**Communication Scaling:**
- Communication volume scales with interface area
- Multigrid reduces communication overhead vs. single-level
- Optimal for large 3D problems with geometric multigrid

### 8.2 Solver Integration

**Parallel Multigrid Operations:**
- Local smoothing operations on each process
- Inter-process communication for restriction/prolongation
- Coarse grid solves coordinated by master process
- Efficient scaling for large processor counts

---

## 9. Implementation Details

### 9.1 Key Files and Modules

| File/Module | Purpose |
|-------------|---------|
| `applications/*/app_init.f90` | Main initialization sequence |
| `source/src_mesh/mesh_refine.f90` | Mesh refinement and connectivity |
| `source/src_mpi/parentcomm.f90` | Master-worker communication setup |
| `source/src_mpi/pp3d_mpi.f90` | MPI data structures and utilities |
| `source/src_quadLS/QuadSc_mpi.f90` | Q2/P1 parallel communication |
| `partitioner/` | METIS-based mesh partitioning tool |

### 9.2 Configuration Parameters

**Key Parameters:**
- `nSubCoarseMesh`: Controls initial mesh subdivision
- `NLMIN`, `NLMAX`: Multigrid level range per process
- `MedLev`: Target finest level for computation
- Communication buffer sizes automatically determined

---

## 10. Performance Characteristics

### 10.1 Scalability Features

**Strong Scaling**: Fixed problem size, increasing processors
- Communication overhead grows with interface/volume ratio
- Optimal for moderately large processor counts (10-1000 cores)
- METIS partitioning minimizes communication volume

**Weak Scaling**: Problem size scales with processors  
- Communication patterns remain consistent
- Memory usage per process remains bounded
- Suitable for very large 3D problems

### 10.2 Optimization Strategies

**Communication Minimization:**
- METIS partitioning reduces edge cuts
- Face-based communication eliminates redundant exchanges
- Optimized message scheduling via `ExtraxtParallelPattern()`

**Load Balancing:**
- Master-worker role differentiation
- Dynamic level assignment based on computational load
- Geometric partitioning accounts for element distribution

---

## 11. Usage and Best Practices

### 11.1 Typical Execution

```bash
# 1. Generate partitioned mesh
./partitioner 8 1 8 "channel_mesh" "project.dat"

# 2. Parallel execution  
mpirun -np 8 ./applications/q2p1_devel/q2p1_devel
```

### 11.2 Performance Tuning

**Mesh Partitioning:**
- Choose `nSubCoarseMesh` ≈ number of MPI processes
- Balance partition sizes for uniform load distribution
- Consider communication topology of target hardware

**Multigrid Levels:**
- Set `MedLev` based on problem size and available memory
- More levels reduce communication but increase setup cost
- Optimal range typically 3-6 levels for 3D problems

**Memory Management:**
- Monitor per-process memory usage
- Adjust refinement levels if memory becomes limiting
- Consider coarser partitioning for memory-constrained systems

---

## 12. Summary

FeatFloWer's domain decomposition strategy represents a sophisticated approach to parallel CFD that effectively combines:

1. **METIS-based graph partitioning** for optimal load balancing
2. **Hierarchical mesh refinement** with automatic connectivity generation
3. **Master-worker coordination** with balanced computational distribution
4. **Comprehensive communication structures** optimized for Q2/P1 elements
5. **Multigrid-aware parallelization** for efficient iterative solvers

This design enables FeatFloWer to handle large-scale 3D fluid flow simulations with complex geometries while maintaining the accuracy benefits of higher-order Q2/P1 finite element discretizations. The approach scales effectively on modern HPC systems and provides a robust foundation for advanced CFD applications including fluid-structure interaction and dense particle suspensions.