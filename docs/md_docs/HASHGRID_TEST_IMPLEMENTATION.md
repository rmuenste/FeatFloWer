# HashGrid Acceleration Test - Implementation Summary

## Overview

This document summarizes the implementation of the HashGrid acceleration test application `q2p1_hashgrid_test` for verifying the correctness of accelerated point-inside queries.

**Latest Update (2026-02-06):** Extended with support for mixed shape types (spheres, boxes, capsules) and triangle mesh testing (cone.obj). The test infrastructure now supports comprehensive validation of HashGrid acceleration across all geometry primitives supported by the PE library.

## Extended Testing Capabilities (2026-02-06)

### Packing Methods

The test application now supports multiple packing configurations via JSON configuration files:

#### 1. Grid Packing (Original)
**Config:** `packingMethod_: "grid"`
- 1000 spheres in 10×10×10 grid
- Sphere radius: 0.01
- Spacing: 0.1
- Domain: [0,1]³

#### 2. Mixed Grid Packing (NEW)
**Config:** `packingMethod_: "mixed_grid"`
**Configuration file:** `pe_user_config.json`

Cycles through three shape types with uniform bounding radius 0.02:
- **Spheres:** radius = 0.02
- **Boxes:** cube side ≈ 0.0231 (fits in bounding sphere)
- **Capsules:** radius = 0.01, length = 0.02

**Distribution:** 10×10×10 grid = 1000 objects
- ~333 spheres + ~333 boxes + ~334 capsules
- Shapes cycle: sphere → box → capsule → sphere → ...

**Purpose:** Validates HashGrid acceleration with different primitive containment algorithms

#### 3. Triangle Mesh Grid Packing (NEW)
**Config:** `packingMethod_: "tm_grid"`
**Configuration file:** `pe_user_config_tm.json`

Loads triangle mesh from `cone.obj` file:
- 1000 cone meshes in 10×10×10 grid
- Each cone fits within bounding radius 0.02
- Tests plane-based containment algorithm (without CGAL)
- Tests AABB early-out optimization

**Purpose:** Validates HashGrid with complex geometries and verifies triangle mesh containment fallback algorithms

### Configuration Files

**File:** `build-ninja-accelerated/applications/q2p1_hashgrid_test/pe_user_config.json`
```json
{
    "timesteps_": 3,
    "stepsize_": 0.001,
    "packingMethod_": "mixed_grid",
    "particleDensity_": 1.0,
    "fluidViscosity_": 0.01,
    "fluidDensity_": 1.0,
    "gravity_": [0.0, 0.0, 0.0]
}
```

**File:** `build-ninja-accelerated/applications/q2p1_hashgrid_test/pe_user_config_tm.json`
```json
{
    "timesteps_": 3,
    "stepsize_": 0.001,
    "packingMethod_": "tm_grid",
    "particleDensity_": 1.0,
    "fluidViscosity_": 0.01,
    "fluidDensity_": 1.0,
    "gravity_": [0.0, 0.0, 0.0]
}
```

### Key Implementation Changes

#### SimulationConfig Extensions
**File:** `libs/pe/pe/config/SimulationConfig.h`
```cpp
enum PackingMethod {
    Grid,             //!< Regular grid packing (spheres only)
    MixedGrid,        //!< Mixed shape grid packing (spheres, boxes, capsules)
    TriangleMeshGrid, //!< Triangle mesh grid packing (cones from cone.obj)
    External,         //!< Load from external file
    None              //!< No initial packing
};
```

#### Box Containment Fix
**File:** `libs/pe/src/interface/object_queries.cpp`

Added missing `boxType` handling in three locations:
- `pointInsideParticles()` - local bodies section
- `pointInsideParticles()` - shadow copies section
- `pointInsideParticlesAccelerated()` - HashGrid accelerated path

**Impact:** Boxes are now correctly detected by CFD mesh queries

#### Triangle Mesh Containment Algorithm
**File:** `libs/pe/pe/core/rigidbody/trianglemeshtrait/Default.h`

Implemented plane-based containment fallback (O(N)) when CGAL is unavailable:
```cpp
// Separating plane algorithm:
// Point is inside if it's behind all triangle planes
for each triangle face:
    compute normal = (B - A) × (C - A)
    if (point - A) · normal > 0:
        return false  // Point in front of plane → outside
return true  // Behind all planes → inside
```

Added AABB early-out optimization (O(1)):
```cpp
// Fast rejection before expensive containment check
if (!getAABB().contains(point)):
    return false  // Outside bounding box
```

**Performance:** AABB rejects ~90%+ of queries instantly for sparse grids

### Testing Workflow

**Test all three packing methods:**

```bash
# 1. Test spheres only (original)
cd build-ninja-accelerated/applications/q2p1_hashgrid_test
mpirun -np 4 ./q2p1_hashgrid_test

# 2. Test mixed shapes (spheres, boxes, capsules)
cp pe_user_config.json ./  # packingMethod: "mixed_grid"
mpirun -np 4 ./q2p1_hashgrid_test

# 3. Test triangle meshes
cp pe_user_config_tm.json pe_user_config.json
cp /path/to/cone.obj ./
mpirun -np 4 ./q2p1_hashgrid_test
```

**Expected verification output:**
```
==================================================
HashGrid Verification & Performance Summary
==================================================

--- Correctness ---
Total queries:         17869
Result:                ✓ Perfect match (0 mismatches)

--- Performance Comparison ---
Accelerated (HashGrid):    0.002407 s
  Per query:                  0.135 μs

Baseline (Linear):        0.221373 s
  Per query:                 12.389 μs

Speedup:                      91.97x
==================================================
```

## Files Created/Modified

### 1. PE Library - Setup Function
**File**: `libs/pe/pe/interface/sim_setup_serial.h`
**Changes**: Added `setupHashGridTest(int cfd_rank)` function (lines ~1712-1920)

**Purpose**: Creates test particles based on JSON configuration
**Key Features**:
- **JSON-driven configuration:** Loads `pe_user_config.json` for parameters
- **Multiple packing methods:** Grid (spheres), MixedGrid (spheres/boxes/capsules), TriangleMeshGrid (cones)
- **Configurable parameters:** Timesteps, stepsize, fluid properties, gravity
- **Iron material:** PE standard material for consistency
- **10×10×10 grid:** 1000 objects in [0,1]³ domain with 0.1 spacing
- **Diagnostic output:** Detailed setup summary from rank 1

### 2. PE Library - C/Fortran Interface
**File**: `libs/pe/src/interface/c2f_interface.cpp`
**Changes**: Added `commf2c_hashgrid_test_()` function (lines ~417-420)

**Purpose**: Routes Fortran call to C++ PE setup function
**Function signature**:
```cpp
extern "C" void commf2c_hashgrid_test_(int *Fcomm, int *FcommEx0, int *remoteRank)
```

### 3. Test Application - Main Program
**File**: `applications/q2p1_hashgrid_test/q2p1_hashgrid_test.f90`
**Status**: NEW

**Purpose**: Minimal main program that:
- Initializes MPI and mesh
- Calls single alpha field computation
- Reports results
- Handles --version flag

**Key Functions**:
- `init_hashgrid_test()` - Initialization
- `updateFBMGeometry()` - Alpha field computation (point queries)

### 4. Test Application - Initialization
**File**: `applications/q2p1_hashgrid_test/app_init.f90`
**Status**: NEW

**Purpose**: Simplified initialization based on q2p1_ATC
**Keeps**:
- MPI initialization
- Mesh loading and partition reading
- Q2/P1 structure setup
- PE initialization (calls `commf2c_hashgrid_test`)
- FBM particle retrieval

**Removes**:
- NS solver setup
- Multigrid initialization
- Transport equation setup
- Mesh deformation (umbrella smoothing)
- HYPRE solver setup

### 5. Test Application - CMake Configuration
**File**: `applications/q2p1_hashgrid_test/CMakeLists.txt`
**Status**: NEW

**Purpose**: Build configuration for test application
**Features**:
- Links to FF_APPLICATION_LIBS
- Version header generation
- Out-of-source build support

### 6. Build System Registration
**File**: `applications/CMakeLists.txt`
**Changes**: Added `add_subdirectory(q2p1_hashgrid_test)` (line 21)

**Purpose**: Register test application with CMake build system

### 7. Documentation
**Files**:
- `applications/q2p1_hashgrid_test/README.md`
- `HASHGRID_TEST_IMPLEMENTATION.md` (this file)

**Purpose**: Usage instructions and implementation notes

## Build Configurations

### Configuration 1: Baseline (No Acceleration)
```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_PE=ON \
      -DUSE_PE_SERIAL_MODE=ON \
      -DUSE_ACCELERATED_POINT_QUERY=OFF \
      -DBUILD_APPLICATIONS=ON \
      ..
```

### Configuration 2: Accelerated
```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_PE=ON \
      -DUSE_PE_SERIAL_MODE=ON \
      -DUSE_ACCELERATED_POINT_QUERY=ON \
      -DBUILD_APPLICATIONS=ON \
      ..
```

## Test Flow

1. **MPI Initialization**
   - Each MPI rank initializes
   - Rank 0 = master (no computation)
   - Ranks 1-N = workers (mesh + particles)

2. **Mesh Loading**
   - Read coarse mesh from `.tri` file
   - Refine to target level (NLMAX)
   - Set up parallel communication structures

3. **PE Setup**
   - Create MPI communicator excluding master
   - Call `commf2c_hashgrid_test()` on worker ranks
   - PE library creates 1000 spheres in 10×10×10 grid

4. **Particle Retrieval**
   - `FBM_GetParticles()` queries PE for particle list
   - `FBM_ScatterParticles()` distributes to worker ranks

5. **Alpha Field Computation**
   - `updateFBMGeometry()` called once
   - For each mesh node, checks if inside any sphere
   - Uses either linear search (baseline) or HashGrid (accelerated)
   - Outputs "Total dofs inside" count

6. **Result Validation**
   - Compare "Total dofs inside" between baseline and accelerated
   - Counts must be identical for test to pass

## Key Code Paths

### Point Query Path (Simplified)
```
updateFBMGeometry() [QuadSc_geometry_utilities.f90]
  └─> QuadScalar_FictKnpr() [QuadSc_boundary.f90]
      └─> FBM_IsNodeInsideParticle() [cinterface.f90]
          └─> isPointInsideParticle() [object_queries.cpp]
              └─> [IF ACCELERATED]
                  CollisionSystem::pointQuery() [c_interface_particle_fbm.h]
                    └─> HashGrid fine collision detection
              └─> [IF BASELINE]
                  pointInsideParticleLinearSearch() [object_queries.cpp]
                    └─> Loop over all particles
```

## Verification Points

### 1. PE Setup Executed
Look for in output:
```
--HASHGRID TEST SETUP
Domain: [0, 1]^3
Spheres: 1000 (10x10x10 grid)
...
```

### 2. Particles Loaded
Look for:
```
Number of particles: 1000
```

### 3. Alpha Field Computed
Look for:
```
> Total dofs inside: XXXXX
> Dofs per Particle: YYY
```

### 4. Counts Match
```bash
grep "Total dofs inside" output_baseline.txt
grep "Total dofs inside" output_accelerated.txt
# Values should be identical
```

## Dependencies

### External Libraries
- MPI (required for parallel execution)
- PE library (physics engine, built with USE_PE=ON)

### Internal Modules
- `def_FEAT` - FeatFlow definitions
- `PP3D_MPI` - MPI utilities
- `MESH_Structures` - Mesh data structures
- `var_QuadScalar` - Scalar field variables
- `cinterface` - C/Fortran interface

### PE Components
- `pe::World` - Physics world
- `pe::CollisionSystem` - Collision detection
- `pe::Material` - Material properties
- `pe::createSphere()` - Sphere primitive creation

## Implementation Notes

### Why Serial PE Mode?
Serial PE mode is used because:
1. Each CFD domain runs independent PE instance
2. All particles are "local" to each domain
3. No MPI communication within PE library
4. Forces synchronized via CFD's MPI layer
5. Simpler setup for large particles (< 20 bodies)

### Why 1000 Objects?
1. Large enough to trigger HashGrid activation (threshold ~100)
2. Small enough for reasonable runtime and verification overhead
3. Sufficient spatial coverage for testing (fills unit cube with 0.1 spacing)
4. 10×10×10 grid is easy to reason about and debug
5. Tests HashGrid performance with realistic particle counts

### Why Iron Material?
Iron is a standard PE material with well-defined properties used in many benchmarks. Using it ensures consistency with other tests.

### Why [0,1]³ Domain?
1. Simple unit cube for easy verification
2. Spheres fit comfortably with 0.1 spacing
3. No boundary effects from domain edges
4. Standard domain size in CFD tests

## Troubleshooting Guide

### Problem: Build Fails with PE Errors
**Solution**: Ensure two-step CMake configuration:
```bash
cmake -DUSE_PE=ON ..
cmake -DUSE_PE_SERIAL_MODE=ON ..
```

### Problem: Different Counts Between Runs
**Diagnosis**:
1. Verify both use PE_SERIAL_MODE
2. Check CMake cache for USE_ACCELERATED_POINT_QUERY
3. Add debug output to point query functions
4. Test with fewer spheres (5×5×5)

### Problem: No Particles Created
**Diagnosis**:
1. Check PE setup output for "Spheres: 1000"
2. Verify `commf2c_hashgrid_test()` is called
3. Check PE library linked correctly

### Problem: Mesh Loading Fails
**Solution**:
1. Ensure mesh file exists in `_data/`
2. Check partition files match mesh
3. Verify NLMAX in parameters

## Completed Enhancements (2026-02-06)

### Implemented Features ✅
1. **Multiple Geometries**: ✅ Spheres, boxes, capsules, triangle meshes all supported
2. **JSON Configuration**: ✅ Parameterized via pe_user_config.json
3. **Performance Metrics**: ✅ Detailed timing and speedup reporting in verification mode
4. **Mixed Shape Testing**: ✅ MixedGrid packing validates all primitive types
5. **Triangle Mesh Support**: ✅ TriangleMeshGrid with cone.obj loading

### Future Enhancements

### Potential Improvements
1. **Parameterized Grid**: Allow command-line override of JSON parameters
2. **Automated Testing**: Add to CI/CD pipeline with regression tests
3. **Visual Output**: Generate VTK for particle positions and alpha field
4. **Benchmark Suite**: Standardized performance benchmarks across platforms

### Extended Tests
1. **Below Threshold**: Test with < 100 objects (verify HashGrid activation threshold)
2. **Mixed Sizes**: Varying particle radii and scales
3. **Non-uniform Distribution**: Random particle positions
4. **Boundary Cases**: Particles at domain edges and overlapping boundaries
5. **Complex Meshes**: Test with higher-resolution triangle meshes (> 1000 triangles)

## References

### Related Files
- `libs/pe/pe/interface/sim_setup_serial.h` - PE setup functions
- `libs/pe/src/interface/c2f_interface.cpp` - C/Fortran interface
- `libs/pe/src/interface/object_queries.cpp` - Point query implementations
- `source/src_quadLS/QuadSc_boundary.f90` - Alpha field computation
- `source/src_quadLS/QuadSc_geometry_utilities.f90` - FBM geometry update

### Plan Document
See the original plan in the session transcript for detailed specification.

## Completion Checklist

### Original Implementation
- [x] Created `setupHashGridTest()` in `sim_setup_serial.h`
- [x] Added `commf2c_hashgrid_test_()` in `c2f_interface.cpp`
- [x] Created `q2p1_hashgrid_test.f90` main program
- [x] Created `app_init.f90` initialization
- [x] Created `CMakeLists.txt` build configuration
- [x] Registered application in `applications/CMakeLists.txt`
- [x] Created README with usage instructions
- [x] Created implementation summary (this document)
- [x] Built and tested baseline version
- [x] Built and tested accelerated version
- [x] Verified correctness with sphere grid

### Extended Testing (2026-02-06)
- [x] Added MixedGrid and TriangleMeshGrid to SimulationConfig
- [x] Implemented mixed shape grid packing (spheres/boxes/capsules)
- [x] Implemented triangle mesh grid packing (cone.obj)
- [x] Fixed box containment detection in object_queries.cpp
- [x] Implemented plane-based triangle mesh containment algorithm
- [x] Added AABB early-out optimization for triangle meshes
- [x] Created pe_user_config.json for mixed grid testing
- [x] Created pe_user_config_tm.json for triangle mesh testing
- [x] Verified correctness with all shape types (spheres, boxes, capsules, meshes)
- [x] Documented deferred acceleration strategy (hashgrid_initialization_fix.md)
- [x] Committed PE library changes and main repo updates

## Next Steps

1. **Build Baseline**:
   ```bash
   cd /data/warehouse17/rmuenste/code/FF-ATC-NEW/ff-accel
   mkdir build_test_baseline
   cd build_test_baseline
   cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PE=ON -DBUILD_APPLICATIONS=ON ..
   cmake -DUSE_PE_SERIAL_MODE=ON -DUSE_ACCELERATED_POINT_QUERY=OFF ..
   make -j8 q2p1_hashgrid_test
   ```

2. **Build Accelerated**:
   ```bash
   cd /data/warehouse17/rmuenste/code/FF-ATC-NEW/ff-accel
   mkdir build_test_accelerated
   cd build_test_accelerated
   cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PE=ON -DBUILD_APPLICATIONS=ON ..
   cmake -DUSE_PE_SERIAL_MODE=ON -DUSE_ACCELERATED_POINT_QUERY=ON ..
   make -j8 q2p1_hashgrid_test
   ```

3. **Prepare Test Data**:
   - Copy mesh file to `applications/q2p1_hashgrid_test/_data/`
   - Ensure partition files are available

4. **Run Tests**:
   ```bash
   cd build_test_baseline/applications/q2p1_hashgrid_test
   mpirun -np 4 ./q2p1_hashgrid_test > output_baseline.txt 2>&1

   cd build_test_accelerated/applications/q2p1_hashgrid_test
   mpirun -np 4 ./q2p1_hashgrid_test > output_accelerated.txt 2>&1
   ```

5. **Compare Results**:
   ```bash
   grep "Total dofs inside" output_baseline.txt
   grep "Total dofs inside" output_accelerated.txt
   ```

If counts match → HashGrid implementation is correct!
If counts differ → Debug point query implementation.
