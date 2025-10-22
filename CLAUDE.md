# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

FeatFloWer uses CMake (minimum version 3.18) as the primary build system for this parallel multigrid Q2/P1 computational fluid dynamics framework.

### Basic Build Commands

```bash
# Initialize submodules (required)
git submodule update --init --recursive

# Standard build process
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON ..
make -j8

# Show available build configurations
cmake -DSHOW_BUILD_IDS=ON ..
```

### Essential CMake Options

- `-DCMAKE_BUILD_TYPE=Release|Debug`: Build type (default: Release)
- `-DBUILD_APPLICATIONS=ON|OFF`: Build application suite (default: ON)
- `-DUSE_CGAL=ON`: Enable computational geometry features
- `-DUSE_MUMPS=ON`: Enable parallel direct solver (Intel compiler required)
- `-DUSE_HYPRE=ON`: Enable scalable linear solvers
- `-DUSE_PE=ON`: Enable rigid body physics engine
- `-DUSE_PE_SERIAL_MODE=ON`: Use serial PE mode for large particles (requires `-DUSE_PE=ON`)
- `-DUSE_BOOST=ON`: Enable Boost C++ libraries
- `-DUSE_OPENMESH=ON`: Enable mesh processing capabilities

### Build IDs and Compiler Configurations

The build system supports compiler-specific optimizations through build IDs:
- Intel: `nehalem-linux-intel-release`, `xeon-linux-intel-release`
- GCC: Various architecture-specific configurations
- Use `./configure --build-ids` to list all available build IDs

### Dependencies

**Required:**
- C/C++ and Fortran compilers
- OpenMP (always required)
- Git 2.17+ for submodule management
- MPI (for parallel execution which is the only realistic way to run a 3D cfd simulation)

**Optional but Important:**
- CGAL (computational geometry)
- MUMPS (requires Intel compiler)
- Hypre (algebraic multigrid solvers)
- Boost C++ (≥1.68.0, required for CGAL)

## Code Architecture

### Main Components

**Core Navier-Stokes Solver (`source/src_pp3d/`):**
- Q2/P1 finite element discretization
- Multigrid preconditioning and parallel solvers
- MPI parallelization with domain decomposition

**Transport Equations:**
- `src_LinSc/`: Linear scalar transport
- `src_quadLS/`: Quadratic/linear scalar systems
- `src_PLin/`: Particle-based transport methods

**Fluid-Structure Interaction:**
- `src_fbm/`: Fictitious Boundary Method
- `src_particles/`: Lagrangian particle methods
- `FullC0ntact/`: Rigid body dynamics engine

**Specialized Physics:**
- `src_visco/`: Viscoelastic and non-Newtonian flows
- `src_mesh/`: Mesh adaptation and geometry processing

### Application Structure (`applications/`)

Applications follow the pattern `q2p1_*` for Navier-Stokes solvers:

**Core Flow Applications:**
- `q2p1_devel/`: General development solver
- `q2p1_cc/`: Coupled convection-diffusion
- `q2p1_die/`: Die swell extrusion simulation
- `q2p1_sse/`: Single screw extruder

**Fluid-Structure Interaction:**
- `q2p1_fc_ext/`: Extended fictitious boundary method
- `q2p1_fsi_*`: Various FSI benchmarks and applications
- `q2p1_particles/`: Eulerian-Lagrangian coupling

**Specialized Applications:**
- `heat/`: Heat equation solver
- `q1_scalar/`: General scalar transport
- Various rheological and biomedical flow applications

### FullC0ntact Physics Engine

Independent rigid body dynamics framework (`FullC0ntact/`) with:
- Collision detection and response systems
- DEM (Discrete Element Method) capabilities  
- GPU acceleration support
- Integration with FeatFloWer through fictitious boundary methods

See `libs/pe/CLAUDE.md` for detailed PE library information.

### PE Integration Modes

FeatFloWer supports three modes for rigid body physics, configured at build time with CMake options:

#### No PE (`USE_PE=OFF`, default)
- No rigid body physics
- Pure fluid simulation only
- Smallest binary size and build time

#### Parallel PE (`USE_PE=ON`, default when PE enabled)
- Full MPI parallelization within PE library
- Distributed shadow copies across domains
- **Best for:** Many small/medium particles distributed across domains
- **Characteristics:**
  - Particles managed by nearest domain
  - Shadow copies for particles near domain boundaries
  - Standard domain decomposition approach

#### Serial PE (`USE_PE=ON` + `USE_PE_SERIAL_MODE=ON`)
- Each CFD domain runs independent serial PE instance
- Forces synchronized via CFD's MPI layer (bypasses PE MPI)
- **Best for:** Few large particles (< 20) that span multiple domains
- **Characteristics:**
  - All domains maintain full particle information
  - No shadow copies or distant process registration
  - Each domain checks all particles for α field computation
  - Avoids "distant process registration" MPI complexity
- **Build command (recommended two-step approach):**
  ```bash
  # Step 1: Enable PE library
  cmake -DUSE_PE=ON ..

  # Step 2: Enable serial mode
  cmake -DUSE_PE_SERIAL_MODE=ON ..

  # Build
  make -j8
  ```

  **Why two steps?** The `cmake_dependent_option` mechanism ensures `USE_PE_SERIAL_MODE` is only available when `USE_PE=ON`. The two-step configuration guarantees proper dependency resolution and prevents the PE library from being built with MPI support in serial mode.

**When to use Serial Mode:**
- Particles larger than domain size (e.g., die geometry in extrusion)
- Fewer than ~20 rigid bodies total
- Very fine mesh with small domain sizes
- Benchmark configurations with large objects

**Implementation details:**
- Preprocessor flag: `-DPE_SERIAL_MODE`
- Forces still synchronized via `COMM_SUMMN` in CFD layer
- Alpha field computation checks all particles (efficient for small counts)
- Deterministic serial PE ensures consistency across domains

## Development Workflow

### Running Applications

```bash
# Sequential execution
./applications/q2p1_devel/q2p1_devel

# Parallel execution with MPI
mpirun -np 4 ./applications/q2p1_devel/q2p1_devel

# Mesh partitioning for parallel runs
./partitioner mesh.tri 4
```

### Configuration Files

- Parameter files: `_data/q2p1_param.dat` (solver configuration)
- Particle configs: `_data/config_particles.dat`
- Mesh files: Various formats (.tri, .off, .obj)
- Physics setups: JSON files for rigid body configurations

### Testing and Validation

```bash
# Enable testing
cmake -DBUILD_TESTING=ON ..
make test

# Individual application testing
cd applications/q2p1_devel
ctest
```

### Common Build Patterns

```bash
# Full-featured build for research/development (Parallel PE)
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_CGAL=ON \
      -DUSE_MUMPS=ON \
      -DUSE_HYPRE=ON \
      -DUSE_PE=ON \
      -DBUILD_APPLICATIONS=ON ..

# Build with Serial PE mode (for large particles)
# Two-step configuration recommended:
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_PE=ON \
      -DBUILD_APPLICATIONS=ON ..
cmake -DUSE_PE_SERIAL_MODE=ON ..

# Minimal CFD-only build (no rigid bodies)
cmake -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_APPLICATIONS=ON ..
```

## Key File Types and Locations

- **Source files**: Fortran 90/95 (`.f90`), Fortran 77 (`.f`), C++ (`.cpp`, `.h`)
- **Mesh formats**: `.tri` (triangular), `.off` (object file format), `.obj`
- **Output formats**: VTU/VTK for ParaView visualization
- **Configuration**: `.dat` parameter files, `.json` for physics setups

## Submodule Management

This repository uses git submodules extensively:

```bash
# Initialize all submodules (essential)
git submodule update --init --recursive

# Update submodules to latest versions
git submodule update --remote --recursive
```

Major submodules include external libraries (LAPACK, METIS, CGAL, etc.) in `extern/libraries/`.

## Special Notes

- **MUMPS solver**: Requires Intel compiler and MKL
- **CGAL integration**: Automatically enables Boost dependency
- **MPI builds**: Use `MPI_VENDOR="openmpi"` for OpenMPI
- **PE library**: Has independent build system, see `libs/pe/CLAUDE.md`
- **Memory requirements**: Large 3D problems may need significant RAM for multigrid hierarchy