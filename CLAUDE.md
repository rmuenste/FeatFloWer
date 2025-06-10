# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

FeatFloWer uses CMake as its primary build system:

```bash
mkdir build
cd build
cmake ..
make -j8
```

### Key CMake Options
- `USE_PE=ON/OFF`: Enable PE (Physics Engine) library integration for rigid body dynamics
- `USE_CGAL=ON/OFF`: Enable CGAL library for computational geometry
- `USE_MUMPS=ON/OFF`: Enable MUMPS parallel direct solver
- `USE_HYPRE=ON/OFF`: Enable Hypre library for solvers
- `USE_OPENMESH=ON/OFF`: Enable OpenMesh library for mesh processing
- `BUILD_APPLICATIONS=ON/OFF`: Build application executables (default: ON)
- `BUILD_METIS=ON/OFF`: Build METIS library for graph partitioning (default: ON)

### Common Build Commands
```bash
# Basic build
make -j8

# Build specific application
make q2p1_fc2

# Clean build
make clean
```

## Project Architecture

FeatFloWer is a computational fluid dynamics framework with fluid-structure interaction capabilities.

### Core Components

**source/**: Main source code directory
- `src_quadLS/`: Q2P1 finite element solver with level-set methods
  - `QuadSc_force.f90`: Force computation for fictitious boundary method
  - `QuadSc_main.f90`: Main Q2P1 solver routines
- `src_fbm/`: Fictitious Boundary Method (FBM) implementation
- `src_particles/`: Particle simulation and tracking
- `src_mpi/`: MPI parallelization support
- `src_pp3d/`: Post-processing and visualization

**applications/**: Application-specific solvers
- `q2p1_fc2/`: Q2P1 fluid-structure interaction solver with fictitious boundary
- `q2p1_particles/`: Particle-laden flow simulations
- `q2p1_die/`: Die swell flow simulations
- `heat/`: Heat transfer applications
- `laplace/`: Laplace equation solver

**FullC0ntact/**: Rigid body dynamics library (C++)
- Physics engine for particle-particle and particle-wall interactions
- GPU acceleration support via CUDA

### Key Solver Types
- **Q2P1**: Taylor-Hood finite elements (Q2 velocity, P1 pressure)
- **LinSc**: Linear scalar transport equations
- **QuadSc**: Quadratic scalar transport with level-set methods

### Parallelization
- MPI-based domain decomposition
- OpenMP threading for shared memory parallelism
- CUDA support for GPU acceleration (in FullC0ntact)

## Development Workflow

### Adding New Applications
1. Create directory in `applications/`
2. Add `CMakeLists.txt` with application target
3. Include `app_init.f90` for initialization
4. Add main solver file (e.g., `your_app.f90`)

### Fictitious Boundary Method (FBM)
- Uses level-set functions to represent immersed boundaries
- `ALPHA` array tracks fluid (0) vs solid (particle ID) regions
- Force computation in `QuadSc_force.f90` integrates stress over particle surfaces

### Common Development Tasks

**Force Computations**: Modify `source/src_quadLS/QuadSc_force.f90` for custom force models

**Boundary Conditions**: Add to `fbm_vel_bc_include.h` and corresponding source files

**Post-processing**: Extend routines in `source/postprocessing/`

**Particle Integration**: Interface with FullC0ntact library via `dem_query` module

### File Organization Conventions
- `.f90`: Modern Fortran source files
- `.f`: Legacy Fortran 77 source files  
- `_def.f90`: Variable and type definitions
- `_main.f90`: Main program entry points
- `app_init.f90`: Application-specific initialization

## Testing and Validation

Applications include test cases in subdirectories. Many solvers have corresponding Python scripts (`.py`) for parameter setup and post-processing.

### Running Applications
```bash
# Navigate to application directory
cd applications/q2p1_fc2

# Run the solver (typically requires parameter files)
./q2p1_fc2

# Or use provided Python scripts
python q2p1_fc_ext_start.py
```

## Integration Points

- **PE Library**: For rigid body dynamics (when `USE_PE=ON`)
- **CGAL**: Computational geometry operations
- **MUMPS/HYPRE**: Advanced linear solvers
- **MPI**: Distributed memory parallelization
- **OpenMesh**: Mesh processing and manipulation

The codebase supports both academic research and industrial CFD applications with particular strength in fluid-structure interaction and particle-laden flows.