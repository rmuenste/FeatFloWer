# PE Library Initialization in FeatFloWer Applications

## Overview

This document describes how FeatFloWer CFD applications integrate with the PE (Physics Engine) library for rigid body dynamics simulations. The PE library provides collision detection, contact resolution, and rigid body motion integration for fluid-structure interaction (FSI) problems.

## Architecture

The integration between FeatFloWer (Fortran) and the PE library (C++) occurs through a carefully designed C/Fortran interface layer that handles:

- MPI communicator conversion between Fortran and C
- Physics world initialization
- Rigid body geometry setup
- Material property configuration
- Simulation parameter loading

## Initialization Flow

### 1. CFD Application Side (Fortran)

The initialization occurs in `applications/*/app_init.f90` within the `General_init_ext()` subroutine:

```fortran
! applications/q2p1_fc_ext/app_init.f90:299-302
call init_fc_rigid_body(myid)
call FBM_GetParticles()
CALL FBM_ScatterParticles()
```

After mesh setup and parallel communication structures are established, the application creates an MPI subcommunicator excluding rank 0 (the master process):

```fortran
! applications/q2p1_fc_ext/app_init.f90:453-464
processRanks(1) = 0
CALL MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_W0, error_indicator)
CALL MPI_GROUP_EXCL(MPI_W0, 1, processRanks, MPI_EX0, error_indicator)
CALL MPI_COMM_CREATE(MPI_COMM_WORLD, MPI_EX0, MPI_Comm_EX0, error_indicator)

#ifdef HAVE_PE
  if (myid .ne. 0) then
    call commf2c_fsi(MPI_COMM_WORLD, MPI_Comm_Ex0, myid)
  end if
#endif
```

**Key Points:**
- Only non-master processes (myid ≠ 0) participate in PE simulation
- Master process (myid = 0) handles I/O and coarse grid operations
- Different interface functions exist for different simulation types

### 2. C++/Fortran Interface Layer

The interface layer (`libs/pe/src/interface/c2f_interface.cpp`) provides multiple entry points based on simulation type:

| Interface Function | Simulation Type | Setup Called |
|-------------------|-----------------|--------------|
| `commf2c_` | Particle sedimentation | `setupParticleBench[Serial]()` |
| `commf2c_fsi_` | FSI benchmark | `setupSpan[Serial]()` → `setupFSIBench[Serial]()` |
| `commf2c_drill_` | Drilling simulation | `setupDrill[Serial]()` |
| `commf2c_lubrication_lab_` | Lubrication experiments | `setupLubricationLab[Serial]()` |
| `commf2c_archimedes_` | Buoyancy tests | `setupArchimedes[Serial]()` |
| `commf2c_kroupa_` | Kroupa benchmark | `setupKroupa[Serial]()` |
| `commf2c_creep_` | Creep flow | `setupCreep[Serial]()` |
| `commf2c_dkt_` | Draft-Kiss-Tumble | `setupDraftKissTumb[Serial]()` |

The interface layer has **conditional compilation** based on PE mode:

#### Parallel PE Mode (HAVE_MPI)

```cpp
// libs/pe/src/interface/c2f_interface.cpp:120-140
extern "C" void commf2c_fsi_(MPI_Fint *Fcomm, MPI_Fint *FcommEx0, int *remoteRank)
{
  int remRank = *remoteRank;

  if(remRank != 0) {
    int rank, size;

    MPI_Comm CcommEx0 = MPI_Comm_f2c(*FcommEx0);
    MPI_Comm_rank (CcommEx0, &rank);
    MPI_Comm_size (CcommEx0, &size);

    if (rank == 1) {
      printf( "%d> C) Configuration FSI bench with %d processes.\n", remRank, size );
    }

    setupSpan(CcommEx0);  // Calls parallel setup from sim_setup.h
  }
}
```

#### Serial PE Mode (PE_SERIAL_MODE)

```cpp
// libs/pe/src/interface/c2f_interface.cpp:326-329
extern "C" void commf2c_fsi_(int *Fcomm, int *FcommEx0, int *remoteRank) {
  // Serial PE mode: Each CFD domain runs independent PE instance
  // Pass CFD rank for unique log filenames (pe<rank>.log)
  pe::setupFSIBenchSerial(*remoteRank);
}
```

## PE Modes: Parallel vs Serial

### Parallel PE Mode (Default)

**Build Configuration:**
```bash
cmake -DUSE_PE=ON ..
make -j8
```

**Characteristics:**
- PE library uses MPI for domain decomposition
- Particles distributed across processes with shadow copies
- Full MPI communication within PE for particle synchronization
- Efficient for many small/medium particles

**Use Cases:**
- Particle-laden flows with many particles (> 20)
- Particles smaller than domain size
- Distributed particle systems

### Serial PE Mode

**Build Configuration:**
```bash
# Step 1: Enable PE library
cmake -DUSE_PE=ON ..

# Step 2: Enable serial mode (recommended two-step approach)
cmake -DUSE_PE_SERIAL_MODE=ON ..

# Build
make -j8
```

**Characteristics:**
- Each CFD domain runs independent serial PE instance
- No MPI within PE library
- All domains maintain full particle information
- Forces synchronized via CFD's MPI layer (bypasses PE MPI)
- No shadow copies or distant process registration

**Use Cases:**
- Few large particles (< 20) that span multiple CFD domains
- Particles larger than CFD domain size (e.g., die geometry)
- Very fine CFD meshes with small domain sizes
- Benchmark configurations with large objects

**Technical Details:**
- Preprocessor flag: `-DPE_SERIAL_MODE`
- Forces synchronized via `COMM_SUMMN` in CFD layer
- Alpha field computation checks all particles (efficient for small counts)
- Deterministic serial PE ensures consistency across domains

## Setup Functions

### Serial PE Setup Functions

Located in `libs/pe/pe/interface/sim_setup_serial.h` as **inline header-only functions**.

#### Example: setupFSIBenchSerial()

```cpp
// libs/pe/pe/interface/sim_setup_serial.h:111-228
inline void setupFSIBenchSerial(int cfd_rank) {
  // 1. Set custom rank for PE logger (unique log file per CFD domain)
  pe::logging::Logger::setCustomRank(cfd_rank);

  // 2. Load configuration from JSON file
  SimulationConfig::loadFromFile("example.json");
  auto &config = SimulationConfig::getInstance();
  config.setCfdRank(cfd_rank);
  const bool isRepresentative = (config.getCfdRank() == 1);

  WorldID world = theWorld();

  // 3. Set gravity from configuration
  world->setGravity(config.getGravity());

  // 4. Configure fluid properties
  real simViscosity(config.getFluidViscosity());
  real simRho(config.getFluidDensity());

  world->setLiquidSolid(true);
  world->setLiquidDensity(simRho);
  world->setViscosity(simViscosity);
  world->setDamping(1.0);

  // 5. CRITICAL: Enable automatic force reset for substepping
  world->setAutoForceReset(true);

  // 6. Configure VTK visualization (if enabled)
  if (isRepresentative && config.getVtk()) {
    unsigned int effectiveVisspacing = config.getVisspacing() * config.getSubsteps();
    vtk::WriterID vtk = vtk::activateWriter(
      "./paraview", effectiveVisspacing, 0,
      config.getTimesteps() * config.getSubsteps(),
      false);
  }

  // 7. Create ground planes
  int idx = 0;
  PlaneID plane = createPlane(idx++, 0.0, 0.0, 1.0, 0.0, granite);
  createPlane(idx++, 0.0, 1.0, 0.0,  0.0, granite);  // +y
  createPlane(idx++, 0.0,-1.0, 0.0, -0.02, granite); // -y

  // 8. Create rigid body (chip/particle)
  real radBench = config.getBenchRadius();
  real rhoParticle(config.getParticleDensity());
  Vec3 position(1.0, 0.01, 0.1275);
  std::string fileName = std::string("span_cm.obj");

  TimeStep::stepsize(config.getStepsize());
  MaterialID chipMat = createMaterial("chip", rhoParticle, 0.01, 0.05, 0.05, 0.2, 80, 100, 10, 11);

  Vec3 chipPos = position;
  TriangleMeshID chip = createTriangleMesh(++idx, chipPos, fileName, chipMat, false, true);

  // 9. Enable DistanceMap acceleration (if CGAL available)
#ifdef PE_USE_CGAL
  chip->enableDistanceMapAcceleration(64, 3);  // spacing, resolution
  const bool distanceMapEnabled = chip->hasDistanceMap();
  const DistanceMap* dm = distanceMapEnabled ? chip->getDistanceMap() : nullptr;
  if (!distanceMapEnabled && isRepresentative) {
    std::cerr << "WARNING: DistanceMap acceleration failed to initialize for chip" << std::endl;
  }
#endif

  // 10. Print simulation setup information
  if (isRepresentative) {
    std::cout << "\n--SIMULATION SETUP--------------------------------------------------------------\n"
              << " Simulation stepsize dt                  = " << TimeStep::size() << "\n"
              << " Fluid viscosity                         = " << simViscosity << "\n"
              << " Fluid density                           = " << simRho << "\n"
              << " Gravity                                 = " << world->getGravity() << "\n"
              << " Triangle mesh file                      = " << fileName << "\n"
              << " VTK output                              = " << (config.getVtk() ? "enabled" : "disabled") << "\n"
              << " Distance map enabled                    = " << (distanceMapEnabled ? "yes" : "no") << "\n";

    if (distanceMapEnabled && dm) {
      std::cout << " Distance map grid size                  = "
                << dm->getNx() << " x " << dm->getNy() << " x " << dm->getNz() << "\n"
                << " Distance map origin                     = ("
                << dm->getOrigin()[0] << ", " << dm->getOrigin()[1] << ", " << dm->getOrigin()[2] << ")\n"
                << " Distance map spacing                    = " << dm->getSpacing() << "\n";
    }

    std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
  }
}
```

#### Available Setup Functions

| Function | Purpose | Source Location |
|----------|---------|-----------------|
| `setupDrillSerial()` | Drill simulation with rotating tool mesh and workpiece | `sim_setup_serial.h:461-599` |
| `setupFSIBenchSerial()` | FSI benchmark with chip geometry | `sim_setup_serial.h:111-228` |
| `setupLubricationLabSerial()` | Lubrication experiments with sphere | `sim_setup_serial.h:388-451` |
| `setupParticleBenchSerial()` | Basic particle sedimentation | `sim_setup_serial.h:24-79` |
| `setupArchimedesSerial()` | Buoyancy test setup | `sim_setup_serial.h:235-254` |
| `setupKroupaSerial()` | Kroupa benchmark setup | `sim_setup_serial.h:261-278` |
| `setupCreepSerial()` | Creep flow configuration | `sim_setup_serial.h:285-302` |
| `setupDraftKissTumbSerial()` | Draft-Kiss-Tumble benchmark | `sim_setup_serial.h:309-326` |
| `setup2x2x2Serial()` | 2x2x2 domain decomposition test | `sim_setup_serial.h:333-343` |
| `setupCylSerial()` | Cylinder flow setup | `sim_setup_serial.h:350-360` |
| `setupDCAVSerial()` | DCAV (Double Concentric Annular Viscometer) | `sim_setup_serial.h:606-616` |

### Parallel PE Setup Functions

Located in `libs/pe/src/interface/sim_setup.cpp` as **compiled functions**.

Similar structure but include MPI-specific domain decomposition and process registration.

## Configuration

### JSON Configuration Files

PE setup functions load simulation parameters from JSON files (typically `example.json`):

```json
{
  "gravity": [0.0, 0.0, -9.807],
  "fluid_viscosity": 1.0e-3,
  "fluid_density": 1000.0,
  "particle_density": 2500.0,
  "bench_radius": 0.0635,
  "stepsize": 0.001,
  "substeps": 1,
  "timesteps": 1000,
  "vtk": true,
  "visspacing": 10,
  "use_checkpointer": false,
  "checkpoint_path": "./checkpoints",
  "pointerspacing": 100
}
```

### SimulationConfig Singleton

Configuration accessed via the singleton pattern:

```cpp
auto &config = SimulationConfig::getInstance();

// Get parameters
Vec3 gravity = config.getGravity();
real viscosity = config.getFluidViscosity();
real density = config.getFluidDensity();
real radius = config.getBenchRadius();
real dt = config.getStepsize();
int substeps = config.getSubsteps();
bool vtkEnabled = config.getVtk();
```

## Time Stepping

### Serial PE Time Stepping

Time stepping with substepping support (`sim_setup_serial.h:634-723`):

```cpp
inline void stepSimulationSerial() {
  WorldID world = theWorld();
  auto &config = SimulationConfig::getInstance();

  static int timestep = 0;

  // Substepping configuration
  int substeps = config.getSubsteps();
  real fullStepSize = config.getStepsize();
  real substepSize = fullStepSize / static_cast<real>(substeps);

  // CRITICAL: Apply external forces (fluid drag) with FULL timestep BEFORE substepping
  for (auto it = theCollisionSystem()->getBodyStorage().begin();
       it != theCollisionSystem()->getBodyStorage().end(); ++it) {
    BodyID body = *it;
    body->applyFluidForces(fullStepSize);  // Uses full dt, not substep dt
  }

  // Set global timestep to substep size for physics integration
  TimeStep::stepsize(substepSize);

  // Execute substeps (for collision handling and gravity)
  for (int istep = 0; istep < substeps; ++istep) {
    // Collision detection and response
    // world->simulationStep(substepSize);  // Commented out in current implementation
  }

  // Restore original timestep size
  TimeStep::stepsize(fullStepSize);

  // Particle diagnostics output (once per main timestep)
  // ... output code ...

  timestep++;
}
```

**Key Points:**
- Fluid forces applied **once** per CFD timestep with full dt
- Physics integration uses substeps for stability
- VTK output frequency adjusted for substepping: `visspacing * substeps`
- Forces reset automatically if `setAutoForceReset(true)`

## Example Application Flow

### q2p1_drill Application

1. **Fortran Initialization** (`applications/q2p1_drill/app_init.f90`):
   ```fortran
   call General_init_ext(79, log_unit)
   ! ... mesh setup, MPI communicators ...
   call init_fc_rigid_body(myid)
   call FBM_GetParticles()
   call FBM_ScatterParticles()
   ! ... later in General_init_ext ...
   if (myid .ne. 0) then
     call commf2c_drill(MPI_COMM_WORLD, MPI_Comm_Ex0, myid)
   end if
   ```

2. **C++ Interface** (`libs/pe/src/interface/c2f_interface.cpp:356-359`):
   ```cpp
   extern "C" void commf2c_drill_(int *Fcomm, int *FcommEx0, int *remoteRank) {
     pe::setupDrillSerial(*remoteRank);
   }
   ```

3. **Setup Function** (`libs/pe/pe/interface/sim_setup_serial.h:461-599`):
   - Loads `example.json` configuration
   - Creates workpiece mesh (`span2_scaled.obj`)
   - Creates drill mesh (`tool_scaled.obj`)
   - Enables DistanceMap acceleration for both meshes
   - Configures materials, gravity, fluid properties
   - Sets up VTK output and timestep

4. **Simulation Loop** (in CFD main loop):
   - CFD computes fluid forces on particles
   - `stepSimulationSerial()` advances rigid body positions
   - Particle positions/velocities sent back to CFD for mesh deformation

## Design Rationale

### Why Header-Only for Serial Mode?

Serial PE setup functions are inline header-only to:
1. Avoid MPI symbol conflicts (serial PE compiled without MPI)
2. Enable compiler optimizations across translation units
3. Simplify linking when PE library built in serial mode
4. Allow easy customization per application without rebuilding PE library

### Why Exclude Rank 0?

The master process (rank 0) is excluded from PE simulation because:
1. It handles coarse grid operations in FeatFloWer's multigrid hierarchy
2. It manages I/O and global synchronization
3. PE particles are only needed on worker processes with fine grids
4. Reduces memory footprint on master process

### Representative Process

One process (typically rank 1) is designated as "representative" to:
- Print diagnostic output once (avoid duplicate messages)
- Generate VTK output files
- Log simulation progress

```cpp
const bool isRepresentative = (config.getCfdRank() == 1);
if (isRepresentative) {
  std::cout << "Simulation setup information..." << std::endl;
}
```

## Troubleshooting

### Parallel Mode Errors for Large Particles

**Error:**
```
ERROR: commf2c_drill_() is not implemented in parallel PE mode
```

**Solution:**
The simulation requires PE_SERIAL_MODE. Rebuild:
```bash
cd build
cmake -DUSE_PE=ON ..
cmake -DUSE_PE_SERIAL_MODE=ON ..
make -j8
```

### Missing DistanceMap Acceleration

**Warning:**
```
WARNING: DistanceMap acceleration failed to initialize for chip
```

**Causes:**
1. PE not built with CGAL support
2. Triangle mesh has issues (non-manifold, degenerate triangles)

**Solution:**
```bash
# Rebuild PE library with CGAL
cd libs/pe/build
cmake -DCGAL=ON ..
make -j8
```

### Configuration File Not Found

**Error:**
```
Could not open config file: example.json
```

**Solution:**
Ensure `example.json` exists in the application run directory. Check `SimulationConfig::loadFromFile()` path.

## Source Code Reference

| Component | File | Description |
|-----------|------|-------------|
| Fortran interface calls | `applications/*/app_init.f90` | CFD application initialization |
| C/Fortran interface layer | `libs/pe/src/interface/c2f_interface.cpp` | MPI communicator conversion and setup dispatch |
| Serial PE setup functions | `libs/pe/pe/interface/sim_setup_serial.h` | Header-only inline setup functions |
| Parallel PE setup functions | `libs/pe/src/interface/sim_setup.cpp` | Compiled MPI-based setup functions |
| Serial PE time stepping | `libs/pe/pe/interface/sim_setup_serial.h:634-723` | `stepSimulationSerial()` function |
| Configuration singleton | `libs/pe/pe/config/SimulationConfig.h` | JSON configuration loader |
| PE interface header | `libs/pe/pe/interface/sim_setup.h` | Parallel mode function declarations |

## Related Documentation

- `libs/pe/CLAUDE.md` - PE library build system and architecture
- `CLAUDE.md` - FeatFloWer build system and PE integration modes
- `docs/rigid_body_dynamics/` - Rigid body dynamics theory and implementation
