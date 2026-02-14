# PE Library Usage with FeatFloWer

## Purpose
The PE (Physics Engine) library provides rigid-body motion, collision detection, and contact
resolution for particle- or body-laden CFD setups. FeatFloWer couples to PE to advance immersed
particles and feed forces back into the fluid solver.

## When to Enable PE
- Particle sedimentation or slurry transport
- Rigid obstacles that move/deform the flow (FSI benchmarks, drilling, lubrication lab, Archimedes)
- Draft-Kiss-Tumble, Kroupa, creep, or custom particle cases
- Not needed for purely single-phase CFD

## Build and Configuration
- Enable PE: `cmake -DUSE_PE=ON ..`
- Serial mode for large or few particles spanning domains: `cmake -DUSE_PE_SERIAL_MODE=ON ..`
- Flags injected: `-DHAVE_PE` (always), `-DPE_SERIAL_MODE` (serial builds); links `pe_static`
- Keep `BUILD_APPLICATIONS=ON` to compile the PE-enabled front ends
- Optional: `-DUSE_CGAL=ON` to accelerate distance maps for complex meshes

## Runtime Coupling Flow
1. Fortran app sets up CFD mesh and MPI, excludes rank 0 from PE, and creates a worker
   subcommunicator.
2. Application calls the appropriate `commf2c_*` hook (e.g., `commf2c_fsi`, `commf2c_drill`) after
   fluid data structures are ready.
3. The C/Fortran interface converts communicators and dispatches to PE setup routines
   (parallel: `sim_setup.cpp`; serial: `sim_setup_serial.h`).
4. CFD loop supplies fluid forces to PE; PE advances rigid bodies and returns updated positions and
   velocities for mesh deformation or force accumulation.

## Picking a Mode
- Parallel (default): many small/medium particles; PE uses MPI domain decomposition and shadow
  copies; best for particle counts > O(20).
- Serial (`USE_PE_SERIAL_MODE=ON`): few large particles or objects that span domains; each CFD rank
  runs its own PE instance and synchronizes forces via CFD MPI; avoids PE-internal MPI.

## Configuration Files
- PE loads JSON configs (e.g., `example.json`) at setup time via `SimulationConfig`.
- Typical keys: gravity, fluid density/viscosity, particle density/size, timestep/substeps, VTK
  output cadence, mesh filenames for rigid bodies.
- Place the JSON in the run directory or supply via your launcher script.

## Running Applications (examples)
- FSI benchmark: `mpirun -np <ranks> ./applications/q2p1_fsi_bench/q2p1_fsi_bench <param>`
- Lubrication lab: `mpirun -np <ranks> ./applications/q2p1_lubrication_lab/q2p1_lubrication_lab <param>`
- Drill: `mpirun -np <ranks> ./applications/q2p1_drill/q2p1_drill <param>`
- Ensure the parameter file points to the particle JSON and meshes, and that partitions match your
  rank count.

## Data Exchange and Output
- Forces: CFD computes fluid loads; PE applies them once per CFD timestep (full dt) before PE
  substepping.
- Kinematics: PE returns positions/velocities each CFD step; used for mesh motion and diagnostics.
- Logging: one representative rank (usually CFD rank 1) prints setup info; serial mode tags logs
  with CFD rank for uniqueness.
- Visualization: VTK output controlled by JSON `vtk` and `visspacing` settings.

## Troubleshooting Checklist
- Build mismatch: ensure `USE_PE`/`USE_PE_SERIAL_MODE` flags align with the intended mode.
- Communicator errors: verify subcommunicator excludes rank 0 and is non-null before calling
  `commf2c_*`.
- Missing config: confirm the JSON path is reachable from the run directory.
- Distance map warnings: enable CGAL and check rigid mesh quality if acceleration fails.
- Parallel-only hooks: some setups require serial mode (e.g., large-particle drill); rebuild with
  `USE_PE_SERIAL_MODE=ON` if you hit “not implemented in parallel PE mode”.

## Further Reading
- `docs/md_docs/pe_initialization.md` — detailed call flow and per-hook notes
- `libs/pe/CLAUDE.md` — PE build system and architecture
