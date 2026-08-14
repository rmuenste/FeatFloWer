# q2p1_sse Application Index

This index is a navigation layer for the historically named `q2p1_sse`
application family. The executable is used for SSE/TSE-style extrusion cases and
also supports DIE/gendie workflows through setup choices and launcher flags.

The goal of this document tree is to make future restructuring work cheaper:
start here, then open only the subsystem page and source files relevant to the
change.

## Source Boundary

- Main executable: `applications/q2p1_sse/q2p1_sse.f90`
- Application initialization: `applications/q2p1_sse/app_init.f90`
- Temperature executable: `applications/q2p1_sse_temp/q2p1_sse_temp.f90`
- Temperature initialization: `applications/q2p1_sse_temp/app_init.f90`
- Build and staging: `applications/q2p1_sse/CMakeLists.txt`
- Main launcher: `tools/e3d_scripts/e3d_start.py`
- YAML simulation-plan launcher: `tools/e3d_scripts/e3d_start_yaml.py`
- YAML plan templates: `tools/e3d_scripts/e3d_*.yaml`
- Legacy launcher: `applications/q2p1_sse/q2p1_sse_start.py`
- E3D/setup parser: `source/src_util/ReadExtrud3DParameters.f90`
- Shared extrusion state: `source/src_quadLS/QuadSc_Sigma_User.f90`
- Momentum solver entry: `source/src_quadLS/QuadSc_transport_extensions.f90`
- Temperature solver entry: `source/src_LinSc/LinSc_transport_extensions.f90`
- Rheology model entry: `source/src_quadLS/QuadSc_user.f90`
- Output entry: `source/postprocessing/solution_io.f90`

## Documentation Tree

- [Launch And Staging](q2p1_sse/launch_and_staging.md)
  - Python driver
  - `setup.e3d` to `_data/Extrud3D.dat` staging
  - mesher and partitioner steps
  - angle, DIE, mesh-reduction, and temperature modes
- [YAML Simulation Plans](q2p1_sse/yaml_simulation_plans.md)
  - `e3d_start_yaml.py`
  - `options` and `stages`
  - init/main/final execution plans
  - momentum, heat, and material-distribution solver steps
- [Initialization Flow](q2p1_sse/initialization_flow.md)
  - `init_q2p1_ext`
  - `General_init_ext`
  - MPI, parameter, mesh, communication, CGAL/FBM, and restart setup
- [Configuration And State](q2p1_sse/configuration_and_state.md)
  - `ReadS3Dfile`
  - `Sigma_User` global state
  - `mySigma`, `myProcess`, `mySetup`, `myMultiMat`, `myOutput`
- [Geometry Handling](q2p1_sse/geometry_handling.md)
  - SSE/TSE/DIE/XSE machine type selection
  - segment/object parsing
  - STL/OFF indexing
  - boundary parametrization and mesh deformation hooks
- [Rheology](q2p1_sse/rheology.md)
  - material/rheology parsing contract
  - `AlphaViscosityMatModel`
  - automatic timestep estimation
- [Solver Loop](q2p1_sse/solver_loop.md)
  - `Q2P1_SSE` main loop
  - `Transport_q2p1_UxyzP_sse`
  - scalar transport and convergence checks
- [Temperature Extension](q2p1_sse/temperature_extension.md)
  - `q2p1_sse_temp`
  - dumped velocity/mesh/temperature snapshot loading
  - transient field interpolation
  - Q1 heat transport, heat fluxes, sensors, and PID state
- [Output And Finalization](q2p1_sse/output_and_finalization.md)
  - `postprocessing_sse`
  - VTK/GMV and protocol output
  - `_1D`, `_hist`, `_RTD`, `_prot*`, `_dump`
  - final status and exit codes

## High-Level Flow

```text
e3d_start.py
  -> validate command line and project folder
  -> copy setup.e3d or Extrud3D.dat to _data/Extrud3D_0.dat
  -> choose q2p1_param.dat template from _data_BU
  -> copy OFF assets into the run directory
  -> run s3d_mesher or copy project meshDir
  -> partition _data/meshDir/file.prj for numProcessors - 1 compute ranks
  -> for each angle:
       copy _data/Extrud3D_0.dat to _data/Extrud3D.dat
       append angle and shared E3DSimulationSettings
       launch q2p1_sse with mpirun or srun

e3d_start_yaml.py
  -> same staging basis as e3d_start.py
  -> require -y/--yaml plan file
  -> read YAML options and staged solver sequence
  -> validate referenced parameter files
  -> execute init/main/final stages:
       MomentumEquation -> q2p1_sse
       HeatEquation -> q2p1_sse_temp
       MaterialDistribution -> q1_scalar_multimat
  -> archive per-stage protocol output

q2p1_sse.f90
  -> parse -a/--angle, -v/--version, -h/--help
  -> init_q2p1_ext
       -> General_init_ext
            -> GDATNEW("SimPar")
            -> partition reader include
            -> Init_QuadScalar
            -> ReadS3Dfile("_data/Extrud3D.dat")
            -> Setup_STL_Segments
            -> mesh read/refine/communication setup
            -> CGAL/FBM particle and parametrization setup
       -> Init_QuadScalar_Structures_sse
       -> linear scalar or multi-material scalar setup
       -> restart or initial-condition path
  -> time loop:
       Transport_q2p1_UxyzP_sse
       optional Transport_LinScalar
       postprocessing_sse
       print_time and handle_statistics
  -> DetermineIfGoalsWereReached
  -> sim_finalize_sse

q2p1_sse_temp.f90
  -> init_q2p1_ext
       -> General_init_ext
            -> GDATNEW("SimPar") using _data/q2p1_paramT.dat
            -> ReadS3Dfile("_data/Extrud3D_0.dat")
            -> mesh, communication, FBM/CGAL, and parametrization setup
       -> Init_QuadScalar_Stuctures and Init_LinScalar
       -> load dumped q2p1_sse fields by angle:
            velocity, coordinates, screw/shell distances, temperature,
            optional segment indicators and multi-material alpha fields
       -> initialize temperature field and optional PID state
  -> nested rotation/time-level/substep loop:
       TemporalFieldInterpolator
       optional Assemble_LinScOperators_XSE
       Transport_LinScalar_XSE
       print_time and handle_statistics
       angle-level dump/Viz/PID output
  -> MPI_Finalize
```

## Main Architectural Boundaries

- The Python launcher owns run-directory preparation, setup staging, meshing,
  partitioning, and repeated angle launches.
- `ReadS3Dfile` owns interpretation of E3D/setup sections and populates the
  shared extrusion state in `Sigma_User`.
- `app_init.f90` connects the generic FeatFloWer initialization machinery with
  the extrusion-specific setup, mesh deformation, CGAL/FBM geometry, and restart
  modes.
- `q2p1_sse.f90` is intentionally thin: it drives the lifecycle but delegates
  almost all domain behavior to shared modules.
- `q2p1_sse_temp` is the heat/temperature extension of the same extrusion
  workflow. It depends on staged velocity/geometry dumps from the momentum run
  and is launched by `e3d_start.py` temperature mode or YAML `HeatEquation`
  steps.
- Solver, rheology, and output are not application-local implementations; they
  are shared `src_quadLS`, transport, and postprocessing mechanisms configured
  by the SSE state.

## Initial Refactoring Hotspots

- `applications/q2p1_sse/app_init.f90` mixes MPI setup, parameter reading,
  mesh construction, geometry registration, restart handling, and timestep
  estimation.
- `source/src_util/ReadExtrud3DParameters.f90` is the central parser and is large
  enough to deserve section-level documentation before behavioral changes.
- `source/src_quadLS/QuadSc_Sigma_User.f90` is both the state contract and a
  large geometry implementation container.
- `tools/e3d_scripts/e3d_start.py` controls important production behavior and
  should be treated as part of the application, not as a simple helper script.
- `applications/q2p1_sse_temp/app_init.f90` contains the temperature-specific
  coupling to dumped momentum fields and should be treated as part of the SSE
  application family during refactoring.
