# heat Application Index

This index is a navigation layer for the standalone `applications/heat`
application. The executable solves an EWIKON-style Q1 linear scalar heat problem
on geometry and material data read from `_data/heat.s3d`.

The goal of this document tree is to save search time and tokens: start here,
then open only the focused subsystem page and source files relevant to the
change.

## Source Boundary

- Main executable: `applications/heat/heat.f90`
- Application initialization: `applications/heat/app_init.f90`
- Build, install, and staging: `applications/heat/CMakeLists.txt`
- Python launcher: `applications/heat/heat_start.py`
- Sensor extraction helper: `applications/heat/sensor_temperature_extraction.sh`
- Default case input: `applications/heat/_data/heat.s3d`
- EWIKON/S3D parser: `source/src_util/ReadExtrud3DParameters.f90`
- Shared heat state types: `source/src_util/types.f90`
- Shared extrusion/heat state container: `source/src_quadLS/QuadSc_Sigma_User.f90`
- Linear scalar initialization and geometry update: `source/src_LinSc/LinSc_main.f90`
- Linear scalar EWIKON solver entry: `source/src_LinSc/LinSc_transport_extensions.f90`
- Heat boundary/source callbacks: `source/src_LinSc/LinSc_user.f90`
- Heat distance functions: `source/src_mesh/geometry_processing.f90`
- Restart loading: `source/initialization/app_initialization.f90`,
  `source/OutputProfiles.f90`
- Output entry: `source/postprocessing/solution_io.f90`

## Documentation Tree

- [Launch And Staging](heat/launch_and_staging.md)
  - `heat_start.py`
  - case folder staging
  - meshing, partitioning, MPI/srun launch, install-time runtime files
- [Initialization Flow](heat/initialization_flow.md)
  - `init_q2p1_ext`
  - `General_init_ext`
  - MPI, parameter, mesh, communication, restart, and scalar setup
- [Parsing And State](heat/parsing_and_state.md)
  - `_data/heat.s3d`
  - `ReadEWIKONfile`
  - segments, materials, process values, heat run modes
- [Geometry Handling](heat/geometry_handling.md)
  - OFF/STL segment files
  - mesher contract
  - `Setup_STL_Segments`
  - `calcDistanceFunction_heat` and sensor distance fields
- [Solver Loop](heat/solver_loop.md)
  - `HEAT` main loop
  - `Transport_LinScalar_EWIKON`
  - source terms, boundary fluxes, nonlinear scalar solver
- [Postprocessing](heat/postprocessing.md)
  - `postprocessing_app_heat`
  - profile output, sensor output, integrated heat quantities, dumps
- [Convergence Evaluation](heat/convergence_evaluation.md)
  - scalar defect convergence
  - `DivergedSolution`
  - `ConvergedSolution`
  - heater/sensor and melt-monitor convergence mechanisms

## High-Level Flow

```text
heat_start.py
  -> parse -f/--in-folder, -n/--num-processors, optional --use-srun
  -> copy OFF/STL assets from the case folder into the working directory
  -> copy <case>/heat.s3d to _data/heat.s3d
  -> run s3d_mesher -a heat
  -> use generated _data/meshDir or copy <case>/meshDir fallback
  -> partition _data/meshDir/file.prj for numProcessors - 1 worker ranks
  -> launch heat through mpirun or srun

heat.f90
  -> init_q2p1_ext
       -> General_init_ext
            -> INIT_MPI
            -> GDATNEW("SimPar")
            -> partition reader include
            -> Init_QuadScalar
            -> Init_Die_Handlers
            -> ReadEWIKONfile("_data/heat.s3d")
            -> Setup_STL_Segments
            -> mesh read/refine/communication setup
       -> Init_QuadScalar_Stuctures
       -> Init_LinScalar
       -> istart path:
            0: InitHeatObjects, InitMeshDeform, InitLinearOperators,
               InitCond_LinScalar_EWIKON
            1: InitHeatObjects, load same-level restart,
               InitLinearOperators, SetTracerToLoadedTemperatue
            2/3: currently rejected
  -> time loop:
       update time counters
       Transport_LinScalar_EWIKON
       postprocessing_app_heat
       print_time and handle_statistics
       exit on final time, divergence, or convergence
  -> print final status
  -> MPI_Finalize
```

## Main Architectural Boundaries

- The Python launcher owns run-directory staging, meshing, partitioning, and
  process launch. It does not interpret the heat physics beyond locating
  `heat.s3d`.
- `ReadEWIKONfile` owns interpretation of `_data/heat.s3d` and populates
  `mySigma`, `myMaterials`, `myProcess`, and `mySetup`.
- `app_init.f90` connects the generic FeatFloWer initialization machinery with
  the heat-specific parser, mesh construction, geometry classification, restart
  handling, and scalar setup.
- `heat.f90` is intentionally thin: it drives the lifecycle and delegates heat
  assembly/solve behavior to the shared Q1 linear scalar modules.
- Geometry classification is shared machinery: parsed segment files are
  registered with `Setup_STL_Segments`, then distance/object/sensor fields are
  built in `updateHeatGeometry`.
- Solver, source, boundary, convergence, and output are not local
  implementations; they are shared `src_LinSc`, `src_mesh`, and
  postprocessing mechanisms configured by EWIKON state.

## Initial Refactoring Hotspots

- `applications/heat/app_init.f90` mixes MPI setup, parameter reading, mesh
  construction, EWIKON parsing, geometry registration, restart handling, and
  scalar initialization.
- `source/src_util/ReadExtrud3DParameters.f90` contains the large
  `ReadEWIKONfile` parser and should be treated as the authoritative input
  contract for `_data/heat.s3d`.
- `source/src_LinSc/LinSc_main.f90` contains heat-object allocation and geometry
  classification; changes there can affect boundary/source assembly.
- `source/src_LinSc/LinSc_transport_extensions.f90` contains both scalar solver
  convergence and heat run-mode convergence mechanisms.
- `applications/heat/heat_start.py` is production behavior for staging and
  partitioning, not a disposable helper.
