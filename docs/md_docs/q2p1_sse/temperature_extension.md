# q2p1_sse Temperature Extension

Primary files:

- `applications/q2p1_sse_temp/q2p1_sse_temp.f90`
- `applications/q2p1_sse_temp/app_init.f90`
- `applications/q2p1_sse_temp/CMakeLists.txt`
- `source/src_LinSc/LinSc_transport_extensions.f90`
- `source/src_LinSc/LinSc_main.f90`
- `source/src_LinSc/LinSc_user.f90`
- `source/src_quadLS/QuadSc_transport_extensions.f90`
- `source/OutputProfiles.f90`
- `tools/e3d_scripts/e3d_start.py`
- `tools/e3d_scripts/e3d_start_yaml.py`

## Role

`q2p1_sse_temp` is the heat/temperature extension of the `q2p1_sse` extrusion
workflow. It is not a standalone geometry or momentum application. It reuses the
same staged case setup and consumes dumped fields produced by the momentum run.

Launcher entry points:

- `e3d_start.py --do-temperature` runs velocity/temperature cycles.
- `e3d_start_yaml.py` runs `q2p1_sse_temp` for YAML `HeatEquation` steps.

The executable is installed into both `bin/q2p1_sse` and `bin/q2p1_gendie`.

## High-Level Flow

```text
q2p1_sse or q2p1_sse YAML MomentumEquation step
  -> writes angle-indexed dumped fields

q2p1_sse_temp
  -> reads _data/q2p1_paramT.dat
  -> reads _data/Extrud3D_0.dat
  -> loads velocity/mesh/temperature snapshots by angle
  -> interpolates transient fields over substeps
  -> solves Q1 scalar heat equation
  -> writes temperature dumps, visualization, sensor output, and PID state
```

## Initialization

`q2p1_sse_temp.f90` calls:

```text
init_q2p1_ext(ufile)
```

The temperature `init_q2p1_ext` does the following:

- Sets `ApplicationString` to the SSE temperature module.
- Calls temperature-local `General_init_ext`.
- Initializes Q2 structures with `Init_QuadScalar_Stuctures`.
- Initializes the Q1 scalar transport with `Init_LinScalar`.
- Sets export names to `tmp` and `temp`.
- Allocates temperature-coupled fields such as `Screw`, `Shell`, `Shearrate`,
  `Viscosity`, and optionally `mySegmentIndicator`.
- Allocates `myTransientSolution` arrays for velocity, coordinates, screw/shell
  distance fields, temperature, and optional segment indicators.
- Loads dumped fields from the momentum run for each available angle.
- Expands periodic copies when `myProcess%Periodicity > 1`.
- Optionally subtracts rotational-frame velocity components when
  `mySetup%bRotationalFramOfReference` is enabled.
- Initializes the Q1 scalar temperature field.
- Restores previous PID data via `READ_PID_DATA` for restart runs.

## General_init_ext Differences

The temperature `General_init_ext` is similar in shape to the `q2p1_sse`
initialization, but has several important differences:

- It sets `myDataFile = '_data/q2p1_paramT.dat'` before `GDATNEW("SimPar", 0)`.
- It reads `_data/Extrud3D_0.dat`, not the per-angle `_data/Extrud3D.dat`.
- It does not call `Setup_STL_Segments` in the local file.
- It uses generic `ProlongateCoordinates`, while the momentum app uses
  `ProlongateCoordinates_SSE`.
- It skips rigid-body initialization when restarting from MPI dump format.

The result is a mesh/communication/parametrization setup compatible with the
dumped momentum fields and the Q1 scalar heat equation.

## Transient Field Preparation

The temperature app depends on dumped fields from the momentum solve. Depending
on `myTransientSolution%DumpFormat`, initialization loads:

- velocity fields `v`
- mesh coordinates/deformation fields `d` and `x`
- temperature field `t`
- screw/shell or segment fields `s`
- optional multi-material alpha fields `q`
- viscosity/shear-related field `y`

The data is copied into `myTransientSolution%Velo`, `%Coor`, `%Dist`, `%Shell`,
`%Temp`, and optional `%iSeg`.

This stage is a key coupling boundary: if the momentum app changes dump names,
field contents, angle numbering, or periodicity conventions, `q2p1_sse_temp`
must be updated with it.

## Main Solution Loop

The temperature main loop is nested:

```text
for iRot = 1, nitns:
  for iStep = 0, myProcess%nTimeLevels - 1:
    for iSubStep = 1, myTransientSolution%nTimeSubStep:
      update time state
      TemporalFieldInterpolator(iStep, iSubStep)
      optional Assemble_LinScOperators_XSE
      Transport_LinScalar_XSE(...)
      print_time
      handle_statistics
    angle-level output and dump handling
```

Time setup:

- `dPeriod = 60 / myProcess%Umdr`
- `dTimeStep = dPeriod / myProcess%nTimeLevels`
- `tstep = dTimeStep / myTransientSolution%nTimeSubStep`
- `dtgmv = dPeriod / myProcess%nTimeLevels`

If `mySetup%bConstantMesh` is true, `Assemble_LinScOperators_XSE` is called once
before the loop. Otherwise it is called after each interpolation.

## Field Interpolation And ALE Correction

`TemporalFieldInterpolator(iL, iS)` lives in
`source/src_quadLS/QuadSc_transport_extensions.f90`.

It interpolates between two angle/time-level snapshots:

- `QuadSc%ValU`, `%ValV`, `%ValW`
- mesh coordinates in `mg_mesh%level(NLMAX)%dcorvg`
- `Screw`
- `Shell`
- optional `mySegmentIndicator`

After interpolation it rebuilds dump structures, transfers coordinates to the
coarse/master representation, recreates the density/mass matrix, recomputes
non-Newtonian viscosity, and subtracts ALE mesh velocity from the interpolated
fluid velocity.

## Heat Transport

The heat solve is:

```text
Transport_LinScalar_XSE(Boundary_LinSc_Val_XSE, AddSource_XSE, ufile, inonln_t)
```

Important work inside `Transport_LinScalar_XSE`:

- Builds boundary-condition markers through `Create_Knpr`.
- Sets `thstep = 0.5 * tstep`.
- Stores old temperature solution.
- Assembles scalar advection/heat equation defect through
  `Matdef_LinScalar_XSE`.
- Adds air-cooling, boundary heat flux, and volumetric heat flux contributions.
- Applies Dirichlet boundary conditions.
- Solves with `Solve_General_LinScalar`.
- Updates global `Temperature`.
- Checks for NaN/non-finite temperature state and sets `DivergedSolution`.
- Outputs sensor temperatures.
- Integrates heat rate and outflow temperature.

The operator setup is handled by `Assemble_LinScOperators_XSE`, which creates
rho-cp convection, mass, lumped mass, lambda diffusion, and scalar mass matrices.

## Sensors And PID

Temperature sensor and PID handling is split between:

- `OutputSensorTemperatures` in `source/src_LinSc/LinSc_transport_extensions.f90`
- `PID_controller` in `source/PID.f90`
- `OUTPUT_PID_DATA` and `READ_PID_DATA` in `source/OutputProfiles.f90`

Sensor handling computes average temperature and volume around configured DIE
sensors, updates PID-controlled heat sources, and handles `+STOP`/`-STOP`
sensor logic for heat-source and constant-temperature switching.

PID state is persisted in:

```text
_dump/PID.dmp
```

## Postprocessing And Output

The temperature executable does not call `postprocessing_sse`. Its output path is
more direct:

- At angle `0`, it calls `viz_output_fields_Simple` for `tmp/temp` output.
- For selected periodic angles, it releases temperature dumps:
  - `Release_ListFiles_General(angle, 't')` for list-file dumps.
  - `ReleaseMPIDumpFiles(angle, 't')` for MPI dump format.
- It calls `OUTPUT_PID_DATA` at angle-level output points.
- It writes timing/statistics through `print_time` and `handle_statistics`.

This means temperature postprocessing is coupled to dump/restart conventions and
simple visualization output, rather than the `postprocessing_sse` mechanism used
by the momentum executable.

## Additional Refactoring-Relevant Stages

Besides initialization, main loop, and output, the temperature app has two
stages worth treating explicitly:

- Transient field preparation: loading and arranging dumped momentum fields by
  angle, periodicity, and substep interpolation.
- Heat-control feedback: sensor/PID state updates that modify heat-source or
  constant-temperature behavior during the solve.

