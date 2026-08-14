# heat Convergence Evaluation

## Entry Points

- `applications/heat/heat.f90`
- `source/src_LinSc/LinSc_transport_extensions.f90`
- `source/src_quadLS/QuadSc_var.f90`
- `source/src_util/ReadExtrud3DParameters.f90`
- `source/postprocessing/solution_io.f90`

## Main Exit Flags

The heat application uses shared flags from `var_QuadScalar`:

- `DivergedSolution`
- `ConvergedSolution`

`heat.f90` exits the time loop when either flag is true. It returns process
status `1` only for `DivergedSolution`.

## Scalar Solver Convergence

Inside `Transport_LinScalar_EWIKON`, each timestep computes a defect criterion:

```text
DefTempCrit = max(RhsTemp * Tracer%prm%defCrit, Tracer%prm%MinDef)
```

The nonlinear scalar loop calls `Solve_General_LinScalar`, rebuilds the defect,
and stops once:

```text
DefTemp <= DefTempCrit
and INL >= Tracer%prm%NLmin
```

The loop is capped by `Tracer%prm%NLmax`. Solver progress is written through
`Protocol_linScalar`.

## Divergence Detection

The linear scalar transport layer marks `DivergedSolution` when scalar defects
become NaN or non-finite. The main program then prints the diverged-solution
status, finalizes MPI, and stops with exit code `1`.

## Heat Run-Mode Convergence

`ReadEWIKONfile` selects the heat run mode:

- PID regulation: `HEAT_RUN_MODE_PID`
- fixed power: `HEAT_RUN_MODE_FIXED`
- no regulation: `HEAT_RUN_MODE_NONE`

Later convergence logic in `LinSc_transport_extensions.f90` sets
`ConvergedSolution` based on the active mode. For fixed-power mode the parser
enables `mySigma%MeltMonitor`, and convergence is taken from
`mySigma%MeltMonitor%Detector%Converged`.

The melt monitor uses:

- `E3DSimulationSettings/MeltMonitorTolerance`
- `E3DSimulationSettings/MeltMonitorLimit`

Invalid or missing values fall back to tolerance `1e-3` and limit `250`.

## Wire Sensor Convergence

Wire segments can define their own convergence detector and regulation limits:

- `ConvergenceCondition`
- `ConvergenceLimit`
- `TemperatureSensorMinRegValue`
- `TemperatureSensorMaxRegValue`
- PID set value and P/I/D constants when `Regulation=PID`

Sensor outputs and regulation state are updated from the heat geometry/sensor
classification created by `updateHeatGeometry`.

## Postprocessing Interaction

`postprocessing_app_heat` checks `ConvergedSolution`. When true, it writes a
special profile output with index `7199` before normal final status handling.

When changing convergence behavior, check both the solver-side flag assignment
and this output convention.
