# heat Postprocessing

## Entry Points

- `applications/heat/heat.f90`
- `source/postprocessing/solution_io.f90`
- `source/src_LinSc/LinSc_transport_extensions.f90`
- `applications/heat/sensor_temperature_extraction.sh`

## Per-Step Postprocessing

`heat.f90` calls:

```text
postprocessing_app_heat(dout, inonln_u, inonln_t, ufile)
print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)
handle_statistics(tt0, itns)
```

`postprocessing_app_heat` handles profile output timing and selected dump hooks.

## Profile Output

`postprocessing_app_heat`:

- records output timing statistics on the first timestep
- writes `Output_Profiles(7199)` when `ConvergedSolution` is true
- writes normal `Output_Profiles(iXgmv)` when `dout <= timens`
- enlarges `tstep` and `dtgmv` by `dTimeStepEnlargmentFactor` after scheduled
  output

Some dump calls are intentionally commented out for production mode, including
initial profile output and intermediate heat solution dumps.

## Solver-Embedded Output

`Transport_LinScalar_EWIKON` performs additional heat-specific output before
returning to the main program:

- `IntegrateOutputQuantities(mfile)`
- `CreateSensorOutputs(mfile)`
- protocol lines for conductive, convective, radiative heat fluxes and heat
  source power
- updates `Temperature`
- accumulates `Temperature_AVG` after the first 10 percent of timesteps or time

Sensor and integrated heat quantities are therefore part of the solver step, not
only the postprocessing routine.

## Sensor Extraction Helper

`applications/heat/sensor_temperature_extraction.sh` is installed beside the
application. It extracts sensor-related temperature information from runtime
outputs and `_data/heat.s3d`; check this script when changing protocol formats
or sensor output names.

## Finalization

The heat main program prints one of three final statuses:

- successfully finished before the timestep limit when `ConvergedSolution` is
  true
- successfully finished when final time or timestep limit is reached
- stopped due to diverged solution when `DivergedSolution` is true

The program calls `MPI_Finalize` directly and exits with status `1` on
divergence.
