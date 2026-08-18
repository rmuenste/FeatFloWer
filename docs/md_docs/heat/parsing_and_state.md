# heat Parsing And State

## Entry Points

- `applications/heat/_data/heat.s3d`
- `source/src_util/ReadExtrud3DParameters.f90`
- `source/src_util/types.f90`
- `source/src_quadLS/QuadSc_Sigma_User.f90`

## Parser

`ReadEWIKONfile(cE3Dfile)` is the heat input parser. It uses `iniparser` and
expects `E3DGeometryData/Machine/Type` to be `HEAT`.

The parser fills shared state:

- `mySigma`: machine type, diameter, segment count, segment geometry files,
  heat-object types, sensor definitions, regulation mode, convergence detector,
  heat run mode
- `myMaterials`: heat conductivity, heat capacity, density, slopes, and derived
  thermal diffusivity-like `Alpha`
- `myProcess`: ambient/far-field/cooling/melt-inflow temperatures, heat transfer
  coefficients, workbench thickness
- `mySetup`: convergence estimator flag and optional box mesher controls

## Geometry Sections

The canonical input shape is:

```text
[E3DGeometryData/Machine]
Type=heat
Unit=mm|cm|dm|m
BarrelDiameter=...
NoOfElements=...

[E3DGeometryData/Machine/Element_i]
Type=STL
ObjectType=block|wire|melt
MaterialIndex=...
InitialTemperature=...
VolumetricHeatSourceMin=...
VolumetricHeatSourceMax=...
TemperatureBC=constant|fullconstant|flux|no
Unit=...
ScrewOFF(1)=...
```

Lengths are converted to centimeters internally. Machine and segment units are
validated against `MM`, `CM`, `DM`, and `M`.

## Segment Semantics

`ObjectType` drives later geometry and source handling:

- `BLOCK`: passive solid/block region, optionally with thermal boundary
  conditions
- `WIRE`: heating wire region, with heat-source limits, sensor definition,
  regulation, and convergence detector
- `MELT`: melt/channel region, can receive source handling in
  `AddSource_EWIKON`

`Type=STL` currently reads one or more `screwOFF` entries into
`mySigma%mySegment(i)%OFFfiles`.

## Sensors And Regulation

Wire segments may define:

- `SensorType=COOR`, `OFF`, or `STL`
- `TemperatureSensorCoor`
- `TemperatureSensorRadius`
- `sensorOFF` files for OFF/STL sensors
- `Regulation=SIMPLE`, `PID`, `NONE`, or `FIXED`

`SIMPLE` uses `TemperatureSensorMinRegValue` and
`TemperatureSensorMaxRegValue`. `PID` reads set value and P/I/D constants and
sets `mySigma%bHasPIDRegulation`. `NONE` and `FIXED` set
`mySigma%bHasFixedPowerRegulation`.

The parser derives:

- `HEAT_RUN_MODE_PID` when PID regulation exists
- `HEAT_RUN_MODE_FIXED` when fixed power regulation exists
- `HEAT_RUN_MODE_NONE` otherwise

Melt-monitor convergence is enabled for fixed-power mode and disabled for PID
and no-regulation modes.

## Material Sections

Material sections are indexed from zero:

```text
[E3DMaterialParameters]
NoOfMaterials=...

[E3DMaterialParameters/Mat_i]
HeatConductivity=...
HeatConductivitySlope=...
HeatCapacity=...
HeatCapacitySlope=...
Density=...
DensitySlope=...
```

Segments reference these materials through `MaterialIndex`.

## Simulation Settings

Optional `E3DSimulationSettings` values include:

- `ConvergenceEstimator`
- `SensorPositions`
- `SensorRadius`
- `MeltMonitorTolerance`
- `MeltMonitorLimit`
- `HexMesher`
- `BoxMesherELems`
- `BoxMesherNumberOfELems`
- `MeshResolution`
- `BoxMesherUnit`
- `BoxMesherX`, `BoxMesherY`, `BoxMesherZ`

These settings affect convergence reporting, external sensor candidates, and
box-mesher setup.
