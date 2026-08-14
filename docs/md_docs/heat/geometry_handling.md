# heat Geometry Handling

## Entry Points

- `applications/heat/heat_start.py`
- `applications/heat/app_init.f90`
- `source/src_util/ReadExtrud3DParameters.f90`
- `source/src_LinSc/LinSc_main.f90`
- `source/src_mesh/geometry_processing.f90`
- `source/src_mesh/umbrella_smoother.f90`

## Geometry Assets

The launcher copies `*.off` and `*.OFF` files from the selected case folder into
the working directory. `_data/heat.s3d` points each heat segment at one or more
OFF/STL-style geometry files through `ScrewOFF(...)`.

The build/install layout may also provide default heat cases under
`_ianus/HEAT`.

## Parser To Geometry Registration

`General_init_ext` calls:

```text
ReadEWIKONfile("_data/heat.s3d")
Setup_STL_Segments
```

`ReadEWIKONfile` records segment object type, material index, temperature
values, heat-source limits, sensor definitions, and OFF files. Then
`Setup_STL_Segments` registers the parsed geometry assets for distance queries.

## Mesh Source

`heat_start.py` runs:

```text
./s3d_mesher -a heat
```

The solver expects `_data/meshDir/file.prj` and the generated mesh directory. If
the mesher does not create `_data/meshDir`, the launcher copies a `meshDir`
fallback from the case folder.

Inside the solver, `General_init_ext` reads and refines the partitioned coarse
mesh:

```text
readTriCoarse
refineMesh
PARENTCOMM
CREATECOMM
E013_CreateComm_coarse
Create_GlobalNumbering
```

## Heat Object Fields

`InitHeatObjects` allocates nodal heat classification arrays:

- `myHeatObjects%Block`
- `myHeatObjects%Wire`
- `myHeatObjects%Channel`
- `myHeatObjects%Sensor`
- `myHeatObjects%Segment`

`myHeatObjects%Segment = 0` means air/background. Positive segment ids map back
to `mySigma%mySegment(iSeg)`.

## Distance And Classification Update

`InitLinearOperators` calls `updateHeatGeometry(mfile)`.

For `mySigma%cType == "HEAT"`, `updateHeatGeometry` calls:

```text
calcDistanceFunction_heat(...)
calcDistanceFunction_sensor(...)
```

These routines classify fine-level mesh nodes/elements against parsed heat
segments and sensor regions. The classification is later consumed by boundary
conditions, material-dependent matrix assembly, source terms, sensor output, and
convergence checks.

## Important Coupling Points

- `Boundary_LinSc_Val_EWIKON` uses `Tracer%knpr` and
  `myHeatObjects%Segment(i)` to assign boundary temperatures.
- `AddSource_EWIKON` uses `myHeatObjects%Segment(i)` to map a node to segment
  object type, material index, heat-source power, and segment volume.
- `Matdef_LinScalar_EWIKON` assembles material-dependent scalar matrices using
  the EWIKON heat/material state.
- `CreateSensorOutputs` uses sensor classifications to evaluate temperature
  feedback and output.
