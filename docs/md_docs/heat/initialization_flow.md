# heat Initialization Flow

## Entry Points

- `applications/heat/heat.f90`
- `applications/heat/app_init.f90`
- `source/initialization/app_initialization.f90`
- `source/OutputProfiles.f90`
- `source/src_LinSc/LinSc_main.f90`

## Main Program Initialization

`heat.f90` calls:

```text
init_q2p1_ext(ufile)
ZTIME
set first output time dout
```

No command-line parsing happens in the Fortran executable. The Python launcher
prepares `_data/heat.s3d`, mesh files, and partition files before `heat` starts.

## `init_q2p1_ext`

`applications/heat/app_init.f90:init_q2p1_ext` performs the heat application
setup:

```text
General_init_ext(79, log_unit)
Init_QuadScalar_Stuctures
Init_LinScalar
select istart path
```

The linear scalar is the temperature unknown. The Q2/P1 structures are still
initialized because shared infrastructure, mesh geometry storage, auxiliary
fields, and postprocessing paths use the broader FeatFloWer state.

## `General_init_ext`

`General_init_ext` performs the generic parallel and mesh setup plus the
heat-specific parser hook:

```text
INIT_MPI
GDATNEW("SimPar")
optional FindNodes
include partition reader
Init_QuadScalar
Init_Die_Handlers
if _data/heat.s3d exists:
    ReadEWIKONfile("_data/heat.s3d")
    Setup_STL_Segments
readTriCoarse
refineMesh
PARENTCOMM and CREATECOMM
E013_CreateComm_coarse
Create_GlobalNumbering
```

`ReadEWIKONfile` must run before heat object geometry is classified because it
allocates and fills the segment, material, process, and heat run-mode state.

## Start Modes

`istart = 0` starts from the initial EWIKON configuration:

```text
CreateDumpStructures
InitHeatObjects
InitMeshDeform
InitLinearOperators
InitCond_LinScalar_EWIKON
```

`InitLinearOperators` calls `updateHeatGeometry`, which classifies mesh nodes
against heat segments and sensors before boundary/source assembly.

`istart = 1` restarts from a same-level solution:

```text
InitHeatObjects
init_sol_same_level_heat(CSTART)
InitLinearOperators
SetTracerToLoadedTemperatue
```

`init_sol_same_level_heat` loads the restart through `SolFromFile_heat`, restores
coordinates from auxiliary Q2 arrays, and exchanges coarse-level node values.

`istart = 2` and `istart = 3` currently stop with a not-supported message.

## Scalar Initial Conditions

`InitCond_LinScalar_EWIKON`:

- temporarily lifts `NLMAX` to access the finest scalar level
- creates scalar boundary markers with `Create_Knpr(LinSc_Knpr)`
- calls `LinSc_InitCond_EWIKON`
- applies `Boundary_LinSc_Val_EWIKON`

For restarts, `SetTracerToLoadedTemperatue` copies loaded `Temperature` into
`Tracer%Val(NLMAX)%x` and reapplies the EWIKON boundary callback.
