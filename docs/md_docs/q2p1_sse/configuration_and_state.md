# q2p1_sse Configuration And State

Primary files:

- `source/src_util/ReadExtrud3DParameters.f90`
- `source/src_quadLS/QuadSc_Sigma_User.f90`
- `tools/e3d_scripts/e3d_start.py`
- `applications/q2p1_sse/_data_BU/`

## Setup File Path

The common path is:

```text
case folder/setup.e3d
  -> _data/Extrud3D_0.dat
  -> _data/Extrud3D.dat
  -> ReadS3Dfile("_data/Extrud3D.dat")
```

The parser name is historical: `ReadS3Dfile` reads the staged E3D-style file.

## Parser Role

`ReadS3Dfile` reads the INI-like setup through `iniparser` and populates global
state from `Sigma_User` plus related `var_QuadScalar` state.

Important target state:

- `mySigma`: machine type, geometry dimensions, segment definitions, material
  count, sensor settings, thermal run mode, and geometry flags.
- `myProcess`: angle, rotation direction, speed, temperature, extrusion speed,
  process phase, and other runtime process values.
- `mySetup`: numerical and workflow options such as automatic timestep control.
- `myMultiMat`: material list and rheology descriptions.
- `myThermodyn` and `myMaterials`: thermal/material data.
- `myOutput`: output-related options.
- `myTransientSolution`: restart/dump format behavior.

## Machine-Type Gate

`ReadS3Dfile` accepts machine types including:

- `SSE`
- `TSE`
- `DIE`
- `XSE`
- `NETZSCH`

The `SoftwareRelease` value in `Sigma_User` is checked against the configured
machine type unless the release is `XSE`.

## Parameter Sources

There are two major parameter channels:

- `_data/q2p1_param.dat` is selected from `_data_BU` by the launcher and read
  by `GDATNEW("SimPar", 0)`.
- `_data/Extrud3D.dat` is read by `ReadS3Dfile` and describes machine geometry,
  process settings, materials, rheology, output, and E3D simulation settings.

This split is important for refactoring: numerical solver parameters and
extruder/domain parameters are staged and parsed by different mechanisms.

