# q2p1_sse Launch And Staging

Primary files:

- `tools/e3d_scripts/e3d_start.py`
- `tools/e3d_scripts/e3d_start_yaml.py`
- `applications/q2p1_sse/q2p1_sse_start.py`
- `applications/q2p1_sse/CMakeLists.txt`
- `applications/q2p1_sse/_data_BU/`
- `applications/q2p1_sse/_ianus/`

## Role

The launcher prepares the run directory before the Fortran executable starts.
For the documented workflow, `e3d_start.py` is the main entry point copied into
the build directory by the `q2p1_sse` CMake target.

`e3d_start_yaml.py` is the staged-plan variant. It uses the same setup basis but
requires `-y/--yaml` and executes solver steps from YAML files in
`tools/e3d_scripts`. See [YAML Simulation Plans](yaml_simulation_plans.md).

`e3d_start_yaml.py` also separates the installation folder (executables,
`partitioner`, default templates and plans) from the case folder (`-C/--case`,
default: current directory), so a simulation folder only has to contain the
case-specific inputs. The shared logic lives in `tools/e3d_scripts/e3d_layout.py`;
`e3d_start.py` still assumes that it runs inside the installation folder.

## Main Responsibilities

- Parse command-line options such as `-n`, `-f`, `-a`, `-d`, `-t`, `--short-test`,
  `--die-simulation`, `--do-temperature`, `--mesh-reduction`, `--skip-setup`,
  and `--only-mesh-creation`.
- Locate `setup.e3d` in the project folder, falling back to `Extrud3D.dat`.
- Copy the chosen setup file to `_data/Extrud3D_0.dat`.
- Select the correct `_data/q2p1_param.dat` template from `_data_BU`.
- Copy project OFF/OFF-like geometry files into the run directory.
- Remove stale `_data/meshDir`.
- Run `s3d_mesher`; if no mesh is generated, copy `meshDir` from the case folder.
- Partition `_data/meshDir/file.prj` for `numProcessors - 1` compute ranks.
- Create `_data/Extrud3D.dat` for each angle by appending the active angle to
  `_data/Extrud3D_0.dat`.
- Launch `q2p1_sse`, `q2p1_sse_temp`, or `q2p1_sse_mesh` depending on mode.
- In YAML mode, execute staged `MomentumEquation`, `HeatEquation`, and
  `MaterialDistribution` steps from the selected plan file.

## Data Flow

```text
case folder/setup.e3d
  -> _data/Extrud3D_0.dat
  -> _data/Extrud3D.dat
       appended:
         [E3DSimulationSettings]
         dAlpha
         Periodicity
         nSolutions
         Angle

_data_BU/q2p1_param*.dat
  -> _data/q2p1_param.dat

case folder/*.off
  -> working directory/*.off

s3d_mesher or case folder/meshDir
  -> _data/meshDir/file.prj
  -> partitioner.partition(...)
```

## Modes

- Normal velocity mode launches `q2p1_sse`.
- DIE mode sets `singleAngle = 0`, selects DIE parameter templates, and launches
  the same `q2p1_sse` executable with DIE-specific setup.
- Temperature mode alternates velocity and heat solver runs using
  `q2p1_sse_temp`; see [Temperature Extension](temperature_extension.md).
- Mesh-reduction mode launches `q2p1_sse_mesh` and copies `ReducedMeshDir` back
  into the case folder.
- Short-test mode selects reduced iteration parameter templates.
- YAML mode launches through `e3d_start_yaml.py` and delegates the sequence of
  solver runs to `e3d_xse*.yaml` or `e3d_die*.yaml` templates.

## Build/Staging Contract

`applications/q2p1_sse/CMakeLists.txt` defines the `q2p1_sse` executable and a
`q2p1_sse_stage` target. The staging logic creates runtime directories such as
`_data`, `_mesh`, `_vtk`, `_dump`, `_1D`, `_hist`, `_prot0`, `_prot1`, and
`_RTD`, copies `_data_BU`, `_ianus`, `e3d_start.py`, `conv_check.sh`,
`MG.dat`, `s3d_mesher`, and partitioner support files.

The install logic also includes `e3d_start_yaml.py` and YAML templates matching
`tools/e3d_scripts/e3d_xse*.yaml` for the `q2p1_sse` install and
`tools/e3d_scripts/e3d_die*.yaml` for the `q2p1_gendie` install.
