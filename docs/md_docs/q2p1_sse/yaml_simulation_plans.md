# q2p1_sse YAML Simulation Plans

Primary files:

- `tools/e3d_scripts/e3d_start_yaml.py`
- `tools/e3d_scripts/e3d_xse.yaml`
- `tools/e3d_scripts/e3d_xse_temp.yaml`
- `tools/e3d_scripts/e3d_die.yaml`
- `tools/e3d_scripts/e3d_die_temp.yaml`
- `tools/e3d_scripts/e3d_die_multi.yaml`
- `tools/e3d_scripts/e3d_die_multi_temp.yaml`

## Role

`e3d_start_yaml.py` is an alternative launcher to `e3d_start.py`. It keeps the
same basic responsibilities for setup staging, meshing, partitioning, MPI launch,
and protocol monitoring, but replaces the fixed simulation loop with a YAML
execution plan.

The YAML mode is useful when a run needs multiple solver phases, for example a
momentum solve followed by heat or material-distribution solves, with different
parameter files per phase.

## Command-Line Contract

The YAML launcher requires the normal project and MPI options plus a YAML plan:

```bash
python ./e3d_start_yaml.py -f _ianus/TSE/Conv -n 5 -y e3d_xse.yaml
python ./e3d_start_yaml.py -f _ianus/DIE/RH_LOC -n 5 -y e3d_die.yaml
```

Inherited options include `-f`, `-n`, `-p`, `-d`, `-a`, `-t`,
`--die-simulation`, `--skip-setup`, `--skip-simulation`,
`--only-mesh-creation`, `--retry-deformation`, and `--use-srun`.

YAML mode adds:

- `-y`, `--yaml`: required path to the YAML execution file.

## Plan Structure

Plans are mappings with optional `options` and required `stages`.

```yaml
options:
  periodicity: 2
  ResolutionLevel: 3

stages:
  init:
    startangle: 0.0
    steps:
      - solver: MomentumEquation
        param_file: _data_BU/q2p1_paramV_0.dat

  main:
    loop: 1
    steps:
      - solver: MomentumEquation
        param_file: _data_BU/q2p1_paramV_1.dat

  final:
    steps:
      - solver: MomentumEquation
        param_file: _data_BU/q2p1_paramV_1.dat
```

Supported stage names:

- `init`
- `main`
- `final`

Supported solver names:

- `MomentumEquation`
- `HeatEquation`
- `MaterialDistribution`

## YAML Options

`apply_yaml_options` reads options from either `options` or `settings`. Some
options may also be supplied as top-level keys for compatibility.

Recognized YAML options include:

- `die-simulation` or `dieSimulation`
- `periodicity`
- `delta-angle` or `deltaAngle`
- `ResolutionLevel`, `resolution-level`, or `resolutionLevel`

Command-line values for periodicity and delta angle take precedence when they
were explicitly provided.

`ResolutionLevel` is applied to staged parameter files. The launcher rewrites
matching entries such as `SimPar@MaxMeshLevel`, `Temper@MGMaxLev`,
`Temper@MGMinLev`, and `Temper@MGMedLev`.

## Stage Execution

`execute_yaml_plan` processes stages in this fixed order:

```text
init -> main -> final
```

For each stage:

- `loop` controls how many times the stage is repeated; default is `1`.
- `steps` must be a list of solver steps.
- `init` may define `startangle` or `start-angle`.
- `startangle` is rejected for DIE simulations and single-angle runs.
- Empty stages are skipped with a warning.

Before execution, `validate_yaml_stage_inputs` verifies that every referenced
`param_file` exists.

## Solver Step Mapping

`MomentumEquation`:

```text
param_file -> _data/q2p1_param.dat
apply ResolutionLevel override
simLoopVelocity -> q2p1_sse
```

`HeatEquation`:

```text
param_file -> _data/q2p1_paramT.dat
apply ResolutionLevel override
launch q2p1_sse_temp
```

The `q2p1_sse_temp` executable is documented in
[Temperature Extension](temperature_extension.md).

`MaterialDistribution`:

```text
param_file -> _data/q2p1_paramAlpha.dat
apply ResolutionLevel override
copy _data_BU/DIE/mesh_names.offs to mesh_names.offs when available
launch q1_scalar_multimat
copy _data/prot.txt to _data/prot_ALPHA.txt
```

## E3D Settings Handling

`write_shared_e3d_settings` rewrites the `[E3DSimulationSettings]` section in
the staged E3D file. It writes `dAlpha`, `Periodicity`, `nSolutions`, and
optional `StartAngle`. For each momentum angle, the usual active angle is still
written into `_data/Extrud3D.dat`.

## Progress And Protocol Output

YAML mode adds YAML-specific status lines to `e3d.log`, including
`YamlStageProgress`, `YamlInStageProgress`, `CurrentYamlSolver`,
`CurrentYamlSolverProgress`, and mirrored solver progress values.

After each stage, `archive_stage_outputs` moves protocol files from `_data` into
`_prot0`, `_prot1`, and so on. The final stage is archived into `_prot_final`.
When a heat step is present in the stage, `_data/prot.txt` is archived as
`prot_TEMP.out`.

## Provided Plan Templates

- `e3d_xse.yaml`: isothermal XSE/SSE/TSE-style momentum plan.
- `e3d_xse_temp.yaml`: XSE/SSE/TSE-style momentum plus heat plan.
- `e3d_die.yaml`: isothermal DIE momentum plan.
- `e3d_die_temp.yaml`: DIE momentum plus heat plan.
- `e3d_die_multi.yaml`: DIE momentum plus material-distribution plan.
- `e3d_die_multi_temp.yaml`: DIE momentum, material-distribution, and heat plan.
