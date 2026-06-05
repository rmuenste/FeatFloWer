# SSE YAML Execution Update for Cluster Automation

## Purpose

This note collects the recent SSE execution updates that should be relevant for Dirk when embedding the new
mechanisms into a cluster-side automation sequence.

The main functional change is the transition from the older `e3d_start.py` launcher to the new
`e3d_start_yaml.py` driver with YAML help-scripts.

## Background

Previously, the production-oriented launcher logic was centered around `e3d_start.py`.
That path was mainly preconfigured for momentum-equation CFD runs.
For XSE/screw cases it could also handle temperature runs, but this was not yet arranged as a clean production
automation path.

This legacy path remains available.
The original `e3d_start.py` keeps its previous functionality, and the previously used parameter files are still
present so the established legacy SSE/gendie track can continue to be used where backward compatibility is needed.

The new setup introduces YAML-driven staged execution plans.
This makes the solver sequence explicit while keeping the runtime call itself a one-liner.
In the FeatFloWer installation, the exported launcher and YAML helper files are available in `tools/e3d_scripts/`.

## Main Improvement

The user-facing execution is now reduced to a single command of the form:

```bash
python3 e3d_start_yaml.py -f <project-folder> -n <ranks> -y <yaml-plan>
```

If the cluster workflow prefers `srun`, the same driver also supports the script option
`-u` / `--use-srun`.

Important clarification for cluster integration:

- In `e3d_start.py` and `e3d_start_yaml.py`, `-u` and `--use-srun` are two spellings for the same
  script flag.
- This flag only switches the launcher from `mpirun ...` to `srun ...`.
- It is not the same as Slurm's own `srun -u` option.
- The driver does not append `-u` to the generated `srun` command; it only prefixes the solver call
  with `srun`.
- When `--use-srun` is active, the script still keeps its internal `-n <ranks>` value for logging
  and partitioning, but it does not translate that value into `srun -n <ranks>`. The actual task
  count must therefore come from the Slurm allocation or outer launch context.

Example:

```bash
python3 e3d_start_yaml.py -f <project-folder> -n 16 -y e3d_xse.yaml -u
```

This means "use the script's `srun` launch path". It does not mean that the script will generate
`srun -u ...`, and it does not guarantee that `16` is forwarded as an `srun -n 16` argument.

The important operational simplification is that the old manual hide/seek around the active
`q2p1_param*.dat` files now happens behind the curtains inside the YAML driver.
The YAML stages point to the backing parameter files, and the driver copies them into the active `_data/`
locations required by the solvers.
This is easier to automate than the older production variant where the surrounding script logic had to manage
the parameter-file switching explicitly.

## DIE Variants to Expose in Automation

All four DIE variants are now available as predefined YAML plans.

### 1. Pure die run

- YAML plan: `e3d_die.yaml`
- Command pattern:

```bash
python3 e3d_start_yaml.py -f <project-folder> -n <ranks> -y e3d_die.yaml
```

- Solver sequence:
  - `init`: `MomentumEquation` with `q2p1_paramV_DIE_0.dat`
  - `main`: `MomentumEquation` with `q2p1_paramV_DIE_2.dat`
  - `final`: `MomentumEquation` with `q2p1_paramV_DIE_2.dat`

This is the isothermal single-material DIE path.

### 2. Die with multimat

- YAML plan: `e3d_die_multi.yaml`
- Command pattern:

```bash
python3 e3d_start_yaml.py -f <project-folder> -n <ranks> -y e3d_die_multi.yaml
```

- Solver sequence:
  - `init`: `MomentumEquation` + `MaterialDistribution`
  - `main`: `MomentumEquation` + `MaterialDistribution`
  - `final`: `MomentumEquation`

This is the isothermal DIE path with alpha/material-distribution transport.
It also resolves the previously used multimaterial track with the long-criticized
parameter-file "switching game" (`Huetchenspiel`), because the required parameter activation is now handled by the
YAML driver instead of by fragile external script logic.
Compared with the legacy track, the material-distribution step is no longer used only once to influence the CFD.
In the current YAML setup it is applied two times, once in `init` and once in `main`.
This gives the material distribution a stronger and more sustained coupling into the flow solution and allows runs
to reach a broader material-filled portion of the flow domain.

### 3. Die with temperature

- YAML plan: `e3d_die_temp.yaml`
- Command pattern:

```bash
python3 e3d_start_yaml.py -f <project-folder> -n <ranks> -y e3d_die_temp.yaml
```

- Solver sequence:
  - `init`: `MomentumEquation` + `HeatEquation`
  - `main`: `MomentumEquation` + `HeatEquation`
  - `final`: `MomentumEquation`

This is the non-isothermal single-material DIE path.

### 4. Die with temperature and multimat

- YAML plan: `e3d_die_multi_temp.yaml`
- Command pattern:

```bash
python3 e3d_start_yaml.py -f <project-folder> -n <ranks> -y e3d_die_multi_temp.yaml
```

- Solver sequence:
  - `init`: `MomentumEquation` + `MaterialDistribution` + `HeatEquation`
  - `main`: `MomentumEquation` + `MaterialDistribution` + `HeatEquation`
  - `final`: `MomentumEquation`

This is the non-isothermal multi-material DIE path.
An important practical point is that the heat solver is no longer limited to the pure CFD path.
With the staged YAML setup it can be combined both with the CFD solver alone and with the multi-material module in
one normalized execution pipeline.
As in the isothermal multi-material YAML track, the material-distribution step is applied two times, once in
`init` and once in `main`, rather than influencing the CFD only once as in the older legacy handling.

## Why This Matters for Dirk's Cluster Automation

The relevant benefit is not only that these combinations exist, but that they are now normalized behind one
entrypoint.

For automation this means:

- The selection of the physical run variant is reduced to choosing the YAML file.
- The concrete solver ordering is encoded in the YAML stage definition instead of in external ad hoc shell logic.
- The parameter-file switching is internalized in the driver.
- The same launcher syntax can be used across the different DIE variants.
- The driver performs a YAML preflight check and stops early if referenced parameter files are missing.
- Stage outputs are archived into legacy-compatible `_prot<stage-index>` folders, which may be useful for structured post-run
  collection on cluster jobs.

In short, the new path is more suitable for robust production sequencing than the older `e3d_start.py`
workflow.
At the same time, it should be seen as an additive upgrade rather than a removal of the old path, since the legacy
launcher and legacy parameterization remain available.

## YAML Stage Progress in `e3d.log`

The YAML driver now also exposes the staged execution state through `e3d.log`.
This is intended for cluster-side monitoring tools that already inspect the `[SimulationStatus]` section and should
not have to parse terminal output to understand where a staged YAML run currently is.

The status block now uses YAML-specific keys:

```ini
YamlStageProgress=<current-stage>/<total-stages>
YamlInStageProgress=<current-solver-in-stage>/<solvers-in-stage>
CurrentYamlSolver=<MomentumEquation|HeatEquation|MaterialDistribution>
CurrentYamlSolverProgress=<solver-local-current>/<solver-local-maximum>
CurrentYamlStatus=<status-text>
```

When the YAML driver is active, transient status keys are intentionally named with `Yaml` in the keyword.
For example, `CurrentStatus` is emitted as `CurrentYamlStatus`, and the final timestamp is emitted as
`YamlFinishingTime`.
Legacy inner-iteration keys such as `CurrentInnerIteration` and `MaxInnerIteration` are not repeated in the YAML
status block because the same information is represented by `CurrentYamlSolverProgress`.

The total stage count is computed from the active YAML plan:

- `init` contributes one stage if it is present and has steps.
- `main` contributes its configured `loop` count.
- `final` contributes one stage if it is present and has steps.

For example, a plan with:

```yaml
stages:
  init:
    steps:
      - solver: MomentumEquation
      - solver: MaterialDistribution
      - solver: HeatEquation

  main:
    loop: 2
    steps:
      - solver: MomentumEquation
      - solver: MaterialDistribution
      - solver: HeatEquation

  final:
    steps:
      - solver: MomentumEquation
```

has four YAML stages in total:

```text
1/4 = init
2/4 = main loop 1
3/4 = main loop 2
4/4 = final
```

Inside each stage, `YamlInStageProgress` tracks normalized stage work, not just the raw solver index.
This matters for XSE/screw runs because one `MomentumEquation` step can expand to multiple angular positions.
If a momentum solve needs `L` angular positions, it contributes `L` in-stage progress units.
`HeatEquation` and `MaterialDistribution` each contribute one in-stage progress unit.

For an `init` or `main` stage with `MomentumEquation + MaterialDistribution + HeatEquation`, the denominator is
therefore `L + 2`.
If an XSE run needs `18` angular positions, the momentum solve advances through `1/20`, `2/20`, ... `18/20`,
the material-distribution step is `19/20`, and the heat step is `20/20`.
For a DIE single-angle run, `L = 1`, so the same stage naturally becomes `1/3`, `2/3`, and `3/3`.

The solver-local progress is sampled from the solver protocol file where possible.
For momentum and heat solves, the driver watches `_data/prot.txt` and parses the existing `itns: current/max`
lines.
For the material-distribution solve, the driver watches the Q1 scalar output and maps the reported
`VolumeFractions[%]` filling degree to an integer progress value against `100`.

Representative `e3d.log` status snapshots for the example pipeline are shown below.
The surrounding fixed header lines such as `PathToE3DFile`, `NumOfCpus`, and `YamlStartingTime` are omitted here for
brevity.
In YAML mode, the legacy header entries `NumOfAnglePositions` and `Periodicity` are not emitted; angular progress is
represented through `YamlInStageProgress`, `CurrentYamlAngleIteration`, and `MaxYamlAngleIteration`.

During the `init` momentum solve:

```ini
YamlStageProgress=1/4
YamlInStageProgress=4/20
CurrentYamlSolver=MomentumEquation
CurrentYamlSolverProgress=4/12
CurrentYamlAngleIteration=4
MaxYamlAngleIteration=18
CurrentYamlStatus=running Momentum Solver
```

During the first `main` material-distribution solve:

```ini
YamlStageProgress=2/4
YamlInStageProgress=19/20
CurrentYamlSolver=MaterialDistribution
CurrentYamlSolverProgress=67/100
CurrentYamlMaterialFillingDegree=67
MaxYamlMaterialFillingDegree=100
CurrentYamlStatus=running Material Distribution Solver
```

During the second `main` heat solve:

```ini
YamlStageProgress=3/4
YamlInStageProgress=20/20
CurrentYamlSolver=HeatEquation
CurrentYamlSolverProgress=6/10
CurrentYamlStatus=running Heat Solver
```

At the end of the `final` momentum stage:

```ini
YamlStageProgress=4/4
YamlInStageProgress=18/18
CurrentYamlSolver=MomentumEquation
CurrentYamlSolverProgress=12/12
CurrentYamlStatus=finished
YamlFinishingTime=...
```

This keeps YAML-mode monitoring unambiguous: automation can restrict itself to keys containing `Yaml` and still
distinguish between `init`, repeated `main` stages, and `final`.

## Constant-Mesh Heat Optimization

There is one further runtime-relevant improvement that should be visible to automation users.

If the setup provides:

```ini
[E3DSimulationSettings]
ConstantMesh = yes
```

then the code detects that the mesh is constant and treats the heat transport accordingly.
In that case the matrix set for the heat solve is built once initially and then reused, instead of being rebuilt
repeatedly during the run.

Operationally, this matters because constant-mesh SSE-TEMP style runs avoid unnecessary matrix regeneration while
still using the same staged execution framework.

## Relation to Screw/XSE Cases

The same YAML concept also exists for screw/XSE cases:

- `e3d_xse.yaml`
- `e3d_xse_temp.yaml`

So the DIE rollout is not a one-off mechanism.
It fits into a broader staged-launch concept for SSE-related executions.

## Optional Extension: ROUND Geometry Keyword

There is a second, independent update that may be relevant for automation as an optional feature flag.

If the setup contains:

```ini
[E3DGeometryData/Machine]
GeometryType = ROUND
```

then the new-generation mesher reacts to this keyword.
That mesh generation step happens in the separate meshing module, not inside FFer itself.

The intended workflow is:

1. The geometry setup advertises `GeometryType = ROUND`.
2. The new mesher creates the corresponding round mesh.
3. FFer then prolongates that mesh to higher levels using cylindrical parametrization for the higher-level
   nodes/vertices.

This means the keyword is not just descriptive.
It activates a geometry-dependent path spanning meshing and later in-solver prolongation.

## Boundary Between Mesher and FFer

To avoid confusion in the handoff:

- Mesh creation for the round setup is handled by the separate, recently rolled-out meshing module.
- The higher-level coordinate prolongation is handled on the FFer side.
- The FFer prolongation uses round/cylindrical information when building the higher mesh levels.

So Dirk's automation may need to treat `GeometryType = ROUND` as a keyword-level switch that affects both the
meshing stage and the subsequent simulation startup path.

## Suggested Handoff Message

The practical message to Dirk is:

The SSE execution path now supports four predefined YAML-driven DIE pipelines
(`die`, `die_multi`, `die_temp`, `die_multi_temp`) through a single launcher,
`e3d_start_yaml.py`.
This reduces automation to selecting the appropriate YAML plan while the driver handles internal parameter-file
activation and stage sequencing, including the former multimat parameter-file switching track.
The heat solver can now be embedded both in pure CFD runs and in combined CFD plus multi-material runs through the
same mechanism.
If `ConstantMesh = yes` is provided, the heat path also avoids repeated matrix rebuilds by reusing the initially
assembled transport matrices.
In addition, `GeometryType = ROUND` can optionally activate the newer round-mesh plus cylindrical-prolongation
workflow, with mesh creation in the separate mesher and higher-level coordinate prolongation in FFer.
The original `e3d_start.py` path and the older parameter files remain available as a legacy-compatible fallback for
the traditional SSE/gendie workflow.
