# heat Launch And Staging

## Entry Points

- `applications/heat/heat_start.py`
- `applications/heat/CMakeLists.txt`
- `applications/heat/_data/heat.s3d`
- `applications/heat/sensor_temperature_extraction.sh`

## Launcher Contract

`heat_start.py` takes an input folder through `-f/--in-folder`, a total MPI
rank count (at least two) through `-n/--num-processors`, and an optional case
directory through `-C/--case` (default: invocation directory).
The input folder must contain `heat.s3d` and its case-specific
`sampleRigidBody.xml` at the same level. It may also contain a `meshDir`
fallback. Geometry is read directly from the input directory or absolute paths.

The installation is the launcher directory, optionally overridden by
`FF_HEAT_HOME`. It holds executables, the partitioner, and shipped defaults.
The case holds runtime files and outputs. Relative `-C` paths are interpreted
against the invocation directory. `-f` is resolved before changing directories:
invocation directory first, then case, then installation. Absolute paths are
used directly. Relative scheduler host/rank-file paths are also made absolute
before changing directories.

```bash
python3 /opt/ff/bin/heat/heat_start.py \
  -C /scratch/heat_01 -f /scratch/heat_01/input -n 4
# Equivalent from the case:
cd /scratch/heat_01
python3 /opt/ff/bin/heat/heat_start.py -f input -n 4
```

Missing runtime directories are created. Missing `_data/MG.dat`,
`_data/q2p1_param.dat`, and `start/data.TXT` are seeded from the installation.
Existing case copies of those defaults are preserved; parameter selection is
reported. The input `sampleRigidBody.xml` is case-specific and is staged to
`<case>/start/sampleRigidBody.xml` on every launch, replacing a stale runtime
copy, with absolute boundary geometry references. To tune a case, edit its
`_data/q2p1_param.dat`.

The launcher performs these steps:

```text
validate paths and rank count; prepare and enter the case directory
generate _data/heat.s3d with absolute segment and sensor geometry references
generate start/sampleRigidBody.xml with absolute boundary geometry references
leave input configurations and geometry untouched
remove old _data/meshDir
run <installation>/s3d_mesher -a heat
if no generated meshDir exists, copy <input>/meshDir
require _data/meshDir/file.prj
partition _data/meshDir/file.prj for numProcessors - 1 worker ranks
launch <installation>/heat with mpirun or srun
retain geometry, configuration, and outputs even on solver failure
```

Rank 0 is the master/control rank, so the partitioner receives
`numProcessors - 1` worker partitions.

## Partition Format

The launcher reads `_data/q2p1_param.dat` and checks
`SimPar@PartitionFormat`. Valid values are `legacy` and `json`; missing or
invalid values fall back to `legacy`.

Node grouping is inferred from scheduler environment variables first
(`SLURM_STEP_NUM_NODES`, `SLURM_JOB_NUM_NODES`, `SLURM_JOB_NODELIST`) and then
from host/rank files. If no node count is found, one node group is assumed.

## Launch Modes

- Default: `mpirun -np <numProcessors> <installation>/heat`
- `-u/--use-srun`: `srun <installation>/heat`

The launcher propagates the solver exit code. A nonzero mesher exit is reported;
the existing mesh/fallback check still runs. A generated directory without
`file.prj` is rejected. Runtime symlinks and input/runtime mesh overlap are
rejected before mesh removal; keep supplied meshes under `input/meshDir`.

Geometry is never copied or deleted. Relative references in input configurations
are relative to the input directory (`-f`), not the runtime case (`-C`). For
example, `steel.off` becomes `/path/to/input/steel.off` in the generated file.
Nested relative references and existing absolute paths work the same way.
Both generated configurations are refreshed each run; originals remain unchanged.
Use a separate input directory, not the generated `_data/heat.s3d` as input.
All MPI ranks must retain access to the input geometry. See
[Geometry Handling](geometry_handling.md) for native path length/character limits.

Migration: references formerly relative to the runtime case must become
input-relative or absolute. Old OFF copies already in a run directory are left
untouched, but are no longer selected by bare input-relative references.

## Build And Install Staging

`applications/heat/CMakeLists.txt`:

- builds the `heat` executable from `app_init.f90` and `heat.f90`
- depends on `s3d_mesher`
- copies `heat_start.py` and `sensor_temperature_extraction.sh` into the build
  directory
- installs/stages the shared `e3d_layout.py` helper
- installs `start`, `_data`, `_mesh`, `_vtk`, `_dump`, the partitioner Python
  package, scripts, and `_data/MG.dat`
- installs the repository test case from `applications/heat/_ianus/HEAT` as
  `bin/heat/_ianus/HEAT` and includes it in the `heat_stage` layout
- optionally installs `_adc/EWIKON_201912` into `bin/heat/_ianus/HEAT` when it
  exists in the build directory

`cmake --build <build> --target heat_stage -j4` builds heat and refreshes the
launcher, helper, partitioner, mesher, and METIS library beside the executable.
It seeds only missing runtime defaults, preserving parameter edits on repeat
staging. Installation continues to use `cmake --install <build> --prefix <prefix>`.

Regression tests: `python3 -m pytest applications/heat/tests tools/e3d_scripts/tests -q`.

For a two-step, four-rank EWIKON smoke test, set `HEAT_SMOKE_INSTALL` to the
installed `bin/heat` directory and run the same tests. This requires MPI and
uses the installed `_ianus/HEAT` test case directly without copying its geometry;
the test creates an isolated case
with mesh level 1 and retains a `heat.log` under pytest's temporary directory.
Two steps are needed because heat suppresses visualization output on step one.
This checks launcher integration and source-file preservation, not physical
convergence. At mesh level 1, the bundled case can have zero-volume sensors
and NaN PID diagnostics; the previous copying launcher shows the same behavior.

## Common Failure Points

- Missing `-f` input folder, `<input>/heat.s3d`, or
  `<input>/sampleRigidBody.xml`.
- `s3d_mesher -a heat` fails and the case folder has no `meshDir` fallback.
- Geometry paths are missing, exceed native reader limits, or are not valid
  input-relative or absolute paths.
- MPI rank count is too small: the application expects rank 0 plus at least one
  worker partition.
