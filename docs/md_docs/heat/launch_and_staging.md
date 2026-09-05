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
`_data/q2p1_param.dat`, `start/sampleRigidBody.xml`, and `start/data.TXT` are
seeded from the installation. Existing case copies are preserved; parameter
selection is reported. To tune a case, edit its `_data/q2p1_param.dat`.

The launcher performs these steps:

```text
validate paths and rank count; prepare and enter the case directory
copy <input>/heat.s3d unchanged to _data/heat.s3d (skip same-file copy)
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

Geometry is never copied or cleaned up. Relative references in `heat.s3d` are
relative to the case directory, not the input folder. Absolute references allow
shared geometry. For `input/steel.off`, reference `input/steel.off` explicitly.
Existing cases relying on automatic OFF copying must update their references.

## Build And Install Staging

`applications/heat/CMakeLists.txt`:

- builds the `heat` executable from `app_init.f90` and `heat.f90`
- depends on `s3d_mesher`
- copies `heat_start.py` and `sensor_temperature_extraction.sh` into the build
  directory
- installs/stages the shared `e3d_layout.py` helper
- installs `start`, `_data`, `_mesh`, `_vtk`, `_dump`, the partitioner Python
  package, scripts, and `_data/MG.dat`
- optionally installs `_adc/EWIKON_201912` into `bin/heat/_ianus/HEAT` when it
  exists in the build directory

`cmake --build <build> --target heat_stage -j4` builds heat and refreshes the
launcher, helper, partitioner, mesher, and METIS library beside the executable.
It seeds only missing runtime defaults, preserving parameter edits on repeat
staging. Installation continues to use `cmake --install <build> --prefix <prefix>`.

Regression tests: `python3 -m pytest applications/heat/tests tools/e3d_scripts/tests -q`.

For a two-step, four-rank EWIKON smoke test, set `HEAT_SMOKE_INSTALL` to the
installed `bin/heat` directory and run the same tests. This requires MPI and
the optional installed `EWIKON_201912` example; the test creates an isolated
case with mesh level 1 and retains a `heat.log` under pytest's temporary directory.
Two steps are needed because heat suppresses visualization output on step one.

## Common Failure Points

- Missing `-f` case folder or missing `<case>/heat.s3d`.
- `s3d_mesher -a heat` fails and the case folder has no `meshDir` fallback.
- Geometry paths inside `heat.s3d` are not valid case-relative or absolute paths.
- MPI rank count is too small: the application expects rank 0 plus at least one
  worker partition.
