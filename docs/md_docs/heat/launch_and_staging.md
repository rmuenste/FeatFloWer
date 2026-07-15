# heat Launch And Staging

## Entry Points

- `applications/heat/heat_start.py`
- `applications/heat/CMakeLists.txt`
- `applications/heat/_data/heat.s3d`
- `applications/heat/sensor_temperature_extraction.sh`

## Launcher Contract

`heat_start.py` expects a case folder through `-f/--in-folder` and a total MPI
rank count through `-n/--num-processors`.

The launcher performs these steps:

```text
parse command line
copy *.off/*.OFF from case folder into working directory
copy <case>/heat.s3d to _data/heat.s3d
remove old _data/meshDir
run ./s3d_mesher -a heat
if no generated meshDir exists, copy <case>/meshDir
partition _data/meshDir/file.prj for numProcessors - 1 worker ranks
launch ./heat with mpirun or srun
remove copied OFF files from working directory
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

- Default: `mpirun -np <numProcessors> ./heat`
- `-u/--use-srun`: `srun ./heat`

The launcher propagates the solver exit code after cleanup.

## Build And Install Staging

`applications/heat/CMakeLists.txt`:

- builds the `heat` executable from `app_init.f90` and `heat.f90`
- depends on `s3d_mesher`
- copies `heat_start.py` and `sensor_temperature_extraction.sh` into the build
  directory
- installs `start`, `_data`, `_mesh`, `_vtk`, `_dump`, the partitioner Python
  package, scripts, and `_data/MG.dat`
- optionally installs `_adc/EWIKON_201912` into `bin/heat/_ianus/HEAT` when it
  exists in the build directory

## Common Failure Points

- Missing `-f` case folder or missing `<case>/heat.s3d`.
- `s3d_mesher -a heat` fails and the case folder has no `meshDir` fallback.
- OFF/STL paths inside `heat.s3d` do not match the copied case assets or the
  installed `_ianus/HEAT` layout.
- MPI rank count is too small: the application expects rank 0 plus at least one
  worker partition.
