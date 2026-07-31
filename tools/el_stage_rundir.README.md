# el_stage_rundir.sh — independent run directories for EL validation cases

## Why

All Euler-Lagrange validation runs used to execute inside one shared staged
directory (`build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow/`). Two
concurrent runs (or agents) editing `_data/q2p1_param.dat` / `example.json`
there clobber each other. This script stages any validation case into an
arbitrary, self-contained run directory so runs can proceed in parallel.

## Usage

```sh
tools/el_stage_rundir.sh [-f] <case-name> <target-dir>
```

- `<case-name>`: a subdirectory of
  `applications/q2p1_el_pipeflow/validation_cases/`
  (currently `pipe_hp_check`, `v3_ss_frozen`, `v1b_tencate_settling`).
  Inputs are always taken from the repository source tree, never from a
  shared staged build directory.
- `<target-dir>`: created if missing. An existing non-empty directory is
  refused unless `-f` is given.
- `-f`: stage into a non-empty directory anyway (existing mesh folders of
  the same name are replaced; other files are left alone).

What gets staged into `<target-dir>`:

- `_data/q2p1_param.dat` — from the case's `q2p1_param.dat`
- `_data/MG.dat` — from the repo-level `_data/MG.dat`
- `example.json`, `particles.xyz` — from the case directory
- each coarse project folder under `<case>/mesh/*` except `partitions`
  (e.g. `ten_cate_box/`, `pipe_ogrid_z/`) — symlinked meshes are
  dereferenced (`v3_ss_frozen/mesh` is a symlink to `../pipe_hp_check/mesh`)
- `<case>/mesh/partitions/*` into `_mesh/` (e.g. `_mesh/TENCATE27/sub0001…`)
- empty output dirs: `_vtk _gmv _ns _dump _sol _sol/1 solution testresults
  protocols`
- an `_adc` symlink replicated from the reference staged dir, if present
  there

The binary is **not** copied; the script prints the exact `mpirun` line at
the end, e.g.:

```sh
BIN=/home/user/rmuenste/nobackup/code/FF-EL/FeatFloWer/build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
mpirun --oversubscribe --wdir <target-dir> -np 28 "$BIN"
```

## Environment overrides

- `EL_BIN` — binary path used in the printed mpirun hint
- `EL_REF_DIR` — reference staged dir whose `_adc` symlink is mirrored
  (default: the shared gcc14 staged dir)

## Typical smoke-run environment

```sh
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
export OMPI_MCA_rmaps_base_oversubscribe=1
MPIPREFIX=$(ldd "$BIN" | grep -m1 'libmpi\.so' | awk '{print $3}' | sed 's#/lib/libmpi.*##')
export PATH="$MPIPREFIX/bin:$PATH" LD_LIBRARY_PATH="$MPIPREFIX/lib:${LD_LIBRARY_PATH:-}"
```

Notes:

- To shorten a run, edit `SimPar@MaxNumStep` / `SimPar@MaxSimTime` in the
  NEW dir's `_data/q2p1_param.dat` and `timesteps_` in the NEW dir's
  `example.json`. Keep PE `stepsize_` equal to `SimPar@TimeStep` — a
  startup assertion checks `|dt - stepsize_|`.
- Never edit or run inside the shared staged directory again; stage a fresh
  run dir instead.
