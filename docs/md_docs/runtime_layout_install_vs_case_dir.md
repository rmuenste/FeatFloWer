# gendie: separating the installation from the case folder

## 1. Request

A gendie user asked for binaries in one place and per-simulation folders that
contain only configuration. In the follow-up the requirement was sharpened:

* `e3d_start_yaml.py` is the launcher of interest. It lives in the install
  (`bin/q2p1_gendie`), not in the case folder.
* The launcher is told where the case folder is (OpenFOAM style). Where it is
  invoked from must not matter. The e3d input folder is always passed
  explicitly.
* The case folder holds only what differs, or could differ, between cases.
  Binaries, launcher scripts, `SCALEXA`, `_ianus` and the like are not copied
  per case.
* "Once, properly" is preferred over a quick fix.

## 2. What the install currently contains, classified

Everything below is installed into `bin/q2p1_gendie` by
`applications/q2p1_sse/CMakeLists.txt` and the sibling application CMake files.

**Static, belongs to the installation only**

| Item | Note |
|---|---|
| `q2p1_sse`, `q2p1_sse_temp`, `s3d_mesher`, `q1_scalar_multimat`, `autoparam`, `stream`, `STLvsTRI`, `meshref` | executables, rpath `$ORIGIN/../lib` |
| `partitioner/` | Python package, imported by the launchers |
| `e3d_start.py`, `e3d_start_yaml.py`, `conv_check.sh`, `computeAreas.py`, `RunnerGenDIE.sh`, `RunnerScript.sh`, `RunSCALEXA.sh`, `RankFileGenerator.sh`, `slurm_Veka_*.sh` | scripts |
| `e3d_die*.yaml` | default execution plans |
| `_data_BU/*.dat` | default solver parameter templates |
| `_data/MG.dat` | static lookup table read by `CreateDumpStructures` (`source/OutputProfiles.f90:1803`), needed by the umbrella/dump code of `q2p1_sse` |
| `start/sampleRigidBody.xml`, `start/data.TXT` | static rigid-body world and FullC0ntact parameters read by `init_fc_rigid_body` (`FullC0ntact/inshape3dcore/fortrancppinterface/init_func.hpp`) from the working directory, for every SSE run |
| `PATCHES/`, `param_meshref.cfg`, `param_vtu_mesher.cfg` | inputs of the preprocessing tools, read relative to their working directory |
| `SCALEXA/` (81 MB), `_ianus/` (12 MB), `VEKA_S/` | tooling and example cases, never read by the gendie solvers |
| `_mesh _vtk _dump _1D _hist _prot0 _prot1 _RTD` | empty output skeleton |

**Per case, belongs to the case folder**

| Item | Note |
|---|---|
| e3d input folder (`setup.e3d` or `Extrud3D.dat`, `meshDir/`, `*.off`) | passed with `-f`, may live anywhere |
| execution plan yaml | passed with `-y`; the shipped default or a case-specific copy |
| optional overrides of `_data_BU/*.dat` | only if the case needs tuned templates |
| generated at run time: `_data/{q2p1_param.dat, Extrud3D_0.dat, Extrud3D.dat, meshDir, prot*}`, `_mesh/`, `_vtk/`, `_dump/`, `_1D/`, `_hist/`, `_prot*/`, `_RTD/`, `e3d.log`, `mesh_names.offs`, copied `*.off` | all relative to the working directory |

## 3. Solver contract (unchanged)

The Fortran executables address every file relative to the **working
directory** with fixed prefixes (`_data/`, `_mesh/`, `_vtk/`, ...). They never
look up their own location or a sibling executable, and `mpirun` and `srun`
propagate the working directory to all ranks. The contract is therefore:

> working directory of the solver == case folder

This is a stable convention across all FeatFloWer applications and is the
anchor of the design below. Nothing in Fortran changes.

## 4. Design

### 4.1 Command-line contract of `e3d_start_yaml.py`

```
e3d_start_yaml.py -c <case_dir> -f <e3d_input_dir> -y <plan.yaml> -n <ranks> [...]
```

* `-c/--case DIR` (new): the case folder. Default: the current working
  directory, which keeps every existing invocation valid.
* `-f` and `-y` are resolved relative to the **invocation** directory (normal
  CLI semantics) and made absolute before anything else happens. A bare plan
  name (`-y e3d_die.yaml`) is looked up in the case folder and then in the
  install.
* Everything else is unchanged.

Typical use:

```bash
# case folder contains only e3d_input/
cd /scratch/cases/veka_01
python3 /opt/ff/bin/q2p1_gendie/e3d_start_yaml.py -f e3d_input -y e3d_die.yaml -n 64

# or from anywhere
python3 /opt/ff/bin/q2p1_gendie/e3d_start_yaml.py -c /scratch/cases/veka_01 -f /scratch/cases/veka_01/e3d_input -y e3d_die.yaml -n 64
```

SLURM: `#SBATCH --chdir=<case>` or `cd "$SLURM_SUBMIT_DIR"`, then the line above.

### 4.2 One shared layout module

`tools/e3d_scripts/e3d_layout.py`, installed next to the launchers. It owns all
knowledge about the two locations, so `e3d_start.py` can adopt it later and the
two launchers cannot drift apart on this point.

```python
class RunLayout:
    install_dir: Path      # Path(__file__).resolve().parent, override FF_GENDIE_HOME
    case_dir: Path         # from -c, default Path.cwd()

    def exe(self, name) -> str
        # install_dir / name, absolute

    def resolve_input(self, p) -> Path
        # absolute -> as is
        # case_dir / p if it exists
        # install_dir / p if it exists
        # else FileNotFoundError naming both candidates
        # prints "<p>: using <case|install> copy" once per file

    def prepare_case(self)
        # mkdir -p the output skeleton
        # seed _data/MG.dat from the install if missing
        # nothing else is copied

    def enter_case(self)
        # os.chdir(case_dir) once, after -f/-y are absolute
```

### 4.3 Changes inside `e3d_start_yaml.py`

1. Parse `-c`, build `RunLayout`, make `-f`/`-y` absolute, call
   `prepare_case()` and `enter_case()` at the top of `main()`.
2. Executables: the five launch sites (`s3d_mesher`, `q2p1_sse`,
   `q2p1_sse_temp`, `q1_scalar_multimat`, and the `srun` variants) use
   `layout.exe(name)` instead of `"./name"` or `os.getcwd() + "/name"`.
3. Templates: the `_data_BU/...` lookups in `folderSetup`,
   `simLoopTemperatureXSE`, the DIE final step and
   `_data_BU/DIE/mesh_names.offs` go through `layout.resolve_input`.
4. Plans: `validate_yaml_stage_inputs` and `run_*_step` resolve `param_file`
   through `layout.resolve_input`. The shipped plans keep their
   `_data_BU/...` entries unchanged; a case-local override with the same
   relative name wins.
5. `import partitioner` needs no change: Python puts the script's own directory
   on `sys.path`.

Roughly 25 changed lines in the launcher plus the ~120 line module. Every path
literal that remains in the script is relative to the working directory and is
now, by construction, relative to the case folder.

### 4.4 Why a single `chdir` rather than threading `case_dir` through every path

The script has about 60 working-directory-relative path literals spread over
duplicated loops, and the Fortran side requires the working directory to be the
case folder anyway. Changing the working directory once, explicitly, at a
single well-named point makes the script honour exactly the contract the
solvers already impose. Threading a `case_dir` argument through every function
would touch most of the file, change nothing for the solvers, and leave the
`mpirun` launch needing the same `chdir` regardless.

### 4.5 Preprocessing (`STLvsTRI`, `meshref`)

These read `param_meshref.cfg` and `PATCHES/` from their working directory and
take the case input with `-f`. With the recently added `bin/preprocessing`
install they already work as a fixed tool folder: run them with the working
directory in `bin/preprocessing` and an absolute `-f <case>/e3d_input`. No
change needed now. A later improvement is a `-p <template dir>` option for
`meshref` so it can be run from the case folder as well.

### 4.6 CMake

No change is required. Optional tidy-up once nothing copies the install
anymore: stop installing `_ianus`, `VEKA_S`, `start` and the empty output
directories into `bin/q2p1_gendie`, or move the examples to `share/gendie`.
`_data/MG.dat` and `_data_BU/` stay in the install as the seed and default
templates.

### 4.7 Implementation status (branch `feature/gendie-case-dir`)

Implemented as described above, with two deviations found during testing:

* The short option for the case folder is `-C` (long `--case`), because `-c`
  was already taken by `--host-conf`.
* `start/sampleRigidBody.xml` and `start/data.TXT` are read by the C++
  rigid-body initialisation of every SSE run (not visible to a Fortran-only
  grep). They are seeded into the case folder like `MG.dat`.
* The `partitioner` package loaded `libmetis.so` from the working directory
  (then `../lib64`, then the system path). With the case folder as working
  directory this failed. `tools/partitioner/part.py` now also searches next to
  the package itself (`bin/q2p1_gendie/libmetis.so`, where the install puts a
  copy) and the install prefix's `lib64`/`lib`, before falling back to the
  system path. The working-directory lookup stays first, so the classic
  layout is unaffected.

Files touched: `tools/e3d_scripts/e3d_layout.py` (new),
`tools/e3d_scripts/e3d_start_yaml.py`, `tools/partitioner/part.py`,
`applications/q2p1_sse/CMakeLists.txt` (install and stage the module),
`tools/e3d_scripts/tests/test_e3d_layout.py` (new),
`docs/md_docs/q2p1_sse/{launch_and_staging,yaml_simulation_plans}.md`.

Tests: `python3 -m pytest tools/e3d_scripts/tests -q` runs unit tests of the
layout module plus two end-to-end runs of the real launcher against a fake
installation whose executables are stub scripts (case outside the install,
and the classic invocation inside the install). They need a built
`libmetis.so` in any `build*` directory for the partitioner.

## 5. Compatibility and risks

* No `-c` given and launched from inside the install folder: install dir and
  case dir coincide, behaviour is identical to today.
* Fortran, parameter formats, mesh formats and the install tree are unchanged.
* New failure mode: a stale case-local `_data_BU` override silently shadows an
  updated install template. Mitigated by the one-line "using case/install
  copy" message for every resolved input.
* Cluster runs need the install prefix visible on all nodes. Already true,
  since users copy the binaries out of that prefix today.
* `e3d_start.py` keeps working unchanged; it can adopt the module in a second
  step.

## 6. Validation

1. Install to a prefix. Create an empty folder with only `e3d_input/` (copy of
   `_ianus/DIE/RH_LOC`). Run the launcher with `-n 4 --short-test` from that
   folder and from an unrelated folder with `-c`. Expect all outputs in the
   case folder and nothing written into the prefix.
2. Run inside `bin/q2p1_gendie` as today and compare `_data/q2p1_param.dat`
   and `_1D` outputs against a pre-change run.
3. Add a case-local `_data_BU/q2p1_paramV_DIE_1.dat`, rerun, confirm the
   override is reported and used.
