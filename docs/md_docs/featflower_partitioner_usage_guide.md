# FeatFloWer Partitioner: Installation and Usage Guide

This guide covers installation and practical use of the `featflower-partitioner` CLI
tool — an installable Python package wrapping the METIS 5 mesh partitioner used by the
FeatFloWer CFD framework.

- Package location: `tools/featflower_partitioner/`
- Source it wraps: `tools/partitioner/` (left untouched, still used at runtime via PYTHONPATH)
- Native dependency: `libmetis.so` built as part of the FeatFloWer CMake build

## 1) Prerequisites

- Python 3.8+
- `pip`
- CMake ≥ 3.18 and a C/C++ compiler (needed to build METIS)
- Git 2.17+ with submodules initialised

## 2) Installation

### Step 1 — Clone and initialise submodules

```bash
git clone <repo-url> FeatFloWer
cd FeatFloWer
git submodule update --init --recursive
```

The METIS 5.1.0 source lives in `extern/libraries/metis-5.1.0/` and is pulled in by
this step.

### Step 2 — Create a virtual environment

It is strongly recommended to install into a virtual environment rather than the
system or user Python, to keep dependencies isolated and avoid conflicts.

```bash
python3 -m venv ~/.venvs/featflower
source ~/.venvs/featflower/bin/activate
```

> **Tip:** Add the `source` line to your shell's startup file (`.bashrc`, `.zshrc`,
> etc.) so the environment is activated automatically in every new session.
>
> You can place the venv anywhere — `~/.venvs/featflower` is just a convention.
> Some teams keep it inside the repo at `.venv/` (add `.venv/` to `.gitignore` in
> that case).

Once activated, `python` and `pip` refer to the venv and all packages installed
there are fully isolated from the system Python.

### Step 4 — Build METIS

The partitioner requires `libmetis.so`. The `BUILD_METIS` CMake option is `ON` by
default, so a minimal CMake build is enough:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=OFF ..
make -j$(nproc) metis
cd ..
```

> `-DBUILD_APPLICATIONS=OFF` skips the full Fortran solver build. If you are building
> FeatFloWer anyway, a standard `make -j$(nproc)` will produce `libmetis.so` as a
> by-product — no extra step needed.

After this step the shared library is at:

```
build/extern/libraries/metis-5.1.0/libmetis/libmetis.so
```

### Step 5 — Install the Python package

Set `FEATFLOWER_BUILD_DIR` to your build directory and install:

```bash
export FEATFLOWER_BUILD_DIR=$(pwd)/build
pip install -e tools/featflower_partitioner/
```

During installation the `setup.py` hook copies `libmetis.so` from the build tree
directly into the package directory. Once installed, **no environment variable is
needed at runtime** — the library travels with the Python environment.

> The `-e` flag installs in editable mode: the package reflects any future changes to
> the source without reinstalling.

### Step 6 — Verify

```bash
featflower-partition --help
```

Expected output:

```
usage: featflower-partition [-h] NPart Method NSubPart MeshName ProjectFile

Partition a FeatFloWer mesh for parallel execution (METIS 5).

positional arguments:
  NPart        Number of partitions (>=1)
  Method       Partitioning method (1=recursive, 2/3=kway, negative=axis-based)
  NSubPart     Number of subgrids (>=1)
  MeshName     Output mesh name (used as directory key under _mesh/)
  ProjectFile  Path to .prj project file
```

If `~/.local/bin` is not in your `$PATH`, invoke via:

```bash
python -m featflower_partitioner --help
```

### Summary of install commands

```bash
git clone <repo-url> FeatFloWer && cd FeatFloWer
git submodule update --init --recursive

python3 -m venv ~/.venvs/featflower
source ~/.venvs/featflower/bin/activate

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=OFF ..
make -j$(nproc) metis
cd ..

export FEATFLOWER_BUILD_DIR=$(pwd)/build
pip install -e tools/featflower_partitioner/
```

## 3) Command-Line Usage

```
featflower-partition <NPart> <Method> <NSubPart> <MeshName> <ProjectFile>
```

| Argument | Type | Description |
|---|---|---|
| `NPart` | integer ≥ 1 | Total number of partitions (MPI ranks) |
| `Method` | integer | Partitioning method (see table below) |
| `NSubPart` | integer ≥ 1 | Number of subgrids the mesh is first split into |
| `MeshName` | string | Output directory key — result written to `_mesh/<MeshName>/` |
| `ProjectFile` | path | Path to the `.prj` project file for the mesh |

### Partitioning methods

| Method | Algorithm |
|--------|-----------|
| `1` | METIS recursive bisection (`METIS_PartGraphRecursive`) |
| `2` | METIS k-way partitioning (`METIS_PartGraphKway`) |
| `3` | METIS k-way partitioning (alias for `2`) |
| `11`, `12`, `13` | As above, but subgrids are written in reversed order |
| negative (e.g. `-1`, `-12`, `-123`) | Axis-based splitting along coordinate directions 1, 2, or 3 |

### Example: 4-way partition of the 2D FAC benchmark

```bash
featflower-partition 4 1 1 NEWFAC testmeshes/2D_FAC/2Dbench.prj
```

Expected output:

```
The projekt folder consists of the following files:
- Grid File: testmeshes/2D_FAC/ffff.tri
- Boundary Files:
  * testmeshes/2D_FAC/bottom.par
  ...
Grid input file: '_mesh/NEWFAC/sub0001/GRID.tri'
...
Calling Metis...
26 edges were cut by Metis.
32
33
32
33
The partitioning was successful!
```

The four numbers (32, 33, 32, 33) are the element counts in each partition. The
partitioned files are written to:

```
_mesh/NEWFAC/
  GRID.tri          # original grid (copy)
  GRID.prj          # original project file (copy)
  sub0001/
    GRID.tri        # full mesh for subgrid 1
    GRID.tri        # ...
    <boundary>.par  # boundary files
  sub0001/
    GRID0001.tri    # partition 1 of subgrid 1
    GRID0002.tri    # partition 2 of subgrid 1
    ...
```

### Example: 2-subgrid, 8-total-partition run

```bash
featflower-partition 8 1 2 MESH3D path/to/mesh.prj
```

This first splits the mesh into 2 subgrids (using METIS recursive bisection), then
partitions each subgrid into 4 parts, giving 8 total partitions.

## 4) How libmetis.so Is Located

The package searches for `libmetis.so` in the following order at runtime:

| Priority | Location | Notes |
|----------|----------|-------|
| 1 | `<package_dir>/libmetis.so` | Bundled during `pip install` — used when `FEATFLOWER_BUILD_DIR` was set at install time |
| 2 | `$FEATFLOWER_BUILD_DIR/extern/libraries/metis-5.1.0/libmetis/libmetis.so` | Build tree fallback if env var is set |
| 3 | `./libmetis.so` | Legacy: current working directory |
| 4 | `../lib64/libmetis.so` | Legacy: relative path |
| 5 | `libmetis.so` | System library path via `ldconfig` |

If no location succeeds, a `RuntimeError` is raised listing every path that was tried,
along with a hint to set `FEATFLOWER_BUILD_DIR`.

## 5) Reinstalling After a Rebuild

If you rebuild METIS (e.g. after a clean build), activate your venv and reinstall the
package to pick up the new library:

```bash
source ~/.venvs/featflower/bin/activate
export FEATFLOWER_BUILD_DIR=$(pwd)/build
pip install -e tools/featflower_partitioner/
```

This re-runs the `setup.py` hook and overwrites the bundled `libmetis.so` with the
freshly built copy.

## 6) Python API

The package also exposes a Python API for use in scripts:

```python
from featflower_partitioner import partition

partition(
    NPart=4,
    PartMethod=1,
    NSubPart=1,
    MeshName="NEWFAC",
    ProjektFile="testmeshes/2D_FAC/2Dbench.prj",
)
```

The `partition()` function creates the `_mesh/` output directory if it does not exist
and then runs the full partitioning pipeline.

Full public API:

| Symbol | Description |
|--------|-------------|
| `partition(NPart, PartMethod, NSubPart, MeshName, ProjektFile)` | Top-level entry point — creates `_mesh/` and calls `MainProcess` |
| `MainProcess(nnPart, pMethod, nSubMesh, MeshName, ProjektFile)` | Core orchestration routine |
| `checkParameters(params)` | Validate a `sys.argv`-style parameter list; raises `sys.exit` on error |
| `mkdir(dir)` | Create a directory only if it does not already exist |

## 7) Relationship to Legacy Components

The package wraps the code in `tools/partitioner/` without modifying it. The original
files remain in place and continue to be used by `e3d_start.py` via `PYTHONPATH`.

| Legacy component | Role of `featflower-partitioner` |
|-----------------|----------------------------------|
| `tools/partitioner/part.py` | Adapted copy lives in the package; METIS loading replaced by the prioritised `_libmetis.py` loader |
| `tools/partitioner/part_main.py` | Verbatim copy in the package |
| `tools/PyPartitioner.py` | Thin `sys.argv` wrapper — still works unchanged |
| `tools/partpy/` | METIS 4 variant with `.pls` and Cartesian `-4` features — separate future migration, not included here |
