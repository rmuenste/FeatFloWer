# Repository Guidelines

## Project Structure & Module Organization
Core solvers live in `source/`, with submodules such as `src_pp3d` (Q2/P1 Navier–Stokes), `src_visco` (non-Newtonian rheology), and `src_fbm` (fictitious boundary coupling). Executable front-ends are under `applications/` and follow the `q2p1_*` naming scheme. The `libs/` and `extern/` trees host third-party components (PE, contact dynamics, CGAL helpers). Canonical parameter files reside in `_data/`, while benchmark meshes and setups are collected in `testcases/`. Extended notes and derivations are maintained in `docs/md_docs/`, and build recipes in `build_guide.md`.

## Build, Test, and Development Commands
Initialize dependencies before first build:
```bash
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON
cmake --build build -- -j8
```
Discover tuned toolchains with `cmake -S . -B build -DSHOW_BUILD_IDS=ON`. Run test targets from an existing build directory:
```bash
cmake --build build --target test
ctest --test-dir build --output-on-failure
```
Use `cmake --build build --target install` for staged packaging. For parallel application runs, partition meshes via `./partitioner mesh.tri 4` and launch solvers with `mpirun -np <ranks> ./applications/q2p1_devel/q2p1_devel`.

## Coding Style & Naming Conventions
Fortran sources adopt two-space indentation, uppercase module and subroutine names, and aligned `::` declarations (see `source/Init.f90`). Keep lines under 100 characters and group related declarations. Follow existing directory prefixes (`src_`, `q2p1_`) when introducing new physics or applications. CMake files prefer lowercase commands and cached options prefixed with `USE_` or `BUILD_`. Document new input parameters in `_data` samples or `docs/md_docs/`.

## Testing Guidelines
Primary regression coverage uses CTest; enable via `-DBUILD_TESTING=ON` and invoke `ctest`. Add focused verification by scripting scenario runs in `testcases/`—mirror the directory structure when adding new benchmarks. Tests should assert stability for both sequential and MPI runs; mention required ranks in accompanying documentation.

## Commit & Pull Request Guidelines
Commits follow short, present-tense summaries (e.g., "Update FullC0ntact submodule"), optionally suffixing domain scope (`src_pp3d:`) for clarity. Squash incidental build artifacts before review. PRs must describe motivation, list reproducible build/test commands, and link to issues or benchmark cases. Include before/after metrics or screenshots when altering solver outputs or documentation plots.

## Simulation Workflow Notes
Keep configuration under version control by copying templates from `_data/`. When introducing new numerical options, echo them in `Description.txt` and cross-link relevant theory notes in `docs/md_docs/`. For large-scale jobs, verify partition quality with `tools/` diagnostics before submitting to clusters.

## PE Library Integration
- Enable with `-DUSE_PE=ON`; add `-DUSE_PE_SERIAL_MODE=ON` to bypass PE’s MPI for large particles. CMake adds `-DHAVE_PE` (and `-DPE_SERIAL_MODE`) and links `pe_static` (`CMakeLists.txt` + `cmake/modules/GenerateLinkerFlags.cmake`).
- Fortran apps exclude rank 0, build a subcommunicator, then call `commf2c_*` from `app_init.f90` (e.g., `commf2c_fsi`, `commf2c_drill`, `commf2c_lubrication_lab`; see `applications/*/app_init.f90`).
- Entry points live in `libs/pe/src/interface/c2f_interface.cpp`: parallel mode converts Fortran MPI communicators and dispatches to `setup*` functions in `libs/pe/src/interface/sim_setup.cpp`; serial mode calls header-only `pe::setup*Serial` helpers in `libs/pe/pe/interface/sim_setup_serial.h`.
- Config and runtime parameters flow through `SimulationConfig` (JSON) inside PE; typical files like `example.json` are loaded from the run directory.
- Detailed flow, mode differences, and setup hooks are documented in `docs/md_docs/pe_initialization.md` (also see `libs/pe/CLAUDE.md` for the PE build).
