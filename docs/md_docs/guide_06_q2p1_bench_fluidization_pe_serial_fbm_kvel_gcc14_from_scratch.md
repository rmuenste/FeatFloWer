# FeatFloWer Guide 06: Configure and Build `q2p1_bench_fluidization` with GCC 14 + PE Serial + JSON + Eigen + FBM/KVEL Acceleration

This guide is written as a direct execution checklist for a smaller coding model.
Follow the commands in order. Do not switch targets to `q2p1_ATC`; this guide is for:

- Application directory: `applications/q2p1_bench_fluidization`
- CMake target: `q2p1_bench_fluidization`
- Executable: `applications/q2p1_bench_fluidization/q2p1_bench_fluidization`
- PE runtime config: `applications/q2p1_bench_fluidization/example.json`
- Solver runtime config: `applications/q2p1_bench_fluidization/_data/q2p1_param.dat`

The validated configuration is:

- `-DUSE_PE=ON`
- `-DUSE_PE_SERIAL_MODE=ON`
- `-DUSE_JSON=ON`
- `-DEIGEN=ON`
- `-DENABLE_FBM_ACCELERATION=ON` (HashGrid + KVEL force acceleration)
- `-DUSE_CGAL=OFF` (default; CGAL is not needed for `q2p1_bench_fluidization`)
- GCC/G++/GFortran 14.3.0 through the OpenMPI 4.1.6 GCC 14 wrapper compilers

## 1) Environment Setup

From the repository root, load the GCC 14 compiler suite and matching OpenMPI wrapper stack:

```bash
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v14 openmpi/options/interface/ethernet openmpi/4.1.6
module list
```

Check that the wrapper compilers resolve to GCC 14:

```bash
which gcc
gcc --version | head -1
which mpicc
mpicc --version | head -1
which mpicxx
mpicxx --version | head -1
which mpifort
mpifort --version | head -1
```

Expected compiler version indicator:

```text
gcc (GCC) 14.3.0
g++ (GCC) 14.3.0
GNU Fortran (GCC) 14.3.0
```

On some cluster or sandbox shells, the OpenMPI wrappers may print:

```text
opal_ifinit: socket() failed with errno=1
```

This warning did not stop configure or build in the validated run.

## 2) Initialize Submodules

From the repository root:

```bash
git submodule update --init --recursive
```

Important for this build:

- PE is built from `libs/pe`.
- FullC0ntact is configured as a subproject.
- Eigen is required because this guide uses `-DEIGEN=ON`.

Quick check:

```bash
ls FullC0ntact/libs/eigen
```

If Eigen is not already available locally, CMake may download it during configure.
If the download fails with a DNS or network error, rerun configure in an environment with network access.

## 3) Configure a New Build Directory

Use a fresh build directory for this workflow:

```bash
cmake -S . -B build-q2p1-bench-fluidization-gcc14-pe-serial-fbm -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DUSE_PE=ON \
  -DUSE_PE_SERIAL_MODE=ON \
  -DUSE_JSON=ON \
  -DEIGEN=ON \
  -DENABLE_FBM_ACCELERATION=ON \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort
```

Do not add `-DUSE_CGAL=ON` for this application. In the validated run, configure reported:

```text
CGAL support is disabled in pe.
```

Expected configure indicators:

- `C compiler identification is GNU 14.3.0`
- `CXX compiler identification is GNU 14.3.0`
- `Fortran compiler identification is GNU 14.3.0`
- `PE SERIAL MODE ENABLED`
- `Accelerated point query (HashGrid spatial hashing) is ENABLED`
- `Eigen support enabled`
- `CGAL support is disabled in pe.`
- `FBM acceleration enabled (HashGrid + KVEL). Control at runtime via q2p1_param.dat`
- `Configuring done`
- `Generating done`
- `Build files have been written to: .../build-q2p1-bench-fluidization-gcc14-pe-serial-fbm`

## 4) Build Only `q2p1_bench_fluidization`

Build the application target:

```bash
cmake --build build-q2p1-bench-fluidization-gcc14-pe-serial-fbm \
  --target q2p1_bench_fluidization -- -j8
```

This is a cold build, so Ninja may compile thousands of shared library objects first.
That is expected. The validated run built `4752` steps and finished by linking:

```text
Linking Fortran executable applications/q2p1_bench_fluidization/q2p1_bench_fluidization
```

## 5) Verify the Build Output

Check the executable:

```bash
ls -lh build-q2p1-bench-fluidization-gcc14-pe-serial-fbm/applications/q2p1_bench_fluidization/q2p1_bench_fluidization
```

Validated result:

```text
-rwxr-xr-x ... 47M ... build-q2p1-bench-fluidization-gcc14-pe-serial-fbm/applications/q2p1_bench_fluidization/q2p1_bench_fluidization
```

Verify the important CMake cache options:

```bash
rg -n "^(CMAKE_(C|CXX|Fortran)_COMPILER|EIGEN|USE_PE|USE_PE_SERIAL_MODE|USE_JSON|USE_CGAL|ENABLE_FBM_ACCELERATION):" \
  build-q2p1-bench-fluidization-gcc14-pe-serial-fbm/CMakeCache.txt
```

Expected values:

- `CMAKE_C_COMPILER:UNINITIALIZED=mpicc`
- `CMAKE_CXX_COMPILER:UNINITIALIZED=mpicxx`
- `CMAKE_Fortran_COMPILER:UNINITIALIZED=mpifort`
- `EIGEN:BOOL=ON`
- `ENABLE_FBM_ACCELERATION:BOOL=ON`
- `USE_CGAL:BOOL=OFF`
- `USE_JSON:BOOL=ON`
- `USE_PE:BOOL=ON`
- `USE_PE_SERIAL_MODE:BOOL=ON`

## 6) Optional Runtime Staging Target

The application defines a staging target:

```bash
cmake --build build-q2p1-bench-fluidization-gcc14-pe-serial-fbm \
  --target q2p1_bench_fluidization_stage -- -j8
```

This target depends on `q2p1_bench_fluidization` and `metis`, copies runtime files, and partitions the mesh with:

- `applications/q2p1_bench_fluidization/example.json` to `example.json`
- `applications/q2p1_bench_fluidization/_data/q2p1_param.dat` to `_data/q2p1_param.dat`
- repository `_data/MG.dat` to `_data/MG.dat`
- `libmetis.so` to the application build directory

```text
tools/PyPartitioner.py ${Q2P1_BENCH_FLUIDIZATION_PARTS} 1 1 NEWFAC _adc/benchSym/bench.prj
```

The default partition count is controlled by:

```text
Q2P1_BENCH_FLUIDIZATION_PARTS=31
```

This guide validates compile/link only. Run the staging target when you also need runtime files and mesh partitions.

## 7) Known Warnings Seen in the Validated Build

The successful GCC 14 build emitted warnings. They did not stop the final executable from linking.

OpenMPI wrapper warning:

```text
opal_ifinit: socket() failed with errno=1
```

Fortran warnings:

- C++-only flags passed to Fortran commands, such as `-std=c++17`.
- COMMON block size warnings for `rparm` and `iparm`.
- Unused variable and unused dummy argument warnings.

PE/C++ warnings:

- Deprecated `std::auto_ptr`.
- `requires` is a C++20 keyword compatibility warning.
- Hidden virtual overload warnings around `containsPoint`.
- Unused variable warnings in VTK/utility code.

Linker warnings:

```text
size of symbol `triad_' changed
size of symbol `triaa_' changed
size of symbol `mgtrd_' changed
size of symbol `mgtra_' changed
```

In the validated build, these were warnings only.

## 8) Minimal Success Criteria

The workflow is successful when all of the following are true:

- Configure exits with code `0`.
- Build exits with code `0`.
- The executable exists at:

```text
build-q2p1-bench-fluidization-gcc14-pe-serial-fbm/applications/q2p1_bench_fluidization/q2p1_bench_fluidization
```

- `CMakeCache.txt` contains `USE_CGAL:BOOL=OFF`, `USE_PE_SERIAL_MODE:BOOL=ON`, and `ENABLE_FBM_ACCELERATION:BOOL=ON`.
