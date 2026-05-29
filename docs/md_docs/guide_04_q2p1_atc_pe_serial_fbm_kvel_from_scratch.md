# FeatFloWer Guide 04: Configure and Build `q2p1_ATC` with PE Serial + CGAL + JSON + Eigen + FBM/KVEL Acceleration

This guide is the next installment in the from-scratch series.
It documents the exact configuration that was validated for `q2p1_ATC`:

- `-DUSE_PE=ON`
- `-DUSE_PE_SERIAL_MODE=ON`
- `-DUSE_JSON=ON`
- `-DUSE_CGAL=ON`
- `-DEIGEN=ON`
- `-DENABLE_FBM_ACCELERATION=ON` (HashGrid + KVEL force acceleration)

## 1) Environment Setup (GCC 13 validated stack)

On RHEL 9.7, load the same compiler/MPI modules used in Guide 02:

```bash
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
module list
```

## 1.1) Optional Environment Setup with GCC 14

`q2p1_ATC` can also be configured with the GCC 14 compiler suite.
Keep the GCC 13 workflow above for reproducing older validated builds; use this stack when you explicitly want GCC 14:

```bash
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v14 openmpi/options/interface/ethernet openmpi/4.1.6
module list
```

Quick compiler check:

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

Expected GCC 14 version indicator:

```text
gcc (GCC) 14.3.0
g++ (GCC) 14.3.0
GNU Fortran (GCC) 14.3.0
```

## 2) Initialize Submodules (required)

From repository root:

```bash
git submodule update --init --recursive
```

Important for this build:
- `FullC0ntact/libs/eigen/` must contain Eigen headers (`Eigen/`, `unsupported/`).
- If this folder is empty, `q2p1_ATC` build fails with `Eigen/...` include errors.

Quick check:

```bash
ls FullC0ntact/libs/eigen
```

## 3) Configure a Dedicated Build Directory

Use the exact configure command that worked:

```bash
cmake -S . -B build-atc-ninja-release-eigen -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DUSE_PE=ON \
  -DUSE_PE_SERIAL_MODE=ON \
  -DUSE_JSON=ON \
  -DUSE_CGAL=ON \
  -DEIGEN=ON \
  -DENABLE_FBM_ACCELERATION=ON \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort
```

For a GCC 14 build, use the GCC 14 module stack from section 1.1 and a separate build directory:

```bash
cmake -S . -B build-atc-ninja-release-eigen-gcc14 -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DUSE_PE=ON \
  -DUSE_PE_SERIAL_MODE=ON \
  -DUSE_JSON=ON \
  -DUSE_CGAL=ON \
  -DEIGEN=ON \
  -DENABLE_FBM_ACCELERATION=ON \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort
```

Expected configure indicators:
- `PE SERIAL MODE ENABLED`
- JSON dependency fetch/availability
- FBM acceleration enabled (`HashGrid + KVEL`)

## 4) Build `q2p1_ATC`

Recommended first-time sequence:

```bash
cmake --build build-atc-ninja-release-eigen --target cgal -- -j8
cmake --build build-atc-ninja-release-eigen --target q2p1_ATC -- -j8
```

For the GCC 14 build directory:

```bash
cmake --build build-atc-ninja-release-eigen-gcc14 --target cgal -- -j8
cmake --build build-atc-ninja-release-eigen-gcc14 --target q2p1_ATC -- -j8
```

The first command avoids first-build ordering issues where PE sources may compile before CGAL headers are staged.

To keep a full log:

```bash
cmake --build build-atc-ninja-release-eigen --target q2p1_ATC -- -j8 \
  > build-atc-ninja-release-eigen/build_q2p1_ATC.log 2>&1
```

## 5) Verify the Result

Executable location:

```bash
ls -lh build-atc-ninja-release-eigen/applications/q2p1_ATC/q2p1_ATC
```

For the GCC 14 build:

```bash
ls -lh build-atc-ninja-release-eigen-gcc14/applications/q2p1_ATC/q2p1_ATC
```

Verify CMake cache options:

```bash
rg -n "^(EIGEN|USE_PE|USE_PE_SERIAL_MODE|USE_JSON|USE_CGAL|ENABLE_FBM_ACCELERATION):" \
  build-atc-ninja-release-eigen/CMakeCache.txt
```

For the GCC 14 build:

```bash
rg -n "^(EIGEN|USE_PE|USE_PE_SERIAL_MODE|USE_JSON|USE_CGAL|ENABLE_FBM_ACCELERATION):" \
  build-atc-ninja-release-eigen-gcc14/CMakeCache.txt
```

Expected values:
- `EIGEN:BOOL=ON`
- `USE_PE:BOOL=ON`
- `USE_PE_SERIAL_MODE:BOOL=ON`
- `USE_JSON:BOOL=ON`
- `USE_CGAL:BOOL=ON`
- `ENABLE_FBM_ACCELERATION:BOOL=ON`

## 6) Notes About FBM/KVEL Acceleration

`-DENABLE_FBM_ACCELERATION=ON` enables compile-time acceleration code paths (HashGrid alpha + KVEL force).
Runtime behavior is controlled via the application parameter file (`q2p1_param.dat`), as indicated by configure output.

## 7) Known Warnings Seen in This Configuration

During the successful build, Fortran/C linker warnings were observed (for example COMMON block size and symbol-size warnings).
In this setup they were warnings only; the final `q2p1_ATC` executable linked successfully.
