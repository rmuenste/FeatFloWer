# FeatFloWer (Q2P1) Build Guide

This document provides instructions for configuring and building the FeatFloWer (internally referred to as Q2P1) software project using CMake.

## From-Scratch Guide Series

For end-to-end, command-by-command workflows (configure/build/stage/run), see:

* `docs/md_docs/guide_01_q2p1_fc_ext_cylinder_benchmark_from_scratch.md`
* `docs/md_docs/guide_02_q2p1_bench_sedimentation_pe_serial_from_scratch.md`
* `docs/md_docs/guide_03_q2p1_sse_tse_gendie_from_scratch.md`

## 1. Prerequisites

Before you begin, ensure you have the following installed on your system:

*   **CMake**: Version 3.18 or higher.
*   **C, C++, and Fortran Compilers**: A compatible set (e.g., GCC, Intel, Clang). The build system has mechanisms to select compiler sets (see Build IDs).
*   **Git**: Required if you plan to build CGAL from its repository (i.e., not using a local installation). Version 2.17 or higher is checked for.
*   **OpenMP**: This is a required dependency. Ensure your compilers support OpenMP.
*   **librt**: On Linux systems, the `librt` library is required (often part of `glibc-devel` or `libc6-dev`).
*   **(Optional) MPI Implementation**: If you plan to use features requiring MPI (like PE or MUMPS), an MPI implementation (e.g., OpenMPI, MPICH, Intel MPI) must be available.
*   **(Optional) Boost Libraries**: Required if `USE_BOOST` is enabled (e.g., by `USE_CGAL`). Version 1.68.0 or newer.
*   **(Optional) MKL (Math Kernel Library)**: If `USE_MUMPS` is enabled with an Intel compiler.

## 2. Getting Started: Basic Build

It is **strongly recommended** to perform an out-of-source build. In-source builds are explicitly disallowed by the build system.

1.  **Create a build directory:**
    ```bash
    mkdir build
    cd build
    ```

2.  **Configure using CMake:**
    The simplest configuration points CMake to the source directory (parent directory `..` in this case):
    ```bash
    cmake ..
    ```
    This will use default settings, including building bundled versions of dependencies like BLAS/LAPACK and METIS, and disabling optional libraries like CGAL, PE, etc.

3.  **Compile the project:**
    ```bash
    cmake --build . -- -j<N>
    # OR, if you used "Unix Makefiles" generator:
    # make -j<N>
    ```
    Replace `<N>` with the number of parallel jobs you want to use for compilation (e.g., `make -j8`).

## 3. Key CMake Configuration Variables

You can customize the build by passing `-D<VARIABLE_NAME>=<VALUE>` options to the `cmake` command during the configuration step.

*   **`CMAKE_INSTALL_PREFIX`**:
    Specifies the directory where the project will be installed when you run `cmake --build . --target install` (or `make install`).
    ```bash
    cmake -DCMAKE_INSTALL_PREFIX=/path/to/your/install/dir ..
    ```
    If not specified, it defaults to a system-dependent location (e.g., `/usr/local`).

*   **`CMAKE_BUILD_TYPE`**:
    Determines the build configuration and associated compiler flags. Common values:
    *   `Release`: (Default) Optimized for performance, no debug symbols.
    *   `Debug`: Includes debug symbols, fewer optimizations.
    *   `RelWithDebInfo`: Optimized build with debug symbols.
    *   `MinSizeRel`: Optimized for smallest binary size.
    ```bash
    cmake -DCMAKE_BUILD_TYPE=Debug ..
    ```

*   **`-G "<GeneratorName>"`**:
    Specifies the build system generator. Examples:
    *   `"Unix Makefiles"` (Default on Linux/macOS)
    *   `"Ninja"`
    ```bash
    cmake -G "Ninja" ..
    ```

*   **`CMAKE_C_COMPILER`, `CMAKE_CXX_COMPILER`, `CMAKE_Fortran_COMPILER`**:
    Specify the C, C++, and Fortran compilers respectively.
    ```bash
    cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran ..
    ```

*   **`Q2P1_COMPILER` (Advanced)**:
    This project-specific variable can be used to select a predefined compiler suite/toolchain via internal build ID mechanisms. It defaults to `gcc` if not set.
    ```bash
    cmake -DQ2P1_COMPILER=intel ..
    ```
    This often works in conjunction with "Build IDs" (see below).

## 4. Available Build Options (`-DOPTION_NAME=ON/OFF`)

These options control which components and features are enabled or disabled.

| Option                          | Default | Description                                                                                                | Notes                                                                                                 |
| :------------------------------ | :------ | :--------------------------------------------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------- |
| `USE_OPENMESH`                  | `OFF`   | Use the OpenMesh library.                                                                                  | Automatically enabled if `BUILD_BOUNDARY_LAYER_TOOLS` is `ON`.                                        |
| `USE_MUMPS`                     | `OFF`   | Use the MUMPS MPI parallel direct coarse grid solver.                                                      | **Requires an Intel compiler toolchain**. Adds `MumpsSolver.f90`.                                    |
| `USE_ODE`                       | `OFF`   | Use the ODE (Open Dynamics Engine) library.                                                                | If `BUILD_APPLICATIONS` is `ON`, adds `FullC0ntact/libs/ode-cmake`.                                  |
| `USE_CGAL`                      | `OFF`   | Use the CGAL (Computational Geometry Algorithms Library).                                                  | See "Building with CGAL" section. Automatically enables `USE_BOOST`.                                 |
| `USE_CGAL_LOCAL`                | `OFF`   | If `USE_CGAL` is `ON`, try to find a local system installation of CGAL instead of building it from source. |                                                                                                       |
| `USE_BOOST`                     | `OFF`   | Use Boost C++ libraries.                                                                                   | Automatically enabled if `USE_CGAL` is `ON`. Requires version 1.68.0+.                                |
| `USE_SYSTEM_BLASLAPACK`         | `OFF`   | Use system-provided BLAS/LAPACK libraries if found.                                                        | If `OFF` or not found, bundled `lapack-3.6.1` is built.                                                 |
| `BUILD_METIS`                   | `ON`    | Enable build of the bundled METIS 4.0.3 library.                                                           |                                                                                                       |
| `BUILD_APPLICATIONS`            | `ON`    | Enable build of applications.                                                                              | Requires the `FullC0ntact` subdirectory.                                                              |
| `BUILD_BOUNDARY_LAYER_TOOLS`    | `OFF`   | Build the boundary layer tools.                                                                            | Forces `USE_OPENMESH=ON`.                                                                             |
| `USE_HYPRE`                     | `OFF`   | Use the HYPRE (High Performance Preconditioners) library.                                                  | Builds bundled HYPRE. Adds `HypreSolver.f90`.                                                        |
| `USE_PE`                        | `OFF`   | Use the PE (rigid body physics engine) library.                                                            | Builds `libs/pe`. Enables `MPI` option.                                                              |
| `USE_PE_SERIAL_MODE`            | `OFF`   | Use serial PE mode for large particles that span multiple domains.                                         | Only available when `USE_PE=ON`. See section 5.2.                                                    |
| `SED_BENCH`                     | `OFF`   | Enable sedimentation benchmark output (force/position/velocity per timestep).                              | Adds `-DSED_BENCH` preprocessor flag; activates `SED_BENCH_FORCE/POS/VEL` log lines.                |
| `VERIFY_HASHGRID`               | `OFF`   | Run HashGrid alpha acceleration alongside brute-force baseline for correctness verification.               | **Significant performance overhead.** For debugging only; do not use in production runs.             |
| `SHOW_BUILD_IDS`                | `OFF`   | Display all available build IDs (compiler/flag configurations) and exit.                                   | Useful for developers to see predefined toolchain settings.                                           |

## 5. Specific Feature Configurations

### 5.1. Building with CGAL

CGAL is a powerful library for computational geometry.

*   **Building CGAL from Source (Default with `USE_CGAL=ON`)**:
    When `USE_CGAL=ON` and `USE_CGAL_LOCAL=OFF` (default), the build system will download and compile CGAL (version v5.3.2) as an external project. This also automatically enables `USE_BOOST=ON`.
    ```bash
    cmake -DUSE_CGAL=ON ..
    ```

*   **Using a Local/System CGAL Installation**:
    If you have CGAL already installed on your system:
    ```bash
    cmake -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=ON ..
    ```
    The build system will attempt to find CGAL in standard system locations or via the `CGAL_DIR` environment variable.

*   **Cluster-specific note**: Linking can fail due to static/dynamic Boost conflicts. If so, try:
    ```bash
    cmake -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=ON -DBoost_USE_STATIC_LIBS=ON ..
    ```

### 5.2. Building with PE (Physics Engine)

The PE library enables rigid body dynamics coupled to the CFD solver via the Fictitious Boundary Method (FBM). Two modes are supported:

#### Parallel PE (default when `USE_PE=ON`)

Full MPI parallelization within the PE library. Best for many small/medium particles distributed across domains.

```bash
cmake -DUSE_PE=ON ..
```

Ensure a working MPI implementation is installed and findable by CMake. You can specify the MPI vendor via `-DMPI_VENDOR=openmpi`.

#### Serial PE (`USE_PE=ON` + `USE_PE_SERIAL_MODE=ON`)

Each CFD domain runs an independent serial PE instance. Forces are synchronized via the CFD MPI layer (bypasses PE's own MPI). Best for a small number of large particles (< ~20) that span multiple subdomains.

**Use the two-step configuration** to ensure correct dependency resolution:

```bash
# Step 1: Enable PE library
cmake -DUSE_PE=ON ..

# Step 2: Enable serial mode
cmake -DUSE_PE_SERIAL_MODE=ON ..

# Build
make -j8
```

`USE_PE_SERIAL_MODE` is only visible to CMake after `USE_PE=ON` is set, due to the `cmake_dependent_option` mechanism.

### 5.3. Building with PE and CGAL

```bash
cmake -DUSE_PE=ON -DUSE_CGAL=ON ..
# Add -DUSE_CGAL_LOCAL=ON if you have a local CGAL installation
```

### 5.4. Building with MUMPS

MUMPS is a parallel direct solver that **requires an Intel compiler toolchain**.

```bash
cmake -DUSE_MUMPS=ON -DQ2P1_COMPILER=intel ..
```

Ensure MKL (Intel Math Kernel Library) is also available and findable by CMake.

### 5.5. Building with HYPRE

```bash
cmake -DUSE_HYPRE=ON ..
```

This will build the bundled version of HYPRE.

### 5.6. Applications (`BUILD_APPLICATIONS=ON`)

By default, applications are built. This requires the `FullC0ntact` subdirectory to be present in the source tree.

### 5.7. BLAS/LAPACK Configuration

*   **Bundled (default)**: The project builds `lapack-3.6.1` from `extern/libraries/lapack-3.6.1`.
*   **System**: To use BLAS/LAPACK already installed on your system:
    ```bash
    cmake -DUSE_SYSTEM_BLASLAPACK=ON ..
    ```

## 6. Acceleration Techniques (PE Serial Mode)

When using PE serial mode (`USE_PE=ON` + `USE_PE_SERIAL_MODE=ON`), two independent acceleration layers reduce the cost of the FBM coupling each timestep.

### 6.1. Runtime Control via Parameter File

Both accelerations can be controlled at runtime via `_data/q2p1_param.dat`:

```
SimPar@UseHashGridAccel = Yes
SimPar@UseKVELAccel = Yes
```

Default values (when not specified): both `Yes` (enabled).

To disable an acceleration for testing, set to `No`:
```
SimPar@UseHashGridAccel = No   # Force brute-force alpha field computation
SimPar@UseKVELAccel = No       # Force brute-force force integration
```

### 6.2. Alpha Field Acceleration (HashGrid)

The FBM geometry query `fbm_getFictKnprFC2` — which classifies every mesh DOF as fluid or solid — is accelerated via a spatial HashGrid built by the PE collision pipeline.

*   **Runtime control:** `SimPar@UseHashGridAccel` in `q2p1_param.dat`
*   On timestep 1 the brute-force baseline (`verifyAllParticles`) is used because the HashGrid is not yet built.
*   From timestep 2 onwards, `checkAllParticles` (HashGrid lookup) is called if enabled.
*   To verify correctness, build with `-DVERIFY_HASHGRID=ON`, which runs both paths in parallel and reports any mismatches (at significant cost).

### 6.3. Force Integration Acceleration (KVEL/KEEL/KAAL candidate elements)

The hydrodynamic force integral over boundary elements is accelerated by restricting the element loop to a candidate set: only elements topologically adjacent to DOFs inside the particle.

*   **Runtime control:** `SimPar@UseKVELAccel` in `q2p1_param.dat`
*   When disabled, the full brute-force element loop is used — identical physics, useful for correctness comparison.
*   Can achieve 10,000x+ reduction in element evaluations for many-particle systems.

### 6.4. Pure Baseline Builds (Compile-Time Disable)

For final correctness verification, you can compile out all acceleration code:

```bash
cmake -DUSE_PE=ON -DUSE_PE_SERIAL_MODE=ON -DENABLE_FBM_ACCELERATION=OFF ..
make -j8
```

This produces a pure brute-force binary where the acceleration code paths are never compiled. The runtime flags in `q2p1_param.dat` have no effect in this build.

### 6.5. Summary

| Layer | Runtime Control | Compile-Time Control | Default (runtime) | Default (compile) |
|---|---|---|---|---|
| Alpha / HashGrid | `SimPar@UseHashGridAccel` | `-DENABLE_FBM_ACCELERATION` | ON | ON |
| Force / KVEL | `SimPar@UseKVELAccel` | `-DENABLE_FBM_ACCELERATION` | ON | ON |
| HashGrid verification | N/A (debug only) | `-DVERIFY_HASHGRID=ON` | N/A | OFF |

**Required CMake flags:** `-DUSE_PE=ON -DUSE_PE_SERIAL_MODE=ON`

**Recommended testing workflow:**
1. Build with `-DENABLE_FBM_ACCELERATION=ON` (default)
2. Run with both flags `Yes` (accelerated)
3. Run with both flags `No` (brute-force, same binary)
4. Compare results to verify correctness
5. For final verification: rebuild with `-DENABLE_FBM_ACCELERATION=OFF` and compare

## 7. Build IDs and Compiler Selection (`Q2P1_BUILD_ID`)

The project includes a system for managing "Build IDs" corresponding to specific compiler versions, flags, and environments.

*   **To list available Build IDs**:
    ```bash
    cmake -DSHOW_BUILD_IDS=ON ..
    ```

*   **To use a specific Build ID**:
    ```bash
    cmake -DQ2P1_COMPILER=intel ..
    ```

## 8. Installation

After a successful build:
```bash
cmake --build . --target install
# OR:
# make install
```

Installed targets:
*   `tools/partitioner` → `<CMAKE_INSTALL_PREFIX>/bin`
*   `tools/PyPartitioner.py` → `<CMAKE_INSTALL_PREFIX>/bin`

## 9. Packaging (CPack)

```bash
cpack
```

Run from the build directory. Generates Debian packages (`.deb`) as configured by `cmake/modules/GeneratePackageDeb.cmake`.

## 10. Example Build Commands

*   **Basic Release Build:**
    ```bash
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=~/featflower_install ..
    make -j8
    make install
    ```

*   **PE Serial Mode (large particle sedimentation):**
    ```bash
    mkdir build_pe_serial && cd build_pe_serial
    cmake -DUSE_PE=ON ..
    cmake -DUSE_PE_SERIAL_MODE=ON ..
    make -j8
    ```

*   **PE Serial Mode with sedimentation benchmark output:**
    ```bash
    cmake -DUSE_PE=ON ..
    cmake -DUSE_PE_SERIAL_MODE=ON -DSED_BENCH=ON ..
    make -j8
    ```

*   **Debug Build with CGAL and PE, using Ninja:**
    ```bash
    mkdir build_debug && cd build_debug
    cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Debug -DUSE_CGAL=ON -DUSE_PE=ON ..
    ninja
    ```

*   **Build with MUMPS (requires Intel toolchain):**
    ```bash
    mkdir build_intel_mumps && cd build_intel_mumps
    cmake -DQ2P1_COMPILER=intel -DUSE_MUMPS=ON ..
    make -j8
    ```
