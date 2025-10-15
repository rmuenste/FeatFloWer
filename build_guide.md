# FeatFloWer (Q2P1) Build Guide

This document provides instructions for configuring and building the FeatFloWer (internally referred to as Q2P1) software project using CMake. The following sections describe a FeatFloWer version that is >=

commit 61644df53ca3eb534626f6ec1a4e1fca01f46377 (HEAD -> master, origin/master)
Author: superman <raphael.muenster@math.tu-dortmund.de>
Date:   Wed Apr 23 17:15:48 2025 +0200

    Forwarded submodule pe

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
    *   `"Visual Studio 17 2022"` (Windows)
    ```bash
    cmake -G "Ninja" ..
    ```
    Using generators other than "Unix Makefiles" is at your own risk, though common ones like Ninja should work.

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
| `USE_MUMPS`                     | `OFF`   | Use the MUMPS MPI parallel direct coarse grid solver.                                                      | **Requires an Intel compiler toolchain** (selected via `Q2P1_BUILD_ID` or `Q2P1_COMPILER`). Adds `MumpsSolver.f90`. |
| `USE_ODE`                       | `OFF`   | Use the ODE (Open Dynamics Engine) library.                                                                | If `BUILD_APPLICATIONS` is `ON`, adds `FullC0ntact/libs/ode-cmake`.                                  |
| `USE_CGAL`                      | `OFF`   | Use the CGAL (Computational Geometry Algorithms Library).                                                  | See "Building with CGAL" section. Automatically enables `USE_BOOST`.                                 |
| `USE_CGAL_LOCAL`                | `OFF`   | If `USE_CGAL` is `ON`, try to find a local system installation of CGAL instead of building it from source. |                                                                                                       |
| `USE_BOOST`                     | `OFF`   | Use Boost C++ libraries.                                                                                   | Automatically enabled if `USE_CGAL` is `ON`. Requires version 1.68.0+.                                |
| `USE_SYSTEM_BLASLAPACK`         | `OFF`   | Use system-provided BLAS/LAPACK libraries if found.                                                        | If `OFF` or not found, bundled `lapack-3.6.1` is built.                                                 |
| `BUILD_METIS`                   | `ON`    | Enable build of the bundled METIS 4.0.3 library.                                                           |                                                                                                       |
| `BUILD_APPLICATIONS`            | `ON`    | Enable build of applications.                                                                              | Requires the `FullC0ntact` subdirectory.                                                     |
| `BUILD_BOUNDARY_LAYER_TOOLS`    | `OFF`   | Build the boundary layer tools.                                                                            | Forces `USE_OPENMESH=ON`.                                                                             |
| `USE_HYPRE`                     | `OFF`   | Use the HYPRE (High Performance Preconditioners) library.                                                  | Builds bundled HYPRE. Adds `HypreSolver.f90`.                                                        |
| `USE_PE`                        | `OFF`   | Use the PE (Parallel Environment) library.                                                                 | Builds `libs/pe`. Enables `MPI` option.                                                              |
| `SHOW_BUILD_IDS`                | `OFF`   | Display all available build IDs (compiler/flag configurations) and exit.                                   | Useful for developers to see predefined toolchain settings.                                           |

## 5. Specific Feature Configurations

### 5.1. Building with CGAL

CGAL is a powerful library for computational geometry.

*   **Building CGAL from Source (Default with `USE_CGAL=ON`)**:
    When `USE_CGAL=ON` and `USE_CGAL_LOCAL=OFF` (default for `USE_CGAL_LOCAL`), the build system will download and compile CGAL (version v5.3.2 specified) as an external project. This also automatically enables `USE_BOOST=ON`.
    ```bash
    cmake -DUSE_CGAL=ON ..
    ```

*   **Using a Local/System CGAL Installation**:
    If you have CGAL already installed on your system and want to use it, enable both `USE_CGAL` and `USE_CGAL_LOCAL`:
    ```bash
    cmake -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=ON ..
    ```
    The build system will attempt to find CGAL in standard system locations or via environment variables like `CGAL_DIR`.
*   **Cluster-specific Adjustements**:
    CGAL needs the Boost libraries. Linking on some clusters can fail due to static/dynamic linking problems. In case of said problems try to use:

    ```bash
    cmake -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=ON -DBoost_USE_STATIC_LIBS=ON ..
    ```    

### 5.2. Building with PE (Parallel Environment)

The PE library enables parallel processing capabilities, likely using MPI.
To enable PE:
```bash
cmake -DUSE_PE=ON ..
```
This will automatically enable the `MPI` option. Ensure you have a working MPI implementation installed and findable by CMake (e.g., `mpicc`, `mpif90` in your `PATH`). You can specify the MPI vendor via `MPI_VENDOR` (e.g., `-DMPI_VENDOR=openmpi`).

**Cluster-specific Adjustements**:
CGAL needs the Boost libraries. Linking on some clusters can fail due to static/dynamic linking problems. In case of said problems try to use:
```bash
    cmake -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=ON -DBoost_USE_STATIC_LIBS=ON ..
```    

### 5.3. Building with PE and CGAL

This combination is supported for advanced geometry components.
To enable both PE and CGAL:
```bash
cmake -DUSE_PE=ON -DUSE_CGAL=ON ..
# Optionally add -DUSE_CGAL_LOCAL=ON if you have a local CGAL and want to use it
```
**Cluster-specific Adjustements**:
CGAL needs the Boost libraries. Linking on some clusters can fail due to static/dynamic linking problems. In case of said problems try to use:
```bash
cmake -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=ON -DBoost_USE_STATIC_LIBS=ON ..
```    

### 5.4. Building with MUMPS

MUMPS is a parallel direct solver that, in this project, **requires an Intel compiler toolchain**.
To enable MUMPS (assuming an Intel toolchain is active or will be selected):
```bash
cmake -DUSE_MUMPS=ON -DQ2P1_COMPILER=intel .. # Or appropriate Q2P1_BUILD_ID for Intel
```
If a non-Intel compiler toolchain is active when `USE_MUMPS=ON`, CMake will produce a fatal error. Ensure MKL (Intel Math Kernel Library) is also available and findable by CMake.

### 5.5. Building with HYPRE

HYPRE provides high-performance preconditioners.
To enable HYPRE:
```bash
cmake -DUSE_HYPRE=ON ..
```
This will build the bundled version of HYPRE.

### 5.6. Applications (`BUILD_APPLICATIONS=ON`)

By default, applications are built (`BUILD_APPLICATIONS=ON`). This requires the `FullC0ntact` subdirectory to be present in the source tree.
*   If `USE_OPENMESH=ON`, `FullC0ntact/libs/OpenMesh` is added.
*   If `USE_ODE=ON`, `FullC0ntact/libs/ode-cmake` is added.
*   If `USE_CGAL=ON`, the `inshape3dcore` target (likely part of `FullC0ntact` or `applications`) will have a dependency on the `cgal` external project target.

### 5.7. BLAS/LAPACK Configuration

*   **Bundled BLAS/LAPACK (Default)**:
    The project builds `lapack-3.6.1` from `extern/libraries/lapack-3.6.1`. This is the default behavior.
    ```bash
    cmake ..
    # or explicitly
    # cmake -DUSE_SYSTEM_BLASLAPACK=OFF ..
    ```
*   **System BLAS/LAPACK**:
    To use BLAS/LAPACK libraries already installed on your system:
    ```bash
    cmake -DUSE_SYSTEM_BLASLAPACK=ON ..
    ```
    CMake will attempt to find them. If not found, it falls back to the bundled version.

## 6. Build IDs and Compiler Selection (`Q2P1_BUILD_ID`)

The project includes a system for managing "Build IDs" (e.g., via `GenerateBuildId.cmake`). These IDs likely correspond to specific compiler versions, flags, and environments (e.g., a "gcc" build ID, an "intel" build ID).

*   **To list available Build IDs**:
    ```bash
    cmake -DSHOW_BUILD_IDS=ON ..
    ```
    This will display the IDs and then exit.

*   **To use a specific Build ID**:
    The exact mechanism to select a build ID (e.g., `-DQ2P1_BUILD_ID=<ID_NAME>`) would be detailed by the output of `SHOW_BUILD_IDS` or the `CompilerQuickCheck.cmake` and `GenerateBuildId.cmake` modules. The variable `Q2P1_COMPILER` seems to be a primary way to influence this compiler selection, which in turn selects the build ID. For example:
    ```bash
    cmake -DQ2P1_COMPILER=intel ..
    ```
    This is particularly important for features like MUMPS that depend on a specific compiler toolchain.

## 7. Installation

After a successful build, you can install the project using:
```bash
cmake --build . --target install
# OR, if using Makefiles:
# make install
```
This will install files to the path specified by `CMAKE_INSTALL_PREFIX`.
The following are explicitly installed:
*   `tools/partitioner` directory to `<CMAKE_INSTALL_PREFIX>/bin`
*   `tools/PyPartitioner.py` script to `<CMAKE_INSTALL_PREFIX>/bin` with execute permissions.

## 8. Packaging (CPack)

The project is configured to generate Debian packages (`.deb`) using CPack.
After building the project, you can typically generate a package by running:
```bash
cpack
```
This command should be run from the build directory. The specifics of package generation are controlled by `cmake/modules/GeneratePackageDeb.cmake` and CPack variables that might be set elsewhere.

## 9. Example Build Commands

*   **Basic Release Build (Unix Makefiles):**
    ```bash
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=~/featflower_install ..
    make -j8
    make install
    ```

*   **Debug Build with CGAL (from source) and PE, using Ninja:**
    ```bash
    mkdir build_debug_cgal_pe && cd build_debug_cgal_pe
    cmake -G "Ninja" \
          -DCMAKE_BUILD_TYPE=Debug \
          -DUSE_CGAL=ON \
          -DUSE_PE=ON \
          -DCMAKE_INSTALL_PREFIX=~/featflower_install_debug \
          ..
    ninja
    ninja install
    ```

*   **Build with MUMPS (requires Intel toolchain):**
    ```bash
    mkdir build_intel_mumps && cd build_intel_mumps
    cmake -DQ2P1_COMPILER=intel \
          -DUSE_MUMPS=ON \
          -DCMAKE_INSTALL_PREFIX=~/featflower_intel_mumps \
          ..
    make -j8
    ```

Remember to adapt compiler settings (`CMAKE_<LANG>_COMPILER` or `Q2P1_COMPILER`/`Q2P1_BUILD_ID`) as needed for your system and desired configuration.
```