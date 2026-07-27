#=========================================================================
# FFDependencies.cmake
#
# Provider pattern for the solver dependencies of FeatFloWer.
#
# Each dependency is resolved by exactly one provider:
#
#   vendored  - the copy shipped in extern/libraries (legacy default)
#   external  - a site install named by the legacy EXTERNAL_*_DIR variables
#   system    - whatever find_package() locates (fails loudly if absent)
#   fetch     - the pinned tarball, downloaded and built as a sub-project
#   auto      - external if the legacy option is set, else vendored if the
#               vendored copy exists, else fetch (find_package still wins
#               through FIND_PACKAGE_ARGS)
#
# The result of every provider is funnelled into a thin INTERFACE target
# (ff::hypre, ff::umfpack) plus the legacy variables the rest of the build
# system already consumes (HYPRE_LIBRARIES, HYPRE_STRUCTURES_INCLUDE_PATH,
# FF_UMFPACK_LINK_LIBS, FF_SUITESPARSE_INCLUDE_DIR, FF_MKL_PARDISO_LIBS), so
# GenerateLinkerFlags.cmake / GenerateIncludeFlags.cmake stay unchanged.
#
# See docs/md_docs/fetchcontent_dependencies.md.
#=========================================================================

include_guard(GLOBAL)
include(FetchContent)

# Extracted tarball timestamps: use the extraction time so that re-running
# CMake does not needlessly re-configure the sub-projects (CMake >= 3.24).
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

#-------------------------------------------------------------------------
# Pinned upstream releases (URL + SHA256 of the actual release tarball)
#-------------------------------------------------------------------------
set(FF_HYPRE_FETCH_VERSION "2.33.0")
set(FF_HYPRE_FETCH_URL
  "https://github.com/hypre-space/hypre/archive/refs/tags/v2.33.0.tar.gz")
set(FF_HYPRE_FETCH_SHA256
  "0f9103c34bce7a5dcbdb79a502720fc8aab4db9fd0146e0791cde7ec878f27da")

set(FF_SUITESPARSE_FETCH_VERSION "7.12.2")
set(FF_SUITESPARSE_FETCH_URL
  "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.12.2.tar.gz")
set(FF_SUITESPARSE_FETCH_SHA256
  "679412daa5f69af96d6976595c1ac64f252287a56e98cc4a8155d09cc7fd69e8")

# The minimum CMake version that supports FetchContent's FIND_PACKAGE_ARGS.
set(FF_FETCHCONTENT_MIN_VERSION "3.24")

#-------------------------------------------------------------------------
# User-facing options
#-------------------------------------------------------------------------
option(FF_FETCH_DEPENDENCIES
  "Allow FetchContent to download and build missing solver dependencies"
  ON
  )

set(FF_HYPRE_PROVIDER "auto" CACHE STRING
  "How to obtain hypre: auto|vendored|external|system|fetch")
set_property(CACHE FF_HYPRE_PROVIDER PROPERTY STRINGS
  auto vendored external system fetch)

set(FF_SUITESPARSE_PROVIDER "auto" CACHE STRING
  "How to obtain SuiteSparse/UMFPACK: auto|vendored|external|system|fetch")
set_property(CACHE FF_SUITESPARSE_PROVIDER PROPERTY STRINGS
  auto vendored external system fetch)

option(FF_HYPRE_GPU_AWARE_MPI
  "Build/require hypre with GPU-aware MPI (HYPRE_ENABLE_GPU_AWARE_MPI)"
  OFF
  )

# OFF by default on purpose: MKL::MKL moves MKL to the end of the link line,
# which on this toolchain splits the Fortran runtime and hangs the solver at
# startup. See the long comment in ff_setup_mkl_pardiso() below.
option(FF_MKL_USE_CMAKE_CONFIG
  "Link MKL through MKLConfig.cmake (MKL::MKL) instead of the explicit static \
link group. Changes MKL's position on the link line; read the libgfortran note \
in ff_setup_mkl_pardiso() before enabling."
  OFF
  )

#-------------------------------------------------------------------------
# Helpers
#-------------------------------------------------------------------------

# Abort with a helpful message when a fetch is required but not permitted.
macro(_ff_assert_fetch_allowed _dep)
  if(NOT FF_FETCH_DEPENDENCIES)
    message(FATAL_ERROR
      "${_dep} would have to be downloaded (FetchContent), but "
      "FF_FETCH_DEPENDENCIES=OFF. Provide a site install instead "
      "(see docs/md_docs/fetchcontent_dependencies.md).")
  endif()
  if(CMAKE_VERSION VERSION_LESS ${FF_FETCHCONTENT_MIN_VERSION})
    message(FATAL_ERROR
      "Fetching ${_dep} needs CMake >= ${FF_FETCHCONTENT_MIN_VERSION} "
      "(FetchContent FIND_PACKAGE_ARGS); this is CMake ${CMAKE_VERSION}. "
      "Use a site install and set FF_${_dep}_PROVIDER=external/system.")
  endif()
endmacro()

# FeatFloWer ships Find*.cmake modules (FindBLAS, FindLAPACK, FindMKL, ...)
# that are tailored to this project. Third-party sub-projects configured
# through FetchContent must not pick them up instead of CMake's own modules,
# so the module path is emptied for the duration of the sub-configure.
macro(_ff_push_clean_module_path)
  set(_ff_saved_module_path "${CMAKE_MODULE_PATH}")
  set(CMAKE_MODULE_PATH "")
endmacro()

macro(_ff_pop_module_path)
  set(CMAKE_MODULE_PATH "${_ff_saved_module_path}")
  unset(_ff_saved_module_path)
endmacro()

#=========================================================================
# hypre
#
# Sets: HYPRE_LIBRARIES, HYPRE_STRUCTURES_INCLUDE_PATH,
#       FF_HYPRE_PROVIDER_RESOLVED, target ff::hypre
#=========================================================================
macro(ff_setup_hypre)

  set(_ff_hypre_vendored_dir "${CMAKE_SOURCE_DIR}/extern/libraries/hypre/src")

  # ---- resolve the provider ---------------------------------------------
  set(_ff_hypre_provider "${FF_HYPRE_PROVIDER}")
  if(_ff_hypre_provider STREQUAL "auto")
    if(USE_EXTERNAL_HYPRE)
      set(_ff_hypre_provider "external")
    elseif(EXISTS "${_ff_hypre_vendored_dir}/CMakeLists.txt")
      set(_ff_hypre_provider "vendored")
    else()
      set(_ff_hypre_provider "fetch")
    endif()
  elseif(_ff_hypre_provider STREQUAL "external" AND NOT USE_EXTERNAL_HYPRE)
    # Explicitly asking for the external provider implies the legacy option.
    set(USE_EXTERNAL_HYPRE ON CACHE BOOL
      "Link a pre-installed hypre instead of the vendored copy" FORCE)
  endif()

  message(STATUS "hypre provider: ${_ff_hypre_provider}")
  set(FF_HYPRE_PROVIDER_RESOLVED "${_ff_hypre_provider}")

  add_library(ff_hypre INTERFACE)
  add_library(ff::hypre ALIAS ff_hypre)

  # ---- external: unchanged legacy path ----------------------------------
  if(_ff_hypre_provider STREQUAL "external")

    if(NOT EXTERNAL_HYPRE_DIR)
      message(FATAL_ERROR
        "USE_EXTERNAL_HYPRE=ON requires EXTERNAL_HYPRE_DIR to name an install prefix")
    endif()

    message(STATUS "Using external hypre from ${EXTERNAL_HYPRE_DIR}")
    find_package(HYPRE 2.31 CONFIG REQUIRED
      PATHS "${EXTERNAL_HYPRE_DIR}"
      NO_DEFAULT_PATH
      )
    target_link_libraries(ff_hypre INTERFACE HYPRE::HYPRE)
    set(HYPRE_STRUCTURES_INCLUDE_PATH "${EXTERNAL_HYPRE_DIR}/include")

  # ---- system: any find_package-able hypre ------------------------------
  elseif(_ff_hypre_provider STREQUAL "system")

    find_package(HYPRE 2.31 CONFIG REQUIRED)
    target_link_libraries(ff_hypre INTERFACE HYPRE::HYPRE)
    get_target_property(_ff_hypre_inc HYPRE::HYPRE INTERFACE_INCLUDE_DIRECTORIES)
    if(_ff_hypre_inc)
      set(HYPRE_STRUCTURES_INCLUDE_PATH "${_ff_hypre_inc}")
    endif()
    unset(_ff_hypre_inc)

  # ---- vendored: unchanged legacy path ----------------------------------
  elseif(_ff_hypre_provider STREQUAL "vendored")

    message(STATUS "Using repository version of hypre")
    if(USE_HYPRE_CUDA)
      set(HYPRE_WITH_CUDA ON CACHE BOOL "" FORCE)
      set(HYPRE_ENABLE_UNIFIED_MEMORY ON CACHE BOOL "" FORCE)
      set(HYPRE_CUDA_SM "80" CACHE STRING "" FORCE)
    endif()
    add_subdirectory("${_ff_hypre_vendored_dir}")
    target_link_libraries(ff_hypre INTERFACE HYPRE)
    set(HYPRE_STRUCTURES_INCLUDE_PATH "${_ff_hypre_vendored_dir}")

  # ---- fetch ------------------------------------------------------------
  elseif(_ff_hypre_provider STREQUAL "fetch")

    _ff_assert_fetch_allowed(HYPRE)

    ff_configure_fetched_hypre_options()

    FetchContent_Declare(HYPRE
      URL "${FF_HYPRE_FETCH_URL}"
      URL_HASH SHA256=${FF_HYPRE_FETCH_SHA256}
      SOURCE_SUBDIR src
      FIND_PACKAGE_ARGS 2.31 NAMES HYPRE
      )

    _ff_push_clean_module_path()
    FetchContent_MakeAvailable(HYPRE)
    _ff_pop_module_path()

    target_link_libraries(ff_hypre INTERFACE HYPRE::HYPRE)
    if(hypre_SOURCE_DIR)
      # HYPREf.h (the Fortran header included by source/HypreSolver.f90)
      # lives next to the sub-project's CMakeLists.
      set(HYPRE_STRUCTURES_INCLUDE_PATH "${hypre_SOURCE_DIR}/src")
      message(STATUS "Fetched hypre ${FF_HYPRE_FETCH_VERSION}: ${hypre_SOURCE_DIR}")
    else()
      # find_package() satisfied the dependency through FIND_PACKAGE_ARGS.
      get_target_property(_ff_hypre_inc HYPRE::HYPRE INTERFACE_INCLUDE_DIRECTORIES)
      if(_ff_hypre_inc)
        set(HYPRE_STRUCTURES_INCLUDE_PATH "${_ff_hypre_inc}")
      endif()
      unset(_ff_hypre_inc)
      message(STATUS "hypre satisfied by find_package (no download needed)")
    endif()

  else()
    message(FATAL_ERROR
      "FF_HYPRE_PROVIDER='${FF_HYPRE_PROVIDER}' is not one of "
      "auto|vendored|external|system|fetch")
  endif()

  # ---- consistency checks for pre-built hypre ---------------------------
  if(_ff_hypre_provider MATCHES "^(external|system)$")
    if(USE_HYPRE_CUDA)
      if(NOT HYPRE_ENABLE_CUDA)
        message(FATAL_ERROR
          "USE_HYPRE_CUDA=ON, but external hypre was not built with CUDA")
      endif()
      if(NOT HYPRE_ENABLE_UNIFIED_MEMORY)
        message(FATAL_ERROR
          "USE_HYPRE_CUDA=ON requires external hypre with unified memory")
      endif()
    endif()
    if(FF_HYPRE_GPU_AWARE_MPI AND NOT HYPRE_ENABLE_GPU_AWARE_MPI)
      message(WARNING
        "FF_HYPRE_GPU_AWARE_MPI=ON, but the pre-built hypre does not report "
        "HYPRE_ENABLE_GPU_AWARE_MPI=ON. Rebuild hypre with GPU-aware MPI or "
        "switch to FF_HYPRE_PROVIDER=fetch.")
    endif()
  endif()

  set(HYPRE_LIBRARIES ff::hypre)
  unset(_ff_hypre_provider)
  unset(_ff_hypre_vendored_dir)

endmacro()

# Drive hypre's own build options from FeatFloWer-level settings. The values
# mirror the manual reference build in stage1-work/hypre-2.33.0-build.
macro(ff_configure_fetched_hypre_options)

  set(HYPRE_ENABLE_MPI          ON  CACHE BOOL "" FORCE)
  set(HYPRE_ENABLE_HYPRE_BLAS   ON  CACHE BOOL "" FORCE)
  set(HYPRE_ENABLE_HYPRE_LAPACK ON  CACHE BOOL "" FORCE)
  set(HYPRE_ENABLE_BIGINT       OFF CACHE BOOL "" FORCE)
  set(HYPRE_ENABLE_MIXEDINT     OFF CACHE BOOL "" FORCE)
  set(HYPRE_BUILD_EXAMPLES      OFF CACHE BOOL "" FORCE)
  set(HYPRE_BUILD_TESTS         OFF CACHE BOOL "" FORCE)

  # Enabler for the GPU-aware-MPI experiment (see solver_phase3 docs).
  set(HYPRE_ENABLE_GPU_AWARE_MPI ${FF_HYPRE_GPU_AWARE_MPI} CACHE BOOL "" FORCE)

  if(USE_HYPRE_CUDA)
    set(HYPRE_ENABLE_CUDA            ON CACHE BOOL "" FORCE)
    set(HYPRE_ENABLE_UNIFIED_MEMORY  ON CACHE BOOL "" FORCE)
    set(HYPRE_ENABLE_CUBLAS          ON CACHE BOOL "" FORCE)
    set(HYPRE_ENABLE_CURAND          ON CACHE BOOL "" FORCE)
    set(HYPRE_ENABLE_CUSPARSE        ON CACHE BOOL "" FORCE)
    set(HYPRE_ENABLE_CUSOLVER        ON CACHE BOOL "" FORCE)
    set(HYPRE_ENABLE_CUDA_STREAMS    ON CACHE BOOL "" FORCE)
    if(NOT CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES "80" CACHE STRING
        "CUDA architectures for the fetched hypre" FORCE)
    endif()
  else()
    set(HYPRE_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
  endif()

endmacro()

#=========================================================================
# SuiteSparse / UMFPACK
#
# Sets: FF_UMFPACK_LINK_LIBS, FF_SUITESPARSE_INCLUDE_DIR,
#       FF_UMFPACK_WRAPPER_PORT, FF_SUITESPARSE_PROVIDER_RESOLVED,
#       target ff::umfpack
#=========================================================================
macro(ff_setup_suitesparse)

  # ---- resolve the provider ---------------------------------------------
  set(_ff_ss_provider "${FF_SUITESPARSE_PROVIDER}")
  if(_ff_ss_provider STREQUAL "auto")
    if(USE_EXTERNAL_SUITESPARSE)
      set(_ff_ss_provider "external")
    elseif(EXISTS "${CMAKE_SOURCE_DIR}/extern/libraries/umfpack4/CMakeLists.txt")
      set(_ff_ss_provider "vendored")
    else()
      set(_ff_ss_provider "fetch")
    endif()
  elseif(_ff_ss_provider STREQUAL "external" AND NOT USE_EXTERNAL_SUITESPARSE)
    set(USE_EXTERNAL_SUITESPARSE ON CACHE BOOL
      "Use an external SuiteSparse UMFPACK instead of the vendored UMFPACK 4" FORCE)
  endif()

  message(STATUS "SuiteSparse provider: ${_ff_ss_provider}")
  set(FF_SUITESPARSE_PROVIDER_RESOLVED "${_ff_ss_provider}")

  add_library(ff_umfpack INTERFACE)
  add_library(ff::umfpack ALIAS ff_umfpack)

  # The vendored umf4_f77wrapper_port.c is only needed when the SuiteSparse
  # in use is newer than the vendored UMFPACK 4 (it re-implements the umf4*
  # handle table on top of the stable umfpack_di_* API).
  # NOTE: do not pre-set FF_SUITESPARSE_INCLUDE_DIR here. find_path() below
  # writes a *cache* entry, which a normal variable of the same name would
  # shadow -- that silently breaks the external provider.
  set(FF_UMFPACK_WRAPPER_PORT FALSE)

  # ---- external / system: unchanged legacy find-path --------------------
  if(_ff_ss_provider MATCHES "^(external|system)$")

    find_path(FF_SUITESPARSE_INCLUDE_DIR umfpack.h
      HINTS ${EXTERNAL_SUITESPARSE_DIR}/include/suitesparse
            ${EXTERNAL_SUITESPARSE_DIR}/include
      PATHS /usr/include/suitesparse
      )
    find_library(FF_UMFPACK_LIBRARY umfpack
      HINTS ${EXTERNAL_SUITESPARSE_DIR}/lib64 ${EXTERNAL_SUITESPARSE_DIR}/lib
      )
    if(NOT FF_SUITESPARSE_INCLUDE_DIR OR NOT FF_UMFPACK_LIBRARY)
      message(FATAL_ERROR
        "USE_EXTERNAL_SUITESPARSE=ON but umfpack.h or libumfpack was not found "
        "(EXTERNAL_SUITESPARSE_DIR='${EXTERNAL_SUITESPARSE_DIR}')")
    endif()
    message(STATUS "External UMFPACK library: ${FF_UMFPACK_LIBRARY}")
    message(STATUS "External UMFPACK headers: ${FF_SUITESPARSE_INCLUDE_DIR}")
    target_link_libraries(ff_umfpack INTERFACE ${FF_UMFPACK_LIBRARY})
    target_include_directories(ff_umfpack INTERFACE ${FF_SUITESPARSE_INCLUDE_DIR})
    set(FF_UMFPACK_WRAPPER_PORT TRUE)

  # ---- vendored: unchanged legacy path ----------------------------------
  elseif(_ff_ss_provider STREQUAL "vendored")

    add_subdirectory(extern/libraries/amd)
    add_subdirectory(extern/libraries/umfpack4)
    target_link_libraries(ff_umfpack INTERFACE amd umfpack4)

  # ---- fetch ------------------------------------------------------------
  elseif(_ff_ss_provider STREQUAL "fetch")

    _ff_assert_fetch_allowed(SUITESPARSE)

    ff_configure_fetched_suitesparse_options()

    FetchContent_Declare(SuiteSparse
      URL "${FF_SUITESPARSE_FETCH_URL}"
      URL_HASH SHA256=${FF_SUITESPARSE_FETCH_SHA256}
      FIND_PACKAGE_ARGS NAMES UMFPACK
      )

    # SuiteSparse runs find_package(BLAS); it must not see FeatFloWer's
    # vendored BLAS target name, nor FeatFloWer's own FindBLAS.cmake, or it
    # would link the reference netlib BLAS instead of an optimized one.
    set(_ff_saved_blas_libraries "${BLAS_LIBRARIES}")
    set(_ff_saved_lapack_libraries "${LAPACK_LIBRARIES}")
    set(_ff_saved_shared_libs "${BUILD_SHARED_LIBS}")
    unset(BLAS_LIBRARIES)
    unset(LAPACK_LIBRARIES)
    set(BUILD_SHARED_LIBS OFF)

    _ff_push_clean_module_path()
    FetchContent_MakeAvailable(SuiteSparse)
    _ff_pop_module_path()

    set(BUILD_SHARED_LIBS "${_ff_saved_shared_libs}")
    set(BLAS_LIBRARIES "${_ff_saved_blas_libraries}")
    set(LAPACK_LIBRARIES "${_ff_saved_lapack_libraries}")
    unset(_ff_saved_shared_libs)
    unset(_ff_saved_blas_libraries)
    unset(_ff_saved_lapack_libraries)

    if(NOT TARGET SuiteSparse::UMFPACK)
      message(FATAL_ERROR
        "SuiteSparse was made available but SuiteSparse::UMFPACK does not exist")
    endif()
    target_link_libraries(ff_umfpack INTERFACE SuiteSparse::UMFPACK)
    set(FF_UMFPACK_WRAPPER_PORT TRUE)
    if(suitesparse_SOURCE_DIR)
      set(FF_SUITESPARSE_INCLUDE_DIR
        "${suitesparse_SOURCE_DIR}/UMFPACK/Include"
        "${suitesparse_SOURCE_DIR}/SuiteSparse_config"
        "${suitesparse_SOURCE_DIR}/AMD/Include")
      message(STATUS
        "Fetched SuiteSparse ${FF_SUITESPARSE_FETCH_VERSION}: ${suitesparse_SOURCE_DIR}")
    else()
      message(STATUS "SuiteSparse satisfied by find_package (no download needed)")
    endif()

  else()
    message(FATAL_ERROR
      "FF_SUITESPARSE_PROVIDER='${FF_SUITESPARSE_PROVIDER}' is not one of "
      "auto|vendored|external|system|fetch")
  endif()

  set(FF_UMFPACK_LINK_LIBS ff::umfpack)
  unset(_ff_ss_provider)

endmacro()

macro(ff_configure_fetched_suitesparse_options)

  set(SUITESPARSE_ENABLE_PROJECTS
    "suitesparse_config;amd;camd;colamd;ccolamd;cholmod;umfpack"
    CACHE STRING "" FORCE)

  set(SUITESPARSE_USE_CUDA   OFF CACHE BOOL "" FORCE)
  set(SUITESPARSE_USE_PYTHON OFF CACHE BOOL "" FORCE)
  set(SUITESPARSE_DEMOS      OFF CACHE BOOL "" FORCE)
  set(BUILD_STATIC_LIBS      ON  CACHE BOOL "" FORCE)

endmacro()

#=========================================================================
# MKL PARDISO
#
# Sets: FF_MKL_PARDISO_LIBS  (+ the MKL_PARDISO_AVAIL definition)
#=========================================================================
macro(ff_setup_mkl_pardiso)

  if(NOT MKL_PARDISO_DIR)
    message(FATAL_ERROR
      "USE_MKL_PARDISO=ON requires MKL_PARDISO_DIR (MKL install root, "
      "e.g. /sfw/intel/2024.0.1/mkl/latest)")
  endif()

  # ---------------------------------------------------------------------
  # Why the hand-assembled list is the default (do not "clean this up"):
  #
  # FF_MKL_PARDISO_LIBS sits in the middle of FF_DEFAULT_LIBS /
  # FF_APPLICATION_LIBS, right after the vendored BLAS/LAPACK archives.
  # As a list of raw flags and .a paths it lands at exactly that position
  # in the generated link line. As an imported target (MKL::MKL) CMake
  # instead resolves it through the target graph and emits it near the END
  # of the link line.
  #
  # That reordering is not cosmetic here. This toolchain mixes the static
  # libgfortran.a from the /sfw GCC 13.2.0 prefix with the system
  # libgfortran.so.5, so which libgfortran.a members the linker pulls in
  # depends on the set of undefined symbols at the moment -lgfortran is
  # scanned. With MKL moved to the end, the executable ended up with
  # _gfortran_st_open/_gfortran_st_read resolved to statically linked
  # members while _gfortran_st_close resolved to the shared library --
  # i.e. two Fortran unit tables in one process. CLOSE then did not close
  # the unit that OPEN/READ were using, so the next sequential READ past
  # EOF returned iostat=5001 ("READ after EOF") instead of -1, and the
  # parameter-parser loops that exit only on -1 spun forever: every rank
  # deterministically hung at startup in GetPresParameters.
  #
  # Verified by relinking the identical object files with only the MKL
  # group moved back to its original position -- the hang disappears and
  # _gfortran_st_close becomes locally defined again.
  #
  # The MKLConfig.cmake route is kept behind the explicit opt-in
  # FF_MKL_USE_CMAKE_CONFIG (declared with the other options above) for
  # sites whose toolchain does not have the split-libgfortran hazard.
  # ---------------------------------------------------------------------
  set(_ff_mkl_config_dir "${MKL_PARDISO_DIR}/lib/cmake/mkl")
  set(FF_MKL_PARDISO_LIBS "")

  if(FF_MKL_USE_CMAKE_CONFIG AND EXISTS "${_ff_mkl_config_dir}/MKLConfig.cmake")
    set(MKL_LINK      "static"     CACHE STRING "" FORCE)
    set(MKL_INTERFACE "lp64"       CACHE STRING "" FORCE)
    set(MKL_THREADING "sequential" CACHE STRING "" FORCE)
    set(MKL_MPI       "openmpi"    CACHE STRING "" FORCE)

    _ff_push_clean_module_path()
    find_package(MKL CONFIG QUIET
      PATHS "${_ff_mkl_config_dir}"
      NO_DEFAULT_PATH
      )
    _ff_pop_module_path()

    if(TARGET MKL::MKL)
      # This project drives PARDISO from GNU Fortran, so the GNU interface
      # layer (libmkl_gf_*) is mandatory. If MKLConfig.cmake resolved to the
      # Intel interface layer instead, fall back to the explicit list below
      # rather than producing a subtly mis-linked binary.
      get_target_property(_ff_mkl_iface MKL::MKL INTERFACE_LINK_LIBRARIES)
      if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU"
          AND NOT "${_ff_mkl_iface}" MATCHES "mkl_gf")
        message(STATUS
          "MKLConfig.cmake did not select the GNU Fortran interface layer; "
          "falling back to the explicit MKL link group")
      else()
        set(FF_MKL_PARDISO_LIBS MKL::MKL)
        message(STATUS
          "MKL PARDISO enabled via MKLConfig.cmake (${_ff_mkl_config_dir}), "
          "link=${MKL_LINK} interface=${MKL_INTERFACE} threading=${MKL_THREADING}")
        message(WARNING
          "FF_MKL_USE_CMAKE_CONFIG=ON moves MKL to the end of the link line. "
          "On toolchains that mix a static libgfortran.a with the system "
          "libgfortran.so.5 this can split the Fortran runtime and hang the "
          "solver at startup in the parameter parser. Verify that "
          "'nm <app> | grep _gfortran_st_close' reports T (not U).")
      endif()
      unset(_ff_mkl_iface)
    endif()
  endif()

  # Default: the hand-assembled static link group, emitted literally at the
  # FF_MKL_PARDISO_LIBS position in FF_DEFAULT_LIBS / FF_APPLICATION_LIBS.
  if(NOT FF_MKL_PARDISO_LIBS)
    foreach(_mkl_lib mkl_gf_lp64 mkl_sequential mkl_core)
      if(NOT EXISTS ${MKL_PARDISO_DIR}/lib/lib${_mkl_lib}.a)
        message(FATAL_ERROR
          "MKL library not found: ${MKL_PARDISO_DIR}/lib/lib${_mkl_lib}.a")
      endif()
    endforeach()
    set(FF_MKL_PARDISO_LIBS
      -Wl,--start-group
      ${MKL_PARDISO_DIR}/lib/libmkl_gf_lp64.a
      ${MKL_PARDISO_DIR}/lib/libmkl_sequential.a
      ${MKL_PARDISO_DIR}/lib/libmkl_core.a
      -Wl,--end-group
      pthread m dl
      )
    message(STATUS
      "MKL PARDISO enabled from ${MKL_PARDISO_DIR} (static link group)")
  endif()

  add_definitions(-DMKL_PARDISO_AVAIL)
  unset(_ff_mkl_config_dir)

endmacro()
