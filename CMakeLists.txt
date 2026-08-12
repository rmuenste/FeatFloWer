# 3.21 is the floor imposed by the vendored/fetched hypre 2.33.0 sub-project.
# (FetchContent's FIND_PACKAGE_ARGS, used by the fetch providers, needs 3.24 --
# see FF_FETCHCONTENT_MIN_VERSION in cmake/modules/FFDependencies.cmake.)
cmake_minimum_required(VERSION 3.21)


# name of the project, specIFy the used languages
project(Q2P1 C CXX Fortran)

# Here we can set a project version
set(PROJECT_VERSION "1.0.0")

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(FF_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# allow less strict IF-else syntax
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# output the system name
message(STATUS "Configuring Q2P1 for a ${CMAKE_SYSTEM} system")

if(UNIX)
  find_program(LSB_RELEASE lsb_release)
  execute_process(COMMAND ${LSB_RELEASE} -is
    OUTPUT_VARIABLE LSB_RELEASE_ID_SHORT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  message(STATUS "Linux flavor: ${LSB_RELEASE_ID_SHORT}")

  find_program(MODULE_COMMAND modulecmd PATHS $ENV{MODULESHOME}/bin)
endif(UNIX)

enable_testing()

include(CTest)
include(CMakeDependentOption)
include(FetchContent)

#=========================================================================
#                         custom targets  
#=========================================================================
include(CustomTargets)  
include(GenerateBuildIdOutput)  
option(SHOW_BUILD_IDS
  "Display all available build ids and exit"
  OFF
  )

if(SHOW_BUILD_IDS)
  generate_string_list()
  message(FATAL_ERROR "Ending script, this is as intended.")
endif(SHOW_BUILD_IDS)


#=========================================================================
#                         compiler settings
#=========================================================================
message(STATUS "Checking build type...")
include(CompilerQuickCheck)

IF(NOT Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER "gcc")
ELSE(NOT Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER ${Q2P1_COMPILER})
ENDIF(NOT Q2P1_COMPILER)

include(GenerateBuildId)  

#=========================================================================
#                             Options
#=========================================================================

option(USE_OPENMESH
  "Use the OpenMesh library"
  OFF
  )

option(USE_OPTICALTWEEZERS
  "Use the library for simulating optical forces"
  OFF
  )

option(USE_MUMPS
  "Use the MUMPS MPI parallel direct coarse grid solver"
  OFF
  )

option(USE_ODE
  "Use the ODE library"
  OFF
  )

option(USE_CGAL
  "Use CGAL library"
  OFF
  )

option(USE_CGAL_LOCAL
  "This option will try to find a local installation of the CGAL library, instead of cloning the Git repository."
  OFF
  )

# USE_CGAL_LOCAL implies USE_CGAL — no need to set both on the command line.
if(USE_CGAL_LOCAL AND NOT USE_CGAL)
  set(USE_CGAL ON CACHE BOOL "Use CGAL library" FORCE)
  message(STATUS "USE_CGAL_LOCAL=ON implies USE_CGAL=ON (auto-enabled)")
endif()

option(USE_BOOST
  "Use Boost library"
  OFF
  )

option(USE_SYSTEM_BLASLAPACK
  "Use the system BLAS/LAPACK libraries if CMAKE can find them."
  OFF
  )

option(BUILD_METIS
  "Enables build of metis library"
  ON
  )

option(BUILD_APPLICATIONS
  "Enables build of applications"
  ON
  )

option(Q1P1_SSE_LINK_STATIC_FORTRAN
  "Enables static linking of the Fortran runtime library"
  OFF 
  )

option(BUILD_BOUNDARY_LAYER_TOOLS
  "Build the boundary layer tools"
  OFF
  )

option(USE_HYPRE
  "Use hypre library"
  OFF
  )

cmake_dependent_option(
  USE_EXTERNAL_HYPRE
  "Link a pre-installed hypre instead of the vendored copy"
  OFF
  "USE_HYPRE"
  OFF
  )

set(EXTERNAL_HYPRE_DIR "" CACHE PATH
  "Install prefix of an external hypre"
  )

cmake_dependent_option(
  USE_HYPRE_CUDA
  "Enable hypre's CUDA execution path"
  OFF
  "USE_HYPRE"
  OFF
  )

cmake_dependent_option(
  USE_HYPRE_GPU_AMG
  "Use GPU-compatible BoomerAMG parameters"
  OFF
  "USE_HYPRE"
  OFF
  )

cmake_dependent_option(
  ENABLE_HYPRE_TIMING
  "Print machine-readable hypre coarse-solver timings"
  OFF
  "USE_HYPRE"
  OFF
  )

option(ENABLE_SOLVER_TIMING
  "Print machine-readable FF_TIMING records for MG stages and coarse solvers"
  OFF
  )

option(USE_EXTERNAL_SUITESPARSE
  "Use an external SuiteSparse UMFPACK instead of the vendored UMFPACK 4"
  OFF
  )

set(EXTERNAL_SUITESPARSE_DIR "" CACHE PATH
  "Install prefix of an external SuiteSparse (empty: system paths)"
  )

option(USE_MKL_PARDISO
  "Enable MKL PARDISO as pressure coarse solver type 9"
  OFF
  )

set(MKL_PARDISO_DIR "$ENV{MKLROOT}" CACHE PATH
  "MKL install root providing lib/libmkl_*.a (required for USE_MKL_PARDISO)"
  )

cmake_dependent_option(
  USE_EXTERNAL_MUMPS
  "Link a pre-installed MUMPS (works with GCC) instead of the vendored MUMPS 5.0.1 + MKL"
  OFF
  "USE_MUMPS"
  OFF
  )

set(EXTERNAL_MUMPS_DIR "" CACHE PATH
  "Install prefix of an external MUMPS (lib/, include/)"
  )

set(EXTERNAL_SCALAPACK_DIR "" CACHE PATH
  "Install prefix of the ScaLAPACK matching the external MUMPS"
  )

set(EXTERNAL_SCOTCH_DIR "" CACHE PATH
  "Install prefix of the SCOTCH/esmumps matching the external MUMPS"
  )

set(EXTERNAL_PARMETIS_DIR "" CACHE PATH
  "Install prefix of the ParMETIS/METIS matching the external MUMPS"
  )

option(USE_PE
  "Use pe library"
  OFF
  )

cmake_dependent_option(
  USE_PE_SERIAL_MODE
  "Use serial PE mode for large particles (bypasses PE MPI)"
  OFF              # Default: use parallel PE
  "USE_PE"         # Only show this option if USE_PE=ON
  OFF              # Value when USE_PE=OFF
)

cmake_dependent_option(
  PE_USE_MPI                       # Dependent Option name (PE-owned, PE_-prefixed)
  "Enable MPI support in PE"       # Description
  ON                               # Value when condition is true (parallel PE needs MPI)
  "USE_PE;NOT USE_PE_SERIAL_MODE"  # Condition: USE_PE=ON AND USE_PE_SERIAL_MODE=OFF
  OFF                              # Value when condition is false (no PE or serial PE)
)

option(SED_BENCH
  "Enable sedimentation benchmark output (force/position/velocity per timestep)"
  OFF
)

option(ENABLE_FBM_ACCELERATION
  "Enable FBM acceleration techniques (HashGrid alpha + KVEL force). Requires USE_PE."
  OFF
)

if(NOT USE_PE AND ENABLE_FBM_ACCELERATION)
  message(STATUS "ENABLE_FBM_ACCELERATION requested with USE_PE=OFF; forcing ENABLE_FBM_ACCELERATION=OFF.")
  set(ENABLE_FBM_ACCELERATION OFF CACHE BOOL
    "Enable FBM acceleration techniques (HashGrid alpha + KVEL force). Requires USE_PE."
    FORCE
  )
endif()


if(BUILD_BOUNDARY_LAYER_TOOLS)
  SET(USE_OPENMESH ON CACHE BOOL "Build OpenMesh library" FORCE)
endif(BUILD_BOUNDARY_LAYER_TOOLS)

#=========================================================================
#      Solver dependency providers (vendored / external / system / fetch)
#      Declares FF_FETCH_DEPENDENCIES, FF_HYPRE_PROVIDER,
#      FF_SUITESPARSE_PROVIDER, FF_HYPRE_GPU_AWARE_MPI and the
#      ff_setup_* macros used further down.
#      See docs/md_docs/fetchcontent_dependencies.md
#=========================================================================
include(FFDependencies)

#option(FC_CUDA_SUPPORT
#    "Enables use of CUDA extensions"
#    OFF
#)
#
## we need additional libraries to support CUDA
#IF(FC_CUDA_SUPPORT)
#
#  add_definitions(-DFC_CUDA_SUPPORT)
#  message(STATUS "CUDA define set")  
#
#  # check for openGL
#  FIND_PACKAGE(OpenGL REQUIRED)
#  
#  # check for the CUDA libraries
#  FIND_PACKAGE( CUDA )
#  
#  IF(CUDA_FOUND)
#    message(STATUS "Found CUDA")
#  ELSEIF(CUDA_FOUND)
#    message(FATAL_ERROR "The CUDA library was not found on the system.")
#  ENDIF(CUDA_FOUND)
#  
#  # search for cuda libraries
#  FIND_PACKAGE(CUDASDK)
#ENDIF(FC_CUDA_SUPPORT)

find_package(OpenMP REQUIRED)
#=========================================================================
#                       CMake BUILD Type
#=========================================================================
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

#=========================================================================
#                     Check for the 'librt' library
#=========================================================================
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
  FIND_LIBRARY(LIBRT_LIBRARY NAMES rt)
  if(NOT LIBRT_LIBRARY)
    find_file(LIBRT_LIBRARY librt.so  PATHS /usr/lib /usr/lib64/)
    if(NOT LIBRT_LIBRARY)
      message(FATAL_ERROR "You need to install the system library librt.")
    endif()
  endif()
endif()

#=========================================================================
#                       Configure BLAS/LAPACK 
#=========================================================================
if(USE_SYSTEM_BLASLAPACK)

  # Try to find BLAS/LAPACK on the system
  message(STATUS "Trying to find BLAS/LAPACK on the system.")
  find_package(BLAS)
  find_package(LAPACK)

  # Use system libraries if they are found
  if(BLAS_FOUND AND LAPACK_FOUND)
    message(STATUS "Found BLAS/LAPACK on the system.")
    message(STATUS "Blas library: ${BLAS_LIBRARIES}")
    message(STATUS "Lapack library: ${LAPACK_LIBRARIES}")

    # Use the sources from the repository if 
    # no libraries are found on the system
  else(BLAS_FOUND AND LAPACK_FOUND)
    message(STATUS "BLAS/LAPACK not found on the system, using repository version.")
    add_subdirectory(extern/libraries/lapack-3.6.1)
    set(BLAS_LIBRARIES blas)
    set(LAPACK_LIBRARIES lapack)
  endif(BLAS_FOUND AND LAPACK_FOUND)

else(USE_SYSTEM_BLASLAPACK)
  # Use the sources from the repository  
  message(STATUS "Using repository versions of BLAS/LAPACK")
  add_subdirectory(extern/libraries/lapack-3.6.1)
  set(BLAS_LIBRARIES blas)
  set(LAPACK_LIBRARIES lapack)
endif(USE_SYSTEM_BLASLAPACK)

if(BUILD_METIS)
  add_subdirectory(extern/libraries/metis-5.1.0)
endif(BUILD_METIS)

#=========================================================================
# Add the cgal library
#=========================================================================
if(USE_CGAL)
  set(USE_BOOST TRUE)

  # CGAL's imported target must carry these exact arithmetic dependencies to
  # every C++ consumer, including PE and FullC0ntact applications.
  find_package(GMP REQUIRED)
  find_package(MPFR REQUIRED)
  message(STATUS "Found GMP: ${GMP_LIBRARIES}")
  message(STATUS "GMP include directory: ${GMP_INCLUDE_DIR}")
  message(STATUS "Found MPFR: ${MPFR_LIBRARIES}")
  message(STATUS "MPFR include directory: ${MPFR_INCLUDE_DIR}")

  # Prevent CGAL's TweakFindBoost.cmake from silently requesting static Boost.
  if(NOT DEFINED CGAL_Boost_USE_STATIC_LIBS)
    set(CGAL_Boost_USE_STATIC_LIBS OFF CACHE BOOL
      "Tell CGAL's FindBoost tweak to prefer shared Boost libraries" FORCE)
    message(STATUS "CGAL_Boost_USE_STATIC_LIBS not set -- defaulting to OFF (shared Boost)")
  else()
    message(STATUS "CGAL_Boost_USE_STATIC_LIBS explicitly set to ${CGAL_Boost_USE_STATIC_LIBS}")
  endif()

  if(USE_CGAL_LOCAL)
    message(STATUS "Trying to find system installation of CGAL")
    find_package(CGAL 5.3 REQUIRED CONFIG)
  else()
    message(STATUS "Fetching CGAL 5.6.3")
    set(CGAL_BUILD_DOCS OFF CACHE BOOL "Build CGAL documentation" FORCE)
    set(CGAL_BUILD_EXAMPLES OFF CACHE BOOL "Build CGAL examples" FORCE)
    set(CGAL_BUILD_TESTING OFF CACHE BOOL "Build CGAL tests" FORCE)
    set(WITH_CGAL_ImageIO OFF CACHE BOOL "Build CGAL ImageIO support" FORCE)
    set(WITH_CGAL_Qt5 OFF CACHE BOOL "Build CGAL Qt5 support" FORCE)

    FetchContent_Declare(
      CGAL
      GIT_REPOSITORY https://github.com/CGAL/cgal.git
      GIT_TAG v5.6.3
      GIT_SHALLOW TRUE
    )
    set(_ff_build_testing "${BUILD_TESTING}")
    set(BUILD_TESTING OFF)
    FetchContent_MakeAvailable(CGAL)
    set(BUILD_TESTING "${_ff_build_testing}")
    unset(_ff_build_testing)
    find_package(CGAL REQUIRED CONFIG
      PATHS "${cgal_BINARY_DIR}"
      NO_DEFAULT_PATH
    )
  endif()

  if(NOT TARGET CGAL::CGAL)
    message(FATAL_ERROR "CGAL support is enabled, but CGAL::CGAL is unavailable")
  endif()

  message(STATUS "Using CGAL version ${CGAL_VERSION}")
  set(CGAL_LIBRARIES CGAL::CGAL)
  if(CGAL_INCLUDE_DIRS)
    set(CGAL_INCLUDE_DIR ${CGAL_INCLUDE_DIRS})
  else()
    get_target_property(CGAL_INCLUDE_DIR CGAL::CGAL INTERFACE_INCLUDE_DIRECTORIES)
  endif()
endif()

#=========================================================================
#                          Configure BoostC++ 
#=========================================================================
# If there are issues on WINDOWS finding the BoostC++ libraries, try to set these variables
# that help CMake find the correct Boost location.
#   Boost_DEBUG              - Set to ON to enable debug output from FindBoost.
#                              Please enable this before filing any bug report.
#   Boost_USE_MULTITHREADED  - Set to OFF to use the non-multithreaded
#                              libraries ('mt' tag).  Default is ON.
#   Boost_USE_STATIC_LIBS    - Set to ON to force the use of the static
#                              libraries.  Default is OFF.
#   Boost_USE_STATIC_RUNTIME - Set to ON or OFF to specify whether to use
#                              libraries linked statically to the C++ runtime
#                              ('s' tag).  Default is platform dependent.
#   Boost_COMPILER           - Set to the compiler-specific library suffix
#
#   Boost_USE_DEBUG_RUNTIME  - Set to ON or OFF to specify whether to use
#                              libraries linked to the MS debug C++ runtime
#                              ('g' tag).  Default is ON.

# find_package(Boost 1.56 REQUIRED COMPONENTS thread): 
# find_package(Boost [version] [EXACT] [REQUIRED] [COMPONENTS <libs> ...]) 
# Here [version] denotes the minimum version of the library
# In case the libraries are found the following variables are set:
# 
# Boost_INCLUDE_DIRS : - Boost include directories
# Boost_LIBRARY_DIRS : - Boost library directories
# Boost_LIBRARIES    : - Boost libraries to be linked
# Boost version vars : - Boost detailed library version: ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}

if(USE_BOOST)

  set(Boost_MIN_VERSION "1.68.0")

  # Find Boost with the specified components
  find_package(Boost ${Boost_MIN_VERSION} REQUIRED COMPONENTS
      thread
      system
      filesystem
      program_options
      random
  )

  # Check if Boost was found
  if (Boost_FOUND)
      message(STATUS "Boost found: ${Boost_VERSION}")
      message(STATUS "Boost include directories: ${Boost_INCLUDE_DIRS}")
      message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
  else()
      message(FATAL_ERROR "Boost ${Boost_MIN_VERSION} not found.")
  endif()

endif(USE_BOOST)

#=========================================================================
#                      Configure MPI/MKL/MUMPS 
#=========================================================================
find_package(MPI)

SET(MPI_VENDOR "openmpi" CACHE STRING "Variant of the OpenMPI implementation")

message(STATUS "CMake MPI_FORTRAN_INCLUDE_PATH : ${MPI_Fortran_INCLUDE_PATH}")
message(STATUS "CMake MPI_Fortran_LIBRARIES : ${MPI_Fortran_LIBRARIES}")

add_subdirectory(extern/libraries/feat2d)
add_subdirectory(extern/libraries/feat3d)

# SuiteSparse/UMFPACK provider: vendored UMFPACK 4, a site install
# (USE_EXTERNAL_SUITESPARSE + EXTERNAL_SUITESPARSE_DIR) or the pinned
# SuiteSparse tarball (FF_SUITESPARSE_PROVIDER=fetch).
ff_setup_suitesparse()

if(USE_MKL_PARDISO)
  ff_setup_mkl_pardiso()
endif()

if(USE_MUMPS)
  if(USE_EXTERNAL_MUMPS)
    foreach(_mumps_lib dmumps mumps_common pord)
      string(TOUPPER ${_mumps_lib} _MUMPS_LIB_U)
      find_library(FF_${_MUMPS_LIB_U}_LIBRARY ${_mumps_lib}
        HINTS ${EXTERNAL_MUMPS_DIR}/lib ${EXTERNAL_MUMPS_DIR}/lib64)
      if(NOT FF_${_MUMPS_LIB_U}_LIBRARY)
        message(FATAL_ERROR
          "External MUMPS: lib${_mumps_lib} not found under ${EXTERNAL_MUMPS_DIR}")
      endif()
    endforeach()
    if(NOT EXISTS ${EXTERNAL_MUMPS_DIR}/include/dmumps_struc.h)
      message(FATAL_ERROR
        "External MUMPS: dmumps_struc.h not found under ${EXTERNAL_MUMPS_DIR}/include")
    endif()
    find_library(FF_SCALAPACK_LIBRARY scalapack
      HINTS ${EXTERNAL_SCALAPACK_DIR}/lib ${EXTERNAL_SCALAPACK_DIR}/lib64)
    if(NOT FF_SCALAPACK_LIBRARY)
      message(FATAL_ERROR
        "External MUMPS: libscalapack not found under ${EXTERNAL_SCALAPACK_DIR}")
    endif()
    # SCOTCH orderings are optional in principle but required by the site
    # MUMPS builds; collect whatever variants the prefix provides.
    set(FF_SCOTCH_LIBRARIES "")
    foreach(_scotch_lib esmumps ptscotch scotch ptscotcherr scotcherr)
      find_library(FF_SCOTCH_${_scotch_lib}_LIBRARY ${_scotch_lib}
        HINTS ${EXTERNAL_SCOTCH_DIR}/lib ${EXTERNAL_SCOTCH_DIR}/lib64)
      if(FF_SCOTCH_${_scotch_lib}_LIBRARY)
        list(APPEND FF_SCOTCH_LIBRARIES ${FF_SCOTCH_${_scotch_lib}_LIBRARY})
      endif()
    endforeach()
    # The site MUMPS builds also use ParMETIS orderings.
    set(FF_PARMETIS_LIBRARIES "")
    foreach(_pm_lib parmetis metis gklib)
      find_library(FF_PARMETIS_${_pm_lib}_LIBRARY ${_pm_lib}
        HINTS ${EXTERNAL_PARMETIS_DIR}/lib ${EXTERNAL_PARMETIS_DIR}/lib64
        NO_DEFAULT_PATH)
      if(FF_PARMETIS_${_pm_lib}_LIBRARY)
        list(APPEND FF_PARMETIS_LIBRARIES ${FF_PARMETIS_${_pm_lib}_LIBRARY})
      endif()
    endforeach()
    if(NOT FF_PARMETIS_LIBRARIES)
      message(FATAL_ERROR
        "External MUMPS: libparmetis not found under ${EXTERNAL_PARMETIS_DIR}")
    endif()
    find_library(FF_OPENBLAS_LIBRARY openblas)
    if(NOT FF_OPENBLAS_LIBRARY)
      message(FATAL_ERROR "External MUMPS: libopenblas not found")
    endif()
    add_definitions(-DMUMPS_AVAIL)
    set(src_q2p1 ${src_q2p1} ${CMAKE_SOURCE_DIR}/source/MumpsSolver.f90)
    set(MUMPS_STRUCTURES_INCLUDE_PATH ${EXTERNAL_MUMPS_DIR}/include)
    message(STATUS "External MUMPS: ${FF_DMUMPS_LIBRARY}")
    message(STATUS "External MUMPS ScaLAPACK: ${FF_SCALAPACK_LIBRARY}")
    message(STATUS "External MUMPS SCOTCH libs: ${FF_SCOTCH_LIBRARIES}")
  else()
    find_package(MKL)
    if(${Q2P1_BUILD_ID} MATCHES "intel")
      add_definitions(-DMUMPS_AVAIL)
      add_subdirectory(extern/libraries/MUMPS_5.0.1)
      set(src_q2p1 ${src_q2p1} ${CMAKE_SOURCE_DIR}/source/MumpsSolver.f90)
    else(${Q2P1_BUILD_ID} MATCHES "intel")
      message(FATAL_ERROR "The vendored MUMPS library requires the intel compiler. Choose a build-id with the intel compiler, or set USE_EXTERNAL_MUMPS=ON with a GCC-compatible MUMPS install.")
    endif(${Q2P1_BUILD_ID} MATCHES "intel")
  endif()
endif(USE_MUMPS)

find_package(Git 2.17)

#=========================================================================
# This builds our small CGAL mesh validation tool
#=========================================================================
if(NOT WIN32)
  if(USE_CGAL)
    add_subdirectory(tools/check_manifold)
  endif(USE_CGAL)
endif(NOT WIN32)


#=========================================================================
# Add the hypre library
#=========================================================================
if(USE_HYPRE)
  add_definitions(-DHYPRE_AVAIL)

  # hypre provider: vendored submodule, a site install
  # (USE_EXTERNAL_HYPRE + EXTERNAL_HYPRE_DIR) or the pinned hypre tarball
  # (FF_HYPRE_PROVIDER=fetch).
  ff_setup_hypre()

  if(USE_HYPRE_CUDA)
    add_definitions(-DHYPRE_CUDA)
    if(NOT USE_HYPRE_GPU_AMG)
      message(STATUS "USE_HYPRE_CUDA enables USE_HYPRE_GPU_AMG")
      set(USE_HYPRE_GPU_AMG ON CACHE BOOL
        "Use GPU-compatible BoomerAMG parameters" FORCE)
    endif()
  endif()

  if(USE_HYPRE_GPU_AMG)
    add_definitions(-DHYPRE_GPU_AMG)
  endif()

  if(ENABLE_HYPRE_TIMING)
    add_definitions(-DHYPRE_TIMING)
  endif()

  set(src_q2p1 ${src_q2p1} ${CMAKE_SOURCE_DIR}/source/HypreSolver.f90)
endif(USE_HYPRE)

if(ENABLE_SOLVER_TIMING)
  add_definitions(-DFF_SOLVER_TIMING)
endif()

#=========================================================================
#      Set the flags/paths/etc. we have determined for libraries 
#=========================================================================
message(STATUS "Module path was set to ${CMAKE_MODULE_PATH}")
include(GenerateIncludeFlags)  
include(GenerateLinkerFlags)  

#=========================================================================
#                         project directories
#=========================================================================
include(ProjectFiles)  

#=========================================================================
#                     PE Library
#=========================================================================
message(STATUS "Adding pe library to build")
if(USE_PE)
  add_definitions(-DHAVE_PE)

  # Additional flag for serial mode
  if(USE_PE_SERIAL_MODE)
    add_definitions(-DPE_SERIAL_MODE)
    message(STATUS "==============================================")
    message(STATUS "PE SERIAL MODE ENABLED")
    message(STATUS "  - Large particle optimization")
    message(STATUS "  - PE runs serially in each CFD domain")
    message(STATUS "  - Forces synchronized via CFD MPI")
    message(STATUS "  - Bypasses PE distributed shadow copies")
    message(STATUS "==============================================")
  else()
    message(STATUS "PE Parallel Mode (default): Full MPI parallelization")
  endif()

  # Propagate FBM acceleration to PE library's HashGrid
  if(ENABLE_FBM_ACCELERATION)
    set(PE_USE_ACCELERATED_POINT_QUERY ON CACHE BOOL
      "Enable HashGrid spatial hashing in PE (controlled by ENABLE_FBM_ACCELERATION)" FORCE)
  endif()

  # Forward FeatFloWer's CGAL choice to the PE library. PE removed its
  # USE_CGAL deprecated alias; its own option is now PE_USE_CGAL (default OFF),
  # so we must set it explicitly to keep PE's CGAL code paths compiled.
  if(USE_CGAL)
    set(PE_USE_CGAL ON CACHE BOOL
      "Build PE with CGAL support (controlled by FeatFloWer USE_CGAL)" FORCE)
  endif()

  add_subdirectory(libs/pe)
endif(USE_PE)

if(SED_BENCH)
  add_definitions(-DSED_BENCH)
  message(STATUS "SED_BENCH: sedimentation benchmark output enabled")
endif()

if(ENABLE_FBM_ACCELERATION)
  add_definitions(-DENABLE_FBM_ACCELERATION)
  message(STATUS "FBM acceleration enabled (HashGrid + KVEL). Control at runtime via q2p1_param.dat")
else()
  message(STATUS "==============================================")
  message(STATUS "FBM ACCELERATION DISABLED")
  message(STATUS "  - Pure baseline (brute-force) build")
  message(STATUS "  - For correctness verification only")
  message(STATUS "==============================================")
endif()

option(DEBUG_FBM_OPTIMIZATION
  "Enable FBM optimization debug mode: run both standard and KVEL force computation, compare results, write per-timestep debug files. Requires ENABLE_FBM_ACCELERATION."
  OFF
)

if(NOT ENABLE_FBM_ACCELERATION AND DEBUG_FBM_OPTIMIZATION)
  message(STATUS "DEBUG_FBM_OPTIMIZATION requested with ENABLE_FBM_ACCELERATION=OFF; forcing OFF.")
  set(DEBUG_FBM_OPTIMIZATION OFF CACHE BOOL
    "Enable FBM optimization debug mode" FORCE)
endif()

if(DEBUG_FBM_OPTIMIZATION)
  add_definitions(-DDEBUG_FBM_OPTIMIZATION)
  message(STATUS "DEBUG_FBM_OPTIMIZATION: dual force computation + per-timestep debug output enabled")
endif()

#=========================================================================
#                         Add applications
#=========================================================================
SET(OUT_OF_SOURCE_BUILD True)
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
  SET(${OUT_OF_SOURCE_BUILD} False)
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

IF(BUILD_APPLICATIONS)
  # add subdirectories
  IF(EXISTS "${CMAKE_SOURCE_DIR}/FullC0ntact")
    add_definitions(-DFEATFLOWLIB)
    message(STATUS "Found FullC0ntact library directory")  
    add_subdirectory(FullC0ntact)
    if(USE_OPENMESH)
      add_subdirectory(FullC0ntact/libs/OpenMesh)
    endif(USE_OPENMESH)
    if(USE_ODE)
      add_definitions(-DWITH_ODE)
      add_subdirectory(FullC0ntact/libs/ode-cmake)
    endif(USE_ODE)
  else()
    message(FATAL_ERROR "FullC0ntact directory is required for applications.")  
  ENDIF()

  # Include CreateDataDirectories.cmake to make the function available to all applications
  include(CreateDataDirectories)
  
  add_subdirectory(applications)
ENDIF(BUILD_APPLICATIONS)

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(DIRECTORY tools/partitioner DESTINATION bin)
install(FILES tools/PyPartitioner.py DESTINATION bin PERMISSIONS 
  WORLD_EXECUTE OWNER_READ OWNER_WRITE 
  OWNER_EXECUTE GROUP_READ GROUP_WRITE 
  GROUP_EXECUTE)

#=========================================================================
#                       CPACK CONFIGURATION
#=========================================================================
include(GeneratePackageDeb)  

include(CPack)
