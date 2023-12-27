cmake_minimum_required(VERSION 3.1...3.23)

# name of the project, specIFy the used languages
project(Q2P1 C CXX Fortran)

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(FF_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# allow less strict IF-else syntax
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# output the system name
message(STATUS "Configuring Q2P1 for a ${CMAKE_SYSTEM} system")

if(UNIX)
  find_program(LSB_RELEASE lsb_release)
  if(LSB_RELEASE)
    execute_process(COMMAND ${LSB_RELEASE} -is
      OUTPUT_VARIABLE LSB_RELEASE_ID_SHORT
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    message(STATUS "Linux flavor: ${LSB_RELEASE_ID_SHORT}")
  endif()

  find_program(MODULE_COMMAND modulecmd PATHS $ENV{MODULESHOME}/bin)
endif()

enable_testing()

include(CTest)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

#=========================================================================
#                         custom targets  
#=========================================================================
include(${CMAKE_MODULE_PATH}/CustomTargets.cmake)  

#=========================================================================
#                         compiler settings
#=========================================================================
message(STATUS "Checking build type...")

IF(Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER ${Q2P1_COMPILER})
ELSE()
  set(Q2P1_DEFAULT_COMPILER "gcc")
ENDIF()

include(${CMAKE_MODULE_PATH}/GenerateBuildId.cmake)  

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

option(BUILD_TESTS
  "Enables build of tests"
  OFF
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

if(BUILD_BOUNDARY_LAYER_TOOLS)
  SET(USE_OPENMESH ON CACHE BOOL "Build OpenMesh library" FORCE) 
endif(BUILD_BOUNDARY_LAYER_TOOLS)

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

#=========================================================================
#                       CMake BUILD Type
#=========================================================================
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
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
  add_subdirectory(extern/libraries/metis-4.0.3)
endif(BUILD_METIS)

#=========================================================================
# Add the cgal library
#=========================================================================
message(STATUS "Running find_package CGAL.")
# Add CGAL submodule path
if(USE_CGAL)
  set(CGAL_DIR ${CMAKE_SOURCE_DIR}/extern/libraries/cgal-install-dir/lib/cmake/CGAL)
  message(STATUS "CGAL dir: ${CGAL_DIR}")
  find_package(CGAL REQUIRED PATHS "${CMAKE_SOURCE_DIR}/extern/libraries/cgal-install-dir/lib/cmake/CGAL")
endif()

if(CGAL_FOUND)
    message(STATUS "CGAL found.")
    message(STATUS "CGAL version: ${CGAL_VERSION}")
    message(STATUS "CGAL include directory: ${CGAL_INCLUDE_DIRS}")
    set(CGAL_INCLUDE_DIR "${CGAL_INCLUDE_DIRS}")
    # Additional messages or actions if needed
else()
    message(FATAL_ERROR "CGAL not found. Please check your CGAL_DIR setting.")
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
#=========================================================================

# Check if USE_CGAL is ON, and set USE_BOOST accordingly
if(USE_CGAL)
    set(USE_BOOST ON)
endif()

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

      # Set Boost_USE_STATIC_LIBS to ON
      set(Boost_USE_STATIC_LIBS ON)
      set(Boost_USE_MULTITHREADED ON )
  else()
      message(FATAL_ERROR "Boost ${Boost_MIN_VERSION} not found.")
  endif()

endif(USE_BOOST)

#=========================================================================
#                     Check for the 'librt' library
#=========================================================================
if (UNIX AND NOT APPLE)
    find_library(LIBRT_LIBRARY NAMES rt)
    if(NOT LIBRT_LIBRARY)
        find_file(LIBRT_LIBRARY librt.so  PATHS /usr/lib /usr/lib64/)
        if(NOT LIBRT_LIBRARY)
            message(FATAL_ERROR "Please install the system library librt to compile this code.")
        endif()
    endif()
elseif (APPLE)
    # Handle macOS case if needed
endif()

#=========================================================================
#                      Configure MPI/MKL/MUMPS 
#=========================================================================
find_package(MPI)

SET(MPI_VENDOR "openmpi" CACHE STRING "Variant of the OpenMPI implementation")

message(STATUS "CMake MPI_FORTRAN_INCLUDE_PATH : ${MPI_Fortran_INCLUDE_PATH}")
message(STATUS "CMake MPI_Fortran_LIBRARIES : ${MPI_Fortran_LIBRARIES}")

add_subdirectory(extern/libraries/feat2d)
add_subdirectory(extern/libraries/feat3d)
add_subdirectory(extern/libraries/amd)
add_subdirectory(extern/libraries/umfpack4)

if(USE_MUMPS)
  find_package(MKL)
  if(${Q2P1_BUILD_ID} MATCHES "intel")
    add_definitions(-DMUMPS_AVAIL)
    add_subdirectory(extern/libraries/MUMPS_5.0.1)
    set(src_q2p1 ${src_q2p1} ${CMAKE_SOURCE_DIR}/source/MumpsSolver.f90)
  else(${Q2P1_BUILD_ID} MATCHES "intel")
    message(FATAL_ERROR "The MUMPS library requires the intel compiler. Choose a build-id with the intel compiler to resolve this issue.")
  endif(${Q2P1_BUILD_ID} MATCHES "intel")
endif(USE_MUMPS)

find_package(Git 2.17)

if(NOT WIN32)
  if(USE_CGAL)
    add_subdirectory(tools/check_manifold)
  endif(USE_CGAL)
endif(NOT WIN32)


#=========================================================================
# Add the hypre library
#=========================================================================
if(USE_HYPRE)
  # Use the sources from the repository  
  message(STATUS "Using repository versions of hypre")
  add_definitions(-DHYPRE_AVAIL)		
  add_subdirectory(extern/libraries/hypre/src)
  set(src_q2p1 ${src_q2p1} ${CMAKE_SOURCE_DIR}/source/HypreSolver.f90)
  set(HYPRE_LIBRARIES HYPRE)
endif(USE_HYPRE)

#=========================================================================
#      Set the flags/paths/etc. we have determined for libraries 
#=========================================================================
message(STATUS "Module path was set to ${CMAKE_MODULE_PATH}")
include(${CMAKE_SOURCE_DIR}/cmake/modules/GenerateIncludeFlags.cmake)  
include(${CMAKE_SOURCE_DIR}/cmake/modules/GenerateLinkerFlags.cmake)  

#=========================================================================
#                         project directories
#=========================================================================
include(${CMAKE_SOURCE_DIR}/cmake/modules/ProjectFiles.cmake)  

#=========================================================================
#                         Add applications
#=========================================================================
SET(OUT_OF_SOURCE_BUILD True)
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
  SET(${OUT_OF_SOURCE_BUILD} False)
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# we need additional libraries to support CUDA
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

  add_subdirectory(applications)
ENDIF(BUILD_APPLICATIONS)

#=========================================================================
#                     PE Library
#=========================================================================
add_subdirectory(libs/pe)

#=========================================================================
#                    TESTS 
#=========================================================================
if(BUILD_TESTS)
  add_subdirectory(tests/cgal_interface)
endif()


#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(DIRECTORY tools/partitioner DESTINATION bin)
file(GLOB PYPARTITIONER "tools/PyPartitioner.py")
install(FILES ${PYPARTITIONER} DESTINATION bin
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)
#install(FILES ${E3D_SCRIPTS} DESTINATION bin
#    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)

#=========================================================================
#                       CPACK CONFIGURATION
#=========================================================================
include(${FF_MODULE_PATH}/GeneratePackageDeb.cmake)  

include(CPack)

#git@github.com:rmuenste/cgal.git
