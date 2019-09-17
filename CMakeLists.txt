cmake_minimum_required(VERSION 3.0)

# name of the project, specIFy the used languages
project(Q2P1 C CXX Fortran)

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

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

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

#=========================================================================
#                         custom targets  
#=========================================================================
include(${CMAKE_MODULE_PATH}/CustomTargets.cmake)  

#=========================================================================
#                         compiler settings
#=========================================================================
message(STATUS "Checking build type...")

IF(NOT Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER "gcc")
ELSE(NOT Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER ${Q2P1_COMPILER})
ENDIF(NOT Q2P1_COMPILER)

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
  OFF
  )

option(BUILD_APPLICATIONS
  "Enables build of applications"
  ON
  )

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
#                         project directories
#=========================================================================
set(src_pp3d
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/projma.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/coeff.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/indat3d.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/gupwd.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/nsdef.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/bndry.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/conv.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/util.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/avs3d.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/period.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/parq3d.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/xmrout.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/rdparm.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/trsort.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/orsc.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/mgrout.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/inout.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/user.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/optcnl.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/dfkt.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/error.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/bmul.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/diff.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/nu_turb.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/prostp.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/def_feat.f90
  #source/src_pproc/pproc_main.f90
  )

set(postprocessing
  ${CMAKE_SOURCE_DIR}/source/postprocessing/solution_io.f90
  ${CMAKE_SOURCE_DIR}/source/postprocessing/post_utils.f90
  ${CMAKE_SOURCE_DIR}/source/postprocessing/visualization_output.f90
  )

set(Elements
  ${CMAKE_SOURCE_DIR}/source/Elements/e013.f
  ${CMAKE_SOURCE_DIR}/source/Elements/e012.f
  )

set(src_util
  ${CMAKE_SOURCE_DIR}/source/src_util/ReadExtrud3DParameters.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/types.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/OctTreeSearch.f90
  )

set(src_fbm
  ${CMAKE_SOURCE_DIR}/source/src_fbm/fbm_aux.f90
  ${CMAKE_SOURCE_DIR}/source/src_fbm/fbm_main.f90
)

set(src_assemblies
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_laplace.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_diff.f  
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_massrho.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_conv.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_shear.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_stress.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_surftens.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_gravity.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_BMatrix.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_barMmatrix.f
)

set(src_quadLS_app
${src_util}
${src_fbm}
${src_assemblies}
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_solver.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_proj.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_force.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_cylforce.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_torque.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_Sigma_User.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mpi.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_main.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_var.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_user.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mg.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_def.f90
${CMAKE_SOURCE_DIR}/source/initialization/app_initialization.f90
)


set(src_mpi
  ${CMAKE_SOURCE_DIR}/source/src_mpi/smooth_mpi.f
  ${CMAKE_SOURCE_DIR}/source/src_mpi/master_mpi.f
  )

set(src_LinSc
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_conv.f
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_diff.f
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_solver.f
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_def.f90
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_dfkt.f90
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_main.f90
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_mg.f90
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_mpi.f90
  )

set(src_PLin
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_def.f90
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_main.f90
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_matstr.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_mass.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_user.f90
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_conv.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_norm.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_flux.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_intphase.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_intpol.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_limiter.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_mpi.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_density.f90
  )

set(src_visco
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_integration.f
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_solver.f
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_main.f90
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_def.f90
  )

set(src_mesh
  ${CMAKE_SOURCE_DIR}/source/src_mesh/mesh_refine.f90
  ${CMAKE_SOURCE_DIR}/source/src_mesh/umbrella_smoother.f90
  ${CMAKE_SOURCE_DIR}/source/src_mesh/geometry_processing.f90
  )

# source files for standard q2p1 
set(src_q2p1
  ${CMAKE_SOURCE_DIR}/source/OutputProfiles.f90
  ${CMAKE_SOURCE_DIR}/source/tetra.f90
  ${CMAKE_SOURCE_DIR}/source/Init.f90
  ${CMAKE_SOURCE_DIR}/source/Statistics.f90
  ${CMAKE_SOURCE_DIR}/source/Umbrella.f90
  ${CMAKE_SOURCE_DIR}/source/src_mpi/pp3d_mpi.f90
  ${CMAKE_SOURCE_DIR}/source/ProcCtrl.f90
  ${CMAKE_SOURCE_DIR}/source/UmfpackSolver.f90
  ${CMAKE_SOURCE_DIR}/source/inverse.f90
  ${CMAKE_SOURCE_DIR}/source/Parametrization.f90
  ${CMAKE_SOURCE_DIR}/source/src_ini/iniparser.f90
  ${CMAKE_SOURCE_DIR}/source/distance.f90
  ${CMAKE_SOURCE_DIR}/source/3x3EigenV.f
  )


set(src_ini_c
  ${CMAKE_SOURCE_DIR}/source/src_ini/isdirectory.c
  ${CMAKE_SOURCE_DIR}/source/src_ini/mkdir_recursive.c
  )

set(src_particles
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_input.f90
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_step.f90
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_tracer.f90
  )
#=========================================================================
#         This will create source folders in Visual Studio
#=========================================================================
source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_particles FILES ${src_particles})
source_group(Elements FILES ${Elements})
source_group(postprocessing FILES ${postprocessing})

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
# Add the optical tweezers library
#=========================================================================
if(USE_OPTICALTWEEZERS)
  set(OPTICALTWEEZERS_LIBRARY True)

  ExternalProject_Add(OPTICALTWEEZERS_PRJ
    GIT_REPOSITORY ssh://rmuenste@arryn.mathematik.tu-dortmund.de/data/warehouse14/rmuenste/code/GitHub/raytracing-ot/raytracing 
    GIT_TAG mit_octree
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/libraries/opticaltweezers
    PREFIX ${CMAKE_SOURCE_DIR}/extern/libraries/opticaltweezers-dir
    #INSTALL_DIR ${CMAKE_SOURCE_DIR}/extern/libraries/cgal-install-dir 
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/extern/libraries/opticaltweezers-install -DBUILD_SHARED_LIBS=False -DOPTICALTWEEZERS_LIBRARY=True
    )

  set(OPTICALTWEEZERS_LIBRARIES 
      ${CMAKE_BINARY_DIR}/extern/libraries/opticaltweezers-install/lib/libvector.a
      ${CMAKE_BINARY_DIR}/extern/libraries/opticaltweezers-install/lib/libstrahl.a
      ${CMAKE_BINARY_DIR}/extern/libraries/opticaltweezers-install/lib/libot.a
     )

endif(USE_OPTICALTWEEZERS)

#=========================================================================
# Add the cgal library
#=========================================================================
set(CGAL_LIBRARIES 
  ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/lib64/libCGAL_Core.a
  ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/lib64/libCGAL.a
  )

set(CGAL_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/include)
if(USE_CGAL)
  set(USE_BOOST True)
  if(NOT WIN32)
        ExternalProject_Add(CGAL_PRJ
          GIT_REPOSITORY https://github.com/rmuenste/cgal.git 
          GIT_TAG rafa-devel
          SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/libraries/cgal
          PREFIX ${CMAKE_SOURCE_DIR}/extern/libraries/cgal-dir
          #INSTALL_DIR ${CMAKE_SOURCE_DIR}/extern/libraries/cgal-install-dir 
          CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/extern/libraries/cgal-install -DBUILD_SHARED_LIBS=False -DWITH_CGAL_Qt5=false -DWITH_CGAL_ImageIO=false
          )

    if(LSB_RELEASE_ID_SHORT STREQUAL "Ubuntu")
      set(CGAL_LIBRARIES 
        ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/lib/libCGAL_Core.a
        ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/lib/libCGAL.a
        )
    elseif(LSB_RELEASE_ID_SHORT STREQUAL "Scientific")
      set(CGAL_LIBRARIES 
        ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/lib64/libCGAL_Core.a
        ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/lib64/libCGAL.a
        )
      message(STATUS "CGAL path: ${CGAL_LIBRARIES} ")
    else()
      set(CGAL_LIBRARIES 
        ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/lib64/libCGAL_Core.a
        ${CMAKE_BINARY_DIR}/extern/libraries/cgal-install/lib64/libCGAL.a
        )
    endif()

  else(NOT WIN32)
    find_package(CGAL REQUIRED COMPONENTS Core)
    if(${CGAL_FOUND})
      find_package(GMP)
      find_package(MPFR)

      if(NOT ${GMP_FOUND})
        message(WARNING "GMP library not found. Cannot build cgal applications.")
      endif()

      if(NOT ${MPFR_FOUND})
        message(WARNING "MPFR library not found. Cannot build cgal applications.")
      endif()

    endif()	
    message("INCLUDE dirs  ${CGAL_INCLUDE_DIRS} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR}")
    message("LIBRARIES: ${CGAL_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES}") 
    set(CGAL_INCLUDE_DIR ${CGAL_INCLUDE_DIRS})
  endif(NOT WIN32)
endif(USE_CGAL)

if (CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
  FIND_LIBRARY(LIBRT_LIBRARY NAMES rt)
  if(NOT LIBRT_LIBRARY)
    find_file(LIBRT_LIBRARY librt.so  PATHS /usr/lib /usr/lib64/)
    if(NOT LIBRT_LIBRARY)
      message(FATAL_ERROR "You need to install the system library librt.")
    endif()
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
  # UNCOMMENT THE LINE BELOW IF YOU HAVE ISSUES FINDING BOOST
  # set(Boost_DEBUG True)
  find_package(Boost 1.53 REQUIRED COMPONENTS thread)
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

#if(Boost_FOUND)
#  message("Boost library version: ${Boost_VERSION}")
#  message("Boost detailed library version: ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
#  message("Boost include directory: ${Boost_INCLUDE_DIRS}")
#  message("Boost library directory: ${Boost_LIBRARY_DIRS}")
#  message("Boost library directory: ${Boost_LIBRARIES}")
#  message(FATAL_ERROR "Found Boost library")
#else()
#  message(FATAL_ERROR "Did not find Boost library")
#endif()
#
#message(FATAL_ERROR "end")

#=========================================================================
#               Build the C directory handling library  
#=========================================================================
add_library(cdirs ${src_ini_c})

#=========================================================================
#      Set the flags/paths/etc. we have determined for libraries 
#=========================================================================
include(${CMAKE_MODULE_PATH}/GenerateIncludeFlags.cmake)  
include(${CMAKE_MODULE_PATH}/GenerateLinkerFlags.cmake)  

#=========================================================================
#                         Add applications
#=========================================================================
SET(OUT_OF_SOURCE_BUILD True)
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
  SET(${OUT_OF_SOURCE_BUILD} False)
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

## we need additional libraries to support CUDA
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
#                       INSTALL CONFIGURATION
#=========================================================================
#install(DIRECTORY tools/partitioner DESTINATION bin)
file(GLOB E3D_SCRIPTS "tools/e3d_scripts/e3d_*")
install(FILES ${E3D_SCRIPTS} DESTINATION bin
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)

#=========================================================================
#                       CPACK CONFIGURATION
#=========================================================================

include(${CMAKE_MODULE_PATH}/GeneratePackageDeb.cmake)  

include(CPack)

