# cmake version
cmake_minimum_required(VERSION 3.2)

# name of the project, specIFy the used languages
project(Q2P1 C CXX Fortran)

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# allow less strict IF-else syntax
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# output the system name
message(STATUS "Configuring Q2P1 for a ${CMAKE_SYSTEM} system")

enable_testing()

include(CTest)

#-------------------------------------------------------------------------------------------------
#                                               compiler settings
#-------------------------------------------------------------------------------------------------

message(STATUS "Checking build type...")

IF(NOT Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER "gcc")
ELSE(NOT Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER ${Q2P1_COMPILER})
ENDIF(NOT Q2P1_COMPILER)

include(${CMAKE_MODULE_PATH}/GenerateBuildId.cmake)  

#-------------------------------------------------------------------------------------------------
#                                                   Options
#-------------------------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------------------
#                                               project directories
#-------------------------------------------------------------------------------------------------

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
${CMAKE_SOURCE_DIR}/source/postprocessing/solution_output.f90
)

set(Elements
${CMAKE_SOURCE_DIR}/source/Elements/e013.f
${CMAKE_SOURCE_DIR}/source/Elements/e012.f
)

set(src_quadLS
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_solver.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_laplace.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_diff.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_massrho.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_proj.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_conv.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_shear.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_stress.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_surftens.f
#${CMAKE_SOURCE_DIR/}source/src_quadLS/QuadSc_force.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_force.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_gravity.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_BMatrix.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mpi.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_main.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_var.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_user.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mg.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_def.f90
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
${CMAKE_SOURCE_DIR}/source/distance.f90
${CMAKE_SOURCE_DIR}/source/3x3EigenV.f
)

#-------------------------------------------------------------------------------------------------
#                                 This will create source folders in Visual Studio
#-------------------------------------------------------------------------------------------------

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(Elements FILES ${Elements})
source_group(postprocessing FILES ${postprocessing})

#-------------------------------------------------------------------------------------------------
#                                    Configure BLAS/LAPACK 
#-------------------------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------------------
#                           Configure MPI/MKL/MUMPS 
#-------------------------------------------------------------------------------------------------

find_package(MPI)

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

if(BUILD_METIS)
  add_subdirectory(extern/libraries/metis-4.0.3)
endif(BUILD_METIS)

if(USE_OPTICALTWEEZERS)
  add_subdirectory(extern/libraries/opticaltweezers)
endif(USE_OPTICALTWEEZERS)

if (CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
  FIND_LIBRARY(LIBRT_LIBRARY NAMES rt)
  if(NOT LIBRT_LIBRARY)
    find_file(LIBRT_LIBRARY librt.so  PATHS /usr/lib /usr/lib64/)
    if(NOT LIBRT_LIBRARY)
      message(FATAL_ERROR "You need to install the system library librt.")
    endif()
  endif()
endif()

#-------------------------------------------------------------------------------------------------
#                               Configure BoostC++ 
#-------------------------------------------------------------------------------------------------

#find_package(Boost 1.56 COMPONENTS filesystem REQUIRED)
#-------------------------------------------------------------------------------------------------
#                      Set the flags/paths/etc. we have determined for libraries 
#-------------------------------------------------------------------------------------------------

include(${CMAKE_MODULE_PATH}/GenerateIncludeFlags.cmake)  
include(${CMAKE_MODULE_PATH}/GenerateLinkerFlags.cmake)  
#-------------------------------------------------------------------------------------------------
#                                    Configuration of the executable
#-------------------------------------------------------------------------------------------------

#add_executable(q2p1 ${Elements} ${src_q2p1}
#  ${src_pp3d} ${src_mpi} ${src_PLin} ${src_LinSc} 
#  ${src_quadLS} ${src_visco} ${src_mesh}
#  ${postprocessing}
#  ${CMAKE_SOURCE_DIR}/source/src_pp3d/les.f) 
#
#target_link_libraries(q2p1 
#${FF_DEFAULT_LIBS}
#)
#
#target_include_directories(q2p1 PUBLIC ${FF_DEFAULT_INCLUDE_PATH})
#
#target_compile_options(q2p1 PRIVATE ${Fortran_FLAGS})
#-------------------------------------------------------------------------------------------------
#                                  Copy files for out-of-source build
#-------------------------------------------------------------------------------------------------

include(${CMAKE_MODULE_PATH}/CreateWorkingDirectories.cmake)
#-------------------------------------------------------------------------------------------------
#                                        Add applications
#-------------------------------------------------------------------------------------------------

## we need additional libraries to support CUDA
IF(BUILD_APPLICATIONS)
  # add subdirectories
  IF(EXISTS "${CMAKE_SOURCE_DIR}/FullC0ntact")
    add_definitions(-DFEATFLOWLIB)
    message(STATUS "Found FullC0ntact library directory")  
    add_subdirectory(FullC0ntact/inshape3dcore)
    add_subdirectory(FullC0ntact/math)
    add_subdirectory(FullC0ntact/util)
  if(USE_OPENMESH)
    add_subdirectory(FullC0ntact/libs/OpenMesh)
  endif(USE_OPENMESH)
  if(USE_ODE)
    add_subdirectory(FullC0ntact/libs/ode-0.13-cmake)
  endif(USE_ODE)
  else()
    message(FATAL_ERROR "FullC0ntact directory is required for applications.")  
  ENDIF()

  add_subdirectory(applications)
ENDIF(BUILD_APPLICATIONS)
#-------------------------------------------------------------------------------------------------
#                                       INSTALL CONFIGURATION
#-------------------------------------------------------------------------------------------------

#include(${CMAKE_MODULE_PATH}/ConfigureInstallDirectories.cmake)

