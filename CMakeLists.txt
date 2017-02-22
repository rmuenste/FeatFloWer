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

#--------------------------------------------------------------------------------------------------------------------
#                                               compiler settings
#--------------------------------------------------------------------------------------------------------------------

message(STATUS "Checking build type...")

IF(NOT Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER "gcc")
ELSE(NOT Q2P1_COMPILER)
  set(Q2P1_DEFAULT_COMPILER ${Q2P1_COMPILER})
ENDIF(NOT Q2P1_COMPILER)

include(${CMAKE_MODULE_PATH}/GenerateBuildId.cmake)  

#--------------------------------------------------------------------------------------------------------------------
#                                                   Options
#--------------------------------------------------------------------------------------------------------------------

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

#--------------------------------------------------------------------------------------------------------------------
#                                               project directories
#--------------------------------------------------------------------------------------------------------------------

#option(USE_FC_LIB
#    "Enables use FullC0ntact library"
#    OFF
#)


# source files for standad q2p1 
set(src_q2p1
source/Elements/e013.f
source/Elements/e012.f
source/src_pp3d/projma.f
source/src_pp3d/coeff.f
source/src_pp3d/indat3d.f
source/src_pp3d/gupwd.f
source/src_pp3d/nsdef.f
source/src_pp3d/bndry.f
source/src_pp3d/conv.f
source/src_pp3d/util.f
source/src_pp3d/avs3d.f
source/src_pp3d/period.f
source/src_pp3d/parq3d.f
source/src_pp3d/xmrout.f
source/src_pp3d/rdparm.f
source/src_pp3d/les.f
source/src_pp3d/trsort.f
source/src_pp3d/orsc.f
source/src_pp3d/mgrout.f
source/src_pp3d/inout.f
source/src_pp3d/user.f
source/src_pp3d/optcnl.f
source/src_pp3d/dfkt.f
source/src_pp3d/error.f
source/src_pp3d/bmul.f
source/src_pp3d/diff.f
source/src_pp3d/nu_turb.f
source/src_pp3d/prostp.f
#source/src_pproc/pproc_main.f90
source/src_quadLS/QuadSc_solver.f
source/src_quadLS/QuadSc_laplace.f
source/src_quadLS/QuadSc_diff.f
source/src_quadLS/QuadSc_massrho.f
source/src_quadLS/QuadSc_proj.f
source/src_quadLS/QuadSc_conv.f
source/src_quadLS/QuadSc_shear.f
source/src_quadLS/QuadSc_stress.f
source/src_quadLS/QuadSc_force.f90
#source/src_quadLS/QuadSc_force.f
source/src_quadLS/QuadSc_gravity.f
source/src_quadLS/QuadSc_BMatrix.f
source/src_mpi/smooth_mpi.f
source/src_mpi/master_mpi.f
source/src_LinSc/LinSc_conv.f
source/src_LinSc/LinSc_diff.f
source/src_LinSc/LinSc_solver.f
source/src_pp3d/def_feat.f90
source/OutputProfiles.f90
source/src_quadLS/QuadSc_mpi.f90
source/src_quadLS/QuadSc_main.f90
source/src_quadLS/QuadSc_var.f90
source/src_quadLS/QuadSc_user.f90
source/src_quadLS/QuadSc_mg.f90
source/src_quadLS/QuadSc_def.f90
source/tetra.f90
source/Init.f90
source/Statistics.f90
source/Umbrella.f90
source/src_mpi/pp3d_mpi.f90
source/src_LinSc/LinSc_def.f90
source/src_LinSc/LinSc_dfkt.f90
source/src_LinSc/LinSc_main.f90
source/src_PLin/PLinSc_def.f90
source/src_PLin/PLinSc_main.f90
source/src_PLin/PLinSc_matstr.f
source/src_PLin/PLinSc_mass.f
source/src_PLin/PLinSc_user.f90
source/src_PLin/PLinSc_conv.f
source/src_PLin/PLinSc_norm.f
source/src_PLin/PLinSc_flux.f
source/src_PLin/PLinSc_intphase.f
source/src_PLin/PLinSc_intpol.f
source/src_PLin/PLinSc_limiter.f
source/src_PLin/PLinSc_mpi.f
source/src_PLin/PLinSc_density.f90
source/ProcCtrl.f90
source/UmfpackSolver.f90
source/inverse.f90
source/Parametrization.f90
source/distance.f90
source/3x3EigenV.f
source/src_visco/Visco_integration.f
source/src_visco/Visco_solver.f
source/src_visco/Visco_main.f90
source/src_visco/Visco_def.f90
source/src_mesh/mesh_refine.f90
)

#--------------------------------------------------------------------------------------------------------------------
#                                    Configuration of the executable
#--------------------------------------------------------------------------------------------------------------------
add_executable(q2p1 ${src_q2p1}) 

#--------------------------------------------------------------------------------------------------------------------
#                                    search for the libraries and link them
#--------------------------------------------------------------------------------------------------------------------

#find_package(FEATFLOW_FEAT3D_LIBRARY REQUIRED)
#find_package(FEATFLOW_FEAT2D_LIBRARY REQUIRED)
#find_package(FEATFLOW_SYSUTILS_LIBRARY REQUIRED)

find_package(BLAS)
find_package(LAPACK)

if(NOT (${BLAS_FOUND} AND ${LAPACK_FOUND}))
  add_subdirectory(extern/libraries/lapack-3.6.1)
  set(BLAS_LIBRARIES blas)
  set(LAPACK_LIBRARIES lapack)
endif()


find_package(MPI)

message(STATUS "mpi include : ${MPI_Fortran_INCLUDE_PATH}")
message(STATUS "mpi libs : ${MPI_Fortran_LIBRARIES}")

add_subdirectory(extern/libraries/feat2d)
add_subdirectory(extern/libraries/feat3d)
add_subdirectory(extern/libraries/sysutils)
add_subdirectory(extern/libraries/amd)
add_subdirectory(extern/libraries/umfpack4)

# add subdirectories
IF(EXISTS "${CMAKE_SOURCE_DIR}/FullC0ntact")
  add_definitions(-DFEATFLOWLIB)
  message(STATUS "Found FullC0ntact library directory")  
  add_subdirectory(FullC0ntact/inshape3dcore)
  add_subdirectory(FullC0ntact/math)
  add_subdirectory(FullC0ntact/util)
  add_subdirectory(FullC0ntact/libs/OpenMesh)
ENDIF(EXISTS "${CMAKE_SOURCE_DIR}/FullC0ntact")

if (CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
  FIND_LIBRARY(LIBRT_LIBRARY NAMES rt)
  if(NOT LIBRT_LIBRARY)
    message(FATAL_ERROR "You need to install the system library librt.")
  endif()
endif()

target_link_libraries(q2p1 
  amd
  umfpack4
  feat2d
  feat3d
  sysutils 
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  inshape3dcore
  Utility
  Math
  OpenMeshCore 
  OpenMeshTools
  ${LIBRT_LIBRARY}
  ${MPI_Fortran_LIBRARIES}
  )


target_include_directories(q2p1 PUBLIC ${MPI_Fortran_INCLUDE_PATH})

target_compile_options(q2p1 PRIVATE ${Fortran_FLAGS})
#set_property(TARGET q2p1 PROPERTY LINKER_LANGUAGE Fortran)

#--------------------------------------------------------------------------------------------------------------------
#                                               Copy files for out-of-source build
#--------------------------------------------------------------------------------------------------------------------

IF(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  # check whether the $ENV{Q2P1_MESH_DIR} variable is set
  IF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    set(DIRECTORYLINKS meshes)
    message(STATUS "Found mesh directory : $ENV{Q2P1_MESH_DIR}")
    IF(NOT EXISTS ${CMAKE_BINARY_DIR}/_adc)
      message(STATUS "Found mesh directory creating symbolic link to $ENV{Q2P1_MESH_DIR}")
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        $ENV{Q2P1_MESH_DIR}
        ${CMAKE_BINARY_DIR}/_adc)
    ENDIF()
  ELSE(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    set(DIRECTORYLINKS _adc meshes)
  ENDIF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")

  set(FILELINKS partitioner)
  set(DIRECTORYCOPIES _data start)
  set(DEF_DIRECTORIES _dump _gmv _mesh _ns solution testresults _vtk)

  message(STATUS "Configuring out-of source build")

  file(COPY ${CMAKE_SOURCE_DIR}/_data DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/start DESTINATION ${CMAKE_BINARY_DIR})

  file(COPY ${CMAKE_SOURCE_DIR}/testcases DESTINATION ${CMAKE_BINARY_DIR})

  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fac/test1 DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fac_nnewt/test2 DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fallingparticle/test3 DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fac_visco/test4 DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fac_visco_elastic/test5 DESTINATION ${CMAKE_BINARY_DIR})

  # establish the file links
  FOREACH(flink ${FILELINKS})
    IF(NOT EXISTS ${CMAKE_BINARY_DIR}/${flink})
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
        ${CMAKE_SOURCE_DIR}/${flink}
        ${CMAKE_BINARY_DIR}/${flink})
    ENDIF ()    
  ENDFOREACH()

  # create the default directories if they are missing
  FOREACH(dir ${DEF_DIRECTORIES})
    IF(NOT EXISTS ${CMAKE_BINARY_DIR}/${dir})
      file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${dir})
    ENDIF()
  ENDFOREACH()

  # copy the data directories
  FOREACH(dir ${DIRECTORYCOPIES})
    IF(NOT EXISTS ${CMAKE_BINARY_DIR}/${dir})
      message(STATUS "${CMAKE_BINARY_DIR}/${dir} does not exist, copying directory")
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/${dir} ${CMAKE_BINARY_DIR}/${dir})
    ENDIF()
  ENDFOREACH()
ENDIF()

# continue to create the neccessary directory structure:
# we know the _dump directory exists now, we need to create the subdirectories 
set(SUBDIRS 01 02 03 04 05 06 07 08 09 10)
FOREACH(sd ${SUBDIRS})
  IF(NOT EXISTS ${CMAKE_BINARY_DIR}/_dump/${sd})
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_dump/${sd})
  ENDIF()  
ENDFOREACH()


#--------------------------------------------------------------------------------------------------------------------
#                                           Predefined test cases
#--------------------------------------------------------------------------------------------------------------------

add_test(q2p1-fac-newt ./test1)
add_test(q2p1-fac-nnewt ./test2)
add_test(q2p1-fallingparticle ./test3)
add_test(q2p1-visco-elastic ./test4)

