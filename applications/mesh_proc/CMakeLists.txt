# CMakeLists of meshproc application
set(src_cinterface
${CMAKE_SOURCE_DIR}/source/cinterface/cinterface.f90
)

set(src_quadLS_meshproc
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
)

# source files for standard q2p1 
set(src_meshproc
${Elements} 
${src_q2p1} 
${src_pp3d} 
${src_mpi} 
${src_PLin} 
${src_LinSc} 
${src_quadLS_meshproc} 
${src_visco} 
${src_mesh}
${src_cinterface}
${postprocessing}
meshproc_main.f90
meshproc_def.f90
)

# Add the application executable
add_executable(meshproc ${src_meshproc}) 

target_link_libraries(meshproc 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(meshproc PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(meshproc PRIVATE ${Fortran_FLAGS})

set_target_properties(meshproc PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

# file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
# file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh/PAR)
# file(COPY Mesh/mesh.tri DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(COPY param.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


