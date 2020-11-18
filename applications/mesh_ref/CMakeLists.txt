# source files for standard q2p1 
set(src_meshref
${CMAKE_SOURCE_DIR}/source/ProcCtrl.f90
meshref_main.f90
meshref_def.f90
meshref_init.f90
meshref_var.f90
meshref_refine.f90
meshref_output.f90
)

# Add the application executable
add_executable(meshref ${src_meshref}) 

target_link_libraries(meshref 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(meshref PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(meshref PRIVATE ${Fortran_FLAGS})

set_target_properties(meshref PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

# file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/start)
# file(COPY Mesh/mesh.tri DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(COPY param.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY start/sampleRigidBody.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/start)


