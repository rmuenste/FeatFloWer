# source files for standard q2p1 
set(src_meshproc
${CMAKE_SOURCE_DIR}/source/ProcCtrl.f90
meshproc_main.f90
meshproc_def.f90
meshproc_PDE.f90
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
  include(${FF_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

# file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/start)
# file(COPY Mesh/mesh.tri DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(COPY param.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY start/sampleRigidBody.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/start)


