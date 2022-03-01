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

IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/PATCHES_C")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/PATCHES_C" "${CMAKE_CURRENT_BINARY_DIR}/PATCHES_C" )
ENDIF()

IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/PATCHES")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/PATCHES" "${CMAKE_CURRENT_BINARY_DIR}/PATCHES" )
ENDIF()

IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/start")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/start" "${CMAKE_CURRENT_BINARY_DIR}/start" )
ENDIF()

IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/VEKA_S")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/testcases/gendie/VEKA_S" "${CMAKE_CURRENT_BINARY_DIR}/VEKA_S" )
ENDIF()

# file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
# file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/start)
# file(COPY Mesh/mesh.tri DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(COPY param_meshref.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# file(COPY param_BU.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY start/sampleRigidBody.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/start)


