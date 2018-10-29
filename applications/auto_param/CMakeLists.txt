# CMakeLists of autoparam application

# source files for standard q2p1 
set(src_autoparam
#${Elements} 
#${src_q2p1} 
#${src_pp3d} 
#${src_mpi} 
#${src_PLin} 
#${src_LinSc} 
#${src_quadLS_mesher} 
#${src_visco} 
#${src_mesh}
#${src_cinterface}
#${postprocessing}
bndry.f90
meshproc_def.f90
meshproc_main.f90
umbrella.f90
)

# Add the application executable
add_executable(autoparam ${src_autoparam}) 

target_link_libraries(autoparam 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(autoparam PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(autoparam PRIVATE ${Fortran_FLAGS})

set_target_properties(autoparam PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh/PAR)
file(COPY Mesh/mesh.tri DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(COPY param.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


