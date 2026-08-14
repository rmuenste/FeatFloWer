# CMakeLists of q2p1_dkt application (drafting-kissing-tumbling benchmark)

# source files for standad q2p1
set(src_q2p1_dkt
app_init.f90
q2p1_dkt.f90
)

# Add the application executable
add_executable(q2p1_dkt ${src_q2p1_dkt})

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_ext})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_cinterface FILES ${src_cinterface})
source_group(Elements FILES ${Elements})

target_link_libraries(q2p1_dkt 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_dkt PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_dkt PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_dkt PROPERTIES LINKER_LANGUAGE Fortran)

#message(STATUS "q2p1_fc_ext source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
#message(STATUS "q2p1_fc_ext bin dir: ${CMAKE_CURRENT_BINARY_DIR}")

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

set(Q2P1_DKT_PARTS "31" CACHE STRING "Partition count for q2p1_dkt (CFD subdomains)")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/_data/MG.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/libs/pe/pe/interface/example.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_SOURCE_DIR}/tools/partitioner DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY cube.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(NOT Python3_EXECUTABLE)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
endif()

add_custom_target(q2p1_dkt_stage
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_BINARY_DIR}/extern/libraries/metis-5.1.0/libmetis/libmetis.so
          ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_SOURCE_DIR}/libs/pe/pe/interface/example.json
          ${CMAKE_CURRENT_BINARY_DIR}/example.json
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_SOURCE_DIR}/_data/MG.dat
          ${CMAKE_CURRENT_BINARY_DIR}/_data/MG.dat
  COMMAND ${CMAKE_COMMAND} -E env
          "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{LD_LIBRARY_PATH}"
          ${Python3_EXECUTABLE}
          ${CMAKE_SOURCE_DIR}/tools/PyPartitioner.py
          ${Q2P1_DKT_PARTS} 1 1 NEWFAC _adc/benchSym/bench.prj
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_dkt metis
  COMMENT "Staging q2p1_dkt runtime files and partitioning mesh"
)
