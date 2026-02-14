# CMakeLists of q2p1_fc_ext application

# source files for standad q2p1 
set(src_q2p1_fc_ext
app_init.f90
q2p1_bench_fluidization.f90
)

# Add the application executable
add_executable(q2p1_bench_fluidization ${src_q2p1_fc_ext}) 

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

target_link_libraries(q2p1_bench_fluidization 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_bench_fluidization PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_bench_fluidization PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_bench_fluidization PROPERTIES LINKER_LANGUAGE Fortran)

#message(STATUS "q2p1_fc_ext source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
#message(STATUS "q2p1_fc_ext bin dir: ${CMAKE_CURRENT_BINARY_DIR}")

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

set(Q2P1_BENCH_FLUIDIZATION_PARTS "31" CACHE STRING "Partition count for q2p1_bench_fluidization (CFD subdomains)")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/_data/MG.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/libs/pe/pe/interface/example.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/tools/dashboard/q2p1_ctest_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY cube.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(q2p1-bench-fluidization python ./q2p1_ctest_start.py)

if(NOT Python3_EXECUTABLE)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
endif()

add_custom_target(q2p1_bench_fluidization_stage
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_BINARY_DIR}/extern/libraries/metis-4.0.3/Lib/libmetis.so
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
          ${Q2P1_BENCH_FLUIDIZATION_PARTS} 1 1 NEWFAC _adc/benchSym/bench.prj
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_bench_fluidization metis
  COMMENT "Staging q2p1_bench_fluidization runtime files and partitioning mesh"
)
