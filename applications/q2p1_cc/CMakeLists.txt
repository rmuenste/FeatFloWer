
# CC-specific sources. Core FEATFLOW sources are provided by
# FF_APPLICATION_LIBS; compiling them again here creates duplicate COMMON
# blocks and module state in the executable.
set(src_q2p1_cc
app_init.f90
assemblies_cc.f
iso_assemblies.f
postprocessing.f90
lin_transport_cc.f90
q2p1_def_cc.f90
q2p1_mg_cc.f90
q2p1_cc_umfpacksolver.f90
q2p1_cc.f90
q2p1_transport_cc.f90
q2p1_var_newton.f90
)

# Add the application executable
add_executable(q2p1_cc ${src_q2p1_cc}) 
add_dependencies(q2p1_cc metis)
add_custom_command(TARGET q2p1_cc POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          $<TARGET_FILE:metis> ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so)

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_ext})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(Elements FILES ${Elements})

target_link_libraries(q2p1_cc 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_cc PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_cc PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_cc PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/_data/q2p1_param.dat
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/_data/MG.dat
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/tools/partitioner
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/tools/PyPartitioner.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY q2p1_cc.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                      GROUP_READ GROUP_EXECUTE
                      WORLD_READ WORLD_EXECUTE)

# Stationary Re-20 cylinder regression. The tests run in a dedicated
# hermetic working directory whose _adc always points at the vendored
# source-tree case, independent of Q2P1_MESH_DIR, so the force baseline
# cannot drift with an external mesh repository. Interactive runs in the
# application build directory keep the env-preferred _adc behavior.
if(BUILD_TESTING AND NOT WIN32)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    set(_q2p1_cc_python ${Python3_EXECUTABLE})
  else()
    set(_q2p1_cc_python python3)
  endif()

  set(_q2p1_cc_test_dir ${CMAKE_CURRENT_BINARY_DIR}/ctest_run)
  file(MAKE_DIRECTORY ${_q2p1_cc_test_dir})
  createDefaultDirectories(${_q2p1_cc_test_dir} ${CMAKE_CURRENT_SOURCE_DIR})
  # Force the vendored mesh data regardless of Q2P1_MESH_DIR.
  file(REMOVE ${_q2p1_cc_test_dir}/_adc)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                  ${CMAKE_SOURCE_DIR}/_adc ${_q2p1_cc_test_dir}/_adc)
  # Refresh the deck and tools on every configure.
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/_data/q2p1_param.dat
       DESTINATION ${_q2p1_cc_test_dir}/_data)
  file(COPY ${CMAKE_SOURCE_DIR}/_data/MG.dat
       DESTINATION ${_q2p1_cc_test_dir}/_data)
  file(COPY ${CMAKE_SOURCE_DIR}/tools/partitioner
       DESTINATION ${_q2p1_cc_test_dir})
  file(COPY ${CMAKE_SOURCE_DIR}/tools/PyPartitioner.py
       DESTINATION ${_q2p1_cc_test_dir})
  file(COPY q2p1_cc.py DESTINATION ${_q2p1_cc_test_dir}
       FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                        GROUP_READ GROUP_EXECUTE
                        WORLD_READ WORLD_EXECUTE)
  # Executable and METIS library are build products; link them (dangling
  # until the first build, which is fine).
  foreach(_lnk q2p1_cc libmetis.so)
    file(REMOVE ${_q2p1_cc_test_dir}/${_lnk})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${CMAKE_CURRENT_BINARY_DIR}/${_lnk}
                    ${_q2p1_cc_test_dir}/${_lnk})
  endforeach()

  add_test(NAME q2p1_cc_cylinder_stationary
           COMMAND ${_q2p1_cc_python} ${_q2p1_cc_test_dir}/q2p1_cc.py
                   -n 4 --log ctest_q2p1_cc.log
                   --expect-drag 5.5981911 --expect-lift 0.0099684876
           WORKING_DIRECTORY ${_q2p1_cc_test_dir})
  set_tests_properties(q2p1_cc_cylinder_stationary PROPERTIES
                       TIMEOUT 900
                       PROCESSORS 4
                       LABELS "mpi;benchmark")

  # Fail-closed guard: a NaN in the reported forces must be rejected by the
  # launcher's validator, never accepted as a passing benchmark.
  file(WRITE ${_q2p1_cc_test_dir}/ctest_nan_sample.log
       "BenchForce:  0.10000000E+01             NaN             NaN\n")
  add_test(NAME q2p1_cc_nan_guard
           COMMAND ${_q2p1_cc_python} ${_q2p1_cc_test_dir}/q2p1_cc.py
                   --validate-log-only ${_q2p1_cc_test_dir}/ctest_nan_sample.log
                   --expect-drag 5.5981911 --expect-lift 0.0099684876
           WORKING_DIRECTORY ${_q2p1_cc_test_dir})
  set_tests_properties(q2p1_cc_nan_guard PROPERTIES
                       PASS_REGULAR_EXPRESSION "FAILED: non-finite forces"
                       TIMEOUT 60
                       LABELS "unit")
endif()
