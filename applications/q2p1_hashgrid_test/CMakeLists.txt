include(CustomTargets)

# source files for q2p1_hashgrid_test
set(src_q2p1_hashgrid_test
app_init.f90
q2p1_hashgrid_test.f90
)

# Add the application executable
add_executable(q2p1_hashgrid_test ${src_q2p1_hashgrid_test})

# Add a custom command to update version.h before the build starts
add_custom_target(SetVersion_HashGridTest
    COMMAND ${CMAKE_COMMAND} -DGIT_DIR=${CMAKE_SOURCE_DIR} -P "${CMAKE_SOURCE_DIR}/cmake/modules/generate_version_h.cmake"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Updating version.h with the latest git commit hash"
)

add_dependencies(q2p1_hashgrid_test SetVersion_HashGridTest)

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

# Target related settings
target_link_libraries(q2p1_hashgrid_test
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_hashgrid_test PUBLIC ${FF_APPLICATION_INCLUDE_PATH} ${CMAKE_CURRENT_BINARY_DIR})

target_compile_options(q2p1_hashgrid_test PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_hashgrid_test PROPERTIES LINKER_LANGUAGE Fortran)

add_dependencies(q2p1_hashgrid_test metis)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

IF(${OUT_OF_SOURCE_BUILD})

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/partitioner")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/tools/partitioner" "${CMAKE_CURRENT_BINARY_DIR}/partitioner" )
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_data/MG.dat")
    file(COPY "${CMAKE_SOURCE_DIR}/_data/MG.dat" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/_data")
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/RankFileGenerator.sh")
    file(COPY "${CMAKE_SOURCE_DIR}/tools/RankFileGenerator.sh" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
  ENDIF()

ENDIF()
