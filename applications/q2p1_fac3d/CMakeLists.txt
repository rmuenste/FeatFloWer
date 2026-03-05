include(CustomTargets)  

# source files for standard q2p1 
set(src_q2p1_fac3d
app_init.f90
q2p1_fac3d.f90
)

# Add the application executable
add_executable(q2p1_fac3d ${src_q2p1_fac3d}) 

# Add a custom command to update version.h before the build starts
add_custom_target(SetVersion_FAC3D
    COMMAND ${CMAKE_COMMAND} -DGIT_DIR=${CMAKE_SOURCE_DIR} -P "${CMAKE_SOURCE_DIR}/cmake/modules/generate_version_h.cmake"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Updating version.h with the latest git commit hash"
)

add_dependencies(q2p1_fac3d SetVersion_FAC3D)

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
target_link_libraries(q2p1_fac3d 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_fac3d PUBLIC ${FF_APPLICATION_INCLUDE_PATH} ${CMAKE_CURRENT_BINARY_DIR})

target_compile_options(q2p1_fac3d PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_fac3d PROPERTIES LINKER_LANGUAGE Fortran)

add_dependencies(q2p1_fac3d metis)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/q2p1_fac3d_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/tools/dashboard/q2p1_ctest_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/tools/partitioner DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(q2p1-fac3d python ./q2p1_ctest_start.py)

add_custom_target(q2p1_fac3d_stage
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_BINARY_DIR}/extern/libraries/metis-4.0.3/Lib/libmetis.so
          ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_fac3d metis
  COMMENT "Staging q2p1_fac3d runtime files in build directory"
)

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
