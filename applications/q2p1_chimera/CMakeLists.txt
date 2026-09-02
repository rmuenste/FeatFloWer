include(CustomTargets)

# Chimera overlapping-mesh application (chimera-integration-design.md v3,
# section 8).  Same scaffold as q2p1_fc_ext; the vendored Chimera FAC
# case lives in _data/CHIMERA_FAC (channel background mesh from
# tools/chimera_meshgen/channel_tri.py; the coarse annulus shell is the
# committed Phase-2 fixture, staged here at configure time).
set(src_q2p1_chimera
app_init.f90
q2p1_chimera.f90
)

add_executable(q2p1_chimera ${src_q2p1_chimera})

add_custom_target(SetVersionChimera
    COMMAND ${CMAKE_COMMAND} -DGIT_DIR=${CMAKE_SOURCE_DIR} -P "${CMAKE_SOURCE_DIR}/cmake/modules/generate_version_h.cmake"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Updating version.h with the latest git commit hash"
)

add_dependencies(q2p1_chimera SetVersionChimera)

target_link_libraries(q2p1_chimera
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_chimera PUBLIC ${FF_APPLICATION_INCLUDE_PATH} ${CMAKE_CURRENT_BINARY_DIR})

target_compile_options(q2p1_chimera PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_chimera PROPERTIES LINKER_LANGUAGE Fortran)

add_dependencies(q2p1_chimera metis)
add_custom_command(TARGET q2p1_chimera POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          $<TARGET_FILE:metis> ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

# Stage the coarse annulus shell (Phase-2 fixture) next to the vendored
# background mesh; the deck references _data/CHIMERA_FAC/annulus_coarse.tri.
file(COPY ${CMAKE_SOURCE_DIR}/source/src_chimera/tests/fixtures/annulus_coarse.tri
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data/CHIMERA_FAC)
file(COPY ${CMAKE_SOURCE_DIR}/tools/partitioner DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

IF(${OUT_OF_SOURCE_BUILD})

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_data/MG.dat")
    file(COPY "${CMAKE_SOURCE_DIR}/_data/MG.dat" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/_data")
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/RankFileGenerator.sh")
    file(COPY "${CMAKE_SOURCE_DIR}/tools/RankFileGenerator.sh" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
  ENDIF()

ENDIF()

# Runtime staging (METIS next to the binary), as recommended by
# docs/md_docs/guide_05_agent_new_application_workflow.md
add_custom_target(q2p1_chimera_stage
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          $<TARGET_FILE:metis> ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_chimera metis
  COMMENT "Staging q2p1_chimera runtime files in build directory"
)
