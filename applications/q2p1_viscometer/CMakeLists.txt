# CMakeLists of q2p1_viscometer application
# D5.1 numerical viscometer (Searle cell: stationary vessel, rotating bore).
# Cloned from applications/q2p1_dkt (serial PE + HardContactAndFluid +
# external particles.xyz packing).

# source files for standard q2p1
set(src_q2p1_viscometer
app_init.f90
q2p1_viscometer.f90
)

# Add the application executable
add_executable(q2p1_viscometer ${src_q2p1_viscometer})

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

target_link_libraries(q2p1_viscometer
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_viscometer PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_viscometer PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_viscometer PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

# Partition count for the D5.1 mesh (walled case -> METIS is acceptable,
# axis-uniform/Cartesian partitioning is not required here).
set(Q2P1_VISCOMETER_PARTS "108" CACHE STRING "Partition count for q2p1_viscometer (CFD subdomains)")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/_data/MG.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/_data/q2p1_param.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data)

file(COPY ${CMAKE_SOURCE_DIR}/tools/partitioner DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/example.json")
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/example.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()
IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/particles.xyz")
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/particles.xyz DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()

# Runtime staging.
#
# NOTE: the D5.1 mesh (SimPar@MeshFolder = "VISCO_108") is produced by the
# meshing workflow and is intentionally NOT part of this application
# directory.  The stage target therefore only copies the runtime files it
# owns; partition the mesh separately, e.g.
#
#   python3 tools/PyPartitioner.py ${Q2P1_VISCOMETER_PARTS} 1 1 VISCO_108 <path>/visco.prj
#
add_custom_target(q2p1_viscometer_stage
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_BINARY_DIR}/extern/libraries/metis-5.1.0/libmetis/libmetis.so
          ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_SOURCE_DIR}/_data/MG.dat
          ${CMAKE_CURRENT_BINARY_DIR}/_data/MG.dat
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_CURRENT_SOURCE_DIR}/_data/q2p1_param.dat
          ${CMAKE_CURRENT_BINARY_DIR}/_data/q2p1_param.dat
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_CURRENT_SOURCE_DIR}/example.json
          ${CMAKE_CURRENT_BINARY_DIR}/example.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_viscometer metis
  COMMENT "Staging q2p1_viscometer runtime files in build directory"
)
