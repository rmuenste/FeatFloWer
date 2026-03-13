# Add this directory for internal users.
include_directories(.)
# Find sources.
file(GLOB metis_sources *.c)
# Build libmetis.
add_library(metis ${METIS_LIBRARY_TYPE} ${GKlib_sources} ${metis_sources})
if(UNIX)
  target_link_libraries(metis m)
endif()

# Post-build: copy libmetis.so next to the q2p1_fc_ext application binary.
add_custom_command(
  TARGET metis
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory
    ${CMAKE_BINARY_DIR}/applications/q2p1_fc_ext
  COMMAND ${CMAKE_COMMAND} -E copy
    $<TARGET_FILE:metis>
    ${CMAKE_BINARY_DIR}/applications/q2p1_fc_ext
)

if(METIS_INSTALL)
  # Primary install location — RPATH is handled by CMake for this one.
  install(TARGETS metis
    LIBRARY DESTINATION lib64
    RUNTIME DESTINATION lib64
    ARCHIVE DESTINATION lib)

  # Additional install locations for FeatFloWer application directories.
  # Use install(FILES ...) to avoid multiple-RPATH conflicts.
  foreach(_dest bin/q2p1_sse bin/heat bin/q2p1_gendie bin/q2p1_sse_partfil bin/q2p1_die)
    install(
      FILES $<TARGET_FILE:metis>
      DESTINATION ${_dest}
    )
  endforeach()
endif()
