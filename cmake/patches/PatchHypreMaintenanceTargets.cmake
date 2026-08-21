# Patch the pinned upstream hypre 2.33.0 source before it is embedded through
# FetchContent.  Custom target names are global in a CMake build, so hypre's
# maintenance targets must only exist when hypre is the top-level project.

if(NOT DEFINED FF_HYPRE_SOURCE_DIR OR FF_HYPRE_SOURCE_DIR STREQUAL "")
  message(FATAL_ERROR "FF_HYPRE_SOURCE_DIR must name the populated hypre source tree")
endif()

set(_hypre_cmakelists "${FF_HYPRE_SOURCE_DIR}/src/CMakeLists.txt")
if(NOT EXISTS "${_hypre_cmakelists}")
  message(FATAL_ERROR "Cannot patch hypre: ${_hypre_cmakelists} does not exist")
endif()

file(READ "${_hypre_cmakelists}" _hypre_contents)

set(_upstream_block
"# Add targets for tags, distclean, and uninstall
add_hypre_target_tags()
add_hypre_target_distclean()
add_hypre_target_uninstall()")

set(_patched_block
"# These maintenance targets are appropriate only for a standalone hypre build.
# Their global names may collide with targets owned by a superproject.
if(HYPRE_IS_TOP_LEVEL)
  add_hypre_target_tags()
  add_hypre_target_distclean()
  add_hypre_target_uninstall()
endif()")

string(FIND "${_hypre_contents}" "${_patched_block}" _patched_position)
if(NOT _patched_position EQUAL -1)
  message(STATUS "hypre maintenance targets are already restricted to top-level builds")
  return()
endif()

string(FIND "${_hypre_contents}" "${_upstream_block}" _upstream_position)
if(_upstream_position EQUAL -1)
  message(FATAL_ERROR
    "Cannot patch hypre maintenance targets: the expected hypre 2.33.0 block "
    "was not found in ${_hypre_cmakelists}")
endif()

string(REPLACE "${_upstream_block}" "${_patched_block}"
  _hypre_patched_contents "${_hypre_contents}")
file(WRITE "${_hypre_cmakelists}" "${_hypre_patched_contents}")
message(STATUS "Restricted hypre maintenance targets to top-level builds")
