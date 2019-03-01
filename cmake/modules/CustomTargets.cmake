#=========================================================================
#                         custom commands
#=========================================================================
set(CUSTOM_DEPENDENCY "BuildID-Output")
add_custom_command(
#    COMMAND ${CMAKE_COMMAND} -E touch ${CUSTOM_DEPENDENCY} 
    COMMAND python ${CMAKE_SOURCE_DIR}/tools/cmake-py/get_build_ids.py 
    OUTPUT ${CUSTOM_DEPENDENCY}
)

add_custom_target(output-build-ids DEPENDS ${CUSTOM_DEPENDENCY})
