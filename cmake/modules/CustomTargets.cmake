#=========================================================================
#                         custom commands
#=========================================================================
set(CUSTOM_DEPENDENCY "my_file.stamp")
add_custom_command(
    COMMAND ${CMAKE_COMMAND} -E touch ${CUSTOM_DEPENDENCY} 
    OUTPUT ${CUSTOM_DEPENDENCY}
)

add_custom_target(custom-target DEPENDS ${CUSTOM_DEPENDENCY})
