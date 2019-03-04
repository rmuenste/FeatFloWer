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

#=========================================================================
#                         UnZip custom target
#=========================================================================
#add_custom_target( unZip ALL)
#add_custom_command(TARGET unZip PRE_BUILD
#   COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/abc/
#   COMMAND ${CMAKE_COMMAND} -E tar xzf {CMAKE_SOURCE_DIR}/abc.zip
#WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#DEPENDS ${CMAKE_SOURCE_DIR}/abc.zip
#COMMENT "Unpacking abc.zip"
#VERBATIM)