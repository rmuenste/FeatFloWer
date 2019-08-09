#=========================================================================
#                         custom commands
#=========================================================================
#set(CUSTOM_DEPENDENCY "BuildID-Output")
#add_custom_command(
##    COMMAND ${CMAKE_COMMAND} -E touch ${CUSTOM_DEPENDENCY} 
#    COMMAND python ${CMAKE_SOURCE_DIR}/tools/cmake-py/get_build_ids.py 
#    OUTPUT ${CUSTOM_DEPENDENCY}
#)
#
#add_custom_target(output-build-ids DEPENDS ${CUSTOM_DEPENDENCY})

#=========================================================================
#           A function that adds a custom command to a target 
#=========================================================================
function(addCompilerInfo targetName)

  add_custom_command(
    TARGET ${targetName}
    POST_BUILD
    COMMAND module list |& cat > compiler-info.txt 
    COMMENT "Generating compiler-info"
  )

endfunction(addCompilerInfo targetName)

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


#set(CUSTOM_DEPENDENCY "compiler-info.txt")
#
#add_custom_command(
#  COMMAND module list |& cat > compiler-info.txt 
#  COMMENT "Generating compiler-info"
#  OUTPUT ${CUSTOM_DEPENDENCY}
#)
