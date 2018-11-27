# source files for standad q2p1
set(src_iniparserdemonstration
${CMAKE_SOURCE_DIR}/source/src_ini/iniparser.f90
parserdemonstration.f90
)

# Add the application executable
add_executable(parserdemonstration ${src_iniparserdemonstration})
target_link_libraries(parserdemonstration cdirs)


target_compile_options(parserdemonstration PRIVATE ${Fortran_FLAGS})

set_target_properties(parserdemonstration PROPERTIES LINKER_LANGUAGE Fortran)



IF(${OUT_OF_SOURCE_BUILD})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/testcfg.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()
