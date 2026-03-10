include(${CMAKE_SOURCE_DIR}/cmake/modules/CustomTargets.cmake)  

# source files for standard q2p1 
set(src_q2p1_dns_drag
app_init.f90
q2p1_dns_drag.f90
)

#add_custom_target(output-build-fc-ext DEPENDS ${CUSTOM_DEPENDENCY})

# Add the application executable
#add_executable(q2p1_fc_ext ${src_q2p1_fc_ext} ${CUSTOM_DEPENDENCY}) 
add_executable(q2p1_dns_drag ${src_q2p1_dns_drag}) 

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
target_link_libraries(q2p1_dns_drag 
  ${FF_APPLICATION_LIBS}
  )
target_include_directories(q2p1_dns_drag PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_dns_drag PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_dns_drag PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

IF(${OUT_OF_SOURCE_BUILD})
  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/example.json")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/example.json" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
  ENDIF()
ENDIF()

add_custom_target(q2p1_dns_drag_stage
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_BINARY_DIR}/extern/libraries/metis-4.0.3/Lib/libmetis.so
          ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_dns_drag metis
  COMMENT "Staging q2p1_dns_drag runtime files in build directory"
)
