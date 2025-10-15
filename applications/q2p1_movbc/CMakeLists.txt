include(CustomTargets)  

# source files for standad q2p1 
set(src_q2p1_fc_ext
app_init.f90
#postprocessing.f90
q2p1_movbc.f90
)

#add_custom_target(output-build-fc-ext DEPENDS ${CUSTOM_DEPENDENCY})

# Add the application executable
#add_executable(q2p1_movbc ${src_q2p1_fc_ext} ${CUSTOM_DEPENDENCY}) 
add_executable(q2p1_movbc ${src_q2p1_fc_ext}) 

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
target_link_libraries(q2p1_movbc 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_movbc PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_movbc PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_movbc PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/q2p1_movbc_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
