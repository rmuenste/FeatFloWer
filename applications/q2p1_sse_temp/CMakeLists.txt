# CMakeLists of q2p1_sse_temp application

# source files for standard q2p1 
set(src_q2p1_sse_temp
app_init.f90
q2p1_sse_temp.f90
)


# Add the application executable
add_executable(q2p1_sse_temp ${src_q2p1_sse_temp}) 

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

if(Q1P1_SSE_LINK_STATIC_FORTRAN)
  target_link_libraries(q2p1_sse_temp 
    ${FF_APPLICATION_LIBS} -static-libgfortran
    )
else(Q1P1_SSE_LINK_STATIC_FORTRAN)
  target_link_libraries(q2p1_sse_temp 
    ${FF_APPLICATION_LIBS}
    )
endif(Q1P1_SSE_LINK_STATIC_FORTRAN)

target_include_directories(q2p1_sse_temp PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_sse_temp PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_sse_temp PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY q2p1_sse_temp.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(TARGETS q2p1_sse_temp
    RUNTIME DESTINATION bin/q2p1_sse
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(TARGETS q2p1_sse_temp
    RUNTIME DESTINATION bin/q2p1_gendie
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
