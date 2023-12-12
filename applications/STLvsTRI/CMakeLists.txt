# CMakeLists of STLvsTRI application

# source files for standard q2p1 
set(src_STLvsTRI
STLvsTRI_def.f90
STLvsTRI_main.f90
)

# Add the application executable
add_executable(STLvsTRI ${src_STLvsTRI}) 

target_link_libraries(STLvsTRI 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(STLvsTRI PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(STLvsTRI PRIVATE ${Fortran_FLAGS})

set_target_properties(STLvsTRI PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${FF_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/IN)
file(COPY IN/mesh.tri DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/IN)
file(COPY IN/surface.off DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/IN)
file(COPY param_STLvsTRI.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

