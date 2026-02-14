add_executable(fem_sandbox
  ${CMAKE_CURRENT_SOURCE_DIR}/fem_sandbox.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/feat_stubs.f90
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/def_feat.f90
  ${CMAKE_SOURCE_DIR}/source/Elements/e013.f
)

target_compile_options(fem_sandbox PRIVATE ${Fortran_FLAGS})

target_include_directories(fem_sandbox PRIVATE ${FF_APPLICATION_INCLUDE_PATH})

set_target_properties(fem_sandbox PROPERTIES LINKER_LANGUAGE Fortran)
