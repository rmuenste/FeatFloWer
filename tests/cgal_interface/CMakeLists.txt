set(TEST_LIBS
    stdc++
    stdc++fs
    inshape3dcore
    Utility
    Math
)

add_executable(off_processing off_processing.f90)
target_link_libraries(off_processing ${TEST_LIBS})

target_compile_options(off_processing PRIVATE ${Fortran_FLAGS})

set_target_properties(off_processing PROPERTIES LINKER_LANGUAGE Fortran)
