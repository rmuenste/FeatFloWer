cmake_minimum_required(VERSION 3.5)

project(UMFPACK C Fortran)

# allow less strict IF-else syntax
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

file(GLOB UMF_SRC ./src/*.c)

file(GLOB UMF_MOD_SRC ./modsrc/*.c)

add_library(umfpack4 ${UMF_SRC} ${UMF_MOD_SRC})

set(AMD_INCLUDE_DIR ../amd/include ../amd/src)

target_include_directories(umfpack4 PRIVATE ./include ./src ${AMD_INCLUDE_DIR})

target_compile_options(umfpack4 PRIVATE)

target_compile_definitions(umfpack4 PRIVATE INT=1 LP64=1 USE_NO_BLAS=1 DNTIMER=1)

target_link_libraries(umfpack4 PUBLIC amd)

