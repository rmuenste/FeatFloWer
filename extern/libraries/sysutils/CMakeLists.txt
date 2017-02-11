cmake_minimum_required(VERSION 3.2)

project(SYSUTILS Fortran)

# allow less strict IF-else syntax
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

set(SYSUTILS_SRC ./src/ztime.f)

add_library(sysutils ${SYSUTILS_SRC})

target_compile_options(sysutils PRIVATE -O3 )

