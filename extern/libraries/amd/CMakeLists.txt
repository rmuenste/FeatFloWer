cmake_minimum_required(VERSION 3.2)

project(AMD C Fortran)

# allow less strict IF-else syntax
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)


file(GLOB AMD_SRC ./src/*.c)

add_library(amd ${AMD_SRC})

target_include_directories(amd PRIVATE ./include)

target_compile_options(amd PRIVATE -O3 )

target_compile_definitions(amd PRIVATE INT=1 LP64=1)

