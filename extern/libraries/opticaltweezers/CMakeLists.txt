cmake_minimum_required(VERSION 3.2)

# name of the project, specIFy the used languages
project(OpticalTweezers)

# allow easier IF-else syntax
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

add_subdirectory(vector)
add_subdirectory(strahl)
#add_subdirectory(ot)

include_directories(./vector ./strahl .)

set(ottest_sources
ottest.cpp
PStrahl.cpp
trace.cpp
surfArray.cpp
  )

#add_executable(${PROJECT_NAME} ${ottest_sources})
add_library(opticaltweezers ${ottest_sources})

target_link_libraries(opticaltweezers strahl vector)
