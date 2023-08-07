# generate_version_h.cmake

# Find Git to get the commit hash
find_package(Git)

if(NOT GIT_FOUND)
    message(FATAL_ERROR "Git is required to build this software. Please install Git.")
endif()

message("Source dir: ${GIT_DIR}")

# Custom command to extract git commit hash and store it in a variable
execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${GIT_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Set additional version information
set(PROJECT_VERSION "1.0.0")

# Write the version.h content using the template
# C/C++ Version
#file(WRITE ${CMAKE_BINARY_DIR}/version.h "#ifndef VERSION_H\n")
#file(APPEND ${CMAKE_BINARY_DIR}/version.h "#define VERSION_H\n")
#file(APPEND ${CMAKE_BINARY_DIR}/version.h "#define PROJECT_VERSION \"${PROJECT_VERSION}\"\n")
#file(APPEND ${CMAKE_BINARY_DIR}/version.h "#define GIT_COMMIT_HASH \"${GIT_COMMIT_HASH}\"\n")
#file(APPEND ${CMAKE_BINARY_DIR}/version.h "#endif\n")

# Fortran90 Version
file(WRITE ${CMAKE_BINARY_DIR}/version.h "  character(len=100), parameter :: PROJECT_VERSION = '${PROJECT_VERSION}'\n")
file(APPEND ${CMAKE_BINARY_DIR}/version.h "  character(len=100), parameter :: GIT_COMMIT_HASH = '${GIT_COMMIT_HASH}'\n")