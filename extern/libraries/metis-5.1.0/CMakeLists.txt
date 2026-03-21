cmake_minimum_required(VERSION 2.8.12...3.31)
project(METIS VERSION 5.1.0)

# GKlib is bundled in the source tree; use an absolute path so that
# the include() below resolves correctly when built as a subdirectory
# of a larger CMake project.
set(GKLIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/GKlib" CACHE PATH "path to GKlib")

# FeatFloWer always needs a shared library (loaded via ctypes from PyPartitioner).
set(METIS_LIBRARY_TYPE SHARED)
set(METIS_INSTALL TRUE)

include(${GKLIB_PATH}/GKlibSystem.cmake)

# Include directories for GKlib internals and the public METIS header.
include_directories(${GKLIB_PATH})
include_directories(include)

add_subdirectory("include")
add_subdirectory("libmetis")
# Note: "programs" subdirectory (gpmetis, ndmetis, …) is intentionally
# omitted — we only need the shared library.
# Post-build copy and install rules live in libmetis/CMakeLists.txt where the
# metis target is defined (required by CMake's add_custom_command scoping rules).
