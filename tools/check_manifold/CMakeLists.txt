cmake_minimum_required(VERSION 3.12)

# name of the project, specIFy the used languages
project(CHECK_MANIFOLD C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(COPY "cube.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "fail.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "fail.stl" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "wurfel.stl" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "check_manifold.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "blender_check_manifold.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy test driver files
file(COPY "run_tests.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "test_validator.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy test files directory
file(COPY "test_files" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# add a binary
add_executable(check_manifold cgal_test.cpp)

# Create executable
add_executable(validate_triangle_mesh validate_triangle_mesh.cpp)

get_target_property(MYTARGET check_manifold NAME)

# TODO:
# This modern syntax can already be used with the USE_CGAL_LOCAL=ON option.
# With the ExternalProject_add mechanism this syntax will fail as the CGAL::CGAL target is not
# available at configure time.
#target_link_libraries(${MYTARGET} PRIVATE CGAL::CGAL CGAL::CGAL_Core)
#target_link_libraries(validate_triangle_mesh PRIVATE CGAL::CGAL CGAL::CGAL_Core)

# add include directories
target_include_directories(${MYTARGET} PUBLIC ${CGAL_INCLUDE_DIR})
target_include_directories(validate_triangle_mesh PUBLIC ${CGAL_INCLUDE_DIR})
