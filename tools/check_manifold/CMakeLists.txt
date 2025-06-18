cmake_minimum_required(VERSION 3.8)

# name of the project, specIFy the used languages
project(CHECK_MANIFOLD C CXX)

file(COPY "cube.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "fail.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "fail.stl" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "wurfel.stl" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "check_manifold.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "blender_check_manifold.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# add a binary
add_executable(check_manifold cgal_test.cpp)
add_executable(validate_triangle_mesh validate_triangle_mesh.cpp)

get_target_property(MYTARGET check_manifold NAME)

target_link_libraries(${MYTARGET} ${CGAL_LIBRARIES}) 
target_link_libraries(validate_triangle_mesh ${CGAL_LIBRARIES}) 

# add include directories
target_include_directories(${MYTARGET} PUBLIC ${CGAL_INCLUDE_DIR})
target_include_directories(validate_triangle_mesh PUBLIC ${CGAL_INCLUDE_DIR})


