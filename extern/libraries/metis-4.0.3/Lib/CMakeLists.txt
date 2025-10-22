
include_directories(.)

file(GLOB metis_SRC "*.h" "*.c")

add_library(metis SHARED ${metis_SRC})

add_custom_command(
  TARGET metis 
  POST_BUILD 
  COMMAND ${CMAKE_COMMAND} -E copy 
 ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so
 ${CMAKE_BINARY_DIR}/applications/q2p1_fc_ext
)

install(TARGETS metis
    RUNTIME DESTINATION bin/q2p1_sse
    LIBRARY DESTINATION bin/q2p1_sse
    ARCHIVE DESTINATION lib
    )

install(TARGETS metis
    RUNTIME DESTINATION bin/heat
    LIBRARY DESTINATION bin/heat
    ARCHIVE DESTINATION lib
    )

install(TARGETS metis
    RUNTIME DESTINATION bin/q2p1_gendie
    LIBRARY DESTINATION bin/q2p1_gendie
    ARCHIVE DESTINATION lib
    )

install(TARGETS metis
    RUNTIME DESTINATION bin/q2p1_sse_partfil
    LIBRARY DESTINATION bin/q2p1_sse_partfil
    ARCHIVE DESTINATION lib
    )

install(TARGETS metis
    RUNTIME DESTINATION bin/q2p1_die
    LIBRARY DESTINATION bin/q2p1_die
    ARCHIVE DESTINATION lib
    )
    
install(TARGETS metis
    RUNTIME DESTINATION lib64 
    LIBRARY DESTINATION lib64 
    ARCHIVE DESTINATION lib
    )
