
include_directories(.)

file(GLOB metis_SRC "*.h" "*.c")

add_library(metis SHARED ${metis_SRC})

install(TARGETS metis
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION bin
    ARCHIVE DESTINATION lib
    )
