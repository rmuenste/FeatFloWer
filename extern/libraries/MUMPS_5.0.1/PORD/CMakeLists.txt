
file(GLOB PORD_SRC ./lib/*.c)

add_library(pord ${PORD_SRC})

target_include_directories(pord PRIVATE ./include)

target_compile_options(pord PRIVATE -O)

target_compile_definitions(pord PRIVATE Add_=1 pord=1 MUMPS_ARITH=MUMPS_Arith_d)
