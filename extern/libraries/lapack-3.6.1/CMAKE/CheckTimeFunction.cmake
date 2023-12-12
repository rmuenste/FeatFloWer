# - Check if the Fortran function exists.
# CHECK_TIME_FUNCTION(FUNCTION VARIABLE TYPE)
# - macro which checks if the Fortran function exists
#  FUNCTION - the name of the Fortran function
#  VARIABLE - variable to store the result
#

macro(CHECK_TIME_FUNCTION FUNCTION VARIABLE)

    set(OUTPUT "")  # Initialize OUTPUT_VARIABLE to an empty string

    try_compile(RES
    ${PROJECT_BINARY_DIR}/INSTALL
    ${PROJECT_SOURCE_DIR}/INSTALL
    TIMING secondtst_${FUNCTION}
    OUTPUT_VARIABLE OUTPUT)

    if(RES)
      #set(${VARIABLE} ${FUNCTION} CACHE INTERNAL "Have Fortran function ${FUNCTION}")
      set(${VARIABLE} ${FUNCTION})
      message(STATUS "Looking for Fortran ${FUNCTION} - found - function name: ${VARIABLE} ${FUNCTION}")
      file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
        "Fortran ${FUNCTION} exists. ${OUTPUT} \n\n")
      message(STATUS "Fortran ${FUNCTION} exists. ${OUTPUT}")
    else(RES)
      message(STATUS "Looking for Fortran ${FUNCTION} - not found")
      file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Fortran ${FUNCTION} does not exist. \n ${OUTPUT} \n")
      message(STATUS "Fortran ${FUNCTION} does not exist. ${OUTPUT}")
    endif(RES)
endmacro(CHECK_TIME_FUNCTION)


