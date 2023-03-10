# find lib2ds includes and library
#
# LIB2DS_INCLUDE_DIR - where the lib2ds directory containing the headers can be
#                      found
# LIB2DS_LIBRARY     - full path to the lib2ds library
# LIB2DS_FOUND       - TRUE if lib2ds was found

FIND_LIBRARY(FEATFLOW_FEAT2D_LIBRARY NAMES feat2d
             PATHS  $ENV{FEATFLOW1_LIBDIR} )


IF(FEATFLOW_FEAT2D_LIBRARY)
    MESSAGE(STATUS "Found libfeat2d library: ${FEATFLOW_FEAT2D_LIBRARY}")
ELSE(FEATFLOW_FEAT2D_LIBRARY)
    MESSAGE(STATUS "Could NOT find libfeat2d library.")
ENDIF(FEATFLOW_FEAT2D_LIBRARY)


IF(FEATFLOW_FEAT2D_LIBRARY)
    SET(FEATFLOW_FEAT2D_FOUND TRUE)
ELSE(FEATFLOW_FEAT2D_LIBRARY)
    SET(FEATFLOW_FEAT2D_FOUND FALSE)
    IF(FEATFLOW_FEAT2D_LIBRARY_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Featflow libfeat2d not found")
    ENDIF(FEATFLOW_FEAT2D_LIBRARY_FIND_REQUIRED)
ENDIF(FEATFLOW_FEAT2D_LIBRARY)

