# find lib3ds includes and library
#
# LIB3DS_INCLUDE_DIR - where the lib3ds directory containing the headers can be
#                      found
# LIB3DS_LIBRARY     - full path to the lib3ds library
# LIB3DS_FOUND       - TRUE if lib3ds was found

FIND_LIBRARY(FEATFLOW_FEAT3D_LIBRARY NAMES libfeat3d feat3d
             PATHS $ENV{FEATFLOW1_LIBDIR} )


IF(FEATFLOW_FEAT3D_LIBRARY)
    MESSAGE(STATUS "Found libfeat3d library: ${FEATFLOW_FEAT3D_LIBRARY}")
ELSE(FEATFLOW_FEAT3D_LIBRARY)
    MESSAGE(STATUS "libfeat3d not found in: $ENV{FEATFLOW1_LIBDIR}")
    MESSAGE(STATUS "Could NOT find libfeat3d library.")
ENDIF(FEATFLOW_FEAT3D_LIBRARY)


IF(FEATFLOW_FEAT3D_LIBRARY)
    SET(FEATFLOW_FEAT3D_FOUND TRUE)
ELSE(FEATFLOW_FEAT3D_LIBRARY)
    SET(FEATFLOW_FEAT3D_FOUND FALSE)
    IF(FEATFLOW_FEAT3D_LIBRARY_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Featflow libfeat3d not found")
    ENDIF(FEATFLOW_FEAT3D_LIBRARY_FIND_REQUIRED)
ENDIF(FEATFLOW_FEAT3D_LIBRARY)

# vim: et sw=4 ts=4
