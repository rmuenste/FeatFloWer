set(CPACK_GENERATOR "DEB;TGZ")

# Common package information
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "Extrud3D is a FEM-based CFD-Package for the simulation of extrusion processes.")

#SET(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://github.com/google/flatbuffers")

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Raphael Muenster <raphael.muenster@math.tu-dortmund.de>")

set(VERSION_MAJOR 0)
set(VERSION_MINOR 2)
set(VERSION_PATCH 3)

SET(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
SET(CPACK_PACKAGE_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CPACK_DEBIAN_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}")

# Derive architecture
IF(NOT CPACK_DEBIAN_PACKAGE_ARCHITECTURE)
  FIND_PROGRAM(DPKG_CMD dpkg)
  IF(NOT DPKG_CMD)
    MESSAGE(STATUS "Can not find dpkg in your path, default to i386.")
    SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
  ENDIF(NOT DPKG_CMD)
  EXECUTE_PROCESS(COMMAND "${DPKG_CMD}" --print-architecture
    OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
ENDIF(NOT CPACK_DEBIAN_PACKAGE_ARCHITECTURE)

# Package name
SET(CPACK_DEBIAN_PACKAGE_NAME "e3d")

# Add a license file
#SET(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/LICENSE.txt)
SET(CPACK_PACKAGE_FILE_NAME
        "${CPACK_DEBIAN_PACKAGE_NAME}_${CPACK_DEBIAN_PACKAGE_VERSION}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")

