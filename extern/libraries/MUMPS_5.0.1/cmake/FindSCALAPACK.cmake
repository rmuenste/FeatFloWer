# - Try to find SCALAPACK
#

if (SCALAPACK_INCLUDES AND SCALAPACK_LIBRARIES)
  set(SCALAPACK_FIND_QUIETLY TRUE)
endif (SCALAPACK_INCLUDES AND SCALAPACK_LIBRARIES)

find_path (SCALAPACK_INCLUDES NAMES lib/libscalapack.a PATHS $ENV{HOME}/scalapack DOC "ScaLapack Directory")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCALAPACK DEFAULT_MSG
                                  SCALAPACK_INCLUDES SCALAPACK_LIBRARIES)

mark_as_advanced(SCALAPACK_INCLUDES SCALAPACK_LIBRARIES)
