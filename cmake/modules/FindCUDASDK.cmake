message(STATUS "Looking for CUDA SDK")
# Look for the SDK stuff
find_path(CUDA_SDK_ROOT_DIR cutil.h
   PATH_SUFFIXES  "common/inc" "C/common/inc"
   "$ENV{NVSDKCUDA_ROOT}"
   "$ENV{NVSDKCOMPUTE_ROOT}"
   "[HKEY_LOCAL_MACHINE\\SOFTWARE\\NVIDIA Corporation\\Installed Products\\NVIDIA SDK 10\\Compute;InstallDir]"
   "/Developer/GPU\ Computing/C"
   )
 
# fallback method for determining CUDA_SDK_ROOT_DIR in case the previous one failed!
if (NOT CUDA_SDK_ROOT_DIR)
  message(FATAL_ERROR "The CUDA SDK COMPUTE_ROOT directory not found.")
else()
  message(STATUS "Found CUDA SDK COMPUTE_ROOT in ${CUDA_SDK_ROOT_DIR}")	
endif()

find_path(SHRUTIL_INCLUDE_DIR shrUtils.h
	${CUDA_SDK_ROOT_DIR}/../shared/inc
	DOC "The path to shrutil.h"
   )

if (NOT SHRUTIL_INCLUDE_DIR)
  message(FATAL_ERROR "The shared utilities were not found in the CUDA SDK.")
else()
  message(STATUS "Found shrUtils include in ${SHRUTIL_INCLUDE_DIR}")	
endif()

IF(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
	# Look for the library.
	FIND_LIBRARY(FREEGLUT_LIBRARY NAMES freeglut PATHS ${CUDA_SDK_ROOT_DIR}/../shared/lib/Win32)

	if (NOT FREEGLUT_LIBRARY)
		message(FATAL_ERROR "Freeglut library not found.")
	else()
		message(STATUS "Freeglut found int ${CUDA_SDK_ROOT_DIR}/../shared/lib/Win32")	
	endif()
ELSE(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
	# Look for the library.
	FIND_LIBRARY(GLUT_LIBRARY NAMES glut glut_x86_64 PATHS /usr/lib/)

	if (NOT GLUT_LIBRARY)
		message(FATAL_ERROR "glut library not found.")
	else()
		message(STATUS "glut found in /usr/lib/ ${GLUT_LIBRARY}")
	endif()
ENDIF(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")

IF(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
	# Look for the library.
	FIND_LIBRARY(GLEW_LIBRARY NAMES glew32 PATHS ${CUDA_SDK_ROOT_DIR}/../shared/lib/Win32)
	if (NOT GLEW_LIBRARY)
		message(FATAL_ERROR "glew32 library not found in the CUDA SDK.")
	else()
		message(STATUS "glew32 found int ${CUDA_SDK_ROOT_DIR}/../shared/lib/Win32")	
	endif()
ELSE(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
	# Look for the library.
	FIND_LIBRARY(GLEW_LIBRARY NAMES GLEW_x86_64 GLEW PATHS ${CUDA_SDK_ROOT_DIR}/../shared/lib/linux)
	if (NOT GLEW_LIBRARY)
		message(FATAL_ERROR "GLEW library not found in the CUDA SDK.")
	else()
		message(STATUS "GLEW found int ${CUDA_SDK_ROOT_DIR}/../shared/lib/linux ${GLEW_LIBRARY}")
	endif()
ENDIF(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")

IF(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
	FIND_LIBRARY(CUTIL_LIBRARY NAMES CUTIL cutil32 PATHS ${CUDA_SDK_ROOT_DIR}/common/lib/Win32
					DOC "The CUTIL library")
	if (NOT CUTIL_LIBRARY)
		message(FATAL_ERROR "cutil32 library not found in the CUDA SDK.")
	else()
		message(STATUS "cutil32 found int ${CUDA_SDK_ROOT_DIR}/common/lib/Win32")	
	endif()
ELSE(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
	FIND_LIBRARY(CUTIL_LIBRARY NAMES CUTIL cutil cutil_x86_64 PATHS ${CUDA_SDK_ROOT_DIR}/lib DOC "The CUTIL library")
	if (NOT CUTIL_LIBRARY)
		message(FATAL_ERROR "cutil library not found in ${CUDA_SDK_ROOT_DIR}/lib.")
	else()
		message(STATUS "cutil found in ${CUDA_SDK_ROOT_DIR}/lib")	
	endif()
ENDIF(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")

IF(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
	FIND_LIBRARY( SHRUTIL_LIBRARY
	NAMES SHRUTIL shrUtils32D
	PATHS
	${CUDA_SDK_ROOT_DIR}/../shared/lib/Win32
	DOC "The SHRUTIL library")
	if (NOT SHRUTIL_LIBRARY)
		message(FATAL_ERROR "shrUtils32 library not found in the CUDA SDK.")
	else()
		message(STATUS "shrUtils32 found int ${CUDA_SDK_ROOT_DIR}/../shared/lib/Win32")	
	endif()
ELSE(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
	FIND_LIBRARY( SHRUTIL_LIBRARY
	NAMES shrutil_x86_64 shrutil
	PATHS
	${CUDA_SDK_ROOT_DIR}/../shared/lib/
	DOC "The SHRUTIL library")
	if (NOT SHRUTIL_LIBRARY)
		message(FATAL_ERROR "shrutil library not found in the CUDA SDK.")
	else()
		message(STATUS "shrutil_x86_64 found in ${SHRUTIL_LIBRARY}")	
	endif()
ENDIF(NOT CMAKE_SYSTEM_NAME MATCHES "Linux")
