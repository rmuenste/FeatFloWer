set(Q2P1_BUILD_ID_FOUND false)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.1.0")
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED YES)
  else()
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED YES)
  endif()
else()
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED YES)
endif()

#===============================================================================================================
#                                              Intel builds
#===============================================================================================================

IF(Q2P1_BUILD_ID STREQUAL "nehalem-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -assume buffered_io -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "core2duo-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -assume buffered_io -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "i7-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -assume buffered_io -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "xeon-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "xeon-linux-intel-release-checks")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div -traceback -check -fpe0)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -fp-model precise -no-prec-div -fpp -traceback -check all,noarg_temp_created -fpe0)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "broadwell-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "xeongold-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xCORE-AVX2  -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xCORE-AVX2 -funroll-loops -assume underscore -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "opteron-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -assume buffered_io -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "opteronx2-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -assume buffered_io -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "phenomIIx4-linux-intel-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -xhost -funroll-loops -fp-model precise -no-prec-div)
  SET(Fortran_FLAGS -xhost -funroll-loops -assume underscore -assume buffered_io -fp-model precise -no-prec-div -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

#===============================================================================================================
#                                              GCC builds
#===============================================================================================================

IF(Q2P1_BUILD_ID STREQUAL "phenomIIx4-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "epyc16core-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "opteronx2-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "skylake-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "broadwell-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "nehalem-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "xeon-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "core2duo-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "i7-linux-gcc-release")
  SET(CMAKE_BUILD_TYPE "Release")
  SET(CXX_FLAGS_FC -march=native)
  SET(Fortran_FLAGS -march=native -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

#===============================================================================================================
#                                              Debug builds
#===============================================================================================================

IF(Q2P1_BUILD_ID STREQUAL "opteron-linux-intel-debug")
  SET(CMAKE_BUILD_TYPE "Debug")
  SET(CXX_FLAGS_FC)
  SET(Fortran_FLAGS -assume underscore -fpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "xeon-linux-intel-debug")
  SET(CMAKE_BUILD_TYPE "Debug")
  SET(CXX_FLAGS_FC -traceback -fpe0)
  SET(Fortran_FLAGS -assume underscore -fpp -traceback -check all,noarg_temp_created -fpe0)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "xeon-linux-gcc-debug")
  SET(CMAKE_BUILD_TYPE "Debug")
  SET(CXX_FLAGS_FC)
  SET(Fortran_FLAGS -finit-local-zero -ffixed-line-length-none -ffree-line-length-none -Wall -cpp)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

IF(Q2P1_BUILD_ID STREQUAL "opteronx2-linux-intel-debug")
  SET(CMAKE_BUILD_TYPE "Debug")
  SET(CXX_FLAGS_FC -traceback -check -fpe0)
  SET(Fortran_FLAGS -assume underscore -fpp -traceback -check all,noarg_temp_created -fpe0)
  SET(Q2P1_BUILD_ID_FOUND true)
ENDIF()

