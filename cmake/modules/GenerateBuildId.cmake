# 
#
# 
# 
# 
# 
IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
  IF(NOT Q2P1_BUILD_ID)
    # A build id is composed of: 
    # ${cpu_type}-${os}-${compiler}-${build_type}
    # Example: nehalem-linux-intel-release

    # We do not have a user defined build id
    # so we select a default build
    set(CMAKE_BUILD_TYPE "Release")
    message(STATUS "No build id selected... trying to determine machine type and compiler settings")
    set(_vendor_id)
    set(_cpu_family)
    set(_cpu_model)
    set(_cpu_flags)

    set(Q2P1_CPU_TYPE)
    set(Q2P1_OS)
    set(Q2P1_DEFAULT_BUILD "release")

      set(Q2P1_OS "linux") 
      file(READ "/proc/cpuinfo" _cpuinfo)
      string(REGEX REPLACE ".*vendor_id[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _vendor_id "${_cpuinfo}")
      string(REGEX REPLACE ".*cpu family[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_family "${_cpuinfo}")
      string(REGEX REPLACE ".*model[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_model "${_cpuinfo}")
      string(REGEX REPLACE ".*flags[ \t]*:[ \t]+([ \ta-zA-Z0-9_-]+)(.*)" "\\1" _cpu_flags "${_cpuinfo}")

    message(STATUS "vendor_id ${_vendor_id}")
    message(STATUS "family ${_cpu_family}")
    message(STATUS "model ${_cpu_model}")
    message(STATUS "flags ${_cpu_flags}")

    IF(${_vendor_id} MATCHES "AuthenticAMD")
      IF(${_cpu_family} EQUAL 16)
        IF(${_cpu_model} EQUAL 4)
          set(Q2P1_CPU_TYPE "phenomIIx4")
        ELSEIF(${_cpu_model} EQUAL 2) 
          set(Q2P1_CPU_TYPE "opteron")      
        ENDIF()
      ELSEIF(${_cpu_family} EQUAL 15)
        IF(${_cpu_model} EQUAL 33)
          set(Q2P1_CPU_TYPE "opteron")
        ELSEIF(${_cpu_model} EQUAL 65 OR ${_cpu_model} EQUAL 66 )
          set(Q2P1_CPU_TYPE "opteronx2")
        ENDIF()
      ELSEIF(${_cpu_family} EQUAL 23)
        IF(${_cpu_model} EQUAL 1)
          set(Q2P1_CPU_TYPE "epyc16core")
        ENDIF()
        IF(${_cpu_model} EQUAL 49)
          set(Q2P1_CPU_TYPE "epyc32core")
        ENDIF()
      ELSEIF(${_cpu_family} EQUAL 25)
        IF(${_cpu_model} EQUAL 80)
          set(Q2P1_CPU_TYPE "zen3")
        ENDIF()
      ELSE()
        message(FATAL_ERROR "Unknown CPU model, cannot set up default configuration")
      ENDIF()    
    ELSEIF(${_vendor_id} MATCHES "GenuineIntel")
      IF(${_cpu_family} EQUAL 6)
        IF("15 21 22 23 29" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "core2duo")
        ELSEIF("26 30 37 44 45 46" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "nehalem")
        ELSEIF("42" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "i7")
        ELSEIF("58 60" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "xeon")
        ELSEIF("63" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "sandy")
        ELSEIF("62" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "haswell")
        ELSEIF("94" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "skylake")
        ELSEIF("85" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "xeongold")
        ELSEIF("79" MATCHES ${_cpu_model})
          set(Q2P1_CPU_TYPE "broadwell")
        ELSE()
          message(FATAL_ERROR "Unknown CPU model, cannot set up default configuration")
        ENDIF()
      ENDIF()
    ELSE()
      message(FATAL_ERROR "Unknown CPU vendor, cannot set up default configuration")
    ENDIF()

    message(STATUS "Processor ${CMAKE_SYSTEM_PROCESSOR}")
    message(STATUS "CMAKE_Fortran_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID}")
    message(STATUS "CMAKE_CXX_COMPILER_ID : " ${CMAKE_CXX_COMPILER_ID} )
    message(STATUS "CMAKE_CXX_COMPILER_VERSION : " ${CMAKE_CXX_COMPILER_VERSION} )

    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(Q2P1_DEFAULT_COMPILER "gcc")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      set(Q2P1_DEFAULT_COMPILER "intel")
    endif()

    set(Q2P1_BUILD_ID "${Q2P1_CPU_TYPE}-${Q2P1_OS}-${Q2P1_DEFAULT_COMPILER}-${Q2P1_DEFAULT_BUILD}")

    include(${CMAKE_MODULE_PATH}/SetFlagsForID.cmake)  

    IF(NOT Q2P1_BUILD_ID_FOUND)
      message(FATAL_ERROR "Build id:<${Q2P1_BUILD_ID}> was not found.")
    ENDIF(NOT Q2P1_BUILD_ID_FOUND)

    message(STATUS "Configuring for build id:<${Q2P1_BUILD_ID}>")


  ELSE(NOT Q2P1_BUILD_ID)
    # set the build id and search for it
    set(Q2P1_BUILD_ID_USER ${Q2P1_BUILD_ID})
    include(${CMAKE_MODULE_PATH}/SetFlagsForID.cmake)  
    IF(NOT Q2P1_BUILD_ID_FOUND)
      message(FATAL_ERROR "User requested build id:<${Q2P1_BUILD_ID_USER}> was not found.")
    ENDIF(NOT Q2P1_BUILD_ID_FOUND)
    message(STATUS "Configuring for build id:<${Q2P1_BUILD_ID}>")
  ENDIF(NOT Q2P1_BUILD_ID)

else(CMAKE_SYSTEM_NAME MATCHES "Linux")
  IF(WIN32)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /heap-arrays0")
    message(STATUS "Set windows flags to ${CMAKE_Fortran_FLAGS}")
  ENDIF(WIN32)
ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")

