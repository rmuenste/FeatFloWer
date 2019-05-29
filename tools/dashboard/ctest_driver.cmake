#=========================================================================
# This is the CMake driver script for the FeatFloWer Test Suite
# 
# The Test Suite is intented to be run under a Linux system.
# In order to run this CMake script from a command line
# a syntax of the following kind is needed:
#
# ctest -S ctest_driver.cmake [parameters]
#
# The required parameters to the script are:
# [-DSRC_DIR=]: the source directory 
# [-DBIN_DIR=]: the binary directory 
# 
# The optional parameters to the script are:
# [-DBUILD_STRING=]: the build string, i.e.: xeon-linux-gcc-release
# [-DCLEAN_BIN=]: Set this to <True> to clean the bin dir before building
#
# Example start command:
# ctest -S ctest_driver.cmake \
# -DBUILD_STRING=xeon-linux-gcc-release \
# -DCONSTRAINT=epyx \
# -DPROCS=16 \
# -DSRC_DIR=$(pwd)/Feat_FloWer \
# -DBIN_DIR=$(pwd)/bin-bttf \
# --verbose
#
#=========================================================================
# allow less strict IF-else syntax
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

#========================================================================
# Save the input parameters in a CMake variable
#========================================================================
set(CTEST_SOURCE_DIRECTORY ${SRC_DIR})
set(CTEST_BINARY_DIRECTORY ${BIN_DIR})
set(BENCH_CLEAN_BIN ${CLEAN_BIN})
set(BENCH_PROCS ${PROCS})

set(USER_BUILD_NAME ${BUILD_STRING})
set(SLURM_CONSTRAINT ${CONSTRAINT})

if("${BENCH_PROCS}" STREQUAL "")
  message(FATAL_ERROR "Required parameter PROCS is not set. Set it using the -DPROCS=<number> option.")
else("${BENCH_PROCS}" STREQUAL "")
  message(STATUS "Using ${PROCS} processors for the benchmark.")
endif("${BENCH_PROCS}" STREQUAL "")

message(STATUS "${CTEST_SOURCE_DIRECTORY} ${CTEST_BINARY_DIRECTORY}")

cmake_host_system_information(RESULT HNAME QUERY HOSTNAME)

#========================================================================
# Set the maximum wall time for the Test Suite
#========================================================================
set(CTEST_TEST_TIMEOUT 7200) 


#========================================================================
# Here we configure the CTest variables that
# are required to submit a CTest result to a CDash dashboard server
#========================================================================
set(CTEST_SITE ${HNAME})

set(CTEST_DROP_METHOD "http")

set(CTEST_BUILD_NAME "${BUILD_STRING}-${CONSTRAINT}")

set(CTEST_DROP_SITE "localhost/CDash/public")

set(CTEST_DROP_LOCATION "/submit.php?project=Feat_FloWer")

#========================================================================
# Configuration of the underlying CMake generator
#========================================================================
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

set(CTEST_BUILD_CONFIGURATION "Profiling")

set(CTEST_BUILD_OPTIONS "-DWITH_SSH1=ON -WITH_SFTP=ON -DWITH_SERVER=ON -DWITH_ZLIB=ON -DWITH_PCAP=ON -DWITH_GCRYPT=OFF")

include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CTEST_BUILD_FLAGS -j7)
  set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()

message(STATUS "Drop site: ${CTEST_DROP_SITE}")
message(STATUS "Hostname: ${HNAME}")
message(STATUS "Build flags: ${CTEST_BUILD_FLAGS}")

#========================================================================
# Configure memory checking
#========================================================================
set(WITH_MEMCHECK false)
set(WITH_COVERAGE false)

#========================================================================
# How should the binary directory be treated
#========================================================================
set(CLEAN_BIN_DIR "")

if(CLEAN_BIN)
  message(STATUS "Preparing to delete the binary directory.")
  if(EXISTS ${CTEST_BINARY_DIRECTORY})
    message(STATUS "Deleting the binary directory.")
    ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
  endif()
else(CLEAN_BIN)
  message(STATUS "Using the current contents of the binary directory")
endif(CLEAN_BIN)

#========================================================================
# Configure the git checkout command
#========================================================================
find_program(CTEST_GIT_COMMAND NAMES git)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)

#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/tests/valgrind.supp)

#========================================================================
# If the source directory does not exist, get it from the repo
#========================================================================
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  message(FATAL_ERROR "CTest source not found: ${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone --recursive ssh://rmuenste@arryn.mathematik.tu-dortmund.de/home/user/git/Feat_FloWer.git")
endif()

set(CTEST_UPDATE_COMMAND "git")

set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION} -DPROCS=${BENCH_PROCS}")

set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTING:BOOL=ON ${CTEST_BUILD_OPTIONS}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

#========================================================================
# Call the CTest routines for invoking the particular tests
#========================================================================
#set(CTEST_SOURCE_DIRECTORY "/home/user/rmuenste/nobackup/code/ode_intermediate/Feat_Flower")
#set(CTEST_SOURCE_DIRECTORY "/data/warehouse14/rmuenste/code/ode_intermediate/Feat_Flower")
#message(FATAL_ERROR "${CTEST_SOURCE_DIRECTORY}")

ctest_start("Experimental")

#ctest_update()

ctest_configure()

ctest_build()

ctest_test(START 1 END 3 STRIDE 1)

#ctest_test()

if (WITH_COVERAGE AND CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif (WITH_COVERAGE AND CTEST_COVERAGE_COMMAND)
if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)

#========================================================================
# Set the name of the note file
# The note file is a json-file that stores the simulation
# output in json format, so it can be proccessed in the 
# CDash website.
#========================================================================
file(GLOB FF_NOTES ${CTEST_BINARY_DIRECTORY}/*-bench.json)
set(CTEST_NOTES_FILES ${FF_NOTES})

ctest_submit()
