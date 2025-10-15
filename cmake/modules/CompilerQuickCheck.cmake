# What happens here:
# attempts to set the C++ standard to C++17.​
#
# ensures that CMake will produce an error if the requested standard is not supported.​
#
# disables compiler-specific extensions, promoting portability.​
#
# The check_cxx_compiler_flag function checks if the compiler supports the -std=c++17 flag. If not, the script falls back to C++14 and issues a warning.

#==============================================================================
# Set the desired C++ standard to 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Check if the compiler supports C++17
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)

if(NOT COMPILER_SUPPORTS_CXX17)
  # Fallback to C++14 if C++17 is not supported
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED YES)
  message(WARNING "Compiler does not support C++17; falling back to C++14.")
endif()
#==============================================================================

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # Additional settings specific to GCC, if necessary
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # Additional settings specific to Intel compiler, if necessary
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # Additional settings specific to Microsoft Visual Studio compiler, if necessary
endif()