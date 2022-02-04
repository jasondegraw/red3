cmake_minimum_required(VERSION 2.8.11)

project(red3)

if(CMAKE_COMPILER_IS_GNUCXX)
  # C++ 20 is allegedly available beginning in GCC v8
  # https://gcc.gnu.org/projects/cxx-status.html#cxx20
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "8.0.0")
    message(FATAL_ERROR "g++ versions earlier than 8.0.0 are not supported")
  endif()
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "10.0.0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
  endif()
elseif(MSVC)
  # https://devblogs.microsoft.com/cppblog/msvc-cpp20-and-the-std-cpp20-switch/
  if(${CMAKE_C_COMPILER_VERSION} VERSION_LESS "19.29.30040.0")
    message(FATAL_ERROR "Visual Studio earlier than VS2019 16.11 is not supported")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest") # Switch to latest to get format, need to figure out what that means for the version
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Need to check for clang version
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# JSON Dependency
#include_directories("dependencies/jsoncpp")
# Eigen Dependency
#include_directories(dependencies/eigen-eigen-5a0156e40feb/)
# fmt
#add_subdirectory(dependencies/fmt-6.0.0)
#include_directories(dependencies/fmt-6.0.0/include)

include_directories(include)

add_subdirectory(src)
add_subdirectory(laminar-channel)
#add_subdirectory(simple2d)
#add_subdirectory(doc)

add_subdirectory(test)

