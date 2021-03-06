cmake_minimum_required(VERSION 2.8.11)

project(red3)

if(CMAKE_COMPILER_IS_GNUCXX)
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.6.0")
    message(FATAL_ERROR "g++ versions earlier than 4.6.0 are not supported")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
elseif(MSVC)
  if(${CMAKE_C_COMPILER_VERSION} VERSION_LESS "15.2.26430.04")
    message(FATAL_ERROR "Visual Studio earlier than VS2017 is not supported")
  endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -stdlib=libc++")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# JSON Dependency
#include_directories("dependencies/jsoncpp")
# Eigen Dependency
include_directories(dependencies/eigen-eigen-5a0156e40feb/)

add_subdirectory(src)
add_subdirectory(laminar-channel)
add_subdirectory(simple2d)
#add_subdirectory(doc)

add_subdirectory(test)

