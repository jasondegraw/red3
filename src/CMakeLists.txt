cmake_minimum_required(VERSION 2.8.11)

project(lib)

set(DEP_SRCS ../dependencies/jsoncpp/jsoncpp.cpp)
include_directories(../dependencies/jsoncpp)

set(SRCS array.cpp
         grid.cpp
         util.cpp
         staggeredgrid.cpp
         upwind.cpp
         util.cpp)

set(HDRS red3api.hpp
         array.hpp
         arrayops.hpp
         def.hpp
         grid.hpp
         red3.hpp
         staggeredgrid.hpp
         upwind.hpp
         util.hpp)

add_library(red3_core SHARED ${HDRS} ${SRCS} ${DEP_SRCS})

