cmake_minimum_required(VERSION 2.8.11)

project(lib)

set(DEP_SRCS ../dependencies/jsoncpp/jsoncpp.cpp)
include_directories(../dependencies/jsoncpp)

set(SRCS block.cpp
         geom.cpp
         grid.cpp
         util.cpp
         main.cpp
         staggeredgrid.cpp
         util.cpp)

set(HDRS red3api.hpp
         arrayops.hpp
         block.hpp
         def.hpp
         geom.hpp
         red3.hpp
         staggeredgrid.hpp
         util.hpp)

add_library(red3_core SHARED ${HDRS} ${SRCS} ${DEP_SRCS})
