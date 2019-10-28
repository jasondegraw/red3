cmake_minimum_required(VERSION 2.8.11)

project(lib)

#set(DEP_SRCS ../dependencies/jsoncpp/jsoncpp.cpp)
#include_directories(../dependencies/jsoncpp) 

set(SRCS array.cpp
         approxfrac.cpp
         grid.cpp
         util.cpp
         staggeredgrid.cpp
         upwind.cpp
         util.cpp
         utilities.cpp)

set(HDRS red3api.hpp
         approxfrac.hpp
         array.hpp
         arrayops.hpp
         def.hpp
         grid.hpp
         red3.hpp
         staggeredgrid.hpp
         upwind.hpp
         util.hpp
         utilities.hpp)

add_library(red3_core SHARED ${HDRS} ${SRCS} ${DEP_SRCS})
target_link_libraries(red3_core fmt::fmt)

