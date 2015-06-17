cmake_minimum_required(VERSION 2.8.11)

project(lib)

set(DEP_SRCS ../dependencies/jsoncpp/jsoncpp.cpp)
include_directories(../dependencies/jsoncpp)

set(SRCS block.cpp
         geom.cpp
         grid.cpp
         util.cpp
         main.cpp
         staggeredmodel.cpp
         util.cpp)

set(HDRS red3api.h
         block.h
         def.h
         geom.h
         red3.h
         staggeredmodel.h
         util.h)

add_library(red3_core SHARED ${HDRS} ${SRCS} ${DEP_SRCS})
