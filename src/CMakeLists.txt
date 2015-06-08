cmake_minimum_required(VERSION 2.8.11)

project(lib)

#set(DEP_SRCS ../dependencies/jsoncpp/jsoncpp.cpp)

set(SRCS block.cpp
)

set(HDRS red3api.h
         block.hi)

add_library(red3_core SHARED ${HDRS} ${SRCS} ${DEP_SRCS})
