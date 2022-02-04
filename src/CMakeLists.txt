cmake_minimum_required(VERSION 2.8.11)

project(lib)

#set(DEP_SRCS ../dependencies/jsoncpp/jsoncpp.cpp)
#include_directories(../dependencies/jsoncpp) 

set(SRCS #
         approxfrac.cpp
         grid.cpp
         #util.cpp
         staggeredgrid.cpp
         #upwind.cpp
         utilities.cpp
)

set(HDRS ../include/red3/defs.hpp
         ../include/red3/approxfrac.hpp
         ../include/red3/array.hpp
         ../include/red3/arrayops.hpp
         #def.hpp
         ../include/red3/grid.hpp
         #red3.hpp
         ../include/red3/staggeredgrid.hpp
         #upwind.hpp
         ../include/red3/utilities.hpp
)

add_library(red3_core SHARED ${HDRS} ${SRCS} ${DEP_SRCS})

