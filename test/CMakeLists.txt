project(test)

enable_testing(true)

set(SRCS approxfractests.cpp
         arraytests.cpp
         staggeredgridtests.cpp
         upwindtests.cpp
         utilitiestests.cpp
)

include_directories(../src) # ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(red3tests catch.hpp ${SRCS})
add_test(red3tests red3tests)
target_link_libraries(red3tests red3_core)
