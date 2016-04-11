project(test)

enable_testing(true)

set(SRCS arraytests.cpp
         staggeredgridtests.cpp
)

include_directories(../src ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(red3tests ${SRCS})
add_test(red3tests red3tests)
target_link_libraries(red3tests red3_core gtest gtest_main)
