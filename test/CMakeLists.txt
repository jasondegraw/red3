project(test)

enable_testing(true)

set(SRCS main.cpp
         #approxfractests.cpp
         arraytests.cpp
         #gridtests.cpp
         #staggeredgridtests.cpp
         #upwindtests.cpp
         #utilitiestests.cpp
)

include_directories(../src)
include_directories(../dependencies)

add_executable(red3tests ${SRCS})
add_test(red3tests red3tests)
#target_link_libraries(red3tests red3_core)
