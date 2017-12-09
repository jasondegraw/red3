cmake_minimum_required(VERSION 2.8.11)

project(simple2d)

include_directories(../src)

add_executable(simple2d geom.cpp)
target_link_libraries(simple2d red3_core)

