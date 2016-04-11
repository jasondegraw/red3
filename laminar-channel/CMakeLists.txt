cmake_minimum_required(VERSION 2.8.11)

project(laminar-channel)

include_directories(../src)

add_executable(laminar-channel geom.cpp)
target_link_libraries(laminar-channel red3_core)

