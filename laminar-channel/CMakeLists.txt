cmake_minimum_required(VERSION 2.8.11)

project(laminar-channel)

include_directories(../include)

add_executable("entrance-length" "entrance-length.cpp")
target_link_libraries("entrance-length" red3_core)

