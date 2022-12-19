cmake_minimum_required(VERSION 3.6)
project(task3)
include_directories(/usr/include/local)
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(hyperbolic_equation hyperbolic_equation.cpp hyperbolic_equation.h)
