cmake_minimum_required(VERSION 3.10)
project(c__)

set(CMAKE_C_COMPILER /usr/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/bin/mpic++)

set(CMAKE_CXX_FLAGS "-g -pg -O3 -ftree-vectorize -funroll-loops")
set(CMAKE_CXX_STANDARD 11)

add_executable(c__ main.cpp)