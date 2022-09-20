cmake_minimum_required(VERSION 3.22)
project(ring_mapper_analysis)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wno-c++11-extensions -std=gnu++1y -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O4")

add_executable(ring_mapper_analysis main.cpp functions.h)
include_directories(/Users/jyesselm/CLionProjects/ring_mapper_analysis)
