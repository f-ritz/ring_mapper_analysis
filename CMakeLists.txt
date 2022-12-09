cmake_minimum_required(VERSION 3.22)
project(ring_mapper_analysis)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wno-c++11-extensions -std=gnu++1y -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O4")

include_directories(/Users/jyesselm/CLionProjects/ring_mapper_analysis/src/)
include_directories(/Users/jyesselm/CLionProjects/ring_mapper_analysis/src/base)
include_directories(/Users/jyesselm/CLionProjects/ring_mapper_analysis/src/unittests)
include_directories(/Users/jyesselm/CLionProjects/ring_mapper_analysis/src/unittests/all_unittests)
include_directories(/Users/jyesselm/CLionProjects/ring_mapper_analysis/src/unittest_resources)

add_executable(ring_mapper_analysis src/main.cpp src/base/functions.h src/base/ThreeDInfoVector.cpp src/base/ThreeDInfoVector.h)
add_executable(ring_mapper_unittest src/unittests/all_unittests/ring_mapper_unittest.cpp src/unittests/doctest.h src/unittests/common.hpp)

add_test(ring_mapper_unittest ring_mapper_unittest)
