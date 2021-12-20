cmake_minimum_required(VERSION 3.14)

project(UnionFindCPP)


if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)
set(LINK_WHAT_YOU_USE TRUE)
set(BUILD_RPATH_USE_ORIGIN TRUE)

option(BUILD_EXECUTABLES "Build C++ executables" OFF)
option(BUILD_TESTS "Build C++ tests" OFF)


add_subdirectory(union_find/src)

if(BUILD_TESTS)
	enable_testing()
endif()


if (BUILD_EXECUTABLES)
	add_subdirectory(examples)
endif()

