cmake_minimum_required(VERSION 3.9)

project(UnionFindCPP)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/externals)
include_directories(${CMAKE_SOURCE_DIR}/robin-map/include)

add_subdirectory(pybind11)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(LINK_WHAT_YOU_USE TRUE)
set(BUILD_RPATH_USE_ORIGIN TRUE)

option(BUILD_EXECUTABLE "Build C++ executables" OFF)

if (BUILD_EXECUTABLE)
	find_package (Eigen3 3.3 REQUIRED NO_MODULE)
endif()

set(JSON_URL "https://github.com/nlohmann/json/releases/download/v3.7.3/json.hpp")
set(JSON_DOWNLOAD_DIR "${CMAKE_SOURCE_DIR}/externals/nlohmann")
set(JSON_DOWNLOAD_PATH "${JSON_DOWNLOAD_DIR}/json.hpp")

set(CATCH_URL "https://github.com/catchorg/Catch2/releases/download/v2.11.1/catch.hpp")
set(CATCH_DOWNLOAD_PATH "${CMAKE_SOURCE_DIR}/externals/catch.hpp")


if (NOT EXISTS "${JSON_DOWNLOAD_PATH}")
	file(MAKE_DIRECTORY "${JSON_DOWNLOAD_DIR}")
	file(DOWNLOAD "${JSON_URL}" "${JSON_DOWNLOAD_PATH}")
endif()

if (NOT EXISTS "${CATCH_DOWNLOAD_PATH}")
	file(DOWNLOAD "${CATCH_URL}" "${CATCH_DOWNLOAD_PATH}")
endif()


if (BUILD_EXECUTABLE)
	add_executable(TestGenerateSyndromes "Test/TestGenerateSyndromes.cpp" "src/utility.cpp" "src/cpp_utils.cpp")

	add_executable(RunUnionFind2D "src/RunUnionFind2D.cpp" "src/utility.cpp" "src/cpp_utils.cpp")
	target_link_libraries(RunUnionFind2D Eigen3::Eigen)

	add_executable(RunUnionFind2D_bitflip "src/RunUnionFind2D_bitflip.cpp" "src/utility.cpp" "src/cpp_utils.cpp")
	target_link_libraries(RunUnionFind2D_bitflip Eigen3::Eigen)

	add_executable(RunUnionFind3D "src/RunUnionFind3D.cpp" "src/utility.cpp" "src/cpp_utils.cpp")
	target_link_libraries(RunUnionFind3D Eigen3::Eigen)
endif()

pybind11_add_module(union_find "src/UnionFindPy.cpp" "src/utility.cpp")
