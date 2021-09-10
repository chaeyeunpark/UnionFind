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

option(BUILD_EXECUTABLES "Build C++ executables" OFF)

if(BUILD_EXECUTABLES)
	find_package(Eigen3 3.3 REQUIRED NO_MODULE)

	find_package(MPI)

	if(MPI_FOUND)
		message("MPI found. Build MPI executables.")
	endif()
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


if (BUILD_EXECUTABLES)
	add_library(toric_utils STATIC "src/error_utils.cpp" "src/toric_utils.cpp" "src/utility.cpp")

	add_executable(TestGenerateSyndromes "Test/TestGenerateSyndromes.cpp")
	target_link_libraries(TestGenerateSyndromes PRIVATE toric_utils)

	# 2D Bitflip noise
	add_executable(run_uf_2d_bitflip "src/run_uf_2d_bitflip.cpp")
	target_link_libraries(run_uf_2d_bitflip PRIVATE toric_utils Eigen3::Eigen)

	add_executable(run_uf_2d_bitflip_lazy "src/run_uf_2d_bitflip.cpp")
	target_compile_definitions(run_uf_2d_bitflip_lazy PUBLIC USE_LAZY)
	target_link_libraries(run_uf_2d_bitflip_lazy PRIVATE toric_utils Eigen3::Eigen)

	# 2D Depolarzing noise
	add_executable(run_uf_2d_depol "src/run_uf_2d_depol.cpp")
	target_link_libraries(run_uf_2d_depol PRIVATE toric_utils Eigen3::Eigen)

	add_executable(run_uf_2d_depol_lazy "src/run_uf_2d_depol.cpp")
	target_compile_definitions(run_uf_2d_depol_lazy PUBLIC USE_LAZY)
	target_link_libraries(run_uf_2d_depol_lazy PRIVATE toric_utils Eigen3::Eigen)

	# 3D Bitflip noise
	add_executable(run_uf_3d_bitflip "src/run_uf_3d_bitflip.cpp")
	target_link_libraries(run_uf_3d_bitflip PRIVATE toric_utils Eigen3::Eigen)

	add_executable(run_uf_3d_bitflip_lazy "src/run_uf_3d_bitflip.cpp")
	target_compile_definitions(run_uf_2d_bitflip_lazy PUBLIC USE_LAZY)
	target_link_libraries(run_uf_3d_bitflip_lazy PRIVATE toric_utils Eigen3::Eigen)

	# 3D Depolarizing noise
	add_executable(run_uf_3d_depol "src/run_uf_3d_depol.cpp")
	target_link_libraries(run_uf_3d_depol PRIVATE toric_utils Eigen3::Eigen)

	add_executable(run_uf_3d_depol_lazy "src/run_uf_3d_depol.cpp")
	target_compile_definitions(run_uf_2d_depol_lazy PUBLIC USE_LAZY)
	target_link_libraries(run_uf_3d_depol_lazy PRIVATE toric_utils Eigen3::Eigen)

	if (MPI_FOUND)
		# 2D Bitflip noise with MPI
		add_executable(run_uf_2d_bitflip_mpi "src/run_uf_2d_bitflip.cpp")
		target_compile_definitions(run_uf_2d_bitflip_mpi PUBLIC USE_MPI)
		target_link_libraries(run_uf_2d_bitflip_mpi PRIVATE toric_utils Eigen3::Eigen MPI::MPI_CXX)

		add_executable(run_uf_2d_bitflip_lazy_mpi "src/run_uf_2d_bitflip.cpp")
		target_compile_definitions(run_uf_2d_bitflip_lazy_mpi PUBLIC USE_MPI
			USE_LAZY)
		target_link_libraries(run_uf_2d_bitflip_lazy_mpi PRIVATE toric_utils Eigen3::Eigen MPI::MPI_CXX)

		# 2D Depolarzing noise with MPI
		add_executable(run_uf_2d_depol_mpi "src/run_uf_2d_depol.cpp")
		target_compile_definitions(run_uf_2d_depol_mpi PUBLIC USE_MPI)
		target_link_libraries(run_uf_2d_depol_mpi PRIVATE toric_utils Eigen3::Eigen MPI::MPI_CXX)

		add_executable(run_uf_2d_depol_lazy_mpi "src/run_uf_2d_depol.cpp")
		target_compile_definitions(run_uf_2d_depol_lazy_mpi PUBLIC USE_MPI
			USE_LAZY)
		target_link_libraries(run_uf_2d_depol_lazy_mpi PRIVATE toric_utils Eigen3::Eigen MPI::MPI_CXX)
	endif ()
endif()


pybind11_add_module(union_find "src/UnionFindPy.cpp" "src/utility.cpp")
