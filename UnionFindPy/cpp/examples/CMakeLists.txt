project(UnionFindCPP_examples)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(fmt EXCLUDE_FROM_ALL)

find_package(MPI)
if(MPI_FOUND)
	message("MPI found. Build MPI executables.")
endif()

add_library(example_utils STATIC "error_utils.cpp" "toric_utils.cpp" "runner_utils.cpp")
target_link_libraries(example_utils PRIVATE union_find_cpp_dependency Eigen3::Eigen)
target_link_libraries(example_utils PRIVATE fmt::fmt-header-only)

# 2D Bitflip noise
add_executable(run_uf_2d_bitflip "run_uf_2d_bitflip.cpp")
target_link_libraries(run_uf_2d_bitflip PRIVATE example_utils union_find_cpp_dependency Eigen3::Eigen)

add_executable(run_uf_2d_bitflip_lazy "run_uf_2d_bitflip.cpp")
target_compile_definitions(run_uf_2d_bitflip_lazy PUBLIC USE_LAZY)
target_link_libraries(run_uf_2d_bitflip_lazy PRIVATE example_utils union_find_cpp_dependency Eigen3::Eigen)

# 3D Bitflip noise
add_executable(run_uf_3d_bitflip "run_uf_3d_bitflip.cpp")
target_link_libraries(run_uf_3d_bitflip PRIVATE example_utils union_find_cpp_dependency Eigen3::Eigen)

add_executable(run_uf_3d_bitflip_lazy "run_uf_3d_bitflip.cpp")
target_compile_definitions(run_uf_2d_bitflip_lazy PUBLIC USE_LAZY)
target_link_libraries(run_uf_3d_bitflip_lazy PRIVATE example_utils union_find_cpp_dependency Eigen3::Eigen)

if (MPI_FOUND)
	# 2D Bitflip noise with MPI
	add_executable(run_uf_2d_bitflip_mpi "run_uf_2d_bitflip.cpp")
	target_compile_definitions(run_uf_2d_bitflip_mpi PUBLIC USE_MPI)
	target_link_libraries(run_uf_2d_bitflip_mpi PRIVATE example_utils union_find_cpp_dependency Eigen3::Eigen MPI::MPI_CXX)

	add_executable(run_uf_2d_bitflip_lazy_mpi "run_uf_2d_bitflip.cpp")
	target_compile_definitions(run_uf_2d_bitflip_lazy_mpi PUBLIC USE_MPI
		USE_LAZY)
	target_link_libraries(run_uf_2d_bitflip_lazy_mpi PRIVATE example_utils union_find_cpp_dependency Eigen3::Eigen MPI::MPI_CXX)

	# 3D Bitflip noise with MPI
	add_executable(run_uf_3d_bitflip_mpi "run_uf_3d_bitflip.cpp")
	target_compile_definitions(run_uf_3d_bitflip_mpi PUBLIC USE_MPI)
	target_link_libraries(run_uf_3d_bitflip_mpi PRIVATE example_utils union_find_cpp_dependency Eigen3::Eigen MPI::MPI_CXX)

	add_executable(run_uf_3d_bitflip_lazy_mpi "run_uf_3d_bitflip.cpp")
	target_compile_definitions(run_uf_3d_bitflip_lazy_mpi PUBLIC USE_MPI
		USE_LAZY)
	target_link_libraries(run_uf_3d_bitflip_lazy_mpi PRIVATE example_utils union_find_cpp_dependency Eigen3::Eigen MPI::MPI_CXX)
endif ()
