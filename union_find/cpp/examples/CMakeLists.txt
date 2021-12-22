project(UnionFindCPP_examples)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(MPI)

if(MPI_FOUND)
	message("MPI found. Build MPI executables.")
endif()

add_library(toric_utils STATIC "src/error_utils.cpp" "src/toric_utils.cpp" "src/utility.cpp")
target_link_libraries(toric_utils PRIVATE Eigen3::Eigen)

add_executable(TestGenerateSyndromes "Test/TestGenerateSyndromes.cpp")
target_link_libraries(TestGenerateSyndromes PRIVATE toric_utils Eigen3::Eigen)

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

	# 3D Bitflip noise with MPI
	add_executable(run_uf_3d_bitflip_mpi "src/run_uf_3d_bitflip.cpp")
	target_compile_definitions(run_uf_3d_bitflip_mpi PUBLIC USE_MPI)
	target_link_libraries(run_uf_3d_bitflip_mpi PRIVATE toric_utils Eigen3::Eigen MPI::MPI_CXX)

	add_executable(run_uf_3d_bitflip_lazy_mpi "src/run_uf_3d_bitflip.cpp")
	target_compile_definitions(run_uf_3d_bitflip_lazy_mpi PUBLIC USE_MPI
		USE_LAZY)
	target_link_libraries(run_uf_3d_bitflip_lazy_mpi PRIVATE toric_utils Eigen3::Eigen MPI::MPI_CXX)


	# 3D Depolarzing noise with MPI
	add_executable(run_uf_3d_depol_mpi "src/run_uf_3d_depol.cpp")
	target_compile_definitions(run_uf_3d_depol_mpi PUBLIC USE_MPI)
	target_link_libraries(run_uf_3d_depol_mpi PRIVATE toric_utils Eigen3::Eigen MPI::MPI_CXX)

	add_executable(run_uf_3d_depol_lazy_mpi "src/run_uf_3d_depol.cpp")
	target_compile_definitions(run_uf_3d_depol_lazy_mpi PUBLIC USE_MPI
		USE_LAZY)
	target_link_libraries(run_uf_3d_depol_lazy_mpi PRIVATE toric_utils Eigen3::Eigen MPI::MPI_CXX)
endif ()
