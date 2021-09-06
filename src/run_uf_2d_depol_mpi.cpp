#include <random>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <atomic>
#include <chrono>

#include <Eigen/Dense>

#include <sstream>

#include <mpi.h>

#include "error_utils.hpp"
#include "Lattice2D.hpp"
#include "UnionFind.hpp"

int main(int argc, char* argv[])
{
	namespace chrono = std::chrono;

	int mpi_rank, mpi_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	std::random_device rd;
	std::default_random_engine re{rd()};

	auto noise_type = NoiseType::Depolarizing;

	if(argc != 3)
	{
		printf("Usage: %s [L] [p]\n", argv[0]);
		return 1;
	}

	const uint32_t n_iter = 1'000'000;
	
	int L;
	double p;

	sscanf(argv[1], "%d", &L);
	sscanf(argv[2], "%lf", &p);

	if(L < 0)
	{
		printf("Error: L must be positive\n");
		return 1;
	}

	if(!((p > 0.0) && (p < 1.0)))
	{
		printf("Error: p must be in between 0.0 and 1.0\n");
		return 1;
	}

	fprintf(stderr, "Processing at rank = %d, size = %d\n", mpi_rank, mpi_size);

	uint32_t n_success = 0u;
	chrono::microseconds node_dur;
	UnionFindDecoder<Lattice2D> decoder(L);
	for(uint32_t k = mpi_rank; k < n_iter; k += mpi_size)
	{
		auto [error_x, error_z] = 
			generate_errors(2*L*L, L, p, re, noise_type);

		auto synd_x = calc_syndromes(lattice, error_x, ErrorType::X);
		auto synd_z = calc_syndromes(lattice, error_z, ErrorType::Z);

		Eigen::ArrayXi error_total_x = error_x.col(L-1);
		Eigen::ArrayXi error_total_z = error_z.col(L-1);

		const auto [measurement_error_x, measurement_error_z]
			= create_measurement_errors(re, L*L, L, p, noise_type);
		
		add_measurement_noise(L, re, synd_x, measurement_error_x);
		add_measurement_noise(L, re, synd_z, measurement_error_z);

		layer_syndrome_diff(L, synd_x);
		layer_syndrome_diff(L, synd_z);

		auto start = chrono::high_resolution_clock::now();
		decoder.clear();
		auto corrections_x = decoder.decode(synd_x);
		decoder.clear();
		auto corrections_z = decoder.decode(synd_z);
		auto end = chrono::high_resolution_clock::now();

		if((!has_logical_error(L, error_total_x, corrections_x, ErrorType::X)) &&
				(!has_logical_error(L, error_total_z, corrections_z, ErrorType::Z)))
		{
			++n_success ;
		}

		auto dur = chrono::duration_cast<chrono::microseconds> (end-start);
		node_dur += dur;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	int64_t dur_in_microseconds = node_dur.count();
	int64_t total_dur_in_microseconds = 0;
	MPI_Allreduce(&dur_in_microseconds, &total_dur_in_microseconds,
			1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	uint32_t total_success = 0;
	MPI_Allreduce(&n_success, &total_success, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);


	if(mpi_rank == 0)
	{
		char filename[255];
		sprintf(filename, "out_L%d_%06d.json", L, int(p*100000+0.5));
		std::ofstream out_data(filename);
		nlohmann::json out_j;
		out_j["L"] = L;
		out_j["average_microseconds"] = double(total_dur_in_microseconds)/n_iter;
		out_j["p"] = p;
		out_j["accuracy"] = double(total_success)/n_iter;

		out_data << out_j.dump(0);
	}

	MPI_Finalize();
	return 0;
}
