#include <random>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <chrono>

#include <Eigen/Dense>

#ifdef USE_MPI
#pragma message ("Build with MPI")
#include <mpi.h>
#endif

#include "Lattice2D.hpp"
#include "UnionFind.hpp"
#include "error_utils.hpp"
#include "toric_utils.hpp"

#ifdef USE_LAZY
#pragma message ("Build with Lazy decoder")
#include "LazyDecoder.hpp"
#endif

int main(int argc, char* argv[])
{
	namespace chrono = std::chrono;

	const auto noise_type = NoiseType::Depolarizing;
	const uint32_t n_iter = 1'000'000;

	int mpi_rank = 0;
	int mpi_size = 1;

#ifdef USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	std::random_device rd;
	std::default_random_engine re{rd()};


	if(argc != 3)
	{
		printf("Usage: %s [L] [p]\n", argv[0]);
		return 1;
	}

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

#ifdef USE_MPI
	fprintf(stderr, "Processing at rank = %d, size = %d\n", mpi_rank, mpi_size);
#endif

	uint32_t n_success = 0u;
	chrono::microseconds node_dur{};
#ifdef USE_LAZY
	LazyDecoder<Lattice2D> lazy_decoder(L);
#endif
	
	UnionFindDecoder<Lattice2D> decoder(L);
	for(uint32_t k = mpi_rank; k < n_iter; k += mpi_size)
	{
		auto [x_errors, z_errors] = create_errors(re, decoder.num_edges(),
				p, noise_type);

		auto synd_x = errors_to_syndromes(L, x_errors, ErrorType::X);
		auto synd_z = errors_to_syndromes(L, z_errors, ErrorType::Z);

		auto start = chrono::high_resolution_clock::now();
#ifdef USE_LAZY
		// Process lazy decoder
		auto [success_x, decoding_x] = lazy_decoder.decode(synd_x);
		if (!success_x)
		{
			decoder.clear();
			auto decoding_uf =  decoder.decode(synd_x);
			decoding_x.insert(decoding_x.end(), decoding_uf.begin(), decoding_uf.end());
		}

		auto [success_z, decoding_z] = lazy_decoder.decode(synd_z);
		if (!success_z)
		{
			decoder.clear();
			auto decoding_uf =  decoder.decode(synd_z);
			decoding_z.insert(decoding_z.end(), decoding_uf.begin(), decoding_uf.end());
		}

#else
		decoder.clear();
		auto decoding_x = decoder.decode(synd_x);
		decoder.clear();
		auto decoding_z = decoder.decode(synd_z);
#endif

		add_corrections(L, decoding_x, x_errors, ErrorType::X);
		add_corrections(L, decoding_z, z_errors, ErrorType::Z);
		auto end = chrono::high_resolution_clock::now();

		if((!logical_error(L, x_errors, ErrorType::X)) && 
				!(logical_error(L, z_errors, ErrorType::Z)))
		{
			++n_success;
		}

		auto dur = chrono::duration_cast<chrono::microseconds> (end-start);
		node_dur += dur;
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	
	int64_t dur_in_microseconds = node_dur.count();
	int64_t total_dur_in_microseconds = 0;
	MPI_Allreduce(&dur_in_microseconds, &total_dur_in_microseconds,
			1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	printf("rank: %d, dur_in_microseconds: %ld, total_dur_in_microsecond: %ld\n", 
			mpi_rank, dur_in_microseconds, total_dur_in_microseconds);

	uint32_t total_success = 0;
	MPI_Allreduce(&n_success, &total_success, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
	uint64_t total_dur_in_microseconds = node_dur.count();
	uint32_t total_success = n_success;
#endif

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

#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
