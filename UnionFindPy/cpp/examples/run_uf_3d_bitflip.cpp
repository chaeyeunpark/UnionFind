// Copyright (C) 2021 UnionFind++ authors
//
// This file is part of UnionFind++.
//
// UnionFind++ is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// UnionFind++ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with UnionFind++.  If not, see <https://www.gnu.org/licenses/>.
#include "Decoder.hpp"
#include "LatticeCubic.hpp"
#include "LazyDecoder.hpp"
#include "error_utils.hpp"
#include "runner_utils.hpp"
#include "toric_utils.hpp"

#include <Eigen/Dense>
#ifdef USE_MPI
#pragma message("Build with MPI")
#include <mpi.h>
#endif
#include <nlohmann/json.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

auto main(int argc, char* argv[]) -> int
{
	namespace chrono = std::chrono;
	using UnionFindCPP::ArrayXu, UnionFindCPP::Decoder, UnionFindCPP::ErrorType,
		UnionFindCPP::LatticeCubic, UnionFindCPP::LazyDecoder, UnionFindCPP::NoiseType,
		UnionFindCPP::add_measurement_noise, UnionFindCPP::add_corrections,
		UnionFindCPP::layer_syndrome_diff;

	const auto noise_type = NoiseType::X;
	const uint32_t n_iter = 1'000'000;

	int mpi_rank = 0;
	int mpi_size = 1;

#ifdef USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	uint32_t L = 0;
	double p = 0.0;
	try
	{
		std::tie(L, p) = parse_args(argc, argv);
	}
	catch(std::exception& e)
	{
		std::cout << e.what() << std::endl;
		return 1;
	}

	std::random_device rd;
	std::default_random_engine re{rd()};

#ifdef USE_MPI
	fmt::print(stderr, "Processing at rank = {}, size = {}\n", mpi_rank, mpi_size);
#endif

	unsigned int n_success = 0U;
	chrono::microseconds node_dur{};
#ifdef USE_LAZY
	LazyDecoder<LatticeCubic> lazy_decoder(L);
#endif
	LatticeCubic lattice(L);
	Decoder<LatticeCubic> decoder(L);
	for(uint32_t k = mpi_rank; k < n_iter; k += mpi_size)
	{
		auto [error_x, error_z] = generate_errors(2 * L * L, L, p, re, noise_type);

		auto synd_x = calc_syndromes(lattice, error_x, ErrorType::X);

		ArrayXu error_total_x = error_x.col(L - 1);

		const auto [measurement_error_x, measurement_error_z]
			= create_measurement_errors(re, L * L, L, p, noise_type);

		add_measurement_noise(L, synd_x, measurement_error_x);

		layer_syndrome_diff(L, synd_x);

		auto start = chrono::high_resolution_clock::now();
#ifdef USE_LAZY
		// Process lazy decoder
		auto [success_x, decoding_x] = lazy_decoder.decode(synd_x);
		if(!success_x)
		{
			decoder.clear();
			auto decoding_uf = decoder.decode(synd_x);
			decoding_x.insert(decoding_x.end(), decoding_uf.begin(), decoding_uf.end());
		}

#else
		decoder.clear();
		auto decoding_x = decoder.decode(synd_x);
#endif
		auto end = chrono::high_resolution_clock::now();

		if(!has_logical_error(L, error_total_x, decoding_x, ErrorType::X))
		{
			++n_success;
		}

		auto dur = chrono::duration_cast<chrono::microseconds>(end - start);
		node_dur += dur;
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	unsigned long dur_in_microseconds = node_dur.count();
	unsigned long total_dur_in_microseconds = 0;
	MPI_Allreduce(&dur_in_microseconds, &total_dur_in_microseconds, 1, MPI_UNSIGNED_LONG,
				  MPI_SUM, MPI_COMM_WORLD);

	fmt::print("rank: {}, dur_in_microseconds: {}, total_dur_in_microsecond: {}\n",
			   mpi_rank, dur_in_microseconds, total_dur_in_microseconds);

	unsigned int total_success = 0;
	MPI_Allreduce(&n_success, &total_success, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
	unsigned long total_dur_in_microseconds = node_dur.count();
	unsigned int total_success = n_success;
#endif

	if(mpi_rank == 0)
	{
		save_to_json(L, p, static_cast<double>(total_dur_in_microseconds) / n_iter,
					 static_cast<double>(total_success) / n_iter);
	}

#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
