#include <random>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <atomic>
#include <chrono>

#include <Eigen/Dense>

#include <sstream>

#include <mpi.h>

#include "ErrorGenerator.hpp"
#include "Lattice2D.hpp"
#include "LatticeCubic.hpp"
#include "UnionFind.hpp"

void print_qubit_errors(const int L, const Eigen::ArrayXXi& qubit_errors)
{
	for(int h = 0; h < L; ++h)
	{
		Eigen::ArrayXi layer_error = qubit_errors.col(h);
		std::cerr << "Layer: " << h << "\n";
		auto m = Eigen::Map<Eigen::MatrixXi>(layer_error.data(), 2*L, L);
		std::cerr << m.transpose() << "\n";
	}
}

template<class RandomEngine>
std::pair<std::vector<int>, std::vector<int>>
noisy_syndromes(const LatticeCubic& lattice, double p, RandomEngine& re)
{
	const int L = lattice.getL();
	Eigen::ArrayXXi qubit_errors = Eigen::ArrayXXi::Zero(2*L*L, L); //repetetions L
	auto syndromes = std::vector<int>(lattice.num_vertices(), 0); //length L^3

	std::uniform_real_distribution<> dist(0., 1.);

	for(int h = 0; h < L; ++h)
	{
		for(int eidx = 0; eidx < 2*L*L; ++eidx)
		{
			qubit_errors(eidx, h) = int(dist(re) < p);
		}
	}

	//cumsum
	for(int h = 1; h < L; ++h)
	{
		qubit_errors.col(h) += qubit_errors.col(h-1);
	}
	qubit_errors = qubit_errors.unaryExpr([](int x){ return x % 2;});

	std::vector<int> error_total(2*L*L);
	for(size_t n = 0; n < 2*L*L; ++n)
	{
		error_total[n] = qubit_errors(n, L-1);
	}

	for(int h = 0; h < L; ++h)
	{
		auto layer_error = std::vector<int>(qubit_errors.data() + h*2*L*L, 
				qubit_errors.data() + (h+1)*2*L*L);
		auto layer_syndrome = z_error_to_syndrome_x(L, layer_error);
		std::copy(layer_syndrome.begin(), layer_syndrome.end(),
				syndromes.begin() + h*L*L);
	}

	//add measurement noise
	Eigen::ArrayXXi measurement_errors = Eigen::ArrayXXi::Zero(L*L, L); 	
	for(int h = 0; h < L-1; ++h)//perfect measurements in the last round
	{
		for(int eidx = 0; eidx < L*L; ++eidx)
		{
			measurement_errors(eidx, h) = int(dist(re) < p);
		}
	}
	Eigen::Map<Eigen::ArrayXXi> syndromes_map(syndromes.data(), L*L, L);
	syndromes_map += measurement_errors;

	for(int h = L-1; h >= 1; --h)
	{
		syndromes_map.col(h) -= syndromes_map.col(h-1);
	}
	syndromes_map = syndromes_map.unaryExpr([](int x){ return (x+2)%2;});

	return std::make_pair(error_total, syndromes);
}

bool has_logical_error(int L, std::vector<int>& error_total, const std::vector<Edge>& corrections)
{
	for(auto edge: corrections)
	{
		if(edge.u/(L*L) == edge.v/(L*L)) //edge is timelike
		{
			auto corr_edge = Edge{edge.u % (L*L), edge.v % (L*L)};
			int corr_qubit = decoder_edge_to_qubit_idx(L, corr_edge, ErrorType::Z);
			error_total[corr_qubit] += 1;
		}
	}

	for(int& e: error_total)
		e %= 2;

	return logical_error(L, error_total, ErrorType::Z);
}


void run_client()
{
	LatticeCubic lattice{L};
	uint32_t n_success = 0u;
	chrono::microseconds node_dur;
	UnionFindDecoder decoder(lattice);

	for(uint32_t k = mpi_rank; k < n_iter; k += mpi_size)
	{
		auto [error_total, syndromes] = noisy_syndromes(lattice, p, re);
		decoder.clear();
		auto start = chrono::high_resolution_clock::now();
		auto corrections = decoder.decode(syndromes);
		auto end = chrono::high_resolution_clock::now();

		auto dur = chrono::duration_cast<chrono::microseconds>(end - start);
		node_dur += dur;

		if(!has_logical_error(L, error_total, corrections))
		{
			++n_success ;
		}
	}
}


int main(int argc, char* argv[])
{
	namespace chrono = std::chrono;

	int mpi_rank, mpi_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	std::random_device rd;
	std::default_random_engine re{rd()};


	if(argc != 3)
	{
		printf("Usage: %s [L] [p]\n", argv[0]);
		return 1;
	}

	const uint32_t n_iter = 10'000'000;
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

	fprintf(stderr, "Processing at rank = %d, size = %d", mpi_rank, mpi_size);

	if(rank == 0)
	{
	}
	else
	{
		run_client();
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
