#include <random>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>

#include <chrono>

#include <Eigen/Dense>

#include "LatticeCubic.hpp"
#include "UnionFind.hpp"

#include "cpp_utils.hpp"

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
std::pair<Eigen::ArrayXXi, Eigen::ArrayXXi>
generate_errors(const LatticeCubic& lattice, double p, RandomEngine& re,
		NoiseType noise_type)
{
	const int L = lattice.getL();
	Eigen::ArrayXXi qubit_errors_x = Eigen::ArrayXXi::Zero(2*L*L, L); //repetetions L
	Eigen::ArrayXXi qubit_errors_z = Eigen::ArrayXXi::Zero(2*L*L, L); //repetetions L

	for(int r = 0; r < L; ++r)
	{
		auto [layer_error_x, layer_error_z] = create_errors(re, 2*L*L, p, noise_type);
		qubit_errors_x.col(r) = layer_error_x;
		qubit_errors_z.col(r) = layer_error_z;
	}

	//cumsum
	for(int h = 1; h < L; ++h)
	{
		qubit_errors_x.col(h) += qubit_errors_x.col(h-1);
		qubit_errors_z.col(h) += qubit_errors_z.col(h-1);
	}
	qubit_errors_x = qubit_errors_x.unaryExpr([](int x){ return x % 2;});
	qubit_errors_z = qubit_errors_z.unaryExpr([](int x){ return x % 2;});

	return std::make_pair(qubit_errors_x, qubit_errors_z);
}

std::vector<int>
calc_syndromes(const LatticeCubic& lattice, 
		const Eigen::ArrayXXi& errors, ErrorType error_type)
{
	const int L = lattice.getL();
	std::vector<int> syndromes(lattice.num_vertices());
	for(int h = 0; h < L; ++h)
	{
		Eigen::ArrayXi layer_error = errors.col(h);
		auto layer_syndrome = [L, error_type, &layer_error]
		{
			switch(error_type)
			{
			case ErrorType::Z:
				return z_error_to_syndrome_x(L, layer_error);
			case ErrorType::X:
				return x_error_to_syndrome_z(L, layer_error);
			}
		}();
		std::copy(layer_syndrome.begin(), layer_syndrome.end(),
				syndromes.begin() + h*L*L);
	}

	return syndromes;
}

template<typename RandomEngine>
void add_measurement_noise(const int L, RandomEngine& re, 
		std::vector<int>& syndromes, const Eigen::ArrayXXi& measurement_error)
{
	Eigen::Map<Eigen::ArrayXXi> syndromes_map(syndromes.data(), L*L, L);
	syndromes_map += measurement_error;
}


void layer_syndrome_diff(const int L, std::vector<int>& syndromes)
{
	Eigen::Map<Eigen::ArrayXXi> syndromes_map(syndromes.data(), L*L, L);
	for(int h = L-1; h >= 1; --h)
	{
		syndromes_map.col(h) -= syndromes_map.col(h-1);
	}
	syndromes_map = syndromes_map.unaryExpr([](int x){ return (x+2)%2;});
}

bool has_logical_error(int L, Eigen::ArrayXi& error_total, 
		const std::vector<Edge>& corrections, ErrorType error_type)
{
	for(auto edge: corrections)
	{
		if(edge.u/(L*L) == edge.v/(L*L)) //edge is timelike
		{
			auto corr_edge = Edge{edge.u % (L*L), edge.v % (L*L)};
			int corr_qubit = decoder_edge_to_qubit_idx(L, corr_edge, error_type);
			error_total[corr_qubit] += 1;

#ifndef NDEBUG
			std::cerr << "[" << corr_edge.u << ", " << corr_edge.v << "], qubit: "
				<< corr_qubit << std::endl;
#endif
		}
	}

	return logical_error(L, error_total, error_type);
}

int main(int argc, char* argv[])
{
	namespace chrono = std::chrono;
	std::random_device rd;
	std::default_random_engine re{rd()};
	auto noise_type = NoiseType::Depolarizing;

	if(argc != 3)
	{
		printf("Usage: %s [L] [p]\n", argv[0]);
		return 1;
	}

	const uint32_t n_iter = 100'000;
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

	LatticeCubic lattice{L};
	UnionFindDecoder<LatticeCubic> decoder(L);
	uint32_t n_success = 0;
	chrono::microseconds total_dur{0};
	for(int n = 0; n < n_iter; ++n)
	{
		std::cerr << "Processing " << n << std::endl;
		auto [error_x, error_z] = 
			generate_errors(lattice, p, re, noise_type);

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
		total_dur += dur;
	}
	char filename[255];
	sprintf(filename, "out_L%d_%05d.json", L, int(p*10000+0.5));
	std::ofstream out_data(filename);
	nlohmann::json out_j;
	out_j["L"] = L;
	out_j["total_dur"] = total_dur.count();
	out_j["p"] = p;
	out_j["accuracy"] = double(n_success)/n_iter;

	out_data << out_j.dump(0);
	return 0;
}
