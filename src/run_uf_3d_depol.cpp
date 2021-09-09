#include <random>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>

#include <chrono>

#include <Eigen/Dense>

#include "LatticeCubic.hpp"
#include "UnionFind.hpp"
#include "error_utils.hpp"
#include "toric_utils.hpp"

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

	const uint32_t n_iter = 1'000;
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
		total_dur += dur;
	}
	char filename[255];
	sprintf(filename, "out_L%d_%05d.json", L, int(p*10000+0.5));
	std::ofstream out_data(filename);
	nlohmann::json out_j;
	out_j["L"] = L;
	out_j["total_dur"] = total_dur.count();
	out_j["average_microseconds"] = double(total_dur.count())/n_iter;
	out_j["p"] = p;
	out_j["accuracy"] = double(n_success)/n_iter;

	out_data << out_j.dump(0);
	return 0;
}
