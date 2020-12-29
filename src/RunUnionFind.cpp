#include <random>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>

#include <chrono>

#include "ErrorGenerator.hpp"
#include "UnionFind.hpp"

int main(int argc, char* argv[])
{
	namespace chrono = std::chrono;
	std::random_device rd;
	ErrorGenerator<std::default_random_engine> gen;


	if(argc != 3)
	{
		printf("Usage: %s [L]\n", argv[0]);
		return 1;
	}

	uint32_t L;
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

	chrono::microseconds total_dur{0};

	const uint32_t n_iter = 1'000;

	fprintf(stderr, "#L = %d, p = %f\n", L, p);

	int acc = 0;
	
	UnionFindDecoder decoder(L);
	for(int n = 0; n < n_iter; ++n)
	{
		auto [x_errors, z_errors] = gen.get_errors(L, p, NoiseType::Depolarizing);

		auto synd_x = errors_to_syndromes(L, x_errors, ErrorType::X);
		auto synd_z = errors_to_syndromes(L, z_errors, ErrorType::Z);


		auto start = chrono::high_resolution_clock::now();
		decoder.clear();
		auto decoding_x = decoder.decode(synd_x);
		decoder.clear();
		auto decoding_z = decoder.decode(synd_z);
		auto end = chrono::high_resolution_clock::now();

		add_corrections(L, decoding_x, x_errors, ErrorType::X);
		add_corrections(L, decoding_z, z_errors, ErrorType::Z);


		if((!logical_error(L, x_errors, ErrorType::X)) && !(logical_error(L, z_errors, ErrorType::Z)))
		{
			++acc;
		}

		printf("Processing L=%d p=%f n=%d\n", L, p, n);
		fflush(stdout);

		auto dur = chrono::duration_cast<chrono::microseconds> (end-start);
		total_dur += dur;
	}

	char filename[255];
	sprintf(filename, "out_L%2d_%04d.json", L, int(p*1000+0.5));
	std::ofstream out_data(filename);
	nlohmann::json out_j;
	out_j["L"] = L;
	out_j["total_dur"] = total_dur.count();
	out_j["p"] = p;
	out_j["accuracy"] = double(acc)/n_iter;

	out_data << out_j.dump(0);

	return 0;
}
