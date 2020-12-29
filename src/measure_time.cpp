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
	using nlohmann::json;
	std::random_device rd;
	ErrorGenerator<std::default_random_engine> gen;


	if(argc != 2)
	{
		printf("Usage: %s [p]\n", argv[0]);
		return 1;
	}

	double p;

	sscanf(argv[1], "%lf", &p);

	if(!((p > 0.0) && (p < 1.0)))
	{
		printf("Error: p must be in between 0.0 and 1.0\n");
		return 1;
	}

	const uint32_t n_iter = 100'000;
	json out_j =  json::array();

	int old_L = 0;

	for(unsigned int n_qubit = 200; n_qubit <= 40000; n_qubit += 200)
	{
		int L = int(sqrt(n_qubit/2.0)+0.5);
		if(L == old_L)
			continue;
		chrono::microseconds total_dur{0};
		fprintf(stderr, "#L = %d, p = %f\n", L, p);

		int acc = 0;
		UnionFindDecoder decoder(L);
		for(int n = 0; n < n_iter; ++n)
		{
			auto [x_errors, z_errors] = gen.get_errors(L, p, NoiseType::Z);

			auto synd_z = errors_to_syndromes(L, z_errors, ErrorType::Z);
			

			auto start = chrono::high_resolution_clock::now();
			decoder.clear();
			auto decoding_z = decoder.decode(synd_z);
			auto end = chrono::high_resolution_clock::now();

			add_corrections(L, decoding_z, z_errors, ErrorType::Z);

			if(!logical_error(L, z_errors, ErrorType::Z))
			{
				++acc;
			}

			printf("Processing L=%d p=%f n=%d\n", L, p, n);
			fflush(stdout);

			auto dur = chrono::duration_cast<chrono::microseconds> (end-start);
			total_dur += dur;
		}

		json j = { {"L", L}, {"duration", total_dur.count()}};
		out_j.push_back(j);

		old_L = L;

	}
	char filename[255];
	sprintf(filename, "out_%04d.json", int(p*1000+0.5));
	std::ofstream out_data(filename);


	out_data << out_j.dump(0);


	return 0;
}
