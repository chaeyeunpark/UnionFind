#include <random>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>

#include <chrono>

#include "Lattice2D.hpp"
#include "UnionFind.hpp"

#include "utility.hpp"
#include "cpp_utils.hpp"


void add_corrections(const int L, const std::vector<Edge>& corrections, 
		Eigen::ArrayXi& error, ErrorType error_type)
{
	for(auto e: corrections)
	{
		auto idx = decoder_edge_to_qubit_idx(L, e, error_type);
		error[idx] += 1;
	}
}

template<typename Lattice>
class LazyDecoder
{
private:
	const Lattice lattice_;
	std::vector<Edge> all_edges_;

public:
	template<typename ...Args>
	LazyDecoder(Args&&... args)
		: lattice_{args...}
	{
		const auto num_edges = lattice_.num_edges();

		for(int i = 0; i < num_edges; ++i)
		{
			all_edges_.emplace_back(lattice_.to_edge(i));
		}
	}

	std::pair<bool, std::vector<Edge>> decode(std::vector<int>& syndromes)
	{
		std::vector<Edge> corrections;

		for(Edge edge: all_edges_)
		{
			if (syndromes[edge.u] == 1 && syndromes[edge.v] == 1)
			{
				corrections.emplace_back(std::move(edge));
			}
		}

		for(const auto& edge: corrections)
		{
			syndromes[edge.u] ^= 1;
			syndromes[edge.v] ^= 1;
		}

		bool success = true;

		for(const auto syndrome: syndromes)
		{
			if(syndrome == 1)
			{
				success = false;
				break;
			}
		}

		return std::make_pair(success, corrections);
	}
};

int main(int argc, char* argv[])
{
	namespace chrono = std::chrono;
	std::random_device rd;
	std::default_random_engine re{rd()};

	if(argc != 3)
	{
		printf("Usage: %s [L]\n", argv[0]);
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

	chrono::microseconds total_dur{0};

	const uint32_t n_iter = 1'000'000;

	fprintf(stderr, "#L = %d, p = %f\n", L, p);

	int acc = 0;
	
	LazyDecoder<Lattice2D> lazy_decoder(L);
	UnionFindDecoder<Lattice2D> decoder(L);
	for(int n = 0; n < n_iter; ++n)
	{
		auto [x_errors, z_errors] = create_errors(re, decoder.num_edges(),
				p, NoiseType::X);

		auto synd_x = errors_to_syndromes(L, x_errors, ErrorType::X);

		auto start = chrono::high_resolution_clock::now();
		decoder.clear();

		auto [success, decoding] = lazy_decoder.decode(synd_x);
		if (!success)
		{
			auto decoding_uf =  decoder.decode(synd_x);
			decoding.insert(decoding.end(), decoding_uf.begin(), decoding_uf.end());
		}
		auto end = chrono::high_resolution_clock::now();

		add_corrections(L, decoding, x_errors, ErrorType::X);
		//add_corrections(L, decoding_z, z_errors, ErrorType::Z);


		if((!logical_error(L, x_errors, ErrorType::X)))
		{
			++acc;
		}

		printf("Processing L=%d p=%f n=%d\n", L, p, n);
		fflush(stdout);

		auto dur = chrono::duration_cast<chrono::microseconds> (end-start);
		total_dur += dur;
	}

	char filename[255];
	sprintf(filename, "out_lazy_L%02d_%04d.json", L, int(p*1000+0.5));
	std::ofstream out_data(filename);
	nlohmann::json out_j;
	out_j["L"] = L;
	out_j["total_dur"] = total_dur.count();
	out_j["average_dur"] = double(total_dur.count()) / n_iter;
	out_j["p"] = p;
	out_j["accuracy"] = double(acc)/n_iter;

	out_data << out_j.dump(0);

	return 0;
}
