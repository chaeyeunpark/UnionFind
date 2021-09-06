#pragma once
#include <vector>
#include "utility.hpp"

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
