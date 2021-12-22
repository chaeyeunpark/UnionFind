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
#pragma once
#include "utility.hpp"
#include <vector>

namespace UnionFindCPP
{
template<typename Lattice> class LazyDecoder
{
private:
	const Lattice lattice_;
	std::vector<Edge> all_edges_;

public:
	template<typename... Args> LazyDecoder(Args&&... args) : lattice_{args...}
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

		for(Edge edge : all_edges_)
		{
			if(syndromes[edge.u] == 1 && syndromes[edge.v] == 1)
			{
				corrections.emplace_back(std::move(edge));
			}
		}

		for(const auto& edge : corrections)
		{
			syndromes[edge.u] ^= 1;
			syndromes[edge.v] ^= 1;
		}

		bool success = true;

		for(const auto syndrome : syndromes)
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
} // namespace UnionFindCPP
