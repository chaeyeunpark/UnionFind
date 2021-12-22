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
#include "toric_utils.hpp"
#include "LatticeCubic.hpp"

namespace UnionFindCPP
{
auto z_error_to_syndrome_x(const int L, const Eigen::ArrayXi& z_error) -> std::vector<int>
{
	std::vector<int> syndromes_array(L * L, 0);
	for(int n = 0; n < 2 * L * L; ++n)
	{
		if(z_error[n] == 0) { continue; }
		Edge e = to_edge(L, n);
		syndromes_array[e.u] += 1;
		syndromes_array[e.v] += 1;
	}
	for(auto& u : syndromes_array) { u %= 2; }

	return syndromes_array;
}

auto x_error_to_syndrome_z(const int L, const Eigen::ArrayXi& x_error) -> std::vector<int>
{
	std::vector<int> syndromes_array(L * L, 0);
	for(int n = 0; n < 2 * L * L; ++n)
	{
		if(x_error[n] == 0) { continue; }
		Edge e = to_edge(L, n);

		if(is_horizontal(L, e))
		{
			auto u = left(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			syndromes_array[u] += 1;
			syndromes_array[to_vertex_index(L, row - 1, col)] += 1;
		}
		else
		{
			auto u = lower(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			syndromes_array[u] += 1;
			syndromes_array[to_vertex_index(L, row, col - 1)] += 1;
		}
	}
	for(auto& u : syndromes_array) { u %= 2; }
	return syndromes_array;
}

auto decoder_edge_to_qubit_idx(const int L, Edge e, ErrorType error_type) -> int
{
	switch(error_type)
	{
	case ErrorType::X:
		// each edge is a qubit in the dual lattice
		if(is_horizontal(L, e))
		{
			auto u = left(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			return L * ((row + 1) % L) + ((col + 1) % L);
		}
		else
		{
			auto u = upper(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			return L * (row % L) + col + L * L;
		}
	case ErrorType::Z:
		return to_edge_idx(L, e);
	}
	__builtin_unreachable();
}

auto to_edge(const int L, int edge_index) -> Edge
{
	int row = edge_index / L; // % L is done in to_vertex_index
	int col = edge_index % L;
	if((edge_index / (L * L)) == 0) // vertical edge
	{
		return Edge(to_vertex_index(L, row, col), to_vertex_index(L, row - 1, col));
	}
	// horizontal edge
	return Edge(to_vertex_index(L, row, col), to_vertex_index(L, row, col + 1));
}

auto calc_syndromes(const LatticeCubic& lattice, const Eigen::ArrayXXi& errors,
					ErrorType error_type) -> std::vector<int>
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
			__builtin_unreachable();
		}();
		std::copy(layer_syndrome.begin(), layer_syndrome.end(),
				  syndromes.begin() + h * L * L);
	}

	return syndromes;
}

void add_corrections(const int L, const std::vector<Edge>& corrections,
					 Eigen::ArrayXi& error, ErrorType error_type)
{
	for(auto e : corrections)
	{
		auto idx = decoder_edge_to_qubit_idx(L, e, error_type);
		error[idx] += 1;
	}
}

auto logical_error(const int L, const Eigen::ArrayXi& error, ErrorType error_type) -> bool
{
	// one may use a counting algorithm for testing logical error
	int sum1 = 0;
	int sum2 = 0;
	switch(error_type)
	{
	case ErrorType::X:
		// need to think in a dual lattice
		for(int u = 0; u < L * L; u += L) { sum1 += error[u]; }
		for(int u = L * L; u < L * L + L; ++u) { sum2 += error[u]; }
		break;
	case ErrorType::Z:
		for(int u = 0; u < L; ++u) { sum1 += error[u]; }
		for(int u = L * L; u < 2 * L * L; u += L) { sum2 += error[u]; }
		break;
	}

	return (sum1 % 2 == 1) || (sum2 % 2 == 1);
}

auto has_logical_error(int L, Eigen::ArrayXi& error_total,
					   const std::vector<Edge>& corrections, ErrorType error_type) -> bool
{
	for(auto edge : corrections)
	{
		if(edge.u / (L * L) == edge.v / (L * L)) // edge is spacelike
		{
			auto corr_edge = Edge{edge.u % (L * L), edge.v % (L * L)};
			int corr_qubit = decoder_edge_to_qubit_idx(L, corr_edge, error_type);
			error_total[corr_qubit] += 1;
		}
	}
	return logical_error(L, error_total, error_type);
}
} // namespace UnionFindCPP
