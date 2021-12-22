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

#include "Lattice2D.hpp"
#include "toric_utils.hpp"
#include "utility.hpp"

namespace UnionFindCPP
{
class LatticeCubic
{
private:
	const int L_;

	inline int to_vertex_index(int row, int col, int h) const
	{
		return UnionFindCPP::to_vertex_index(L_, row, col) + h * (L_ * L_);
	}

public:
	using Vertex = int;

	explicit LatticeCubic(int L) : L_{L} { }

	inline int getL() const { return L_; }

	constexpr int vertex_connection_count(Vertex v) const
	{
		int L = L_;

		int h = v / (L * L);

		if(h == L - 1 || h == 0) return 5;

		return 6;
	}

	std::vector<Vertex> vertex_connections(Vertex v) const
	{
		int L = L_;

		int h = v / (L * L);
		int row = (v / L) % L;
		int col = v % L;

		auto res = std::vector<Vertex>{
			to_vertex_index(row - 1, col, h),
			to_vertex_index(row + 1, col, h),
			to_vertex_index(row, col - 1, h),
			to_vertex_index(row, col + 1, h),
		};
		if(h < L - 1) res.emplace_back(to_vertex_index(row, col, h + 1));
		if(h > 0) res.emplace_back(to_vertex_index(row, col, h - 1));

		return res;
	}

	inline int num_vertices() const { return L_ * L_ * L_; }

	inline int num_edges() const { return 3 * L_ * L_ * L_ - L_ * L_; }

	inline int edge_idx(const Edge& edge) const
	{
		int L = L_;
		int uh = edge.u / (L * L);

		if((edge.u / (L * L)) == (edge.v / (L * L))) // edge is spacelike
		{
			return to_edge_idx(L, Edge{edge.u % (L * L), edge.v % (L * L)})
				   + 3 * L * L * uh;
		}

		{ // edge is in the time direction
			int row = (edge.u / L) % L;
			int col = edge.u % L;
			return 3 * L * L * uh + 2 * L * L + L * row + col;
		}
	}

	Edge to_edge(int edge_index) const
	{
		int L = L_;
		int h = edge_index / (3 * L * L);
		int layer_idx = edge_index % (3 * L * L);
		if(layer_idx >= 2 * L * L) // edge is in the time direction
		{
			auto [row, col] = UnionFindCPP::vertex_to_coord(L, layer_idx - 2 * L * L);
			return Edge{to_vertex_index(row, col, h), to_vertex_index(row, col, h + 1)};
		}
		else // edge is spacelike
		{
			Edge e_2d = UnionFindCPP::to_edge(L, layer_idx);
			return Edge{e_2d.u + h * L * L, e_2d.v + h * L * L};
		}
	}
};
} // namespace UnionFindCPP
