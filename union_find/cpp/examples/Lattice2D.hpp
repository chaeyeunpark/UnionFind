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

#include "toric_utils.hpp"
#include "utility.hpp"

namespace UnionFindCPP
{
class Lattice2D
{
private:
	int L_;

public:
	using Vertex = int;

	explicit Lattice2D(int L) : L_{L} { }

	inline int to_vertex_index(int row, int col) const
	{
		return UnionFindCPP::to_vertex_index(L_, row,
											 col); // call function from the parent scope
	}

	constexpr int vertex_connection_count(Vertex v) const { return 4; }

	std::array<Vertex, 4> vertex_connections(Vertex v) const
	{
		int row = v / L_;
		int col = v % L_;

		return {
			to_vertex_index(row - 1, col),
			to_vertex_index(row + 1, col),
			to_vertex_index(row, col - 1),
			to_vertex_index(row, col + 1),
		};
	}

	inline int num_vertices() const { return L_ * L_; }

	inline int num_edges() const { return 2 * L_ * L_; }

	inline int edge_idx(const Edge& edge) const { return to_edge_idx(L_, edge); }

	inline Edge to_edge(int edge_index) const
	{
		return UnionFindCPP::to_edge(L_, edge_index);
	}
};
} // namespace UnionFindCPP
