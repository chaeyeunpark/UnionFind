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
	uint32_t L_;

public:
	using Vertex = uint32_t;

	explicit Lattice2D(uint32_t L) : L_{L} { }

	[[nodiscard]] inline auto to_vertex_index(uint32_t row, uint32_t col) const
		-> uint32_t
	{
		return UnionFindCPP::to_vertex_index(L_, row,
											 col); // call function from the parent scope
	}

	[[nodiscard]] constexpr static auto vertex_connection_count(Vertex /*v*/) -> uint32_t
	{
		return 4;
	}

	[[nodiscard]] auto vertex_connections(Vertex v) const -> std::array<Vertex, 4>
	{
		uint32_t row = v / L_;
		uint32_t col = v % L_;

		return {
			to_vertex_index(row - 1, col),
			to_vertex_index(row + 1, col),
			to_vertex_index(row, col - 1),
			to_vertex_index(row, col + 1),
		};
	}

	[[nodiscard]] inline auto num_vertices() const -> uint32_t { return L_ * L_; }

	[[nodiscard]] inline auto num_edges() const -> uint32_t { return 2 * L_ * L_; }

	[[nodiscard]] inline auto edge_idx(const Edge& edge) const -> uint32_t
	{
		return to_edge_idx(L_, edge);
	}

	[[nodiscard]] inline auto to_edge(uint32_t edge_index) const -> Edge
	{
		return UnionFindCPP::to_edge(L_, edge_index);
	}
};
} // namespace UnionFindCPP
