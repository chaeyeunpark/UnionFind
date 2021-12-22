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
#include "LatticeFromParity.hpp"

#include <random>
#include <set>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using UnionFindCPP::LatticeFromParity;
/**
 * The Lattice2D from v0.1 uses the simplest ordering.
 *
 * For example, for Lx=4, Ly=2 lattice, qubits are labeled as:
 *
 *   	─┼─8─┼─9─┼───┼─
 *   	 4   5   6   7
 *   	─┼─0─┼─1─┼─2─┼3
 *
 * and the X stabiliers are labeles as:
 *
 *   	─4───5───6───7─
 *   	 │   │   │   │
 *   	─0───1───2───3─
 *
 * which means that $P_0 = X_0 X_3 X_4 X_{12}, \cdots$.
 */
auto toric_x_stabilizers_qubits_old(const int Lx, const int Ly, int vertex)
	-> std::set<int>
{
	int row = vertex / Lx;
	int col = vertex % Lx;

	return std::set<int>{2 * row * Lx + col, 2 * row * Lx + (col + Lx),
						 2 * row * Lx + (col - 1 + Lx) % Lx,
						 2 * ((row - 1 + Ly) % Ly) * Lx + col + Lx};
}

auto has_same_elts(std::vector<int> v1, std::vector<int> v2) -> bool
{
	std::set<int> s1(std::make_move_iterator(v1.begin()),
					 std::make_move_iterator(v1.end()));
	std::set<int> s2(std::make_move_iterator(v2.begin()),
					 std::make_move_iterator(v2.end()));

	return s1 == s2;
}

TEST_CASE("Test internal functions", "[internal]")
{
	// Test using Lx=4, Ly=2
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 0) == std::set<int>{0, 3, 4, 12});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 1) == std::set<int>{0, 1, 5, 13});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 2) == std::set<int>{1, 2, 6, 14});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 3) == std::set<int>{2, 3, 7, 15});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 4) == std::set<int>{4, 8, 11, 12});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 5) == std::set<int>{5, 8, 9, 13});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 6) == std::set<int>{6, 9, 10, 14});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 7) == std::set<int>{7, 10, 11, 15});
}
/**
 * LatticeFromParity class only get a list of parity operators (in terms of sites) as an
 * input. It does not know any information in advance, but constructs a graph where
 * vertices are the index of the parity operator and edges are qubits.
 */
TEST_CASE("Test whether the constructed lattice is correct with a toric code",
		  "[LatticeFromParity]")
{
	std::random_device rd;
	std::default_random_engine re{rd()};
	{
		const int Lx = 4;
		const int Ly = 2;

		std::vector<int> col_indices;
		std::vector<int> indptr;

		indptr.push_back(0);

		int n_parity = 0;

		for(int nx = 0; nx < Lx; ++nx)
		{
			for(int ny = 0; ny < Ly; ++ny)
			{
				auto m = toric_x_stabilizers_qubits_old(Lx, Ly, nx * Ly + ny);
				indptr.push_back(indptr.back() + m.size());
				col_indices.insert(col_indices.end(), std::make_move_iterator(m.begin()),
								   std::make_move_iterator(m.end()));
			}
		}

		auto lattice = LatticeFromParity(Lx * Ly, 2 * Lx * Ly, 4 * Lx * Ly,
										 col_indices.data(), indptr.data());
		REQUIRE(has_same_elts(lattice.vertex_connections(0), {1, 3, 4}));
		REQUIRE(has_same_elts(lattice.vertex_connections(1), {0, 2, 5}));
		REQUIRE(has_same_elts(lattice.vertex_connections(2), {1, 3, 6}));
		REQUIRE(has_same_elts(lattice.vertex_connections(3), {0, 2, 7}));
		REQUIRE(has_same_elts(lattice.vertex_connections(4), {0, 5, 7}));
		REQUIRE(has_same_elts(lattice.vertex_connections(5), {1, 4, 6}));
		REQUIRE(has_same_elts(lattice.vertex_connections(6), {2, 5, 7}));
		REQUIRE(has_same_elts(lattice.vertex_connections(7), {3, 4, 6}));
	}

	for(const int L : {7, 15, 31, 63, 127})
	{
		auto to_vertex = [L](int row, int col)
		{
			return ((row + L) % L) * L + ((col + L) % L);
		};

		std::vector<int> col_indices;
		std::vector<int> indptr;

		indptr.push_back(0);
		for(int nx = 0; nx < L; ++nx)
		{
			for(int ny = 0; ny < L; ++ny)
			{
				auto m = toric_x_stabilizers_qubits_old(L, L, nx * L + ny);
				indptr.push_back(indptr.back() + m.size());
				col_indices.insert(col_indices.end(), std::make_move_iterator(m.begin()),
								   std::make_move_iterator(m.end()));
			}
		}
		auto lattice = LatticeFromParity(L * L, 2 * L * L, 4 * L * L, col_indices.data(),
										 indptr.data());

		std::uniform_int_distribution<int> col_dist(0, L - 1);
		for(int row = 0; row < L; ++row)
		{
			auto col = col_dist(re);

			int vertex = row * L + col;
			REQUIRE(lattice.vertex_connection_count(vertex) == 4);

			REQUIRE(has_same_elts(lattice.vertex_connections(vertex),
								  {to_vertex(row - 1, col), to_vertex(row + 1, col),
								   to_vertex(row, col + 1), to_vertex(row, col - 1)}));
			REQUIRE(lattice.edge_idx({to_vertex(row, col), to_vertex(row, col + 1)})
					== 2 * row * L + col);
			REQUIRE(lattice.edge_idx({to_vertex(row, col), to_vertex(row + 1, col)})
					== 2 * row * L + col + L);
		}
	}
}
