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
#include "../examples/Lattice2D.hpp"
#include "../examples/LatticeCubic.hpp"
#include "LatticeConcept.hpp"
#include "LatticeFromParity.hpp"

#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

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
auto toric_x_stabilizers_qubits_old(const uint32_t Lx, const uint32_t Ly, uint32_t vertex)
	-> std::set<uint32_t>
{
	uint32_t row = vertex / Lx;
	uint32_t col = vertex % Lx;

	return std::set<uint32_t>{2 * row * Lx + col, 2 * row * Lx + (col + Lx),
							  2 * row * Lx + (col - 1 + Lx) % Lx,
							  2 * ((row - 1 + Ly) % Ly) * Lx + col + Lx};
}

template<class ContainerType1, class ContainerType2>
auto have_same_elts(const ContainerType1& c1, const ContainerType2& c2) -> bool
{
	static_assert(std::is_same_v<typename ContainerType1::value_type,
								 typename ContainerType2::value_type>,
				  "Two container types must have the same value_type");
	using T = typename ContainerType1::value_type;
	std::set<T> s1(c1.begin(), c1.end());
	std::set<T> s2(c2.begin(), c2.end());

	return s1 == s2;
}

template<UnionFindCPP::LatticeConcept LatticeType1,
		 UnionFindCPP::LatticeConcept LatticeType2>
void check_same_lattice(const LatticeType1& lattice1, const LatticeType2& lattice2)
{
	REQUIRE(lattice1.num_vertices() == lattice2.num_vertices());
	REQUIRE(lattice1.num_edges() == lattice2.num_edges());

	for(uint32_t v = 0; v < lattice1.num_vertices(); ++v)
	{
		REQUIRE(lattice1.vertex_connection_count(v)
				== lattice2.vertex_connection_count(v));
		REQUIRE(have_same_elts(lattice1.vertex_connections(v),
							   lattice2.vertex_connections(v)));
	}

	for(uint32_t v = 0; v < lattice1.num_vertices(); ++v)
	{
		for(uint32_t v2 : lattice2.vertex_connections(v))
		{
			REQUIRE(lattice1.edge_idx({v, v2}) == lattice2.edge_idx({v, v2}));
		}
	}
}

/**
 * @brief generate a parity matrix for a prepetition code
 */
auto repetition_code(uint32_t L) -> Eigen::SparseMatrix<uint32_t, Eigen::RowMajor>
{
	Eigen::SparseMatrix<uint32_t, Eigen::RowMajor> m(L, L);
	m.reserve(2 * L);
	for(uint32_t l = 0; l < L; ++l)
	{
		m.insert(l, l) = 1;
		m.insert(l, (l + 1) % L) = 1;
	}
	m.makeCompressed();
	return m;
}

auto toric_x_stabilizers_qubits_new(const uint32_t L)
{
	using SpMatu = Eigen::SparseMatrix<uint32_t, Eigen::RowMajor>;
	SpMatu Id(L, L);
	Id.setIdentity();
	SpMatu Hr = repetition_code(L);
	SpMatu HrId = Eigen::kroneckerProduct(Hr, Id);
	SpMatu IdHr = Eigen::kroneckerProduct(Id, Hr.transpose());

	SpMatu H(L * L, 2 * L * L);

	for(int k = 0; k < HrId.outerSize(); ++k)
	{
		for(SpMatu::InnerIterator it(HrId, k); it; ++it)
		{
			H.coeffRef(it.row(), it.col()) = it.value();
		}
	}

	for(int k = 0; k < IdHr.outerSize(); ++k)
	{
		for(SpMatu::InnerIterator it(IdHr, k); it; ++it)
		{
			H.coeffRef(it.row(), it.col() + L * L) = it.value();
		}
	}

	H.makeCompressed();

	return H;
}

TEST_CASE("Test internal functions", "[internal]")
{
	// Test using Lx=4, Ly=2
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 0) == std::set<uint32_t>{0, 3, 4, 12});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 1) == std::set<uint32_t>{0, 1, 5, 13});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 2) == std::set<uint32_t>{1, 2, 6, 14});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 3) == std::set<uint32_t>{2, 3, 7, 15});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 4) == std::set<uint32_t>{4, 8, 11, 12});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 5) == std::set<uint32_t>{5, 8, 9, 13});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 6) == std::set<uint32_t>{6, 9, 10, 14});
	REQUIRE(toric_x_stabilizers_qubits_old(4, 2, 7) == std::set<uint32_t>{7, 10, 11, 15});
}
/**
 * LatticeFromParity class only get a list of parity operators (in terms of sites) as an
 * input. It does not know any information in advance, but constructs a graph where
 * vertices are the index of the parity operator and edges are qubits.
 */
TEST_CASE("Test whether the constructed lattice is correct for a toric code (2D)",
		  "[LatticeFromParity]")
{
	std::random_device rd;
	std::default_random_engine re{rd()};

	SECTION("Lx=4, Ly=2, old indexing")
	{
		const uint32_t Lx = 4;
		const uint32_t Ly = 2;

		std::vector<int> col_indices;
		std::vector<int> indptr;

		indptr.push_back(0);

		for(uint32_t nx = 0; nx < Lx; ++nx)
		{
			for(uint32_t ny = 0; ny < Ly; ++ny)
			{
				auto m = toric_x_stabilizers_qubits_old(Lx, Ly, nx * Ly + ny);
				indptr.push_back(indptr.back() + m.size());
				col_indices.insert(col_indices.end(), std::make_move_iterator(m.begin()),
								   std::make_move_iterator(m.end()));
			}
		}

		auto lattice
			= LatticeFromParity(Lx * Ly, 2 * Lx * Ly, col_indices.data(), indptr.data());
		REQUIRE(have_same_elts(lattice.vertex_connections(0),
							   std::vector<uint32_t>{1, 3, 4}));
		REQUIRE(have_same_elts(lattice.vertex_connections(1),
							   std::vector<uint32_t>{0, 2, 5}));
		REQUIRE(have_same_elts(lattice.vertex_connections(2),
							   std::vector<uint32_t>{1, 3, 6}));
		REQUIRE(have_same_elts(lattice.vertex_connections(3),
							   std::vector<uint32_t>{0, 2, 7}));
		REQUIRE(have_same_elts(lattice.vertex_connections(4),
							   std::vector<uint32_t>{0, 5, 7}));
		REQUIRE(have_same_elts(lattice.vertex_connections(5),
							   std::vector<uint32_t>{1, 4, 6}));
		REQUIRE(have_same_elts(lattice.vertex_connections(6),
							   std::vector<uint32_t>{2, 5, 7}));
		REQUIRE(have_same_elts(lattice.vertex_connections(7),
							   std::vector<uint32_t>{3, 4, 6}));
	}

	SECTION("L=7 to 127, old indexing")
	{
		for(const uint32_t L : {7, 15, 31, 63, 127})
		{
			auto to_vertex = [L](uint32_t row, uint32_t col)
			{
				return ((row + L) % L) * L + ((col + L) % L);
			};

			std::vector<int> col_indices;
			std::vector<int> indptr;

			indptr.push_back(0);
			for(uint32_t nx = 0; nx < L; ++nx)
			{
				for(uint32_t ny = 0; ny < L; ++ny)
				{
					auto m = toric_x_stabilizers_qubits_old(L, L, nx * L + ny);
					indptr.push_back(indptr.back() + m.size());
					col_indices.insert(col_indices.end(),
									   std::make_move_iterator(m.begin()),
									   std::make_move_iterator(m.end()));
				}
			}
			auto lattice
				= LatticeFromParity(L * L, 2 * L * L, col_indices.data(), indptr.data());

			std::uniform_int_distribution<int> col_dist(0, L - 1);
			for(uint32_t row = 0; row < L; ++row)
			{
				auto col = col_dist(re);

				uint32_t vertex = row * L + col;
				REQUIRE(lattice.vertex_connection_count(vertex) == 4);

				REQUIRE(have_same_elts(lattice.vertex_connections(vertex),
									   std::vector<uint32_t>{to_vertex(row - 1, col),
															 to_vertex(row + 1, col),
															 to_vertex(row, col + 1),
															 to_vertex(row, col - 1)}));
				REQUIRE(lattice.edge_idx({to_vertex(row, col), to_vertex(row, col + 1)})
						== 2 * row * L + col);
				REQUIRE(lattice.edge_idx({to_vertex(row, col), to_vertex(row + 1, col)})
						== 2 * row * L + col + L);
			}
		}
	}

	SECTION("L=7 to 127, new indexing")
	{
		using UnionFindCPP::Lattice2D;
		for(uint32_t L : {7, 15, 31, 63})
		{
			auto H = toric_x_stabilizers_qubits_new(L);
			auto lattice_from_H = LatticeFromParity(H.rows(), H.cols(), H.innerIndexPtr(),
													H.outerIndexPtr());
			auto lattice2d = Lattice2D(L);

			check_same_lattice(lattice_from_H, lattice2d);
		}
	}
}

TEST_CASE("Test whether the constructed lattice is correct for a toric code (3D)",
		  "[LatticeFromParity]")
{
	SECTION("L=7 to 127, new indexing")
	{
		using UnionFindCPP::LatticeCubic;
		for(uint32_t L : {3, 7, 15})
		{
			auto H = toric_x_stabilizers_qubits_new(L);
			auto lattice_from_H = LatticeFromParity(H.rows(), H.cols(), H.innerIndexPtr(),
													H.outerIndexPtr(),
													/*repetitions = */ L);
			auto lattice3d = LatticeCubic(L);

			check_same_lattice(lattice_from_H, lattice3d);
		}
	}
}
