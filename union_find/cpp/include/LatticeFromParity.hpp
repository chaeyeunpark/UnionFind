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

#include "tsl/robin_map.h"

#include <array>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include <iostream>

namespace UnionFindCPP
{
/**
 * This class implements LatticeConcept using the input of sparse parity matrix.
 */
class LatticeFromParity
{
private:
	int num_vertices_;
	int num_edges_;

	/* key: qubit_idx, value: vector of all connected verticies (parities). Length
	 * num_qubits */
	std::vector<std::vector<int>> qubit_parities_;
	/* connectivity of vertices (parities). Length num_parities */
	std::vector<std::vector<int>> vertex_connections_;

	tsl::robin_map<Edge, int> edge_idx_;

	void construct_connections(size_t num_parities, size_t num_qubits,
							   int* col_indices, int* indptr)
	{
		qubit_parities_.resize(num_qubits);
		for(size_t p_idx = 0; p_idx < num_parities; ++p_idx)
		{
			for(int idx = indptr[p_idx]; idx < indptr[p_idx + 1]; ++idx)
			{
				qubit_parities_[col_indices[idx]].push_back(p_idx);
			}
		}

		vertex_connections_.resize(num_parities);
		for(size_t q_idx = 0; q_idx < num_qubits; ++q_idx)
		{
			const auto& q_parities = qubit_parities_[q_idx];
			if(q_parities.size() != 2)
			{
				throw std::invalid_argument(
					"Each qubit must appear exactly twice in the parity.");
			}

			vertex_connections_[q_parities[0]].push_back(q_parities[1]);
			vertex_connections_[q_parities[1]].push_back(q_parities[0]);

			auto edge = Edge{q_parities[0], q_parities[1]};
			auto iter = edge_idx_.find(edge);

			if(iter == edge_idx_.end()) // first appear
			{
				edge_idx_.emplace(edge, q_idx);
			}
		}
	}

public:
	/**
	 * @brief construct a Lattice class from a given parity matrix (CSR format)
	 *
	 * @param num_parities total number of parities. Same as the number of rows of the
	 * matrix.
	 * @param num_qubits total number of qubits. Same as the number of columns of the
	 * matrix
	 * @param col_indices indices[idx] indicate the column index of the element data[idx]
	 * @param indptr indptr[row+1]-indptr[row] indicate the number of elements in the row
	 */
	LatticeFromParity(int num_parities, int num_qubits, int* col_indices,
					  int* indptr)
		: num_vertices_{num_parities}, num_edges_{num_qubits}
	{
		construct_connections(num_parities, num_qubits, col_indices, indptr);
	}

	const std::vector<int>& vertex_connections(int v) const
	{
		return vertex_connections_[v];
	}

	int vertex_connection_count(int vertex) const
	{
		return vertex_connections_[vertex].size();
	}

	inline int edge_idx(const Edge& edge) const { return edge_idx_.at(edge); }

	int num_edges() const { return num_edges_; }

	int num_vertices() const { return num_vertices_; }
};
} // namespace UnionFindCPP
