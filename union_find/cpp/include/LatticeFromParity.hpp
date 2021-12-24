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
	uint32_t num_vertices_;
	uint32_t num_edges_;

	/* connectivity of vertices (parities). Length num_parities */
	std::vector<std::vector<uint32_t>> vertex_connections_;

	tsl::robin_map<Edge, uint32_t> edge_idx_;

	static auto construct_qubit_associated_parities(uint32_t num_parities,
													uint32_t num_qubits, int* col_indices,
													const int* indptr)
		-> std::vector<std::vector<uint32_t>>
	{
		/**
		 * key: qubit_idx, value: vector of all connected verticies (parities). Length
		 * num_qubits
		 */
		std::vector<std::vector<uint32_t>> qubit_associated_parities;
		qubit_associated_parities.resize(num_qubits);
		for(uint32_t p_idx = 0; p_idx < num_parities; ++p_idx)
		{
			// NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
			for(uint32_t idx = indptr[p_idx]; idx < indptr[p_idx + 1]; ++idx)
			{
				// NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
				qubit_associated_parities[col_indices[idx]].push_back(p_idx);
			}
		}
		return qubit_associated_parities;
	}

	static auto construct_edge_idx(
		uint32_t num_qubits,
		const std::vector<std::vector<uint32_t>>& qubit_associated_parities)
		-> tsl::robin_map<Edge, uint32_t>
	{
		tsl::robin_map<Edge, uint32_t> edge_idx;
		for(uint32_t q_idx = 0; q_idx < num_qubits; ++q_idx)
		{
			const auto& q_parities = qubit_associated_parities[q_idx];
			auto edge = Edge{q_parities[0], q_parities[1]};
			auto iter = edge_idx.find(edge);

			if(iter == edge_idx.end()) // first appear
			{
				edge_idx.emplace(edge, q_idx);
			}
		}
		return edge_idx;
	}

	/**
	 * @brief Use edge_idx_ to construct vertex_connections_
	 */
	void construct_vertex_connections_from_edges()
	{
		vertex_connections_.resize(num_vertices_);
		for(const auto& [edge, qidx] : edge_idx_)
		{
			const auto [u, v] = edge;
			vertex_connections_[u].emplace_back(v);
			vertex_connections_[v].emplace_back(u);
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
	LatticeFromParity(uint32_t num_parities, uint32_t num_qubits, int* col_indices,
					  int* indptr)
		: num_vertices_{num_parities}, num_edges_{num_qubits}
	{
		auto qubit_associated_parities = construct_qubit_associated_parities(
			num_parities, num_qubits, col_indices, indptr);

		edge_idx_ = construct_edge_idx(num_qubits, qubit_associated_parities);
		construct_vertex_connections_from_edges();
	}

	LatticeFromParity(uint32_t layer_vertex_size, uint32_t layer_num_qubits,
					  int* col_indices, int* indptr, uint32_t repetitions)
		: num_vertices_{layer_vertex_size * repetitions},
		  num_edges_{layer_num_qubits * repetitions
					 + layer_vertex_size * (repetitions - 1)}
	{
		if(repetitions < 2)
		{
			throw std::invalid_argument("Repetition must be greater than or equal to 2.");
		}

		auto qubit_associated_parities = construct_qubit_associated_parities(
			layer_vertex_size, layer_num_qubits, col_indices, indptr);
		auto layer_edge_idx
			= construct_edge_idx(layer_num_qubits, qubit_associated_parities);

		// Construct edge_idx_
		edge_idx_.reserve(num_edges_);
		// Add spacelike edges
		for(uint32_t depth = 0; depth < repetitions; ++depth)
		{
			for(const auto& [layer_edge, q_idx] : layer_edge_idx)
			{
				auto edge = Edge{layer_edge.u + depth * layer_vertex_size,
								 layer_edge.v + depth * layer_vertex_size};

				edge_idx_.emplace(edge,
								  q_idx + depth * (layer_vertex_size + layer_num_qubits));
			}
		}

		for(uint32_t depth = 0; depth < repetitions - 1; ++depth)
		{
			// Add timelike edges
			for(int layer_vertex_idx = 0; layer_vertex_idx < layer_vertex_size;
				++layer_vertex_idx)
			{
				auto edge = Edge{depth * layer_vertex_size + layer_vertex_idx,
								 (depth + 1) * layer_vertex_size + layer_vertex_idx};
				edge_idx_.emplace(edge,
								  layer_vertex_idx + layer_num_qubits * (depth + 1));
			}
		}

		// Construct vertex_connections_
		construct_vertex_connections_from_edges();
	}

	[[nodiscard]] auto vertex_connections(uint32_t v) const
		-> const std::vector<uint32_t>&
	{
		return vertex_connections_[v];
	}

	[[nodiscard]] auto vertex_connection_count(int vertex) const -> uint32_t
	{
		return static_cast<int>(vertex_connections_[vertex].size());
	}

	[[nodiscard]] inline auto edge_idx(const Edge& edge) const -> uint32_t
	{
		return edge_idx_.at(edge);
	}

	[[nodiscard]] auto num_edges() const -> uint32_t { return num_edges_; }

	[[nodiscard]] auto num_vertices() const -> uint32_t { return num_vertices_; }

	[[nodiscard]] inline auto
	edge_idx_all() const& -> const tsl::robin_map<Edge, uint32_t>&
	{
		return edge_idx_;
	}
	[[nodiscard]] inline auto edge_idx_all() && -> tsl::robin_map<Edge, uint32_t>
	{
		return edge_idx_;
	}
};
} // namespace UnionFindCPP
