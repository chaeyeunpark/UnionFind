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

#include "typedefs.hpp"
#include "utility.hpp"

#include <Eigen/Dense>

namespace UnionFindCPP
{
inline auto is_horizontal(uint32_t L, Edge e) -> bool
{
	return ((e.v - e.u) == 1) || ((e.v - e.u) == (L - 1));
}
inline auto is_vertical(uint32_t L, Edge e) -> bool
{
	return !is_horizontal(L, e);
}
inline auto lower(uint32_t L, Edge e) -> uint32_t // works only when vertical
{
	if((e.v - e.u) == L) { return e.u; }
	// else
	return e.v;
}
inline auto upper(uint32_t L, Edge e) -> uint32_t
{
	if((e.v - e.u) == L) { return e.v; }
	// else
	return e.u;
}

inline auto left(uint32_t /*L*/, Edge e) -> uint32_t // works only when horizontal
{
	if((e.v - e.u) == 1) { return e.u; }
	// else
	return e.v;
}
inline auto right(uint32_t /*L*/, Edge e) -> uint32_t // works only when horizontal
{
	if((e.v - e.u) == 1) { return e.v; }
	// else
	return e.u;
}

/**
 * This file contains functions utility functions for the Toric code.
 */

auto z_error_to_syndrome_x(uint32_t L, const ArrayXu& z_error) -> std::vector<uint32_t>;
auto x_error_to_syndrome_z(uint32_t L, const ArrayXu& x_error) -> std::vector<uint32_t>;

/*
 * if error_type is X, the output is error locations in the dual lattice.
 * @param	error	An array of length 2*L*L where element 1 indicates the error at the
 * given index
 * */
inline auto errors_to_syndromes(const uint32_t L, const ArrayXu& error,
								ErrorType error_type) -> std::vector<uint32_t>
{
	switch(error_type)
	{
	case ErrorType::X:
		return x_error_to_syndrome_z(L, error);
	case ErrorType::Z:
		return z_error_to_syndrome_x(L, error);
	}
	return {};
}

constexpr auto to_vertex_index(uint32_t L, const uint32_t row, const uint32_t col)
	-> uint32_t
{
	return (((row + L) % L)) * L + (col + L) % L;
}

inline auto vertex_to_coord(const uint32_t L, const uint32_t vidx)
	-> std::tuple<uint32_t, uint32_t>
{
	return std::make_tuple(vidx / L, vidx % L);
}

inline auto to_edge_idx(const uint32_t L, const Edge e) -> uint32_t
{
	if(is_horizontal(L, e))
	{
		auto u = left(L, e);
		const auto [row, col] = vertex_to_coord(L, u);
		return L * row + col + L * L;
	}
	// else
	{
		auto u = upper(L, e);
		const auto [row, col] = vertex_to_coord(L, u);
		return L * row + col;
	}
}

auto decoder_edge_to_qubit_idx(uint32_t L, Edge e, ErrorType error_type) -> uint32_t;

auto to_edge(uint32_t L, uint32_t edge_index) -> Edge;

class LatticeCubic; // forward def

auto calc_syndromes(const LatticeCubic& lattice, const ArrayXXu& errors,
					ErrorType error_type) -> std::vector<uint32_t>;

void add_corrections(uint32_t L, const std::vector<Edge>& corrections, ArrayXu& error,
					 ErrorType error_type);

auto logical_error(uint32_t L, const ArrayXu& error, ErrorType error_type) -> bool;

auto has_logical_error(uint32_t L, ArrayXu& error_total,
					   const std::vector<Edge>& corrections, ErrorType error_type)
	-> bool;
} // namespace UnionFindCPP
