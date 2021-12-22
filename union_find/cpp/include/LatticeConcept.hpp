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

#include <array>
#include <concepts>
#include <vector>

namespace UnionFindCPP
{
namespace detail
{
	template<typename T> struct is_std_array : std::false_type
	{
	};

	template<typename T, std::size_t N>
	struct is_std_array<std::array<T, N>> : std::true_type
	{
	};

	template<typename T> concept std_array = is_std_array<T>::value;
} // namespace detail

template<typename T>
concept vertex_connections_result
	= std::convertible_to<T, std::vector<int>> || detail::std_array<T>;

template<typename T>
concept LatticeConcept = requires(const T lattice, int vertex, Edge e)
{
	{
		lattice.num_vertices()
		} -> std::convertible_to<int>;
	{
		lattice.num_edges()
		} -> std::convertible_to<int>;
	{
		lattice.vertex_connections(vertex)
		} -> vertex_connections_result;
	{
		lattice.vertex_connection_count(vertex)
		} -> std::convertible_to<int>;
	{
		lattice.edge_idx(e)
		} -> std::convertible_to<int>;
};
}
