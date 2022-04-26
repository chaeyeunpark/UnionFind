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
	= std::convertible_to<T, std::vector<uint32_t>> || detail::std_array<T>;

/**
 * @brief Define LatticeConcept that custom Lattice classes should follow.
 * */
template<typename T>
concept LatticeConcept = requires(const T lattice, uint32_t vertex, Edge e)
{
	{
		lattice.num_vertices()
		} -> std::convertible_to<uint32_t>;
	{
		lattice.num_edges()
		} -> std::convertible_to<uint32_t>;
	{
		lattice.vertex_connections(vertex)
		} -> vertex_connections_result;
	{
		lattice.vertex_connection_count(vertex)
		} -> std::convertible_to<uint32_t>;
	{
		lattice.edge_idx(e)
		} -> std::convertible_to<uint32_t>;
};
} // namespace UnionFindCPP
