#pragma once
#include "utility.hpp"

#include <array>
#include <concepts>
#include <vector>

namespace detail 
{
template<typename T>
struct is_std_array: std::false_type {};

template<typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template<typename T>
concept std_array = is_std_array<T>::value;
} // namespace detail

template<typename T>
concept vertex_connections_result =
    std::same_as<T, std::vector<int>> || detail::std_array<T>;

template <typename T>
concept LatticeConcept = requires (const T lattice, int vertex, Edge e)
{
	{lattice.num_vertices()} -> std::convertible_to<int>;
	{lattice.num_edges()} -> std::convertible_to<int>;
	{lattice.vertex_connections(vertex)} -> vertex_connections_result;
	{lattice.vertex_connection_count(vertex)} -> std::convertible_to<int>;
	{lattice.edge_idx(e)} -> std::convertible_to<int>;
};
