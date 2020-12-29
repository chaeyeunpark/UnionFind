#pragma once
#include <cstdint>
#include <algorithm>
#include <functional>
#include <vector>

#include <nlohmann/json.hpp>
#include <tsl/robin_hash.h>

#ifndef NDEBUG
#define DEBUG
#endif

enum class ErrorType
{
	X, Z
};

struct Edge
{
	int u;
	int v;

	Edge(int ul, int vl)
	{
		u = std::min(ul,vl);
		v = std::max(ul,vl);
	}

	inline bool operator==(const Edge& rhs) const
	{
		return (u == rhs.u) && (v == rhs.v);
	}
};



void to_json(nlohmann::json& j, const Edge& e);
void from_json(const nlohmann::json& j, Edge& e);


inline bool is_horizontal(int L, Edge e)
{
	return ((e.v - e.u) == 1) || ((e.v - e.u) == (L-1));
}
inline bool is_vertical(int L, Edge e)
{
	return !is_horizontal(L, e);
}

int lower(int L, Edge e); // works only when vertical
int left(int L, Edge e); // works only when horizontal

constexpr int to_vertex_index(int L, const int row, const int col)
{
	return (((row + L) % L)) * L + (col + L) % L;
}

inline std::tuple<int, int> vertex_to_coord(const int L, const int vidx)
{
	return std::make_tuple(vidx / L, vidx % L);
}

Edge to_edge(int L, int edge_index);
std::vector<int> syndrome_locations(const int L, 
		const std::vector<int>& syndromes_array);


/*
 * if error_type is X, the output is error locations in the dual lattice.
 * */
std::vector<int> errors_to_syndromes(const int L, const std::vector<int>& error,
		ErrorType error_type);

int decoder_edge_to_qubit_idx(const int L, Edge e, ErrorType error_type);

void add_corrections(const int L, const std::vector<Edge>& corrections, 
		std::vector<int>& error, ErrorType error_type);

bool logical_error(const int L, const std::vector<int>& error, ErrorType error_type);
