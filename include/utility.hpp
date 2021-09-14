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

template <>
struct std::hash<Edge>
{
	std::size_t operator()(const Edge& e) const noexcept
	{
		auto h1 = std::hash<int>()(e.u);
		auto h2 = std::hash<int>()(e.v);
		return h1 ^ (h2 << 1);
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

inline int lower(int L, Edge e) // works only when vertical
{
	if((e.v - e.u) == L)
		return e.u;
	else
		return e.v;
}

inline int upper(int L, Edge e)
{
	if((e.v - e.u) == L)
		return e.v;
	else
		return e.u;
}

inline int left(int L, Edge e) // works only when horizontal
{
	if((e.v - e.u) == 1)
		return e.u;
	else
		return e.v;
}
inline int right(int L, Edge e) // works only when horizontal
{
	if((e.v - e.u) == 1)
		return e.v;
	else
		return e.u;
}


