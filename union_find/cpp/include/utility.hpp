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
//
#pragma once
#include <algorithm>
#include <cstdint>
#include <functional>
#include <ostream>
#include <vector>

#include <nlohmann/json.hpp>

#ifndef NDEBUG
#define DEBUG
#endif

namespace UnionFindCPP
{

enum class ErrorType
{
	X,
	Z
};

struct Edge
{
	int u;
	int v;

	Edge(int ul, int vl)
	{
		u = std::min(ul, vl);
		v = std::max(ul, vl);
	}

	inline bool operator==(const Edge& rhs) const { return (u == rhs.u) && (v == rhs.v); }
};

void to_json(nlohmann::json& j, const Edge& e);
void from_json(const nlohmann::json& j, Edge& e);

inline bool is_horizontal(int L, Edge e)
{
	return ((e.v - e.u) == 1) || ((e.v - e.u) == (L - 1));
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

std::ostream& operator<<(std::ostream&, const UnionFindCPP::Edge& e);
} // namespace UnionFindCPP

template<> struct std::hash<UnionFindCPP::Edge>
{
	std::size_t operator()(const UnionFindCPP::Edge& e) const noexcept
	{
		auto h1 = std::hash<int>()(e.u);
		auto h2 = std::hash<int>()(e.v);
		return h1 ^ (h2 << 1);
	}
};
