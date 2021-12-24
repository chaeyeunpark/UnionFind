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
#include "utility.hpp"

#include <ostream>
namespace UnionFindCPP
{
void to_json(nlohmann::json& j, const Edge& e)
{
	j = nlohmann::json{e.u, e.v};
}
void from_json(const nlohmann::json& j, Edge& e)
{
	e = Edge(j[0], j[1]);
}

auto operator<<(std::ostream& os, const Edge& e) -> std::ostream&
{
	os << "[" << e.u << ", " << e.v << "]";
	return os;
}
} // namespace UnionFindCPP
