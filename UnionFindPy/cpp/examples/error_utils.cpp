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
#include <Eigen/Dense>

#include "error_utils.hpp"
#include "utility.hpp"

namespace UnionFindCPP
{
void add_measurement_noise(const uint32_t L, std::vector<uint32_t>& syndromes,
						   const ArrayXXu& measurement_error)
{
	Eigen::Map<ArrayXXu> syndromes_map(syndromes.data(), L * L, L);
	syndromes_map += measurement_error;
}

void layer_syndrome_diff(const uint32_t L, std::vector<uint32_t>& syndromes)
{
	Eigen::Map<ArrayXXu> syndromes_map(syndromes.data(), L * L, L);
	for(int h = L - 1; h >= 1; --h) { syndromes_map.col(h) += syndromes_map.col(h - 1); }
	syndromes_map = syndromes_map.unaryExpr([](uint32_t x) { return x % 2; });
}
} // namespace UnionFindCPP
