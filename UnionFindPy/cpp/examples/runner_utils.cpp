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
#include "runner_utils.hpp"

#include <fmt/core.h>
#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <span>
#include <stdexcept>

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
auto parse_args(int argc, const char* const argv[]) -> std::pair<uint32_t, double>
{
	auto args = std::span(argv, size_t(argc));
	if(argc != 3) { throw std::invalid_argument(fmt::format("Usage: {} L p", args[0])); }

	uint32_t L = 0;
	double p = 0.0;

	try
	{
		unsigned long Ll = std::stoul(args[1]);
		if(Ll > std::numeric_limits<uint32_t>::max())
		{
			throw std::invalid_argument("L must be in 32 bit unsigned integer");
		}
		L = Ll;
		p = std::stod(args[2]);
	}
	catch(std::exception& e)
	{
		throw e;
	}

	if(!((p > 0.0) && (p < 1.0)))
	{
		throw std::invalid_argument("Error: p must be in between 0.0 and 1.0\n");
	}
	return std::make_pair(L, p);
}

void save_to_json(uint32_t L, double p, double avg_dur_in_microseconds,
				  double avg_success)
{
	constexpr static int p_precision = 5;
	auto p_format = [](double p) -> long
	{
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
		return std::lround(p * std::pow(10, p_precision));
	};

	std::string filename
		= fmt::format("out_L{:d}_P{:0{}d}.json", L, p_format(p), p_precision + 1);
	std::ofstream out_data(filename);
	nlohmann::json out_j;
	out_j["L"] = L;
	out_j["average_microseconds"] = double(avg_dur_in_microseconds);
	out_j["p"] = p;
	out_j["accuracy"] = double(avg_success);

	out_data << out_j.dump(0);
}
