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
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "typedefs.hpp"
#include "utility.hpp"

/**
 * This file contains <em>lattice independent</em> functions for generating
 * qubit and syndrome measurement errors.
 */

namespace UnionFindCPP
{
enum class NoiseType
{
	Depolarizing = 0,
	Independent,
	X,
	Z
};

template<class RandomEngine>
auto create_errors(RandomEngine& re, const uint32_t num_qubits, const double p,
				   NoiseType noise_type) -> std::pair<ArrayXu, ArrayXu>
{
	ArrayXu x_error = ArrayXu::Zero(num_qubits);
	ArrayXu z_error = ArrayXu::Zero(num_qubits);
	std::uniform_real_distribution<> urd(0.0, 1.0);

	switch(noise_type)
	{
	case NoiseType::Depolarizing:
	{
		std::uniform_int_distribution<> uid(0, 2);
		for(int i = 0; i < num_qubits; ++i)
		{
			if(!(urd(re) < p)) { continue; }
			switch(uid(re))
			{
			case 0:
				x_error[i] = 1;
				break;
			case 1:
				z_error[i] = 1;
				break;
			case 2:
				x_error[i] = 1;
				z_error[i] = 1;
				break;
			}
		}
	}
	break;

	case NoiseType::Independent:
		for(int i = 0; i < num_qubits; i++) { x_error[i] = (urd(re) < p) ? 1 : 0; }
		for(int i = 0; i < num_qubits; i++) { z_error[i] = (urd(re) < p) ? 1 : 0; }
		break;

	case NoiseType::X:
		for(int i = 0; i < num_qubits; i++) { x_error[i] = (urd(re) < p) ? 1 : 0; }
		break;

	case NoiseType::Z:
		for(int i = 0; i < num_qubits; i++) { z_error[i] = (urd(re) < p) ? 1 : 0; }
		break;
	}

	return std::make_pair(x_error, z_error);
}

template<class RandomEngine>
auto create_measurement_errors(RandomEngine& re, const uint32_t num_qubits,
							   const uint32_t repetition, const double p,
							   NoiseType noise_type)
{
	ArrayXXu measurement_error_x = ArrayXXu::Zero(num_qubits, repetition);
	ArrayXXu measurement_error_z = ArrayXXu::Zero(num_qubits, repetition);
	std::uniform_real_distribution<> urd(0.0, 1.0);

	for(int h = 0; h < repetition - 1; ++h) // perfect measurement in the last round
	{
		const auto [layer_error_x, layer_error_z]
			= create_errors(re, num_qubits, p, noise_type);
		measurement_error_x.col(h) = layer_error_x;
		measurement_error_z.col(h) = layer_error_z;
	}

	return std::make_pair(measurement_error_x, measurement_error_z);
}

template<class RandomEngine>
auto generate_errors(const uint32_t num_qubits, const uint32_t repetition, double p,
					 RandomEngine& re, NoiseType noise_type)
{
	ArrayXXu qubit_errors_x = ArrayXXu::Zero(num_qubits, repetition);
	ArrayXXu qubit_errors_z = ArrayXXu::Zero(num_qubits, repetition);

	for(int r = 0; r < repetition; ++r)
	{
		auto [layer_error_x, layer_error_z]
			= create_errors(re, num_qubits, p, noise_type);
		qubit_errors_x.col(r) = layer_error_x;
		qubit_errors_z.col(r) = layer_error_z;
	}

	// cumsum
	for(int h = 1; h < repetition; ++h)
	{
		qubit_errors_x.col(h) += qubit_errors_x.col(h - 1);
		qubit_errors_z.col(h) += qubit_errors_z.col(h - 1);
	}
	qubit_errors_x = qubit_errors_x.unaryExpr([](uint32_t x) { return x % 2; });
	qubit_errors_z = qubit_errors_z.unaryExpr([](uint32_t x) { return x % 2; });

	return std::make_pair(qubit_errors_x, qubit_errors_z);
}

void add_measurement_noise(uint32_t L, std::vector<uint32_t>& syndromes,
						   const ArrayXXu& measurement_error);

void layer_syndrome_diff(uint32_t L, std::vector<uint32_t>& syndromes);

} // namespace UnionFindCPP
