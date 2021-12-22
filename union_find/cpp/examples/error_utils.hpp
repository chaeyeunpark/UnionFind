#pragma once
#include <vector>
#include <random>

#include <Eigen/Dense>

#include "utility.hpp"

/**
 * This file contains <em>lattice independent</em> functions for generating
 * qubit and syndrome measurement errors.
 */

enum class NoiseType
{
	Depolarizing = 0, Independent, X, Z
};


template<class RandomEngine>
auto create_errors(RandomEngine& re, const int num_qubits,
		const double p, NoiseType noise_type)
{
	Eigen::ArrayXi x_error = Eigen::ArrayXi::Zero(num_qubits);
	Eigen::ArrayXi z_error = Eigen::ArrayXi::Zero(num_qubits);
	std::uniform_real_distribution<> urd(0.0,1.0);

	switch(noise_type)
	{
	case NoiseType::Depolarizing:
		{
		std::uniform_int_distribution<> uid(0, 2);
		for(int i = 0; i < num_qubits; ++i)
		{
			if(!(urd(re) < p))
				continue;
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
		for(int i = 0; i < num_qubits; i++)
		{
			x_error[i] = (urd(re) < p)?1:0;
		}
		for(int i = 0; i < num_qubits; i++)
		{
			z_error[i] = (urd(re) < p)?1:0;
		}
		break;

	case NoiseType::X:
		for(int i = 0; i < num_qubits; i++)
		{
			x_error[i] = (urd(re) < p)?1:0;
		}
		break;

	case NoiseType::Z:
		for(int i = 0; i < num_qubits; i++)
		{
			z_error[i] = (urd(re) < p)?1:0;
		}
		break;
	}

	return std::make_pair(x_error, z_error);
}

template<class RandomEngine>
auto create_measurement_errors(RandomEngine& re, const int num_qubits, 
		const int repetition, const double p, NoiseType noise_type)
{
	Eigen::ArrayXXi measurement_error_x = Eigen::ArrayXXi::Zero(num_qubits, repetition);
	Eigen::ArrayXXi measurement_error_z = Eigen::ArrayXXi::Zero(num_qubits, repetition);
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
std::pair<Eigen::ArrayXXi, Eigen::ArrayXXi>
generate_errors(const int num_qubits, const int repetition,	double p, 
		RandomEngine& re, NoiseType noise_type)
{
	Eigen::ArrayXXi qubit_errors_x = Eigen::ArrayXXi::Zero(num_qubits, repetition);
	Eigen::ArrayXXi qubit_errors_z = Eigen::ArrayXXi::Zero(num_qubits, repetition);

	for(int r = 0; r < repetition; ++r)
	{
		auto [layer_error_x, layer_error_z] = create_errors(re, num_qubits, p, noise_type);
		qubit_errors_x.col(r) = layer_error_x;
		qubit_errors_z.col(r) = layer_error_z;
	}

	//cumsum
	for(int h = 1; h < repetition; ++h)
	{
		qubit_errors_x.col(h) += qubit_errors_x.col(h-1);
		qubit_errors_z.col(h) += qubit_errors_z.col(h-1);
	}
	qubit_errors_x = qubit_errors_x.unaryExpr([](int x){ return x % 2;});
	qubit_errors_z = qubit_errors_z.unaryExpr([](int x){ return x % 2;});

	return std::make_pair(qubit_errors_x, qubit_errors_z);
}


template<typename RandomEngine>
void add_measurement_noise(const int L, RandomEngine& re, 
		std::vector<int>& syndromes, const Eigen::ArrayXXi& measurement_error)
{
	Eigen::Map<Eigen::ArrayXXi> syndromes_map(syndromes.data(), L*L, L);
	syndromes_map += measurement_error;
}

void layer_syndrome_diff(const int L, std::vector<int>& syndromes);
