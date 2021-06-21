#pragma once
#include <vector>
#include <random>

#include <Eigen/Dense>

#include "utility.hpp"

enum class NoiseType
{
	Depolarizing = 0, Independent, X, Z
};

bool logical_error(const int L, const Eigen::ArrayXi& error, ErrorType error_type);

std::vector<int> z_error_to_syndrome_x(const int L, const Eigen::ArrayXi& z_error);
std::vector<int> x_error_to_syndrome_z(const int L, const Eigen::ArrayXi& x_error);

/*
 * if error_type is X, the output is error locations in the dual lattice.
 * @param	error	An array of length 2*L*L where element 1 indicates the error at the given index
 * */
inline std::vector<int> errors_to_syndromes(const int L, 
		const Eigen::ArrayXi& error, ErrorType error_type) 
{
	switch(error_type)
	{
	case ErrorType::X:
		return x_error_to_syndrome_z(L, error);
	case ErrorType::Z:
		return z_error_to_syndrome_x(L, error);
	}
	return {};
}

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
auto create_measurement_errors(RandomEngine& re, 
		const int num_qubits, const int repetition,
		const double p, NoiseType noise_type)
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

