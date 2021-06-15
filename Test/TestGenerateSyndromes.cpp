#include <random>
#include <nlohmann/json.hpp>
#include <iostream>

#include "ErrorGenerator.hpp"
#include "utility.hpp"

int main()
{
	std::random_device rd;
	ErrorGenerator<std::default_random_engine> gen;

	const uint32_t L = 11;

	auto [x_errors, z_errors] = gen.get_errors(L, 0.06, NoiseType::Depolarizing);

	auto synd_x = errors_to_syndromes(L, x_errors, ErrorType::X);
	auto synd_z = errors_to_syndromes(L, z_errors, ErrorType::Z);

	nlohmann::json err_j = {};
	err_j["L"] = L;
	err_j["x_errors"] = x_errors;
	err_j["z_errors"] = z_errors;

	err_j["syndromes_for_x"] = synd_x;
	err_j["syndromes_for_z"] = synd_z;

	std::cout << err_j << std::endl;

	return 0;
}
