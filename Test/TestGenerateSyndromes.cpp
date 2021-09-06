#include <random>
#include <nlohmann/json.hpp>
#include <iostream>

#include "error_utils.hpp"
#include "utility.hpp"

int main()
{
	std::random_device rd;
	std::default_random_engine re{rd()};

	const uint32_t L = 5;

	auto [x_errors, z_errors] = create_errors(re, 2*L*L, 0.06, NoiseType::Depolarizing);

	auto synd_x = errors_to_syndromes(L, x_errors, ErrorType::X);
	auto synd_z = errors_to_syndromes(L, z_errors, ErrorType::Z);

	nlohmann::json err_j = {};
	err_j["L"] = L;
	err_j["x_errors"] = std::vector<int>(x_errors.data(), x_errors.data() + 2*L*L);
	err_j["z_errors"] = std::vector<int>(z_errors.data(), z_errors.data() + 2*L*L);

	err_j["syndromes_for_x"] = synd_x;
	err_j["syndromes_for_z"] = synd_z;

	std::cout << err_j << std::endl;

	return 0;
}
