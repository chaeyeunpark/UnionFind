#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

#include "utility.hpp"

int main(int argc, char* argv[])
{
	using nlohmann::json;
	std::ifstream fin(argv[1]);

	json err_j;
	fin >> err_j;

	std::vector<int> z_err_corr = err_j["z_corrected"];

	std::cout << logical_error(err_j["L"].get<int>(), z_err_corr, ErrorType::Z)
		<< std::endl;

	return 0;
}
