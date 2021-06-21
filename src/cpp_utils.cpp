#include <Eigen/Dense>
#include "utility.hpp"
#include "cpp_utils.hpp"


bool logical_error(const int L, const Eigen::ArrayXi& error, ErrorType error_type)
{
	//one may use a counting algorithm for testing logical error
	int sum1 = 0;
	int sum2 = 0;
	switch(error_type)
	{
	case ErrorType::X:
		//need to think in a dual lattice
		for(int u = 0; u < L*L; u += L)
		{
			sum1 += error[u];
		}
		for(int u = L*L; u < L*L+L; ++u)
		{
			sum2 += error[u];
		}
		break;
	case ErrorType::Z:
		for(int u = 0; u < L; ++u)
		{
			sum1 += error[u];
		}
		for(int u = L*L; u < 2*L*L; u += L)
		{
			sum2 += error[u];
		}
		break;
	}

	return (sum1 % 2 == 1 ) || (sum2 % 2 == 1 );
}


std::vector<int> z_error_to_syndrome_x(const int L, const Eigen::ArrayXi& z_error) 
{
	std::vector<int> syndromes_array(L*L, 0u);
	for(int n = 0; n < 2*L*L; ++n)
	{
		if(z_error[n] == 0)
			continue;
		Edge e = to_edge(L, n);
		syndromes_array[e.u] += 1;
		syndromes_array[e.v] += 1;
	}
	for(auto& u : syndromes_array)
		u %= 2;

	return syndromes_array;
}

std::vector<int> x_error_to_syndrome_z(const int L, const Eigen::ArrayXi& x_error) 
{
	std::vector<int> syndromes_array(L*L, 0u);
	for(int n = 0; n < 2*L*L; ++n)
	{
		if(x_error[n] == 0)
			continue;
		Edge e = to_edge(L, n);

		if(is_horizontal(L, e))
		{
			auto u = left(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			syndromes_array[u] += 1;
			syndromes_array[::to_vertex_index(L, row-1, col)] += 1;
		}
		else
		{
			auto u = lower(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			syndromes_array[u] += 1;
			syndromes_array[::to_vertex_index(L, row, col-1)] += 1;
		}

	}
	for(auto& u : syndromes_array)
		u %= 2;
	return syndromes_array;
}
