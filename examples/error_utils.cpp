#include <Eigen/Dense>

#include "error_utils.hpp"
#include "utility.hpp"


void layer_syndrome_diff(const int L, std::vector<int>& syndromes)
{
	Eigen::Map<Eigen::ArrayXXi> syndromes_map(syndromes.data(), L*L, L);
	for(int h = L-1; h >= 1; --h)
	{
		syndromes_map.col(h) -= syndromes_map.col(h-1);
	}
	syndromes_map = syndromes_map.unaryExpr([](int x){ return (x+2)%2;});
}

