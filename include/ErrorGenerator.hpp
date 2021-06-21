#pragma once
#include <vector>
#include <cstdint>
#include <random>
#include <Eigen/Dense>

#include "utility.hpp"



class ErrorGenerator
{
private:
	const int L_;
public:

	explicit ErrorGenerator(int L)
		: L_{L}
	{
	}


};

