#pragma once
#include <vector>
#include <cstdint>
#include <random>

#include "utility.hpp"

enum class NoiseType
{
	Depolarizing = 0, Independent, X, Z
};

template<class RandomEngine>
class ErrorGenerator
{
private:
	RandomEngine re_;

public:
	ErrorGenerator()
	{
		std::random_device rd;
		re_.seed(rd());
	}

	auto get_errors(int L, const double p, NoiseType noise_type = NoiseType::Depolarizing)
	{
		std::vector<int> x_error(2*L*L, 0u);
		std::vector<int> z_error(2*L*L, 0u);
		std::uniform_real_distribution<> urd(0.0,1.0);

		switch(noise_type)
		{
		case NoiseType::Depolarizing:
			{
			std::uniform_int_distribution<> uid(0, 2);
			for(int i = 0; i < 2*L*L; ++i)
			{
				if(!(urd(re_) < p))
					continue;
				switch(uid(re_))
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
			for(int i = 0; i < 2*L*L; i++)
			{
				x_error[i] = (urd(re_) < p)?1:0;
			}
			for(int i = 0; i < 2*L*L; i++)
			{
				z_error[i] = (urd(re_) < p)?1:0;
			}
			break;

		case NoiseType::X:
			for(int i = 0; i < 2*L*L; i++)
			{
				x_error[i] = (urd(re_) < p)?1:0;
			}
			break;

		case NoiseType::Z:
			for(int i = 0; i < 2*L*L; i++)
			{
				z_error[i] = (urd(re_) < p)?1:0;
			}
			break;
		}

		return std::make_pair(x_error, z_error);
	}
};

