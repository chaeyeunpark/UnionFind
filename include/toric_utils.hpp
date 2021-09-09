#pragma once

#include <Eigen/Dense>

#include "utility.hpp"
#include "LatticeCubic.hpp"

/**
 * This file contains functions utility functions for the Toric code.
 */

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


int decoder_edge_to_qubit_idx(const int L, Edge e, ErrorType error_type);

constexpr int to_vertex_index(int L, const int row, const int col)
{
	return (((row + L) % L)) * L + (col + L) % L;
}

inline std::tuple<int, int> vertex_to_coord(const int L, const int vidx)
{
	return std::make_tuple(vidx / L, vidx % L);
}


Edge to_edge(const int L, int edge_index);

std::vector<int>
calc_syndromes(const LatticeCubic& lattice, const Eigen::ArrayXXi& errors, ErrorType error_type);

void add_corrections(const int L, const std::vector<Edge>& corrections, 
		Eigen::ArrayXi& error, ErrorType error_type);


bool logical_error(const int L, const Eigen::ArrayXi& error, ErrorType error_type);

bool has_logical_error(int L, Eigen::ArrayXi& error_total, 
		const std::vector<Edge>& corrections, ErrorType error_type);