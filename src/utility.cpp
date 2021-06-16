#include "utility.hpp"


void to_json(nlohmann::json& j, const Edge& e)
{
	j = nlohmann::json{e.u, e.v};
}
void from_json(const nlohmann::json& j, Edge& e)
{
	e = Edge(j[0], j[1]);
}

inline bool is_horizontal(int L, Edge e)
{
	return ((e.v - e.u) == 1) || ((e.v - e.u) == (L-1));
}
inline bool is_vertical(int L, Edge e)
{
	return !is_horizontal(L, e);
}

inline int lower(int L, Edge e) // works only when vertical
{
	if((e.v - e.u) == L)
		return e.u;
	else
		return e.v;
}

inline int upper(int L, Edge e)
{
	if((e.v - e.u) == L)
		return e.v;
	else
		return e.u;
}

inline int left(int L, Edge e) // works only when horizontal
{
	if((e.v - e.u) == 1)
		return e.u;
	else
		return e.v;
}
inline int right(int L, Edge e) // works only when horizontal
{
	if((e.v - e.u) == 1)
		return e.v;
	else
		return e.u;
}

int decoder_edge_to_qubit_idx(const int L, Edge e, ErrorType error_type)
{
	int idx = 0;
	switch(error_type)
	{
	case ErrorType::X:
		//each edge is a qubit in the dual lattice
		if(is_horizontal(L, e))
		{
			auto u = left(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			idx = L*row + ((col+1) % L);
		}
		else
		{
			auto u = upper(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			idx = L*((row+1) % L) + col + L*L;
		}
		break;
	case ErrorType::Z:
		//each edge in correction is a actual qubit
		if(is_horizontal(L, e))
		{
			auto u = left(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			idx = L*row + col + L*L;
		}
		else
		{
			auto u = upper(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			idx = L*row + col;
		}
		break;
	}
	return idx;
}


void add_corrections(const int L, const std::vector<Edge>& corrections, 
		std::vector<int>& error, ErrorType error_type)
{
	for(auto e: corrections)
	{
		auto idx = decoder_edge_to_qubit_idx(L, e, error_type);
		error[idx] += 1;
	}
}

bool logical_error(const int L, const std::vector<int>& error, ErrorType error_type)
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

Edge to_edge(const int L, int edge_index) 
{
	int row = edge_index / L; // % L is done in to_vertex_index
	int col = edge_index % L;
	if((edge_index / (L*L)) == 0) //vertical edge
	{
		return Edge(to_vertex_index(L, row, col), to_vertex_index(L, row-1, col));
	}
	else //horizontal edge
	{
		return Edge(to_vertex_index(L, row, col), to_vertex_index(L, row, col+1));
	}
}


std::vector<int> z_error_to_syndrome_x(const int L, const std::vector<int>& z_error) 
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

std::vector<int> x_error_to_syndrome_z(const int L, const std::vector<int>& x_error) 
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


