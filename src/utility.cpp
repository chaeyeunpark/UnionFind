#include "utility.hpp"


void to_json(nlohmann::json& j, const Edge& e)
{
	j = nlohmann::json{e.u, e.v};
}
void from_json(const nlohmann::json& j, Edge& e)
{
	e = Edge(j[0], j[1]);
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
			idx = L*((row+1) % L) + ((col+1) % L);
		}
		else
		{
			auto u = upper(L, e);
			const auto [row, col] = vertex_to_coord(L, u);
			idx = L*(row % L) + col + L*L;
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

