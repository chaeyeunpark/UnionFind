#pragma once

#include "utility.hpp"

class Lattice2D
{
private:
	int L_;

public:
	using Vertex = int;


	explicit Lattice2D(int L)
		: L_{L}
	{
	}

	inline int to_vertex_index(int row, int col) const
	{
		return ::to_vertex_index(L_, row, col);
	}

	constexpr int vertex_connection_count(Vertex v) const
	{
		return 4;
	}

	std::array<Vertex, 4> vertex_connections(Vertex v) const
	{
		int row = v / L_;
		int col = v % L_;

		return {
			to_vertex_index(row-1, col),
			to_vertex_index(row+1,col),
			to_vertex_index(row,col-1),
			to_vertex_index(row,col+1),
		};
	}
	
	inline int num_vertices() const
	{
		return L_*L_;
	}

	inline int num_edges() const
	{
		return 2*L_*L_;
	}

	inline int edge_idx(const Edge& edge) const
	{
		return decoder_edge_to_qubit_idx(L_, edge, ErrorType::Z);
	}


	inline Edge to_edge(int edge_index) const
	{
		return ::to_edge(L_, edge_index);
	}
};
