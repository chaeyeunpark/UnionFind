// Copyright (C) 2021 UnionFind++ authors
//
// This file is part of UnionFind++.
// 
// UnionFind++ is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// UnionFind++ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with UnionFind++.  If not, see <https://www.gnu.org/licenses/>.

#include "Decoder.hpp"
#include "LatticFromParity.hpp"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"

#include <stdexcept>

void free_int_arr(void* p)
{
	int* p_int = reinterpret_cast<int*>(p);
	delete[] p_int;
}

namespace py = pybind11;

PYBIND11_MODULE(_union_find_py, m)
{
using UnionFindFromParity = UnionFindCPP::Decoder<UnionFindCPP::LatticeFromParity>;
py::class_<UnionFindFromParity>(m, "DecoderFromParity")
	.def(py::init(
		[](int num_qubits, int num_parity, int nnz, py::array_t<int> col_indices,
			py::array_t<int> indptr)
	{
		// check the given dimension is correct
		return UnionFindFromParity(num_qubits, num_parity, nnz,
				static_cast<int*>(col_indices.request().ptr), 
				static_cast<int*>(indptr.request().ptr));
	}))
	.def("clear", &UnionFindFromParity::clear, "Clear decoder's internal data")
	.def_property_readonly("num_edges", &UnionFindFromParity::num_edges, "Get total number of edges (qubits) of the decoder")
	.def_property_readonly("num_vertices", &UnionFindFromParity::num_vertices, "Get total number of vertices (parity operators) of the decoder")
	.def("decode", 
		[](UnionFindFromParity& decoder, std::vector<int>& syndromes) -> py::array_t<int>
	{
		if(decoder.num_vertices() != syndromes.size())
		{
			throw std::invalid_argument("Size of syndromes should be the same as "
					"the size of vertices");
		}
		auto res = decoder.decode(syndromes);

		int* corrections = new int[decoder.num_edges()];
		memset(corrections, 0, sizeof(int)*decoder.num_edges());
		for(size_t i = 0; i < res.size(); ++i)
		{
			corrections[decoder.edge_idx(res[i])] = 1;
		}
		py::capsule free_when_done(corrections, free_int_arr);

		return py::array_t<int>(
			{(int64_t)decoder.num_edges()},
			corrections, free_when_done);
	}, "Decode the given syndroms");
}
