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
#include "LatticeFromParity.hpp"

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <stdexcept>

// NOLINTBEGIN(cppcoreguidelines-*)
void free_arr_uint32(void* p)
{
	uint32_t* p_int = reinterpret_cast<uint32_t*>(p);
	delete[] p_int;
}
// NOLINTEND(cppcoreguidelines-*)

namespace py = pybind11;

// NOLINTNEXTLINE(cppcoreguidelines-*)
PYBIND11_MODULE(_union_find_py, m)
{
	using UnionFindFromParity = UnionFindCPP::Decoder<UnionFindCPP::LatticeFromParity>;
	py::class_<UnionFindFromParity>(m, "DecoderFromParity")
		.def(py::init(
			[](int num_parities, int num_qubits,
			   py::array_t<int, py::array::c_style | py::array::forcecast> col_indices,
			   py::array_t<int, py::array::c_style | py::array::forcecast> indptr)
			{
				if(num_parities <= 0)
				{
					throw std::invalid_argument(
						"Number of partiy operators must be larger than 0");
				}
				if(num_qubits <= 0)
				{
					throw std::invalid_argument("Number of qubits must be larger than 0");
				}
				// check the given dimension is correct
				return UnionFindFromParity(static_cast<uint32_t>(num_parities),
										   static_cast<uint32_t>(num_qubits),
										   static_cast<int*>(col_indices.request().ptr),
										   static_cast<int*>(indptr.request().ptr));
			}))
		.def(py::init(
			[](int num_parities, int num_qubits,
			   py::array_t<int, py::array::c_style | py::array::forcecast> col_indices,
			   py::array_t<int, py::array::c_style | py::array::forcecast> indptr,
			   int repetitions)
			{
				if(num_parities <= 0)
				{
					throw std::invalid_argument(
						"Number of partiy operators must be larger than 0");
				}
				if(num_qubits <= 0)
				{
					throw std::invalid_argument("Number of qubits must be larger than 0");
				}
				if(repetitions <= 1)
				{
					throw std::invalid_argument("Repetitions must be larger than 1");
				}
				// check the given dimension is correct
				return UnionFindFromParity(static_cast<uint32_t>(num_parities),
										   static_cast<uint32_t>(num_qubits),
										   static_cast<int*>(col_indices.request().ptr),
										   static_cast<int*>(indptr.request().ptr),
										   static_cast<uint32_t>(repetitions));
			}))
		.def("clear", &UnionFindFromParity::clear, "Clear decoder's internal data")
		.def_property_readonly("num_edges", &UnionFindFromParity::num_edges,
							   "Get total number of edges (qubits) of the decoder")
		.def_property_readonly(
			"num_vertices", &UnionFindFromParity::num_vertices,
			"Get total number of vertices (parity operators) of the decoder")
		.def(
			"decode",
			[](UnionFindFromParity& decoder,
			   std::vector<uint32_t> syndromes) -> py::array_t<uint32_t>
			{
				if(decoder.num_vertices() != syndromes.size())
				{
					throw std::invalid_argument("Size of syndromes should be the same as "
												"the size of vertices");
				}
				auto res = decoder.decode(syndromes);

				uint32_t* corrections = new uint32_t[decoder.num_edges()];
				memset(corrections, 0, sizeof(uint32_t) * decoder.num_edges());
				for(size_t i = 0; i < res.size(); ++i)
				{
					corrections[decoder.edge_idx(res[i])] = 1;
				}
				py::capsule free_when_done(corrections, free_arr_uint32);

				return py::array_t<uint32_t>({(int64_t)decoder.num_edges()}, corrections,
											 free_when_done);
			},
			"Decode the given syndroms");
}
