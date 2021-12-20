#include "UnionFind.hpp"
#include "Lattice2D.hpp"
#include "LatticeCubic.hpp"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"

#include <stdexcept>

namespace py = pybind11;

void free_int_arr(void* p)
{
	int* p_int = reinterpret_cast<int*>(p);
	delete[] p_int;
}

PYBIND11_MODULE(union_find_py, m)
{
using UnionFindToric = UnionFindDecoder<Lattice2D>;
py::class_<UnionFindToric>(m, "UnionFindToric")
	.def(py::init(
		[](int L)
	{
		return UnionFindToric(L);
	}))
	.def("clear", &UnionFindToric::clear)
	.def_property_readonly("num_edges", &UnionFindToric::num_edges)
	.def_property_readonly("num_vertices", &UnionFindToric::num_vertices)
	.def("decode", 
		[](UnionFindToric& decoder, std::vector<int>& syndromes) -> py::array_t<int>
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
	});

using UnionFind3D = UnionFindDecoder<LatticeCubic>;
py::class_<UnionFind3D>(m, "UnionFind3D")
	.def(py::init(
		[](int L)
	{
		return UnionFind3D(L);
	}))
	.def("clear", &UnionFind3D::clear)
	.def_property_readonly("num_edges", &UnionFind3D::num_edges)
	.def_property_readonly("num_vertices", &UnionFind3D::num_vertices)
	.def("decode", 
		[](UnionFind3D& decoder, std::vector<int>& syndromes) -> py::array_t<int>
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
	});
}
