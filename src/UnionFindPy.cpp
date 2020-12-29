#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "UnionFind.hpp"

namespace py = pybind11;

PYBIND11_MODULE(union_find, m)
{
	py::class_<UnionFindDecoder>(m, "UnionFindDecoder")
		.def(py::init<int>())
		.def("clear", &UnionFindDecoder::clear)
		.def("decode", [](UnionFindDecoder& decoder, std::vector<int>& syndromes)
		{
			std::vector<std::pair<int, int>> corrections;
			auto res = decoder.decode(syndromes);
			for(auto edge: res)
			{
				corrections.emplace_back(edge.u, edge.v);
			}
			return corrections;
		});
}
