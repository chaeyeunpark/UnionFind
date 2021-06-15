Installation
===============

.. highlight:: console

You can install this project by clonning the source code tree and compiling it.
To clone the source tree, use `git clone` as ::

	git clone https://github.com/chaeyeunpark/UnionFind.git
	cd UnionFind
	git submodule update --init --recursive

Then you can compile the code using ::

	mkdir build && cd build
	cmake ..
	make union_find

.. highlight:: python

The compiled library file can be found under name `union_find.cpython-[some extra string].so`. You can now use it ::

	from union_find import UnionFindToric
	decoder = UnionFindToric(lattice_size) # L = lattice_size
	decoder.decode(syndromes) # syndromes is a list of size L^2
	decoder.clear() # should be called for reuse

This project does not support PyPI now. 
