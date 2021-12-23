Implementing Custom Lattice
============================

Even though one can use our python binding for any possible lattices, implementing a lattice directly in C++ can boost decoding speeds.

As our ``UnionFind`` is a template class, you can add your own custom lattice without much difficulty. 
A custom lattice class should implement 5 methods (implementing ``LatticeConcept``).

.. code-block:: c++

    class CustomLattice
    {
    public:
        using Vertex = int;

        int num_vertices(); //return number of all vertices in the lattice
        int num_edges(); //return number of all edges (qubits) in the lattice

        std::vector<int> vertex_connections(Vertex v);  //return nearest neighbor vertices

        int edge_idx(Edge edge); //return index of edge for a given edge
        int vertex_connection_count(Vertex v); //return the number of nearest neighbor vertices
    };

Then you can use our ``UnionFind`` template class in your C++ code as

.. code-block:: c++

	#include <UnionFind.hpp>
	auto decoder = UnionFind<CustomLattice>(args...);

All parameters of the ``UnionFind`` constructor are perfectly forwarded to the constructor of ``CustomLattice`` class.

We are planning to support a Python interface to generate a custom lattice.


