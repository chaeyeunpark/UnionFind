#pragma once

#include "LatticeConcept.hpp"
#include "RootManager.hpp"
#include "utility.hpp"

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <queue>
#include <set>
#include <vector>

namespace UnionFindCPP
{
template<LatticeConcept Lattice> class Decoder
{
public:
	using Vertex = uint32_t;
	using RootIterator = tsl::robin_set<Vertex>::const_iterator;

private:
	const Lattice lattice_;

	/* index: vertex */
	std::vector<Vertex> connection_counts_;

	/* index: edge index */
	std::vector<uint32_t> support_;
	std::deque<Edge> fuse_list_;

	/* index: vertex */
	std::vector<Vertex> root_of_vertex_; // root of vertex

	RootManager mgr_;
	/* key: root, value: borders */
	tsl::robin_map<Vertex, tsl::robin_set<Vertex>> border_vertices_;

	/* Data for peeling */
	std::deque<Edge> peeling_edges_;

	void init_cluster(const std::vector<uint32_t>& roots)
	{
		connection_counts_ = std::vector<Vertex>(lattice_.num_vertices(), 0);
		support_ = std::vector<uint32_t>(lattice_.num_edges(), 0);
		mgr_.initialize_roots(roots);
		for(auto root : roots) { border_vertices_[root].emplace(root); }

		root_of_vertex_.resize(lattice_.num_vertices());

		for(uint32_t u = 0; u < lattice_.num_vertices(); ++u) { root_of_vertex_[u] = u; }
	}

	void grow(Vertex root)
	{
		for(auto border_vertex : border_vertices_[root])
		{
			for(auto v : lattice_.vertex_connections(border_vertex))
			{
				auto edge = Edge(border_vertex, v);

				auto& elt = support_[lattice_.edge_idx(edge)];
				if(elt == 2) { continue; }
				if(++elt == 2)
				{
					connection_counts_[edge.u]++;
					connection_counts_[edge.v]++;
					fuse_list_.emplace_back(edge);
				}
			}
		}
	}

	auto find_root(Vertex vertex) -> Vertex
	{
		Vertex tmp = root_of_vertex_[vertex];
		if(tmp == vertex) { return vertex; }

		std::vector<Vertex> path;
		Vertex root{};
		do {
			root = tmp;
			path.emplace_back(root);
			tmp = root_of_vertex_[root];
		} while(tmp != root);

		// now root == (tmp = root_of_vertex_[root])

		for(const auto v : path) { root_of_vertex_[v] = root; }
		return root;
	}

	void merge_boundary(Vertex root1, Vertex root2)
	{
		border_vertices_[root1].insert(border_vertices_[root2].cbegin(),
									   border_vertices_[root2].cend());

		for(auto vertex : border_vertices_[root2])
		{
			if(connection_counts_[vertex] == lattice_.vertex_connection_count(vertex))
			{
				border_vertices_[root1].erase(vertex);
			}
		}
		border_vertices_.erase(root2);
	}

	void fusion()
	{
		while(!fuse_list_.empty())
		{
			Edge fuse_edge = fuse_list_.front();
			fuse_list_.pop_front();
			auto root1 = find_root(fuse_edge.u);
			auto root2 = find_root(fuse_edge.v);

			if(root1 == root2)
			{
				continue; // do nothing
			}

			peeling_edges_.push_back(fuse_edge);

			// let the size of the cluster of root1 be larger than that of root2
			if(mgr_.size(root1) < mgr_.size(root2)) { std::swap(root1, root2); }

			root_of_vertex_[root2] = root1;

			if(!mgr_.is_root(root2)) // if merging one is a single vertex
			{
				++mgr_.size(root1);
				border_vertices_[root1].emplace(root2);
			}
			else
			{
				mgr_.merge(root1, root2);
				merge_boundary(root1, root2);
			}
		}
	}

	auto peeling(std::vector<uint32_t>& syndromes) -> std::vector<Edge>
	{
		std::vector<Edge> corrections;

		// vs vector?
		tsl::robin_map<Vertex, int> vertex_count;

		for(Edge edge : peeling_edges_)
		{
			++vertex_count[edge.u];
			++vertex_count[edge.v];
		}

		while(!peeling_edges_.empty())
		{
			Edge leaf_edge = peeling_edges_.back();
			peeling_edges_.pop_back();
			auto u = Vertex{};
			auto v = Vertex{};
			if(vertex_count[leaf_edge.u] == 1)
			{
				u = leaf_edge.u;
				v = leaf_edge.v;
			}
			else if(vertex_count[leaf_edge.v] == 1)
			{
				u = leaf_edge.v;
				v = leaf_edge.u;
			}
			else // not a leaf
			{
				peeling_edges_.push_front(leaf_edge);
				continue;
			}

			--vertex_count[u];
			--vertex_count[v];

			if(syndromes[u] == 1)
			{
				corrections.emplace_back(leaf_edge);
				syndromes[u] = 0;
				syndromes[v] ^= 1U;
			}
		}
		return corrections;
	}

public:
	template<typename... Args> explicit Decoder(Args&&... args) : lattice_{args...} { }

	auto decode(std::vector<uint32_t>& syndromes) -> std::vector<Edge>
	{
		assert(syndromes.size() == lattice_.num_vertices());
		std::vector<Vertex> syndrome_vertices;
		for(uint32_t n = 0; n < syndromes.size(); ++n)
		{
			if((syndromes[n] % 2) != 0) { syndrome_vertices.emplace_back(n); }
		}

		init_cluster(syndrome_vertices);

		while(!mgr_.isempty_odd_root())
		{
			for(auto root : mgr_.odd_roots()) { grow(root); }
			fusion();
		}

		return peeling(syndromes);
	}

	[[nodiscard]] inline auto num_vertices() const -> int
	{
		return lattice_.num_vertices();
	}

	[[nodiscard]] inline auto num_edges() const -> int { return lattice_.num_edges(); }

	[[nodiscard]] inline auto edge_idx(const Edge& edge) const -> int
	{
		return lattice_.edge_idx(edge);
	}

	void clear()
	{
		std::deque<Edge>().swap(fuse_list_);

		mgr_.clear();

		tsl::robin_map<Vertex, tsl::robin_set<Vertex>>().swap(border_vertices_);

		std::deque<Edge>().swap(peeling_edges_);
	}
};
} // namespace UnionFindCPP
