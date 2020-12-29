#pragma once
#include <cstdint>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <queue>
#include <algorithm>
#include <nlohmann/json.hpp>

#include <iostream>

#include <tsl/robin_set.h>
#include <tsl/robin_map.h>

#include "utility.hpp"
#include "RootManager.hpp"


class UnionFindDecoder
{
public:
	using Vertex = int;
	using RootIterator = tsl::robin_set<Vertex>::const_iterator;

private:
	int L_;

	/* index: vertex */
	std::vector<Vertex> connection_counts_;
	std::vector<int> support_;
	std::deque<Edge> fuse_list_;

	/* index: vertex */
	std::vector<Vertex> root_of_vertex_; // root of vertex

	RootManager mgr_;
	/* key: root, value: borders */
	tsl::robin_map<Vertex, tsl::robin_set<Vertex>> border_vertices_; 

	/* Data for peeling */
	std::deque<Edge> peeling_edges_;


	inline int to_vertex_index(int row, int col) const
	{
		return ::to_vertex_index(L_, row, col);
	}

	void init_cluster(const std::vector<int>& roots)
	{
		connection_counts_ = std::vector<Vertex>(L_*L_, 0);
		support_ = std::vector<int>(2*L_*L_, 0);
		mgr_.initialize_roots(roots);
		for(auto root: roots)
		{
			border_vertices_[root].emplace(root);
		}

		root_of_vertex_.resize(L_*L_);

		for(int u = 0; u < L_*L_; ++u)
		{
			root_of_vertex_[u] = u;
		}
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


	void grow(Vertex root)
	{
		for(auto border_vertex: border_vertices_[root])
		{
			for(auto v: vertex_connections(border_vertex))
			{
				auto edge = Edge(border_vertex, v);

				int& elt = support_[decoder_edge_to_qubit_idx(L_, edge, ErrorType::Z)];
				++elt;
				
				if(elt == 2)
				{
					connection_counts_[edge.u] ++;
					connection_counts_[edge.v] ++;
					fuse_list_.emplace_back(edge);
				}
			}
		}
	}

	Vertex find_root(Vertex vertex)
	{
		Vertex tmp = root_of_vertex_[vertex];
		if(tmp == vertex)
			return vertex;

		std::vector<Vertex> path;
		Vertex root;
		do
		{
			root = tmp;
			path.emplace_back(root);
			tmp = root_of_vertex_[root];
		}while(tmp != root);
		
		// now root == (tmp = root_of_vertex_[root])

		for(const auto v: path)
		{
			root_of_vertex_[v] = root;
		}
		return root;
	}

	void merge_boundary(Vertex root1, Vertex root2)
	{
		border_vertices_[root1].insert(
			border_vertices_[root2].cbegin(), border_vertices_[root2].cend());
		
		for(auto vertex: border_vertices_[root2])
		{
			if(connection_counts_[vertex] == 4)
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
				continue; //do nothing
			}
			
			peeling_edges_.push_back(fuse_edge);


			// let the size of the cluster of root1 be larger than that of root2
			if(mgr_.size(root1) < mgr_.size(root2))
				std::swap(root1, root2);

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

	std::vector<Edge> peeling(std::vector<int>& syndromes)
	{
		std::vector<Edge> corrections;
		tsl::robin_map<Vertex, int> vertex_count;

		for(Edge edge: peeling_edges_)
		{
			++vertex_count[edge.u];
			++vertex_count[edge.v];
		}

		while(!peeling_edges_.empty())
		{
			Edge leaf_edge = peeling_edges_.back();
			peeling_edges_.pop_back();
			Vertex u, v;
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
				--syndromes[u];
				syndromes[v] = 1-syndromes[v];
			}
		}
		return corrections;
	}


public:
	UnionFindDecoder(int L)
		: L_{L}
	{
	}

	std::vector<Edge> decode(std::vector<int>& syndromes)
	{
		std::vector<Vertex> syndrome_vertices = 
			syndrome_locations(L_, syndromes);

		init_cluster(syndrome_vertices);

		while(!mgr_.isempty_odd_root())
		{
			for(auto root: mgr_.odd_roots())
			{
				grow(root);
			}
			fusion();
		}

		return peeling(syndromes);
	}

	void clear()
	{
		std::deque<Edge>().swap(fuse_list_);

		mgr_.clear();

		tsl::robin_map<Vertex, tsl::robin_set<Vertex> >().swap(
				border_vertices_);

		std::deque<Edge>().swap(peeling_edges_);
	}

	nlohmann::json clusters()
	{
		std::map<Vertex, std::vector<Vertex>> cluster;
		for(Vertex v = 0; v < L_*L_; ++v)
		{
			Vertex root = find_root(v);
			if(mgr_.size(root) != 0)
			{
				cluster[root].emplace_back(v);
			}
		}
		return cluster;
	}
};
