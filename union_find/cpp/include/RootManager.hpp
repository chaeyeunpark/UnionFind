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
#pragma once
#include <vector>

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

#include <nlohmann/json.hpp>

class RootManager
{
public:
	using Vertex = int;

private:
	/* set of roots */
	tsl::robin_set<Vertex> roots_;
	/* set of roots with odd parity */
	tsl::robin_set<Vertex> odd_roots_; // root size comb
	/* key: root, value: size of the cluster */
	tsl::robin_map<Vertex, int> size_;
	/* key: root, value: parity of the cluster */
	tsl::robin_map<Vertex, int> parity_;

public:
	class SizeProxy
	{
	private:
		RootManager& mgr_;
		Vertex root_;

	public:
		SizeProxy(RootManager& mgr, Vertex root) : mgr_{mgr}, root_{root} { }

		operator int() const
		{
			if(!mgr_.is_root(root_)) return 0;
			const auto it = mgr_.size_.find(root_);
			return it->second;
		}

		SizeProxy& operator=(int new_size)
		{
			mgr_.size_[root_] = new_size;
			return *this;
		}

		SizeProxy& operator++()
		{
			auto it = mgr_.size_.find(root_);
			++it.value();
			return *this;
		}
	};

	friend class SizeProxy;

	void initialize_roots(const std::vector<Vertex>& roots)
	{
		const auto n_reserve = 2 * roots.size();
		roots_.reserve(n_reserve);
		odd_roots_.reserve(n_reserve);
		size_.reserve(n_reserve);
		parity_.reserve(n_reserve);
		for(const auto root : roots)
		{
			roots_.emplace(root);
			odd_roots_.emplace(root);
			size_.emplace(root, 1);
			parity_.emplace(root, 1);
		}
	}

	inline SizeProxy size(Vertex root) { return SizeProxy(*this, root); }

	inline int size(Vertex root) const
	{
		if(!is_root(root)) return 0;
		const auto it = size_.find(root);
		return it->second;
	}

	inline int parity(Vertex root) const
	{
		auto it = parity_.find(root);
		if(it == parity_.end()) return 0;
		return it->second;
	}

	inline bool is_root(Vertex v) const { return roots_.count(v) == 1; }

	inline bool is_odd_root(Vertex v) const { return odd_roots_.count(v) == 1; }

	// size of the cluster of root1 is larger than that of root2
	void merge(Vertex root1, Vertex root2)
	{
		const auto new_parity = parity(root1) + parity(root2);

		if((new_parity % 2) == 1) { odd_roots_.emplace(root1); }
		else
		{
			odd_roots_.erase(root1);
		}

		size_[root1] += size(root2);
		parity_[root1] = new_parity;

		odd_roots_.erase(root2);

		size_.erase(root2);
		parity_.erase(root2);
		roots_.erase(root2);
	}

	void remove(Vertex root)
	{
		odd_roots_.erase(root);
		roots_.erase(root);
		size_.erase(root);
		parity_.erase(root);
	}

	bool isempty_odd_root() const { return odd_roots_.empty(); }

	void clear()
	{
		tsl::robin_set<Vertex>().swap(roots_);
		tsl::robin_set<Vertex>().swap(odd_roots_);
		tsl::robin_map<Vertex, int>().swap(size_);
		tsl::robin_map<Vertex, int>().swap(parity_);
	}

	const tsl::robin_set<Vertex>& odd_roots() const& { return odd_roots_; }
	tsl::robin_set<Vertex> odd_roots() && { return odd_roots_; }

	void print(std::ostream& os) const
	{
		nlohmann::json p;

		auto roots_j = nlohmann::json::array();
		for(auto root : roots_)
		{
			roots_j.emplace_back(
				nlohmann::json::array({root, size_.at(root), parity_.at(root)}));
		}
		p["roots"] = roots_j;
		p["odd_roots"] = odd_roots_;

		os << p << std::endl;
	}
};
