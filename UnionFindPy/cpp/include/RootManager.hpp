#pragma once
#include <vector>

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

#include <nlohmann/json.hpp>

class RootManager
{
public:
	using Vertex = uint32_t;

private:
	/* set of roots */
	tsl::robin_set<Vertex> roots_;
	/* set of roots with odd parity */
	tsl::robin_set<Vertex> odd_roots_; // root size comb
	/* key: root, value: size of the cluster */
	tsl::robin_map<Vertex, uint32_t> size_;
	/* key: root, value: parity of the cluster */
	tsl::robin_map<Vertex, uint32_t> parity_;

public:
	class SizeProxy
	{
	private:
		RootManager& mgr_;
		Vertex root_;

	public:
		SizeProxy(RootManager& mgr, Vertex root) : mgr_{mgr}, root_{root} { }

		// NOLINTNEXTLINE(google-explicit-constructor,hicpp-explicit-conversions)
		operator uint32_t() const
		{
			if(!mgr_.is_root(root_)) { return 0; }
			const auto it = mgr_.size_.find(root_);
			return it->second;
		}

		auto operator=(int new_size) -> SizeProxy&
		{
			mgr_.size_[root_] = new_size;
			return *this;
		}

		auto operator++() -> SizeProxy&
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

	inline auto size(Vertex root) -> SizeProxy { return SizeProxy(*this, root); }

	[[nodiscard]] inline auto size(Vertex root) const -> uint32_t
	{
		if(!is_root(root)) { return 0; }
		const auto it = size_.find(root);
		return it->second;
	}

	[[nodiscard]] inline auto parity(Vertex root) const -> uint32_t
	{
		auto it = parity_.find(root);
		if(it == parity_.end()) { return 0; }
		return it->second;
	}

	[[nodiscard]] inline auto is_root(Vertex v) const -> bool
	{
		return roots_.count(v) == 1;
	}

	[[nodiscard]] inline auto is_odd_root(Vertex v) const -> bool
	{
		return odd_roots_.count(v) == 1;
	}

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

	[[nodiscard]] auto isempty_odd_root() const -> bool { return odd_roots_.empty(); }

	void clear()
	{
		tsl::robin_set<Vertex>().swap(roots_);
		tsl::robin_set<Vertex>().swap(odd_roots_);
		tsl::robin_map<Vertex, uint32_t>().swap(size_);
		tsl::robin_map<Vertex, uint32_t>().swap(parity_);
	}

	[[nodiscard]] auto odd_roots() const& -> const tsl::robin_set<Vertex>&
	{
		return odd_roots_;
	}
	[[nodiscard]] auto odd_roots() && -> tsl::robin_set<Vertex> { return odd_roots_; }

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
