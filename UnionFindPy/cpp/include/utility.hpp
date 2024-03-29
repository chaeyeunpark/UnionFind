#pragma once
#include <algorithm>
#include <climits>
#include <cstdint>
#include <functional>
#include <ostream>
#include <vector>

#include <nlohmann/json.hpp>

#ifndef NDEBUG
#define DEBUG
#endif

namespace UnionFindCPP
{

enum class ErrorType
{
	X,
	Z
};

struct Edge
{
	uint32_t u;
	uint32_t v;

	Edge(uint32_t ul, uint32_t vl)
	{
		u = std::min(ul, vl);
		v = std::max(ul, vl);
	}

	inline auto operator==(const Edge& rhs) const -> bool
	{
		return (u == rhs.u) && (v == rhs.v);
	}
};

void to_json(nlohmann::json& j, const Edge& e);
void from_json(const nlohmann::json& j, Edge& e);
auto operator<<(std::ostream& os, const UnionFindCPP::Edge& e) -> std::ostream&;
} // namespace UnionFindCPP

template<> struct std::hash<UnionFindCPP::Edge>
{
	auto operator()(const UnionFindCPP::Edge& e) const noexcept -> std::size_t
	{
		auto h1 = std::hash<uint32_t>()(e.u);
		auto h2 = std::hash<uint32_t>()(e.v);
		return h1 ^ (h2 << 1U);
	}
};
