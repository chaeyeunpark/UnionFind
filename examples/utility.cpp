#include "utility.hpp"

void to_json(nlohmann::json& j, const Edge& e)
{
	j = nlohmann::json{e.u, e.v};
}
void from_json(const nlohmann::json& j, Edge& e)
{
	e = Edge(j[0], j[1]);
}

