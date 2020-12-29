#include <nlohmann/json.hpp>
#include <fstream>
#include "UnionFind.hpp"

int main(int argc, char* argv[])
{
	using nlohmann::json;
	std::ifstream ifs(argv[1]);
	json jf = json::parse(ifs);
	UnionFindDecoder decoder(jf["L"]);

	std::vector<int> syndromes; 
	jf.at("syndromes_for_z").get_to<std::vector<int>>(syndromes);
	auto decoding = decoder.decode(syndromes);

	json to_print;
	to_print["cluster"] = decoder.clusters();
	to_print["decoding"] = decoding;
	std::ofstream fout("out.json");
	fout << to_print;

	return 0;
}
