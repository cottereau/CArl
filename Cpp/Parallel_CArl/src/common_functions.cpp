#include "common_functions.h"

libMesh::Real kronecker_delta(unsigned int i,
				   unsigned int j)
{
	return i == j ? 1. : 0.;
};

void carl::invert_index_unordered_map(
		const std::unordered_map<int,int>& input_map,
		std::unordered_map<int,int>& output_map)
{
	int map_length = input_map.size();
	output_map.reserve(map_length);

	std::unordered_map<int,int>::const_iterator mapIt = input_map.begin();
	std::unordered_map<int,int>::const_iterator end_mapIt = input_map.end();

	for( ; mapIt != end_mapIt; ++mapIt)
	{
		output_map[mapIt->second] = mapIt->first;
	}
}

void clear_line()
{
	std::cout << '\r' << "                                                                               " << "\r";
};
