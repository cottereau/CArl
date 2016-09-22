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

int carl::voigt_index_converter(int aaa, int bbb)
{
	if(aaa == bbb)
	{
		// 00 -> 0, 11 -> 1, 22 -> 2
		return aaa;
	}
	else
	{
		// 12, 21 -> 3
		// 02, 20 -> 4
		// 01, 10 -> 5
		int sum = aaa + bbb;

		if(sum == 3)
		{
			return 3;
		}
		else if(sum == 2)
		{
			return 4;
		}
		else if(sum == 1)
		{
			return 5;
		}
	}

	std::cerr << "Bad indexes! " << aaa << " " << bbb << std::endl;
	homemade_error_msg(" You shouldn't be here! (voigt_index_converter)");
	return -1;
};
