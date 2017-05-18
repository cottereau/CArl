#include "common_functions.h"

int kronecker_delta(unsigned int i,
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

void carl::print_stats_to_file(std::vector<double>& vec_data, const std::string filename)
{
	std::ofstream output_stream(filename,std::ofstream::app);
	libMesh::StatisticsVector<double> statistics_vec(vec_data.size(),0);
	for(unsigned int iii = 0; iii < vec_data.size(); ++iii)
	{
		statistics_vec[iii] = vec_data[iii];
	}

	output_stream 	<< statistics_vec.minimum() << " "
					<< statistics_vec.maximum() << " "
					<< statistics_vec.mean() << " "
					<< statistics_vec.median() << " "
					<< statistics_vec.stddev() << std::endl;

	output_stream.close();
};
