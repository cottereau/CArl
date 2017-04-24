#include "common_functions.h"

int kronecker_delta(unsigned int i,
				   unsigned int j)
{
	return i == j ? 1. : 0.;
};

std::string carl::ClusterSchedulerType_to_string(ClusterSchedulerType input)
{
	switch (input)
	{
		case ClusterSchedulerType::PBS :	return "PBS";
						break;

		case ClusterSchedulerType::SLURM :	return "SLURM";
						break;
	}

	homemade_error_msg("Invalid enumerate argument!");
	return "";
};

std::string carl::BaseCGPrecondType_to_string(BaseCGPrecondType input)
{
	switch (input)
	{
		case BaseCGPrecondType::NO_PRECONDITIONER :	return "NONE";
						break;

		case BaseCGPrecondType::COUPLING_OPERATOR :	return "Coupling_operator";
						break;

		case BaseCGPrecondType::COUPLING_JACOBI :	return "Coupling_operator_jacobi";
						break;
	}

	homemade_error_msg("Invalid enumerate argument!");
	return "";
};

std::string carl::ExtSolverType_to_string(ExtSolverType input)
{
	switch (input)
	{
		case ExtSolverType::LIBMESH_LINEAR :	return "LIBMESH_LINEAR";
						break;

		case ExtSolverType::DUMMY :	return "DUMMY";
						break;
	}

	homemade_error_msg("Invalid enumerate argument!");
	return "";
};

std::string carl::exec_command(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != NULL)
            result += buffer.data();
    }
    return result;
}

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
