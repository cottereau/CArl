#include "mpi_carl_tools.h"

void carl::broadcast_index_unordered_map(
		std::unordered_map<int,int>& index_map,
		const libMesh::Parallel::Communicator& CommComm,
		int origin_rank)
{
	int rank = CommComm.rank();
	int dummy_vector_size = -1;
	std::vector<int> dummy_vector;

	if(rank == origin_rank)
	{
		dummy_vector_size = index_map.size();
		dummy_vector.resize(2*dummy_vector_size);

		std::unordered_map<int,int>::iterator mapIt = index_map.begin();
		std::unordered_map<int,int>::iterator end_mapIt = index_map.end();
		int dummy_iii = 0;
		for( ; mapIt != end_mapIt; ++mapIt)
		{
			dummy_vector[2*dummy_iii] = mapIt->first;
			dummy_vector[2*dummy_iii + 1] = mapIt->second;
			++dummy_iii;
		}
	}

	CommComm.barrier();
	CommComm.broadcast(dummy_vector_size,origin_rank);

	if( rank != origin_rank )
	{
		dummy_vector.resize(2*dummy_vector_size);
	}
	CommComm.broadcast(dummy_vector,origin_rank);

	if( rank != origin_rank)
	{
		index_map.reserve(dummy_vector_size);
		for(int iii = 0; iii < dummy_vector_size; ++iii)
		{
			index_map[dummy_vector[2*iii]] = dummy_vector[2*iii+1];
		}
	}
};
