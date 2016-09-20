/*
 * common_functions.h
 *
 *  Created on: Jan 26, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_COMMON_FUNCTIONS_H_
#define COMMON_COMMON_FUNCTIONS_H_

#include "common_header_libmesh.h"
#include "common_header.h"

libMesh::Real kronecker_delta(unsigned int i,
				   unsigned int j);

void clear_line();

/*
 *  A little structure used to save the intersection data
 */

namespace carl
{

void invert_index_unordered_map(
		const std::unordered_map<int,int>& input_map,
		std::unordered_map<int,int>& output_map);

template<typename T>
void jump_lines(T& filestream, unsigned int numberOfLines = 1)
{
	std::string dummy;
	for(unsigned int iii = 0; iii < numberOfLines; ++iii)
		std::getline(filestream,dummy);
};

struct IntersectionData
{
	int InterMeshIdx;
	int AMeshIdx;
	int BMeshIdx;
	int IntersectionID;
};

// Point hash function
struct PointHash_3D {
	std::size_t operator()(const std::vector<long>& k) const
	{
		long prime0 = 73856093;
		long prime1 = 19349669;
		long prime2 = 83492791;
		long primeN = 2038074743;

		return ( ( k[0] * prime0 ) ^ ( k[1] * prime1 ) ^ ( k[2] * prime2 ) ) % primeN;
	}
};
 
struct PointHash_3D_Equal {
	bool operator()(const std::vector<long>& lhs, const std::vector<long>& rhs) const
	{
		return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
	}
};

}

#endif /* COMMON_COMMON_FUNCTIONS_H_ */
