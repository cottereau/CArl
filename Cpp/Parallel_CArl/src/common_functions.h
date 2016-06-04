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
	for(int iii = 0; iii < numberOfLines; ++iii)
		std::getline(filestream,dummy);
};

struct IntersectionData
{
	int InterMeshIdx;
	int AMeshIdx;
	int BMeshIdx;
	int IntersectionID;
};
}

#endif /* COMMON_COMMON_FUNCTIONS_H_ */
