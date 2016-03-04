/*
 * mpi_carl_tools.h
 *
 *  Created on: Feb 22, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_MPI_CARL_TOOLS_H_
#define COMMON_LIBMESH_CODE_MPI_CARL_TOOLS_H_

#include "common_header.h"
#include "common_header_libmesh.h"

namespace carl
{

void broadcast_index_unordered_map(
		std::unordered_map<int,int>& index_map,
		const libMesh::Parallel::Communicator& CommComm,
		int origin_rank = 0);

}
#endif /* COMMON_LIBMESH_CODE_MPI_CARL_TOOLS_H_ */
