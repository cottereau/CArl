/*
 * mpi_carl_tools.h
 *
 *  Created on: Feb 22, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_MPI_CARL_TOOLS_H_
#define COMMON_LIBMESH_CODE_MPI_CARL_TOOLS_H_

#include "carl_headers.h"

namespace carl
{

void broadcast_index_unordered_map(
		std::unordered_map<int,int>& index_map,
		const libMesh::Parallel::Communicator& CommComm,
		int origin_rank = 0);

template <typename T>
inline void MPI_reduce_vector(std::vector<T> & r, int root, const libMesh::Parallel::Communicator& Comm)
{
	if (Comm.size() > 1 && !r.empty())
		{
			libmesh_assert(Comm.verify(r.size()));

			std::vector<T> temp(r);
			MPI_Reduce(&temp[0], &r[0], libMesh::cast_int<int>(r.size()),
						libMesh::Parallel::StandardType<T>(&temp[0]), MPI_SUM, root,
						Comm.get());
		}
	}
}
#endif /* COMMON_LIBMESH_CODE_MPI_CARL_TOOLS_H_ */
