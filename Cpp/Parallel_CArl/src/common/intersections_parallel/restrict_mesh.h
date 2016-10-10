/*
 * restrict_mesh.h
 *
 *  Created on: May 17, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_RESTRICT_MESH_H_
#define COMMON_INTERSECTIONS_PARALLEL_RESTRICT_MESH_H_

#include "carl_headers.h"
#include "mesh_tables.h"

#include "intersection_tools.h"
#include "patch_construction.h"

namespace carl
{

/*
 * 		Mesh_restriction class
 *
 * 			This class contains the methods needed to restrict a mesh to the
 * 		region defined by the coupling mesh. This class has a lot in common
 * 		with the Patch_Construction class: essentially, the mesh restriction
 * 		process is identical to the patch building looped over all the coupling
 * 		mesh elements. Hence, it is a derived class.
 *
 */
class Mesh_restriction: public Patch_construction
{
public:
	Mesh_restriction(libMesh::Mesh & mesh, const libMesh::Parallel::Communicator& local_comm, bool debugOutput = false)
	: Patch_construction(mesh, local_comm, debugOutput)
	{

	};

	void BuildRestriction(const libMesh::SerialMesh 	& Coupling_mesh);

	libMesh::SerialMesh & restricted_mesh();

	void export_restriction_mesh(const std::string & filename_base);
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_RESTRICT_MESH_H_ */
