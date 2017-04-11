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


// **********************
// Mesh_restriction class
// **********************

/** \brief Class used to build a restriction of a parent mesh to the coupling region.
 *  
 *		This class is derived from the carl::Patch_construction class, and it contains 
 *	the methods needed to restrict a mesh to the region defined by the coupling mesh.
 */

class Mesh_restriction: public Patch_construction
{
public:

	// Constructors
	/// Constructor with a pre-defined parent mesh and a local communicator.
	Mesh_restriction(libMesh::Mesh & mesh, const libMesh::Parallel::Communicator& local_comm, bool debugOutput = false)
	: Patch_construction(mesh, local_comm, debugOutput)
	{

	};

	/// Build the restriction of the parent mesh to the coupling region defined by Coupling_mesh.
	void BuildRestriction(const libMesh::ReplicatedMesh 	& Coupling_mesh);

	/// Returns the restricted mesh.
	libMesh::ReplicatedMesh & restricted_mesh();

	///	Export the restricted mesh to a file.
	void export_restriction_mesh(const std::string & filename_base);

	/// Build the restriction of the parent mesh from a given element set. This version is useful if the intersection search was already done.
	void BuildRestrictionFromSet(const std::unordered_set<unsigned int> * restricted_mesh_set);
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_RESTRICT_MESH_H_ */
