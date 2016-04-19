/*
 * intersection_mesh_libmesh.h
 *
 *  Created on: Apr 11, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_
#define COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_

#include "common_header.h"
#include "common_header_libmesh.h"
#include "mesh_tables.h"

#include "CGAL_typedefs.h"

namespace carl
{

/*
 * 		Mesh_Intersection class
 *
 * 			This class is an interface for the construction of the intersection
 * 		meshes - hence why the mesh member is not constant.
 *
 */

class	Mesh_Intersection
{
protected:
	libMesh::Mesh&				   m_libMesh_Mesh;
	const libMesh::Parallel::Communicator& m_comm;

public:

	Mesh_Intersection(libMesh::Mesh & mesh) :
		m_libMesh_Mesh { mesh },
		m_comm { m_libMesh_Mesh.comm() }

	{
	};

	const libMesh::Mesh & mesh()
	{
		return m_libMesh_Mesh;
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_ */
