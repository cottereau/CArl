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
	const libMesh::Parallel::Communicator& 	m_comm;
	int 									m_rank;

	libMesh::SerialMesh&				   	m_libMesh_Mesh;
	libMesh::SerialMesh					   	m_libMesh_PolyhedronMesh;

	libMesh::TetGenMeshInterface 			m_TetGenInterface;

	Kernel_to_ExactKernel 		   ConvertInexactToExact;
	ExactKernel_to_Kernel 		   ConvertExactToInexact;

public:

	Mesh_Intersection(libMesh::SerialMesh & mesh) :
		m_comm { mesh.comm() },
		m_libMesh_Mesh { mesh },
		m_libMesh_PolyhedronMesh { libMesh::SerialMesh(m_comm) },
		m_TetGenInterface { libMesh::TetGenMeshInterface(m_libMesh_PolyhedronMesh) }
	{
		m_rank = m_comm.rank();
	};

	const libMesh::SerialMesh & mesh()
	{
		return m_libMesh_Mesh;
	}

	void triangulate_intersection(std::set<Point_3> & input_points)
	{
		m_libMesh_PolyhedronMesh.clear();
		libMesh::Point	dummy_point;

		for(std::set<Point_3>::const_iterator it_set = input_points.begin();
				it_set != input_points.end();
				++it_set)
		{
			dummy_point = libMesh::Point(it_set->x(),it_set->y(),it_set->z());
			m_libMesh_PolyhedronMesh.add_point(dummy_point);
		}
		m_TetGenInterface.triangulate_pointset();
	}

	double get_intersection_volume(std::set<Point_3> & input_points)
	{
		double volume = 0;

		triangulate_intersection(input_points);

		for(libMesh::SerialMesh::element_iterator it_elem = m_libMesh_PolyhedronMesh.elements_begin();
				it_elem != m_libMesh_PolyhedronMesh.elements_end();
				++it_elem)
		{
			const libMesh::Elem * elem = * it_elem;
			volume += elem->volume();
		}

		return volume;
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_ */
