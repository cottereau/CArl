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

	DelaunayTri_3			   	m_temp_triangulation;
	Tetrahedron                	m_dummyTetrahedron;

	Kernel_to_ExactKernel 		   ConvertInexactToExact;
	ExactKernel_to_Kernel 		   ConvertExactToInexact;

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

	void triangulate_polyhedron(Polyhedron & input_poly)
	{
		m_temp_triangulation.clear();
		m_temp_triangulation.insert(input_poly.points_begin(),input_poly.points_end());
	}

	double get_polyhedron_volume(Polyhedron & input_poly)
	{
		double volume = 0;
		triangulate_polyhedron(input_poly);
//		if(!input_poly.is_empty())
//			std::cout << " Empty polyhedron ..." << std::endl;

		for(DelaunayTri_3::Finite_cells_iterator it_cell = m_temp_triangulation.finite_cells_begin();
			it_cell != m_temp_triangulation.finite_cells_end();
			++it_cell)
		{
			m_dummyTetrahedron = Tetrahedron(	it_cell->vertex(0)->point(),
					it_cell->vertex(1)->point(),
					it_cell->vertex(2)->point(),
					it_cell->vertex(3)->point());
			volume += std::abs(m_dummyTetrahedron.volume());
		}

		return volume;
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_ */
