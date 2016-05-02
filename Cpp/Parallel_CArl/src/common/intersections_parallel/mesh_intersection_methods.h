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

	//	Data structure and variables associated to the "grid" used to collapse
	//  elements that are too small.

	double		m_eps;

	std::vector<long int > m_GridN;

	libMesh::Point m_Grid_MinPoint;
	libMesh::Point m_Grid_MaxPoint;

	std::unordered_map<long int, unsigned int>	m_Grid_to_mesh_vertex_idx;

	int m_nb_of_elements;
	int m_nb_of_vertices;


public:

	Mesh_Intersection(libMesh::SerialMesh & mesh, double average_length = 1E-6, int map_preallocation = 1E6) :
		m_comm { mesh.comm() },
		m_libMesh_Mesh { mesh },
		m_libMesh_PolyhedronMesh { libMesh::SerialMesh(m_comm) },
		m_TetGenInterface { libMesh::TetGenMeshInterface(m_libMesh_PolyhedronMesh) },
		m_eps { average_length * 1E-6 },
		m_GridN { std::vector<long int >(3,-1) },
		m_Grid_MinPoint { libMesh::Point(0,0,0) },
		m_nb_of_elements { 0 },
		m_nb_of_vertices { 0 }
	{
		m_rank = m_comm.rank();
		m_Grid_to_mesh_vertex_idx.reserve(map_preallocation);
	};

	Mesh_Intersection(libMesh::SerialMesh & mesh, const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B, double average_length = 1E-6, int map_preallocation = 1E6) :
		m_comm { mesh.comm() },
		m_libMesh_Mesh { mesh },
		m_libMesh_PolyhedronMesh { libMesh::SerialMesh(m_comm) },
		m_TetGenInterface { libMesh::TetGenMeshInterface(m_libMesh_PolyhedronMesh) },
		m_eps { average_length * 1E-6 },
		m_GridN { std::vector<long int >(3,-1) },
		m_Grid_MinPoint { libMesh::Point(0,0,0) },
		m_nb_of_elements { 0 },
		m_nb_of_vertices { 0 }
	{
		m_rank = m_comm.rank();
		m_Grid_to_mesh_vertex_idx.reserve(map_preallocation);

		set_grid_contraints(mesh_A,mesh_B);
	};

	void preallocate_grid( int map_preallocation )
	{
		m_Grid_to_mesh_vertex_idx.reserve(map_preallocation);
	}

	void set_grid_contraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B)
	{
		libMesh::MeshTools::BoundingBox bbox_A = libMesh::MeshTools::bounding_box(mesh_A);
		libMesh::MeshTools::BoundingBox bbox_B = libMesh::MeshTools::bounding_box(mesh_B);

		// Just to be sure, test if the bboxes intersect!
		homemade_assert_msg(bbox_A.intersect(bbox_B),"Meshes' bounding boxes do not intersect!\n");

		// Set the (future) intersection bbox corners
		for(unsigned int iii = 0; iii < 3; ++iii)
		{
			m_Grid_MinPoint(iii) = std::max(bbox_A.min()(iii),bbox_B.min()(iii)) - 2*m_eps;
			m_Grid_MaxPoint(iii) = std::min(bbox_A.max()(iii),bbox_B.max()(iii)) + 2*m_eps;
			m_GridN[iii] = (m_Grid_MaxPoint(iii) - m_Grid_MinPoint(iii)) / m_eps + 1;
		}
	}

	const libMesh::SerialMesh & mesh()
	{
		return m_libMesh_Mesh;
	}

	void triangulate_intersection(const std::set<libMesh::Point> & input_points)
	{
		m_libMesh_PolyhedronMesh.clear();

		for(std::set<libMesh::Point>::const_iterator it_set = input_points.begin();
				it_set != input_points.end();
				++it_set)
		{
			m_libMesh_PolyhedronMesh.add_point(*it_set);
		}
		m_TetGenInterface.triangulate_pointset();
	}

	void increase_intersection_mesh(	const std::set<libMesh::Point> & input_points,
										std::vector<unsigned int> & points_indexes)
	{
		// 	First, add the points to the grid-to-mesh map, and create new
		// vertices if needed.
		update_intersection_vertices(input_points,points_indexes);

		//	Second, triagulate the point set
		triangulate_intersection(input_points);

		// 	The sets conserve the order, so there is no need of an equivalence
		// table between the m_libMesh_PolyhedronMesh and m_libMesh_Mesh
		// vertices.




	}

	void update_intersection_vertices(	const std::set<libMesh::Point> & input_points,
										std::vector<unsigned int> & points_indexes)
	{
		long int dummy_long_int;
		points_indexes.clear();
		for(	std::set<libMesh::Point>::const_iterator it_points = input_points.begin();
				it_points != input_points.end();
				++it_points)
		{
			dummy_long_int = convert_to_grid(*it_points);
			if(m_Grid_to_mesh_vertex_idx.find(dummy_long_int)!=m_Grid_to_mesh_vertex_idx.end())
			{
				// New vertex! Add it to the mesh
				m_Grid_to_mesh_vertex_idx[dummy_long_int] = m_nb_of_vertices;
				points_indexes.push_back(m_nb_of_vertices);
				m_libMesh_Mesh.add_point(*it_points);
				++m_nb_of_vertices;
			}
			else
			{
				// Recover the index corresponding to the point
				points_indexes.push_back(m_Grid_to_mesh_vertex_idx[dummy_long_int]);
			}
		}
	}

	double get_intersection_volume(std::set<libMesh::Point> & input_points)
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

	long int convert_to_grid(const libMesh::Point iPoint)
	{
		return    lround( (iPoint(0) -  m_Grid_MinPoint(0))/m_eps) * m_GridN[1]*m_GridN[2]
				+ lround( (iPoint(1) -  m_Grid_MinPoint(1))/m_eps) * m_GridN[1]
				+ lround( (iPoint(2) -  m_Grid_MinPoint(2))/m_eps);
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_ */
