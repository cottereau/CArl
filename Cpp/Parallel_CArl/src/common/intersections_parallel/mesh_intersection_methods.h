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

#include "algorithm"

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
	// Communicators and parallel data
	/*
	 * 		The global variants are used for the full system
	 */
	const libMesh::Parallel::Communicator& 	m_comm;
	const unsigned int						m_nodes;
	const unsigned int 						m_rank;

	const libMesh::Parallel::Communicator& 	m_global_comm;
	const unsigned int						m_global_nodes;
	const unsigned int 						m_global_rank;

	// Address of the intersection mesh
	libMesh::SerialMesh&				   	m_libMesh_Mesh;

	// Mesh guarding the tetrahedrization of the intersection polyhedron
	libMesh::SerialMesh					   	m_libMesh_PolyhedronMesh;

	// TetGen interface
	libMesh::TetGenMeshInterface 			m_TetGenInterface;

	// CGAL Kernel converters
	Kernel_to_ExactKernel 		   ConvertInexactToExact;
	ExactKernel_to_Kernel 		   ConvertExactToInexact;

	// Data structures used to save the data needed to build the intersection meshes
	unsigned int 															 m_nb_of_intersections;
	std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > m_intersection_pairs;
	std::unordered_map<unsigned int, unsigned int>  						 m_intersection_couplings;
	std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > m_intersection_element_range;

	// Data structure and variables associated to the "grid" used to collapse
	// elements that are too small.

	// Precision
	double		m_eps;

	// Minimum volume
	double		m_vol_tol;

	// Number of integers over each dimension of the grid
	std::vector<long > m_GridN;
	std::vector<long > m_dummy_discrete_point;

	// Min and max points of the grid
	libMesh::Point m_Grid_MinPoint;
	libMesh::Point m_Grid_MaxPoint;

	// Map between the indexes and the grid positions
//	std::unordered_map<long, unsigned int>	m_Grid_to_mesh_vertex_idx;
	std::unordered_map<std::vector<long>, unsigned int, PointHash_3D, PointHash_3D_Equal >
		m_discrete_vertices;

	// Minimal of integers over each dimension of the grid
	long m_GridN_min;

	// Parameters of the final mesh
	unsigned int m_nb_of_elements;
	unsigned int m_nb_of_vertices;

	std::vector<unsigned int> 	m_intersection_point_indexes;
	unsigned int 				m_nb_of_points;

	//	Boolean controlling if the intersection mesh was finished or not
	bool m_bMeshFinalized;

	// Perflog and debug variables
	bool m_bPrintDebug;

	// PROTECTED constructor
	Mesh_Intersection();

	// PROTECTED methods
	/*
	 * 		Add new points to the intersection data structures
	 */
	// void update_intersection_vertices(	const std::set<libMesh::Point> & input_points);

	/*
	 * 		Triangulate an intersection defined by a set of libMesh::Points
	 */
	void triangulate_intersection(const std::set<libMesh::Point> & input_points);

	/*
	 * 		Update the final mesh with the polyhedron mesh data
	 */
	void update_intersection_mesh();

	/*
	 * 		Insert a new intersection pair
	 */
	void update_intersection_pairs(unsigned int elem_idx_A, unsigned int elem_idx_B, unsigned int inter_id);

	/*
	 * 		Update the intersection range of the intersection no. inter_id
	 */
	void update_intersection_element_range(unsigned int range_start, unsigned int range_end, unsigned int inter_id);

public:

	// Constructors
	Mesh_Intersection(	libMesh::SerialMesh & mesh, const libMesh::Mesh & mesh_A,
						const libMesh::Mesh & mesh_B,
						int map_preallocation = 1E6, long grid_n_min = static_cast<long>(1E9), bool debugOutput = false) :
		m_comm { mesh.comm() },
		m_nodes { m_comm.size() },
		m_rank { m_comm.rank() },
		m_global_comm { mesh_A.comm() },
		m_global_nodes { m_global_comm.size() },
		m_global_rank { m_global_comm.rank() },
		m_libMesh_Mesh { mesh },
		m_libMesh_PolyhedronMesh { libMesh::SerialMesh(m_comm) },
		m_TetGenInterface { libMesh::TetGenMeshInterface(m_libMesh_PolyhedronMesh) },
		m_nb_of_intersections { 0 },
		m_eps { -1 },
		m_GridN { std::vector<long >(3,-1) },
		m_dummy_discrete_point { std::vector<long >(3,-1) },
		m_Grid_MinPoint { libMesh::Point(0,0,0) },
		m_Grid_MaxPoint { libMesh::Point(1,1,1) },
		m_GridN_min { grid_n_min },
		m_nb_of_elements { 0 },
		m_nb_of_vertices { 0 },
		m_nb_of_points { 0 },
		m_bMeshFinalized { false },
		m_bPrintDebug { debugOutput }
	{
//		m_Grid_to_mesh_vertex_idx.reserve(map_preallocation);
		m_discrete_vertices.reserve(map_preallocation);
		m_intersection_point_indexes.resize(24);

		set_grid_constraints(mesh_A,mesh_B);

		m_libMesh_Mesh.reserve_nodes(map_preallocation);
		m_libMesh_Mesh.reserve_elem(10*map_preallocation);

		m_intersection_pairs.reserve(map_preallocation);
		m_intersection_couplings.reserve(map_preallocation);
		m_intersection_element_range.reserve(map_preallocation);

		m_libMesh_Mesh.allow_renumbering(false);
	};

	// Getters
	const libMesh::SerialMesh & mesh();
	libMesh::Point & min_point();
	libMesh::Point & max_point();
	double eps();
	double min_vol();
	std::vector<long> & grid_sizes();
	long grid_min_size();

	// PUBLIC methods
	/*
	 * 		Restart/initialize the data structures
	 */
	void initialize();

	/*
	 * 		Preallocate the discrete points grid
	 */
	void preallocate_grid( int map_preallocation );

	/*
	 * 		Set the boundaries of the discrete points grid
	 */
	void set_grid_constraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B, double vol_tol = -1);

	/*
	 *		Increase the intersection mesh and update its data structures
	 */
	void increase_intersection_mesh(	const std::set<libMesh::Point> & input_points,
										unsigned int elem_idx_A,
										unsigned int elem_idx_B,
										unsigned int elem_idx_C);

	/*
	 *		Export both the intersection mesh its data structures
	 */
	void export_intersection_data(const std::string & filename_base, const std::string & mesh_format = std::string(".e"));

	/*
	 * 		Prepare the mesh for use
	 */
	void prepare_for_use();

	/*
	 * 		Calculate the mesh's total volume
	 */
	double get_total_volume();

	/*
	 * 		Calculate the volume of an polyhedron defined by the point set
	 */
	double get_intersection_volume(std::set<libMesh::Point> & input_points);

	/*
	 * 		Convert a point to an grid index
	 */
	void convert_to_discrete(const libMesh::Point& iPoint, std::vector<long>& oPoint);
//	long convert_to_grid(const libMesh::Point iPoint);
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_ */
