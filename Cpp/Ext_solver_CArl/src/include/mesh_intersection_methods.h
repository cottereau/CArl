/*
 * intersection_mesh_libmesh.h
 *
 *  Created on: Apr 11, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_
#define COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_

#include "carl_headers.h"
#include "mesh_tables.h"

#include "algorithm"

namespace carl
{

// ***********************
// Mesh_Intersection class
// ***********************

/** \brief Class used to construct the intersection meshes.
 *
 *		It contains methods to tetrahedrize polyhedrons (using either CGAL or Tetgen) 
 *  and to increase the intersection mesh. To avoid rounding error problems during the
 *	mesh constructions, each vertex point is associated to a discrete point inside a 
 *	grid. This grid is constructed taking into account the dimensions of the meshes 
 *  for which we want to find the intersections. 
 */

class	Mesh_Intersection
{
protected:
	// Communicators and parallel data
	/*
	 * 		The global variants are used for the full system
	 */
	const libMesh::Parallel::Communicator& 	m_comm;		///< Intersection mesh communicator.
	const unsigned int						m_nodes;	///< Intersection mesh number of processors (**useless!**).
	const unsigned int 						m_rank;		///< Intersection mesh processor rank (**useless!**).

	const libMesh::Parallel::Communicator& 	m_global_comm;	///< Full system communicator.
	const unsigned int						m_global_nodes;	///< Full system number of processors.
	const unsigned int 						m_global_rank;	///< Full system processor rank.

	libMesh::ReplicatedMesh&				   	m_libMesh_Mesh;	///< Address of the intersection mesh.

	/// libMesh Mesh containing the tetrahedrization of an intersection polyhedron (the polyhedron mesh).
	libMesh::Mesh					   	    m_libMesh_PolyhedronMesh;

	/// Auxiliary CGAL mesh, used if CGAL is chosen for the tetrahedrization.
	DT_3		 							m_CGAL_PolyhedronMesh;

	/// TetGen interface.
	libMesh::TetGenMeshInterface 			m_TetGenInterface;

	// CGAL Kernel converters
	Kernel_to_ExactKernel 			ConvertInexactToExact;	///< Convert inexact CGAL constructs to exact constructs.
	ExactKernel_to_Kernel 			ConvertExactToInexact;	///< Convert exact CGAL constructs to inexact constructs.

	unsigned int 															 m_nb_of_intersections;

	/// Unordered map containing the element pair associated to each intersection (key).
	std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > m_intersection_pairs;

	/// Unordered map associating each intersection (key) to a coupling element (**useless!**)
	std::unordered_map<unsigned int, unsigned int>  						 m_intersection_couplings;

	/// Unordered map associating each intersection (key) to a range of elements in the intersection mesh.
 	std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > m_intersection_element_range;

	// Data structure and variables associated to the "grid" used to collapse
	// elements that are too small.

	/// Precision of the discrete point grid.
	double		m_eps;

	/// Minimal volume for the intersection polyhedrons.
	double		m_vol_tol;

	/// Dimensions of the discrete point grid.
	std::vector<long > m_GridN;

	std::vector<long > m_dummy_discrete_point;

	/// Minimal point of the discrete grid.
	libMesh::Point m_Grid_MinPoint;
	
	/// Maximum point of the discrete grid.
	libMesh::Point m_Grid_MaxPoint;

	/// Unordered map (or hash) for the discrete points (key).
	std::unordered_map<std::vector<long>, unsigned int, PointHash_3D, PointHash_3D_Equal >
		m_discrete_vertices;

	/// Minimal number of integers over each dimension of the discrete grid.
	long m_GridN_min;

	/// Number of elements of the intersection mesh.
	unsigned int m_nb_of_elements;

	/// Number of vertices of the intersection mesh.
	unsigned int m_nb_of_vertices;

	/// (**Unused!**)
	std::vector<unsigned int> 	m_intersection_point_indexes;

	/// (**Unused!**)
	unsigned int 				m_nb_of_points;

	///	Boolean controlling if the intersection mesh was finished or not.
	bool m_bMeshFinalized;

	/// Print debug information? *Default:* false.
	bool m_bPrintDebug;

	/// Choice of meshing algorithm.
	/** 	Can be either IntersectionMeshingMethod::LIBMESH_TETGEN or IntersectionMeshingMethod:: CGAL. 
	 *  *Default:* IntersectionMeshingMethod:: CGAL.
	 */
	IntersectionMeshingMethod m_MeshingMethod;

	/// (Protected) default constructor
 	Mesh_Intersection();

	// PROTECTED methods
	/// Update the intersection mesh with the current polyhedron mesh data.
	void update_intersection_mesh();

	/// Update the intersection pair map with a new pair associated to the intersection no. inter_id.
	void update_intersection_pairs(unsigned int elem_idx_A, unsigned int elem_idx_B, unsigned int inter_id);

	/// Update the intersection element range map with a new range associated to the intersection no. inter_id
	void update_intersection_element_range(unsigned int range_start, unsigned int range_end, unsigned int inter_id);

public:

	// Constructors
	Mesh_Intersection(	libMesh::ReplicatedMesh & mesh, const libMesh::Mesh & mesh_A,
						const libMesh::Mesh & mesh_B, IntersectionMeshingMethod MeshingMethod = IntersectionMeshingMethod::CGAL,
						int map_preallocation = 1E6, long grid_n_min = static_cast<long>(1E9), bool debugOutput = false) :
		m_comm { mesh.comm() },
		m_nodes { m_comm.size() },
		m_rank { m_comm.rank() },
		m_global_comm { mesh_A.comm() },
		m_global_nodes { m_global_comm.size() },
		m_global_rank { m_global_comm.rank() },
		m_libMesh_Mesh { mesh },
		m_libMesh_PolyhedronMesh { libMesh::Mesh(m_comm) },
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
		m_bPrintDebug { debugOutput },
		m_MeshingMethod { MeshingMethod }
	{
		m_discrete_vertices.reserve(map_preallocation);
		m_intersection_point_indexes.resize(24);

		set_grid_constraints(mesh_A,mesh_B);

		m_libMesh_Mesh.reserve_nodes(map_preallocation);
		m_libMesh_Mesh.reserve_elem(10*map_preallocation);

		m_intersection_pairs.reserve(map_preallocation);
		m_intersection_couplings.reserve(map_preallocation);
		m_intersection_element_range.reserve(map_preallocation);

		m_libMesh_Mesh.allow_renumbering(false);

		switch(m_MeshingMethod)
		{
		case IntersectionMeshingMethod::LIBMESH_TETGEN :
			std::cout << " -> Using TetGen to generate mesh" << std::endl;
			break;
		case IntersectionMeshingMethod::CGAL :
			std::cout << " -> Using CGAL to generate mesh" << std::endl;
			break;
		}
	};

	// Getters
	const libMesh::ReplicatedMesh & mesh();		///< Returns the address of the intersection mesh.
	libMesh::Point & min_point();				///< Returns the minimal point of the discrete grid.
	libMesh::Point & max_point();				///< Returns the maximum point of the discrete grid.
	double eps();						///< Returns the discrete grid step.
	double min_vol();					///< Returns the minimal polyhedron volume.
	std::vector<long> & grid_sizes();	///< Returns the discrete grid sizes.
	long grid_min_size();				///< Returns the minimum discrete grid size.

	// PUBLIC methods
	/// Initialize the data structures.
	void initialize();

	/// Preallocate the discrete points grid.
	void preallocate_grid( int map_preallocation );

	/// Set the boundaries of the discrete points grid, taking into account the geometry of meshes mesh_A and mesh_B.
	void set_grid_constraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B, double vol_tol = -1);

	/// Triangulate an intersection defined by a set of libMesh::Points, using either CGAL or Tetgen.
	void triangulate_intersection(const std::set<libMesh::Point> & input_points);

	/// Increase the intersection mesh with the current polyhedron mesh and update its data structures.
	void increase_intersection_mesh(	const std::set<libMesh::Point> & input_points,
										unsigned int elem_idx_A,
										unsigned int elem_idx_B,
										unsigned int elem_idx_C);

	/// Export both the intersection mesh its data structures.
	void export_intersection_data(const std::string & filename_base, const std::string & mesh_format = std::string(".e"));

	/// Prepare the intersection mesh for use
	void prepare_for_use();

	/// Calculate the intersection mesh's total volume
	double get_total_volume();

	/// Calculate the volume of an polyhedron defined by the point set.
	double get_intersection_volume(std::set<libMesh::Point> & input_points);

	/// Convert a real valued point to a discrete point.
	void convert_to_discrete(const libMesh::Point& iPoint, std::vector<long>& oPoint);
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_MESH_LIBMESH_H_ */
