/*
 * stitch_meshes.h
 *
 *  Created on: Jun 2, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_STITCH_MESHES_H_
#define COMMON_INTERSECTIONS_PARALLEL_STITCH_MESHES_H_

#include "carl_headers.h"
#include "mesh_tables.h"
#include "mesh_intersection_methods.h"

namespace carl
{
class	Stitch_Intersection_Meshes
{
protected:
	// Communicators and parallel data
	const libMesh::Parallel::Communicator& 	m_world_comm;
	const unsigned int						m_nodes;
	const unsigned int 						m_rank;

	// Address of the intersection mesh
	libMesh::ReplicatedMesh&				   	m_Stitched_mesh;

	// Vectors containing the input files
	std::vector<std::string>				m_mesh_filenames;
	std::vector<std::string>				m_table_filenames;
	unsigned int 							m_nb_files;

	// Output file
	std::string								m_base_output;
	std::string								m_mesh_output;
	std::string								m_table_output;

	// Vector containing all the intersection pairs ...
	std::vector<std::pair<unsigned int, unsigned int> > m_intersection_pairs;

	// Unordered sets containing the elements from each mesh involved
	// in the intersection
	std::unordered_set<unsigned int> m_restriction_set_first;
	std::unordered_set<unsigned int> m_restriction_set_second;

	// ... and the number of elements inside each intersection
	std::vector<unsigned int> m_intersection_nb_of_elements;

	// Data structure and variables associated to the "grid" used to collapse
	// elements that are too small - reproduces the Mesh_Intersection class!

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
	// std::unordered_map<long, unsigned int>	m_Grid_to_mesh_vertex_idx;
	std::unordered_map<std::vector<long>, unsigned int, PointHash_3D, PointHash_3D_Equal >
		m_discrete_vertices;

	// Minimal of integers over each dimension of the grid
	long m_GridN_min;

	// Parameters of the final mesh
	unsigned int m_nb_of_intersections;
	unsigned int m_nb_of_elements;
	unsigned int m_nb_of_nodes;
	unsigned int m_maximum_nb_of_nodes;

	// Boolean set when the filename base is also set
	bool m_bFilenamesSet;

	// Boolean set when the grid is defined
	bool m_bGridDefined;
	bool m_bGridPreallocated;

	// Boolean controlling if the stitched intersection mesh was finished
	bool m_bMeshFinalized;

	// Perflog and debug variables
	bool m_bPrintDebug;

	// PROTECTED constructor
	Stitch_Intersection_Meshes();

	// PROTECTED methods

public:

	// Constructors
	Stitch_Intersection_Meshes(	libMesh::Mesh& output_mesh, const std::string output_filename = "test_stitched",  long grid_n_min = static_cast<long>(1E9), bool debugOutput = true) :
		m_world_comm { output_mesh.comm() },
		m_nodes { m_world_comm.size() },
		m_rank { m_world_comm.rank() },
		m_Stitched_mesh { output_mesh },
		m_nb_files { 0 },
		m_base_output { output_filename },
		m_eps { -1 },
		m_vol_tol { -1 },
		m_GridN { std::vector<long> (3,-1) },
		m_dummy_discrete_point { std::vector<long> (3,-1) },
		m_Grid_MinPoint { libMesh::Point(0,0,0) },
		m_Grid_MaxPoint { libMesh::Point(1,1,1) },
		m_GridN_min { grid_n_min },
		m_nb_of_intersections { 0 },
		m_nb_of_elements { 0 },
		m_nb_of_nodes { 0 },
		m_maximum_nb_of_nodes { 0 },
		m_bFilenamesSet { false },
		m_bGridDefined { false },
		m_bGridPreallocated { false },
		m_bMeshFinalized { false },
		m_bPrintDebug { debugOutput }
	{
		m_Stitched_mesh.allow_renumbering(false);
	};

	// Getters
	const libMesh::ReplicatedMesh & mesh();

	// PUBLIC methods
	/*
	 * 		Set filenames - and prepare the grid preallocations!
	 */
	void set_base_filenames(std::vector<std::string> & mesh_filenames, std::vector<std::string> & table_filenames);
	void set_base_filenames(const std::string & filename_base = std::string("test_r_"), const std::string & mesh_format = std::string(".msh"), unsigned int nb_of_files = 0);

	/*
	 * 		Preallocate the discrete points grid
	 */
	void preallocate_grid(int map_preallocation);

	/*
	 * 		Set the boundaries of the discrete points grid - either by taking as
	 * 	an input the original meshes, either by copying the data from an
	 * 	Mesh_Intersection object.
	 */
	void set_grid_constraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B, double vol_tol = -1);
	void set_grid_constraints(Mesh_Intersection & mesh_inter_obj);

	/*
	 * 		Stitch the meshes!
	 */
	void join_tables();
	void stitch_meshes();

	/*
	 * 		Convert a point to an grid index
	 */
	void convert_to_discrete(const libMesh::Point& iPoint, std::vector<long>& oPoint);

	/*
	 * 		Get the intersecting element sets
	 */
	const std::unordered_set<unsigned int>* get_restricted_set_pointer_first();
	const std::unordered_set<unsigned int>* get_restricted_set_pointer_second();
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_STITCH_MESHES_H_ */
