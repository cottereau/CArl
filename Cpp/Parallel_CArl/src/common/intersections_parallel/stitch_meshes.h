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

// *******************
// Stitch_Meshes class
// *******************

/** \brief Class used to stitch together different meshes.
 *
 *		This class supposes that the meshes being stitched have compatible interfaces
 *	- in other words, that the stitched faces have the same vertices and surface elements.
 *  For the intersection search algorithms, it also joins the local intersection tables 
 *	into a single, global table.
 *
 * 		To avoid rounding error problems during the mesh constructions, each vertex point 
 *  is associated to a discrete point inside a grid. This grid is constructed taking into
 *  account the dimensions of the mesh, in a similar fashion to what is done for the 
 *  carl::Mesh_Intersection class. (**idea**: create a grid class, create a more primitive class ignoring the intersections).
 */

class	Stitch_Meshes
{
protected:
	// Communicators and parallel data.
	const libMesh::Parallel::Communicator& 	m_world_comm;	///< MPI Communicator.
	const unsigned int						m_nodes;		///< Number of processors.
	const unsigned int 						m_rank;			///< Processor rank.

	/// Final, stitched mesh.
	libMesh::ReplicatedMesh&				   	m_Stitched_mesh;

	/// File names of the meshes to be stitched.
	std::vector<std::string>				m_mesh_filenames;

	/// File names of intersection tables to joined.
	std::vector<std::string>				m_table_filenames;
	unsigned int 							m_nb_files;

	// Output files
	std::string								m_base_output;	//< Base output filename
	std::string								m_mesh_output;	//< Mesh output filename
	std::string								m_table_output;	//< Intersection table output filename

	/// Vector containing all the intersection pairs
	std::vector<std::pair<unsigned int, unsigned int> > m_intersection_pairs;

	// Unordered sets containing the elements from each mesh involved
	// in the intersection
	std::unordered_set<unsigned int> m_restriction_set_first;	//< Set of elements used for the restriction of the first mesh
	std::unordered_set<unsigned int> m_restriction_set_second;	//< Set of elements used for the restriction of the second mesh

	/// nNmber of elements inside each intersection
	std::vector<unsigned int> m_intersection_nb_of_elements;

	// Data structure and variables associated to the "grid" used to collapse
	// elements that are too small - reproduces the Mesh_Intersection class!

	/// Precision of the discrete point grid.
	double		m_eps;

	/// Grid minimum volume
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

	// Parameters of the final mesh
	unsigned int m_nb_of_intersections;	//< Final mesh's number of intersections.
	unsigned int m_nb_of_elements;		//< Final mesh's number of elements.
	unsigned int m_nb_of_nodes;			//< Final mesh's number of nodes.
	unsigned int m_maximum_nb_of_nodes;	//< Upper limit for the final mesh's number of nodes.

	// Boolean set when the filename base is also set
	bool m_bFilenamesSet; //< Have the filenames been set?

	// Boolean set when the grid is defined
	bool m_bGridDefined;		//< Is the grid defined?
	bool m_bGridPreallocated; 	//< Is the grid preallocated?

	// Boolean controlling if the stitched intersection mesh was finished
	bool m_bMeshFinalized;		//< Is final mesh finalized?

	// Perflog and debug variables
	bool m_bPrintDebug;		//< Print debug information? *Default:* false.

	// PROTECTED constructor
	Stitch_Meshes();

	// PROTECTED methods

public:

	// Constructors
	Stitch_Meshes(	libMesh::Mesh& output_mesh, const std::string output_filename = "test_stitched",  long grid_n_min = static_cast<long>(1E9), bool debugOutput = true) :
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

	/// Returns the stitched mesh
	const libMesh::ReplicatedMesh & mesh();

	/// Returns the pointer to the set with the elements used to form the first restricted mesh.
	const std::unordered_set<unsigned int>* get_restricted_set_pointer_first();

	/// Returns the pointer to the set with the elements used to form the second restricted mesh.
	const std::unordered_set<unsigned int>* get_restricted_set_pointer_second();

	// PUBLIC methods
	/// Set by hand the filenames used to build the stitched mesh.
	void set_base_filenames(std::vector<std::string> & mesh_filenames, std::vector<std::string> & table_filenames);

	/// Set the filenames used to build the stitched mesh from a common filename base.
	void set_base_filenames(const std::string & filename_base = std::string("test_r_"), const std::string & mesh_format = std::string(".msh"), unsigned int nb_of_files = 0);

	/// Preallocate the discrete points grid.
	void preallocate_grid(int map_preallocation);

	///Set the boundaries of the discrete points grid, using the intersected meshes as a base.
	void set_grid_constraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B, double vol_tol = -1);

	///Copy the boundaries of the discrete points grid from a carl::Mesh_Intersection object.
	void set_grid_constraints(Mesh_Intersection & mesh_inter_obj);

	/// Join the intersection tables.
	void join_tables();

	/// Stitch the meshes.
	void stitch_meshes();

	/// Convert a real valued point to a discrete point.
	void convert_to_discrete(const libMesh::Point& iPoint, std::vector<long>& oPoint);


};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_STITCH_MESHES_H_ */
