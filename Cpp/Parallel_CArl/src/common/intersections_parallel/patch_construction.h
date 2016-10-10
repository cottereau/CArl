/*
 * patch_construction.h
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_
#define COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_

#include "carl_headers.h"
#include "mesh_tables.h"

#include "intersection_tools.h"

namespace carl
{

/*
 * 		Patch_Construction class
 *
 * 			This class contains the methods needed to build the patches used by
 * 		both the global intersection search algorithm and Gander's algorithm.
 *
 */

class	Patch_construction
{
protected:

	// Communicators and parallel data
	/*
	 * 		The local variants are used for the patch meshes
	 */
	const libMesh::Parallel::Communicator&  m_comm;
	const unsigned int 						m_nodes;
	const unsigned int 						m_rank;

	const libMesh::Parallel::Communicator&  m_local_comm;

	// Meshes and point locators
	libMesh::SerialMesh&						   	m_Mesh;
	libMesh::SerialMesh								m_Mesh_patch;
	std::unique_ptr<libMesh::PointLocatorBase>      m_Patch_Point_Locator;

	// Object used to do the intersection tests
	Intersection_Tools								m_Intersection_Test;

	// Data structures used to build the patch
	std::unordered_set<unsigned int> 									m_Patch_Elem_indexes;
	std::unordered_set<unsigned int>									m_Patch_Node_indexes;
	std::unordered_map<unsigned int,std::unordered_set<unsigned int> >	m_Patch_Elem_Neighbours;

	// Data structures used by the advancing front intersection search
	unsigned int 							m_working_element_id;
	bool									m_bTestNeighsForNewPairs;
	std::deque<int> 						m_element_intersection_queue;
	std::deque<int> 						m_element_test_queue;

	std::unordered_map<unsigned int,int>  	m_element_already_treated;
	std::unordered_map<unsigned int,int> 	m_element_inside_intersection_queue;
	std::unordered_set<unsigned int>		m_element_neighbours_to_search;

	std::unordered_map<unsigned int,int>				m_node_map_Global_Output;
	std::unordered_map<unsigned int,int>				m_node_map_Output_Global;
	std::unordered_map<unsigned int,int>				m_elem_map_Global_Output;
	std::unordered_map<unsigned int,int>				m_elem_map_Output_Global;

	// Perflog and debug variables
	bool m_bPrintDebug;

	// PROTECTED constructor
	Patch_construction();

	// PROTECTED methods
	/*
	 * 		Insert element inside thepatch data structures
	 */
	void insert_patch_element(const libMesh::Elem		* Patch_elem);

public:

	// Constructors
	Patch_construction(libMesh::Mesh & mesh, const libMesh::Parallel::Communicator& local_comm, bool debugOutput = false) :
		m_comm { mesh.comm() },
		m_nodes { m_comm.size() },
		m_rank { m_comm.rank() },
		m_local_comm {  local_comm },
		m_Mesh { mesh },
		m_Mesh_patch { libMesh::SerialMesh(m_local_comm) },

		m_bPrintDebug { debugOutput }
	{
		m_Patch_Point_Locator = m_Mesh.sub_point_locator();

		// 		Instruction needed to avoid the code from crashing if a query
		//	is outside the mesh
		m_Patch_Point_Locator->enable_out_of_mesh_mode();

		m_working_element_id = 0;
		m_bTestNeighsForNewPairs = true;

		m_Patch_Elem_indexes.reserve(m_Mesh.n_elem());
		m_Patch_Node_indexes.reserve(m_Mesh.n_nodes());
	};

	// Getters
	libMesh::SerialMesh & mesh();
	libMesh::SerialMesh & patch_mesh();
	std::unordered_set<unsigned int> & elem_indexes();
	std::unordered_set<unsigned int> & node_indexes();
	unsigned int size();
	const libMesh::Elem * elem(unsigned int idx);
	std::unordered_map<unsigned int,std::unordered_set<unsigned int> > & patch_elem_neighbours();

	// PUBLIC methods
	/*
	 * 		Implementation of the patch construction algorithm without any
	 * 	neighboring information concerning the query element.
	 */
	void BuildPatch(const libMesh::Elem 	* Query_elem);

	/*
	 * 		Index converters between the global mesh and the patch, and
	 * 	vice-versa.
	 */
	unsigned int convert_global_to_patch_elem_id(unsigned int input);
	unsigned int convert_patch_to_global_elem_id(unsigned int input);

	/*
	 * 		Build a patch mesh from the patch data structures
	 */
	void build_patch_mesh();

	/*
	 * 		Export the patch mesh to an file
	 */
	void export_patch_mesh(std::string & filename_base);

	/*
	 * 		Current element information
	 */
	unsigned int current_elem_id();
	const libMesh::Elem * current_elem_pointer();

	/*
	 *		The following methods are all used by the advancing front algorithm,
	 *	found inside the "intersection_search.h" files
	 */
	void intersection_queue_push_back(unsigned int elem_id);
	void set_elem_as_treated(unsigned int elem_id);
	void set_elem_as_inside_queue(unsigned int elem_id);
	bool intersection_queue_empty();
	bool test_queue_empty();
	unsigned int intersection_queue_extract_front_elem();
	unsigned int test_queue_extract_front_elem();
	std::unordered_set<unsigned int> & neighbors_to_search_next_pair();
	bool set_neighbors_to_search_next_pairs();
	void add_neighbors_to_test_list();

	void FrontSearch_initialize();
	void FrontSearch_reset();
	unsigned int FrontSearch_prepare_for_probed_test();
};

}

#endif /* COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_ */
