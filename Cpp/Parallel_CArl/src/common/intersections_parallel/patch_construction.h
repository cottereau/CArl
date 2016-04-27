/*
 * patch_construction.h
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_
#define COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_

#include "common_header.h"
#include "common_header_libmesh.h"
#include "mesh_tables.h"

#include "CGAL_typedefs.h"

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
	const libMesh::Mesh&				   			m_Mesh;
	std::unique_ptr<libMesh::PointLocatorBase>      m_Patch_Point_Locator;
	Intersection_Tools								m_Intersection_Test;

	const libMesh::Parallel::Communicator& m_comm;

	std::unordered_set<unsigned int> 									m_Patch_Indexes;
	std::unordered_map<unsigned int,std::unordered_set<unsigned int> >	m_Patch_Neighbours;

	// Data structures used by the advancing front intersection search
	unsigned int 		m_working_element_id;
	bool				m_bTestNeighsForNewPairs;
	std::deque<int> 	m_element_intersection_queue;
	std::deque<int> 	m_element_test_queue;
	std::unordered_map<unsigned int,int>  	m_element_already_treated;
	std::unordered_map<unsigned int,int> 	m_element_inside_intersection_queue;
	std::unordered_set<unsigned int>		m_element_neighbours_to_search;

public:

	Patch_construction(const libMesh::Mesh & mesh) :
		m_Mesh { mesh },
		m_comm { m_Mesh.comm() }
	{
		m_Patch_Point_Locator = m_Mesh.sub_point_locator();

		// 		Instruction needed to avoid the code from crashing if a query
		//	is outside the mesh
		m_Patch_Point_Locator->enable_out_of_mesh_mode();

		m_working_element_id = 0;
		m_bTestNeighsForNewPairs = true;
	};

	const libMesh::Mesh & mesh()
	{
		return m_Mesh;
	}

	std::unordered_set<unsigned int> & indexes()
	{
		return m_Patch_Indexes;
	};

	unsigned int size()
	{
		return m_Patch_Indexes.size();
	}

	const libMesh::Elem * elem(unsigned int idx)
	{
		return m_Mesh.elem(idx);
	}

	std::unordered_map<unsigned int,std::unordered_set<unsigned int> > & patch_elem_neighbours()
	{
		return m_Patch_Neighbours;
	}

	/*
	 * 			Implementation of the patch construction algorithm without any
	 * 		neighboring information concerning the query element.
	 *
	 */
	void BuildPatch(const libMesh::Elem 	* Query_elem)
	{
		bool bDoIntersect = false;
		const libMesh::Elem		* First_Patch_elem = m_Intersection_Test.FindFirstIntersection(Query_elem,m_Patch_Point_Locator);
		m_Patch_Indexes.clear();
		m_Patch_Neighbours.clear();

		// Deque containing the indices of the elements to test
		std::deque<int> Patch_Test_Queue;

		// Unordered set, used to avoid double testing elements
		std::unordered_set<int> Treated_From_Mesh(m_Mesh.n_elem());

		// Index and pointer to element being tested right now
		unsigned int	Tested_idx;

		// Candidate index
		unsigned int 	Candidate_idx;

		// First element is ok!
		Treated_From_Mesh.insert(First_Patch_elem->id());
		m_Patch_Indexes.insert(First_Patch_elem->id());

		libMesh::Elem * elem_candidate;
		for(unsigned int iii = 0; iii < First_Patch_elem->n_neighbors(); ++iii)
		{
			elem_candidate = First_Patch_elem->neighbor(iii);
			if(elem_candidate != NULL)
			{
				Patch_Test_Queue.push_back(elem_candidate->id());
				Treated_From_Mesh.insert(elem_candidate->id());
			}
		}

		// Debug vars
		int nbOfTests = 1;
		int nbOfPositiveTests = 1;

		while(!Patch_Test_Queue.empty())
		{
			// Extract element from the list
			Tested_idx = Patch_Test_Queue[0];
			Patch_Test_Queue.pop_front();
			const libMesh::Elem 	* Tested_elem = m_Mesh.elem(Tested_idx);

			// Test it
			bDoIntersect = m_Intersection_Test.libMesh_exact_do_intersect(Query_elem,Tested_elem);
			++nbOfTests;

			// If it does intersect ...
			if(bDoIntersect)
			{
				++nbOfPositiveTests;

				// Add it to the output list ...
				m_Patch_Indexes.insert(Tested_idx);

				// ... And add its neighbours (if they weren't tested yet)
				for(unsigned int iii = 0; iii < Tested_elem->n_neighbors(); ++iii)
				{
					elem_candidate = Tested_elem->neighbor(iii);
					if(elem_candidate != NULL)
					{
						Candidate_idx = elem_candidate->id();
						if(Treated_From_Mesh.find(Candidate_idx)==Treated_From_Mesh.end())
						{
							Patch_Test_Queue.push_back(Candidate_idx);
							Treated_From_Mesh.insert(Candidate_idx);
						}
					}
				}
			}
		}

		// Now, let us build the patch's neighbour table
		std::unordered_set<unsigned int> dummy_neighbour_set(12);

		std::unordered_set<unsigned int>::iterator it_set = m_Patch_Indexes.begin();
		std::unordered_set<unsigned int>::iterator it_set_end = m_Patch_Indexes.end();
		m_Patch_Neighbours.reserve(m_Patch_Indexes.size());

		libMesh::Elem * elem_neighbour;
		unsigned int elem_neighbour_id;

		for( ; it_set != it_set_end; ++it_set)
		{
			const libMesh::Elem 	* elem = m_Mesh.elem(*it_set);
			dummy_neighbour_set.clear();

			for(unsigned int iii = 0; iii < elem->n_neighbors(); ++iii)
			{
				elem_neighbour = elem->neighbor(iii);

				if(elem_neighbour != NULL)
				{
					// Then test if it is inside the patch
					elem_neighbour_id = elem_neighbour->id();
					if(m_Patch_Indexes.find(elem_neighbour_id) != m_Patch_Indexes.end())
					{
						// Then the neighbor is inside the patch, insert it
						dummy_neighbour_set.insert(elem_neighbour->id());
					}
				}
			}
			m_Patch_Neighbours.insert(std::pair<unsigned int,std::unordered_set<unsigned int> >(*it_set,dummy_neighbour_set));
		}

		std::cout << "    DEBUG: patch search results" << std::endl;
		std::cout << " -> Nb. of intersections found : " << m_Patch_Indexes.size() << std::endl << std::endl;

		std::cout << " -> Nb. of mesh elements       : " << m_Mesh.n_elem() << std::endl;
		std::cout << " -> Patch size %               : " << 100.*m_Patch_Indexes.size()/m_Mesh.n_elem() << " %" << std::endl << std::endl;

		std::cout << " -> Nb. of tests               : " << nbOfTests << std::endl;
		std::cout << " -> Nb. of positive tests      : " << nbOfPositiveTests << std::endl;
		std::cout << " -> Positive %                 : " << 100.*nbOfPositiveTests/nbOfTests << " %" << std::endl << std::endl;
	}

	/*
	 *
	 * 		Methods and getters associated to the advancing front search method.
	 *
	 */

	/*
	 * 	Getters, setters and extractors
	 */
	void intersection_queue_push_back(unsigned int elem_id)
	{
		m_element_intersection_queue.push_back(elem_id);
	}

	void set_elem_as_treated(unsigned int elem_id)
	{
		m_element_already_treated[elem_id] = 1;
	}

	void set_elem_as_inside_queue(unsigned int elem_id)
	{
		m_element_inside_intersection_queue[elem_id] = 1;
	}

	bool intersection_queue_empty()
	{
		return m_element_intersection_queue.empty();
	}

	bool test_queue_empty()
	{
		return m_element_test_queue.empty();
	}

	unsigned int intersection_queue_extract_front_elem()
	{
		m_working_element_id = m_element_intersection_queue[0];
		m_element_intersection_queue.pop_front();
		return m_working_element_id;
	}

	unsigned int test_queue_extract_front_elem()
	{
		m_working_element_id = m_element_test_queue[0];
		m_element_test_queue.pop_front();
		return m_working_element_id;
	};

	unsigned int current_elem_id()
	{
		return m_working_element_id;
	};

	const libMesh::Elem * current_elem_pointer()
	{
		return m_Mesh.elem(m_working_element_id);
	}

	std::unordered_set<unsigned int> & neighbors_to_search_next_pair()
	{
		return m_element_neighbours_to_search;
	}

	// Determinate which neighbors need to be tested
	bool set_neighbors_to_search_next_pairs()
	{
		m_bTestNeighsForNewPairs = true;

		// Iterator over all the neighbors
		std::unordered_set<unsigned int>::iterator it_neigh, it_neigh_end;

		it_neigh 		= m_Patch_Neighbours[m_working_element_id].begin();
		it_neigh_end	= m_Patch_Neighbours[m_working_element_id].end();

		// Clear the "to test" set
		m_element_neighbours_to_search.clear();

		// And fill it, if needed
		for( ; it_neigh != it_neigh_end; ++ it_neigh)
		{
			if( m_element_inside_intersection_queue[*it_neigh] == 0 )
			{
				// This element does not have an initial intersection pair,
				// add it to the list
				m_element_neighbours_to_search.insert(*it_neigh);
			}
		}

		if(m_element_neighbours_to_search.empty())
		{
			// Then do not do the intersection test!
			m_bTestNeighsForNewPairs = false;
		}

		return m_bTestNeighsForNewPairs;
	}

	void add_neighbors_to_test_list()
	{
		std::unordered_set<unsigned int>::iterator it_neigh, it_neigh_end;

		it_neigh 		= m_Patch_Neighbours[m_working_element_id].begin();
		it_neigh_end	= m_Patch_Neighbours[m_working_element_id].end();

		for( ; it_neigh != it_neigh_end; ++ it_neigh)
		{
			if( m_element_already_treated[*it_neigh] == 0 )
			{
				// New element!
				m_element_test_queue.push_back(*it_neigh);
				m_element_already_treated[*it_neigh] = 1;
			}
		}
	};

	/*
	 * 	Initialize the deques, vectors and sets
	 */
	void FrontSearch_initialize()
	{
		homemade_assert_msg(!m_Patch_Indexes.empty(),"Patch is empty!");

		m_element_intersection_queue.clear();
		m_element_test_queue.clear();

		m_element_already_treated.clear();
		m_element_inside_intersection_queue.clear();

		m_element_already_treated.reserve(2*m_Patch_Indexes.size());
		m_element_inside_intersection_queue.reserve(2*m_Patch_Indexes.size());

		m_element_neighbours_to_search.clear();
		m_element_neighbours_to_search.reserve(12);

		std::unordered_set<unsigned int>::iterator it_idx, it_idx_end;
		it_idx 		= m_Patch_Indexes.begin();
		it_idx_end	= m_Patch_Indexes.end();

		for( ; it_idx != it_idx_end; ++it_idx)
		{
			m_element_already_treated[*it_idx] = 0;
			m_element_inside_intersection_queue[*it_idx] = 0;
		}

		m_working_element_id = 0;
	}

	/*
	 * 	Reset everything, with the exception of the intersection queue
	 */
	void FrontSearch_reset()
	{
		homemade_assert_msg(!m_Patch_Indexes.empty(),"Patch is empty!");

		m_element_test_queue.clear();

		std::unordered_set<unsigned int>::iterator it_idx, it_idx_end;
		it_idx 		= m_Patch_Indexes.begin();
		it_idx_end	= m_Patch_Indexes.end();

		for( ; it_idx != it_idx_end; ++it_idx)
		{
			m_element_already_treated[*it_idx] = 0;
			m_element_inside_intersection_queue[*it_idx] = 0;
		}

		m_working_element_id = 0;
	}

	/*
	 * 	Prepare the patch for a new round of tests
	 */
	unsigned int FrontSearch_prepare_for_probed_test()
	{
		// Extract the intersection queue's first element
		intersection_queue_extract_front_elem();

		// Add it to the test queue
		m_element_test_queue.clear();
		m_element_test_queue.push_back(m_working_element_id);

		return m_working_element_id;
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_ */
