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

	std::unordered_set<unsigned int> 							m_Patch_Indexes;
	std::unordered_map<unsigned int,std::set<unsigned int> >	m_Patch_Neighbours;

public:

	Patch_construction(const libMesh::Mesh & mesh) :
		m_Mesh { mesh },
		m_comm { m_Mesh.comm() }
	{
		m_Patch_Point_Locator = m_Mesh.sub_point_locator();

		// 		Instruction needed to avoid the code from crashing if a query
		//	is outside the mesh
		m_Patch_Point_Locator->enable_out_of_mesh_mode();
	};

	const libMesh::Mesh & mesh()
	{
		return m_Mesh;
	}

	/*
	 * 			Implementation of the patch construction algorithm without any
	 * 		neighboring information concerning the query element.
	 *
	 */
	void BuildPatch(const libMesh::Elem 	* Query_elem, std::unordered_set<unsigned int> & Patch_Indexes)
	{
		bool bDoIntersect = false;
		const libMesh::Elem		* First_Patch_elem = m_Intersection_Test.FindFirstIntersection(Query_elem,m_Patch_Point_Locator);
		m_Patch_Indexes.clear();

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
		std::set<unsigned int> dummy_neighbour_set;

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
			m_Patch_Neighbours.insert(std::pair<unsigned int,std::set<unsigned int> >(*it_set,dummy_neighbour_set));
		}

		std::cout << "    DEBUG: patch search results" << std::endl;
		std::cout << " -> Nb. of intersections found : " << m_Patch_Indexes.size() << std::endl << std::endl;

		std::cout << " -> Nb. of mesh elements       : " << m_Mesh.n_elem() << std::endl;
		std::cout << " -> Patch size %               : " << 100.*m_Patch_Indexes.size()/m_Mesh.n_elem() << " %" << std::endl << std::endl;

		std::cout << " -> Nb. of tests               : " << nbOfTests << std::endl;
		std::cout << " -> Nb. of positive tests      : " << nbOfPositiveTests << std::endl;
		std::cout << " -> Positive %                 : " << 100.*nbOfPositiveTests/nbOfTests << " %" << std::endl << std::endl;

		Patch_Indexes = m_Patch_Indexes;
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_ */
