/*
 * restrict_mesh.cpp
 *
 *  Created on: May 17, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "restrict_mesh.h"
namespace carl
{

void Mesh_restriction::BuildRestrictionFromSet(const std::unordered_set<unsigned int> * restricted_mesh_set)
{
	// Intersection and distribution work done already by the stitch algorithm.

	// Only do the work over a single processor, to avoid communications.
	if(m_rank == 0)
	{
		m_Patch_Elem_indexes.clear();
		m_Patch_Node_indexes.clear();
		m_Patch_Elem_Neighbours.clear();

		std::unordered_set<unsigned int >::iterator elem_idx_it = restricted_mesh_set->begin();
		std::unordered_set<unsigned int >::iterator elem_idx_it_end = restricted_mesh_set->end();

		for( ; elem_idx_it != elem_idx_it_end; ++ elem_idx_it)
		{
			const libMesh::Elem * elem_to_add = m_Mesh.elem(*elem_idx_it);
			insert_patch_element(elem_to_add);
		}

		build_patch_mesh();
	}
}

void Mesh_restriction::BuildRestriction(const libMesh::ReplicatedMesh 	& Coupling_mesh)
{
	bool bDoIntersect = false;

	// Deque containing the indices of the elements to test
	std::deque<int> Restriction_Test_Queue;

	// Unordered set, used to avoid double testing elements
	unsigned int treated_from_mesh_preallocation = m_Mesh.n_elem();
	std::unordered_set<int> Treated_From_Mesh;

	// Only do the work over a single processor, to avoid communications.
	if(m_rank ==0)
	{
		Treated_From_Mesh.reserve(treated_from_mesh_preallocation);

		std::set<unsigned int> Intersecting_elems;

		// Index and pointer to element being tested right now
		unsigned int	Tested_idx;

		// Candidate index
		unsigned int 	Candidate_idx;

		m_Patch_Elem_indexes.clear();
		m_Patch_Node_indexes.clear();
		m_Patch_Elem_Neighbours.clear();

		std::set<unsigned int>::iterator it_set_start;

		libMesh::Mesh::const_element_iterator it_coupling = Coupling_mesh.elements_begin();
		libMesh::Mesh::const_element_iterator it_coupling_end = Coupling_mesh.elements_end();

		// Debug vars
		int nbOfTests = 1;
		int nbOfPositiveTests = 1;

		for( ; it_coupling != it_coupling_end; ++it_coupling)
		{
			const libMesh::Elem * Query_elem = * it_coupling;

			Treated_From_Mesh.clear();

			// Find elements from the original mesh intersecting the coupling element
			m_Intersection_Test.FindAllIntersection(Query_elem,m_Patch_Point_Locator,Intersecting_elems);

			libMesh::Elem * 	First_Restriction_elems = NULL;
			libMesh::Elem * 	elem_candidate = NULL;
			it_set_start      = Intersecting_elems.begin();

			for( ; it_set_start != Intersecting_elems.end(); ++it_set_start)
			{
				Treated_From_Mesh.insert(*it_set_start);
				First_Restriction_elems = m_Mesh.elem(*it_set_start);
				insert_patch_element(First_Restriction_elems);

				for(unsigned int iii = 0; iii < First_Restriction_elems->n_neighbors(); ++iii)
				{
					elem_candidate = First_Restriction_elems->neighbor(iii);
					if(elem_candidate != NULL && Treated_From_Mesh.find(elem_candidate->id())==Treated_From_Mesh.end())
					{
						Restriction_Test_Queue.push_back(elem_candidate->id());
						Treated_From_Mesh.insert(elem_candidate->id());
					}
				}
			}

			while(!Restriction_Test_Queue.empty())
			{
				// Extract element from the list
				Tested_idx = Restriction_Test_Queue[0];
				Restriction_Test_Queue.pop_front();
				const libMesh::Elem 	* Tested_elem = m_Mesh.elem(Tested_idx);

				// Test it
				bDoIntersect = m_Intersection_Test.libMesh_exact_do_intersect(Query_elem,Tested_elem);
				++nbOfTests;

				// If it does intersect ...
				if(bDoIntersect)
				{
					++nbOfPositiveTests;

					// Add it to the output list ...
					insert_patch_element(Tested_elem);

					// ... And add its neighbours (if they weren't tested yet)
					for(unsigned int iii = 0; iii < Tested_elem->n_neighbors(); ++iii)
					{
						elem_candidate = Tested_elem->neighbor(iii);
						if(elem_candidate != NULL)
						{
							Candidate_idx = elem_candidate->id();
							if(Treated_From_Mesh.find(Candidate_idx)==Treated_From_Mesh.end())
							{
								Restriction_Test_Queue.push_back(Candidate_idx);
								Treated_From_Mesh.insert(Candidate_idx);
							}
						}
					}
				}
			}
		}

		build_patch_mesh();

		if(m_bPrintDebug)
		{
			std::cout << "    DEBUG: Restriction search results" << std::endl;
			std::cout << " -> Nb. of intersections found  : " << m_Patch_Elem_indexes.size() << std::endl << std::endl;

			std::cout << " -> Nb. of mesh elements        : " << m_Mesh.n_elem() << std::endl;
			std::cout << " -> Nb. of Restriction elements : " << m_Patch_Elem_indexes.size() << std::endl;
			std::cout << " -> Patch elem %                : " << 100.*m_Patch_Elem_indexes.size()/m_Mesh.n_elem() << " %" << std::endl << std::endl;

			std::cout << " -> Nb. of mesh nodes           : " << m_Mesh.n_nodes() << std::endl;
			std::cout << " -> Nb. of patch nodes          : " << m_Patch_Node_indexes.size() << std::endl;
			std::cout << " -> Patch node %                : " << 100.*m_Patch_Node_indexes.size()/m_Mesh.n_nodes() << " %" << std::endl << std::endl;

			std::cout << " -> Nb. of tests                : " << nbOfTests << std::endl;
			std::cout << " -> Nb. of positive tests       : " << nbOfPositiveTests << std::endl;
			std::cout << " -> Positive %                  : " << 100.*nbOfPositiveTests/nbOfTests << " %" << std::endl << std::endl;
		}
	}
}

void Mesh_restriction::export_restriction_mesh(const std::string & filename_base)
{
	if(m_rank == 0)
	{
		std::string filename_mesh = filename_base + ".msh";

		// Print mesh
		libMesh::GmshIO output_mesh(m_Mesh_patch);
		output_mesh.binary() = true;
		output_mesh.write(filename_mesh);

		std::string filename_elements = filename_base + "_restrict.dat";

		std::ofstream elems_out(filename_elements);

		elems_out << m_elem_map_Output_Global.size() << std::endl;
		for(unsigned int iii = 0; iii < m_elem_map_Output_Global.size(); ++iii)
		{
			elems_out << iii << " " << m_elem_map_Output_Global[iii] << std::endl;
		}

		elems_out.close();
	}
}

libMesh::ReplicatedMesh &  Mesh_restriction::restricted_mesh()
{
	return m_Mesh_patch;
}
}
