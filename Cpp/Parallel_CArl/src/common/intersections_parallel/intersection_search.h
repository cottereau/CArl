/*
 * intersection_search.h
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_
#define COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_

#include "common_header.h"
#include "common_header_libmesh.h"
#include "mesh_tables.h"

#include "CGAL_typedefs.h"

#include "mesh_intersection_methods.h"
#include "patch_construction.h"

namespace carl
{

enum SearchMethod
{
	BRUTE,
	FRONT,
	BOTH
};

/*
 * 		Intersection_Search class
 *
 * 			This class contains the structure needed to find all the
 * 		intersections between two meshes A and B, and the coupling region mesh
 * 		(C or Coupling).
 *
 */

class	Intersection_Search
{

protected:
	libMesh::Mesh&				   m_Mesh_A;
	libMesh::Mesh&				   m_Mesh_B;
	libMesh::Mesh&				   m_Mesh_Coupling;

	const libMesh::Parallel::Communicator& m_comm;

	Patch_construction					   m_Patch_Constructor_A;
	Patch_construction					   m_Patch_Constructor_B;

	Mesh_Intersection					   m_Mesh_Intersection;

	std::unordered_map<unsigned int,std::pair<unsigned int,unsigned int> > m_Intersection_Pairs;

	std::unordered_map<unsigned int,std::set<unsigned int> > m_DEBUG_Intersection_Pairs;

	Intersection_Tools m_Intersection_test;
	Intersection_Tools m_Intersection_test_neighbors;

	bool MASTER_bPerfLog_intersection_search;

	libMesh::PerfLog m_perf_log;

public:

	Intersection_Search(libMesh::Mesh & mesh_A,
						libMesh::Mesh & mesh_B,
						libMesh::Mesh & mesh_Coupling,
						libMesh::Mesh & mesh_I,
						bool  bDoPerf_log = true) :
		m_Mesh_A { mesh_A },
		m_Mesh_B { mesh_B },
		m_Mesh_Coupling { mesh_Coupling },
		m_comm { m_Mesh_Coupling.comm() },
		m_Patch_Constructor_A { Patch_construction(m_Mesh_A)},
		m_Patch_Constructor_B { Patch_construction(m_Mesh_B)},
		m_Mesh_Intersection { Mesh_Intersection(mesh_I)},
		MASTER_bPerfLog_intersection_search {bDoPerf_log},
		m_perf_log { libMesh::PerfLog("Intersection search", MASTER_bPerfLog_intersection_search) }

	{
		// Reserve space for the unordered sets
		m_Intersection_Pairs.reserve(mesh_A.n_elem()*mesh_B.n_elem());
		m_DEBUG_Intersection_Pairs.reserve(mesh_A.n_elem()*mesh_B.n_elem());
	};

	libMesh::Mesh & mesh_A()
	{
		return m_Mesh_A;
	}

	libMesh::Mesh & mesh_B()
	{
		return m_Mesh_B;
	}

	libMesh::Mesh & mesh_Coupling()
	{
		return m_Mesh_Coupling;
	}

	void BuildCoupledPatches(const libMesh::Elem 	* Query_elem, int patch_counter)
	{
//		libMesh::SerialMesh mesh_out(m_Mesh_A.comm(),3);
//		mesh_out.prepare_for_use();

		m_Patch_Constructor_A.BuildPatch(Query_elem);
		std::string filename = "patch_mesh_A_" + std::to_string(patch_counter);
		m_Patch_Constructor_A.export_patch_mesh(filename);

		m_Patch_Constructor_B.BuildPatch(Query_elem);
		filename = "patch_mesh_B_" + std::to_string(patch_counter);
		m_Patch_Constructor_B.export_patch_mesh(filename);

//		std::unordered_set<unsigned int>::iterator it_idx_A, it_idx_A_end;
//		it_idx_A = m_Patch_Constructor_A.indexes().begin();
//		it_idx_A_end = m_Patch_Constructor_A.indexes().end();
//		for( ; it_idx_A != it_idx_A_end; ++it_idx_A)
//		{
//			std::cout << 749 + *it_idx_A << std::endl;
//		}
	}

	void BuildPatchIntersections_Brute(const libMesh::Elem 	* Query_elem)
	{
		// TODO: For now it is more of a "find" than a build ...

		m_perf_log.push("Preamble","Brute force algorithm");
		// Code for the brute force intersection tests
		m_Intersection_Pairs.clear();

		// Set containing all the intersection points
		std::set<Point_3> intersection_vertices;

		// Intersection_Tools
		m_Intersection_test.libmesh_set_coupling_nef_polyhedron(Query_elem);

		// Boolean: do they intersect?
		bool bDoIntersect = false;

		std::unordered_set<unsigned int> & Patch_Set_A = m_Patch_Constructor_A.elem_indexes();
		std::unordered_set<unsigned int> & Patch_Set_B = m_Patch_Constructor_B.elem_indexes();

		std::unordered_set<unsigned int>::iterator it_patch_A;
		std::unordered_set<unsigned int>::iterator it_patch_B;

		// Debug vars
		int nbOfTests = 0;
		int nbOfPositiveTests = 0;
		m_perf_log.pop("Preamble","Brute force algorithm");

		bool bCreateNewNefForA = true;

		double total_volume = 0;
		for(	it_patch_A =  Patch_Set_A.begin();
				it_patch_A != Patch_Set_A.end();
				++it_patch_A)
		{
			const libMesh::Elem * elem_A = m_Mesh_A.elem(*it_patch_A);
			bCreateNewNefForA = true;

			for(	it_patch_B =  Patch_Set_B.begin();
					it_patch_B != Patch_Set_B.end();
					++it_patch_B)
			{
				m_perf_log.push("Build exact intersections","Brute force algorithm");
				++nbOfTests;
				const libMesh::Elem * elem_B = m_Mesh_B.elem(*it_patch_B);

				intersection_vertices.clear();

//				std::cout << " hop ! " << *it_patch_A << " " << *it_patch_B << std::endl;
				bDoIntersect = m_Intersection_test.libMesh_exact_intersection_inside_coupling(elem_A,elem_B,intersection_vertices,bCreateNewNefForA);
				m_perf_log.pop("Build exact intersections","Brute force algorithm");

				if(bDoIntersect)
				{
					bCreateNewNefForA = false;
				}

				m_perf_log.push("Get volume","Brute force algorithm");
				if(bDoIntersect && intersection_vertices.size() >= 4)
				{
					total_volume += m_Mesh_Intersection.get_intersection_volume(intersection_vertices);
					m_Intersection_Pairs[nbOfPositiveTests] = std::pair<int,int>(*it_patch_A,*it_patch_B);
					m_DEBUG_Intersection_Pairs[*it_patch_A].insert(*it_patch_B);
					++nbOfPositiveTests;
				}
				m_perf_log.pop("Get volume","Brute force algorithm");
			}
		}
		std::cout << "    DEBUG: brute force search results" << std::endl;
		std::cout << " -> Positives / tests             : " << nbOfPositiveTests << " / " << nbOfTests
				  << " (" << 100.*nbOfPositiveTests/nbOfTests << "%)" << std::endl;
		std::cout << " -> Intersection / element volume : " << total_volume << " / " << Query_elem->volume() << std::endl << std::endl;
	}

	/*
	 * 		Find the first intersecting pair from a patch
	 *
	 * 		output: std::pair( guide_elem_id, probed_elem_id )
	 *
	 */
	void FindFirstPair(		Patch_construction * 	Patch_guide,
							Patch_construction * 	Patch_probed,
							std::pair<unsigned int,unsigned int> &		First_intersection)
	{
		// Set up some references for a simpler code
		std::unordered_set<unsigned int> & 	Patch_guide_ids  = Patch_guide->elem_indexes();
		std::unordered_set<unsigned int> & 	Patch_probed_ids = Patch_probed->elem_indexes();
		const libMesh::Mesh&  		Mesh_guide	= Patch_guide->mesh();
		const libMesh::Mesh&		Mesh_probed	= Patch_probed->mesh();

		// Set up locator for the first intersecting pair
		std::unique_ptr<libMesh::PointLocatorBase> Patch_guide_Locator = Mesh_guide.sub_point_locator();

		// Instruction needed to avoid the code from crashing if a query is
		// outside the mesh
		Patch_guide_Locator->enable_out_of_mesh_mode();

		// Find the first pair
		const libMesh::Elem * Patch_probed_first_elem = Mesh_probed.elem(*Patch_probed_ids.begin());
		unsigned int Patch_probed_first_elem_id = Patch_probed_first_elem->id();

		// Set of intersecting terms
		std::set<unsigned int> Intersecting_elems;
		m_Intersection_test.FindAllIntersection(Patch_probed_first_elem,Patch_guide_Locator,Intersecting_elems);

		// Search for one intersecting term inside the guide patch
		unsigned int Patch_guide_first_elem_id;
		bool bFoundFirstInter = false;

		std::cout << Intersecting_elems.size() << " : " << *Intersecting_elems.begin() << " " << Patch_probed_first_elem_id << std::endl;
		for(std::set<unsigned int>::iterator it_inter = Intersecting_elems.begin();
				it_inter != Intersecting_elems.end();
				++it_inter)
		{
			Patch_guide_first_elem_id = (*it_inter);
			if(Patch_guide_ids.find(Patch_guide_first_elem_id) != Patch_guide_ids.end())
			{
				// Found one !
				bFoundFirstInter = true;
				First_intersection.first = Patch_guide_first_elem_id;
				First_intersection.second = Patch_probed_first_elem_id;
				break;
			}
		}

		homemade_assert_msg(bFoundFirstInter, "Couldn't find a first intersecting pair!\n");
	};

	void BuildPatchIntersections_Front(const libMesh::Elem 	* Query_elem)
	{
		// TODO: For now it is more of a "find" than a build ...

		m_perf_log.push("Preamble","Advancing front algorithm");
		// Code for the front intersection tests
		m_Intersection_Pairs.clear();

		// Set containing all the intersection points
		std::set<Point_3> intersection_vertices;
		std::set<Point_3> dummy_intersection_vertices;

		// Intersection_Tools
		m_Intersection_test.libmesh_set_coupling_nef_polyhedron(Query_elem);
		m_Intersection_test_neighbors.libmesh_set_coupling_nef_polyhedron(Query_elem);

		// --- Set up variables needed by Gander's algorithm

		// First, determinate which patch will be the guide, and which will be
		// the probed
		Patch_construction * Patch_guide  = &m_Patch_Constructor_A;
		Patch_construction * Patch_probed = &m_Patch_Constructor_B;

		// Exchange them if the guide is smaller
		bool bGuidedByB = Patch_guide->size() > Patch_probed->size();
		if(bGuidedByB)
		{
			std::cout << "    DEBUG:     Probe: A     |     Guide: B" << std::endl;
			Patch_guide  = &m_Patch_Constructor_B;
			Patch_probed = &m_Patch_Constructor_A;
		}
		else
		{
			std::cout << "    DEBUG:     Probe: B     |     Guide: A" << std::endl;
		}

		// Ids of the working elements
		unsigned int Guide_working_elem_id;
		unsigned int Probed_working_elem_id;

		// Find the first pair, (guide, probed)
		std::pair<unsigned int,unsigned int> First_intersection;
		FindFirstPair(Patch_guide,Patch_probed,First_intersection);

		// Index and iterators guide neighbor elements to be tested
		unsigned int guide_neighbor_index = 0;
		std::unordered_set<unsigned int>::iterator it_neigh, it_neigh_end;

		// Boolean saying if two elements intersect (exactly)
		bool bDoIntersect = false;

		// Boolean saying if a neighbor element intersects (exactly)
		bool bDoNeighborIntersect = false;

		// Boolean saying if the two elements intersect inside the coupling
		bool bDoIntersectInsideCoupling = false;

		// Boolean used to short-circuit the neighbor search test if all
		// the concerned elements already have neighbours.
		bool bDoTestGuideNeighbours = true;

		// --- Preamble - initializations
		// Insert the first elements in the queues
		Patch_guide->FrontSearch_initialize();
		Patch_probed->FrontSearch_initialize();

		Patch_guide->intersection_queue_push_back(First_intersection.first);
		Patch_probed->intersection_queue_push_back(First_intersection.second);
		Patch_guide->set_elem_as_inside_queue(First_intersection.first);

		// Debug vars
		int nbOfGuideElemsTested = 0;
		int nbOfTests = 0;
		int nbOfPositiveTests = 0;
		double total_volume = 0;

		bool bCreateNewNefForGuide = true;

		bool bCreateNewNefForGuideNeigh = true;

		m_perf_log.pop("Preamble","Advancing front algorithm");

		while(!Patch_guide->intersection_queue_empty())
		{
			m_perf_log.push("Set up new guide element","Advancing front algorithm");

			// Pop out the first element from the guide intersection queue
			Guide_working_elem_id = Patch_guide->intersection_queue_extract_front_elem();
			Patch_guide->set_elem_as_treated(Guide_working_elem_id);

			// Get its pointer
			const libMesh::Elem * elem_guide = Patch_guide->current_elem_pointer();

			++nbOfGuideElemsTested;

			// Set up the neighbor test short-circuits
			bDoTestGuideNeighbours = Patch_guide->set_neighbors_to_search_next_pairs();

			// Clear and initialize the probed patch lists
			Patch_probed->FrontSearch_reset();

			// Pop out the first element from the probed intersection queue, and
			// insert it inside the test queue
			Patch_probed->FrontSearch_prepare_for_probed_test();

			bCreateNewNefForGuide = true;

			m_perf_log.pop("Set up new guide element","Advancing front algorithm");

			while(!Patch_probed->test_queue_empty())
			{
				m_perf_log.push("Search probe intersections","Advancing front algorithm");
				// Pop out the first elem to be tested
				Probed_working_elem_id = Patch_probed->test_queue_extract_front_elem();
				Patch_probed->set_elem_as_treated(Probed_working_elem_id);

				const libMesh::Elem * elem_probed = Patch_probed->current_elem_pointer();

				// Test the intersection
				bDoIntersect = m_Intersection_test.libMesh_exact_do_intersect(elem_guide,elem_probed);
				++nbOfTests;
				m_perf_log.pop("Search probe intersections","Advancing front algorithm");

				// If they do intersect, we have to ...
				if(bDoIntersect)
				{
					m_perf_log.push("Build exact intersections","Advancing front algorithm");
					// 1) test if the intersection is inside the coupling
					//    element. If true, process it.
					bDoIntersectInsideCoupling =  m_Intersection_test.libMesh_exact_intersection_inside_coupling(elem_guide,elem_probed,intersection_vertices,bCreateNewNefForGuide,false,false);
					m_perf_log.pop("Build exact intersections","Advancing front algorithm");

					if(bDoIntersectInsideCoupling)
					{
						bCreateNewNefForGuide = false;
					}

					m_perf_log.push("Get volume","Advancing front algorithm");
					if(bDoIntersectInsideCoupling && intersection_vertices.size() >= 4)
					{
						// TODO For now, just add the volume
						total_volume += m_Mesh_Intersection.get_intersection_volume(intersection_vertices);
						m_DEBUG_Intersection_Pairs[Guide_working_elem_id].insert(Probed_working_elem_id);
						++nbOfPositiveTests;
					}
					m_perf_log.pop("Get volume","Advancing front algorithm");

					m_perf_log.push("Update test queue","Advancing front algorithm");
					// 2) add elem_probed's neighbors, if not treated yet
					Patch_probed->add_neighbors_to_test_list();
					m_perf_log.pop("Update test queue","Advancing front algorithm");

					m_perf_log.push("Update intersection queue","Advancing front algorithm");
					// 3) determinate if any of the guide neighbors are good
					//    candidates with the probed element.
					if(bDoTestGuideNeighbours)
					{
						it_neigh = Patch_guide->neighbors_to_search_next_pair().begin();
						it_neigh_end = Patch_guide->neighbors_to_search_next_pair().end();

						bCreateNewNefForGuideNeigh = true;

						for( ; it_neigh != it_neigh_end; ++it_neigh )
						{
							// This element doesn't have an intersecting pair
							// yet. Test it, but only save it if it is inside
							// the coupling region
							guide_neighbor_index = *it_neigh;
							const libMesh::Elem * elem_neigh = Patch_guide->elem(guide_neighbor_index);
							bDoNeighborIntersect = m_Intersection_test_neighbors.libMesh_exact_do_intersect_inside_coupling(elem_probed,elem_neigh,bCreateNewNefForGuideNeigh);

							if(bDoNeighborIntersect)
							{
								bCreateNewNefForGuideNeigh = false;
							}

							if(bDoNeighborIntersect)
							{
								// Now it has an intersecting pair!
								Patch_guide->intersection_queue_push_back(guide_neighbor_index);
								Patch_probed->intersection_queue_push_back(Probed_working_elem_id);

								// Set the guide elements as "already inside the queue"
								Patch_guide->set_elem_as_inside_queue(guide_neighbor_index);
							}
						}

						bDoTestGuideNeighbours = Patch_guide->set_neighbors_to_search_next_pairs();
					}
					m_perf_log.pop("Update intersection queue","Advancing front algorithm");
				}
			}
		}

		std::cout << "    DEBUG: advancing front search results" << std::endl;
		std::cout << " -> Guide elements tested / all   : " << nbOfGuideElemsTested << " / " << Patch_guide->size() << std::endl;
		std::cout << " -> Positives / tests             : " << nbOfPositiveTests << " / " << nbOfTests
				  << " (" << 100.*nbOfPositiveTests/nbOfTests << "%)" << std::endl;
		std::cout << " -> Intersection / element volume : " << total_volume << " / " << Query_elem->volume() << std::endl << std::endl;
	}

	void BuildIntersections(SearchMethod search_type = BRUTE)
	{
		// Prepare iterators
		libMesh::Mesh::const_element_iterator it_coupl = m_Mesh_Coupling.elements_begin();
		libMesh::Mesh::const_element_iterator it_coupl_end = m_Mesh_Coupling.elements_end();

		int patch_counter = 0;

		for( ; it_coupl != it_coupl_end; ++it_coupl)
		{
			const libMesh::Elem * Query_elem = * it_coupl;

			m_perf_log.push("Build patches");
			BuildCoupledPatches(Query_elem,patch_counter);
			m_perf_log.pop("Build patches");

			++patch_counter;

			switch (search_type)
			{
				case BRUTE : 	m_perf_log.push("Find intersections - brute");
								m_DEBUG_Intersection_Pairs.clear();
								BuildPatchIntersections_Brute(Query_elem);
								m_perf_log.pop("Find intersections - brute");

								break;
				case FRONT : 	m_perf_log.push("Find intersections - advancing front");
								m_DEBUG_Intersection_Pairs.clear();
								BuildPatchIntersections_Front(Query_elem);
								m_perf_log.pop("Find intersections - advancing front");
								break;
				case BOTH : 	m_perf_log.push("Find intersections - brute");
								m_DEBUG_Intersection_Pairs.clear();
								BuildPatchIntersections_Brute(Query_elem);
								m_perf_log.pop("Find intersections - brute");

								m_perf_log.push("Find intersections - advancing front");
								m_DEBUG_Intersection_Pairs.clear();
								BuildPatchIntersections_Front(Query_elem);
								m_perf_log.pop("Find intersections - advancing front");
								break;
			}
		}
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_ */
