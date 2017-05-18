/*
 * intersection_search.cpp
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "intersection_search.h"

namespace carl
{
	libMesh::Mesh & Intersection_Search::mesh_A()
	{
		return m_Mesh_A;
	}

	libMesh::Mesh & Intersection_Search::mesh_B()
	{
		return m_Mesh_B;
	}

	libMesh::Mesh & Intersection_Search::mesh_Coupling()
	{
		return m_Mesh_Coupling;
	}

	/*
	 * 		Build both patches associated to the query element
	 */
	void Intersection_Search::BuildCoupledPatches(const libMesh::Elem 	* Query_elem, int patch_counter)
	{
		// Unbreakable Patches!
		m_Patch_Constructor_A.BuildPatch(Query_elem);

		// Trusty Patches!
		m_Patch_Constructor_B.BuildPatch(Query_elem);

		if(m_bPrintDebug)
		{
			std::string filename = "meshes/3D/tests/output/patch_mesh_A_" + std::to_string(patch_counter) + "_" + std::to_string(m_rank);
			m_Patch_Constructor_A.export_patch_mesh(filename);
			filename = "meshes/3D/tests/output/patch_mesh_B_" + std::to_string(patch_counter) + "_" + std::to_string(m_rank);
			m_Patch_Constructor_B.export_patch_mesh(filename);
		}
	}

	/*
	 * 		Find all the intersections between the patches, using a brute force
	 * 	method (all elements from a patch are tested against all the elements
	 * 	from the other patch).
	 */
	void Intersection_Search::FindPatchIntersections_Brute(const libMesh::Elem 	* Query_elem)
	{
		m_perf_log.push("Preamble","Brute force algorithm");

		// Code for the brute force intersection tests
		m_Intersection_Pairs_multimap.clear();

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

		for(	it_patch_A =  Patch_Set_A.begin();
				it_patch_A != Patch_Set_A.end();
				++it_patch_A)
		{
			const libMesh::Elem * elem_A = m_Mesh_A.elem(*it_patch_A);

			for(	it_patch_B =  Patch_Set_B.begin();
					it_patch_B != Patch_Set_B.end();
					++it_patch_B)
			{
				m_perf_log.push("Find intersection","Brute force algorithm");
				++nbOfTests;
				const libMesh::Elem * elem_B = m_Mesh_B.elem(*it_patch_B);

				bDoIntersect = m_Intersection_test.libMesh_exact_do_intersect(elem_A,elem_B);
				m_perf_log.pop("Find intersection","Brute force algorithm");

				if(bDoIntersect)
				{
					if(m_bSaveInterData)
					{
						m_Intersection_Pairs_multimap.insert(std::pair<unsigned int, unsigned int>(*it_patch_A,*it_patch_B));
					}
					++nbOfPositiveTests;
				}
			}
		}

		// Partitioning vars
		unsigned int query_idx = Query_elem->id();
		m_Nb_Of_Intersections_Elem_C[query_idx] = nbOfPositiveTests;

		if(m_bPrintDebug)
		{
			std::cout << "    DEBUG: brute force search results" << std::endl;
			std::cout << " -> Positives / tests             : " << nbOfPositiveTests << " / " << nbOfTests
					  << " (" << 100.*nbOfPositiveTests/nbOfTests << "%)" << std::endl  << std::endl;
		}
	}

	/*
	 * 		Find the first intersecting pair from a patch. Do so by :
	 *
	 * 		a) locating the elements from the guide patch that contain the
	 * 		vertices of an element from the probed patch.
	 * 		b) for each one of these guide elements, find the one that
	 * 		intersects the probed element inside the coupling.
	 * 		c) if this fails, use a brute force search method.
	 */
	void Intersection_Search::FindFirstPair(		Patch_construction * 	Patch_guide,
							Patch_construction * 	Patch_probed,
							std::pair<unsigned int,unsigned int> &		First_intersection)
	{
		// Set up some references for a simpler code
		libMesh::ReplicatedMesh&  	Mesh_patch_guide	= Patch_guide->patch_mesh();
		libMesh::ReplicatedMesh&		Mesh_patch_probed	= Patch_probed->patch_mesh();

		// Set up locator for the first intersecting pair
		std::unique_ptr<libMesh::PointLocatorBase> Patch_guide_Locator = Mesh_patch_guide.sub_point_locator();

		// Instruction needed to avoid the code from crashing if a query is
		// outside the mesh
		Patch_guide_Locator->enable_out_of_mesh_mode();

		// Find the first pair
		libMesh::Elem * Patch_probed_first_elem = NULL;
		unsigned int Patch_probed_first_elem_id = 0;

		// Set of intersecting terms:
		std::set<unsigned int> Intersecting_guide_elems;
		bool bFoundIntersection = false;
		for(libMesh::ReplicatedMesh::element_iterator element_probed_it = Mesh_patch_probed.elements_begin();
				 element_probed_it != Mesh_patch_probed.elements_end();
				 ++element_probed_it)
		{
			Patch_probed_first_elem = * element_probed_it;
			Patch_probed_first_elem_id = Patch_probed->convert_patch_to_parent_elem_id(Patch_probed_first_elem->id());
			bFoundIntersection = m_Intersection_test.FindAllIntersection(Patch_probed_first_elem,Patch_guide_Locator,Intersecting_guide_elems);
			if(bFoundIntersection)
			{
				break;
			}
		}

		// Search for one intersecting term inside the guide patch
		libMesh::Elem * Patch_guide_first_elem = NULL;
		unsigned int Patch_guide_first_elem_id;
		bool bFoundFirstInter = false;

		for(std::set<unsigned int>::iterator it_inter = Intersecting_guide_elems.begin();
				it_inter != Intersecting_guide_elems.end();
				++it_inter)
		{
			Patch_guide_first_elem = Mesh_patch_guide.elem(*it_inter);
			Patch_guide_first_elem_id = Patch_guide->convert_patch_to_parent_elem_id(*it_inter);

			// Must test if the intersection is inside the coupling region
			bFoundFirstInter = m_Intersection_test_neighbors.libMesh_exact_do_intersect_inside_coupling(Patch_probed_first_elem,Patch_guide_first_elem);

			if(bFoundFirstInter)
			{
				First_intersection.first = Patch_guide_first_elem_id;
				First_intersection.second = Patch_probed_first_elem_id;
				break;
			}
		}

		if(!bFoundFirstInter)
		{
			// Couldn't find the first intersection using a point search ...
			// Will have to do the things the hard way ...
			BruteForce_FindFirstPair(Patch_guide,Patch_probed,First_intersection);
		}
	};

	/*
	 * 		Find the first intersecting pair from a patch, doing a full scan of
	 * 	the patches (essentially, a brute force algorithm set to stop after the
	 * 	first positive test).
	 */
	void Intersection_Search::BruteForce_FindFirstPair(	Patch_construction * 	Patch_guide,
									Patch_construction * 	Patch_probed,
									std::pair<unsigned int,unsigned int> &		First_intersection)
	{
		bool bFoundFirstInter = false;
		std::unordered_set<unsigned int> & Patch_Set_Guide = Patch_guide->elem_indexes();
		std::unordered_set<unsigned int> & Patch_Set_Probed = Patch_probed->elem_indexes();

		libMesh::ReplicatedMesh&  		Mesh_patch_guide	= Patch_guide->patch_mesh();
		libMesh::ReplicatedMesh&		Mesh_patch_probed	= Patch_probed->patch_mesh();

		std::unordered_set<unsigned int>::iterator it_patch_Guide;
		std::unordered_set<unsigned int>::iterator it_patch_Probed;

		for(	it_patch_Guide =  Patch_Set_Guide.begin();
				it_patch_Guide != Patch_Set_Guide.end();
				++it_patch_Guide)
		{
			const libMesh::Elem * elem_Guide = Mesh_patch_guide.elem(*it_patch_Guide);

			for(	it_patch_Probed =  Patch_Set_Probed.begin();
					it_patch_Probed != Patch_Set_Probed.end();
					++it_patch_Probed)
			{
				const libMesh::Elem * elem_Probed = Mesh_patch_probed.elem(*it_patch_Probed);

				bFoundFirstInter = m_Intersection_test.libMesh_exact_do_intersect_inside_coupling(elem_Guide,elem_Probed);

				if(bFoundFirstInter)
				{
					First_intersection.first = *it_patch_Guide;
					First_intersection.second = *it_patch_Probed;
					return;
				}
			}
		}

		homemade_assert_msg(bFoundFirstInter,"Could't find a first intersecting pair. Are you sure that the patches do intersect?\n");
	}

	/*
	 * 		Find all the intersections between the patches, using an advancing
	 * 	front method.
	 */
	void Intersection_Search::FindPatchIntersections_Front(const libMesh::Elem 	* Query_elem)
	{
		m_perf_log.push("Preamble","Advancing front algorithm");
		// Code for the front intersection tests
		m_Intersection_Pairs_multimap.clear();

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
			if(m_bPrintDebug)
			{
				std::cout << "    DEBUG:     Probe: A     |     Guide: B" << std::endl;
			}
			Patch_guide  = &m_Patch_Constructor_B;
			Patch_probed = &m_Patch_Constructor_A;
		}
		else
		{
			if(m_bPrintDebug)
			{
				std::cout << "    DEBUG:     Probe: B     |     Guide: A" << std::endl;
			}
		}

		// Ids of the working elements
		unsigned int Guide_working_elem_id;
		unsigned int Probed_working_elem_id;

		// Find the first pair, (guide, probed)
		std::pair<unsigned int,unsigned int> First_intersection;
		BruteForce_FindFirstPair(Patch_guide,Patch_probed,First_intersection);

		// Index and iterators guide neighbor elements to be tested
		unsigned int guide_neighbor_index = 0;
		std::unordered_set<unsigned int>::iterator it_neigh, it_neigh_end;

		// Boolean saying if two elements intersect (exactly)
		bool bDoIntersect = false;

		// Boolean saying if a neighbor element intersects (exactly)
		bool bDoNeighborIntersect = false;

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

			m_perf_log.pop("Set up new guide element","Advancing front algorithm");

			while(!Patch_probed->test_queue_empty())
			{
				m_perf_log.push("Find intersection","Advancing front algorithm");
				// Pop out the first elem to be tested
				Probed_working_elem_id = Patch_probed->test_queue_extract_front_elem();
				Patch_probed->set_elem_as_treated(Probed_working_elem_id);

				const libMesh::Elem * elem_probed = Patch_probed->current_elem_pointer();

				// Test the intersection
				bDoIntersect = m_Intersection_test.libMesh_exact_do_intersect(elem_guide,elem_probed);
				++nbOfTests;
				m_perf_log.pop("Find intersection","Advancing front algorithm");

				// If they do intersect, we have to ...
				if(bDoIntersect)
				{
					// 1) add elements to intersection multimap, watching out
					//    for the element order
					if(m_bSaveInterData)
					{
						if(bGuidedByB)
						{
							m_Intersection_Pairs_multimap.insert(std::pair<unsigned int, unsigned int>(Probed_working_elem_id,Guide_working_elem_id));
						}
						else
						{
							m_Intersection_Pairs_multimap.insert(std::pair<unsigned int, unsigned int>(Guide_working_elem_id,Probed_working_elem_id));
						}
					}
					++nbOfPositiveTests;

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

		// Partitioning vars
		unsigned int query_idx = Query_elem->id();
		m_Nb_Of_Intersections_Elem_C[query_idx] = nbOfPositiveTests;

		if(m_bPrintDebug)
		{
			std::cout << "    DEBUG: advancing front search results" << std::endl;
			std::cout << " -> Guide elements tested / all   : " << nbOfGuideElemsTested << " / " << Patch_guide->size() << std::endl;
			std::cout << " -> Positives / tests             : " << nbOfPositiveTests << " / " << nbOfTests
					  << " (" << 100.*nbOfPositiveTests/nbOfTests << "%)" << std::endl  << std::endl;
		}
	}

	/*
	 * 		For each coupling element, build the patches and find their
	 * 	intersections, using the brute force method.
	 */
	void Intersection_Search::FindAndBuildIntersections_Brute()
	{
		// Prepare iterators
		libMesh::Mesh::const_element_iterator it_coupl = m_Mesh_Coupling.local_elements_begin();
		libMesh::Mesh::const_element_iterator it_coupl_end = m_Mesh_Coupling.local_elements_end();

		int patch_counter = 0;
		double real_volume = 0;

		m_Mesh_Intersection.initialize();

		for( ; it_coupl != it_coupl_end; ++it_coupl)
		{
			const libMesh::Elem * Query_elem = * it_coupl;

			m_perf_log.push("Build patches","Brute force");
			BuildCoupledPatches(Query_elem,patch_counter);
			m_perf_log.pop("Build patches","Brute force");

			++patch_counter;

			if(m_bPrintDebug)
			{
				real_volume += Query_elem->volume();
			}

			m_perf_log.push("Find intersections","Brute force");
			FindPatchIntersections_Brute(Query_elem);
			m_perf_log.pop("Find intersections","Brute force");
			if(!m_bSkipIntersectionConstruction)
			{
				m_perf_log.push("Build intersections","Brute force");
				UpdateCouplingIntersection(Query_elem);
				m_perf_log.pop("Build intersections","Brute force");
			}
		};

		if(!m_bSkipIntersectionConstruction)
		{
			m_perf_log.push("Prepare mesh","Brute force");
			m_Mesh_Intersection.prepare_for_use();
			m_perf_log.pop("Prepare mesh","Brute force");
		}

		if(m_bPrintTimingData)
		{
			std::vector<double> timing_find_intersections(m_nodes,0);
			std::vector<double> timing_build_intersections(m_nodes,0);

			libMesh::PerfData performance_data = m_perf_log.get_perf_data("Find intersections","Brute force");
			timing_find_intersections[m_rank] = performance_data.tot_time_incl_sub;
			performance_data = m_perf_log.get_perf_data("Build intersections","Brute force");
			timing_build_intersections[m_rank] = performance_data.tot_time_incl_sub;

			m_comm.sum(timing_find_intersections);
			m_comm.sum(timing_build_intersections);

			if(m_rank == 0)
			{
				print_stats_to_file(timing_find_intersections,m_timing_data_file_base + "_search_brute.dat");
				if(!m_bSkipIntersectionConstruction)
				{
					print_stats_to_file(timing_build_intersections,m_timing_data_file_base + "_build_brute.dat");
				}
			}
		}

		if(m_bPrintDebug && !m_bSkipIntersectionConstruction)
		{
			m_perf_log.push("Calculate volume","Brute force");
			double total_volume = m_Mesh_Intersection.get_total_volume();
			m_perf_log.push("Calculate volume","Brute force");

			std::cout << "    DEBUG volume on proc. " << m_rank << "(BRUTE)" << std::endl;
			std::cout << " -> Mesh elems, nodes             : " << m_Mesh_Intersection.mesh().n_elem() << " , " << m_Mesh_Intersection.mesh().n_nodes() << std::endl;
			std::cout << " -> Intersection volume / real    : " << total_volume << " / " << real_volume << std::endl << std::endl;
		}
	}

	void Intersection_Search::FindIntersections_Brute()
	{
		// Prepare iterators
		libMesh::Mesh::const_element_iterator it_coupl = m_Mesh_Coupling.local_elements_begin();
		libMesh::Mesh::const_element_iterator it_coupl_end = m_Mesh_Coupling.local_elements_end();

		int patch_counter = 0;

		for( ; it_coupl != it_coupl_end; ++it_coupl)
		{
			const libMesh::Elem * Query_elem = * it_coupl;

			m_perf_log.push("Build patches","Brute force (preamble run)");
			BuildCoupledPatches(Query_elem,patch_counter);
			m_perf_log.pop("Build patches","Brute force (preamble run)");

			++patch_counter;

			m_perf_log.push("Find intersections","Brute force (preamble run)");
			FindPatchIntersections_Brute(Query_elem);
			m_perf_log.pop("Find intersections","Brute force (preamble run)");
		};
	}

	/*
	 * 		For each coupling element, build the patches and find their
	 * 	intersections, using the advancing front method.
	 */
	void Intersection_Search::FindAndBuildIntersections_Front()
	{
		// Prepare iterators
		libMesh::Mesh::const_element_iterator it_coupl = m_Mesh_Coupling.local_elements_begin();
		libMesh::Mesh::const_element_iterator it_coupl_end = m_Mesh_Coupling.local_elements_end();

		int patch_counter = 0;
		double real_volume = 0;

		m_Mesh_Intersection.initialize();

		for( ; it_coupl != it_coupl_end; ++it_coupl)
		{
			const libMesh::Elem * Query_elem = * it_coupl;

			m_perf_log.push("Build patches","Advancing front");
			BuildCoupledPatches(Query_elem,patch_counter);
			m_perf_log.pop("Build patches","Advancing front");

			++patch_counter;

			if(m_bPrintDebug)
			{
				real_volume += Query_elem->volume();
			}
			m_perf_log.push("Find intersections","Advancing front");
			FindPatchIntersections_Front(Query_elem);
			m_perf_log.pop("Find intersections","Advancing front");
			if(!m_bSkipIntersectionConstruction)
			{
				m_perf_log.push("Build intersections","Advancing front");
				UpdateCouplingIntersection(Query_elem);
				m_perf_log.pop("Build intersections","Advancing front");
			}
		};

		if(!m_bSkipIntersectionConstruction)
		{
			m_perf_log.push("Prepare mesh","Advancing front");
			m_Mesh_Intersection.prepare_for_use();
			m_perf_log.pop("Prepare mesh","Advancing front");
		}

		if(m_bPrintTimingData)
		{
			std::vector<double> timing_find_intersections(m_nodes,0);
			std::vector<double> timing_build_intersections(m_nodes,0);

			libMesh::PerfData performance_data = m_perf_log.get_perf_data("Find intersections","Advancing front");
			timing_find_intersections[m_rank] = performance_data.tot_time_incl_sub;
			performance_data = m_perf_log.get_perf_data("Build intersections","Advancing front");
			timing_build_intersections[m_rank] = performance_data.tot_time_incl_sub;

			m_comm.sum(timing_find_intersections);
			m_comm.sum(timing_build_intersections);

			if(m_rank == 0)
			{
				print_stats_to_file(timing_find_intersections,m_timing_data_file_base + "_search_advancing.dat");
				if(!m_bSkipIntersectionConstruction)
				{
					print_stats_to_file(timing_build_intersections,m_timing_data_file_base + "_build_advancing.dat");
				}
			}
		}

		if(m_bPrintDebug && !m_bSkipIntersectionConstruction)
		{
			m_perf_log.push("Calculate volume","Advancing front");
			double total_volume = m_Mesh_Intersection.get_total_volume();
			m_perf_log.push("Calculate volume","Advancing front");

			std::cout << "    DEBUG volume on proc. " << m_rank << "(FRONT)" << std::endl;
			std::cout << " -> Mesh elems, nodes             : " << m_Mesh_Intersection.mesh().n_elem() << " , " << m_Mesh_Intersection.mesh().n_nodes() << std::endl;
			std::cout << " -> Intersection volume / real    : " << total_volume << " / " << real_volume << std::endl << std::endl;
		}
	}

	void Intersection_Search::FindIntersections_Front()
	{
		// Prepare iterators
		libMesh::Mesh::const_element_iterator it_coupl = m_Mesh_Coupling.local_elements_begin();
		libMesh::Mesh::const_element_iterator it_coupl_end = m_Mesh_Coupling.local_elements_end();

		int patch_counter = 0;

		for( ; it_coupl != it_coupl_end; ++it_coupl)
		{
			const libMesh::Elem * Query_elem = * it_coupl;

			m_perf_log.push("Build patches (preamble)","Advancing front (preamble run)");
			BuildCoupledPatches(Query_elem,patch_counter);
			m_perf_log.pop("Build patches","Advancing front (preamble run)");

			++patch_counter;

			m_perf_log.push("Find intersections","Advancing front (preamble run)");
			FindPatchIntersections_Front(Query_elem);
			m_perf_log.pop("Find intersections","Advancing front (preamble run)");
		};
	}

	/*
	 * 		Preallocate run. It essentially does the intersection run, but
	 * 	without saving the data or building the intersections themselves.
	 *
	 */
	void Intersection_Search::PreparePreallocationAndLoad(SearchMethod search_type)
	{
		m_bSaveInterData = false;
		switch (search_type)
		{
			case SearchMethod::BRUTE :	FindIntersections_Brute();
							break;

			case SearchMethod::FRONT :	FindIntersections_Front();
							break;

			case SearchMethod::BOTH :		FindIntersections_Brute();
							FindIntersections_Front();
							break;
		}

		m_comm.sum(m_Nb_Of_Intersections_Elem_C);

		m_bDidPreliminarySearch = true;
	}

	void Intersection_Search::PreallocateAndPartitionCoupling()
	{
		if(m_bDidPreliminarySearch)
		{
			// Redo the partitioning
			m_coupling_weights.resize(m_Nb_Of_Intersections_Elem_C.size());

			for(unsigned int iii = 0; iii < m_Nb_Of_Intersections_Elem_C.size(); ++iii)
			{
				m_coupling_weights[iii] = m_Nb_Of_Intersections_Elem_C[iii];
			}

			if(!m_bSkipIntersectionPartitioning)
			{
				libMesh::Partitioner * dummy_partitioner = m_Mesh_Coupling.partitioner().get();
				dummy_partitioner->attach_weights(&m_coupling_weights);
				m_Mesh_Coupling.partition(m_nodes);
			}

			// Allocate the intersection maps ...
			std::vector<unsigned int> nb_of_inters_per_rank(m_nodes,0);
			unsigned int dummy_nb_of_inters = 0;
			libMesh::Mesh::const_element_iterator it_local = m_Mesh_Coupling.local_elements_begin();
			libMesh::Mesh::const_element_iterator it_local_end = m_Mesh_Coupling.local_elements_end();
			for( ; it_local != it_local_end; ++ it_local)
			{
				const libMesh::Elem * dummy_elem = * it_local;
				nb_of_inters_per_rank[m_rank] += m_coupling_weights[dummy_elem->id()];
			}
			m_comm.sum(nb_of_inters_per_rank);
			m_Intersection_Pairs_multimap.reserve(2*dummy_nb_of_inters);

			// ... and the intersection grid
			m_Mesh_Intersection.preallocate_grid(192*dummy_nb_of_inters);
			m_bHavePreallocData = true;

			if(m_bPrintIntersectionsPerPartData)
			{
				if(m_rank == 0)
				{
					std::string intersection_statistics_filename;
					if(!m_bSkipIntersectionPartitioning)
					{
						intersection_statistics_filename = m_timing_data_file_base + "_inters_per_partition__new_partitioning.dat";
					}
					else
					{
						intersection_statistics_filename = m_timing_data_file_base + "_inters_per_partition__original.dat";
					}

					libMesh::StatisticsVector<int> intersection_statistics(m_nodes,0);
					for(unsigned int iii = 0; iii < m_nodes; ++iii)
					{
						intersection_statistics[iii] = nb_of_inters_per_rank[iii];
					}

					std::ofstream inter_statistics_output(intersection_statistics_filename,std::ofstream::app);

					inter_statistics_output 	<< intersection_statistics.minimum() << " "
									<< intersection_statistics.maximum() << " "
									<< intersection_statistics.mean() << " "
									<< intersection_statistics.median() << " "
									<< intersection_statistics.stddev() << std::endl;

					inter_statistics_output.close();
				}

			}
		}
	}

	/*
	 * 		Interface for the user to build the intersections. By default, it
	 * 	uses the brute force algorithm, but the argument can be changed to
	 * 	carl::FRONT to use the advancing front method, of to carl::BOTH to use
	 * 	both methods (useful for benchmarking).
	 */
	void Intersection_Search::BuildIntersections(SearchMethod search_type)
	{
		if(!m_bHavePreallocData)
		{
			// Must do an preeeety expensive preallocation. Ouch ...
			m_Intersection_Pairs_multimap.reserve(m_Mesh_A.n_elem()*m_Mesh_B.n_elem());
		}

		m_bSaveInterData = true;
		m_bIntersectionsBuilt = false;

		switch (search_type)
		{
			case SearchMethod::BRUTE :	FindAndBuildIntersections_Brute();
							break;

			case SearchMethod::FRONT :	FindAndBuildIntersections_Front();
							break;

			case SearchMethod::BOTH :		FindAndBuildIntersections_Brute();
							FindAndBuildIntersections_Front();
							break;
		}

		if(!m_bSkipIntersectionConstruction)
		{
			m_Mesh_Intersection.export_intersection_data(m_Output_filename_base);
			m_bIntersectionsBuilt = true;
		}
	}

	/*
	 * 		Calculate the volume over all the processors
	 */
	void Intersection_Search::CalculateGlobalVolume()
	{
		// Use a barrier to guarantee that all procs are in the same position
		m_comm.barrier();

		// Calculate the volume of the intersection on each processor
		double global_volume = m_Mesh_Intersection.get_total_volume();

		// Add it!
		m_comm.sum(global_volume);

		// Calculate the volume of the coupling region on each processor
		libMesh::Mesh::const_element_iterator it_coupl = m_Mesh_Coupling.local_elements_begin();
		libMesh::Mesh::const_element_iterator it_coupl_end = m_Mesh_Coupling.local_elements_end();
		double real_volume = 0;

		for( ; it_coupl != it_coupl_end; ++it_coupl)
		{
			const libMesh::Elem * Query_elem = * it_coupl;
			real_volume += Query_elem->volume();
		}

		// Add it!
		m_comm.sum(real_volume);

		std::cout << "    TOTAL volume, proc. " << m_rank << std::endl;
		std::cout << " -> Intersection volume / real    : " << global_volume << " / " << real_volume << std::endl << std::endl;

	}
	/*
	 * 		Legacy function, used to calculate the volume of the intersections
	 * 	without updating the intersection mesh.
	 */
	void Intersection_Search::CalculateIntersectionVolume(const libMesh::Elem 	* Query_elem)
	{
		std::unordered_set<unsigned int> & Patch_Set_A = m_Patch_Constructor_A.elem_indexes();

		std::unordered_set<unsigned int>::iterator it_patch_A;

		bool bDoIntersect = true;
		bool bCreateNewNefForA = true;

		m_Intersection_test.libmesh_set_coupling_nef_polyhedron(Query_elem);

		std::set<libMesh::Point> points_out;
		double total_volume = 0;

		// Debug vars
		int nbOfTests = 0;
		int nbOfPositiveTests = 0;
		for(	it_patch_A =  Patch_Set_A.begin();
				it_patch_A != Patch_Set_A.end();
				++it_patch_A)
		{
			auto iterator_pair = m_Intersection_Pairs_multimap.equal_range(*it_patch_A);

			const libMesh::Elem * elem_A = m_Mesh_A.elem(*it_patch_A);
			bCreateNewNefForA = true;

			for(auto it_patch_B = iterator_pair.first ;
					 it_patch_B != iterator_pair.second ;
					 ++it_patch_B )
			{
				const libMesh::Elem * elem_B = m_Mesh_B.elem(it_patch_B->second);
				points_out.clear();

				m_perf_log.push("Build intersection polyhedron","Build intersections");
				bDoIntersect = m_Intersection_test.libMesh_exact_intersection_inside_coupling(elem_A,elem_B,points_out,bCreateNewNefForA,true,false);
				m_perf_log.pop("Build intersection polyhedron","Build intersections");
				++nbOfTests;

				if(bDoIntersect)
				{
					bCreateNewNefForA = false;
					m_perf_log.push("Calculate volume","Build intersections");
					total_volume += m_Mesh_Intersection.get_intersection_volume(points_out);
					m_perf_log.pop("Calculate volume","Build intersections");
					++nbOfPositiveTests;
				}
			}
		}

		if(m_bPrintDebug)
		{
			std::cout << "    DEBUG: calculate the volume" << std::endl;
			std::cout << " -> Intersection volume / real    : " << total_volume << " / " << Query_elem->volume() << std::endl << std::endl;
			std::cout << " -> Positives / tests             : " << nbOfPositiveTests << " / " << nbOfTests
					  << " (" << 100.*nbOfPositiveTests/nbOfTests << "%)" << std::endl  << std::endl;
		}
	}

	/*
	 * 		Take the intersection tables info and update the intersection mesh.
	 */
	void Intersection_Search::UpdateCouplingIntersection(const libMesh::Elem 	* Query_elem)
	{
		// Preamble
		std::unordered_set<unsigned int> & Patch_Set_A = m_Patch_Constructor_A.elem_indexes();
		std::unordered_set<unsigned int>::iterator it_patch_A;

		bool bDoIntersect = true;
		bool bCreateNewNefForA = true;
		m_Intersection_test.libmesh_set_coupling_nef_polyhedron(Query_elem);

		std::set<libMesh::Point> points_out;

		// Debug vars
		int nbOfTests = 0;
		int nbOfPositiveTests = 0;

		// For each element inside patch A ...
		for(	it_patch_A =  Patch_Set_A.begin();
				it_patch_A != Patch_Set_A.end();
				++it_patch_A)
		{
			auto iterator_pair = m_Intersection_Pairs_multimap.equal_range(*it_patch_A);

			const libMesh::Elem * elem_A = m_Mesh_A.elem(*it_patch_A);
			bCreateNewNefForA = true;

			// ... get the intersections with patch B ...
			for(auto it_patch_B = iterator_pair.first ;
					 it_patch_B != iterator_pair.second ;
					 ++it_patch_B )
			{
				const libMesh::Elem * elem_B = m_Mesh_B.elem(it_patch_B->second);
				points_out.clear();

				// ... test if they intersect inside the coupling region ...
				m_perf_log.push("Build intersection polyhedron","Build intersections");
				bDoIntersect = m_Intersection_test.libMesh_exact_intersection_inside_coupling(elem_A,elem_B,points_out,bCreateNewNefForA,true,false);
				m_perf_log.pop("Build intersection polyhedron","Build intersections");
				++nbOfTests;
				
				if(bDoIntersect)
				{
					// ... and, if they do, update the intersection mesh!
					bCreateNewNefForA = false;
					m_perf_log.push("Update intersection","Build intersections");
					m_Mesh_Intersection.increase_intersection_mesh(points_out,*it_patch_A,it_patch_B->second,Query_elem->id());
					m_perf_log.pop("Update intersection","Build intersections");
					++nbOfPositiveTests;
				}
			}
		}

		if(m_bPrintDebug)
		{
			std::cout << "    DEBUG: update intersection mesh" << std::endl;
			std::cout << " -> Positives / tests             : " << nbOfPositiveTests << " / " << nbOfTests
					  << " (" << 100.*nbOfPositiveTests/nbOfTests << "%)" << std::endl  << std::endl;
		}
	}

	void Intersection_Search::SetScalingFiles(const std::string& timing_data_file_base)
	{
		m_bPrintTimingData = true;
		m_bPrintIntersectionsPerPartData = true;
		m_timing_data_file_base = timing_data_file_base;
	}
}


