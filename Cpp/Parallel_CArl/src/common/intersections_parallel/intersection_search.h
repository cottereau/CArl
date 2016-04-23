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
	FRONT
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
	const libMesh::Mesh&				   m_Mesh_A;
	const libMesh::Mesh&				   m_Mesh_B;
	const libMesh::Mesh&				   m_Mesh_Coupling;

	const libMesh::Parallel::Communicator& m_comm;

	Patch_construction					   m_Patch_Constructor_A;
	Patch_construction					   m_Patch_Constructor_B;

	std::unordered_set<unsigned int>	   m_Patch_Set_A;
	std::unordered_set<unsigned int>	   m_Patch_Set_B;

	Mesh_Intersection					   m_Mesh_Intersection;

	std::unordered_map<unsigned int,std::pair<unsigned int,unsigned int> > m_Intersection_Pairs;

	Intersection_Tools m_Intersection_test;

public:

	Intersection_Search(const libMesh::Mesh & mesh_A,
						const libMesh::Mesh & mesh_B,
						const libMesh::Mesh & mesh_Coupling,
						libMesh::Mesh & mesh_I) :
		m_Mesh_A { mesh_A },
		m_Mesh_B { mesh_B },
		m_Mesh_Coupling { mesh_Coupling },
		m_comm { m_Mesh_Coupling.comm() },
		m_Patch_Constructor_A { Patch_construction(m_Mesh_A)},
		m_Patch_Constructor_B { Patch_construction(m_Mesh_B)},
		m_Mesh_Intersection { Mesh_Intersection(mesh_I) }
	{
		// Reserve space for the unordered sets
		m_Patch_Set_A.reserve(mesh_A.n_elem());
		m_Patch_Set_B.reserve(mesh_B.n_elem());
		m_Intersection_Pairs.reserve(mesh_A.n_elem()*mesh_B.n_elem());
	};

	const libMesh::Mesh & mesh_A()
	{
		return m_Mesh_A;
	}

	const libMesh::Mesh & mesh_B()
	{
		return m_Mesh_B;
	}

	const libMesh::Mesh & mesh_Coupling()
	{
		return m_Mesh_Coupling;
	}

	void BuildCoupledPatches(const libMesh::Elem 	* Query_elem)
	{
		m_Patch_Constructor_A.BuildPatch(Query_elem,m_Patch_Set_A);
		m_Patch_Constructor_B.BuildPatch(Query_elem,m_Patch_Set_B);
	}

	void BuildPatchIntersections_Brute(const libMesh::Elem 	* Query_elem)
	{
		// TODO: For now it is more of a "find" than a build ...

		// Code for the brute force intersection tests
		m_Intersection_Pairs.clear();

		// Set containing all the intersection points
		std::set<Point_3> intersection_vertices;

		// Intersection_Tools
		m_Intersection_test.libmesh_set_coupling_nef_polyhedron(Query_elem);

		// Boolean: do they intersect?
		bool bDoIntersect = false;

		// Debug vars
		int nbOfTests = 0;
		int nbOfPositiveTests = 0;

		std::unordered_set<unsigned int>::iterator it_patch_A;
		std::unordered_set<unsigned int>::iterator it_patch_B;

		double total_volume = 0;
		for(	it_patch_A =  m_Patch_Set_A.begin();
				it_patch_A != m_Patch_Set_A.end();
				++it_patch_A)
		{
			const libMesh::Elem * elem_A = m_Mesh_A.elem(*it_patch_A);

			for(	it_patch_B =  m_Patch_Set_B.begin();
					it_patch_B != m_Patch_Set_B.end();
					++it_patch_B)
			{
				++nbOfTests;
				const libMesh::Elem * elem_B = m_Mesh_B.elem(*it_patch_B);

				intersection_vertices.clear();
				bDoIntersect = m_Intersection_test.libMesh_exact_intersection_inside_coupling(elem_A,elem_B,intersection_vertices);

				if(bDoIntersect && intersection_vertices.size() >= 4)
				{
					total_volume += m_Mesh_Intersection.get_intersection_volume(intersection_vertices);
					m_Intersection_Pairs[nbOfPositiveTests] = std::pair<int,int>(*it_patch_A,*it_patch_B);
					++nbOfPositiveTests;
				}
			}
		}
		std::cout << "    DEBUG: brute force search results" << std::endl;
		std::cout << " -> Positives / tests             : " << nbOfPositiveTests << " / " << nbOfTests
				  << " (" << 100.*nbOfPositiveTests/nbOfTests << "%)" << std::endl;
		std::cout << " -> Intersection / element volume : " << total_volume << " / " << Query_elem->volume() << std::endl << std::endl;
	}

	void FindFirstPair(		std::unordered_set<unsigned int> & 	Patch_guide,
							std::unordered_set<unsigned int> & 	Patch_probed,
							const libMesh::Mesh *  		Mesh_guide,
							const libMesh::Mesh *		Mesh_probed,
							std::pair<unsigned int,unsigned int> &		First_intersection)
	{
		// Set up locator for the first intersecting pair
		std::unique_ptr<libMesh::PointLocatorBase> Patch_guide_Locator = Mesh_guide->sub_point_locator();

		// Instruction needed to avoid the code from crashing if a query is
		// outside the mesh
		Patch_guide_Locator->enable_out_of_mesh_mode();

		// Find the first pair
		const libMesh::Elem * Patch_probed_first_elem = Mesh_probed->elem(*Patch_probed.begin());
		unsigned int Patch_probed_first_elem_id = Patch_probed_first_elem->id();

		// Set of intersecting terms
		std::set<const libMesh::Elem *> Intersecting_elems;
		m_Intersection_test.FindAllIntersection(Patch_probed_first_elem,Patch_guide_Locator,Intersecting_elems);

		// Search for one intersecting term inside the guide patch
		unsigned int Patch_guide_first_elem_id;
		bool bFoundFirstInter = false;
		for(std::set<const libMesh::Elem *>::iterator it_inter = Intersecting_elems.begin();
				it_inter != Intersecting_elems.end();
				++it_inter)
		{
			Patch_guide_first_elem_id = (*it_inter)->id();
			if(Patch_guide.find(Patch_guide_first_elem_id) != Patch_guide.end())
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

		// Code for the front intersection tests
		m_Intersection_Pairs.clear();

		// Set containing all the intersection points
		std::set<Point_3> intersection_vertices;

		// Intersection_Tools
		m_Intersection_test.libmesh_set_coupling_nef_polyhedron(Query_elem);

		// Boolean: do they intersect?
		bool bDoIntersect = false;

		// Debug vars
		int nbOfTests = 0;
		int nbOfPositiveTests = 0;

		double total_volume = 0;

		// --- Set up variables needed by Gander's algorithm

		// First, determinate which patch will be the guide, and which will be
		// the probed
		std::unordered_set<unsigned int> & Patch_guide	= m_Patch_Set_A;
		std::unordered_set<unsigned int> & Patch_probed 	= m_Patch_Set_B;
		const libMesh::Mesh *	  Mesh_guide	= &m_Mesh_A;
		const libMesh::Mesh *	  Mesh_probed	= &m_Mesh_B;

		// Exchange them if B is smaller
		if(m_Patch_Set_A.size() > m_Patch_Set_B.size())
		{
			std::cout << "    DEBUG: Using A as probe and B as guide" << std::endl;
			Patch_guide = m_Patch_Set_B;
			Patch_probed = m_Patch_Set_A;
			Mesh_guide	= &m_Mesh_B;
			Mesh_probed	= &m_Mesh_A;
		}
		else
		{
			std::cout << "    DEBUG: Using B as probe and A as guide" << std::endl;
		}

		// Find the first pair
		std::pair<unsigned int,unsigned int> First_intersection;
		FindFirstPair(Patch_guide,Patch_probed,Mesh_guide,Mesh_probed,First_intersection);

		std::cout << " > " << First_intersection.first << " " << First_intersection.second << std::endl << std::endl;

		// Deques containing the lists of the next tetrahedrons to be treated.
		std::deque<int> Patch_guide_Queue;
		std::deque<int> Patch_probed_Queue;

		// Marker vector, used to indicate if a tetrahedron from the guide was
		// already treated (=1) or not (=0).
		std::vector<int> Treated_From_guide(Patch_guide.size(),0);

		// Marker unordered set,used to indicate if a tetrahedron from the probe
		// was already treated (=1) or not (=0).
		std::unordered_set<unsigned int> Treated_From_probed(Patch_probed.size());

		std::cout << " > I do nothing!" << std::endl << std::endl;

		bDoIntersect = m_Intersection_test.libMesh_exact_intersection_inside_coupling(
				Mesh_guide->elem(First_intersection.first),Mesh_probed->elem(First_intersection.second),intersection_vertices);

		if(bDoIntersect && intersection_vertices.size() >= 4)
		{
			std::cout << m_Mesh_Intersection.get_intersection_volume(intersection_vertices) << std::endl;
		}
	}

	void BuildIntersections(SearchMethod search_type = BRUTE)
	{
		// Prepare iterators
		libMesh::Mesh::const_element_iterator it_coupl = m_Mesh_Coupling.elements_begin();
		libMesh::Mesh::const_element_iterator it_coupl_end = m_Mesh_Coupling.elements_end();

		for( ; it_coupl != it_coupl_end; ++it_coupl)
		{
			const libMesh::Elem * Query_elem = * it_coupl;

			BuildCoupledPatches(Query_elem);
			switch (search_type)
			{
				case BRUTE : 	BuildPatchIntersections_Brute(Query_elem);
								break;
				case FRONT : 	BuildPatchIntersections_Front(Query_elem);
								break;
			}
		}
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_ */
