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

	std::unordered_set<int>				   m_Patch_Set_A;
	std::unordered_set<int>				   m_Patch_Set_B;

	Mesh_Intersection					   m_Mesh_Intersection;

	std::unordered_map<int,std::pair<int,int> > m_Intersection_Pairs;

public:

	enum SearchMethod
	{
		BRUTE,
		FRONT
	};

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
		Intersection_Tools intersection_test(Query_elem);

		// Boolean: do they intersect?
		bool bDoIntersect = false;

		// Debug vars
		int nbOfTests = 0;
		int nbOfPositiveTests = 0;

		std::unordered_set<int>::iterator it_patch_A;
		std::unordered_set<int>::iterator it_patch_B;

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
				bDoIntersect = intersection_test.libMesh_exact_intersection_inside_coupling(elem_A,elem_B,intersection_vertices);

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

	void BuildPatchIntersections_Front()
	{

	}

	void BuildIntersections(SearchMethod = BRUTE)
	{
		// Prepare iterators
		libMesh::Mesh::const_element_iterator it_coupl = m_Mesh_Coupling.elements_begin();
		libMesh::Mesh::const_element_iterator it_coupl_end = m_Mesh_Coupling.elements_end();

		for( ; it_coupl != it_coupl_end; ++it_coupl)
		{
			const libMesh::Elem * Query_elem = * it_coupl;
			BuildCoupledPatches(Query_elem);
			BuildPatchIntersections_Brute(Query_elem);
		}
	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_ */
