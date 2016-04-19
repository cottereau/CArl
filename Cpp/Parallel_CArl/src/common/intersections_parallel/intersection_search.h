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

	void BuildPatchIntersections_Brute()
	{
		// Code for the brute force intersection tests
		m_Intersection_Pairs.clear();

		// Dummy polyhedron containing the intersection
		ExactPolyhedron dummy_poly;

		// Intersection_Tools
		Intersection_Tools intersection_test;

		std::unordered_set<int>::iterator it_patch_A;
		std::unordered_set<int>::iterator it_patch_B;
		for(	it_patch_A =  m_Patch_Set_A.begin();
				it_patch_A != m_Patch_Set_A.end();
				++it_patch_A)
		{
			for(	it_patch_B =  m_Patch_Set_B.begin();
					it_patch_B != m_Patch_Set_B.end();
					++it_patch_B)
			{
				// For now I do nothing ...
			}
		}
	}

	void BuildPatchIntersections_Front()
	{

	}
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_ */
