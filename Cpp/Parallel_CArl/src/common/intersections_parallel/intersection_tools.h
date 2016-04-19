/*
 * intersection_tools.h
 *
 *  Created on: Apr 17, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_TOOLS_H_
#define COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_TOOLS_H_

#include "common_header.h"
#include "common_header_libmesh.h"
#include "mesh_tables.h"

#include "CGAL_typedefs.h"

namespace carl
{

/*
 * 		Intersection_tools
 *
 * 			This class contains several tools that help to convert libMesh's
 * 		data structures and vice-versa. A class is used instead of a collection
 * 		of functions to reduce the amount of calls to constructors of CGAL data
 * 		exact data structures.
 *
 *		TODO : 	optimize the code by converting this class to one that saves
 *				all the exact geometry of a mesh for reuse.
 */

class	Intersection_Tools
{
protected:
	Nef_Polyhedron m_nef_A;
	Nef_Polyhedron m_nef_B;
	Nef_Polyhedron m_nef_I;

	std::vector<ExactPoint_3> m_exact_points_A;
	std::vector<ExactPoint_3> m_exact_points_B;

	ExactPolyhedron m_dummyPoly;

	CGAL::Convex_hull_traits_3<ExactKernel> ExactHullTraits;

	Kernel_to_ExactKernel ConvertInexactToExact;
	ExactKernel_to_Kernel ConvertExactToInexact;

public:
	Intersection_Tools()
	{
		m_exact_points_A.resize(8);
		m_exact_points_B.resize(8);

		m_dummyPoly.reserve(8,18,12);
	};

	/*
	 * 		Test if two elements intersect
	 */
	bool libMesh_exact_do_intersect(const libMesh::Elem * elem_A,
									const libMesh::Elem * elem_B)
	{
		bool bElemIntersect = false;

		unsigned int n_nodes_A = elem_A->n_nodes();
		unsigned int n_nodes_B = elem_B->n_nodes();

		// First, convert both elements to CGAL exact point vectors
		convert_elem_to_exact_points(elem_A,m_exact_points_A);
		convert_elem_to_exact_points(elem_B,m_exact_points_B);

		// Fast check: bounding boxes
		std::vector<ExactPoint_3>::const_iterator exact_points_A_begin = m_exact_points_A.begin();
		std::vector<ExactPoint_3>::const_iterator exact_points_B_begin = m_exact_points_B.begin();

		Bbox_3 exact_bbox_A = CGAL::bbox_3(exact_points_A_begin,exact_points_A_begin + n_nodes_A);
		Bbox_3 exact_bbox_B = CGAL::bbox_3(exact_points_B_begin,exact_points_B_begin + n_nodes_B);

		bool bBboxIntersect = CGAL::do_intersect(exact_bbox_A,exact_bbox_B);

		// If they do intersect, we'll need to do a better test
		if(bBboxIntersect)
		{
			// Convert to Nef polyhedron ... ouch
			convert_exact_points_to_Nef(	exact_points_A_begin,
											exact_points_A_begin + n_nodes_A,
											m_nef_A);

			convert_exact_points_to_Nef(	exact_points_B_begin,
											exact_points_B_begin + n_nodes_B,
											m_nef_B);

			m_nef_I = m_nef_A*m_nef_B;
			if(!m_nef_I.is_empty())
			{
				bElemIntersect = true;
			}
			else
			{
				bElemIntersect = false;
			}
		}
		else
		{
			bElemIntersect = false;
		}

		return bElemIntersect;
	}

	/*
	 * 		Build two elements intersection
	 */
	bool libMesh_exact_intersection(const libMesh::Elem * elem_A,
									const libMesh::Elem * elem_B,
									ExactPolyhedron poly_out)
	{
		// Test the intersection and build the Nef polyhedron (if true)
		bool bElemIntersect = libMesh_exact_do_intersect(elem_A,elem_B);

		if(bElemIntersect)
		{
			Nef_Polyhedron::Volume_const_iterator itVol = ++m_nef_I.volumes_begin();
			m_nef_I.convert_inner_shell_to_polyhedron(itVol->shells_begin(), poly_out);
		}

		return bElemIntersect;
	}

	/*
	 * 		Convert an libMesh element to list of exact points
	 */
	void convert_elem_to_exact_points(	const libMesh::Elem       *	elem_input,
										std::vector<ExactPoint_3> &	points_output)
	{
		libMesh::Point dummyPoint;
		for(unsigned int iii = 0; iii < elem_input->n_nodes(); ++iii)
		{
			dummyPoint = elem_input->point(iii);
			points_output[iii] = ExactPoint_3(	dummyPoint(0),
												dummyPoint(1),
												dummyPoint(2));
		}
	}

	/*
	 * 		Convert an exact point set to a Nef polyhedron
	 */
	void convert_exact_points_to_Nef(std::vector<ExactPoint_3>::const_iterator it_begin,
									 std::vector<ExactPoint_3>::const_iterator it_end,
									 Nef_Polyhedron & nef_out)
	{
		CGAL::convex_hull_3(it_begin,it_end,m_dummyPoly,ExactHullTraits);
		nef_out = Nef_Polyhedron(m_dummyPoly);
	}
};
}



#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_TOOLS_H_ */
