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

// Can be used to convert polyhedron from exact to inexact and vice-versa
template <class Polyhedron_input,
class Polyhedron_output>
struct Copy_polyhedron_to
        : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
        Copy_polyhedron_to(const Polyhedron_input& in_poly)
                : in_poly(in_poly) {}

        void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
        {
                typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
//                typedef typename Polyhedron_input::HalfedgeDS Input_HDS;

                CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

                typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
                typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
                typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

                builder.begin_surface(in_poly.size_of_vertices(),
                        in_poly.size_of_facets(),
                        in_poly.size_of_halfedges());

                for(Vertex_const_iterator
                        vi = in_poly.vertices_begin(), end = in_poly.vertices_end();
                        vi != end ; ++vi)
                {
                        typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
                                ::CGAL::to_double( vi->point().y()),
                                ::CGAL::to_double( vi->point().z()));
                        builder.add_vertex(p);
                }

                typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
                Index index( in_poly.vertices_begin(), in_poly.vertices_end());

                for(Facet_const_iterator
                        fi = in_poly.facets_begin(), end = in_poly.facets_end();
                        fi != end; ++fi)
                {
                        HFCC hc = fi->facet_begin();
                        HFCC hc_end = hc;
                        builder.begin_facet ();
                        do {
                                builder.add_vertex_to_facet(index[hc->vertex()]);
                                ++hc;
                        } while( hc != hc_end);
                        builder.end_facet();
                }
                builder.end_surface();
        } // end operator()(..)
private:
        const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_B, class Poly_A>
void poly_copy(Poly_B& poly_b, const Poly_A& poly_a)
{
        poly_b.clear();
        Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
        poly_b.delegate(modifier);
}


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
	Nef_Polyhedron m_nef_C;

	std::vector<ExactPoint_3> m_exact_points_A;
	std::vector<ExactPoint_3> m_exact_points_B;
	std::vector<ExactPoint_3> m_exact_points_C;

	ExactPolyhedron m_dummyPoly;

	CGAL::Convex_hull_traits_3<ExactKernel> ExactHullTraits;

	Kernel_to_ExactKernel ConvertInexactToExact;
	ExactKernel_to_Kernel ConvertExactToInexact;

public:
	Intersection_Tools(const libMesh::Elem * elem_C)
	{
		m_exact_points_A.resize(8);
		m_exact_points_B.resize(8);
		m_exact_points_C.resize(8);

		m_dummyPoly.reserve(8,18,12);

		libmesh_set_coupling_nef_polyhedron(elem_C);
	};

	Intersection_Tools()
	{
		m_exact_points_A.resize(8);
		m_exact_points_B.resize(8);

		m_dummyPoly.reserve(8,18,12);

		m_nef_C.clear(Nef_Polyhedron::EMPTY);
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
									ExactPolyhedron & poly_out)
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
	 * 		Set the coupling element Nef polyhedron
	 */
	void libmesh_set_coupling_nef_polyhedron(const libMesh::Elem * elem_C)
	{
		unsigned int n_nodes_C = elem_C->n_nodes();

		// First, convert the element to a CGAL exact point vector
		convert_elem_to_exact_points(elem_C,m_exact_points_C);

		// Then convert to Nef
		std::vector<ExactPoint_3>::const_iterator exact_points_C_begin = m_exact_points_C.begin();
		convert_exact_points_to_Nef(	exact_points_C_begin,
										exact_points_C_begin + n_nodes_C,
										m_nef_C);
	}

	/*
	 * 		Test if two elements intersect inside the coupling region
	 */
	bool libMesh_exact_do_intersect_inside_coupling(const libMesh::Elem * elem_A,
													const libMesh::Elem * elem_B)
	{
		homemade_assert_msg(!m_nef_C.is_empty(), "Coupling restriction element was not set yet!\n");

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

			m_nef_I = m_nef_A*m_nef_B*m_nef_C;
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
	 * 		Build two elements intersection inside the coupling region
	 */
	bool libMesh_exact_intersection_inside_coupling(const libMesh::Elem * elem_A,
													const libMesh::Elem * elem_B,
													std::set<Point_3> & points_out)
	{
		// Test the intersection and build the Nef polyhedron (if true)
		bool bElemIntersect = libMesh_exact_do_intersect_inside_coupling(elem_A,elem_B);

		if(bElemIntersect && m_nef_I.number_of_volumes() > 1)
		{
			for(Nef_Polyhedron::Vertex_const_iterator it_vertex = m_nef_I.vertices_begin();
					it_vertex != m_nef_I.vertices_end();
					++it_vertex)
			{
				points_out.insert(ConvertExactToInexact(it_vertex->point()));
			}
		}
		else
		{
			bElemIntersect = false;
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
