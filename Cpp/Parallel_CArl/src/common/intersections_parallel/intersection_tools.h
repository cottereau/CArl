/*
 * intersection_tools.h
 *
 *  Created on: Apr 17, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_TOOLS_H_
#define COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_TOOLS_H_

#include "carl_headers.h"
#include "mesh_tables.h"

namespace carl
{

// ************************
// Intersection_Tools class
// ************************

/** \brief Class with a series of methods to find the intersections between 
 *  libMesh's elements, using CGAL internally.
 */

class	Intersection_Tools
{
protected:

	// Nef polyhedrons
	Nef_Polyhedron m_nef_A;		///< Nef polyhedrons for the elements.
	Nef_Polyhedron m_nef_B;
	Nef_Polyhedron m_nef_I;
	Nef_Polyhedron m_nef_C;

	/**  \brief Nef polyhedron used to save the intersection between m_nef_A and
	 * 	 m_nef_C, reducing the number of intersection operations.
	 */
	Nef_Polyhedron m_nef_I_AC;

	// Exact point vectors
	std::vector<ExactPoint_3> m_exact_points_A; ///< CGAL Exact point vectors.
	std::vector<ExactPoint_3> m_exact_points_B;
	std::vector<ExactPoint_3> m_exact_points_C;

	ExactPolyhedron m_dummyPoly; ///< Intersection polyhedron.

	CGAL::Convex_hull_traits_3<ExactKernel> ExactHullTraits; ///< CGAL ExactKernel's convex hull traits.

	// CGAL Kernel converters
	Kernel_to_ExactKernel ConvertInexactToExact;	///< Convert inexact CGAL constructs to exact constructs.
	ExactKernel_to_Kernel ConvertExactToInexact;	///< Convert exact CGAL constructs to inexact constructs.

	// Dummy exact tetrahedron and triangle
	ExactTetrahedron 	m_test_tetra;
	ExactTriangle_3 	m_test_triangle;

	double m_Min_Inter_Volume;	///<  Minimal volume cutoff for the intersections.

	// Perflog and debug variables
	bool MASTER_bPerfLog_intersection_tools; ///< Do performance log? *Default:* false.
	libMesh::PerfLog m_perf_log;			 ///< libMesh::PerfLog object.

	// Data structures
	/*
	 * 		CGAL's "do_intersect" methods only works with tetrahedrons and triangles, at most.
	 *  The vectors below map an TET or HEX element to its decomposition with tetrahedrons and
	 *  triangles (which is trivial for the TET elements).
	 */
	std::vector<std::vector<unsigned int> > m_TET_tetrahedrons;	///< Vertex indices for a TET* element tetrahedron
	std::vector<std::vector<unsigned int> > m_TET_triangles;	///< Vertex indices for a TET* element surface triangles
	std::vector<std::vector<unsigned int> > m_TET_edges;		///< Vertex indices for a TET* element edges (**not used!**)

	std::vector<std::vector<unsigned int> > m_HEX_tetrahedrons;	///< Vertex indices for a HEX* element tetrahedron decomposition
	std::vector<std::vector<unsigned int> > m_HEX_triangles; 	///< Vertex indices for a HEX* element surface triangles decomposition
	std::vector<std::vector<unsigned int> > m_HEX_edges; 		///< Vertex indices for a HEX* element edges (**not used!**)

	std::vector<std::vector<unsigned int> > * m_elem_C_tetrahedrons;	///< Pointer for the appropriate tetrahedron vertex index list (**useless!**)
	std::vector<std::vector<unsigned int> > * m_elem_C_triangles;		///< Pointer for the appropriate triangle vertex index list (**useless!**)

	// PROTECTED methods
	/** \brief Build the tetrahedron and triangle decomposition mappings.
	 * 
	 * 		CGAL's "do_intersect" methods are not compatible with generic polyhedrons - at most,
	 *  we can determinate if a triangle intersects a tetrahedron. To determinate if two elements 
	 *  A and B intersect, we must then decompose them into these geometrical objects, and test 
	 *  the intersections pairwise. The Intersection_Tools::set_element_indexes method fills the
	 *  vectors "m_XYZ_tetrahedrons" and "m_XYZ_triangles" with a mapping between an TET or HEX
	 *  element and these decompositions. 
	 */
	void set_element_indexes();

	/** \brief Determinate if two elements intersect, using their points and decompositions.
	 * 
	 * 		CGAL's "do_intersect" methods are not compatible with generic polyhedrons - at most,
	 *  we can determinate if a triangle intersects a tetrahedron. To determinate if two elements 
	 *  A and B intersect, we must then decompose them into these geometrical objects, and test 
	 *  the intersections pairwise. This method takes the elements' points and decompositions into 
	 *  tetrahedrons and triangles and test each tetrahedron / triangle pair for intersections.
	 */
	bool elements_do_intersect(
			std::vector<ExactPoint_3> & elem_C_points,
			std::vector<std::vector<unsigned int> > & elem_C_tetras,
			std::vector<std::vector<unsigned int> > & elem_C_triangles,
			std::vector<ExactPoint_3> & elem_D_points,
			std::vector<std::vector<unsigned int> > & elem_D_tetras,
			std::vector<std::vector<unsigned int> > & elem_D_triangles);

public:

	// Constructors
	/** \brief Constructor with a pre-defined coupling element. It preallocates most objects
	 *	and build the decomposition mappings.
	 */
	Intersection_Tools(const libMesh::Elem * elem_C, double Min_Inter_Volume = 1E-21, bool bDoPerf_log = false) :
		m_Min_Inter_Volume { Min_Inter_Volume },
		MASTER_bPerfLog_intersection_tools {bDoPerf_log},
		m_perf_log { libMesh::PerfLog("Intersection tools", MASTER_bPerfLog_intersection_tools) }
	{
		m_exact_points_A.resize(8);
		m_exact_points_B.resize(8);
		m_exact_points_C.resize(8);

		m_dummyPoly.reserve(8,18,12);

		this->libmesh_set_coupling_nef_polyhedron(elem_C);

		this->set_element_indexes();

		m_elem_C_tetrahedrons = NULL;
		m_elem_C_triangles = NULL;
	};

	/** \brief Constructor preallocating most objects and building the decomposition mappings.
	 */
	Intersection_Tools(double Min_Inter_Volume = 1E-21, bool bDoPerf_log = false) :
		m_Min_Inter_Volume { Min_Inter_Volume },
		MASTER_bPerfLog_intersection_tools {bDoPerf_log},
		m_perf_log { libMesh::PerfLog("Intersection tools", MASTER_bPerfLog_intersection_tools) }
	{
		m_exact_points_A.resize(8);
		m_exact_points_B.resize(8);
		m_exact_points_C.resize(8);

		m_dummyPoly.reserve(8,18,12);

		m_nef_C.clear(Nef_Polyhedron::EMPTY);

		this->set_element_indexes();

		m_elem_C_tetrahedrons = NULL;
		m_elem_C_triangles = NULL;
	};

	// PUBLIC methods
	/** \brief Find an element from a mesh intersecting Query_elem, using
	 *	the mesh's libMesh::PointLocatorBase.
	 * 			
	 * 		By default, a test is done to be sure that the query element does
	 * 	indeed intersect the tested mesh. The test can be bypassed using the flag
	 *	bGuaranteeQueryIsInMesh.
	 */
	const libMesh::Elem * FindFirstIntersection(	const libMesh::Elem * Query_elem,
								std::unique_ptr<libMesh::PointLocatorBase> & point_locator,
								bool				bGuaranteeQueryIsInMesh = false);

	/** \brief Find all elements from the mesh intersecting Query_elem, using
	 *	the mesh's libMesh::PointLocatorBase.
	 */
	bool FindAllIntersection(	const libMesh::Elem * Query_elem,
								std::unique_ptr<libMesh::PointLocatorBase> & point_locator,
								std::set<unsigned int>	&	Intersecting_elems);

	/** \brief Determinate if two elements intersect.
	 *
	 *		Serves as a frontend for the private Intersection_Tools::elements_do_intersect method.
	 */
	bool libMesh_exact_do_intersect(const libMesh::Elem * elem_A,
									const libMesh::Elem * elem_B);

	/** \brief Build the intersection between two elements.
	 */
	bool libMesh_exact_intersection(const libMesh::Elem * elem_A,
									const libMesh::Elem * elem_B,
									std::set<Point_3> & points_out,
									bool bCreateNewNefForA = true,
									bool bConvertPoints = true,
									bool bTestNeeded = true);

	/** \brief Determinate if two elements intersect inside the coupling region.
	 * 
	 *		This method uses CGAL's Nef polyhedrons. The boolean flag "bCreateNewNefForA" 
     * 	is used to avoid re-building the elem_A's Nef polyhedron. 
	 */
	bool libMesh_exact_do_intersect_inside_coupling(const libMesh::Elem * elem_A,
													const libMesh::Elem * elem_B,
													bool bCreateNewNefForA = true);

	/** \brief Build the intersection between two elements intersect inside the coupling region.
	 * 
	 *		This method uses CGAL's Nef polyhedrons. The boolean flag "bCreateNewNefForA" 
     * 	is used to avoid re-building the elem_A's Nef polyhedron. If the boolean flag "bConvertPoints"
	 *	is set to true, the elements points are converted to CGAL's exact points. When the boolean
	 *  flag "bTestNeeded" is set to true, a preliminary intersection test is done.
	 */
	bool libMesh_exact_intersection_inside_coupling(const libMesh::Elem * elem_A,
													const libMesh::Elem * elem_B,
													std::set<libMesh::Point> & points_out,
													bool bCreateNewNefForA = true,
													bool bConvertPoints = true,
													bool bTestNeeded = false);

	/** \brief Set the Nef polyhedron for the coupling mesh element 
	 */
	void libmesh_set_coupling_nef_polyhedron(const libMesh::Elem * elem_C);

	/** \brief Convert an libMesh element to a list of CGAL exact points
	 */
	void convert_elem_to_exact_points(	const libMesh::Elem       *	elem_input,
										std::vector<ExactPoint_3> &	points_output);

	/** \brief Convert an CGAL exact point set to a CGAL Nef polyhedron
	 */
	void convert_exact_points_to_Nef(std::vector<ExactPoint_3>::const_iterator it_begin,
									 std::vector<ExactPoint_3>::const_iterator it_end,
									 Nef_Polyhedron & nef_out);
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_TOOLS_H_ */
