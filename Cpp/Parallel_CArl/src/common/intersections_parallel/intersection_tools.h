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

/*
 * 		Intersection_tools
 *
 * 			This class contains several tools that help to convert libMesh's
 * 		data structures and vice-versa. A class is used instead of a collection
 * 		of functions to reduce the amount of calls to constructors of CGAL data
 * 		exact data structures.
 *
 */

class	Intersection_Tools
{
protected:

	// Nef polyhedrons
	Nef_Polyhedron m_nef_A;
	Nef_Polyhedron m_nef_B;
	Nef_Polyhedron m_nef_I;
	Nef_Polyhedron m_nef_C;

	/*
	 * 		Nef polyhedron representing the intersection between m_nef_A and
	 * 	 m_nef_C, used to reduce the number of intersection operations
	 */
	Nef_Polyhedron m_nef_I_AC;

	// Exact point vectors
	std::vector<ExactPoint_3> m_exact_points_A;
	std::vector<ExactPoint_3> m_exact_points_B;
	std::vector<ExactPoint_3> m_exact_points_C;

	// Intersection polyhedron
	ExactPolyhedron m_dummyPoly;

	// ExactKernel's convex hull traits
	CGAL::Convex_hull_traits_3<ExactKernel> ExactHullTraits;

	// CGAL Kernel converters
	Kernel_to_ExactKernel ConvertInexactToExact;
	ExactKernel_to_Kernel ConvertExactToInexact;

	// Dummy exact tetrahedron and triangle
	ExactTetrahedron 	m_test_tetra;
	ExactTriangle_3 	m_test_triangle;

	// Volume cutoff for the intersections
	double m_Min_Inter_Volume;

	// Perflog and debug variables
	bool MASTER_bPerfLog_intersection_tools;
	libMesh::PerfLog m_perf_log;

	// Data structures
	/*
	 * 		Each vector contains the indexes of the geometric decompositions of
	 * 	either a tetrahedron or an hexahedron into tetrahedrons, triangles and
	 * 	edges. These decompositions are needed because CGAL do_intersect methods
	 * 	only accept the decomposed geometry entities.
	 */
	std::vector<std::vector<unsigned int> > m_tetra_tetrahedrons;
	std::vector<std::vector<unsigned int> > m_tetra_triangles;
	std::vector<std::vector<unsigned int> > m_tetra_edges;

	std::vector<std::vector<unsigned int> > m_hex_tetrahedrons;
	std::vector<std::vector<unsigned int> > m_hex_triangles;
	std::vector<std::vector<unsigned int> > m_hex_edges;

	/*
	 * 		Vector pointers used represent the data structure above
	 */
	std::vector<std::vector<unsigned int> > * m_elem_C_tetras;
	std::vector<std::vector<unsigned int> > * m_elem_C_triangles;

	// PROTECTED methods
	/*
	 * 		Method defining these data structures
	 */
	void set_element_indexes();

	/*
	 * 		Method used to determinate if two elements intersect
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
	Intersection_Tools(const libMesh::Elem * elem_C, double Min_Inter_Volume = 1E-21, bool bDoPerf_log = false) :
		m_Min_Inter_Volume { Min_Inter_Volume },
		MASTER_bPerfLog_intersection_tools {bDoPerf_log},
		m_perf_log { libMesh::PerfLog("Intersection tools", MASTER_bPerfLog_intersection_tools) }
	{
		m_exact_points_A.resize(8);
		m_exact_points_B.resize(8);
		m_exact_points_C.resize(8);

		m_dummyPoly.reserve(8,18,12);

		libmesh_set_coupling_nef_polyhedron(elem_C);

		set_element_indexes();

		m_elem_C_tetras = NULL;
		m_elem_C_triangles = NULL;
	};

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

		set_element_indexes();

		m_elem_C_tetras = NULL;
		m_elem_C_triangles = NULL;
	};

	// PUBLIC methods
	/*
	 * 			Find a element from the mesh intersecting the query element.
	 * 		Does so while doing a test to be sure that the query element does
	 * 		indeed intersect the tested mesh. The test can be bypassed using a
	 * 		boolean.
	 *
	 */
	const libMesh::Elem * FindFirstIntersection(	const libMesh::Elem * Query_elem,
								std::unique_ptr<libMesh::PointLocatorBase> & point_locator,
								bool				bGuaranteeQueryIsInMesh = false);

	/*
	 * 		Find all elements from the mesh intersecting the query element.
	 * 	Does so while doing a test to be sure that the query element does
	 * 	indeed intersect the tested mesh.
	 */
	bool FindAllIntersection(	const libMesh::Elem * Query_elem,
								std::unique_ptr<libMesh::PointLocatorBase> & point_locator,
								std::set<unsigned int>	&	Intersecting_elems);

	/*
	 * 		Test if two elements intersect
	 */
	bool libMesh_exact_do_intersect(const libMesh::Elem * elem_A,
									const libMesh::Elem * elem_B);
	/*
	 * 		Build two elements intersection
	 */
	bool libMesh_exact_intersection(const libMesh::Elem * elem_A,
									const libMesh::Elem * elem_B,
									std::set<Point_3> & points_out,
									bool bCreateNewNefForA = true,
									bool bConvertPoints = true,
									bool bTestNeeded = true);

	/*
	 * 		Test if two elements intersect inside the coupling region
	 */
	bool libMesh_exact_do_intersect_inside_coupling(const libMesh::Elem * elem_A,
													const libMesh::Elem * elem_B,
													bool bCreateNewNefForA = true);

	/*
	 * 		Build two elements intersection inside the coupling region
	 */
	bool libMesh_exact_intersection_inside_coupling(const libMesh::Elem * elem_A,
													const libMesh::Elem * elem_B,
													std::set<libMesh::Point> & points_out,
													bool bCreateNewNefForA = true,
													bool bConvertPoints = true,
													bool bTestNeeded = true);
	/*
	 * 		Set the coupling element Nef polyhedron
	 */
	void libmesh_set_coupling_nef_polyhedron(const libMesh::Elem * elem_C);

	/*
	 * 		Convert an libMesh element to list of exact points
	 */
	void convert_elem_to_exact_points(	const libMesh::Elem       *	elem_input,
										std::vector<ExactPoint_3> &	points_output);

	/*
	 * 		Convert an exact point set to a Nef polyhedron
	 */
	void convert_exact_points_to_Nef(std::vector<ExactPoint_3>::const_iterator it_begin,
									 std::vector<ExactPoint_3>::const_iterator it_end,
									 Nef_Polyhedron & nef_out);
};
}



#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_TOOLS_H_ */
