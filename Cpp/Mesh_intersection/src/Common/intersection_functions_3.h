/*
 * intersection_functions_3.h
 *
 *  Created on: Jul 31, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef INTERSECTION_FUNCTIONS_3_H_
#define INTERSECTION_FUNCTIONS_3_H_

#include "common_header.h"
#include "CGAL_typedefs.h"
#include "triangular_mesh_3.h"

/**
 * 			This file contains the classes and declarations of functions needed
 * 		to do the intersection operations of two 2D meshes. This includes:
 *
 * 		- The intersection methods and the associated functions
 * 		- The intersection "visitors": classes depending of the "Variant" Boost
 * 		  class. They allow to apply different operations on a intersection's
 * 		  result depending on its type, i.e. whenever the intersection is a
 * 		  point, a section, a triangle or a polygon.
 *
 */

// -------- Visitors -------------------

/*
 * 			IntersectionPointsVisitor_3
 *
 * 			Intermediary visitor, used to keep the the points of the
 * 		intersections between two tetrahedrons.
 *
 * 		TODO: adapt the output operator overload (or build a new function) such
 * 		that it is well formated for FreeFem++
 *
 */
class IntersectionPointsVisitor_3
	: public boost::static_visitor<void>
{
public:
	// --- Members

	// Vector containing all the points
	std::vector<Point_3>	IntersectionData;

	// --- Constructors
	// Empty constructor.
	IntersectionPointsVisitor_3()
	{
	};

	// Constructor reserving a space "reserveSize" for the vectors. Notice that
	// this does not resize the vectors, but can optimize memory allocations
	IntersectionPointsVisitor_3(int reserveSize)
	{
		IntersectionData.reserve(reserveSize);
		IntersectionData.reserve(reserveSize);
	};

	// --- Visitor operators

	// Intersection is a point
	void operator()(const Point_3& p)
	{
		IntersectionData.push_back(p);
	};

	// Intersection is a segment
	void operator()(const Segment_3& s)
	{
		IntersectionData.push_back(s.source());
		IntersectionData.push_back(s.target());
	};

	// Intersection is a triangle
	void operator()(const Triangle_3& t)
	{
		IntersectionData.push_back(t.vertex(0));
		IntersectionData.push_back(t.vertex(1));
		IntersectionData.push_back(t.vertex(2));
	};

	// Intersection is a polygon
	void operator()(const std::vector<Point_3>& q)
	{

		for(std::vector<Point_3>::const_iterator itPoly = q.begin(); itPoly != q.end(); ++itPoly)
		{
			IntersectionData.push_back(*itPoly);
		}
	};

	// Return information from the visitor (size of the vectors)
	std::string printStatics()
	{
		return 	std::to_string(IntersectionData.size());
	}

	void add(const Point_3& p)
	{
		IntersectionData.push_back(p);
	}
};

//	Friend function - overloaded operator<< for 'IntersectionPointsVisitor_3'
std::ostream& operator<<(std::ostream& out, IntersectionPointsVisitor_3& input);

// -------- Functions ---------------------------------------------------------

/*
 * 			Brute force intersection finding algorithm for 2D triangulations.
 *
 * 			Runs over the whole meshes dtA and dtB, looking for the
 * 		intersections. Complexity is O(n*m), where n and m are the number of
 * 		triangles in each triangulation.
 *
 */
void BuildMeshIntersections_Brute(
		const Triangular_Mesh_3& dtA,
		const Triangular_Mesh_3& dtB,
		std::vector<Polyhedron>& 	output
		);

/*
 * 			Optimized intersection finding algorithm for 2D triangulations.
 *
 * 			The optimized algorithm runs the external loop on all the
 * 		intersecting triangles of dtB, using a "frontier algorithm" to find the
 * 		intersections with the triangulation dtA. The internal loop runs only
 * 		over the triangles of dtA near a given triangle of dtB. Its complexity
 * 		is, then, O(m) at most cases, where m is the number of triangles inside
 * 		dtB, and O(n*m) ~ O(m^2) in the worst case (all triangles intersect each
 * 		other).
 *
 * 			If the starting triangles are given as an argument (initA and
 * 		initB),	he algorithm itself does not do any search operations, else it
 * 		only does one for the first pair of intersecting triangles.
 *
 * 			Algorithm based on the article
 *
 * 		Gander, M. J. et al., "An Algorithm for Non-Matching Grid Projections
 * 			with Linear Complexity". Domain Decomposition Methods in Science and
 * 			Engineering XVIII, p. 185-192, 2009, Springer.
 *
 *
 * 		TODO: implement version with initA and initB.
 *
 */
void BuildMeshIntersections(
		const Triangular_Mesh_3& 	dtA,
		const Triangular_Mesh_3& 	dtB,
		std::vector<Polyhedron>& 	output
		);

/*
 * 			Find if there is an intersection between two triangles.
 *
 * 			If such an intersection exists, determinate it. The vector
 * 		candidatesAUpdate is needed for BuildMeshIntersections' algorithm
 * 		(detailed in the source file).
 *
 */
void IntersectTetrahedrons(
		const Triangular_Mesh_3& 	dtB,
		const Cell_handle_3& 		workingTetrahedronB,
		const Triangular_Mesh_3& 	dtA,
		const Cell_handle_3& 		workingTetrahedronA,
		std::vector<int>& 			candidatesAUpdate,
		bool& 						queryIntersect,
		bool&						inexactQueryIntersect,
		Polyhedron& 				output
		);

/*
 * 			Find an intersecting pair of triangles between two triangulations.
 *
 */
void FindFirstPair(
		const Triangular_Mesh_3& 	dtA,
		const Triangular_Mesh_3& 	dtB,
		Cell_handle_3&				FirstA,
		Cell_handle_3&				FirstB
		);

#endif /* INTERSECTION_FUNCTIONS_H_ */
