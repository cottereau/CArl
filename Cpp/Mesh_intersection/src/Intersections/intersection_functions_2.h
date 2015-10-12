/*
 * intersection_functions_2.h
 *
 *  Created on: Jun 12, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef INTERSECTION_FUNCTIONS_2_H_
#define INTERSECTION_FUNCTIONS_2_H_

#include "common_header.h"
#include "CGAL_typedefs.h"
#include "triangular_mesh_2.h"
#include "triangular_mesh_intersection_2.h"

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
 * 			IntersectionTallyVisitor
 *
 * 			This visitor simply does a tally of all intersections, divided by
 * 		their type. Output is done with the 'operator<<' overload, and it
 * 		returns	a table with the tallies, formatted for easy (human)
 * 		readability.
 *
 */
class IntersectionTallyVisitor
	: public boost::static_visitor<void>
{
public:
	// --- Members
	// Vector containing the point intersections.
	std::vector<Point_2>				point_intersections;

	// Vector containing the segment intersections.
	std::vector<Segment_2>				seg_intersections;

	// Vector containing the triangle intersections.
	std::vector<Triangle_2>				tri_intersections;

	// Vector containing the polygon intersections.
	std::vector<std::vector<Point_2> >	poly_intersections;

	// --- Constructors
	// Empty constructor.
	IntersectionTallyVisitor()
	{

	};

	// Constructor reserving a space "reserveSize" for the vectors. Notice that
	// this does not resize the vectors, but can optimize memory allocations
	IntersectionTallyVisitor(int reserveSize)
	{
		point_intersections.reserve(reserveSize);
		seg_intersections.reserve(reserveSize);
		tri_intersections.reserve(reserveSize);
		poly_intersections.reserve(reserveSize);
	};

	// --- Visitor operators

	// Push back a point intersection.
	void operator()(const Point_2& p)
	{
		point_intersections.push_back(p);
	};

	// Push back a segment intersection.
	void operator()(const Segment_2& s)
	{
		seg_intersections.push_back(s);
	};

	// Push back a triangle intersection.
	void operator()(const Triangle_2& t)
	{
		tri_intersections.push_back(t);
	};

	// Push back a polygon intersection.
	void operator()(const std::vector<Point_2>& q)
	{
		poly_intersections.push_back(q);
	};

	// Return information from the visitor (size of the vectors)
	std::string printStatics()
	{
		return 	std::to_string(point_intersections.size()) + " " +
				std::to_string(seg_intersections.size()) + " " +
				std::to_string(tri_intersections.size()) + " " +
				std::to_string(poly_intersections.size());
	};

	// Finalize does nothing for this visitor
	void Finalize()
	{

	};
};

//	Friend function - overloaded operator<< for 'IntersectionTallyVisitor
std::ostream& operator<<(std::ostream& out, IntersectionTallyVisitor& input);

/*
 * 			PolyIntersectionVisitor
 *
 * 			This visitor keeps only the informations from the 'triangle' and
 * 		'polygon' intersection results. Output is done with the 'operator<<'
 * 		overload, and is formatted for easy (machine) reading.
 *
 * 		TODO: implement it using CGAL's Nef polyhedrons (for generality and to
 * 			  avoid using the inexact intersection operations)
 */
class PolyIntersectionVisitor
	: public boost::static_visitor<void>
{
public:
	// --- Members

	// Vector containing the triangle intersections.
	std::vector<Triangle_2>				tri_intersections;

	// Vector containing the polygon intersections.
	std::vector<std::vector<Point_2> >	poly_intersections;

	// --- Constructors
	// Empty constructor.
	PolyIntersectionVisitor()
	{

	};

	// Constructor reserving a space "reserveSize" for the vectors. Notice that
	// this does not resize the vectors, but can optimize memory allocations
	PolyIntersectionVisitor(int reserveSize)
	{
		tri_intersections.reserve(reserveSize);
		poly_intersections.reserve(reserveSize);
	};

	// --- Visitor operators

	// Intersection is a point -> do nothing.
	void operator()(const Point_2& p)
	{
	};

	// Intersection is a segment -> do nothing.
	void operator()(const Segment_2& s)
	{
	};

	// Push back a triangle intersection.
	void operator()(const Triangle_2& t)
	{
		tri_intersections.push_back(t);
	};

	// Push back a polygon intersection.
	void operator()(const std::vector<Point_2>& q)
	{
		poly_intersections.push_back(q);
	};

	// Finalize does nothing for this visitor
	void Finalize()
	{

	};

	// Return information from the visitor (size of the vectors)
	std::string printStatics()
	{
		return 	std::to_string(tri_intersections.size()) + " " +
				std::to_string(poly_intersections.size());
	}
};

//	Friend function - overloaded operator<< for 'PolyIntersectionVisitor'
std::ostream& operator<<(std::ostream& out, PolyIntersectionVisitor& input);

/*
 * 			TriangulationIntersectionVisitor
 *
 * 			This visitor takes the informations from the 'triangle' and
 * 		'polygon' intersection results, and convert them to a Triangulation_2.
 * 		NOTE: the triangulation itself is NOT valid for the CGAL standards.
 *
 */
class TriangulationIntersectionVisitor
	: public boost::static_visitor<void>
{
private:
	// --- Constructor

	// Empty constructor - set as private
	TriangulationIntersectionVisitor()
	{
		CharacteristicLength = -1;
		CharacteristicArea   = -1;
	};
public:
	// --- Members

	Intersection_Mesh_2	IntersectionTDS_2;
	double				CharacteristicLength;
	double				CharacteristicArea;

	// --- Constructors
	// Constructor with an characteristic length (and area)
	TriangulationIntersectionVisitor(double iLength,int iVertexMapLength, double proportion = 0.01)
	{
		IntersectionTDS_2.InitializeIntersection(iVertexMapLength);
		CharacteristicLength = iLength*proportion;
		CharacteristicArea   = CharacteristicLength*CharacteristicLength;
	};

	// --- Visitor operators

	// Intersection is a point -> do nothing.
	void operator()(const Point_2& p)
	{
	};

	// Intersection is a segment -> do nothing.
	void operator()(const Segment_2& s)
	{
	};

	// Push back a triangle intersection.
	void operator()(const Triangle_2& t)
	{
		if(std::abs(t.area()) > CharacteristicArea)
		{
			// Then the triangle is "big enough"
			IntersectionTDS_2.AddPolygon(t);
		}
	};

	// Push back a polygon intersection.
	void operator()(const std::vector<Point_2>& q)
	{
		Polygon_2 	tempPolygon(q.begin(),q.end());
		if(std::abs(tempPolygon.area()) > CharacteristicArea)
		{
			// Then the triangle is "big enough"
			IntersectionTDS_2.AddPolygon(tempPolygon,q.size(),CharacteristicArea);
		}
	};

	// Finalize: remove the dummy first face
	void Finalize()
	{
		IntersectionTDS_2.Finalize();
		IntersectionTDS_2.CleanUp();
	};
};

//	Friend function - overloaded operator<< for 'PolyIntersectionVisitor'
std::ostream& operator<<(std::ostream& out, TriangulationIntersectionVisitor& input);

// -------- Functions ---------------------------------------------------------

/*
 * 			Brute force intersection finding algorithm for 2D triangulations.
 *
 * 			Runs over the whole meshes dtA and dtB, looking for the
 * 		intersections. Complexity is O(n*m), where n and m are the number of
 * 		triangles in each triangulation.
 *
 */
template<typename TemplateVisitor>
extern void BuildMeshIntersections_Brute(
		const Triangular_Mesh_2& dtA,
		const Triangular_Mesh_2& dtB,
		TemplateVisitor& output
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
template<typename TemplateVisitor>
extern void BuildMeshIntersections(
		const Triangular_Mesh_2& dtA,
		const Triangular_Mesh_2& dtB,
		TemplateVisitor& output
		);

/*
 * 			Find if there is an intersection between two triangles.
 *
 * 			If such an intersection exists, determinate it. The vector
 * 		candidatesAUpdate is needed for BuildMeshIntersections' algorithm
 * 		(detailed in the source file).
 *
 */
template<typename TemplateVisitor>
extern void IntersectTriangles(
		const Triangular_Mesh_2& 	dtB,
		const Face_handle_2& 		workingTriangleB,
		const Triangular_Mesh_2& 	dtA,
		const Face_handle_2& 		workingTriangleA,
		std::vector<int>& 			candidatesAUpdate,
		bool& 						queryIntersect,
		TemplateVisitor& 	output
		);

/*
 * 			Find an intersecting pair of triangles between two triangulations.
 *
 */
void FindFirstPair(
		const Triangular_Mesh_2& 	dtA,
		const Triangular_Mesh_2& 	dtB,
		Face_handle_2&				FirstA,
		Face_handle_2&				FirstB
		);

#endif /* INTERSECTION_FUNCTIONS_H_ */
