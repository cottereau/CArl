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
#include "triangular_mesh_intersection_3.h"

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
 * 		intersections between two tetrahedrons. It also contains the data
 * 		structures used during the intersection operation, such as the Nef
 * 		polyhedrons (see https://en.wikipedia.org/wiki/Nef_polygon), and the
 * 		exact/inexact data structure converters.
 *
 * 		TODO: adapt the output operator overload (or build a new function) such
 * 		that it is well formated for FreeFem++
 *
 */
class IntersectionPointsVisitor_3
	: public boost::static_visitor<void>
{
protected:
	// --- Constructors
	// Empty constructor.
	IntersectionPointsVisitor_3()
	{
		dummyIndex = 0;
	};

	// --- Members

	int dummyIndex;

	Nef_Polyhedron NefPolyA;
	Nef_Polyhedron NefPolyB;
	Nef_Polyhedron NefPolyI;

	ExactPolyhedron exactDummyPolyA;
	ExactPolyhedron exactDummyPolyB;

	std::vector<ExactPoint_3> exactVerticesA;
	std::vector<ExactPoint_3> exactVerticesB;

	Point_3	inexactPointI;

	Kernel_to_ExactKernel ConvertInexactToExact;
	ExactKernel_to_Kernel ConvertExactToInexact;

	// Vector containing all the intersection points
	std::vector<Point_3>	IntersectionData;

public:

	// --- Constructors
	// Constructor reserving a space "reserveSize" for the vectors.
	IntersectionPointsVisitor_3(int reserveSize)
	{
		IntersectionData.resize(reserveSize);
		exactVerticesA.resize(4);
		exactVerticesB.resize(4);
		dummyIndex = 0;
	};

	// --- Visitor operators
	// Intersection is a point
	void operator()(const Point_3& p)
	{
		add(p);
	};

	// Intersection is a segment
	void operator()(const Segment_3& s)
	{
		add(s.source());
		add(s.target());
	};

	// Intersection is a triangle
	void operator()(const Triangle_3& t)
	{
		add(t.vertex(0));
		add(t.vertex(1));
		add(t.vertex(2));
	};

	// Intersection is a polygon
	void operator()(const std::vector<Point_3>& q)
	{

		for(std::vector<Point_3>::const_iterator itPoly = q.begin(); itPoly != q.end(); ++itPoly)
		{
			add(*itPoly);
		}
	};

	// Return information from the visitor (size of the vectors)
	std::string printStatics()
	{
		return 	std::to_string(IntersectionData.size());
	}

	void convertCellsToExact(Cell_handle_3 cellA, Cell_handle_3 cellB)
	{
		for(int iii = 0;iii < 4; ++iii)
		{
			exactVerticesA[iii] = ConvertInexactToExact(cellA->vertex(iii)->point());
			exactVerticesB[iii] = ConvertInexactToExact(cellB->vertex(iii)->point());
		}

		exactDummyPolyA.clear();
		exactDummyPolyA.make_tetrahedron(	exactVerticesA[0],
											exactVerticesA[1],
											exactVerticesA[2],
											exactVerticesA[3]);

//		std::cout << exactVerticesA[0] << " "
//				<< exactVerticesA[1] << " "
//				<< exactVerticesA[2] << " "
//				<< exactVerticesA[3] << std::endl;

		NefPolyA = Nef_Polyhedron(exactDummyPolyA);

		exactDummyPolyB.clear();
		exactDummyPolyB.make_tetrahedron(	exactVerticesB[0],
											exactVerticesB[1],
											exactVerticesB[2],
											exactVerticesB[3]);

//		std::cout << exactVerticesB[0] << " "
//				<< exactVerticesB[1] << " "
//				<< exactVerticesB[2] << " "
//				<< exactVerticesB[3] << std::endl;

		NefPolyB = Nef_Polyhedron(exactDummyPolyB);
	}

	bool intersect(Polyhedron& polyOut)
	{
		NefPolyI = NefPolyA*NefPolyB;

		if(NefPolyI.number_of_volumes() == 2)
		{
			// Then we have an intersection composed by only one polyhedron
			// (the second volume is the polyhedron's exterior)
			for(Nef_Polyhedron::Vertex_const_iterator
											itVertex = NefPolyI.vertices_begin();
											itVertex != NefPolyI.vertices_end();
											++itVertex)
			{
				inexactPointI = ConvertExactToInexact(itVertex->point());
				add(inexactPointI);
			}

			std::vector<Point_3>::iterator visPointListBegin = IntersectionData.begin();
			std::vector<Point_3>::iterator visPointListEnd = visPointListBegin;
			std::advance(visPointListEnd,size());

			CGAL::convex_hull_3(visPointListBegin,visPointListEnd,polyOut);
			return true;
		}

		return false;
	}

	// Operators to access/change the IntersectionData vector
	std::vector<Point_3>::iterator dataBegin()
	{
		return IntersectionData.begin();
	}

	std::vector<Point_3>::iterator dataEnd()
	{
		return IntersectionData.end();
	}

	void clear()
	{
		dummyIndex = 0;
	}

	int size()
	{
		return dummyIndex;
	}

	void add(const Point_3& p)
	{
		IntersectionData[dummyIndex] = p;
		++dummyIndex;
	}

	void PrintData()
	{
		std::vector<Point_3>::iterator listBegin = IntersectionData.begin();
		std::vector<Point_3>::iterator listEnd = listBegin;
		std::advance(listEnd,size());

		// Test for outlier points
		for(std::vector<Point_3>::iterator 	sillyIt = listBegin;
											sillyIt != listEnd;
											++sillyIt)
		{
			std::cout	<< "(" << *sillyIt << ") ";
		}
		std::cout	<< std::endl << std::endl;
	}
};

//	Friend function - overloaded operator<< for 'IntersectionPointsVisitor_3'
std::ostream& operator<<(std::ostream& out, IntersectionPointsVisitor_3& input);

/*
 * 			TriangulationIntersectionVisitor_3
 *
 * 			This is not really a visitor, but rather a class needed for saving
 * 		the intersection mesh information.
 *
 */
class TriangulationIntersectionVisitor_3
	: public boost::static_visitor<void>
{
private:
	// --- Constructor

	// Empty constructor - set as private
	TriangulationIntersectionVisitor_3()
	{
		CharacteristicLength = -1;
		CharacteristicVolume = -1;
	};
public:
	// --- Members

	Intersection_Mesh_3	IntersectionTDS_3;
	double				CharacteristicLength;
	double				CharacteristicVolume;

	// --- Constructors
	// Constructor with an characteristic length (and area)
	TriangulationIntersectionVisitor_3(double iLength,int iVertexMapLength, double proportion = 0.01)
	{
		IntersectionTDS_3.InitializeIntersection(iVertexMapLength);
		CharacteristicLength = iLength*proportion;
		CharacteristicVolume  = CharacteristicLength*CharacteristicLength*CharacteristicLength;
	};

	void InsertPolyhedron(Polyhedron& poly)
	{
		IntersectionTDS_3.AddPolyhedron(poly,CharacteristicVolume);
	}

	// Finalize: remove the dummy first face
	void Finalize()
	{
		IntersectionTDS_3.Finalize();
		IntersectionTDS_3.CleanUp();
	};
};

//	Friend function - overloaded operator<< for 'PolyIntersectionVisitor'
std::ostream& operator<<(std::ostream& out, TriangulationIntersectionVisitor_3& input);

// -------- Functions ---------------------------------------------------------

/*
 * 			Brute force intersection finding algorithm for 3D triangulations.
 *
 * 			Runs over the whole meshes dtA and dtB, looking for the
 * 		intersections. Complexity is O(n*m), where n and m are the number of
 * 		triangles in each triangulation.
 *
 */
void BuildMeshIntersections_Brute(
		const Triangular_Mesh_3& dtA,
		const Triangular_Mesh_3& dtB,
		TriangulationIntersectionVisitor_3& 	output
		);

/*
 * 			Optimized intersection finding algorithm for 3D triangulations.
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
		TriangulationIntersectionVisitor_3& 	output
		);

/*
 * 			Find if there is an intersection between two tetrahedrons.
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
		Polyhedron& 				output,
		IntersectionPointsVisitor_3&	visPointList
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

/*
 * 			Build intersection between two tetrahedrons (if it exists). This
 * 		function does this by
 *
 * 			a) extracting the (inexact) tetrahedrons from the cells;
 * 			b) converting them to (exact) Nef polyhedrons;
 * 			c) applying a boolean operation between them;
 * 			d) converting the resulting (exact) Nef polyhedron into an
 * 			   (inexact) polyhedron and point list.
 *
 * 			The conversions are needed due to the incompatibility between CGAL's
 * 		intersection operations and the inexact constructions kernels.
 *
 */
bool tetra_intersection(
		const Triangular_Mesh_3& 	dtA,
		const Cell_handle_3 		cellA,
		const Triangular_Mesh_3& 	dtB,
		const Cell_handle_3 		cellB,
		Polyhedron& polyOut,
		IntersectionPointsVisitor_3&	visPointList
		);

/*
 * 			Test the existence of an intersection between two tetrahedrons. This
 * 		function does not use the Nef polyhedrons conversions, since no
 * 		intersections are build, and hence we use exact predicates.
 */
bool tetra_do_intersect(
		const Triangular_Mesh_3& 	dtA,
		const Cell_handle_3 		cellA,
		const Triangular_Mesh_3& 	dtB,
		const Cell_handle_3 		cellB
		);

#endif /* INTERSECTION_FUNCTIONS_H_ */
