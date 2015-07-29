/*
 * intersection_functions_2.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "intersection_functions_2.h"

template<typename TemplateVisitor>
void BuildMeshIntersections_Brute(
		const Triangular_Mesh_2& dtA,
		const Triangular_Mesh_2& dtB,
		TemplateVisitor& output
		)
{

	// --- Input
	/*
	 * 		dtA, dtB	:	2D triangulations to be analyzed
	 * 		output		:	visitor class (see "Boost::variant::") that will be
	 * 						applied to the intersections
	 */

	// --- Preamble

	//	Starting and ending iterators for buth meshes
	Finite_face_iterator_2 finiteFacesBegin_A = dtA.mesh.finite_faces_begin();
	Finite_face_iterator_2 finiteFacesEnd_A = dtA.mesh.finite_faces_end();

	Finite_face_iterator_2 finiteFacesBegin_B = dtB.mesh.finite_faces_begin();
	Finite_face_iterator_2 finiteFacesEnd_B = dtB.mesh.finite_faces_end();

	//	Total number of intersections found
	int nbOfIntersections = 0;

	for(Finite_face_iterator_2 itB = finiteFacesBegin_B; itB != finiteFacesEnd_B; ++itB)
	{
		for(Finite_face_iterator_2 itA = finiteFacesBegin_A; itA != finiteFacesEnd_A; ++itA)
		{

			// Crossing
			Triangle_2_Intersection_Variant result = CGAL::intersection(dtA.mesh.triangle(itA),dtB.mesh.triangle(itB));

			if (result) {
				boost::apply_visitor(output, *result);
				++nbOfIntersections;
			}
		}
	}
};

template<typename TemplateVisitor>
void BuildMeshIntersections(
		const Triangular_Mesh_2& dtA,
		const Triangular_Mesh_2& dtB,
		 TemplateVisitor& output
		)
{

	// --- Input
	/*
	 * 		dtA, dtB	:	2D triangulations to be analyzed
	 * 		output		:	visitor class (see "Boost::variant::") that will be
	 * 						applied to the intersections
	 */

	// --- Notes
	/*
	 * 		1)	This algorithm was built with CGAL's "exact predicates, inexact
	 * 		constructions" kernel, and so the (cheap) "do_intersect" tests,
	 * 		linked to the algorithm's logic, are exact. On the other hand, the
	 * 		"intersection" operations (linked to the output) are inexact.
	 * 		This might result in the algorithm finding an intersection, but not
	 * 		 outputing it. The size of these missed intersections is of the
	 * 		 order of the "double" format rounding error (15 - 17 significant
	 * 		 decimal digits), and can be ignored.
	 *
	 * 		 	One can use CGAL's "exact predicates, exact constructions"
	 * 		 kernel to obtain all the intersections, but this comes at a
	 * 		 considerable resource and time cost.
	 *
	 * 		2)	The "vector_deque" data structure is essentially a vector,
	 * 		but with an auxiliary index, used to access its data, and "length".
	 * 		This allows to reuse the data structure without reallocating the
	 * 		memory (setting the "length" to zero). This data structure does
	 * 		fake but efficient "remove()" operations by moving the value of
	 * 		the last entry to the entry to be removed and reducing the
	 * 		artificial length. The order of the elements is, then, changed after
	 * 		each "remove()" operation.
	 *
	 */

	// --- Preamble - declarations

	// Total number of inexact intersections found (see note 1 above)
	int nbOfIntersections = 0;

	// Handles of the first triangles from dtA and dtB to be tested
	Face_handle_2 FirstA;
	Face_handle_2 FirstB;

	// Find the first pair
	FindFirstPair(dtA,dtB,FirstA,FirstB);

	// Dummy variables for the working triangle handles
	Face_handle_2 workingTriangleA;
	Face_handle_2 workingTriangleB;

	// Dummy variables for the neighboring triangle handles
	Face_handle_2 neighTriangleA;
	Face_handle_2 neighTriangleB;

	// Dummy variables containing the neighboring triangles indexes
	int neighAindex;
	int neighBindex;

	/*
	 *
	 *  	The outer while loop runs over the deque "MeshBQueue", extracting its first element
	 *  at each iteration, until the deque is empty (meaning that all triangles of "dtB" were
	 *  treated). It also extracts the (corresponding) first element from "MeshAQueue", and set it
	 *  as the first element of "MeshAToTest".
	 *
	 *
	 *
	 */

	/*
	 *
	 * 			Vector_deques containing the lists of the next triangles to be
	 * 		treated. If the algorithm is correct, these "vector_deques" are
	 * 		twins: the elements "MeshBQueue[iii]" and "MeshAQueue[iii]"
	 * 		correspond to a pair of triangles from, respectively, dtB and dtA
	 * 		that are guaranteed to intersect (exactly, see note 1 above).
	 * 		"MeshAQueue" and "MeshBQueue" correspond to "bil" and "bl" in
	 * 		Gander's article.
	 *
	 */
	vector_deque<Face_handle_2> MeshAQueue(dtA.get_nb_of_faces() + 1);
	vector_deque<Face_handle_2> MeshBQueue(dtB.get_nb_of_faces() + 1);

	/*
	 * 			Marker vector, used to indicate if a triangle from dtB was
	 * 		already treated (=1) or not (=0). The last term is a marker for ALL
	 * 		the infinite faces.
	 *
	 */
	std::vector<int> treatedFromB(dtB.get_nb_of_faces() + 1,0);

	/*
	 * 			Marker unordered_set, used to indicate if a triangle from dtA
	 * 		was trated (find()==true) or not (find()==false). dtA uses an
	 * 		unordered_set instead of a vector due to the need of efficiently
	 * 		emptying it before the start of the inner loop.
	 *
	 */
	std::unordered_set<int> treatedFromA(dtA.get_nb_of_faces() + 1);
	treatedFromB[dtB.get_nb_of_faces()] = 1;

	// Vector_deque of the elements from dtA that must be tested yet
	vector_deque<Face_handle_2>	MeshAToTest(dtA.get_nb_of_faces() + 1);

	// Handles of the elements that might be added to MeshAQueue
	std::vector<Face_handle_2> candidatesA(3);

	// Marker vectors, used to determinate if a triangle must be added
	// to MeshAQueue
	std::vector<int> candidatesAUpdate(3,0);
	std::vector<int> candidatesAIsOccupied(3,0);

	// Boolean saying if two triangles intersect (exactly)
	bool queryIntersect = false;

	// --- Preamble - initializations
	// Insert the first elements in the queues
	MeshBQueue.add(FirstB);
	MeshAQueue.add(FirstA);

	// Mark the first element from dtB as "treated"
	treatedFromB[FirstB->info().ExtIndex] = 1;

	// --- Preamble - timing variables
	std::chrono::time_point<std::chrono::system_clock> 	timing_start,
														timing_end;

	std::chrono::duration<double> 	elapsed_seconds_total;
	timing_start = std::chrono::system_clock::now();

	// Number of operations per cycle of the external loop
	int DebugNumberOfOperations = 0;

	while(!MeshBQueue.isEmpty())
	{
		// Pop out working triangle from mesh B
		workingTriangleB = MeshBQueue[0];
		MeshBQueue.remove(0);

		// Clear "MeshAToTest" and initialize it with the first element from
		// "MeshAQueue"
		MeshAToTest.clearAll();
		MeshAToTest.add(MeshAQueue[0]);
		MeshAQueue.remove(0);

		// Clear "treatedFromA", mark the first element and the infinite
		// triangles as "Treated"
		treatedFromA.erase(treatedFromA.begin(),treatedFromA.end());
		treatedFromA.insert(dtA.get_nb_of_faces());
		treatedFromA.insert(MeshAToTest[0]->info().ExtIndex);

		// Reset the 'candidatesAIsOccupied' variables
		candidatesAIsOccupied[0] = 0;
		candidatesAIsOccupied[1] = 0;
		candidatesAIsOccupied[2] = 0;

		while(!MeshAToTest.isEmpty())
		{
			// Pop out working triangle from mesh A's test triangles
			workingTriangleA = MeshAToTest[0];
			MeshAToTest.remove(0);

			queryIntersect = false;
			IntersectTriangles(
					dtB,
					workingTriangleB,
					dtA,
					workingTriangleA,
					candidatesAUpdate,
					queryIntersect,
					output
				);

			/* 		If they intersect, we must a) test over workingTriangleA's
			 * 	neighbors, and b) set up candidates for the next crossings
			 */
			if(queryIntersect)
			{
				++nbOfIntersections;
				for(int iii = 0; iii < 3; ++iii)
				{
					// Add neighbors to "MeshAToTest"
					neighTriangleA = workingTriangleA->neighbor(iii);
					neighAindex = neighTriangleA->info().ExtIndex;

					if(treatedFromA.find(neighAindex)==treatedFromA.end())
					{
						// New triangle!
						MeshAToTest.add(neighTriangleA);
						treatedFromA.insert(neighAindex);
					}

					// Set up candidates for the next crossings
					if(candidatesAUpdate[iii]==1)
					{
						candidatesAIsOccupied[iii] = 1;
						candidatesA[iii] = workingTriangleA;
					}
				}
			}
			++DebugNumberOfOperations;
		}
		// End of inner loop

		// Must now update the "MeshBQueue" vector
		for(int iii = 0; iii < 3; ++iii)
		{
			neighTriangleB = workingTriangleB->neighbor(iii);
			neighBindex = neighTriangleB->info().ExtIndex;

			if(treatedFromB[neighBindex]==0 && candidatesAIsOccupied[iii]==1)
			{
				// New triangle!
				MeshBQueue.add(neighTriangleB);
				MeshAQueue.add(candidatesA[iii]);
				treatedFromB[neighBindex]=1;
			}
		}
	}

	// --- Finish - Timing and debug
	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_total = timing_end-timing_start;

	std::cout 	<< dtA.get_nb_of_faces() << " "
				<< dtB.get_nb_of_faces() << " "
				<< nbOfIntersections << " "
				<< DebugNumberOfOperations << " "
				<< double(DebugNumberOfOperations)/dtB.get_nb_of_faces()  << " "
				<< elapsed_seconds_total.count() << std::endl;

};

template<typename TemplateVisitor>
void IntersectTriangles(
		const Triangular_Mesh_2& 	dtB,
		const Face_handle_2& 		workingTriangleB,
		const Triangular_Mesh_2& 	dtA,
		const Face_handle_2& 		workingTriangleA,
		std::vector<int>& 			candidatesAUpdate,
		bool& 						queryIntersect,
		TemplateVisitor& 	output
		)
{
	// Test exact intersection between "workingTriangleB" and "workingTriangleA"
	queryIntersect = CGAL::do_intersect(
										dtA.mesh.triangle(workingTriangleA),
										dtB.mesh.triangle(workingTriangleB)
									);

	// If true, build inexact intersection
	if(queryIntersect)
	{
		Triangle_2_Intersection_Variant result;
		result = CGAL::intersection(
									dtA.mesh.triangle(workingTriangleA),
									dtB.mesh.triangle(workingTriangleB)
								);
		if (result)
		{
			boost::apply_visitor(output, *result);
		}
	}

	// Test exact intersection with "workingTriangleB" 's neighbors
	bool neighTest = false;
	for(int iii = 0; iii < 3; ++iii)
	{
		neighTest = CGAL::do_intersect(dtB.mesh.triangle(workingTriangleB->neighbor(iii)),
				dtA.mesh.triangle(workingTriangleA));

		if(neighTest)
		{
			candidatesAUpdate[iii] = 1;
		}
		else
		{
			candidatesAUpdate[iii] = 0;
		}
	}
}

void FindFirstPair(
		const Triangular_Mesh_2& 	dtA,
		const Triangular_Mesh_2& 	dtB,
		Face_handle_2&				FirstA,
		Face_handle_2&				FirstB
		)
{
	Face_handle_2	tempAFace;
	Face_circulator_2	tempBFace_begin;
	Face_circulator_2	tempBFace;

	FirstB = NULL;

	for(Finite_vertices_iterator_2 itBVertex = dtB.mesh.finite_vertices_begin();
			itBVertex != dtB.mesh.finite_vertices_end(); ++itBVertex)
	{
		tempAFace = dtA.mesh.locate(itBVertex->point());
		if(tempAFace!=NULL && !dtA.mesh.is_infinite(tempAFace))
		{

			// There is a intersection AND it's not with a infinite face
			FirstA = tempAFace;
			tempBFace_begin = dtB.mesh.incident_faces(itBVertex);
			tempBFace = tempBFace_begin;

			if(!dtB.mesh.is_infinite(tempBFace_begin))
			{
				FirstB = tempBFace_begin;
			}
			else
			{
				++tempBFace;
				for( ; tempBFace != tempBFace_begin; ++tempBFace)
				{
					if(!dtB.mesh.is_infinite(tempBFace))
					{
						FirstB = tempBFace;
						break;
					}
				}
			}
			break;
		}
	}
}

std::ostream& operator<<(
		std::ostream& out,
		IntersectionTallyVisitor& input
		)
{
	std::vector<Point_2>::iterator ptItBegin = input.point_intersections.begin();
	std::vector<Segment_2>::iterator segItBegin = input.seg_intersections.begin();
	std::vector<Triangle_2>::iterator triItBegin = input.tri_intersections.begin();
	std::vector<std::vector<Point_2> >::iterator polyItBegin = input.poly_intersections.begin();
	std::vector<Point_2>::iterator vecItBegin;

	out << "Point intersections : " << std::endl;
	for(std::vector<Point_2>::iterator ptIt = ptItBegin;
			ptIt != input.point_intersections.end(); ++ptIt)
	{
		out <<  *ptIt << std::endl;
	}

	out << std::endl;
	out << "Segment intersections : " << std::endl;
	for(std::vector<Segment_2>::iterator segIt = segItBegin;
			segIt != input.seg_intersections.end(); ++segIt)
	{
		out << segIt->source() << " " <<  segIt->target() << std::endl;
	}

	out << std::endl;
	out << "Triangle intersections : " << std::endl;
	for(std::vector<Triangle_2>::iterator triIt = triItBegin;
			triIt != input.tri_intersections.end(); ++triIt)
	{
		out <<  triIt->vertex(0) << " " <<  triIt->vertex(1) << " "
				<<  triIt->vertex(2) << std::endl;
	}

	out << std::endl;
	out << "Polygon intersections : " << std::endl;
	for(std::vector<std::vector<Point_2> >::iterator polyIt = polyItBegin;
			polyIt != input.poly_intersections.end(); ++polyIt)
	{
		vecItBegin = polyIt->begin();
		for(std::vector<Point_2>::iterator vecIt = vecItBegin;
				vecIt != polyIt->end()-1; ++vecIt)
		{
			out <<  *vecIt << " ";
		}
		out <<  *(polyIt->end()-1) << std::endl;
	}
	return out;
};

std::ostream& operator<<(
		std::ostream& out,
		PolyIntersectionVisitor& input
		)
{
	std::vector<Triangle_2>::iterator triItBegin = input.tri_intersections.begin();
	std::vector<std::vector<Point_2> >::iterator polyItBegin = input.poly_intersections.begin();
	std::vector<Point_2>::iterator vecItBegin;

	for(std::vector<Triangle_2>::iterator triIt = triItBegin;
			triIt != input.tri_intersections.end(); ++triIt)
	{
		out <<  triIt->vertex(0) << " " <<  triIt->vertex(1) << " "
				<<  triIt->vertex(2) << std::endl;
	}

	for(std::vector<std::vector<Point_2> >::iterator polyIt = polyItBegin;
			polyIt != input.poly_intersections.end(); ++polyIt)
	{
		vecItBegin = polyIt->begin();
		for(std::vector<Point_2>::iterator vecIt = vecItBegin;
				vecIt != polyIt->end()-1; ++vecIt)
		{
			out <<  *vecIt << " ";
		}
		out <<  *(polyIt->end()-1) << std::endl;
	}
	return out;
};

// --- Template specializations
/*
 * 			Needed to allow template implementations inside a .cpp file. A new
 * 		ensemble of declarations must be done for each new visitor.
 */

// Template visitor == ntersectionTallyVisitor
template void BuildMeshIntersections_Brute<IntersectionTallyVisitor>(
		const Triangular_Mesh_2& dtA,
		const Triangular_Mesh_2& dtB,
		IntersectionTallyVisitor& output
		);

template void BuildMeshIntersections<IntersectionTallyVisitor>(
		const Triangular_Mesh_2& dtA,
		const Triangular_Mesh_2& dtB,
		IntersectionTallyVisitor& output
		);

template void IntersectTriangles<IntersectionTallyVisitor>(
		const Triangular_Mesh_2& 	dtB,
		const Face_handle_2& 		workingTriangleB,
		const Triangular_Mesh_2& 	dtA,
		const Face_handle_2& 		workingTriangleA,
		std::vector<int>& 			candidatesAUpdate,
		bool& 						queryIntersect,
		IntersectionTallyVisitor& 	output
		);

// Template visitor == PolyIntersectionVisitor
template void BuildMeshIntersections_Brute<PolyIntersectionVisitor>(
		const Triangular_Mesh_2& dtA,
		const Triangular_Mesh_2& dtB,
		PolyIntersectionVisitor& output
		);

template void BuildMeshIntersections<PolyIntersectionVisitor>(
		const Triangular_Mesh_2& dtA,
		const Triangular_Mesh_2& dtB,
		PolyIntersectionVisitor& output
		);

template void IntersectTriangles<PolyIntersectionVisitor>(
		const Triangular_Mesh_2& 	dtB,
		const Face_handle_2& 		workingTriangleB,
		const Triangular_Mesh_2& 	dtA,
		const Face_handle_2& 		workingTriangleA,
		std::vector<int>& 			candidatesAUpdate,
		bool& 						queryIntersect,
		PolyIntersectionVisitor& 	output
		);
