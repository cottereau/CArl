/*
 * intersection_functions_2.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "intersection_functions_3.h"

bool tetra_intersection(
						const Triangular_Mesh_3& 	dtA,
						const Cell_handle_3 		cellA,
						const Triangular_Mesh_3& 	dtB,
						const Cell_handle_3 		cellB,
						Polyhedron& polyOut
						)
{
	// Fast intersection test, using bbox
	bool bPreliminary = CGAL::do_intersect(
								dtA.mesh.tetrahedron(cellA).bbox(),
								dtB.mesh.tetrahedron(cellB).bbox()
								);

	IntersectionPointsVisitor_3		visPointList;

	bool bTriangleIntersect = false;
	if(bPreliminary)
	{
		// If true, must do real test
		Triangle_3_Intersection_Variant varTriangleIntersection;

		for(int iii = 0; iii < 4; ++iii)
		{
			for(int jjj =0; jjj < 4; ++jjj)
			{
				varTriangleIntersection = CGAL::intersection(
									dtA.mesh.triangle(cellA,iii),
									dtB.mesh.triangle(cellB,jjj)
									);

				if(varTriangleIntersection)
				{
					bTriangleIntersect = true;
					boost::apply_visitor(visPointList, *varTriangleIntersection);
				}
			}
		}


		if(bTriangleIntersect)
		{
			// --- Add the corner points
			CGAL::Bounded_side pointLocation;

			// The dummy variables are needed for using side_of_cell ...
			Locate_type_3 dummyLT;
			int dummyLI, dummyLJ;

			for(int iii = 0; iii < 4; ++ iii)
			{
				pointLocation = dtA.mesh.side_of_cell(
													dtB.mesh.point(cellB,iii),
													cellA,
													dummyLT,
													dummyLI,
													dummyLJ
													);

				if(pointLocation==CGAL::ON_BOUNDED_SIDE || pointLocation==CGAL::ON_BOUNDARY)
				{
					visPointList.add(dtB.mesh.point(cellB,iii));
				}

				pointLocation = dtB.mesh.side_of_cell(
													dtA.mesh.point(cellA,iii),
													cellB,
													dummyLT,
													dummyLI,
													dummyLJ
													);

				if(pointLocation==CGAL::ON_BOUNDED_SIDE || pointLocation==CGAL::ON_BOUNDARY)
				{
					visPointList.add(dtA.mesh.point(cellA,iii));
				}
			}

			// Build the polyhedron - if there are at least 4 points!
			if(visPointList.IntersectionData.size()>3)
			{
				CGAL::convex_hull_3(
									visPointList.IntersectionData.begin(),
									visPointList.IntersectionData.end(),
									polyOut
									);
			}
			else
			{
				// Well, this was useless ...
				bTriangleIntersect = false;
			}
		}
	}

	return bTriangleIntersect;
};

bool tetra_do_intersect(
						const Triangular_Mesh_3& 	dtA,
						const Cell_handle_3 		cellA,
						const Triangular_Mesh_3& 	dtB,
						const Cell_handle_3 		cellB
						)
{
	// Fast intersection test, using bbox
	bool bPreliminary = CGAL::do_intersect(
								dtA.mesh.tetrahedron(cellA).bbox(),
								dtB.mesh.tetrahedron(cellB).bbox()
								);

	bool bTriangleIntersect = false;
	// If true, must do real test
	if(bPreliminary)
	{
		for(int iii = 0; iii < 4; ++iii)
		{
			bTriangleIntersect =
					CGAL::do_intersect(
									   dtA.mesh.tetrahedron(cellA),
									   dtB.mesh.triangle(cellB,iii)
									   );
			if(bTriangleIntersect)
			{
				break;
			}
		}
	}

	return bTriangleIntersect;
};

void BuildMeshIntersections_Brute(
		const Triangular_Mesh_3& dtA,
		const Triangular_Mesh_3& dtB,
		std::vector<Polyhedron>& 	output
		)
{

	// --- Input
	/*
	 * 		dtA, dtB	:	3D triangulations to be analyzed
	 * 		output		:	visitor class (see "Boost::variant::") that will be
	 * 						applied to the intersections
	 */

	// --- Preamble

	//	Starting and ending iterators for both meshes
	Finite_cells_iterator_3 finiteCellsBegin_A = dtA.mesh.finite_cells_begin();
	Finite_cells_iterator_3 finiteCellsEnd_A = dtA.mesh.finite_cells_end();

	Finite_cells_iterator_3 finiteCellsBegin_B = dtB.mesh.finite_cells_begin();
	Finite_cells_iterator_3 finiteCellsEnd_B = dtB.mesh.finite_cells_end();

	//	Total number of intersections found
	int nbOfIntersections = 0;

	//  Result initialization
	bool result = false;

	// 	Dummy polyhedron
	Polyhedron dummyPoly;

	// Clean up output
	output.clear();

	for(Finite_cells_iterator_3 itB = finiteCellsBegin_B; itB != finiteCellsEnd_B; ++itB)
	{
		for(Finite_cells_iterator_3 itA = finiteCellsBegin_A; itA != finiteCellsEnd_A; ++itA)
		{
			// Crossing
			result = tetra_intersection(dtA,itA,dtB,itB,dummyPoly);

			if (result) {
				output.push_back(dummyPoly);
				++nbOfIntersections;
			}
		}
	}
};

void BuildMeshIntersections(
		const Triangular_Mesh_3& 	dtA,
		const Triangular_Mesh_3& 	dtB,
		std::vector<Polyhedron>& 	output
		)
{

	// --- Input
	/*
	 * 		dtA, dtB	:	3D triangulations to be analyzed
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
	 */

	// --- Preamble - declarations

	// Total number of inexact intersections found (see note 1 above)
	int nbOfIntersections = 0;

	// Handles of the first tetrahedrons from dtA and dtB to be tested
	Cell_handle_3 FirstA;
	Cell_handle_3 FirstB;

	// Find the first pair
	FindFirstPair(dtA,dtB,FirstA,FirstB);

	// Dummy variables for the working tetrahedron handles
	Cell_handle_3 workingTetrahedronA;
	Cell_handle_3 workingTetrahedronB;

	// Dummy variables for the neighboring tetrahedron handles
	Cell_handle_3 neighTetrahedronA;
	Cell_handle_3 neighTetrahedronB;

	// Dummy variables containing the neighboring tetrahedrons indexes
	int neighAindex;
	int neighBindex;

	/*
	 *
	 * 			Deques containing the lists of the next tetrahedrons to be treated.
	 * 		If the algorithm is correct, these "deques" are twins: the elements
	 * 		"MeshBQueue[iii]" and "MeshAQueue[iii]" correspond to a pair of
	 * 		tetrahedrons from, respectively, dtB and dtA that are guaranteed to
	 * 		intersect (exactly, see note 1 above). "MeshAQueue" and "MeshBQueue"
	 * 		correspond to "bil" and "bl" in Gander's article. Deques were used
	 * 		for these data lists due to the need to push_back() and pop_front()
	 * 		elements efficiently.
	 *
	 */
	std::deque<Cell_handle_3> MeshAQueue;
	std::deque<Cell_handle_3> MeshBQueue;

	/*
	 * 			Marker vector, used to indicate if a tetrahedron from dtB was
	 * 		already treated (=1) or not (=0). The last term is a marker for ALL
	 * 		the infinite faces.
	 *
	 */
	std::vector<int> treatedFromB(dtB.get_nb_of_cells() + 1,0);

	/*
	 * 			Marker unordered_set, used to indicate if a tetrahedron from dtA
	 * 		was trated (find()==true) or not (find()==false). dtA uses an
	 * 		unordered_set instead of a vector or deque due to the need of
	 * 		efficient insertion, search and emptying of the data. The number of
	 * 		buckets chosen (dtA.get_nb_of_faces() + 1) guarantees that there
	 * 		will be no collisions.
	 *
	 */
	std::unordered_set<int> treatedFromA(dtA.get_nb_of_cells() + 1);
	treatedFromB[dtB.get_nb_of_cells()] = 1;

	// Deque of the elements from dtA that must be tested yet
	std::deque<Cell_handle_3>	MeshAToTest;

	// Handles of the elements that might be added to MeshAQueue
	std::vector<Cell_handle_3> candidatesA(3);

	// Marker vectors, used to determinate if a tetrahedron must be added
	// to MeshAQueue
	std::vector<int> candidatesAUpdate(3,0);
	std::vector<int> candidatesAIsOccupied(3,0);

	// Boolean saying if two tetrahedrons intersect (exactly)
	bool queryIntersect = false;

	// Boolean saying if two tetrahedrons intersect (inexactly)
	bool inexactQueryIntersect = false;

	// --- Preamble - initializations
	// Insert the first elements in the queues
	MeshAQueue.push_back(FirstA);
	MeshBQueue.push_back(FirstB);

	// Mark the first element from dtB as "treated"
	treatedFromB[FirstB->info().ExtIndex] = 1;

	// --- Preamble - timing variables
	std::chrono::time_point<std::chrono::system_clock> 	timing_start,
														timing_end;

	std::chrono::duration<double> 	elapsed_seconds_total;
	timing_start = std::chrono::system_clock::now();

	// Number of operations per cycle of the external loop
	int DebugNumberOfOperations = 0;

	// Dummy polyhedron
	Polyhedron dummyPoly;

	while(!MeshBQueue.empty())
	{
		// Pop out working tetrahedron from mesh B
		workingTetrahedronB = MeshBQueue[0];
		MeshBQueue.pop_front();

		// Clear "MeshAToTest" and initialize it with the first element from
		// "MeshAQueue"
		MeshAToTest.clear();
		MeshAToTest.push_back(MeshAQueue[0]);
		MeshAQueue.pop_front();

		// Clear "treatedFromA", mark the first element and the infinite
		// tetrahedrons as "Treated"
		treatedFromA.erase(treatedFromA.begin(),treatedFromA.end());
		treatedFromA.insert(dtA.get_nb_of_cells());
		treatedFromA.insert(MeshAToTest[0]->info().ExtIndex);

		// Reset the 'candidatesAIsOccupied' variables
		candidatesAIsOccupied[0] = 0;
		candidatesAIsOccupied[1] = 0;
		candidatesAIsOccupied[2] = 0;

		while(!MeshAToTest.empty())
		{
			// Pop out working tetrahedron from mesh A's test tetrahedrons
			workingTetrahedronA = MeshAToTest[0];
			MeshAToTest.pop_front();

			queryIntersect = false;
			IntersectTetrahedrons(
					dtB,
					workingTetrahedronB,
					dtA,
					workingTetrahedronA,
					candidatesAUpdate,
					queryIntersect,
					inexactQueryIntersect,
					dummyPoly
				);

			/* 		If they intersect, we must a) test over workingTetrahedronA's
			 * 	neighbors, and b) set up candidates for the next crossings
			 */
			if(queryIntersect)
			{
				++nbOfIntersections;
				if(inexactQueryIntersect)
				{
					output.push_back(dummyPoly);
				}

				for(int iii = 0; iii < 4; ++iii)
				{
					// Add neighbors to "MeshAToTest"
					neighTetrahedronA = workingTetrahedronA->neighbor(iii);
					neighAindex = neighTetrahedronA->info().ExtIndex;

					if(treatedFromA.find(neighAindex)==treatedFromA.end())
					{
						// New tetrahedron!
						MeshAToTest.push_back(neighTetrahedronA);
						treatedFromA.insert(neighAindex);
					}

					// Set up candidates for the next crossings
					if(candidatesAUpdate[iii]==1)
					{
						candidatesAIsOccupied[iii] = 1;
						candidatesA[iii] = workingTetrahedronA;
					}
				}
			}
			++DebugNumberOfOperations;
		}
		// End of inner loop

		// Must now update the "MeshBQueue" vector
		for(int iii = 0; iii < 4; ++iii)
		{
			neighTetrahedronB = workingTetrahedronB->neighbor(iii);
			neighBindex = neighTetrahedronB->info().ExtIndex;

			if(treatedFromB[neighBindex]==0 && candidatesAIsOccupied[iii]==1)
			{
				// New tetrahedron!
				MeshBQueue.push_back(neighTetrahedronB);
				MeshAQueue.push_back(candidatesA[iii]);
				treatedFromB[neighBindex]=1;
			}
		}
	}

	// --- Finish - Timing and debug
	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_total = timing_end-timing_start;

	std::cout 	<< dtA.get_nb_of_cells() << " "
				<< dtB.get_nb_of_cells() << " "
				<< nbOfIntersections << " "
				<< DebugNumberOfOperations << " "
				<< double(DebugNumberOfOperations)/dtB.get_nb_of_cells()  << " "
				<< elapsed_seconds_total.count() << std::endl;

};

void IntersectTetrahedrons(
		const Triangular_Mesh_3& 	dtB,
		const Cell_handle_3& 		workingTetrahedronB,
		const Triangular_Mesh_3& 	dtA,
		const Cell_handle_3& 		workingTetrahedronA,
		std::vector<int>& 			candidatesAUpdate,
		bool& 						queryIntersect,
		bool&						inexactQueryIntersect,
		Polyhedron& 				output
		)
{
	// Test exact intersection between "workingTetrahedronB" and "workingTetrahedronA"
	queryIntersect = tetra_do_intersect(
										dtA,
										workingTetrahedronA,
										dtB,
										workingTetrahedronB
									   );

	// If true, build inexact intersection
	if(queryIntersect)
	{
		inexactQueryIntersect = tetra_intersection(
										dtA,
										workingTetrahedronA,
										dtB,
										workingTetrahedronB,
										output
									   );
	}

	// Test exact intersection with "workingTetrahedronB" 's neighbors
	bool neighTest = false;
	for(int iii = 0; iii < 3; ++iii)
	{
		neighTest = tetra_do_intersect(dtB,workingTetrahedronB->neighbor(iii),
				dtA,workingTetrahedronA);

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
		const Triangular_Mesh_3& 	dtA,
		const Triangular_Mesh_3& 	dtB,
		Cell_handle_3&				FirstA,
		Cell_handle_3&				FirstB
		)
{
	Cell_handle_3				tempACell;
	std::vector<Cell_handle_3>	BCellsToTest;

	std::vector<Cell_handle_3>::iterator	BCellsToTest_begin;

	FirstB = NULL;

	for(Finite_vertices_iterator_3 itBVertex = dtB.mesh.finite_vertices_begin();
			itBVertex != dtB.mesh.finite_vertices_end(); ++itBVertex)
	{
		tempACell = dtA.mesh.locate(itBVertex->point());
		if(tempACell!=NULL && !dtA.mesh.is_infinite(tempACell))
		{

			// There is a intersection AND it's not with a infinite face
			FirstA = tempACell;
			dtB.mesh.incident_cells(itBVertex,BCellsToTest.begin());

			BCellsToTest_begin = BCellsToTest.begin();
			for(std::vector<Cell_handle_3>::iterator itVec = BCellsToTest.begin();
					itVec != BCellsToTest.end(); ++itVec)
			{
				if(!dtB.mesh.is_infinite(*itVec))
				{
					FirstB = *itVec;
					break;
				}
			}
			break;
		}
	}
}

std::ostream& operator<<(
		std::ostream& out,
		IntersectionPointsVisitor_3& input
		)
{

	for(std::vector<Point_3>::iterator vecIt = input.IntersectionData.begin();
			vecIt != input.IntersectionData.end(); ++vecIt)
	{
		out <<  *vecIt << " ";
	}
	out << std::endl;
	return out;
};
