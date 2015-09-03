/*
 * main.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"

int main(int argc, char *argv[])
{

	// ************ //
	// Preamble     //
	// ************ //

	// Declare meshes
	std::string filenameA = "data/test_mesh_A.msh";
	std::string filenameB = "data/test_mesh_B.msh";
	std::string filenameOutput = "data/test_mesh_intersection.mesh";

	if(argc == 3)
	{
		filenameA = argv[1];
		filenameB = argv[2];
	}

	if(argc == 4)
	{
		filenameA = argv[1];
		filenameB = argv[2];
		filenameOutput = argv[3];
	}

	// --- Preamble - timing variables
	std::chrono::time_point<std::chrono::system_clock> 	timing_start,
														timing_end,
														timing_start_total,
														timing_end_total;

	std::chrono::duration<double> 	elapsed_seconds_input,
									elapsed_seconds_intersections,
									elapsed_seconds_export,
									elapsed_seconds_total;

	// Import meshes
	timing_start 		= std::chrono::system_clock::now();
	timing_start_total 	= std::chrono::system_clock::now();

	std::string sillyname;
	sillyname = "Testing A";
	Triangular_Mesh_2 dtA(sillyname,filenameA);

	sillyname = "Testing B";
	Triangular_Mesh_2 dtB(sillyname,filenameB);

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_input = timing_end-timing_start;

	// ****************************** //
	// Intersection                   //
	// ****************************** //
	timing_start = std::chrono::system_clock::now();

	TriangulationIntersectionVisitor tallyFastVector(std::min(dtA.LengthOrder(),dtB.LengthOrder()),
														6*std::min(dtA.get_nb_of_faces(),dtB.get_nb_of_faces()));

	BuildMeshIntersections(dtA,dtB,tallyFastVector);

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_intersections = timing_end-timing_start;

	timing_start = std::chrono::system_clock::now();

	// Export the meshes
	tallyFastVector.IntersectionTDS_2.ExportMedit(filenameOutput);

	timing_end   		= std::chrono::system_clock::now();
	timing_end_total	= std::chrono::system_clock::now();

	elapsed_seconds_export = timing_end-timing_start;
	elapsed_seconds_total  = timing_end_total - timing_start_total;


	std::cout 	<< dtA.get_nb_of_faces() << " "
				<< dtB.get_nb_of_faces() << " "
				<< tallyFastVector.IntersectionTDS_2.get_nb_of_faces() << " "
				<< elapsed_seconds_input.count() << " "
				<< elapsed_seconds_intersections.count() << " "
				<< elapsed_seconds_export.count() << " "
				<< elapsed_seconds_total.count() << std::endl;

	return 0;
}
