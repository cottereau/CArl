/*
 * main.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"


int main(int argc, char *argv[])
{
	// Declare meshes
	std::string filenameA = "data/3D/test_inter_A_3D.mesh";
	std::string filenameB = "data/3D/test_inter_B_3D.mesh";
	std::string filenameOutput = "data/3D/test_inter_I_3D.mesh";

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

	std::cout << " ---> Mesh files : " << std::endl;
	std::cout << "           (input) " << filenameA << std::endl;
	std::cout << "                   " << filenameB << std::endl;
	std::cout << "          (output) " << filenameOutput << std::endl;

	// Import meshes
	timing_start 		= std::chrono::system_clock::now();
	timing_start_total 	= std::chrono::system_clock::now();

	std::cout << " ---> Importing mesh ... ";
	std::string sillyname;
	sillyname = "Testing A";
	Triangular_Mesh_3 dtA(sillyname,filenameA);

	sillyname = "Testing B";
	Triangular_Mesh_3 dtB(sillyname,filenameB);
	std::cout << "finished!" << std::endl;

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_input = timing_end-timing_start;

	// ****************************** //
	// Intersection                   //
	// ****************************** //
	timing_start = std::chrono::system_clock::now();

	std::cout << " ---> Building intersections ... ";
	std::cout.flush();
	TriangulationIntersectionVisitor_3 	tallyFastVector(std::min(dtA.LengthOrder(),dtB.LengthOrder()),
															36*std::min(dtA.get_nb_of_cells(),dtB.get_nb_of_cells()));

	BuildMeshIntersections(dtA,dtB,tallyFastVector);
	std::cout << "finished!" << std::endl;

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_intersections = timing_end-timing_start;

	timing_start = std::chrono::system_clock::now();

	// Export the meshes
	std::cout << " ---> Exporting data ... ";
	tallyFastVector.IntersectionTDS_3.ExportMedit(filenameOutput);
	std::cout << "finished!" << std::endl;

	timing_end   		= std::chrono::system_clock::now();
	timing_end_total	= std::chrono::system_clock::now();

	elapsed_seconds_export = timing_end-timing_start;
	elapsed_seconds_total  = timing_end_total - timing_start_total;

	std::cout << " ---> Timing : " << std::endl;
	std::cout 	<< dtA.get_nb_of_cells() << " "
				<< dtB.get_nb_of_cells() << " "
				<< tallyFastVector.IntersectionTDS_3.get_nb_of_cells() << " "
				<< elapsed_seconds_input.count() << " "
				<< elapsed_seconds_intersections.count() << " "
				<< elapsed_seconds_export.count() << " "
				<< elapsed_seconds_total.count() << std::endl;

	return 0;
}
