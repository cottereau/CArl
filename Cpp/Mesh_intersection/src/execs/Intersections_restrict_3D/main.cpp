/*
 * main.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"

#include <tuple>

int main(int argc, char *argv[])
{
	// Declare meshes
	// -> BIG mesh
	std::string filenameA = "meshes/3D/test_restriction_cube.msh";
	// -> small mesh
	std::string filenameB = "meshes/3D/test_restriction_micro.msh";
	// -> restrict region
	std::string filenameRestrict = "meshes/3D/restriction_data.txt";

	std::string filenameOutput = "meshes/3D/output/test_restriction_inter_I_3D.msh";

	std::string filenameTableOutputBase = "meshes/equivalence_tables/intersection_carl";

	if(argc == 3)
	{
		filenameA = argv[1];
		filenameB = argv[2];
		filenameOutput = "meshes/3D/output/test_full_inter_I_3D.msh";
		filenameTableOutputBase = "meshes/equivalence_tables/intersection_carl_full";
	}


	if(argc == 4)
	{
		filenameA = argv[1];
		filenameB = argv[2];
		filenameRestrict = argv[3];
	}

	if(argc == 5)
	{
		filenameA = argv[1];
		filenameB = argv[2];
		filenameRestrict = argv[3];
		filenameOutput = argv[4];
	}

	if(argc == 6)
	{
		filenameA = argv[1];
		filenameB = argv[2];
		filenameRestrict = argv[3];
		filenameOutput = argv[4];
		filenameTableOutputBase = argv[5];
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

	std::cout << " ---> Program call : ./intersection_3D_restrict "
			  << filenameA << " "
			  << filenameB << " "
			  << filenameRestrict << " "
			  << filenameOutput << std::endl;

	std::cout << " ---> Mesh files : " << std::endl;
	std::cout << "           (input) " << filenameA << std::endl;
	std::cout << "                   " << filenameB << std::endl;
	std::cout << "                   " << filenameRestrict << std::endl;
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

	std::cout << " ---> Generating the Nef polyhedron ... " << std::endl;
	Nef_Polyhedron couplingRegion;
	if(argc == 3)
	{
		std::cout << " ---> Using the full space" << std::endl;
		couplingRegion.clear(Nef_Polyhedron::COMPLETE);
	}
	else
	{
		GenerateNefRestrictedRegion(filenameRestrict,couplingRegion);
	}
	std::cout << " ---> Generating the Nef polyhedron ... finished!" << std::endl;

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_input = timing_end-timing_start;

	// ****************************** //
	// Intersection                   //
	// ****************************** //
	timing_start = std::chrono::system_clock::now();

	std::cout << " ---> Building intersections ... ";
	std::cout.flush();
	TriangulationIntersectionVisitor_3 	tallyFastVector(std::min(dtA.LengthOrder(),dtB.LengthOrder()),
														dtA.get_nb_of_cells(),dtB.get_nb_of_cells());

	BuildMeshIntersections(dtA,dtB,tallyFastVector,couplingRegion);
	std::cout << "finished!" << std::endl;

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_intersections = timing_end-timing_start;

	timing_start = std::chrono::system_clock::now();

	// Export the meshes
	std::cout << " ---> Exporting data ... ";
	tallyFastVector.IntersectionTDS_3.ExportMesh_3(filenameOutput);
	tallyFastVector.PrintIntersectionTables(filenameTableOutputBase);
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
