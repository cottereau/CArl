/*
 * main.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *      This program restricts a given mesh to a region defined by a Nef
 *      polyhedron.
 */

#include "main.h"

int main(int argc, char *argv[])
{
	// Declare meshes
	std::string filenameMesh;
	std::string filenamePoints;
	std::string filenameOutput;
	std::string filenameTableOutput;

	filenameMesh = "meshes/3D/test_restriction_cube.msh";
	filenamePoints = "meshes/3D/restriction_data.txt";
	filenameOutput = "meshes/3D/output/test_restriction_output.msh";
	filenameTableOutput = "meshes/equivalence_tables/equivalence_carl_restrict_A.dat";

	if(argc == 3)
	{
		filenameMesh = argv[1];
		filenamePoints = argv[1];
		filenameOutput = argv[2];
	}
	else if(argc == 4)
	{
		filenameMesh = argv[1];
		filenamePoints = argv[2];
		filenameOutput = argv[3];
	}
	else if(argc == 5)
	{
		filenameMesh = argv[1];
		filenamePoints = argv[2];
		filenameOutput = argv[3];
		filenameTableOutput = argv[4];
	}
	else
	{
		std::cout << " >>>> IGNORING ARGUMENTS, WILL USE DEFAULT FILES !!! <<<< " << std::endl;
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

	std::cout << " ---> Files : " << std::endl;
	std::cout << "             (input mesh) " << filenameMesh << std::endl;
	std::cout << "      (input restriction) " << filenamePoints << std::endl;
	std::cout << "                 (output) " << filenameOutput << std::endl;
	std::cout << "           (output table) " << filenameTableOutput << std::endl;

	// Import meshes
	timing_start 		= std::chrono::system_clock::now();
	timing_start_total 	= std::chrono::system_clock::now();

	std::cout << " ---> Importing mesh ... " << std::endl;;
	std::string sillyname;
	sillyname = "Testing A";
	Triangular_Mesh_3 dtMesh(sillyname,filenameMesh);
	Triangular_Mesh_3 dtMeshSilly;

	std::cout << " ---> Importing mesh ... finished!" << std::endl;

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_input = timing_end-timing_start;

	// ****************************** //
	// Restriction                    //
	// ****************************** //

	timing_start = std::chrono::system_clock::now();

	std::cout << " ---> Generating the Nef polyhedron ... " << std::endl;
	Nef_Polyhedron couplingRegion;
	GenerateNefRestrictedRegion(filenamePoints,couplingRegion);
	std::cout << " ---> Generating the Nef polyhedron ... finished!" << std::endl;

	std::cout << " ---> Generating the restricted mesh ... " << std::endl;
	dtMesh.RestrictMesh(couplingRegion,dtMeshSilly,filenameTableOutput);
	std::cout << " ---> Generating the restricted mesh ... finished!" << std::endl;

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_intersections = timing_end-timing_start;

	timing_start = std::chrono::system_clock::now();
	dtMeshSilly.ExportMesh_3(filenameOutput);
	timing_end   		= std::chrono::system_clock::now();

	timing_end_total	= std::chrono::system_clock::now();

	elapsed_seconds_export = timing_end-timing_start;
	elapsed_seconds_total  = timing_end_total - timing_start_total;

	std::cout   << elapsed_seconds_input.count() << " "
				<< elapsed_seconds_intersections.count() << " "
				<< elapsed_seconds_export.count() << " "
				<< elapsed_seconds_total.count() << std::endl;

	return 0;
}
