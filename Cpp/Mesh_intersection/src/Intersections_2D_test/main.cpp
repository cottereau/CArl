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
	std::string filenameOutput = "data/test_mesh_intersection.msh";

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

	std::string sillyname;
	sillyname = "Testing A";
	Triangular_Mesh_2 dtA(sillyname,filenameA);

	sillyname = "Testing B";
	Triangular_Mesh_2 dtB(sillyname,filenameB);

	// Print triangulations in files
	std::ofstream outputF("testeA.cgmesh",std::ios::trunc);
	outputF << dtA.mesh << std::endl;
	outputF.close();

	outputF.open("testeB.cgmesh",std::ios::trunc);
	outputF << dtB.mesh << std::endl;
	outputF.close();

	// ****************************** //
	// Intersection                   //
	// ****************************** //

	TriangulationIntersectionVisitor tallyFastVector(pow(10,-2),6*std::min(dtA.get_nb_of_faces(),dtB.get_nb_of_faces()));

	BuildMeshIntersections(dtA,dtB,tallyFastVector);

	tallyFastVector.IntersectionTDS_2.ExportGmsh(filenameOutput);

	return 0;
}
