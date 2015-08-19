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

	if(argc == 3)
	{
		filenameA = argv[1];
		filenameB = argv[2];
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

	PolyIntersectionVisitor tallyFastVector;

	BuildMeshIntersections(dtA,dtB,tallyFastVector);

	std::string filename = "intersection_poly_" + std::to_string(dtA.get_nb_of_faces())
	+ "_" + std::to_string(dtB.get_nb_of_faces()) + ".dat";

	outputF.open(filename.c_str(),std::ios::trunc);
	outputF << tallyFastVector << std::endl;
	outputF.close();

	CGAL::Geomview_stream gv(CGAL::Bbox_3(-5,-5,-5,5,5,5));

	gv << dtA.mesh;
	gv << dtB.mesh;

	std::cout << "Enter a key to finish" << std::endl;
	std::cin.ignore();

	return 0;
}
