/*
 * main.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"

boost::random::lagged_fibonacci607 m_rng(42);

void ReadInput(Point_2& iA, Point_2& fA, Point_2& iB, Point_2& fB, int& nx, int& ny)
{
	std::cout << "Points of the first square domain: " << std::endl;
	std::cin >> iA >> fA;
	std::cout << "Points of the second square domain: " << std::endl;
	std::cin >> iB >> fB;
	std::cout << "Discretization integers: " << std::endl;
	std::cin >> nx >> ny;
};

void PrintInput(Point_2 iA, Point_2 fA, Point_2 iB, Point_2 fB, int nx, int ny)
{
	std::cout << "Points of the first square domain are " << std::endl;
	std::cout << iA << ", " << fA << std::endl;
	std::cout << "Points of the second square domain are " << std::endl;
	std::cout << iB << ", " << fB << std::endl;
	std::cout << "Discretization integers are " << std::endl;
	std::cout << nx << ", " << ny << std::endl;;
};

int main(int argc, char *argv[])
{

	// ************ //
	// Preamble     //
	// ************ //

	// Declare meshes
	std::string sillyname;
	sillyname = "Testing A";
	Triangular_Mesh_2 dtA(sillyname);

	sillyname = "Testing B";
	Triangular_Mesh_2 dtB(sillyname);

	Point_2 iA, fA, iB, fB;
	int nA, nB;
	if(argc == 11)
	{
		iA = Point_2(atof(argv[1]),atof(argv[2]));
		fA = Point_2(atof(argv[3]),atof(argv[4]));
		iB = Point_2(atof(argv[5]),atof(argv[6]));
		fB = Point_2(atof(argv[7]),atof(argv[8]));

		nA = atoi(argv[9]);
		nB = atoi(argv[10]);
	}
	else
	{
		ReadInput(iA, fA, iB, fB, nA, nB);
	}

	dtA.GenerateTestMeshSquare(iA,fA,nA,nA);
	dtB.GenerateTestMeshSquare(iB,fB,nB,nB);

	// Print triangulations in files
	std::ofstream outputF("testeA.cgmesh",std::ios::trunc);
	outputF << dtA.mesh << std::endl;
	outputF.close();

	outputF.open("testeB.cgmesh",std::ios::trunc);
	outputF << dtB.mesh << std::endl;
	outputF.close();

	// Set indexes
	dtA.set_indexes();
	dtB.set_indexes();

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

	return 0;
}
