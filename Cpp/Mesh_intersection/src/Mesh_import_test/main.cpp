/*
 * main.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"

int main(int argc, char *argv[])
{
	std::string meshTestFilename = "data/test_mesh.msh";

	Triangular_Mesh_2 meshTest;

	meshTest.importGmsh(meshTestFilename);

	return 0;
}


