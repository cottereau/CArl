/*
 * main.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"

int main(int argc, char *argv[])
{

	std::string meshTestFilename = "data/3D/test_pyramid.mesh";
	std::string meshOutputFilename = "data/3D/test_pyramid_out.mesh";

	if(argc == 3)
	{
		meshTestFilename = argv[1];
		meshOutputFilename = argv[2];
	}

	std::string sillyName = "sillyName";

	Triangular_Mesh_3 meshTest(sillyName,meshTestFilename);
	meshTest.ExportMedit(meshOutputFilename);

	return 0;
}
