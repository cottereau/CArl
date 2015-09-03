/*
 * main.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"

int main(int argc, char *argv[])
{
	std::string meshTestFilename = "data/test_mesh.mesh";
	std::string sillyName = "sillyName";
	Triangular_Mesh_2 meshTest(sillyName,meshTestFilename);

	std::string meshOutputFilename = "data/test_mesh_out.msh";
	meshTest.ExportGmsh(meshOutputFilename);
	meshOutputFilename = "data/test_mesh_out.mesh";
	meshTest.ExportMedit(meshOutputFilename);

	return 0;
}
