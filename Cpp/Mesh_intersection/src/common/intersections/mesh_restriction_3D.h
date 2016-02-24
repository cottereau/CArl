/*
 * mesh_restriction_3D.h
 *
 *  Created on: Jan 3, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_MESH_RESTRICTION_MESH_RESTRICTION_3D_H_
#define COMMON_MESH_RESTRICTION_MESH_RESTRICTION_3D_H_

#include "CGAL_typedefs.h"
#include "common_header.h"
#include "triangular_mesh_3.h"
#include "intersection_functions_3.h"

// Inexact to exact converter
Kernel_to_ExactKernel ConvertInexactToExact;

void GenerateNefRestrictedRegion(std::string& filenamePoints, Nef_Polyhedron& outputNef)
{
	// Set fstream and buffer variables
	std::ifstream filePoints(filenamePoints);
	assert(filePoints.good());

	// Buffer string
	std::string bufferLine;
	std::stringstream	dataBuffer;

	// Clean up input Nef polyhedron
	outputNef.clear();

	// And a "holes" Nef polyhedron
	Nef_Polyhedron holesNef;

	// Map and vector to convert points list. The Gmsh indexes are the keys
	// of the map.
	std::map<int,int> exactPointsIndexMap;
	std::vector<ExactPoint_3> exactPointsVector;

	// Declare a few boolean control variables
	bool bHasNodes = false;
	bool bHasOuterShell = false;
	bool bHasHoles = false;

	bool bNodes = false;
	bool bOuterShell = false;
	bool bHoles = false;

	// Read info until the file ends
	while(std::getline(filePoints,bufferLine))
	{
		// Set the booleans
		bNodes = bufferLine.find("$RestrictionNodes")!=std::string::npos;
		bOuterShell = bufferLine.find("$RestrictionOuterShell")!=std::string::npos;
		bHoles = bufferLine.find("$RestrictionHoles")!=std::string::npos;

		if(bNodes)
		{
			bHasNodes = true;
		}

		if(bOuterShell)
		{
			bHasOuterShell = true;
		}

		if(bHoles)
		{
			bHasHoles = true;
		}

		// Search for the nodes tag
		if(bNodes)
		{
			// Declare and read the variables that we'll need
			unsigned int nbOfPoints;
			int bufferNodeIndex;
			double bufferX = 0, bufferY = 0, bufferZ = 0;
			Point_3 dummyInexactPoint;

			filePoints >> nbOfPoints;
			exactPointsVector.resize(nbOfPoints);

			std::getline(filePoints,bufferLine);
			for(unsigned int iii = 0; iii < nbOfPoints; ++iii)
			{
				std::getline(filePoints,bufferLine);

				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				// The "nodes" follow the same structure as Gmsh's element nodes
				dataBuffer >> bufferNodeIndex >> bufferX >> bufferY >> bufferZ;
				dummyInexactPoint = Point_3(bufferX,bufferY,bufferZ);
				exactPointsVector[iii] = ConvertInexactToExact(dummyInexactPoint);

				exactPointsIndexMap.insert(std::pair<int,int>(bufferNodeIndex,iii));
			}
		}

		// Search for an outer shell tag
		if(bOuterShell)
		{
			// Declare and read the variables that we'll need
			unsigned int nbOfOuterShells;
			unsigned int dummyShellIndex;
			unsigned int nbOfPointsInShell;
			unsigned int bufferPointIdx;
			std::vector<ExactPoint_3> nodeList(4);
			ExactPolyhedron dummyExactPolyhedron;
			Nef_Polyhedron dummyNefPolyhedron;

			filePoints >> nbOfOuterShells;

			std::getline(filePoints,bufferLine);
			for(unsigned int iii = 0; iii < nbOfOuterShells; ++iii)
			{
				std::getline(filePoints,bufferLine);

				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				/*
				 * 	The shell lines follow the structure
				 *
				 * 		[dummyIndex] [n] [pointIdx_1] ... [pointIdx_n]
				 *
				 * 	Note that the point index starts at 1, following Gmsh
				 * 	notation!
				 */

				dataBuffer >> dummyShellIndex >> nbOfPointsInShell;
				nodeList.resize(nbOfPointsInShell);

				for(unsigned int jjj = 0; jjj < nbOfPointsInShell; ++jjj)
				{
					dataBuffer >> bufferPointIdx;
					nodeList[jjj] = exactPointsVector[exactPointsIndexMap[bufferPointIdx]];
				}

				// Now build the exact polyhedron and add it to the Nef polyhedron
				CGAL::convex_hull_3(	nodeList.begin(),nodeList.end(),
										dummyExactPolyhedron);

				dummyNefPolyhedron = Nef_Polyhedron(dummyExactPolyhedron);
				outputNef += dummyNefPolyhedron;
			}
		}

		// Search for an holes shell tag
		if(bHoles)
		{
			// Declare and read the variables that we'll need
			unsigned int nbOfHoles;
			unsigned int dummyHoleIndex;
			unsigned int nbOfPointsInHole;
			unsigned int bufferPointIdx;
			std::vector<ExactPoint_3> nodeList(4);
			ExactPolyhedron dummyExactPolyhedron;
			Nef_Polyhedron dummyNefPolyhedron;

			filePoints >> nbOfHoles;

			std::getline(filePoints,bufferLine);
			for(unsigned int iii = 0; iii < nbOfHoles; ++iii)
			{
				std::getline(filePoints,bufferLine);

				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				/*
				 * 	The Hole lines follow the structure
				 *
				 * 		[dummyIndex] [n] [pointIdx_1] ... [pointIdx_n]
				 *
				 * 	Note that the point index starts at 1, following Gmsh
				 * 	notation!
				 */

				dataBuffer >> dummyHoleIndex >> nbOfPointsInHole;
				nodeList.resize(nbOfPointsInHole);

				for(unsigned int jjj = 0; jjj < nbOfPointsInHole; ++jjj)
				{
					dataBuffer >> bufferPointIdx;
					nodeList[jjj] = exactPointsVector[exactPointsIndexMap[bufferPointIdx]];
				}

				// Now build the exact polyhedron and add it to the Nef polyhedron
				CGAL::convex_hull_3(	nodeList.begin(),nodeList.end(),
										dummyExactPolyhedron);

				dummyNefPolyhedron = Nef_Polyhedron(dummyExactPolyhedron);
				holesNef += dummyNefPolyhedron;
			}
		}
	}

	homemade_assert_msg(bHasNodes,"No nodes section!");
	homemade_assert_msg(bHasOuterShell,"No shell section!");
	// The holes are optional

	// Finally: compute (outer_shell - holes), if needed
	if(bHasHoles)
	{
		outputNef -= holesNef;
	}
};

#endif /* COMMON_MESH_RESTRICTION_MESH_RESTRICTION_3D_H_ */
