/*
 * main.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"

using namespace CGAL::parameters;

int main(int argc, char *argv[])
{
    Point_3 pA( 0.0, 0.0, 0.0);
    Point_3 qA( 2.0, 2.0, 2.0);

    Point_3 pB( 0.5, 0.5, 0.5);
    Point_3 qB( 1.5, 1.5, 1.5);

    int nA, nB;
	if(argc == 3)
	{
		nA = atoi(argv[1]);
		nB = atoi(argv[2]);
	}

    Triangular_Mesh_3 dt3A, dt3B;

    dt3A.GenerateTestMeshCube(pA,qA,nA,nA,nA);
    dt3B.GenerateTestMeshCube(pB,qB,nB,nB,nB);

    std::ofstream outputF;

    std::vector<Polyhedron> 	outputBrute;
    std::vector<Polyhedron> 	outputFast;

//    BuildMeshIntersections_Brute(dt3A,dt3B,outputBrute);
//
//    outputF.open("3D_mesh_intersections_brute.dat",std::ios::trunc);
//    for(unsigned int iii = 0; iii < outputBrute.size(); ++iii)
//    {
//    	outputF << outputBrute[iii] << std::endl;
//    }
//    outputF.close();

    BuildMeshIntersections(dt3A,dt3B,outputFast);

//    outputF.open("3D_mesh_intersections_fast.dat",std::ios::trunc);
//    for(unsigned int iii = 0; iii < outputFast.size(); ++iii)
//    {
//    	outputF << outputFast[iii] << std::endl;
//    }
//    outputF.close();

//    CGAL::Geomview_stream gv(CGAL::Bbox_3(-1, -1, -1, 2.5, 2.5, 2.5));
//    gv.set_line_width(4);
//
//    gv.set_bg_color(CGAL::Color(0, 200, 200));
//
//    gv << dt3A.mesh;
//    gv << dt3B.mesh;
//
//    std::cout << "Enter a key to finish" << std::endl;
//    std::cin.ignore();

	return 0;
}


