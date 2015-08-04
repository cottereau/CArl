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
    Point_3 qA( 1.0, 1.0, 1.0);

    Point_3 pB( 0.5, 0.5, 0.5);
    Point_3 qB( 1.5, 1.5, 1.5);

    int nx = 2; int ny = 2; int nz = 2;

    Triangular_Mesh_3 dt3A, dt3B;

    dt3A.GenerateTestMeshCube(pA,qA,nx,ny,nz);
    dt3B.GenerateTestMeshCube(pB,qB,nx,ny,nz);

    std::vector<Polyhedron> 	output;

    BuildMeshIntersections_Brute(dt3A,dt3B,output);

    std::cout << output.size() << std::endl;
    std::cout << "-----------" << std::endl;

    for(unsigned int iii = 0; iii < output.size(); ++iii)
    {
    	std::cout << output[iii] << std::endl;
    }
    CGAL::Geomview_stream gv(CGAL::Bbox_3(-1, -1, -1, 2.5, 2.5, 2.5));
    gv.set_line_width(4);

    gv.set_bg_color(CGAL::Color(0, 200, 200));

    gv << dt3A.mesh;
    gv << dt3B.mesh;

    std::cout << "Enter a key to finish" << std::endl;
    std::cin.ignore();

	return 0;
}


