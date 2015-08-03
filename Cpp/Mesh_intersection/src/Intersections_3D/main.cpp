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
    Point_3 p( 0.0, 0.0, 0.0);
    Point_3 q( 1.0, 1.0, 1.0);

    int nx = 4; int ny = 4; int nz = 4;

    Triangular_Mesh_3 biba;

    biba.GenerateTestMeshCube(p,q,nx,ny,nz);

    CGAL::Geomview_stream gv(CGAL::Bbox_3(-1, -1, -1, 2, 2, 2));
    gv.set_line_width(4);

    gv.set_bg_color(CGAL::Color(0, 200, 200));

    gv << biba.mesh;

    std::cout << "Enter a key to finish" << std::endl;
    char ch;
    std::cin >> ch;

	return 0;
}


