/*
 * triangular_mesh_2.h
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef TRIANGULAR_MESH_2_H_
#define TRIANGULAR_MESH_2_H_

#include "common_header.h"
#include "CGAL_typedefs.h"

class	Triangular_Mesh_2
{
private:
	// Members
	std::string			mName;
	int					mSize_vertices;
	int					mSize_faces;

public:
	// Members
	CDT		mesh;

	// Constructors
	Triangular_Mesh_2()
	{
		mSize_vertices = 0;
		mSize_faces = 0;
	}

	Triangular_Mesh_2(std::string &iName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_faces = 0;
	}

	// Methods
	std::string get_name();

	void set_indexes();

	void GenerateTestMeshSquare(const Point_2& initPoint, const Point_2& finalPoint, int nx, int ny);

	void set_nb_of_faces();

	void set_nb_of_vertices();

	int get_nb_of_faces() const;

	int get_nb_of_vertices() const;
};
#endif /* TRIANGULAR_MESH_2_H_ */
