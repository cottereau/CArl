/*
 * triangular_mesh_3.h
 *
 *  Created on: Jul 31, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef TRIANGULAR_MESH_3_H_
#define TRIANGULAR_MESH_3_H_

#include "common_header.h"
#include "CGAL_typedefs.h"
#include "defines_3D.h"

class	Triangular_Mesh_3
{
protected:
	// Members
	std::string			mName;
	int					mSize_vertices;
	int					mSize_facets;
	int					mSize_cells;

public:
	// Members
	DT_3		mesh;

	// Constructors
	Triangular_Mesh_3()
	{
		mSize_vertices = 0;
		mSize_facets = 0;
		mSize_cells = 0;
	}

	Triangular_Mesh_3(std::string &iName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_facets = 0;
		mSize_cells = 0;
	}

	// Methods
	std::string get_name();

	void set_indexes();

	void GenerateTestMeshCube(const Point_3& initPoint, const Point_3& finalPoint, int nx, int ny, int nz, double amplitude = 0.15);

	void set_nb_of_cells();

	void set_nb_of_facets();

	void set_nb_of_vertices();

	int get_nb_of_cells() const;

	int get_nb_of_facets() const;

	int get_nb_of_vertices() const;
};
#endif /* TRIANGULAR_MESH_3_H_ */
