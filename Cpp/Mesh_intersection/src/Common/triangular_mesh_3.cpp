/*
 * triangular_mesh_3.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "triangular_mesh_3.h"

std::string Triangular_Mesh_3::get_name()
{
	return mName;
};

void Triangular_Mesh_3::set_indexes()
{
	// Set up indexes for the faces. Each finite face receives a different
	// integer index, while the infinite ones receive the same index, equal to
	// mesh.number_of_faces()

	int dummy = 0;
	int infiniteDummy = mesh.number_of_cells();

	for(All_cells_iterator_3	itCell = mesh.all_cells_begin();
			itCell !=  mesh.all_cells_end(); ++itCell)
	{
		if(!mesh.is_infinite(itCell))
		{
			itCell->info().ExtIndex = dummy;
			++dummy;
		}
		else
		{
			itCell->info().ExtIndex = infiniteDummy;
		}
	}
};

void Triangular_Mesh_3::GenerateTestMeshCube(const Point_3& initPoint, const Point_3& finalPoint, int nx, int ny, int nz, double amplitude)
{
	mesh.clear();

	double dx = (finalPoint.x() - initPoint.x())/(nx-1);
	double dy = (finalPoint.y() - initPoint.y())/(ny-1);
	double dz = (finalPoint.z() - initPoint.z())/(nz-1);

	double dummyX;
	double dummyY;
	double dummyZ;

	boost::random::lagged_fibonacci607 m_rng;
	boost::random::uniform_real_distribution<double> shift(-1,1);

	// --- Set points

	// Set corners
	mesh.insert(initPoint);
	mesh.insert(finalPoint);

	mesh.insert(Point_3(initPoint.x(),initPoint.y(),finalPoint.z()));
	mesh.insert(Point_3(finalPoint.x(),finalPoint.y(),initPoint.z()));

	mesh.insert(Point_3(initPoint.x(),finalPoint.y(),initPoint.z()));
	mesh.insert(Point_3(finalPoint.x(),initPoint.y(),finalPoint.z()));

	mesh.insert(Point_3(initPoint.x(),finalPoint.y(),finalPoint.z()));
	mesh.insert(Point_3(finalPoint.x(),initPoint.y(),initPoint.z()));

	// Set edges
	std::vector<double>		constPosA(4,0);
	std::vector<double>		constPosB(4,0);

	// X
	constPosA[0] = initPoint.y();	constPosB[0] = initPoint.z();
	constPosA[1] = finalPoint.y();	constPosB[1] = initPoint.z();
	constPosA[2] = initPoint.y();	constPosB[2] = finalPoint.z();
	constPosA[3] = finalPoint.y();	constPosB[3] = finalPoint.z();

	for(int iii = 1; iii< nx - 1; ++iii)
	{
		for(int jjj = 0; jjj < 4; ++jjj)
		{
			dummyX = initPoint.x() + iii*dx + amplitude*dx*shift(m_rng);
			mesh.insert(Point_3(dummyX,constPosA[jjj],constPosB[jjj]));
		}
	}

	// Y
	constPosA[0] = initPoint.x();	constPosB[0] = initPoint.z();
	constPosA[1] = finalPoint.x();	constPosB[1] = initPoint.z();
	constPosA[2] = initPoint.x();	constPosB[2] = finalPoint.z();
	constPosA[3] = finalPoint.x();	constPosB[3] = finalPoint.z();

	for(int iii = 1; iii< ny - 1; ++iii)
	{
		for(int jjj = 0; jjj < 4; ++jjj)
		{
			dummyY = initPoint.y() + iii*dy + amplitude*dy*shift(m_rng);
			mesh.insert(Point_3(constPosA[jjj],dummyY,constPosB[jjj]));
		}
	}

	// Z
	constPosA[0] = initPoint.x();	constPosB[0] = initPoint.y();
	constPosA[1] = finalPoint.x();	constPosB[1] = initPoint.y();
	constPosA[2] = initPoint.x();	constPosB[2] = finalPoint.y();
	constPosA[3] = finalPoint.x();	constPosB[3] = finalPoint.y();

	for(int iii = 1; iii< nz - 1; ++iii)
	{
		for(int jjj = 0; jjj < 4; ++jjj)
		{
			dummyZ = initPoint.z() + iii*dz + amplitude*dz*shift(m_rng);
			mesh.insert(Point_3(constPosA[jjj],constPosB[jjj],dummyZ));
		}
	}

	// Set faces
	// XY
	constPosA[0] = initPoint.z();
	constPosA[1] = finalPoint.z();
	for(int iii = 1; iii < nx -1; ++iii)
	{
		for(int jjj = 1; jjj < ny-1; ++jjj)
		{
			for(int kkk = 0; kkk < 2; ++kkk)
			{
				dummyX = initPoint.x() + iii*dx + amplitude*dx*shift(m_rng);
				dummyY = initPoint.y() + jjj*dy + amplitude*dy*shift(m_rng);
				mesh.insert(Point_3(dummyX,dummyY,constPosA[kkk]));
			}
		}
	}

	// XZ
	constPosA[0] = initPoint.y();
	constPosA[1] = finalPoint.y();
	for(int iii = 1; iii < nx-1; ++iii)
	{
		for(int jjj = 1; jjj < nz-1; ++jjj)
		{
			for(int kkk = 0; kkk < 2; ++kkk)
			{
				dummyX = initPoint.x() + iii*dx + amplitude*dx*shift(m_rng);
				dummyZ = initPoint.z() + jjj*dz + amplitude*dz*shift(m_rng);
				mesh.insert(Point_3(dummyX,constPosA[kkk],dummyZ));
			}
		}
	}

	// YZ
	constPosA[0] = initPoint.x();
	constPosA[1] = finalPoint.x();
	for(int iii = 1; iii < ny-1; ++iii)
	{
		for(int jjj = 1; jjj < nz-1; ++jjj)
		{
			for(int kkk = 0; kkk < 2; ++kkk)
			{
				dummyY = initPoint.y() + iii*dy + amplitude*dy*shift(m_rng);
				dummyZ = initPoint.z() + jjj*dz + amplitude*dz*shift(m_rng);
				mesh.insert(Point_3(constPosA[kkk],dummyY,dummyZ));
			}
		}
	}

	// Center
	for(int iii = 1; iii < nx-1; ++iii)
	{
		for(int jjj = 1; jjj < ny-1; ++jjj)
		{
			for(int kkk = 1; kkk < nz-1; ++kkk)
			{
				dummyX = initPoint.x() + iii*dx + amplitude*dx*shift(m_rng);
				dummyY = initPoint.y() + jjj*dy + amplitude*dy*shift(m_rng);
				dummyZ = initPoint.z() + kkk*dy + amplitude*dz*shift(m_rng);
				mesh.insert(Point_3(dummyX,dummyY,dummyZ));
			}
		}
	}

	set_indexes();
	set_nb_of_vertices();
	set_nb_of_facets();
	set_nb_of_cells();
};

void Triangular_Mesh_3::set_nb_of_facets()
{
	mSize_facets = mesh.number_of_facets();
};

void Triangular_Mesh_3::set_nb_of_vertices()
{
	mSize_vertices = mesh.number_of_vertices();
};

void Triangular_Mesh_3::set_nb_of_cells()
{
	mSize_cells = mesh.number_of_cells();
};

int Triangular_Mesh_3::get_nb_of_facets() const
{
	return mSize_facets;
};

int Triangular_Mesh_3::get_nb_of_vertices() const
{
	return mSize_vertices;
};

int Triangular_Mesh_3::get_nb_of_cells() const
{
	return mSize_cells;
};
