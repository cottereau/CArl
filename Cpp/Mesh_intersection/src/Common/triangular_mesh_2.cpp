/*
 * triangular_mesh_2.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "triangular_mesh_2.h"

std::string Triangular_Mesh_2::get_name()
{
	return mName;
};

void Triangular_Mesh_2::set_indexes()
{
	// Set up indexes for the faces. Each finite face receives a different
	// integer index, while the infinite ones receive the same index, equal to
	// mesh.number_of_faces()

	int dummy = 0;
	int infiniteDummy = mesh.number_of_faces();

	for(All_Face_iterator_2	itFace = mesh.all_faces_begin();
			itFace !=  mesh.all_faces_end(); ++itFace)
	{
		if(!mesh.is_infinite(itFace))
		{
			itFace->info().ExtIndex = dummy;
			++dummy;
		}
		else
		{
			itFace->info().ExtIndex = infiniteDummy;
		}
	}
};

void Triangular_Mesh_2::GenerateTestMeshSquare(const Point_2& initPoint, const Point_2& finalPoint, int nx, int ny)
{
	mesh.clear();

	//double dx = (CGAL::to_double(finalPoint.x()) - CGAL::to_double(initPoint.x()))/(nx-1);
	//double dy = (CGAL::to_double(finalPoint.y()) - CGAL::to_double(initPoint.y()))/(ny-1);
	double dx = (finalPoint.x() - initPoint.x())/(nx-1);
	double dy = (finalPoint.y() - initPoint.y())/(ny-1);

	double dummyX;
	double dummyY;

	boost::random::uniform_real_distribution<double> shift(-1,1);

	// Set boundaries

	Vertex_handle_2 vA = mesh.insert(initPoint);
	Vertex_handle_2 vB = mesh.insert(Point_2(initPoint.x(),finalPoint.y()));
	Vertex_handle_2 vC = mesh.insert(finalPoint);
	Vertex_handle_2 vD = mesh.insert(Point_2(finalPoint.x(),initPoint.y()));

	mesh.insert_constraint(vA,vB);
	mesh.insert_constraint(vB,vC);
	mesh.insert_constraint(vC,vD);
	mesh.insert_constraint(vD,vA);

//	dummyX = (CGAL::to_double(finalPoint.x()) + CGAL::to_double(initPoint.x()))/2 + 0.2*dx*shift(m_rng);
//	dummyY = (CGAL::to_double(finalPoint.y()) + CGAL::to_double(initPoint.y()))/2 + 0.2*dy*shift(m_rng);

	dummyX = (finalPoint.x() + initPoint.x())/2 + 0.2*dx*shift(m_rng);
	dummyY = (finalPoint.y() + initPoint.y())/2 + 0.2*dy*shift(m_rng);

	mesh.insert(Point_2(dummyX,dummyY));
	CGAL::refine_Delaunay_mesh_2(mesh,Criteria(0.125, std::min(dx,dy)));
	set_nb_of_faces();
	set_nb_of_vertices();
};

//void Triangular_Mesh_2::GenerateTestMeshSquare(const Point_2& initPoint, const Point_2& finalPoint, int nx, int ny)
//{
//	mesh.clear();
//
//	double dx = (finalPoint.x() - initPoint.x())/(nx-1);
//	double dy = (finalPoint.y() - initPoint.y())/(ny-1);
//
//	double dummyX;
//	double dummyY;
//
//	boost::random::uniform_real_distribution<double> shift(-1,1);
//
//	// Set points
//
//	// --- BL corner
//	mesh.insert(initPoint);
//
//	// --- Left Border
//	dummyX = initPoint.x();
//	for(int iii = 1; iii < ny-1; ++iii)
//	{
//		dummyY = initPoint.y() + iii*dy + 0.2*dy*shift(m_rng);
//		mesh.insert(Point_2(dummyX,dummyY));
//	}
//
//	// --- TL corner
//	mesh.insert(Point_2(initPoint.x(),finalPoint.y()));
//
//	// --- Top Border
//	dummyY = finalPoint.y();
//	for(int iii = 1; iii < nx-1; ++iii)
//	{
//		dummyX = initPoint.x() + iii*dx + 0.2*dx*shift(m_rng);
//		mesh.insert(Point_2(dummyX,dummyY));
//	}
//
//	// --- TR corner
//	mesh.insert(finalPoint);
//
//	// --- Right border
//	dummyX = finalPoint.x();
//	for(int iii = 1; iii < ny-1; ++iii)
//	{
//		dummyY = initPoint.y() + iii*dy + 0.2*dy*shift(m_rng);
//		mesh.insert(Point_2(dummyX,dummyY));
//	}
//
//	// --- BR corner
//	mesh.insert(Point_2(finalPoint.x(),initPoint.y()));
//	dummyY = initPoint.y();
//	for(int iii = 1; iii < nx-1; ++iii)
//	{
//		dummyX = initPoint.x() + iii*dx + 0.2*dx*shift(m_rng);
//		mesh.insert(Point_2(dummyX,dummyY));
//	}
//
//	// Middle
//	for(int iii = 1; iii < nx-1; ++iii)
//	{
//		for(int jjj = 1; jjj < ny-1; ++jjj)
//		{
//			dummyX = initPoint.x() + iii*dx + 0.2*dx*shift(m_rng);
//			dummyY = initPoint.y() + jjj*dy + 0.2*dy*shift(m_rng);
//			mesh.insert(Point_2(dummyX,dummyY));
//		}
//	}
//
//	set_nb_of_faces();
//	set_nb_of_vertices();
//};

void Triangular_Mesh_2::set_nb_of_faces()
{
	mSize_faces = mesh.number_of_faces();
};

void Triangular_Mesh_2::set_nb_of_vertices()
{
	mSize_vertices = mesh.number_of_vertices();
};

int Triangular_Mesh_2::get_nb_of_faces() const
{
	return mSize_faces;
};

int Triangular_Mesh_2::get_nb_of_vertices() const
{
	return mSize_vertices;
};
