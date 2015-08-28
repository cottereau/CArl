/*
 * triangular_mesh_intersection_2.cpp
 *
 *  Created on: Aug 25, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "triangular_mesh_intersection_2.h"


void Intersection_Mesh_2::InitializeIntersection(int& iPreallocation)
{
	Initialize();
	mInterVertexDummyIndex = 1;
	mInterFaceDummyIndex = 1;
	mInterVertexHandleIndexMap.reserve(iPreallocation);
	mInterVertexIncidentFaces.reserve(iPreallocation);

	mInterVertexHandle.resize(6);
	mInterVertexIndex.resize(6);
}

void Intersection_Mesh_2::AddPolygon(const Triangle_2& t)
{
	Point_2 dummyPoint;
	std::pair<int, Face_handle_2> dummyPair;

	for(int iii = 0; iii < 3; ++iii)
	{
		dummyPoint = t.vertex(iii);

		UpdateBbox(dummyPoint);
		mInterVertexHandle[iii] = Create_Vertex_2(dummyPoint,mInterVertexDummyIndex);
		mVertexHandleIndexMap[mInterVertexDummyIndex]=mInterVertexHandle[iii];
		mInterVertexIndex[iii] = mInterVertexDummyIndex;
		++mInterVertexDummyIndex;
	}

	Face_handle_2 dummyFace = Create_Face_2(mInterVertexHandle[2],mInterVertexHandle[1],mInterVertexHandle[0],mInterFaceDummyIndex);
	for(int iii = 0; iii < 3; ++iii)
	{
		dummyPair = std::pair<int, Face_handle_2>(mInterVertexIndex[iii],dummyFace);
		mInterVertexIncidentFaces.insert(dummyPair);
	}

	++mInterFaceDummyIndex;
};

void Intersection_Mesh_2::AddPolygon(Polygon_2& t,int nbOfVertices, double CharacteristicArea)
{
	Triangle_2	dummyTriangle;

	Point_2 dummyPoint;
	Face_handle_2 dummyFace;
	std::pair<int, Face_handle_2> dummyPair;

	for(int iii = 0; iii < nbOfVertices; ++iii)
	{
		dummyPoint = t.vertex(iii);

		UpdateBbox(dummyPoint);
		mInterVertexHandle[iii] = Create_Vertex_2(dummyPoint,mInterVertexDummyIndex);
		mVertexHandleIndexMap[mInterVertexDummyIndex]=mInterVertexHandle[iii];
		mInterVertexIndex[iii] = mInterVertexDummyIndex;
		++mInterVertexDummyIndex;
	}

	for(int iii = 2; iii < nbOfVertices; ++iii)
	{
		dummyTriangle = Triangle_2(mInterVertexHandle[0]->point(),mInterVertexHandle[iii-1]->point(),mInterVertexHandle[iii]->point());
		if(std::abs(dummyTriangle.area()) > CharacteristicArea)
		{
			dummyFace = Create_Face_2(mInterVertexHandle[iii],mInterVertexHandle[iii-1],mInterVertexHandle[0],mInterFaceDummyIndex);
			dummyPair = std::pair<int, Face_handle_2>(mInterVertexIndex[0],dummyFace);
			mInterVertexIncidentFaces.insert(dummyPair);
			dummyPair = std::pair<int, Face_handle_2>(mInterVertexIndex[iii-1],dummyFace);
			mInterVertexIncidentFaces.insert(dummyPair);
			dummyPair = std::pair<int, Face_handle_2>(mInterVertexIndex[iii],dummyFace);
			mInterVertexIncidentFaces.insert(dummyPair);

			++mInterFaceDummyIndex;
		}
	}
};

void Intersection_Mesh_2::CleanUp(double CharacteristicArea)
{
	mEps = LengthOrder()*pow(10,-6);
	mGridNx = XLength()/mEps + 1;
	mGridNy = YLength()/mEps + 1;

	// For all the vertices, test if they are valid
	for(Finite_vertices_iterator_2 	itVertex = mesh.finite_vertices_begin();
									itVertex != mesh.finite_vertices_end();
									++itVertex)
	{
		TestVertex(itVertex);
	}

	set_nb_of_faces();
	set_nb_of_vertices();
}

void Intersection_Mesh_2::TestVertex(Vertex_handle_2 itVertex)
{
	long int index = ConvertToIndex(itVertex->point());
	int triangleVertexIdx = -1;
	std::unordered_map<long int,Vertex_handle_2>::const_iterator searchVertex;
	searchVertex = mInterVertexHandleIndexMap.find(index);

	if(searchVertex == mInterVertexHandleIndexMap.end())
	{
		// Vertex not present, insert it in the mapping!
		mInterVertexHandleIndexMap[index] = itVertex;
	}
	else
	{
		// Vertex is already in the map - must remove it!
		faceRangeIteratorPair FaceRange = mInterVertexIncidentFaces.equal_range(itVertex->info().ExtIndex);
		faceRangeIterator	  itRange = FaceRange.first;

		for( ; itRange != FaceRange.second;++itRange)
		{
			// Exchange vertex on the first triangle
			triangleVertexIdx = itRange->second->index(itVertex);
			itRange->second->set_vertex(triangleVertexIdx,searchVertex->second);
			searchVertex->second->set_face(itRange->second);
		}

		mesh.delete_vertex(itVertex);
	}
}

long int Intersection_Mesh_2::ConvertToIndex(Point_2 iPoint)
{
	return lround((iPoint.x() - mPointMin.x())/mEps)*mGridNy + lround((iPoint.y() - mPointMin.y())/mEps);
}

