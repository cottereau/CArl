/*
 * triangular_mesh_intersection_3.cpp
 *
 *  Created on: Aug 25, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "triangular_mesh_intersection_3.h"


void Intersection_Mesh_3::InitializeIntersection(int& iPreallocation)
{
	Initialize();
	mInterVertexDummyIndex = 1;
	mInterCellDummyIndex = 1;
	mInterVertexHandleIndexMap.reserve(iPreallocation);
	mInterVertexIncidentCells.reserve(iPreallocation);

	mInterVertexHandle.resize(6);
	mInterVertexIndex.resize(6);
}

void Intersection_Mesh_3::AddPolyhedron(const Triangle_3& t)
{
//	Point_3 dummyPoint;
	std::pair<int, Cell_handle_3> dummyPair;

	for(int iii = 0; iii < 3; ++iii)
	{
		AddPolyhedronVertex(iii,t);
//		dummyPoint = t.vertex(iii);
//
//		UpdateBbox(dummyPoint);
//		mInterVertexHandle[iii] = Create_Vertex_3(dummyPoint,mInterVertexDummyIndex);
//		mVertexHandleIndexMap[mInterVertexDummyIndex]=mInterVertexHandle[iii];
//		mInterVertexIndex[iii] = mInterVertexDummyIndex;
//		++mInterVertexDummyIndex;
	}

	Cell_handle_3 dummyCell = Create_Cell_3(mInterVertexHandle[2],mInterVertexHandle[1],mInterVertexHandle[0],mInterCellDummyIndex);
	for(int iii = 0; iii < 3; ++iii)
	{
		dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[iii],dummyCell);
		mInterVertexIncidentCells.insert(dummyPair);
	}

	++mInterCellDummyIndex;
};

template<typename GeometryType>
void Intersection_Mesh_3::AddPolyhedronVertex(int iii, GeometryType& t)
{
	Point_3 dummyPoint;
	dummyPoint = t.vertex(iii);

	UpdateBbox(dummyPoint);
	mInterVertexHandle[iii] = Create_Vertex_3(dummyPoint,mInterVertexDummyIndex);
	mVertexHandleIndexMap[mInterVertexDummyIndex]=mInterVertexHandle[iii];
	mInterVertexIndex[iii] = mInterVertexDummyIndex;
	++mInterVertexDummyIndex;
}

void Intersection_Mesh_3::AddPolyhedron(Polyhedron& t,int nbOfVertices, double CharacteristicArea)
{
	Triangle_3	dummyTriangle;


	Cell_handle_3 dummyCell;
	std::pair<int, Cell_handle_3> dummyPair;
	bool createdPrevious = false;

	// Add first vertex
	AddPolyhedronVertex(0,t);

	for(int iii = 2; iii < nbOfVertices; ++iii)
	{
		dummyTriangle = Triangle_3(mInterVertexHandle[0]->point(),t.vertex(iii-1),t.vertex(iii));
		if(std::abs(dummyTriangle.area()) > CharacteristicArea)
		{
			// Create the vertices
			if(!createdPrevious)
			{
				AddPolyhedronVertex(iii-1,t);
			}
			AddPolyhedronVertex(iii,t);

			// Create face
			dummyCell = Create_Cell_3(mInterVertexHandle[iii],mInterVertexHandle[iii-1],mInterVertexHandle[0],mInterCellDummyIndex);
						dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[0],dummyCell);
						mInterVertexIncidentCells.insert(dummyPair);
						dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[iii-1],dummyCell);
						mInterVertexIncidentCells.insert(dummyPair);
						dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[iii],dummyCell);
						mInterVertexIncidentCells.insert(dummyPair);

			++mInterCellDummyIndex;

			createdPrevious = true;
		}
		else
		{
			createdPrevious = false;
		}
	}

//	for(int iii = 0; iii < nbOfVertices; ++iii)
//	{
//		dummyPoint = t.vertex(iii);
//
//		UpdateBbox(dummyPoint);
//		mInterVertexHandle[iii] = Create_Vertex_3(dummyPoint,mInterVertexDummyIndex);
//		mVertexHandleIndexMap[mInterVertexDummyIndex]=mInterVertexHandle[iii];
//		mInterVertexIndex[iii] = mInterVertexDummyIndex;
//		++mInterVertexDummyIndex;
//	}
//
//	for(int iii = 2; iii < nbOfVertices; ++iii)
//	{
//		dummyTriangle = Triangle_3(mInterVertexHandle[0]->point(),mInterVertexHandle[iii-1]->point(),mInterVertexHandle[iii]->point());
//		if(std::abs(dummyTriangle.area()) > CharacteristicArea)
//		{
//			dummyCell = Create_Cell_3(mInterVertexHandle[iii],mInterVertexHandle[iii-1],mInterVertexHandle[0],mInterCellDummyIndex);
//			dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[0],dummyCell);
//			mInterVertexIncidentCells.insert(dummyPair);
//			dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[iii-1],dummyCell);
//			mInterVertexIncidentCells.insert(dummyPair);
//			dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[iii],dummyCell);
//			mInterVertexIncidentCells.insert(dummyPair);
//
//			++mInterCellDummyIndex;
//		}
//	}
};

void Intersection_Mesh_3::CleanUp()
{
	set_nb_of_faces();
	set_nb_of_vertices();

	mEps = LengthOrder()*pow(10,-6);
	mGridNx = XLength()/mEps + 1;
	mGridNy = YLength()/mEps + 1;

	// For all the vertices, test if they are valid
	for(Finite_vertices_iterator_3 	itVertex = mesh.finite_vertices_begin();
									itVertex != mesh.finite_vertices_end();
									++itVertex)
	{
		TestVertex(itVertex);
	}

	set_nb_of_faces();
	set_nb_of_vertices();

	// Now, reorder them
	int dummyIdx = 1;
	for(Finite_vertices_iterator_3 	itVertex = mesh.finite_vertices_begin();
									itVertex != mesh.finite_vertices_end();
									++itVertex)
	{
		itVertex->info().ExtIndex = dummyIdx;
		++dummyIdx;
	}
}

void Intersection_Mesh_3::TestVertex(Vertex_handle_3 itVertex)
{
	long int index = ConvertToIndex(itVertex->point());
	int triangleVertexIdx = -1;
	std::unordered_map<long int,Vertex_handle_3>::const_iterator searchVertex;
	searchVertex = mInterVertexHandleIndexMap.find(index);

	if(searchVertex == mInterVertexHandleIndexMap.end())
	{
		// Vertex not present, insert it in the mapping!
		mInterVertexHandleIndexMap[index] = itVertex;
	}
	else
	{
		// Vertex is already in the map - must remove it!
		faceRangeIteratorPair CellRange = mInterVertexIncidentCells.equal_range(itVertex->info().ExtIndex);
		faceRangeIterator	  itRange = CellRange.first;

		for( ; itRange != CellRange.second;++itRange)
		{
			// Exchange vertex on the first triangle
			triangleVertexIdx = itRange->second->index(itVertex);
			itRange->second->set_vertex(triangleVertexIdx,searchVertex->second);
			searchVertex->second->set_face(itRange->second);
		}

		// Remove it!
		mVertexHandleIndexMap.erase(itVertex->info().ExtIndex);
		mesh.delete_vertex(itVertex);
	}
}

long int Intersection_Mesh_3::ConvertToIndex(Point_3 iPoint)
{
	return lround((iPoint.x() - mPointMin.x())/mEps)*mGridNy + lround((iPoint.y() - mPointMin.y())/mEps);
}

template
void Intersection_Mesh_3::AddPolyhedronVertex<Polyhedron>(int iii, Polyhedron& t);

template
void Intersection_Mesh_3::AddPolyhedronVertex<Triangle_3>(int iii, Triangle_3& t);
