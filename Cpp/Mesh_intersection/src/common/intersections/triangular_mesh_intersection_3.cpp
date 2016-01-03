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

	mInterVertexHandle.resize(36);
	mInterVertexIndex.resize(36);
}

void Intersection_Mesh_3::AddPolyhedron(const Tetrahedron& t)
{
//	Point_3 dummyPoint;
	std::pair<int, Cell_handle_3> dummyPair;

	for(int iii = 0; iii < 4; ++iii)
	{
		AddPolyhedronVertex(iii,t);
	}

	Cell_handle_3 dummyCell = Create_Cell_3(
									mInterVertexHandle[3],
									mInterVertexHandle[2],
									mInterVertexHandle[1],
									mInterVertexHandle[0],
									mInterCellDummyIndex);
	for(int iii = 0; iii < 4; ++iii)
	{
		dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[iii],dummyCell);
		mInterVertexIncidentCells.insert(dummyPair);
	}

	++mInterCellDummyIndex;
};

void Intersection_Mesh_3::AddPolyhedronVertex(int iii, const Tetrahedron& t)
{
	Point_3 dummyPoint;
	dummyPoint = t.vertex(iii);

	UpdateBbox(dummyPoint);
	mInterVertexHandle[iii] = Create_Vertex_3(dummyPoint,mInterVertexDummyIndex);
	mVertexHandleIndexMap[mInterVertexDummyIndex]=mInterVertexHandle[iii];
	mInterVertexIndex[iii] = mInterVertexDummyIndex;
	++mInterVertexDummyIndex;
}

void Intersection_Mesh_3::AddPolyhedronVertex(int iii, Point_3& dummyPoint)
{
	UpdateBbox(dummyPoint);
	mInterVertexHandle[iii] = Create_Vertex_3(dummyPoint,mInterVertexDummyIndex);
	mVertexHandleIndexMap[mInterVertexDummyIndex]=mInterVertexHandle[iii];
	mInterVertexIndex[iii] = mInterVertexDummyIndex;
	++mInterVertexDummyIndex;
}

bool Intersection_Mesh_3::AddPolyhedron(Polyhedron& t, double dummyVolume, int intersectionIdx)
{
	polyTriang.ConvertToTriangulation_3(t);
	Cell_handle_3 dummyCell;
	Vertex_handle_3 bufferHandle0, bufferHandle1, bufferHandle2, bufferHandle3;

	std::pair<int, Cell_handle_3> dummyPair;
	bool bInsertedATetrahedron = false;

	// First, test if the cells should be added
	for(Finite_cells_iterator_3 itSubCells = polyTriang.mesh.finite_cells_begin();
								itSubCells != polyTriang.mesh.finite_cells_end();
								++itSubCells)
	{
		dummyTetrahedron = Tetrahedron(	itSubCells->vertex(0)->point(),
										itSubCells->vertex(1)->point(),
										itSubCells->vertex(2)->point(),
										itSubCells->vertex(3)->point());
		if(std::abs(dummyTetrahedron.volume()) > dummyVolume)
		{
			bInsertedATetrahedron = true;
			itSubCells->info().ToAdd = true;
			for(int iii = 0; iii < 4; ++iii)
			{
				itSubCells->vertex(iii)->info().ToAdd = true;
			}
		}
	}

	if(bInsertedATetrahedron)
	{
		// Insert the vertices
		int dummyIndex = 0;
		for(Finite_vertices_iterator_3 itSubVertices = polyTriang.mesh.finite_vertices_begin();
									itSubVertices != polyTriang.mesh.finite_vertices_end();
									++itSubVertices)
		{
			if(itSubVertices->info().ToAdd==true)
			{
				AddPolyhedronVertex(dummyIndex,itSubVertices->point());
				itSubVertices->info().ExtIndex = dummyIndex;
				++dummyIndex;
			}
		}

		// Insert the cells
		for(Finite_cells_iterator_3 itSubCells = polyTriang.mesh.finite_cells_begin();
									itSubCells != polyTriang.mesh.finite_cells_end();
									++itSubCells)
		{
			if(itSubCells->info().ToAdd==true)
			{
				bufferHandle0 = mInterVertexHandle[itSubCells->vertex(0)->info().ExtIndex];
				bufferHandle1 = mInterVertexHandle[itSubCells->vertex(1)->info().ExtIndex];
				bufferHandle2 = mInterVertexHandle[itSubCells->vertex(2)->info().ExtIndex];
				bufferHandle3 = mInterVertexHandle[itSubCells->vertex(3)->info().ExtIndex];

				// Create the cell
				dummyCell = Create_Cell_3(bufferHandle0,bufferHandle1,bufferHandle2,bufferHandle3,mInterCellDummyIndex);
				dummyCell->info().IntersectionIndex = intersectionIdx;

				// Set up neighbours
				dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[itSubCells->vertex(0)->info().ExtIndex],dummyCell);
				mInterVertexIncidentCells.insert(dummyPair);
				dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[itSubCells->vertex(1)->info().ExtIndex],dummyCell);
				mInterVertexIncidentCells.insert(dummyPair);
				dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[itSubCells->vertex(2)->info().ExtIndex],dummyCell);
				mInterVertexIncidentCells.insert(dummyPair);
				dummyPair = std::pair<int, Cell_handle_3>(mInterVertexIndex[itSubCells->vertex(3)->info().ExtIndex],dummyCell);
				mInterVertexIncidentCells.insert(dummyPair);

				++mInterCellDummyIndex;
			}
		}
	}
	return bInsertedATetrahedron;
};

void Intersection_Mesh_3::CleanUp()
{
	set_nb_of_cells();
	set_nb_of_vertices();

	mEps = LengthOrder()*pow(10,-6);
	mGridNx = XLength()/mEps + 1;
	mGridNy = YLength()/mEps + 1;
	mGridNz = ZLength()/mEps + 1;

	mInterVertexToRemoveMap.reserve(mSize_vertices);

	// For all the vertices, test if they are valid
	for(Finite_vertices_iterator_3 	itVertex = mesh.finite_vertices_begin();
									itVertex != mesh.finite_vertices_end();
									++itVertex)
	{
		TestVertex(itVertex);
	}

	Vertex_handle_3 removableHandle;
	for(std::unordered_set<int>::iterator itRemove = mInterVertexToRemoveMap.begin();
			itRemove != mInterVertexToRemoveMap.end();
			++itRemove)
	{
		removableHandle = mVertexHandleIndexMap[*itRemove];
		mVertexHandleIndexMap.erase(*itRemove);
		mesh.tds().delete_vertex(removableHandle);
	}

	set_nb_of_cells();
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

	// Set data structures
	set_indexes();

//	// Connect tetrahedrons
//	ConnectTetrahedrons_3();
//
//	// Create and connect the INFINITE cells
//	AddInfiniteCells_3();
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
		cellRangeIteratorPair CellRange = mInterVertexIncidentCells.equal_range(itVertex->info().ExtIndex);
		cellRangeIterator	  itRange = CellRange.first;

		for( ; itRange != CellRange.second;++itRange)
		{
			// Exchange vertex on the first triangle
			triangleVertexIdx = itRange->second->index(itVertex);
			itRange->second->set_vertex(triangleVertexIdx,searchVertex->second);
			searchVertex->second->set_cell(itRange->second);
		}

		// Mark for removal!
		mInterVertexToRemoveMap.insert(itVertex->info().ExtIndex);
	}
}

long int Intersection_Mesh_3::ConvertToIndex(Point_3 iPoint)
{
	return lround((iPoint.x() - mPointMin.x())/mEps)*mGridNy*mGridNz
			+ lround((iPoint.y() - mPointMin.y())/mEps)*mGridNz
			+ lround((iPoint.z() - mPointMin.z())/mEps);
}
