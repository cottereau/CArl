/*
 * triangular_mesh_3.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *		Methods of the class "Triangular_Mesh_3". A detailed description of each
 *		method can be found inside the header file.
 *
 */

#include "triangular_mesh_3.h"
// PRIVATE Methods

Vertex_handle_3 Triangular_Mesh_3::Create_Vertex_3(Point_3& pointInput, int indexInput)
{
	Vertex_handle_3 outputHandle;
	outputHandle = mesh.tds().create_vertex();
	outputHandle->set_point(pointInput);
	outputHandle->info().ExtIndex = indexInput;
	return outputHandle;
}

Cell_handle_3 Triangular_Mesh_3::Create_Cell_3(int i0, int i1, int i2, int i3, int idx)
{
	// Create cell
	Cell_handle_3 outputHandle;
	outputHandle = mesh.tds().create_cell(	mVertexHandleIndexMap[i0],
											mVertexHandleIndexMap[i1],
											mVertexHandleIndexMap[i2],
											mVertexHandleIndexMap[i3]);

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(4,0);

	// Set incident cells
	mVertexHandleIndexMap[i0]->set_cell(outputHandle);
	mVertexHandleIndexMap[i1]->set_cell(outputHandle);
	mVertexHandleIndexMap[i2]->set_cell(outputHandle);
	mVertexHandleIndexMap[i3]->set_cell(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;

	return outputHandle;
}

Cell_handle_3 Triangular_Mesh_3::Create_Cell_3(std::vector<int>& vertices, int idx)
{
	// Create cell
	Cell_handle_3 outputHandle;
	outputHandle = mesh.tds().create_cell(
								mVertexHandleIndexMap[vertices[0]],
								mVertexHandleIndexMap[vertices[1]],
								mVertexHandleIndexMap[vertices[2]],
								mVertexHandleIndexMap[vertices[3]]);

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(4,0);

	// Set incident cells
	mVertexHandleIndexMap[vertices[0]]->set_cell(outputHandle);
	mVertexHandleIndexMap[vertices[1]]->set_cell(outputHandle);
	mVertexHandleIndexMap[vertices[2]]->set_cell(outputHandle);
	mVertexHandleIndexMap[vertices[3]]->set_cell(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;

	return outputHandle;
}

Cell_handle_3 Triangular_Mesh_3::Create_Cell_3(std::vector<Vertex_handle_3>& vertices, int idx)
{
	// Create cell
	Cell_handle_3 outputHandle;
	outputHandle = mesh.tds().create_cell(vertices[0],vertices[1],vertices[2],vertices[3]);

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(4,0);

	// Set incident cells
	vertices[0]->set_cell(outputHandle);
	vertices[1]->set_cell(outputHandle);
	vertices[2]->set_cell(outputHandle);
	vertices[3]->set_cell(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;

	return outputHandle;
}

Cell_handle_3 Triangular_Mesh_3::Create_Cell_3(Vertex_handle_3 v0, Vertex_handle_3 v1, Vertex_handle_3 v2, Vertex_handle_3 v3, int idx)
{
	// Create cell
	Cell_handle_3 outputHandle;
	outputHandle = mesh.tds().create_cell(v0,v1,v2,v3);

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(4,0);

	// Set incident cells
	v0->set_cell(outputHandle);
	v1->set_cell(outputHandle);
	v2->set_cell(outputHandle);
	v3->set_cell(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;

	return outputHandle;
}

Cell_handle_3 Triangular_Mesh_3::Create_Infinite_Cell_3(int i0, int i1, int i2)
{
	// Create cell
	Cell_handle_3 outputHandle = mesh.tds().create_cell(
											mVertexHandleIndexMap[i2],
											mVertexHandleIndexMap[i1],
											mVertexHandleIndexMap[i0],
											mesh.infinite_vertex());

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(4,0);

	outputHandle->set_vertices(	mVertexHandleIndexMap[i2],
									mVertexHandleIndexMap[i1],
									mVertexHandleIndexMap[i0],
									mesh.infinite_vertex());

	// Set incident cells
	mVertexHandleIndexMap[i0]->set_cell(outputHandle);
	mVertexHandleIndexMap[i1]->set_cell(outputHandle);
	mVertexHandleIndexMap[i2]->set_cell(outputHandle);
	mesh.infinite_vertex()->set_cell(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = mSize_cells;

	return outputHandle;
}

bool Triangular_Mesh_3::TestForNeighbor_3(Finite_cells_iterator_3& cellHandleA, Finite_cells_iterator_3& cellHandleB,int idxNeigh)
{
	bool output = false;

	int idxIII = mesh.vertex_triple_index(idxNeigh,0);
	int idxJJJ = mesh.vertex_triple_index(idxNeigh,1);
	int idxKKK = mesh.vertex_triple_index(idxNeigh,2);

	//		Indexes used for the reciprocity relation
	int idxBNeigh = -1;
	int idxBIII = -1;
	int idxBJJJ = -1;
	int idxBKKK = -1;

	//		Handles for A's vertices
	Vertex_handle_3 testIII = cellHandleA->vertex(idxIII);
	Vertex_handle_3 testJJJ = cellHandleA->vertex(idxJJJ);
	Vertex_handle_3 testKKK = cellHandleA->vertex(idxKKK);

	// 		Watch out, the order of the indexes is inversed when changing the
	// tetrahedron:  III <-> JJJ
	if(cellHandleB->has_vertex(testIII,idxBJJJ)
			&& cellHandleB->has_vertex(testJJJ,idxBIII)
			&& cellHandleB->has_vertex(testKKK,idxBKKK)
			)
	{
		idxBNeigh	= mesh.next_around_edge(idxBJJJ,idxBIII);

		cellHandleA->set_neighbor(idxNeigh,cellHandleB);
		cellHandleB->set_neighbor(idxBNeigh,cellHandleA);
		
		cellHandleA->info().faceHasNeighbour[idxNeigh] = 1;
		cellHandleB->info().faceHasNeighbour[idxBNeigh] = 1;

		++mNbOfNeighs[cellHandleA->info().ExtIndex];
		++mNbOfNeighs[cellHandleB->info().ExtIndex];

		// Set up booleans
		output = true;
	}
	return output;
}

void Triangular_Mesh_3::ConnectTetrahedrons_3()
{
	mNbOfNeighs.resize(mSize_cells,0);
	Finite_cells_iterator_3 itCompare;
	bool areNeighbors = false;

	// --- Connect the (finite) tetrahedrons
	for(Finite_cells_iterator_3 itCells = mesh.finite_cells_begin(); itCells != mesh.finite_cells_end(); ++itCells)
	{
		// 		The neighbor relations are reciprocal, so the inner loop starts
		// at the current cell iterator + 1;
		Finite_cells_iterator_3 itCompare = itCells;
		++itCompare;

		for( ; itCompare != mesh.finite_cells_end(); ++itCompare)
		{
			// Test for neighbors ...
			for(int iii = 0; iii < 4; ++iii)
			{
				areNeighbors = TestForNeighbor_3(itCells,itCompare,iii);

				if(areNeighbors)
				{
					break;
				}
			}

			if(mNbOfNeighs[itCells->info().ExtIndex]==4)
			{
				// Already found all the 4 neighbors, stop the inner loop
				break;
			}
		}
	}
}

void Triangular_Mesh_3::AddInfiniteCells_3()
{
	int idxIII = -1;
	int idxJJJ = -1;
	int idxKKK = -1;

	int vertexIII = -1;
	int vertexJJJ = -1;
	int vertexKKK = -1;

	Vertex_handle_3 						infiniteVertex = mesh.infinite_vertex();
	std::unordered_map<int,Cell_handle_3>	infiniteCells;
	infiniteCells.reserve(4*mSize_cells);

	int FakeIndex = 0;

	// Build the infinite cells
	Cell_handle_3 itCells;
	for(int iii = 0; iii < mSize_cells; ++iii)
	{
		itCells = mCellHandleIndexMap[iii];

		if(mNbOfNeighs[itCells->info().ExtIndex]!=4)
		{
			for(int idxNeigh = 0; idxNeigh < 4; ++idxNeigh)
			{
				if(itCells->info().faceHasNeighbour[idxNeigh]==0)
				{
					// The cell has an empty neighbor, create an infinite cell
					idxIII = mesh.vertex_triple_index(idxNeigh,0);
					idxJJJ = mesh.vertex_triple_index(idxNeigh,1);
					idxKKK = mesh.vertex_triple_index(idxNeigh,2);

					vertexIII = itCells->vertex(idxIII)->info().ExtIndex;
					vertexJJJ = itCells->vertex(idxJJJ)->info().ExtIndex;
					vertexKKK = itCells->vertex(idxKKK)->info().ExtIndex;

					infiniteCells[FakeIndex] = Create_Infinite_Cell_3(vertexIII,vertexJJJ,vertexKKK);

					// Set as neighbor for the finite cell
					itCells->set_neighbor(idxNeigh,infiniteCells[FakeIndex]);
					itCells->info().faceHasNeighbour[idxNeigh]=1;

					// And the reciprocal relation
					for(int jjj = 0; jjj < 4; ++jjj)
					{
						if(mesh.is_infinite(infiniteCells[FakeIndex]->vertex(jjj)))
						{
							infiniteCells[FakeIndex]->set_neighbor(jjj,itCells);
							infiniteCells[FakeIndex]->info().faceHasNeighbour[jjj]=1;
						}
					}

					++mNbOfNeighs[itCells->info().ExtIndex];
					++FakeIndex;
				}
			}
		}
	}

//	// Set neighboring relations between the infinite cells
	int infiniteVertexIdx = 0;
	std::vector<int> FakeNbOfNeighs(FakeIndex,1);

	int idxBIII = -1;
	int idxBJJJ = -1;
	int idxBKKK = -1;

	for(int iii = 0; iii < FakeIndex; ++iii)
	{
		infiniteVertexIdx = infiniteCells[iii]->index(infiniteVertex);
		idxIII = mesh.vertex_triple_index(infiniteVertexIdx,0);
		idxJJJ = mesh.vertex_triple_index(infiniteVertexIdx,1);
		idxKKK = mesh.vertex_triple_index(infiniteVertexIdx,2);

		vertexIII = infiniteCells[iii]->vertex(idxIII)->info().ExtIndex;
		vertexJJJ = infiniteCells[iii]->vertex(idxJJJ)->info().ExtIndex;
		vertexKKK = infiniteCells[iii]->vertex(idxKKK)->info().ExtIndex;

		for(int jjj = iii + 1; jjj < FakeIndex; ++jjj)
		{
			if(FakeNbOfNeighs[iii]==4)
			{
				break;
			}

			if(infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexIII],idxBIII) &&
					infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexJJJ],idxBJJJ)
					)
			{
				idxBKKK = mesh.next_around_edge(idxBJJJ,idxBIII);
				infiniteCells[iii]->set_neighbor(idxKKK,infiniteCells[jjj]);
				infiniteCells[jjj]->set_neighbor(idxBKKK,infiniteCells[iii]);

				infiniteCells[iii]->info().faceHasNeighbour[idxKKK]=1;
				infiniteCells[jjj]->info().faceHasNeighbour[idxBKKK]=1;

				++FakeNbOfNeighs[iii];
				++FakeNbOfNeighs[jjj];
			}

			if(infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexJJJ],idxBJJJ) &&
					infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexKKK],idxBKKK)
					)
			{
				idxBIII = mesh.next_around_edge(idxBKKK,idxBJJJ);
				infiniteCells[iii]->set_neighbor(idxIII,infiniteCells[jjj]);
				infiniteCells[jjj]->set_neighbor(idxBIII,infiniteCells[iii]);

				infiniteCells[iii]->info().faceHasNeighbour[idxIII]=1;
				infiniteCells[jjj]->info().faceHasNeighbour[idxBIII]=1;

				++FakeNbOfNeighs[iii];
				++FakeNbOfNeighs[jjj];
			}

			if(infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexKKK],idxBKKK) &&
					infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexIII],idxBIII)
					)
			{
				idxBJJJ = mesh.next_around_edge(idxBIII,idxBKKK);
				infiniteCells[iii]->set_neighbor(idxJJJ,infiniteCells[jjj]);
				infiniteCells[jjj]->set_neighbor(idxBJJJ,infiniteCells[iii]);

				infiniteCells[iii]->info().faceHasNeighbour[idxJJJ]=1;
				infiniteCells[jjj]->info().faceHasNeighbour[idxBJJJ]=1;

				++FakeNbOfNeighs[iii];
				++FakeNbOfNeighs[jjj];
			}
		}
	}
}

void Triangular_Mesh_3::PrintDebugInfo()
{
	int debugIII;
	int debugCounter;

	for(All_vertices_iterator_3 itVertex = mesh.all_vertices_begin(); itVertex != mesh.all_vertices_end(); ++itVertex)
	{
		std::cout << "Vertex no. " << itVertex->info().ExtIndex << ": ";
		std::cout << "(";
		if(mesh.is_infinite(itVertex))
		{
			std::cout << "infty";
		}
		else
		{
			std::cout << itVertex->point();
		}
		std::cout << ")" << std::endl;
		if(!itVertex->is_valid())
		{
			std::cout << "   ---> INVALID VERTEX!!!" << std::endl;
		}
	}

	std::cout << " ----------------- " << std::endl;

	for(All_cells_iterator_3 itCells = mesh.all_cells_begin(); itCells != mesh.all_cells_end(); ++itCells)
	{
		debugCounter = 0;
		debugIII = itCells->info().ExtIndex;
		std::cout << "Cell no. " << debugIII << ": " << std::endl;
		std::cout << "   vertices  : ";
		for(int jjj = 0; jjj < 4; ++jjj)
		{
			std::cout << itCells->vertex(jjj)->info().ExtIndex << " ";
		}
		std::cout << std::endl;
//		std::cout << "   neighbors : ";
//		for(int jjj = 0; jjj < 4; ++jjj)
//		{
//			std::cout << "(";
//			if(itCells->info().faceHasNeighbour[jjj]==0)
//			{
//				std::cout << jjj << " , undef";
//			}
//			else
//			{
//				std::cout << jjj << " , " << itCells->neighbor(jjj)->info().ExtIndex;
//				++debugCounter;
//			}
//			std::cout << ") ";
//		}
		if(!itCells->is_valid())
		{
			std::cout << "   ---> INVALID CELL!!!" << std::endl;
		}

		std::cout << std::endl;
	}

	std::cout << " ----------------- " << std::endl;

	if(!mesh.tds().is_valid())
	{
		std::cout << "   ---> INVALID TDS!!!" << std::endl;
	}
	if(!mesh.is_valid(true,1))
	{
		std::cout << "   ---> INVALID MESH!!!" << std::endl;
	}
}

// PUBLIC Methods

// --- Getters
std::string Triangular_Mesh_3::get_name()
{
	return mName;
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

// --- Setters
void Triangular_Mesh_3::set_nb_of_facets()
{
	set_nb_of_vertices();
	mSize_facets = mesh.number_of_facets();
};

void Triangular_Mesh_3::set_nb_of_vertices()
{
	mSize_vertices = mesh.number_of_vertices();
};

void Triangular_Mesh_3::set_nb_of_cells()
{
	mSize_cells = mesh.number_of_finite_cells();
};

void Triangular_Mesh_3::set_indexes()
{
	// Set up indexes for the cells. Each finite cell receives a different
	// integer index, while the infinite ones receive the same index, equal to
	// mesh.number_of_finite_cells()

	int dummy = 0;
	int infiniteDummyCell = mesh.number_of_finite_cells();
	int infiniteDummyVertex = mesh.number_of_vertices() + 1;
	mCellHandleIndexMap.reserve(mesh.number_of_finite_cells() + 3);
	mVertexHandleIndexMap.reserve(mesh.number_of_vertices() + 3);

	for(All_cells_iterator_3	itCell = mesh.all_cells_begin();
			itCell !=  mesh.all_cells_end(); ++itCell)
	{
		if(!mesh.is_infinite(itCell))
		{
			itCell->info().ExtIndex = dummy;
			mCellHandleIndexMap[dummy] = itCell;
			++dummy;
		}
		else
		{
			itCell->info().ExtIndex = infiniteDummyCell;
		}
	}

	dummy = 1;
	for(All_vertices_iterator_3	itVertex = mesh.all_vertices_begin();
			itVertex !=  mesh.all_vertices_end(); ++itVertex)
	{
		if(!mesh.is_infinite(itVertex))
		{
			itVertex ->info().ExtIndex = dummy;
			mVertexHandleIndexMap[dummy] = itVertex;
			++dummy;
		}
		else
		{
			itVertex->info().ExtIndex = infiniteDummyVertex;
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

	set_nb_of_vertices();
	set_nb_of_facets();
	set_nb_of_cells();
	set_indexes();
};

void Triangular_Mesh_3::Initialize()
{
	mesh.clear();
	mesh.tds().set_dimension(3);
	
	mVertexHandleIndexMap.clear();
	mCellHandleIndexMap.clear();
	
	mNbOfNeighs.clear();
	
	mStarterCell = mesh.infinite_cell();
	mesh.infinite_vertex()->info().ExtIndex = -1;
}

void Triangular_Mesh_3::Finalize()
{
	mesh.tds().delete_cell(mStarterCell);
	set_nb_of_vertices();
	set_nb_of_cells();
}

//  --- Import an triangulation from a Gmsh file.
//void Triangular_Mesh_3::ImportGmsh(std::string &ifName)
//{
//	// Open data stream and check if the file exists
//	std::ifstream dataF(ifName);
//	assert(dataF.good());
//
//	// Initialize the mesh
//	Initialize();
//
//	// Set up safety booleans
//	bool hasHeader = false;
//	bool hasNodes = false;
//	bool hasElements = false;
//
//	// Variables needed by gmsh
//	unsigned int		gmshNumberOfNodes = 0;
//	unsigned int 		gmshNumberOfElements = 0;
//	double				bufferX = 1;
//	double				bufferY = 1;
//	int 				bufferNodeIndex = 1;
//	Point_3				bufferPoint;
//
//	int					bufferElementIndex = 1;
//	int					bufferElementType = 1;
//	int					bufferElementTagNumber = 1;
//	std::vector<int>	bufferElementTags(5);
//	std::vector<int>	bufferElementNodes(4);
//
//	double 				dummyZ = 1;
//	int					dummyExtIndex = 0;
//
//	// Variables needed to build the triangulation
//	Vertex_handle_3		workingVertex;
//
//	// Buffer string
//	std::string bufferLine;
//	std::stringstream	dataBuffer;
//
//
//	// Read info until the file ends
//	while(std::getline(dataF,bufferLine))
//	{
//		if(bufferLine.compare("$MeshFormat")==0)
//		{
//			// Read mesh format!
//			hasHeader = true;
//
////			// DEBUG
////			std::cout << bufferLine << std::endl;
////			std::cout << "2.2 0 " << sizeof(double) << std::endl;
////			std::cout << "$EndMeshFormat" << std::endl;
//		}
//
//		if(bufferLine.compare("$Nodes")==0)
//		{
//			// Read nodes!
//			hasNodes = true;
//			dataF >> gmshNumberOfNodes;
//
//			// Reserve space for the nodes
//			mVertexHandleIndexMap.reserve(2*gmshNumberOfNodes);
//
////			// DEBUG
////			std::cout << bufferLine << std::endl;
////			std::cout << gmshNumberOfNodes << std::endl;
//
//			// Line structure (* = to be ignored):
//			// [index] [X] [Y] [Z*]
//			std::getline(dataF,bufferLine);
//			for(unsigned int iii = 0; iii < gmshNumberOfNodes; ++iii)
//			{
//				std::getline(dataF,bufferLine);
//				dataBuffer.str("");
//				dataBuffer.clear();
//				dataBuffer << bufferLine;
//
//				dataBuffer >> bufferNodeIndex >> bufferX >> bufferY >> dummyZ;
//
////				// DEBUG
////				std::cout << bufferNodeIndex << " " << bufferX << " " << bufferY << " " << dummyZ << std::endl;
//
//				// Add point to intersection
//				bufferPoint = Point_3(bufferX,bufferY);
//				UpdateBbox(bufferPoint);
//				mVertexHandleIndexMap[bufferNodeIndex] = Create_Vertex_3(bufferPoint,bufferNodeIndex);
//			}
//
////			// DEBUG
////			std::cout << "$EndNodes" << std::endl;
//		}
//
//		if(bufferLine.compare("$Elements")==0)
//		{
//			// Read the elements!
//			hasElements = true;
//			dataF >> gmshNumberOfElements;
//
////			// DEBUG
////			std::cout << bufferLine << std::endl;
////			std::cout << gmshNumberOfElements << std::endl;
//
//			// Element structures (* = to be ignored):
//			// [index*] [type] [number of tags] [tags*] [nodes]
//
//			// Types :
//			//   2, 9, 20-25	: triangle
//			//   4, 11, 29-31	: tetrahedron
//			// Other types and extra nodes will be ignored (for now)
//			std::getline(dataF,bufferLine);
//			for(unsigned int iii = 0; iii < gmshNumberOfElements; ++iii)
//			{
//				std::getline(dataF,bufferLine);
//				dataBuffer.str("");
//				dataBuffer.clear();
//				dataBuffer << bufferLine;
//
//				dataBuffer 	>> bufferElementIndex
//							>> bufferElementType
//							>> bufferElementTagNumber;
//
////				// DEBUG
////				std::cout 	<< bufferElementIndex << " "
////							<< bufferElementType << " "
////							<< bufferElementTagNumber << " ";
//
//				if(bufferElementTagNumber > bufferElementTags.size())
//				{
//					bufferElementTags.resize(2*bufferElementTagNumber);
//				}
//
//				for(int jjj = 0; jjj < bufferElementTagNumber; ++jjj)
//				{
//					dataBuffer >> bufferElementTags[jjj];
//
////					// DEBUG
////					std::cout << bufferElementTags[jjj] << " ";
//				}
//
//				if	(	bufferElementType==2 ||
//						bufferElementType==9 ||
//						bufferElementType==20 ||
//						bufferElementType==21 ||
//						bufferElementType==22 ||
//						bufferElementType==23 ||
//						bufferElementType==24 ||
//						bufferElementType==25
//					)
//				{
//					// It's a triangle!
//					for(int jjj = 0; jjj < 3; ++jjj)
//					{
//						dataBuffer >> bufferElementNodes[jjj];
//
////						// DEBUG
////						std::cout << bufferElementNodes[jjj] << " ";
//					}
//
//					// Create the triangle
//					Create_Face_3(bufferElementNodes,dummyExtIndex);
//					++dummyExtIndex;
//				}
//
////				// DEBUG
////				std::cout << std::endl;
//			}
////			// DEBUG
////			std::cout << "$EndElements" << std::endl;
//		}
//	}
//
//	dataF.close();
//
//
//
//
//	// 		MUST remove the first triangle! CGAL starts the triangulations with
//	//	an almost empty triangle, connected to the infinite vertex, which is not
//	//	used by this algorithm.
//	Finalize();
//
//	/*
//	 * 		At this point, we built
//	 * 		- The vertices (all finite + 1 infinite)
//	 * 		- The faces (all finite)
//	 *
//	 * 		and we've set up
//	 * 		- Incident relations between vertices and finite faces
//	 * 		- Indexes for the finite faces
//	 *
//	 * 		We must
//	 * 		- Set up neighboring relations between the finite faces
//	 * 		- Build infinite faces
//	 * 		- Set up neighboring relations between the finite and infinite faces
//	 * 		- Set up incident relations with infinite faces
//	 * 		- Set up neighboring relations between the infinite faces
//	 *
//	 */
//
//	// Create neighbor relations between the FINITE faces
//	ConnectTriangles_3();
//
////	// DEBUG
////	PrintDebugInfo();
//}

//  --- Import an triangulation from a Medit file.
void Triangular_Mesh_3::ImportMedit(std::string &ifName)
{
	// Open data stream and check if the file exists
	std::ifstream dataF(ifName);
	assert(dataF.good());

	// Initialize the mesh
	Initialize();

	// Set up safety booleans
	bool hasHeader = false;
	bool hasNodes = false;
	bool hasElements = false;

	// Variables needed by medit
	unsigned int		meditNumberOfNodes = 0;
	unsigned int 		meditNumberOfElements = 0;
	unsigned int 		meditDimension = 0;
	double				bufferX = 1;
	double				bufferY = 1;
	double				bufferZ = 1;

	int 				bufferNodeIndex = 1;
	Point_3				bufferPoint;

	std::vector<int>	bufferElementNodes(4);

	int					dummyDomain = 1;
	int					dummyExtIndex = 0;

	// Variables needed to build the triangulation
	Vertex_handle_3		workingVertex;

	// Buffer string
	std::string bufferLine;
	std::stringstream	dataBuffer;


	// Read info until the file ends
	while(std::getline(dataF,bufferLine))
	{
//		if(bufferLine.compare("MeshVersionFormatted 2")==0)
//		{
//			// Read mesh format!
//			hasHeader = true;
//
//			// DEBUG
//			std::cout << bufferLine << std::endl;
//			std::cout << "2.2 0 " << sizeof(double) << std::endl;
//			std::cout << "$EndMeshFormat" << std::endl;
//		}

		if(bufferLine.find("Dimension")!=std::string::npos)
		{
			dataF >> meditDimension;
		}

		if(bufferLine.find("Vertices")!=std::string::npos)
		{
			assert(meditDimension==3);

			// Read nodes!
			hasNodes = true;
			dataF >> meditNumberOfNodes;

			// Reserve space for the nodes
			mVertexHandleIndexMap.reserve(2*meditNumberOfNodes);

//			// DEBUG
//			std::cout << bufferLine << std::endl;
//			std::cout << meditNumberOfNodes << std::endl;

			// Line structure (* = to be ignored):
			// [X] [Y] [Z] [domain*]
			std::getline(dataF,bufferLine);
			for(unsigned int iii = 0; iii < meditNumberOfNodes; ++iii)
			{
				std::getline(dataF,bufferLine);
				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				dataBuffer >> bufferX >> bufferY >> bufferZ >> dummyDomain;

//				// DEBUG
//				std::cout << bufferNodeIndex << " " << bufferX << " " << bufferY << " " << dummyZ << std::endl;

				// Add point to intersection
				bufferPoint = Point_3(bufferX,bufferY,bufferZ);
				UpdateBbox(bufferPoint);
				mVertexHandleIndexMap[bufferNodeIndex] = Create_Vertex_3(bufferPoint,bufferNodeIndex);
				++bufferNodeIndex;
			}

//			// DEBUG
//			std::cout << "$EndNodes" << std::endl;
		}

		if(bufferLine.find("Tetrahedra")!=std::string::npos)
		{
			// Read the triangles!
			hasElements = true;
			dataF >> meditNumberOfElements;

//			// DEBUG
//			std::cout << bufferLine << std::endl;
//			std::cout << meditNumberOfElements << std::endl;

			// Element structures (* = to be ignored):
			// [nodes] [domain*]

			mCellHandleIndexMap.reserve(2*meditNumberOfElements);

			std::getline(dataF,bufferLine);
			for(unsigned int iii = 0; iii < meditNumberOfElements; ++iii)
			{
				std::getline(dataF,bufferLine);
				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				for(int jjj = 0; jjj < 4; ++jjj)
				{
					dataBuffer >> bufferElementNodes[jjj];

//						// DEBUG
//						std::cout << bufferElementNodes[jjj] << " ";
				}

				// Create the triangle
				mCellHandleIndexMap[dummyExtIndex] = Create_Cell_3(bufferElementNodes,dummyExtIndex);
				++dummyExtIndex;


//				// DEBUG
//				std::cout << std::endl;
			}
//			// DEBUG
//			std::cout << "$EndElements" << std::endl;
		}
	}

	dataF.close();

	// 		MUST remove the first triangle! CGAL starts the triangulations with
	//	an almost empty triangle, connected to the infinite vertex, which is not
	//	used by this algorithm.
	Finalize();

	/*
	 * 		At this point, we built
	 * 		- The vertices (all finite + 1 infinite)
	 * 		- The faces (all finite)
	 *
	 * 		and we've set up
	 * 		- Incident relations between vertices and finite faces
	 * 		- Indexes for the finite faces
	 *
	 * 		We must
	 * 		- Set up neighboring relations between the finite faces
	 * 		- Build infinite faces
	 * 		- Set up neighboring relations between the finite and infinite faces
	 * 		- Set up incident relations with infinite faces
	 * 		- Set up neighboring relations between the infinite faces
	 *
	 */

	// Create neighbor relations between the FINITE faces
	ConnectTetrahedrons_3();

	// Create and connect the INFINITE cells
	AddInfiniteCells_3();
}

//  --- Export an triangulation in a Gmsh file.
//void Triangular_Mesh_3::ExportGmsh(std::string &ofName)
//{
//	std::ofstream dataF(ofName);
//	dataF.precision(15);
//
//	// Small precaution ...
//	set_nb_of_vertices();
//	set_nb_of_faces();
//
//	assert(mSize_faces > 0 && mSize_vertices > 3);
//
//	dataF << "$MeshFormat" << std::endl;
//	dataF << "2.2 0 " << sizeof(double) << std::endl;
//	dataF << "$EndMeshFormat" << std::endl;
//
//	dataF << "$Nodes" << std::endl;
//	dataF << mSize_vertices << std::endl;
//	for(Finite_vertices_iterator_3 	itVertex = mesh.finite_vertices_begin();
//									itVertex != mesh.finite_vertices_end();
//									++itVertex)
//	{
//		dataF << itVertex->info().ExtIndex << " " ;
//		dataF << itVertex->point() << " " << "0" << std::endl;
//	}
//	dataF << "$EndNodes" << std::endl;
//
//	dataF << "$Elements" << std::endl;
//	dataF << mSize_faces << std::endl;
//	for(Finite_faces_iterator_3 	itFace = mesh.finite_faces_begin();
//								itFace != mesh.finite_faces_end();
//								++itFace)
//	{
//		dataF << itFace->info().ExtIndex + 1 << " " ; // Gmsh index starts at 1
//		dataF << "2 2 1 1 "; // dummy tags
//		for(int iii = 0; iii < 3; ++ iii)
//		{
//			dataF << itFace->vertex(iii)->info().ExtIndex << " ";
//		}
//		dataF << std::endl;
//	}
//	dataF << "$EndElements" << std::endl;
//	dataF.close();
//
//}

//  --- Export an triangulation in a Medit file.
void Triangular_Mesh_3::ExportMedit(std::string &ofName)
{
	std::ofstream dataF(ofName);
	dataF.precision(15);

	// Small precaution ...
	set_nb_of_vertices();
	set_nb_of_cells();

	assert(mSize_cells > 0 && mSize_vertices > 4);

	dataF << " MeshVersionFormatted 2" << std::endl;
	dataF << " Dimension" << std::endl;
	dataF << " 3" << std::endl;

	dataF << " Vertices" << std::endl;
	dataF << mSize_vertices << std::endl;
	for(Finite_vertices_iterator_3 	itVertex = mesh.finite_vertices_begin();
									itVertex != mesh.finite_vertices_end();
									++itVertex)
	{
		dataF << "    " ;
		dataF << itVertex->point() << " " << "   1" << std::endl;
	}

//	int dummySize_facets = 0;
//	DT_3::Facet dummyFacet;
//	std::stringstream dummyBuffer;
//	for(Finite_facets_iterator_3	itFacet = mesh.finite_facets_begin();
//									itFacet != mesh.finite_facets_end();
//									++itFacet)
//	{
//		dummyFacet = *itFacet;
//		++dummySize_facets;
//		dummyBuffer << "    ";
//		for(int iii = 0; iii < 3; ++ iii)
//		{
//			dummyBuffer << itFacet->first->vertex((itFacet->second + iii)%4)->info().ExtIndex
//				  << " ";
//		}
//		dummyBuffer << "   1" << std::endl;
//	}
//
//	dataF << " Triangles" << std::endl;
//	dataF << dummySize_facets << std::endl;
//	dataF << dummyBuffer.str();

	dataF << " Tetrahedra" << std::endl;
	dataF << mSize_cells << std::endl;
	for(Finite_cells_iterator_3 	itCell = mesh.finite_cells_begin();
									itCell != mesh.finite_cells_end();
									++itCell)
	{
		dataF << "    ";
		for(int iii = 0; iii < 4; ++ iii)
		{
			dataF << itCell->vertex(iii)->info().ExtIndex << " ";
		}
		dataF << "   1" << std::endl;
	}
	dataF.close();

}
