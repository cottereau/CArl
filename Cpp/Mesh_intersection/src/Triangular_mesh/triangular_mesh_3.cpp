/*
 * triangular_mesh_3.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
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
	outputHandle = mesh.tds().create_cell();

	outputHandle->set_vertices(	mVertexHandleIndexMap[vertices[0]],
								mVertexHandleIndexMap[vertices[1]],
								mVertexHandleIndexMap[vertices[2]],
								mVertexHandleIndexMap[vertices[3]]);

	outputHandle->set_neighbors();

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
	outputHandle = mesh.tds().create_cell();

	outputHandle->set_vertices(vertices[0],vertices[1],vertices[2],vertices[3]);

	outputHandle->set_neighbors();

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
	outputHandle = mesh.tds().create_cell();

	outputHandle->set_vertices(v0,v1,v2,v3);

	outputHandle->set_neighbors();

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
	Cell_handle_3 outputHandle = mesh.tds().create_cell();
	outputHandle->set_neighbors();

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

		++mNbOfNeighs[cellHandleA->info().ExtIndex];
		++mNbOfNeighs[cellHandleB->info().ExtIndex];
		// Set up boolean
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

	// Create and connect the INFINITE cells
	AddInfiniteCells_3();
}

void Triangular_Mesh_3::AddInfiniteCells_3()
{
	int idxIII = -1;
	int idxJJJ = -1;
	int idxKKK = -1;

	int vertexIII = -1;
	int vertexJJJ = -1;
	int vertexKKK = -1;

	Vertex_handle_3 			infiniteVertex = mesh.infinite_vertex();
	std::vector<Cell_handle_3>	infiniteCells(mSize_vertices,Cell_handle_3());
	int FakeIndex = 0;

	// Build the infinite cells
	for(Finite_cells_iterator_3 itCells = mesh.finite_cells_begin(); itCells != mesh.finite_cells_end(); ++itCells)
	{
		if(mNbOfNeighs[itCells->info().ExtIndex]!=4)
		{
			for(int idxNeigh = 0; idxNeigh < 4; ++idxNeigh)
			{
				if(itCells->neighbor(idxNeigh)==Cell_handle_3())
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

					// And the reciprocal relation
					for(int jjj = 0; jjj < 4; ++jjj)
					{
						if(mesh.is_infinite(infiniteCells[FakeIndex]->vertex(jjj)))
						{
							infiniteCells[FakeIndex]->set_neighbor(jjj,itCells);
						}
					}

					++mNbOfNeighs[itCells->info().ExtIndex];
					++FakeIndex;
				}
			}
		}
	}

	// Set neighboring relations between the infinite cells
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
				idxBKKK = mesh.next_around_edge(idxBIII,idxBJJJ);
				infiniteCells[iii]->set_neighbor(idxKKK,infiniteCells[jjj]);
				infiniteCells[jjj]->set_neighbor(idxBKKK,infiniteCells[iii]);

				++FakeNbOfNeighs[iii];
				++FakeNbOfNeighs[jjj];
			}

			if(infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexJJJ],idxBJJJ) &&
					infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexKKK],idxBKKK)
					)
			{
				idxBIII = mesh.next_around_edge(idxBJJJ,idxBKKK);
				infiniteCells[iii]->set_neighbor(idxIII,infiniteCells[jjj]);
				infiniteCells[jjj]->set_neighbor(idxBIII,infiniteCells[iii]);

				++FakeNbOfNeighs[iii];
				++FakeNbOfNeighs[jjj];
			}

			if(infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexKKK],idxBKKK) &&
					infiniteCells[jjj]->has_vertex(mVertexHandleIndexMap[vertexIII],idxBIII)
					)
			{
				idxBJJJ = mesh.next_around_edge(idxBKKK,idxBIII);
				infiniteCells[iii]->set_neighbor(idxJJJ,infiniteCells[jjj]);
				infiniteCells[jjj]->set_neighbor(idxBJJJ,infiniteCells[iii]);

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
		std::cout << "Face no. " << debugIII << ": " << std::endl;
		std::cout << "   vertices  : ";
		for(int jjj = 0; jjj < 4; ++jjj)
		{
			std::cout << itCells->vertex(jjj)->info().ExtIndex << " ";
		}
		std::cout << std::endl;
		std::cout << "   neighbors : ";
		for(int jjj = 0; jjj < 4; ++jjj)
		{
			std::cout << "(";
			if(itCells->neighbor(jjj)==Cell_handle_3())
			{
				std::cout << jjj << " , undef";
			}
			else
			{
				std::cout << jjj << " , " << itCells->neighbor(jjj)->info().ExtIndex;
				++debugCounter;
			}
			std::cout << ") ";
		}
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

void Triangular_Mesh_3::set_indexes()
{
	// Set up indexes for the cells. Each finite cell receives a different
	// integer index, while the infinite ones receive the same index, equal to
	// mesh.number_of_cells()

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

void Triangular_Mesh_3::Initialize()
{
	mesh.clear();
	mesh.tds().set_dimension(3);
	mVertexHandleIndexMap.clear();
	mNbOfNeighs.clear();
	mStarterCell = mesh.infinite_cell();
	mesh.infinite_vertex()->info().ExtIndex = -1;
}

void Triangular_Mesh_3::Finalize()
{
	mesh.tds().delete_cell(mStarterCell);
	set_nb_of_cells();
	set_nb_of_vertices();
}
