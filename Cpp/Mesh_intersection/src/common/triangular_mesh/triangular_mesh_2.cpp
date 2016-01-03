/*
 * triangular_mesh_2.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *		Methods of the class "Triangular_Mesh_2". A detailed description of each
 *		method can be found inside the header file.
 *
 */

#include "triangular_mesh_2.h"
// PRIVATE Methods

Vertex_handle_2 Triangular_Mesh_2::Create_Vertex_2(Point_2& pointInput, int indexInput)
{
	Vertex_handle_2 outputHandle;
	outputHandle = mesh.tds().create_vertex();
	outputHandle->set_point(pointInput);
	outputHandle->info().ExtIndex = indexInput;
	return outputHandle;
}

Face_handle_2 Triangular_Mesh_2::Create_Face_2(int i0, int i1, int i2, int idx)
{
	// Create face
	Face_handle_2 outputHandle;
	outputHandle = mesh.tds().create_face(	mVertexHandleIndexMap[i0],
											mVertexHandleIndexMap[i1],
											mVertexHandleIndexMap[i2]);

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(3,0);
	
	// Set incident faces
	mVertexHandleIndexMap[i0]->set_face(outputHandle);
	mVertexHandleIndexMap[i1]->set_face(outputHandle);
	mVertexHandleIndexMap[i2]->set_face(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;

	return outputHandle;
}

Face_handle_2 Triangular_Mesh_2::Create_Face_2(std::vector<int>& vertices, int idx)
{
	// Create face
	Face_handle_2 outputHandle;
	outputHandle = mesh.tds().create_face(
								mVertexHandleIndexMap[vertices[0]],
								mVertexHandleIndexMap[vertices[1]],
								mVertexHandleIndexMap[vertices[2]]);

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(3,0);
	
	// Set incident faces
	mVertexHandleIndexMap[vertices[0]]->set_face(outputHandle);
	mVertexHandleIndexMap[vertices[1]]->set_face(outputHandle);
	mVertexHandleIndexMap[vertices[2]]->set_face(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;

	return outputHandle;
}

Face_handle_2 Triangular_Mesh_2::Create_Face_2(std::vector<Vertex_handle_2>& vertices, int idx)
{
	// Create face
	Face_handle_2 outputHandle;
	outputHandle = mesh.tds().create_face(vertices[0],vertices[1],vertices[2]);

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(3,0);

	// Set incident faces
	vertices[0]->set_face(outputHandle);
	vertices[1]->set_face(outputHandle);
	vertices[2]->set_face(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;

	return outputHandle;
}

Face_handle_2 Triangular_Mesh_2::Create_Face_2(Vertex_handle_2 v0, Vertex_handle_2 v1, Vertex_handle_2 v2, int idx)
{
	// Create face
	Face_handle_2 outputHandle;
	outputHandle = mesh.tds().create_face(v0,v1,v2);

	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(3,0);

	// Set incident faces
	v0->set_face(outputHandle);
	v1->set_face(outputHandle);
	v2->set_face(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;

	return outputHandle;
}

Face_handle_2 Triangular_Mesh_2::Create_Infinite_Face_2(int i0, int i1)
{
	// Create face
	Face_handle_2 outputHandle = mesh.tds().create_face(
										mVertexHandleIndexMap[i1],
										mVertexHandleIndexMap[i0],
										mesh.infinite_vertex());
										
	outputHandle->set_neighbors();
	outputHandle->info().faceHasNeighbour.resize(3,0);
	
	// Set incident faces
	mVertexHandleIndexMap[i0]->set_face(outputHandle);
	mVertexHandleIndexMap[i1]->set_face(outputHandle);
	mesh.infinite_vertex()->set_face(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = mSize_faces;

	return outputHandle;
}

bool Triangular_Mesh_2::TestForNeighbor_2(Finite_faces_iterator_2& faceHandleA, Finite_faces_iterator_2& faceHandleB,int idxNeigh)
{
	bool output = false;

	int idxIII = mesh.ccw(idxNeigh);
	int idxJJJ = mesh.ccw(idxIII);

	//		Indexes used for the reciprocity relation
	int idxBNeigh = -1;
	int idxBIII = -1;
	int idxBJJJ = -1;

	//		Handles for A's vertices
	Vertex_handle_2 testIII = faceHandleA->vertex(idxIII);
	Vertex_handle_2 testJJJ = faceHandleA->vertex(idxJJJ);

	// 		Watch out, order of the indexes is inversed when changing the
	//	triangle:  III <-> JJJ
	if(faceHandleB->has_vertex(testIII,idxBJJJ) && faceHandleB->has_vertex(testJJJ,idxBIII))
	{
		//	Then we have neighbors, set up indexes
		idxNeigh	= mesh.ccw(idxJJJ);
		idxBNeigh	= mesh.ccw(idxBJJJ);

		faceHandleA->set_neighbor(idxNeigh,faceHandleB);
		faceHandleB->set_neighbor(idxBNeigh,faceHandleA);
		
		faceHandleA->info().faceHasNeighbour[idxNeigh] = 1;
		faceHandleB->info().faceHasNeighbour[idxBNeigh] = 1;

		++mNbOfNeighs[faceHandleA->info().ExtIndex];
		++mNbOfNeighs[faceHandleB->info().ExtIndex];

		// Set up booleans
		output = true;
	}
	return output;
}

void Triangular_Mesh_2::ConnectTriangles_2()
{
	mNbOfNeighs.resize(mSize_faces,0);
	Finite_faces_iterator_2 itCompare;
	bool areNeighbors = false;

	// --- Connect the (finite) triangles
	for(Finite_faces_iterator_2 itFaces = mesh.finite_faces_begin(); itFaces != mesh.finite_faces_end(); ++itFaces)
	{
		// 		The neighbor relations are reciprocal, so the inner loop starts
		// at the current face iterator + 1;
		Finite_faces_iterator_2 itCompare = itFaces;
		++itCompare;

		for( ; itCompare != mesh.finite_faces_end(); ++itCompare)
		{
			// Test for neighbors ...
			for(int iii = 0; iii < 3; ++iii)
			{
				areNeighbors = TestForNeighbor_2(itFaces,itCompare,iii);

				if(areNeighbors)
				{
					break;
				}
			}

			if(mNbOfNeighs[itFaces->info().ExtIndex]==3)
			{
				// Already found all the 3 neighbors, stop the inner loop
				break;
			}
		}
	}

	// Create and connect the INFINITE faces
	AddInfiniteFaces_2();
}

void Triangular_Mesh_2::AddInfiniteFaces_2()
{
	int idxIII = -1;
	int idxJJJ = -1;

	int vertexIII = -1;
	int vertexJJJ = -1;

	Vertex_handle_2 						infiniteVertex = mesh.infinite_vertex();
	std::unordered_map<int,Face_handle_2>	infiniteFaces;
	infiniteFaces.reserve(4*mSize_faces);
	
	int FakeIndex = 0;

	// Build the infinite faces
	Face_handle_2 itFaces;
	for(int iii = 0; iii < mSize_faces; ++iii)
	{
		itFaces = mFaceHandleIndexMap[iii];
		
		if(mNbOfNeighs[itFaces->info().ExtIndex]!=3)
		{
			for(int idxNeigh = 0; idxNeigh < 3; ++idxNeigh)
			{
				if(itFaces->info().faceHasNeighbour[idxNeigh]==0)
				{
					// The face has an empty neighbor, create an infinite face
					idxIII = mesh.ccw(idxNeigh);
					idxJJJ = mesh.ccw(idxIII);

					vertexIII = itFaces->vertex(idxIII)->info().ExtIndex;
					vertexJJJ = itFaces->vertex(idxJJJ)->info().ExtIndex;

					infiniteFaces[FakeIndex] = Create_Infinite_Face_2(vertexIII,vertexJJJ);

					// Set as neighbor for the finite face
					itFaces->set_neighbor(idxNeigh,infiniteFaces[FakeIndex]);
					itFaces->info().faceHasNeighbour[idxNeigh]=1;
					
					// And the reciprocal relation
					for(int jjj = 0; jjj < 3; ++jjj)
					{
						if(mesh.is_infinite(infiniteFaces[FakeIndex]->vertex(jjj)))
						{
							infiniteFaces[FakeIndex]->set_neighbor(jjj,itFaces);
							infiniteFaces[FakeIndex]->info().faceHasNeighbour[jjj]=1;
						}
					}

					++mNbOfNeighs[itFaces->info().ExtIndex];
					++FakeIndex;
				}
			}
		}
	}

	// Set neighboring relations between the infinite faces
	int infiniteVertexIdx = 0;
	std::vector<int> FakeNbOfNeighs(FakeIndex,1);

	int idxBIII = -1;
	int idxBJJJ = -1;

	for(int iii = 0; iii < FakeIndex; ++iii)
	{
		infiniteVertexIdx = infiniteFaces[iii]->index(infiniteVertex);
		idxIII = mesh.ccw(infiniteVertexIdx);
		idxJJJ = mesh.ccw(idxIII);

		vertexIII = infiniteFaces[iii]->vertex(idxIII)->info().ExtIndex;
		vertexJJJ = infiniteFaces[iii]->vertex(idxJJJ)->info().ExtIndex;

		for(int jjj = iii + 1; jjj < FakeIndex; ++jjj)
		{
			if(FakeNbOfNeighs[iii]==3)
			{
				break;
			}

			if(infiniteFaces[jjj]->has_vertex(mVertexHandleIndexMap[vertexIII],idxBIII))
			{
				idxBJJJ = mesh.cw(idxBIII);
				infiniteFaces[iii]->set_neighbor(idxJJJ,infiniteFaces[jjj]);
				infiniteFaces[jjj]->set_neighbor(idxBJJJ,infiniteFaces[iii]);

				infiniteFaces[iii]->info().faceHasNeighbour[idxJJJ]=1;
				infiniteFaces[jjj]->info().faceHasNeighbour[idxBJJJ]=1;
				
				++FakeNbOfNeighs[iii];
				++FakeNbOfNeighs[jjj];
			}

			if(infiniteFaces[jjj]->has_vertex(mVertexHandleIndexMap[vertexJJJ],idxBJJJ))
			{
				idxBIII = mesh.ccw(idxBJJJ);
				infiniteFaces[iii]->set_neighbor(idxIII,infiniteFaces[jjj]);
				infiniteFaces[jjj]->set_neighbor(idxBIII,infiniteFaces[iii]);
				
				infiniteFaces[iii]->info().faceHasNeighbour[idxIII]=1;
				infiniteFaces[jjj]->info().faceHasNeighbour[idxBIII]=1;		

				++FakeNbOfNeighs[iii];
				++FakeNbOfNeighs[jjj];
			}
		}
	}
}

void Triangular_Mesh_2::PrintDebugInfo()
{
	int debugIII;
	int debugCounter;

	for(All_vertices_iterator_2 itVertex = mesh.all_vertices_begin(); itVertex != mesh.all_vertices_end(); ++itVertex)
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

	for(All_faces_iterator_2 itFaces = mesh.all_faces_begin(); itFaces != mesh.all_faces_end(); ++itFaces)
	{
		debugCounter = 0;
		debugIII = itFaces->info().ExtIndex;
		std::cout << "Face no. " << debugIII << ": " << std::endl;
		std::cout << "   vertices  : ";
		for(int jjj = 0; jjj < 3; ++jjj)
		{
			std::cout << itFaces->vertex(jjj)->info().ExtIndex << " ";
		}
		std::cout << std::endl;
		std::cout << "   neighbors : ";
		for(int jjj = 0; jjj < 3; ++jjj)
		{
			std::cout << "(";
			if(itFaces->info().faceHasNeighbour[jjj]==0)
			{
				std::cout << jjj << " , undef";
			}
			else
			{
				std::cout << jjj << " , " << itFaces->neighbor(jjj)->info().ExtIndex;
				++debugCounter;
			}
			std::cout << ") ";
		}
		if(!itFaces->is_valid())
		{
			std::cout << "   ---> INVALID FACE!!!" << std::endl;
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

//  --- Getters
std::string Triangular_Mesh_2::get_name()
{
	return mName;
};

int Triangular_Mesh_2::get_nb_of_faces() const
{
	return mSize_faces;
};

int Triangular_Mesh_2::get_nb_of_vertices() const
{
	return mSize_vertices;
};

//  --- Setters
void Triangular_Mesh_2::set_nb_of_faces()
{
	mSize_faces = mesh.number_of_faces();
};

void Triangular_Mesh_2::set_nb_of_vertices()
{
	mSize_vertices = mesh.number_of_vertices();
};

void Triangular_Mesh_2::set_indexes()
{
	// Set up indexes for the faces. Each finite face receives a different
	// integer index, while the infinite ones receive the same index, equal to
	// mesh.number_of_faces()

	int dummy = 0;
	int infiniteDummy = mesh.number_of_faces();

	for(All_faces_iterator_2	itFace = mesh.all_faces_begin();
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

//  --- DEBUG Generate a dummy square triangulation, used for testing.
void Triangular_Mesh_2::GenerateTestMeshSquare(const Point_2& initPoint, const Point_2& finalPoint, int nx, int ny, double amplitude)
{
	mesh.clear();

	double dx = (finalPoint.x() - initPoint.x())/(nx-1);
	double dy = (finalPoint.y() - initPoint.y())/(ny-1);

	double dummyX;
	double dummyY;

	boost::random::uniform_real_distribution<double> shift(-1,1);

	// Set points

	// --- Corners
	mesh.insert(initPoint);
	mesh.insert(finalPoint);
	mesh.insert(Point_2(initPoint.x(),finalPoint.y()));
	mesh.insert(Point_2(finalPoint.x(),initPoint.y()));

	// --- Left and right border

	for(int iii = 1; iii < ny-1; ++iii)
	{
		dummyX = initPoint.x();
		dummyY = initPoint.y() + iii*dy + amplitude*dy*shift(m_rng);
		mesh.insert(Point_2(dummyX,dummyY));

		dummyX = finalPoint.x();
		dummyY = initPoint.y() + iii*dy + amplitude*dy*shift(m_rng);
		mesh.insert(Point_2(dummyX,dummyY));
	}

	// --- Top and bottom border

	for(int iii = 1; iii < nx-1; ++iii)
	{
		dummyY = initPoint.y();
		dummyX = initPoint.x() + iii*dx + amplitude*dx*shift(m_rng);
		mesh.insert(Point_2(dummyX,dummyY));

		dummyY = finalPoint.y();
		dummyX = initPoint.x() + iii*dx + amplitude*dx*shift(m_rng);
		mesh.insert(Point_2(dummyX,dummyY));
	}

	// Middle
	for(int iii = 1; iii < nx-1; ++iii)
	{
		for(int jjj = 1; jjj < ny-1; ++jjj)
		{
			dummyX = initPoint.x() + iii*dx + amplitude*dx*shift(m_rng);
			dummyY = initPoint.y() + jjj*dy + amplitude*dy*shift(m_rng);
			mesh.insert(Point_2(dummyX,dummyY));
		}
	}

	set_indexes();
	set_nb_of_faces();
	set_nb_of_vertices();
};

void Triangular_Mesh_2::Initialize()
{
	mesh.clear();
	mesh.tds().set_dimension(2);
	
	mVertexHandleIndexMap.clear();
	mFaceHandleIndexMap.clear();
	
	mNbOfNeighs.clear();
	
	mStarterFace = mesh.infinite_face();
	mesh.infinite_vertex()->info().ExtIndex = -1;
}

void Triangular_Mesh_2::Finalize()
{
	mesh.tds().delete_face(mStarterFace);
	set_nb_of_faces();
	set_nb_of_vertices();
}

//  --- Import an triangulation from a Gmsh file.
void Triangular_Mesh_2::ImportGmsh(std::string &ifName)
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

	// Variables needed by gmsh
	unsigned int		gmshNumberOfNodes = 0;
	unsigned int 		gmshNumberOfElements = 0;
	double				bufferX = 1;
	double				bufferY = 1;
	int 				bufferNodeIndex = 1;
	Point_2				bufferPoint;

	int					bufferElementIndex = 1;
	int					bufferElementType = 1;
	int					bufferElementTagNumber = 1;
	std::vector<int>	bufferElementTags(5);
	std::vector<int>	bufferElementNodes(4);

	double 				dummyZ = 1;
	int					dummyExtIndex = 0;

	// Variables needed to build the triangulation
	Vertex_handle_2		workingVertex;

	// Buffer string
	std::string bufferLine;
	std::stringstream	dataBuffer;


	// Read info until the file ends
	while(std::getline(dataF,bufferLine))
	{
		if(bufferLine.compare("$MeshFormat")==0)
		{
			// Read mesh format!
			hasHeader = true;

//			// DEBUG
//			std::cout << bufferLine << std::endl;
//			std::cout << "2.2 0 " << sizeof(double) << std::endl;
//			std::cout << "$EndMeshFormat" << std::endl;
		}

		if(bufferLine.compare("$Nodes")==0)
		{
			// Read nodes!
			hasNodes = true;
			dataF >> gmshNumberOfNodes;

			// Reserve space for the nodes
			mVertexHandleIndexMap.reserve(2*gmshNumberOfNodes);

//			// DEBUG
//			std::cout << bufferLine << std::endl;
//			std::cout << gmshNumberOfNodes << std::endl;

			// Line structure (* = to be ignored):
			// [index] [X] [Y] [Z*]
			std::getline(dataF,bufferLine);
			for(unsigned int iii = 0; iii < gmshNumberOfNodes; ++iii)
			{
				std::getline(dataF,bufferLine);
				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				dataBuffer >> bufferNodeIndex >> bufferX >> bufferY >> dummyZ;

//				// DEBUG
//				std::cout << bufferNodeIndex << " " << bufferX << " " << bufferY << " " << dummyZ << std::endl;

				// Add point to intersection
				bufferPoint = Point_2(bufferX,bufferY);
				UpdateBbox(bufferPoint);
				mVertexHandleIndexMap[bufferNodeIndex] = Create_Vertex_2(bufferPoint,bufferNodeIndex);
			}

//			// DEBUG
//			std::cout << "$EndNodes" << std::endl;
		}

		if(bufferLine.compare("$Elements")==0)
		{
			// Read the elements!
			hasElements = true;
			dataF >> gmshNumberOfElements;

//			// DEBUG
//			std::cout << bufferLine << std::endl;
//			std::cout << gmshNumberOfElements << std::endl;

			// Element structures (* = to be ignored):
			// [index*] [type] [number of tags] [tags*] [nodes]

			// Types :
			//   2, 9, 20-25	: triangle
			//   4, 11, 29-31	: tetrahedron
			// Other types and extra nodes will be ignored (for now)
			std::getline(dataF,bufferLine);
			for(unsigned int iii = 0; iii < gmshNumberOfElements; ++iii)
			{
				std::getline(dataF,bufferLine);
				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				dataBuffer 	>> bufferElementIndex
							>> bufferElementType
							>> bufferElementTagNumber;

//				// DEBUG
//				std::cout 	<< bufferElementIndex << " "
//							<< bufferElementType << " "
//							<< bufferElementTagNumber << " ";

				if(bufferElementTagNumber > bufferElementTags.size())
				{
					bufferElementTags.resize(2*bufferElementTagNumber);
				}

				for(int jjj = 0; jjj < bufferElementTagNumber; ++jjj)
				{
					dataBuffer >> bufferElementTags[jjj];

//					// DEBUG
//					std::cout << bufferElementTags[jjj] << " ";
				}

				if	(	bufferElementType==2 ||
						bufferElementType==9 ||
						bufferElementType==20 ||
						bufferElementType==21 ||
						bufferElementType==22 ||
						bufferElementType==23 ||
						bufferElementType==24 ||
						bufferElementType==25
					)
				{
					// It's a triangle!
					for(int jjj = 0; jjj < 3; ++jjj)
					{
						dataBuffer >> bufferElementNodes[jjj];

//						// DEBUG
//						std::cout << bufferElementNodes[jjj] << " ";
					}

					// Create the triangle
					Create_Face_2(bufferElementNodes,dummyExtIndex);
					++dummyExtIndex;
				}

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
	ConnectTriangles_2();

//	// DEBUG
//	PrintDebugInfo();
}

//  --- Import an triangulation from a Medit file.
void Triangular_Mesh_2::ImportMedit(std::string &ifName)
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
	int 				bufferNodeIndex = 1;
	Point_2				bufferPoint;

	std::vector<int>	bufferElementNodes(4);

	double				dummyZ = 1;
	int					dummyDomain = 1;
	int					dummyExtIndex = 0;

	// Variables needed to build the triangulation
	Vertex_handle_2		workingVertex;

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
			
			// Read nodes!
			hasNodes = true;
			dataF >> meditNumberOfNodes;

			// Reserve space for the nodes
			mVertexHandleIndexMap.reserve(2*meditNumberOfNodes);

//			// DEBUG
//			std::cout << bufferLine << std::endl;
//			std::cout << meditNumberOfNodes << std::endl;

			// Line structure (* = to be ignored):
			// [X] [Y] [domain*]
			std::getline(dataF,bufferLine);
			for(unsigned int iii = 0; iii < meditNumberOfNodes; ++iii)
			{
				std::getline(dataF,bufferLine);
				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				if(meditDimension==2)
				{
					dataBuffer >> bufferX >> bufferY >> dummyDomain;
				}
				else if(meditDimension==3)
				{
					dataBuffer >> bufferX >> bufferY >> dummyZ >> dummyDomain;
				}


//				// DEBUG
//				std::cout << bufferNodeIndex << " " << bufferX << " " << bufferY << " " << dummyZ << std::endl;

				// Add point to intersection
				bufferPoint = Point_2(bufferX,bufferY);
				UpdateBbox(bufferPoint);
				mVertexHandleIndexMap[bufferNodeIndex] = Create_Vertex_2(bufferPoint,bufferNodeIndex);
				++bufferNodeIndex;
			}

//			// DEBUG
//			std::cout << "$EndNodes" << std::endl;
		}

		if(bufferLine.find("Triangles")!=std::string::npos)
		{
			// Read the triangles!
			hasElements = true;
			dataF >> meditNumberOfElements;

//			// DEBUG
//			std::cout << bufferLine << std::endl;
//			std::cout << meditNumberOfElements << std::endl;

			// Element structures (* = to be ignored):
			// [nodes] [domain*]
			
			mFaceHandleIndexMap.reserve(2*meditNumberOfElements);

			std::getline(dataF,bufferLine);
			for(unsigned int iii = 0; iii < meditNumberOfElements; ++iii)
			{
				std::getline(dataF,bufferLine);
				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				for(int jjj = 0; jjj < 3; ++jjj)
				{
					dataBuffer >> bufferElementNodes[jjj];

//						// DEBUG
//						std::cout << bufferElementNodes[jjj] << " ";
				}

				// Create the triangle
				mFaceHandleIndexMap[dummyExtIndex] = Create_Face_2(bufferElementNodes,dummyExtIndex);
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
	ConnectTriangles_2();

//	// DEBUG
//	PrintDebugInfo();
}

//  --- Export an triangulation in a Gmsh file.
void Triangular_Mesh_2::ExportGmsh(std::string &ofName)
{
	std::ofstream dataF(ofName);
	dataF.precision(15);

	// Small precaution ...
	set_nb_of_vertices();
	set_nb_of_faces();

	assert(mSize_faces > 0 && mSize_vertices > 3);

	dataF << "$MeshFormat" << std::endl;
	dataF << "2.2 0 " << sizeof(double) << std::endl;
	dataF << "$EndMeshFormat" << std::endl;

	dataF << "$Nodes" << std::endl;
	dataF << mSize_vertices << std::endl;
	for(Finite_vertices_iterator_2 	itVertex = mesh.finite_vertices_begin();
									itVertex != mesh.finite_vertices_end();
									++itVertex)
	{
		dataF << itVertex->info().ExtIndex << " " ;
		dataF << itVertex->point() << " " << "0" << std::endl;
	}
	dataF << "$EndNodes" << std::endl;

	dataF << "$Elements" << std::endl;
	dataF << mSize_faces << std::endl;
	for(Finite_faces_iterator_2 	itFace = mesh.finite_faces_begin();
								itFace != mesh.finite_faces_end();
								++itFace)
	{
		dataF << itFace->info().ExtIndex + 1 << " " ; // Gmsh index starts at 1
		dataF << "2 2 1 1 "; // dummy tags
		for(int iii = 0; iii < 3; ++ iii)
		{
			dataF << itFace->vertex(iii)->info().ExtIndex << " ";
		}
		dataF << std::endl;
	}
	dataF << "$EndElements" << std::endl;
	dataF.close();

}

//  --- Export an triangulation in a Medit file.
void Triangular_Mesh_2::ExportMedit(std::string &ofName)
{
	std::ofstream dataF(ofName);
	dataF.precision(15);

	// Small precaution ...
	set_nb_of_vertices();
	set_nb_of_faces();

	assert(mSize_faces > 0 && mSize_vertices > 3);

	dataF << " MeshVersionFormatted 2" << std::endl;
	dataF << " Dimension" << std::endl;
	dataF << " 2" << std::endl;

	dataF << " Vertices" << std::endl;
	dataF << mSize_vertices << std::endl;
	for(Finite_vertices_iterator_2 	itVertex = mesh.finite_vertices_begin();
									itVertex != mesh.finite_vertices_end();
									++itVertex)
	{
		dataF << "    " ;
		dataF << itVertex->point() << " " << "   1" << std::endl;
	}

	dataF << " Triangles" << std::endl;
	dataF << mSize_faces << std::endl;
	for(Finite_faces_iterator_2 	itFace = mesh.finite_faces_begin();
								itFace != mesh.finite_faces_end();
								++itFace)
	{
		dataF << "    ";
		for(int iii = 0; iii < 3; ++ iii)
		{
			dataF << itFace->vertex(iii)->info().ExtIndex << " ";
		}
		dataF << "   1" << std::endl;
	}
	dataF.close();

}
