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

//void Triangular_Mesh_2::GenerateTestMeshSquare(const Point_2& initPoint, const Point_2& finalPoint, int nx, int ny)
//{
//	mesh.clear();
//
//	//double dx = (CGAL::to_double(finalPoint.x()) - CGAL::to_double(initPoint.x()))/(nx-1);
//	//double dy = (CGAL::to_double(finalPoint.y()) - CGAL::to_double(initPoint.y()))/(ny-1);
//	double dx = (finalPoint.x() - initPoint.x())/(nx-1);
//	double dy = (finalPoint.y() - initPoint.y())/(ny-1);
//
//	double dummyX;
//	double dummyY;
//
//	boost::random::lagged_fibonacci607 m_rng;
//	boost::random::uniform_real_distribution<double> shift(-1,1);
//
//	// Set boundaries
//
//	Vertex_handle_2 vA = mesh.insert(initPoint);
//	Vertex_handle_2 vB = mesh.insert(Point_2(initPoint.x(),finalPoint.y()));
//	Vertex_handle_2 vC = mesh.insert(finalPoint);
//	Vertex_handle_2 vD = mesh.insert(Point_2(finalPoint.x(),initPoint.y()));
//
//	mesh.insert_constraint(vA,vB);
//	mesh.insert_constraint(vB,vC);
//	mesh.insert_constraint(vC,vD);
//	mesh.insert_constraint(vD,vA);
//
////	dummyX = (CGAL::to_double(finalPoint.x()) + CGAL::to_double(initPoint.x()))/2 + 0.2*dx*shift(m_rng);
////	dummyY = (CGAL::to_double(finalPoint.y()) + CGAL::to_double(initPoint.y()))/2 + 0.2*dy*shift(m_rng);
//
//	dummyX = (finalPoint.x() + initPoint.x())/2 + 0.2*dx*shift(m_rng);
//	dummyY = (finalPoint.y() + initPoint.y())/2 + 0.2*dy*shift(m_rng);
//
//	mesh.insert(Point_2(dummyX,dummyY));
//	CGAL::refine_Delaunay_mesh_2(mesh,Criteria_2(0.125, std::min(dx,dy)));
//	set_nb_of_faces();
//	set_nb_of_vertices();
//};

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

	set_nb_of_faces();
	set_nb_of_vertices();
};

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

void Triangular_Mesh_2::Create_Vertex_2(Point_2& pointInput, int indexInput)
{
	Vertex_handle_2 outputHandle;
	outputHandle = mesh.tds().create_vertex();
	outputHandle->set_point(pointInput);
	outputHandle->info().ExtIndex = indexInput;
	mVertexHandleIndexMap[indexInput] = outputHandle;
}

void Triangular_Mesh_2::Create_Face_2(int i0, int i1, int i2, int idx)
{
	// Create face
	Face_handle_2 outputHandle;
	outputHandle = mesh.tds().create_face(	mVertexHandleIndexMap[i0],
											mVertexHandleIndexMap[i1],
											mVertexHandleIndexMap[i2]);

	outputHandle->set_neighbors();

	// Set incident faces
	mVertexHandleIndexMap[i0]->set_face(outputHandle);
	mVertexHandleIndexMap[i1]->set_face(outputHandle);
	mVertexHandleIndexMap[i2]->set_face(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;
}

void Triangular_Mesh_2::Create_Face_2(std::vector<int>& vertices, int idx)
{
	// Create face
	Face_handle_2 outputHandle;
	outputHandle = mesh.tds().create_face();

	outputHandle->set_vertices(	mVertexHandleIndexMap[vertices[0]],
								mVertexHandleIndexMap[vertices[1]],
								mVertexHandleIndexMap[vertices[2]]);

	outputHandle->set_neighbors();

	// Set incident faces
	mVertexHandleIndexMap[vertices[0]]->set_face(outputHandle);
	mVertexHandleIndexMap[vertices[1]]->set_face(outputHandle);
	mVertexHandleIndexMap[vertices[2]]->set_face(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = idx;
}

void Triangular_Mesh_2::Create_Infinite_Face_2(int i0, int i1, Face_handle_2& outputHandle)
{
	// Create face
	outputHandle = mesh.tds().create_face();
	outputHandle->set_neighbors();

	outputHandle->set_vertices(	mVertexHandleIndexMap[i0],
								mVertexHandleIndexMap[i1],
								mesh.infinite_vertex());

	// Set incident faces
	mVertexHandleIndexMap[i0]->set_face(outputHandle);
	mVertexHandleIndexMap[i1]->set_face(outputHandle);
	mesh.infinite_vertex()->set_face(outputHandle);

	// Set external index
	outputHandle->info().ExtIndex = mSize_faces;
}

void Triangular_Mesh_2::AddInfiniteFaces_2()
{
	int idxIII = -1;
	int idxJJJ = -1;

	int vertexIII = -1;
	int vertexJJJ = -1;

	Vertex_handle_2 			infiniteVertex = mesh.infinite_vertex();
	std::vector<Face_handle_2>	infiniteFaces(mSize_vertices,Face_handle_2());
	int FakeIndex = 0;

	// Build the infinite faces
	for(Finite_face_iterator_2 itFaces = mesh.finite_faces_begin(); itFaces != mesh.finite_faces_end(); ++itFaces)
	{
		if(mNbOfNeighs[itFaces->info().ExtIndex]!=3)
		{
			for(int idxNeigh = 0; idxNeigh < 3; ++idxNeigh)
			{
				if(itFaces->neighbor(idxNeigh)==Face_handle_2())
				{
					// The face has an empty neighbor, create an infinite face
					idxIII = mesh.ccw(idxNeigh);
					idxJJJ = mesh.ccw(idxIII);

					vertexIII = itFaces->vertex(idxIII)->info().ExtIndex;
					vertexJJJ = itFaces->vertex(idxJJJ)->info().ExtIndex;

					Create_Infinite_Face_2(vertexIII,vertexJJJ,infiniteFaces[FakeIndex]);

					// Set as neighbor for the finite face
					itFaces->set_neighbor(idxNeigh,infiniteFaces[FakeIndex]);

					// And the reciprocal relation
					for(int jjj = 0; jjj < 3; ++jjj)
					{
						if(mesh.is_infinite(infiniteFaces[FakeIndex]->vertex(jjj)))
						{
							infiniteFaces[FakeIndex]->set_neighbor(jjj,itFaces);
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

				++FakeNbOfNeighs[iii];
				++FakeNbOfNeighs[jjj];
			}

			if(infiniteFaces[jjj]->has_vertex(mVertexHandleIndexMap[vertexJJJ],idxBJJJ))
			{
				idxBIII = mesh.ccw(idxBJJJ);
				infiniteFaces[iii]->set_neighbor(idxIII,infiniteFaces[jjj]);
				infiniteFaces[jjj]->set_neighbor(idxBIII,infiniteFaces[iii]);

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
	for(All_faces_iterator_2 itFaces = mesh.all_faces_begin(); itFaces != mesh.all_faces_end(); ++itFaces)
	{
		debugCounter = 0;
		debugIII = itFaces->info().ExtIndex;
		std::cout << "Face no. " << debugIII << ": " << std::endl;
		std::cout << "   vertices  : ";
		for(int jjj = 0; jjj < 3; ++jjj)
		{
			std::cout << "(";
			if(mesh.is_infinite(itFaces->vertex(jjj)))
			{
				std::cout << "infty";
			}
			else
			{
				std::cout << itFaces->vertex(jjj)->point();
			}
			std::cout << ")";
		}
		std::cout << std::endl;
		std::cout << "   neighbors : ";
		for(int jjj = 0; jjj < 3; ++jjj)
		{
			std::cout << "(";
			if(itFaces->neighbor(jjj)==Face_handle_2())
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
	if(!mesh.tds().is_valid())
	{
		std::cout << "   ---> INVALID TDS!!!" << std::endl;
	}
	if(!mesh.is_valid(true,1))
	{
		std::cout << "   ---> INVALID MESH!!!" << std::endl;
	}

}

bool Triangular_Mesh_2::TestForNeighbor_2(Finite_face_iterator_2& faceHandleA, Finite_face_iterator_2& faceHandleB,int idxNeigh)
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

		// Set up boolean
		output = true;
	}
	return output;
}

void Triangular_Mesh_2::Connect_Triangles_2()
{
	mNbOfNeighs.resize(mSize_faces,0);
	Finite_face_iterator_2 itCompare;
	bool areNeighbors = false;

	// --- Connect the (finite) triangles
	for(Finite_face_iterator_2 itFaces = mesh.finite_faces_begin(); itFaces != mesh.finite_faces_end(); ++itFaces)
	{
		// 		The neighbor relations are reciprocal, so the inner loop starts
		// at the current face iterator + 1;
		Finite_face_iterator_2 itCompare = itFaces;
		++itCompare;

		for( ; itCompare != mesh.finite_faces_end(); ++itCompare)
		{
			// Test for first neighbor ...
			areNeighbors = TestForNeighbor_2(itFaces,itCompare,0);
			if(!areNeighbors)
			{
				// ... second neighbor ...
				areNeighbors = TestForNeighbor_2(itFaces,itCompare,1);
				if(!areNeighbors)
				{
					// ... and third neighbor
					areNeighbors = TestForNeighbor_2(itFaces,itCompare,2);
				}
			}

			if(areNeighbors)
			{
				++mNbOfNeighs[itFaces->info().ExtIndex];
				++mNbOfNeighs[itCompare->info().ExtIndex];
			}

			if(mNbOfNeighs[itFaces->info().ExtIndex]==3)
			{
				// Already found all the 3 neighbors, stop the inner loop
				break;
			}
		}
	}

	// --- Generate and connect the (finite) triangles
	AddInfiniteFaces_2();
}

void Triangular_Mesh_2::importGmsh(std::string &ifName)
{
	// Clear and set up the mesh
	mesh.clear();
	mesh.tds().set_dimension(2);

	/* 		Must remove the first face of the triangulation, created when the
	 * 	mesh was created.
	 */
	Face_handle_2 starterFace = mesh.infinite_face();
	mesh.infinite_vertex()->info().ExtIndex = -1;

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


	// Open data stream
	std::ifstream dataF(ifName);

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
				Create_Vertex_2(bufferPoint,bufferNodeIndex);
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

	mesh.tds().delete_face(starterFace);

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

	// Set up number of FINITE vertices and faces
	set_nb_of_faces();
	set_nb_of_vertices();

	// Create neighbor relations between the faces
	Connect_Triangles_2();

	// DEBUG
	PrintDebugInfo();
}
