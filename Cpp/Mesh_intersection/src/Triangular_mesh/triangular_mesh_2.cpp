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

void Triangular_Mesh_2::Create_Face_2(int i0, int i1, int i2)
{
	// Create face
	Face_handle_2 outputHandle;
	outputHandle = mesh.tds().create_face(mVertexHandleIndexMap[i0],mVertexHandleIndexMap[i1],mVertexHandleIndexMap[i2]);
	outputHandle->set_neighbors();

	// Set incident faces
	mVertexHandleIndexMap[i0]->set_face(outputHandle);
	mVertexHandleIndexMap[i1]->set_face(outputHandle);
	mVertexHandleIndexMap[i2]->set_face(outputHandle);
}

void Triangular_Mesh_2::Create_Face_2(std::vector<int>& idx)
{
	// Create face
	Face_handle_2 outputHandle;
	outputHandle = mesh.tds().create_face();
	outputHandle->set_neighbors();

	Vertex_handle_2 dummyHandle0  = mVertexHandleIndexMap[idx[0]];
	Vertex_handle_2 dummyHandle1  = mVertexHandleIndexMap[idx[1]];
	Vertex_handle_2 dummyHandle2  = mVertexHandleIndexMap[idx[2]];
	outputHandle->set_vertices(dummyHandle0,dummyHandle1,dummyHandle2);

	// Set incident faces
	mVertexHandleIndexMap[idx[0]]->set_face(outputHandle);
	mVertexHandleIndexMap[idx[1]]->set_face(outputHandle);
	mVertexHandleIndexMap[idx[2]]->set_face(outputHandle);

	std::cout << "Got here?" << std::endl;
}

void Triangular_Mesh_2::Connect_Triangles_2()
{
	std::cout << mesh.tds().number_of_faces() << " " << mesh.number_of_vertices() << std::endl;

	for(Face_iterator_2 itFaces = mesh.faces_begin(); itFaces != mesh.faces_end(); ++itFaces)
	{
		std::cout << "I am a face!" << std::endl;
		if(mesh.is_infinite(itFaces))
		{
			std::cout << "This face is infinite !" << std::endl;
		}
		else
		{
			std::cout << "This face is finite !" << std::endl;
		}
	}

}

void Triangular_Mesh_2::importGmsh(std::string &ifName)
{
	// Clear the mesh
	mesh.clear();

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
			mVertexHandleIndexMap.reserve(gmshNumberOfNodes + 1);

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

						// DEBUG
						std::cout << bufferElementNodes[jjj] << " ";
					}

					// Create the triangle
					Create_Face_2(bufferElementNodes);
				}

//				// DEBUG
//				std::cout << std::endl;
			}
//			// DEBUG
//			std::cout << "$EndElements" << std::endl;
		}
	}

	dataF.close();

	// Now, must create neighbor relations between the faces
	Connect_Triangles_2();
}
