/*
 * triangular_mesh_2.h
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef TRIANGULAR_MESH_2_H_
#define TRIANGULAR_MESH_2_H_

#include "common_header.h"
#include "CGAL_typedefs.h"

class	Triangular_Mesh_2
{
protected:
	// Members
	std::string			mName;
	int					mSize_vertices;
	int					mSize_faces;
	Face_handle_2		mStarterFace;

	Point_2				mPointMax;
	Point_2				mPointMin;

	/*
	 * --- Map linking each vertex index to a vertex handle. Needed because
	 * 		CGAL does not have a search by "info" functionality.
	 */
	std::unordered_map<int,Vertex_handle_2>		mVertexHandleIndexMap;

	/*
	 * --- Vector counting the number of neighbors of each triangle. Used to
	 *		speed up some of the mesh import operations
	 */
	std::vector<int>							mNbOfNeighs;

	// Methods

	/*
	 *  --- Create vertex at point "pointInput", with an index "indexInput"
	 */
	Vertex_handle_2 	Create_Vertex_2(Point_2& pointInput, int indexInput);

	/*
	 *  --- Create a FINITE face using the vertices indexed as "i0", "i1", "i2",
	 *  	and associate to it the index "idx". The output is the corresponding
	 *  	face handle.
	 */
	Face_handle_2		Create_Face_2(int i0, int i1, int i2, int idx);

	/*
	 *  --- Create a FINITE face using the vertex handles "vertices", and
	 *  	associate to it the index "idx". The  output is the corresponding
	 *  	face handle.
	 */
	Face_handle_2		Create_Face_2(std::vector<int>& vertices, int idx);

	Face_handle_2 		Create_Face_2(std::vector<Vertex_handle_2>& vertices, int idx);

	Face_handle_2		Create_Face_2(Vertex_handle_2 v0, Vertex_handle_2 v1, Vertex_handle_2 v2, int idx);

	/*
	 *  --- Create a INFINITE face using the vertices indexed as "i0", "i1". The
	 *  	output is the corresponding face handle.
	 */
	Face_handle_2 		Create_Infinite_Face_2(int i0, int i1);

	/*
	 *  --- Test whenever "faceHandleB" is a neighbor of "faceHandleA" in the
	 *  	"idxNeigh" direction. If true, build the neighboring relations.
	 */
	bool	TestForNeighbor_2(Finite_face_iterator_2& faceHandleA,
							  Finite_face_iterator_2& faceHandleB,int idxNeigh);

	/*
	 *  --- Connect the FINITE triangles and connect them. The current version
	 *  	is O(n^2) in time: the function has two for loops, with the outer
	 *  	one loop running over all the faces, and the inner one running over
	 *  	most of them (even if there are tests to stop it after finding 3
	 *  	neighbors)
	 *
	 *  	TODO	: optimize this algorithm
	 */
	void 	ConnectTriangles_2();

	/*
	 *  --- Generate and connect the infinite faces. Must NOT be called before
	 *  	setting the connections between the finite triangles (it uses the
	 *  	information of which ones among these still have less than three
	 *  	neighbors to determinate where to build the infinite faces).
	 */
	void	AddInfiniteFaces_2();

	/*
	 *  --- Update the corner points, used to crate the bbox
	 */
	void UpdateBbox(Point_2& iPoint)
	{
		if(iPoint.x() < mPointMin.x())
		{
			mPointMin = Point_2(iPoint.x(),mPointMin.y());
		}
		else if(iPoint.x() > mPointMax.x())
		{
			mPointMax = Point_2(iPoint.x(),mPointMax.y());
		}

		if(iPoint.y() < mPointMin.y())
		{
			mPointMin = Point_2(mPointMin.x(),iPoint.y());
		}
		else if(iPoint.y() > mPointMax.y())
		{
			mPointMax = Point_2(mPointMax.x(),iPoint.y());
		}
	}

	/*
	 *  --- DEBUG: print all the vertices and the faces in an human-readable
	 *  	format.
	 */
	void 	PrintDebugInfo();

public:
	// Members
	CDT		mesh;

	// Constructors
	Triangular_Mesh_2()
	{
		mSize_vertices = 0;
		mSize_faces = 0;
		mPointMin = Point_2(0,0);
		mPointMax = Point_2(0,0);
	}

	Triangular_Mesh_2(std::string &iName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_faces = 0;
		mPointMin = Point_2(0,0);
		mPointMax = Point_2(0,0);
	}

	Triangular_Mesh_2(std::string &iName, std::string &fileName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_faces = 0;

		mPointMin = Point_2(0,0);
		mPointMax = Point_2(0,0);

		// TODO : better method to identify the mesh: read the first line
		if(boost::filesystem::path(fileName).extension().string().compare(".msh")==0)
		{
			// Its a Gmsh .msh file
			ImportGmsh(fileName);
		}
		else if(boost::filesystem::path(fileName).extension().string().compare(".mesh")==0)
		{
			// Its a Medit file
			ImportMedit(fileName);
		}
	}

	// Methods
	//  --- Getters
	std::string get_name();
	int get_nb_of_faces() const;
	int get_nb_of_vertices() const;

	/*
	 *  --- Returns the bounding box of the triangulation
	 */
	Bbox_2 bbox()
	{
		return Bbox_2(mPointMin.x(),mPointMin.y(),mPointMax.x(),mPointMax.y());
	}

	/*
	 *  --- Returns a length of the order of the mesh's average distance between
	 *      the vertices
	 */
	double LengthOrder()
	{
		return sqrt(XLength()*YLength()/mSize_faces);
	}

	double XLength()
	{
		return mPointMax.x()-mPointMin.x();
	}

	double YLength()
	{
		return mPointMax.y()-mPointMin.y();
	}

	//  --- Setters
	void set_nb_of_faces();
	void set_nb_of_vertices();
	void set_indexes();

	/*
	 * 	--- DEBUG Generate a dummy square triangulation, used for testing
	 */
	void GenerateTestMeshSquare(const Point_2& initPoint, const Point_2& finalPoint, int nx, int ny, double amplitude = 0.15);

	/*
	 *  --- Initialize the mesh: clean up the data structures and set up the
	 *  	initial dimension.
	 */
	void Initialize();

	/*
	 *  --- Finalize the mesh: remove the dummy initial face
	 */
	void Finalize();

	/*
	 *  --- Import an triangulation from a Gmsh file. ATTENTION: the current
	 *  	implementation was made only to be used with Gander's algorithm,
	 *  	where we only need the geometry information (nodes and triangular
	 *  	elements). As such, information like the tags and second and third
	 *  	order nodes are ignored, and a triangulation imported using this
	 *  	method will not be identical to the original when exported.
	 */
	void ImportGmsh(std::string &ifName);

	/*  --- Import an triangulation from a Medit file. ATTENTION: the current
	 *  	implementation was made only to be used with Gander's algorithm.
	 */
	void ImportMedit(std::string &ifName);

	/*  --- Export an triangulation in a Gmsh file. ATTENTION: the resulting
	 * 		file may not be identical to the original one, read the description
	 * 		of "ImportGmsh"
	 */
	void ExportGmsh(std::string &ofName);

	/*  --- Export an triangulation in a Medit file. ATTENTION: the resulting
	 * 		file may not be identical to the original one.
	 */
	void ExportMedit(std::string &ofName);
};
#endif /* TRIANGULAR_MESH_2_H_ */
