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
private:
	// Members
	std::string			mName;
	int					mSize_vertices;
	int					mSize_faces;

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
	 *  --- Connect the FINITE triangles and connect them
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
	}

	Triangular_Mesh_2(std::string &iName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_faces = 0;
	}

	Triangular_Mesh_2(std::string &iName, std::string &fileName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_faces = 0;
		if(boost::filesystem::path(fileName).extension().string().compare(".msh")==0)
		{
			// Its a Gmsh .msh file
			ImportGmsh(fileName);
		}
	}

	// Methods
	//  --- Getters
	std::string get_name();
	int get_nb_of_faces() const;
	int get_nb_of_vertices() const;

	//  --- Setters
	void set_nb_of_faces();
	void set_nb_of_vertices();
	void set_indexes();

	/*
	 * 	--- DEBUG Generate a dummy square triangulation, used for testing
	 */
	void GenerateTestMeshSquare(const Point_2& initPoint, const Point_2& finalPoint, int nx, int ny, double amplitude = 0.15);

	/*
	 *  --- Import an triangulation from a Gmsh file. ATTENTION: the current
	 *  	implementation was made only to be used with Gander's algorithm,
	 *  	where we only need the geometry information (nodes and triangular
	 *  	elements). As such, information like the tags and second and third
	 *  	order nodes are ignored, and a triangulation imported using this
	 *  	method will not be identical to the original when exported.
	 */
	void ImportGmsh(std::string &ifName);

	/*  --- Export an triangulation in a Gmsh file. ATTENTION: the resulting
	 * 		file may not be identical to the original one, read the description
	 * 		of "ImportGmsh"
	 */
	void ExportGmsh(std::string &ofName);

	/*
	 * 		TODO : implement Medit file import/export operations
	 */
//	void importMedit(std::string &ifName);
//
//	void exportMedit(std::string &ofName);
};
#endif /* TRIANGULAR_MESH_2_H_ */
