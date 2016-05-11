/*
 * triangular_mesh_3.h
 *
 *  Created on: Jul 31, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef TRIANGULAR_MESH_3_H_
#define TRIANGULAR_MESH_3_H_

#include "common_header.h"
#include "CGAL_typedefs.h"

class	Triangular_Mesh_3
{
protected:
	// Members
	std::string			mName;
	int					mSize_vertices;
	int					mSize_facets;
	int					mSize_cells;
	Cell_handle_3		mStarterCell;

	Point_3				mPointMax;
	Point_3				mPointMin;

	// Inexact to exact converter
	Kernel_to_ExactKernel ConvertInexactToExact;

	/*
	 * --- Map linking each vertex index to a vertex handle. Needed because
	 * 		CGAL does not have a search by "info" functionality.
	 */
	std::unordered_map<int,Vertex_handle_3>		mVertexHandleIndexMap;

	/*
	 * --- Map linking each cell index to a cell handle. Needed because
	 * 		CGAL does not have a search by "info" functionality.
	 */
	std::unordered_map<int,Cell_handle_3>		mCellHandleIndexMap;

	/*
	 * --- Vector counting the number of neighbors of each triangle. Used to
	 *		speed up some of the mesh import operations
	 */
	std::vector<int>							mNbOfNeighs;

	// Methods

	/*
	 *  --- Create vertex at point "pointInput", with an index "indexInput"
	 */
	Vertex_handle_3 	Create_Vertex_3(Point_3& pointInput, int indexInput);

	/*
	 *  --- Create a FINITE cell using the vertices indexed as "i0", "i1", "i2",
	 *  	and associate to it the index "idx". The output is the corresponding
	 *  	cell handle.
	 */
	Cell_handle_3		Create_Cell_3(int i0, int i1, int i2, int i3, int idx);

	/*
	 *  --- Create a FINITE cell using the vertex handles "vertices", and
	 *  	associate to it the index "idx". The  output is the corresponding
	 *  	cell handle.
	 */
	Cell_handle_3		Create_Cell_3(std::vector<int>& vertices, int idx);

	Cell_handle_3 		Create_Cell_3(std::vector<Vertex_handle_3>& vertices, int idx);

	Cell_handle_3		Create_Cell_3(Vertex_handle_3 v0, Vertex_handle_3 v1, Vertex_handle_3 v2, Vertex_handle_3 v3, int idx);

	/*
	 *  --- Create a INFINITE cell using the vertices indexed as "i0", "i1". The
	 *  	output is the corresponding cell handle.
	 */
	Cell_handle_3 		Create_Infinite_Cell_3(int i0, int i1, int i2);

	/*
	 *  --- Test whenever "cellHandleB" is a neighbor of "cellHandleA" in the
	 *  	"idxNeigh" direction. If true, build the neighboring relations.
	 */
	bool	TestForNeighbor_3(Finite_cells_iterator_3& cellHandleA,
							  Finite_cells_iterator_3& cellHandleB,int idxNeigh);

	/*
	 *  --- Update the corner points, used to crate the bbox
	 */
	void UpdateBbox(Point_3& iPoint)
	{
		if(iPoint.x() < mPointMin.x())
		{
			mPointMin = Point_3(iPoint.x(),mPointMin.y(),mPointMin.z());
		}
		else if(iPoint.x() > mPointMax.x())
		{
			mPointMax = Point_3(iPoint.x(),mPointMax.y(),mPointMax.z());
		}

		if(iPoint.y() < mPointMin.y())
		{
			mPointMin = Point_3(mPointMin.x(),iPoint.y(),mPointMin.z());
		}
		else if(iPoint.y() > mPointMax.y())
		{
			mPointMax = Point_3(mPointMax.x(),iPoint.y(),mPointMax.z());
		}

		if(iPoint.z() < mPointMin.z())
		{
			mPointMin = Point_3(mPointMin.x(),mPointMin.y(),iPoint.z());
		}
		else if(iPoint.z() > mPointMax.z())
		{
			mPointMax = Point_3(mPointMax.x(),mPointMax.y(),iPoint.z());
		}
	}

public:
	// Members
	DT_3		mesh;

	// Constructors
	Triangular_Mesh_3()
	{
		mSize_vertices = 0;
		mSize_facets = 0;
		mSize_cells = 0;
		mPointMin = Point_3(0,0,0);
		mPointMax = Point_3(0,0,0);
	}

	Triangular_Mesh_3(std::string &iName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_facets = 0;
		mSize_cells = 0;
		mPointMin = Point_3(0,0,0);
		mPointMax = Point_3(0,0,0);
	}

	Triangular_Mesh_3(std::string &iName, std::string &fileName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_cells = 0;

		mPointMin = Point_3(0,0,0);
		mPointMax = Point_3(0,0,0);

		// TODO : better method to identify the mesh: read the first line
		//if(boost::filesystem::path(fileName).extension().string().compare(".msh")==0)
		//{
			// Its a Gmsh .msh file
			ImportGmsh(fileName);
		//}
		//else if(boost::filesystem::path(fileName).extension().string().compare(".mesh")==0)
		//{
//			Its a Medit file
		//	ImportMedit(fileName);
		//}
	}

	// Methods
	std::string get_name();

	void set_indexes();

	void GenerateTestMeshCube(const Point_3& initPoint, const Point_3& finalPoint, int nx, int ny, int nz, double amplitude = 0.15);

	/*
	 *  --- Connect the FINITE triangles and connect them. The current version
	 *  	is O(n^2) in time: the function has two for loops, with the outer
	 *  	one loop running over all the cells, and the inner one running over
	 *  	most of them (even if there are tests to stop it after finding 3
	 *  	neighbors)
	 *
	 *  	TODO	: optimize this algorithm
	 */
	void 	ConnectTetrahedrons_3();

	/*
	 *  --- Generate and connect the infinite cells. Must NOT be called before
	 *  	setting the connections between the finite triangles (it uses the
	 *  	information of which ones among these still have less than three
	 *  	neighbors to determinate where to build the infinite cells).
	 */
	void	AddInfiniteCells_3();

	void set_nb_of_cells();

	void set_nb_of_facets();

	void set_nb_of_vertices();

	int get_nb_of_cells() const;

	int get_nb_of_facets() const;

	int get_nb_of_vertices() const;

	double get_volume()
	{
		double vol = 0;
		Tetrahedron dummyTetrahedron;

		for(Finite_cells_iterator_3	itCell = mesh.finite_cells_begin();
									itCell != mesh.finite_cells_end();
									++itCell)
		{
			dummyTetrahedron = Tetrahedron(	itCell->vertex(0)->point(),
											itCell->vertex(1)->point(),
											itCell->vertex(2)->point(),
											itCell->vertex(3)->point());
			vol += dummyTetrahedron.volume();
		}
		return vol;
	}
	/*
	 *  --- DEBUG: print all the vertices and the faces in an human-readable
	 *  	format.
	 */
	void 	PrintDebugInfo(std::ostream& outStream = std::cout);

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
		return sqrt(XLength()*YLength()*ZLength()/mSize_cells);
	}

	double XLength()
	{
		return mPointMax.x()-mPointMin.x();
	}

	double YLength()
	{
		return mPointMax.y()-mPointMin.y();
	}

	double ZLength()
	{
		return mPointMax.z()-mPointMin.z();
	}

	/*
	 *  --- Initialize the mesh: clean up the data structures and set up the
	 *  	initial dimension.
	 */
	void Initialize(int reserveNbOfVertices = 0, int reserveNbOfCells = 0);

	/*
	 *  --- Finalize the mesh: remove the dummy initial face
	 */
	void Finalize();

	/*
	 *  --- Insert a vertex in the structure
	 */
	Vertex_handle_3 Add_Vertex(Point_3& inputPoint, int inputIdx);

	Cell_handle_3 Add_Cell(std::vector<int>& inputVertexList, int inputIdx);

	Cell_handle_3 Add_Cell(	std::vector<int>& inputVertexList, int inputIdx,
							int bufferElementType, int Ntags,
							std::vector<int>& tags);

	void ConvertToTriangulation_3(Polyhedron& t)
	{
		mesh.clear();
		mesh.insert(t.points_begin(),t.points_end());
		set_nb_of_cells();
		set_nb_of_facets();
		set_nb_of_vertices();

		for(Finite_cells_iterator_3	itCell = mesh.finite_cells_begin();
									itCell != mesh.finite_cells_end();
									++itCell)
		{
			itCell->info().ToAdd = false;
		}

		for(Finite_vertices_iterator_3	itVertex = mesh.finite_vertices_begin();
										itVertex != mesh.finite_vertices_end();
										++itVertex)
		{
			itVertex->info().ToAdd = false;
		}
	}

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

	void ExportMesh_3(std::string &ofName)
	{
		//if(boost::filesystem::path(ofName).extension().string().compare(".msh")==0)
		//{
			// Its a Gmsh .msh file
			ExportGmsh(ofName);
		//}
		//else if(boost::filesystem::path(ofName).extension().string().compare(".mesh")==0)
		//{
			// Its a Medit file
		//	ExportMedit(ofName);
		//}
	}

	/*
	 *  --- Restrict the mesh to only the elements intersecting the Nef
	 *      polyhedron "nefRestriction", and save the output to "outputMesh"
	 *
	 *      SHOULD NOT CHANGE THE CURRENT MESH !!!
	 */
	void RestrictMesh(Nef_Polyhedron& nefRestriction, Triangular_Mesh_3& outputMesh, const std::string tableFilename);
};
#endif /* TRIANGULAR_MESH_3_H_ */
