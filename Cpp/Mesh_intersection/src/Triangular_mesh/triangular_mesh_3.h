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
//#include "defines_3D.h"

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

	/*
	 * --- Map linking each vertex index to a vertex handle. Needed because
	 * 		CGAL does not have a search by "info" functionality.
	 */
	std::unordered_map<int,Vertex_handle_3>		mVertexHandleIndexMap;

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

	/*
	 *  --- DEBUG: print all the vertices and the faces in an human-readable
	 *  	format.
	 */
	void 	PrintDebugInfo();

public:
	// Members
	DT_3		mesh;

	// Constructors
	Triangular_Mesh_3()
	{
		mSize_vertices = 0;
		mSize_facets = 0;
		mSize_cells = 0;
	}

	Triangular_Mesh_3(std::string &iName)
	{
		mName = iName;
		mSize_vertices = 0;
		mSize_facets = 0;
		mSize_cells = 0;
	}

	// Methods
	std::string get_name();

	void set_indexes();

	void GenerateTestMeshCube(const Point_3& initPoint, const Point_3& finalPoint, int nx, int ny, int nz, double amplitude = 0.15);

	void set_nb_of_cells();

	void set_nb_of_facets();

	void set_nb_of_vertices();

	int get_nb_of_cells() const;

	int get_nb_of_facets() const;

	int get_nb_of_vertices() const;

	/*
	 *  --- Initialize the mesh: clean up the data structures and set up the
	 *  	initial dimension.
	 */
	void Initialize();

	/*
	 *  --- Finalize the mesh: remove the dummy initial face
	 */
	void Finalize();
};
#endif /* TRIANGULAR_MESH_3_H_ */
