/*
 * triangular_mesh_intersection_3.h
 *
 *  Created on: Aug 25, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *      	Class used to save the intersection's mesh, and doing cleanup of the
 *      points. Set up an an inherited class to avoid bloating the mesh class
 *      with useless (in that context) spatial search methods. This header file
 *      also contains the specific typedefs needed.
 */

#ifndef INTERSECTIONS_TRIANGULAR_MESH_INTERSECTION_3_H_
#define INTERSECTIONS_TRIANGULAR_MESH_INTERSECTION_3_H_

#include "common_header.h"
#include "CGAL_typedefs.h"
#include "triangular_mesh_3.h"

// Local typedefs
typedef std::unordered_multimap<int, Cell_handle_3>::iterator	cellRangeIterator;
typedef std::pair<cellRangeIterator,cellRangeIterator> 			cellRangeIteratorPair;

class	Intersection_Mesh_3 : public Triangular_Mesh_3
{
protected:

	// Variables used to build a mesh from an intersection
	int												mInterVertexDummyIndex;
	int												mInterCellDummyIndex;
	std::vector<Vertex_handle_2>					mInterVertexHandle;
	std::vector<int>								mInterVertexIndex;
	std::unordered_map<long int,Vertex_handle_2>	mInterVertexHandleIndexMap;
	std::unordered_multimap<int, Cell_handle_3>		mInterVertexIncidentCells;

	double			mEps;

	long int mGridNx;
	long int mGridNy;

public:

	// --- Constructor
	Intersection_Mesh_3()
	{
		mInterVertexDummyIndex = -1;
		mInterCellDummyIndex = -1;
		mEps = -1;

		mGridNx = -1;
		mGridNy = -1;
	}
	// --- Methods
	/*
	 *  --- Initialize the mesh to save the intersection triangulation
	 */
	void InitializeIntersection(int& iPreallocation);

	/*
	 *  --- Triangle / polygon insertion, used with Gander's algorithm
	 */
	void AddPolyhedron(const Triangle_3& t);

	template<typename GeometryType>
	void AddPolyhedronVertex(int iii, GeometryType& t);

	void AddPolyhedron(Polyhedron& t, int nbOfVertices, double CharacteristicArea);

	void CleanUp();

	void TestVertex(Vertex_handle_3 itVertex);

	long int ConvertToIndex(Point_3 iPoint);
};

#endif /* INTERSECTIONS_TRIANGULAR_MESH_INTERSECTION_3_H_ */
