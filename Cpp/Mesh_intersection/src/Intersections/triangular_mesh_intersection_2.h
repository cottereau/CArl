/*
 * triangular_mesh_intersection_2.h
 *
 *  Created on: Aug 25, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *      	Class used to save the intersection's mesh, and doing cleanup of the
 *      points. Set up an an inherited class to avoid bloating the mesh class
 *      with useless (in that context) spatial search methods. This header file
 *      also contains the specific typedefs needed.
 */

#ifndef INTERSECTIONS_TRIANGULAR_MESH_INTERSECTION_2_H_
#define INTERSECTIONS_TRIANGULAR_MESH_INTERSECTION_2_H_

#include "common_header.h"
#include "CGAL_typedefs.h"
#include "triangular_mesh_2.h"

// Local typedefs
typedef std::unordered_multimap<int, Face_handle_2>::iterator	faceRangeIterator;
typedef std::pair<faceRangeIterator,faceRangeIterator> 			faceRangeIteratorPair;

class	Intersection_Mesh_2 : public Triangular_Mesh_2
{
protected:

	// Variables used to build a mesh from an intersection
	int												mInterVertexDummyIndex;
	int												mInterFaceDummyIndex;
	std::vector<Vertex_handle_2>					mInterVertexHandle;
	std::vector<int>								mInterVertexIndex;
	std::unordered_map<long int,Vertex_handle_2>	mInterVertexHandleIndexMap;
	std::unordered_multimap<int, Face_handle_2>		mInterVertexIncidentFaces;

	double			mEps;

	long int mGridNx;
	long int mGridNy;

public:

	// --- Constructor
	Intersection_Mesh_2()
	{
		mInterVertexDummyIndex = -1;
		mInterFaceDummyIndex = -1;
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
	void AddPolygon(const Triangle_2& t);

	void AddPolygonVertex(int iii, Polygon_2& t);

	void AddPolygon(Polygon_2& t, int nbOfVertices, double CharacteristicArea);

	void CleanUp(double CharacteristicArea);

	void TestVertex(Vertex_handle_2 itVertex);

	long int ConvertToIndex(Point_2 iPoint);
};

#endif /* INTERSECTIONS_TRIANGULAR_MESH_INTERSECTION_2_H_ */
