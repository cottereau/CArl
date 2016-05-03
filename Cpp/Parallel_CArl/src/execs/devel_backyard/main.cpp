/*
 * main.cpp
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */
#include "main.h"

/*
 * 		What do I have to do:
 *
 * 		- DONE: create and test a first intersection search algorithm
 * 		- DONE: create and test a patch construction algorithm
 * 		- DONE: create a FindIntersections algorithm
 * 		- TODO: decouple the intersection search from the intersection
 * 		        construction
 * 		- TODO: expand it to build the intersections
 * 		- TODO: assemble these algorithms for a single proc
 *
 * 		- TODO: (P) create a partitioner for \Omega_C
 * 		- TODO: (P) call the assembled algorithm for each part, printing each
 * 		            intersection mesh to a file ( |\Omega_C| files )
 * 		- TODO: (P) stitch the intersection meshes of a processor, printing now
 * 		            n_p files
 * 		- TODO: (P) either read several meshes in the assemble step, or stitch
 * 		            them all together ...
 *
 */
int main(int argc, char *argv[])
{
	// Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);
	libMesh::Parallel::Communicator& WorldComm = init.comm();

	// Read the three meshes
	libMesh::Mesh test_mesh_A(WorldComm);
	libMesh::Mesh test_mesh_B(WorldComm);
	libMesh::Mesh test_mesh_C(WorldComm);
	libMesh::Mesh test_mesh_I(WorldComm);

	test_mesh_A.read("meshes/3D/tests/test_intersection_A_1.msh");
	test_mesh_B.read("meshes/3D/tests/test_intersection_B_1.msh");
	test_mesh_C.read("meshes/3D/tests/test_intersection_C_1.msh");

	// Set up the search
	carl::Intersection_Search search_coupling_intersections(test_mesh_A,test_mesh_B,test_mesh_C,test_mesh_I);

	// Search!
	search_coupling_intersections.BuildIntersections(carl::BOTH);

	return 0;
}


