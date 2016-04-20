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
 * 		- TODO: create a FindIntersections algorithm
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
 * 		Classes that I'll need:
 *
 * 		- patch_construction: methods to build the intersection patch
 * 		                   -> no equivalent in current code
 *
 * 		- mesh_intersection_methods: methods to build the intersection mesh
 * 		  				          -> corresponds to "triangular_mesh_*"
 *
 * 		- intersection_search: methods to search and build the intersections
 * 		                    -> corresponds to an encapsulated version of
 * 		                       "intersection_functions_*"
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

	test_mesh_A.print_info();
	test_mesh_B.print_info();
	test_mesh_C.print_info();

	carl::Intersection_Search search_coupling_intersections(test_mesh_A,test_mesh_B,test_mesh_C,test_mesh_I);
	search_coupling_intersections.BuildIntersections();
//	const libMesh::Elem * query_elem_1 = test_mesh_C.elem(0);
//	search_coupling_intersections.BuildCoupledPatches(query_elem_1);
//
//	const libMesh::Elem * query_elem_2 = test_mesh_C.elem(1);
//	search_coupling_intersections.BuildCoupledPatches(query_elem_2);
//
//	// First, let us test with only one element
//	const libMesh::Elem * query_elem = query_mesh.elem(0);
//
//	// Get intersection!
//	std::unordered_set<int> Patch_Indexes;
//	Patch_Indexes.reserve(test_mesh.n_elem());
//
//	carl::Patch_construction test_patch(test_mesh);
//	test_patch.BuildPatch(query_elem,Patch_Indexes);
//
//	std::cout << " -> Intersection indexes: " << std::endl;
//	for(std::unordered_set<int>::iterator 	it_begin = Patch_Indexes.begin();
//											it_begin != Patch_Indexes.end();
//											++it_begin)
//	{
//		std:: cout << " " << *it_begin + 1;
//	}
//	std::cout << std::endl;
//

	return 0;
}


