/*
 * main.cpp
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */
#include "patch_construction.h"

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

	// Generate a cube mesh
	libMesh::Mesh test_mesh(WorldComm);

	double L = 2.0;
	double t = 0.0;
	libMesh::MeshTools::Generation::build_cube	(	test_mesh,
													2, 2, 2,
													t - L/2, t + L/2,
													t - L/2, t + L/2,
													t - L/2, t + L/2,
													libMesh::TET4);

	// Generate the query mesh
	libMesh::Mesh query_mesh(WorldComm);

	double L_query= 1.5;
	double t_query= 0.0;
	libMesh::MeshTools::Generation::build_cube	(	query_mesh,
													2, 2, 2,
													t_query - L_query/2, t_query + L_query/2,
													t_query - L_query/2, t_query + L_query/2,
													t_query - L_query/2, t_query + L_query/2,
													libMesh::HEX8);

	// Save elements in files
	test_mesh.write("meshes/3D/tests/output/patch_test.msh");
	query_mesh.write("meshes/3D/tests/output/patch_query.msh");

	// First, let us test with only one element
	const libMesh::Elem * query_elem = query_mesh.elem(0);

	// Get intersection!
	std::unordered_set<int> Patch_Indexes;
	Patch_Indexes.reserve(test_mesh.n_elem());

	carl::Patch_construction test_patch(test_mesh);
	test_patch.BuildPatch(query_elem,Patch_Indexes);

	std::cout << " -> Intersection indexes: " << std::endl;
	for(std::unordered_set<int>::iterator 	it_begin = Patch_Indexes.begin();
											it_begin != Patch_Indexes.end();
											++it_begin)
	{
		std:: cout << " " << *it_begin + 1;
	}
	std::cout << std::endl;


	return 0;
}


