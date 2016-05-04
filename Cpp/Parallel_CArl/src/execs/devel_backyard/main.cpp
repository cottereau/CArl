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
 * 		- DONE: decouple the intersection search from the intersection
 * 		        construction
 * 		- DONE: expand it to build the intersections
 * 		- DONE: assemble these algorithms for a single proc
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

struct parallel_intersection_test_params {
	std::string mesh_A;
	std::string mesh_B;
	std::string mesh_C;
	std::string output_base;

	carl::SearchMethod search_type;
};

void get_input_params(GetPot& field_parser,
		parallel_intersection_test_params& input_params) {

	// Set mesh files
	if (field_parser.search(3, "--meshA", "-mA", "MeshA")) {
		input_params.mesh_A = field_parser.next(
				input_params.mesh_A);
	} else {
		input_params.mesh_A = "meshes/3D/tests/test_intersection_A_1.msh";
	}

	if (field_parser.search(3, "--meshB", "-mB", "MeshB")) {
		input_params.mesh_B = field_parser.next(
				input_params.mesh_B);
	} else {
		input_params.mesh_B = "meshes/3D/tests/test_intersection_B_1.msh";
	}

	if (field_parser.search(3, "--meshC", "-mC", "MeshC")) {
		input_params.mesh_C = field_parser.next(
				input_params.mesh_C);
	} else {
		input_params.mesh_C = "meshes/3D/tests/test_intersection_C_1.msh";
	}

	if (field_parser.search(3, "--output", "-mO", "OutputBase")) {
		input_params.output_base = field_parser.next(
				input_params.output_base);
	} else {
		input_params.output_base = "meshes/3D/tests/output/test";
	}

	std::string search_type;
	if (field_parser.search(2, "--searchType", "SearchType")) {
		search_type = field_parser.next(
				search_type);
		if(search_type == "Front" || search_type == "front" || search_type == "FRONT")
		{
			input_params.search_type = carl::FRONT;
		}
		else if(search_type == "Brute" || search_type == "brute" || search_type == "BRUTE")
		{
			input_params.search_type = carl::BRUTE;
		}
		else if(search_type == "Both" || search_type == "both" || search_type == "BOTH")
		{
			input_params.search_type = carl::BOTH;
		}
		else
		{
			input_params.search_type = carl::BRUTE;
		}
	}
	else
	{
		input_params.search_type = carl::BRUTE;
	}
}
;

int main(int argc, char *argv[])
{
	// Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);
	libMesh::Parallel::Communicator& WorldComm = init.comm();

	// Set up inputs
	GetPot command_line(argc, argv);
	GetPot field_parser;
	std::string input_filename;

	if (command_line.search(1, "--inputfile")) {
		input_filename = command_line.next(input_filename);
		field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
	} else {
		field_parser = command_line;
	}

	parallel_intersection_test_params input_params;
	get_input_params(field_parser, input_params);

	// Read the three meshes
	libMesh::Mesh test_mesh_A(WorldComm);
	libMesh::Mesh test_mesh_B(WorldComm);
	libMesh::Mesh test_mesh_C(WorldComm);
	libMesh::Mesh test_mesh_I(WorldComm);

	test_mesh_A.read(input_params.mesh_A);
	test_mesh_B.read(input_params.mesh_B);
	test_mesh_C.read(input_params.mesh_C);

	// Set up the search
	carl::Intersection_Search search_coupling_intersections(test_mesh_A,test_mesh_B,test_mesh_C,test_mesh_I,input_params.output_base);

	// Search!
	search_coupling_intersections.BuildIntersections(input_params.search_type);

	return 0;
}


