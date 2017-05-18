/*
 * main.cpp
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */
#include "main.h"

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
			input_params.search_type = carl::SearchMethod::FRONT;
		}
		else if(search_type == "Brute" || search_type == "brute" || search_type == "BRUTE")
		{
			input_params.search_type = carl::SearchMethod::BRUTE;
		}
		else if(search_type == "Both" || search_type == "both" || search_type == "BOTH")
		{
			input_params.search_type = carl::SearchMethod::BOTH;
		}
		else
		{
			input_params.search_type = carl::SearchMethod::BRUTE;
		}
	}
	else
	{
		input_params.search_type = carl::SearchMethod::BRUTE;
	}
}
;

int main(int argc, char *argv[])
{
	// Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);
	libMesh::Parallel::Communicator& WorldComm = init.comm();
	unsigned int nodes = WorldComm.size();
	unsigned int rank =  WorldComm.rank();
	libMesh::Parallel::Communicator LocalComm;
	WorldComm.split(rank,0,LocalComm);

	libMesh::PerfLog perf_log("Main program");

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

	// And set the (single processor) output mesh
	libMesh::Mesh test_mesh_I(LocalComm);

	test_mesh_A.read(input_params.mesh_A);
	test_mesh_B.read(input_params.mesh_B);
	test_mesh_C.read(input_params.mesh_C);

	// Set the processor ids for the coupling mesh elements
	unsigned int nb_coupling_elems = test_mesh_C.n_elem();
	std::vector<unsigned int> nb_coupling_elements_per_node(nodes,nb_coupling_elems/nodes);
	nb_coupling_elements_per_node[nodes - 1] = nb_coupling_elems - (nodes - 1)*(nb_coupling_elems/nodes);
	unsigned int coupling_counter = 0;
	for(unsigned int iii = 0; iii < nodes; ++iii)
	{
		for(unsigned int jjj = 0; jjj < nb_coupling_elements_per_node[iii]; ++jjj)
		{
			test_mesh_C.elem(coupling_counter)->processor_id(iii);
			++coupling_counter;
		}
	}

	// Debug output - libMesh automatically suppresses output from  processors
	// other the rank == 0.
	std::cout << " -> Nb. of processors                 : " << nodes << std::endl;
	std::cout << " -> Nb. coupling elements (partition) : " << nb_coupling_elems << " ( ";
	for(unsigned int iii = 0; iii < nodes; ++iii)
	{
		std::cout << nb_coupling_elements_per_node[iii] << " ";
	}
	std::cout << ")" << std::endl << std::endl;

	// Set up the search
	perf_log.push("Set up");
	carl::Intersection_Search search_coupling_intersections(test_mesh_A,test_mesh_B,test_mesh_C,test_mesh_I,input_params.output_base);
	perf_log.pop("Set up");

	// Search!
	perf_log.push("Search intersection");
	search_coupling_intersections.BuildIntersections(input_params.search_type);
	search_coupling_intersections.CalculateGlobalVolume();
	perf_log.pop("Search intersection");

	return 0;
}


