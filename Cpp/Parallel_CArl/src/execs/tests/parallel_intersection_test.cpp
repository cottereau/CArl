/*
 * main.cpp
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "parallel_intersection_test.h"

/*
 *
 *  	This program takes as input a configuration file which is parsed by the carl::get_input_params(GetPot& field_parser,
		parallel_intersection_test_params& input_params) function.
 */
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

	carl::parallel_intersection_test_params input_params;
	carl::get_input_params(field_parser, input_params);

	// Read the three meshes
	libMesh::Mesh test_mesh_A(WorldComm);
	libMesh::Mesh test_mesh_B(WorldComm);
	libMesh::Mesh test_mesh_C(WorldComm);

	// And set the (single processor) output mesh
	libMesh::Mesh test_mesh_I(LocalComm);

	test_mesh_A.read(input_params.mesh_A);
	test_mesh_B.read(input_params.mesh_B);
	test_mesh_C.read(input_params.mesh_C);

	test_mesh_A.prepare_for_use();
	test_mesh_B.prepare_for_use();
	test_mesh_C.prepare_for_use();

	test_mesh_C.partition(nodes);

	// Set up the search
	perf_log.push("Set up");
	carl::Intersection_Search search_coupling_intersections(test_mesh_A,test_mesh_B,test_mesh_C,test_mesh_I,input_params.output_base,input_params.inter_meshing_method);
	perf_log.pop("Set up");

	if(input_params.bSkipIntersectionConstruction)
	{
		search_coupling_intersections.SkipIntersectionConstruction(true);
		std::cout << "    ---> Skipping the intersection construction!" << std::endl << std::endl;
	}

	if(input_params.bSkipIntersectionPartitioning)
	{
		search_coupling_intersections.SkipIntersectionPartitioning(true);
		std::cout << "    ---> Skipping the intersection partitioning!" << std::endl << std::endl;
	}

	if(input_params.bExportScalingData)
	{
		search_coupling_intersections.SetScalingFiles(input_params.output_base);
	}

	std::cout << " -> Nb. coupling elements (intersection partition) : " << test_mesh_C.n_partitions() << " ( ";
	for(unsigned int iii = 0; iii < nodes; ++iii)
	{
			std::cout << test_mesh_C.n_elem_on_proc(iii) << " ";
	}
	std::cout << ")" << std::endl << std::endl;

	// Preallocate the data
	perf_log.push("Prepare intersection load");
	search_coupling_intersections.PreparePreallocationAndLoad(input_params.search_type);
	search_coupling_intersections.PreallocateAndPartitionCoupling();
	perf_log.push("Prepare intersection load");

	std::cout << " -> Nb. coupling elements (intersection partition) : " << test_mesh_C.n_partitions() << " ( ";
	for(unsigned int iii = 0; iii < nodes; ++iii)
	{
			std::cout << test_mesh_C.n_elem_on_proc(iii) << " ";
	}
	std::cout << ")" << std::endl << std::endl;

	// Search!
	perf_log.push("Search intersection");
	search_coupling_intersections.BuildIntersections(input_params.search_type);
	if(!input_params.bSkipIntersectionConstruction)
	{
		search_coupling_intersections.CalculateGlobalVolume();
	}
	perf_log.pop("Search intersection");

	// Join the intersection tables, stitch the meshes, build the restrictions!
	if(!input_params.bSkipIntersectionConstruction)
	{
		libMesh::Mesh test_mesh_full_I(LocalComm,3);
		carl::Stitch_Meshes	join_meshes(test_mesh_full_I,input_params.output_base + "_stitched");
		join_meshes.set_grid_constraints(test_mesh_A,test_mesh_B);

		if(rank == 0)
		{
			join_meshes.set_base_filenames(input_params.output_base,".e",nodes);
		}
		perf_log.push("Join intersection tables");
		if(rank == 0)
		{
			join_meshes.join_tables();
		}
		perf_log.pop("Join intersection tables");

		perf_log.push("Stitch intersection meshes");
		if(rank == 0 && !input_params.bSkipMeshStitching)
		{
			join_meshes.stitch_meshes();
		}
		perf_log.pop("Stitch intersection meshes");

		// Restrict the meshes!
		if(!input_params.bSkipRestriction)
		{
			const std::unordered_set<unsigned int> * restrict_set_A_ptr = join_meshes.get_restricted_set_pointer_first();
			const std::unordered_set<unsigned int> * restrict_set_B_ptr = join_meshes.get_restricted_set_pointer_second();

			perf_log.push("Restrict meshes");
			carl::Mesh_restriction restrict_A(test_mesh_A,LocalComm);
			restrict_A.BuildRestrictionFromSet(restrict_set_A_ptr);
			restrict_A.export_restriction_mesh(input_params.output_base + "_A_restriction");

			carl::Mesh_restriction restrict_B(test_mesh_B,LocalComm);
			restrict_B.BuildRestrictionFromSet(restrict_set_B_ptr);
			restrict_B.export_restriction_mesh(input_params.output_base + "_B_restriction");
			perf_log.pop("Restrict meshes");
		}
	}

	return 0;
}


