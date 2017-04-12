/*
 * main.cpp
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "CArl_build_intersections.h"

/*	\brief Main test for the parallel intersection test.
 *
 *  	This program takes as input a configuration file which is parsed by the carl::get_input_params(GetPot& field_parser,
		parallel_intersection_test_params& input_params) function.
 */
int main(int argc, char *argv[])
{
	// --- Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);

	// libMesh's C++ / MPI communicator wrapper
	libMesh::Parallel::Communicator& WorldComm = init.comm();

	// Number of processors and processor rank.
	unsigned int nodes = WorldComm.size();
	unsigned int rank =  WorldComm.rank();

	// Create local communicator
	libMesh::Parallel::Communicator LocalComm;
	WorldComm.split(rank,rank,LocalComm);

	// Main program performance log
	libMesh::PerfLog perf_log("Main program");

	// --- Set up inputs

	// Command line parser
	GetPot command_line(argc, argv);

	// File parser
	GetPot field_parser;

	// If there is an input file, parse it to get the parameters. Else, parse the command line
	std::string input_filename;
	if (command_line.search(1, "--inputfile")) {
		input_filename = command_line.next(input_filename);
		field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
	} else {
		field_parser = command_line;
	}

	carl::parallel_intersection_test_params input_params;
	carl::get_input_params(field_parser, input_params);

	// Declare the three meshes to be intersected
	libMesh::Mesh test_mesh_A(WorldComm);
	libMesh::Mesh test_mesh_B(WorldComm);
	libMesh::Mesh test_mesh_C(WorldComm);

	// Declare the LOCAL output mesh
	libMesh::Mesh test_mesh_I(LocalComm);

	// Read the mesh files and prepare for use
	test_mesh_A.read(input_params.mesh_A);
	test_mesh_B.read(input_params.mesh_B);
	test_mesh_C.read(input_params.mesh_C);

	test_mesh_A.prepare_for_use();
	test_mesh_B.prepare_for_use();
	test_mesh_C.prepare_for_use();

	// --- Set up the Intersection_search object
	/*
	 *	This object is the main userinterface with the intersection search algorithms.
	 */
	perf_log.push("Set up");
	carl::Intersection_Search search_coupling_intersections(
			test_mesh_A,	// Input meshes
			test_mesh_B,	//
			test_mesh_C,	//
			
			test_mesh_I,	// LOCAL output mesh
			
			input_params.output_base,			// Common output filename path
			input_params.inter_meshing_method	// Intersection meshing method
			);
	perf_log.pop("Set up");

	// Skip the intersection construction?
	if(input_params.bSkipIntersectionConstruction)
	{
		search_coupling_intersections.SkipIntersectionConstruction(true);
		std::cout << "    ---> Skipping the intersection construction!" << std::endl << std::endl;
	}

	// Skip the intersection search repartitioning?
	if(input_params.bSkipIntersectionPartitioning)
	{
		search_coupling_intersections.SkipIntersectionPartitioning(true);
		std::cout << "    ---> Skipping the intersection partitioning!" << std::endl << std::endl;
	}

	// Export the scaling data?
	if(input_params.bExportScalingData)
	{
		search_coupling_intersections.SetScalingFiles(input_params.output_base);
	}

	// Print the current partitioning of the intersection search
	std::cout << " -> Nb. coupling elements (intersection partition) : " << test_mesh_C.n_partitions() << " ( ";
	for(unsigned int iii = 0; iii < nodes; ++iii)
	{
			std::cout << test_mesh_C.n_elem_on_proc(iii) << " ";
	}
	std::cout << ")" << std::endl << std::endl;

	// Preallocate the intersection data structures (intersection tables, etc ...)
	perf_log.push("Prepare intersection load");
	search_coupling_intersections.PreparePreallocationAndLoad(input_params.search_type);
	search_coupling_intersections.PreallocateAndPartitionCoupling();
	perf_log.push("Prepare intersection load");

	// Print the current partitioning of the intersection search
	std::cout << " -> Nb. coupling elements (intersection partition) : " << test_mesh_C.n_partitions() << " ( ";
	for(unsigned int iii = 0; iii < nodes; ++iii)
	{
			std::cout << test_mesh_C.n_elem_on_proc(iii) << " ";
	}
	std::cout << ")" << std::endl << std::endl;

	// --- Do the intersection search!
	perf_log.push("Search intersection");
	search_coupling_intersections.BuildIntersections(input_params.search_type);

	// Calculate the total intersection volume (to check for errors)
	if(!input_params.bSkipIntersectionConstruction)
	{
		search_coupling_intersections.CalculateGlobalVolume();
	}
	perf_log.pop("Search intersection");

	// --- Join the intersection tables, stitch the meshes, build the restrictions!
	if(!input_params.bSkipIntersectionConstruction)
	{
		// - Declarations
		// LOCAL stitched mesh
		libMesh::Mesh test_mesh_full_I(LocalComm,3);

		// Object used to stitch the meshes and to join the local intersection tables into a global one
		carl::Stitch_Meshes	join_meshes(test_mesh_full_I,input_params.output_base + "_stitched");
		join_meshes.set_grid_constraints(test_mesh_A,test_mesh_B);

		// Set filenames
		if(rank == 0)
		{
			join_meshes.set_base_filenames(input_params.output_base,".e",nodes);
		}

		// Join the intersection tables
		perf_log.push("Join intersection tables");
		if(rank == 0)
		{
			join_meshes.join_tables();
		}
		perf_log.pop("Join intersection tables");

		// Stitch the intersection meshes
		perf_log.push("Stitch intersection meshes");
		if(rank == 0 && !input_params.bSkipMeshStitching)
		{
			join_meshes.stitch_meshes();
		}
		perf_log.pop("Stitch intersection meshes");

		// Build and save the mesh restrictions
		if(!input_params.bSkipRestriction)
		{
			// Get set of elements forming the restricted mesh
			const std::unordered_set<unsigned int> * restrict_set_A_ptr = join_meshes.get_restricted_set_pointer_first();
			const std::unordered_set<unsigned int> * restrict_set_B_ptr = join_meshes.get_restricted_set_pointer_second();

			// Export the meshes and the element equivalence tables
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


