/*
 * main.cpp
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

 /** \file CArl_build_intersections.cpp
 \brief Implementation of the parallel intersection search
 
 Usage: `./CArl_build_intersections -i [input file]`

This program finds and constructs all the intersections between the meshes provided by the user. 
The input file is parsed by the carl::get_input_params(GetPot& field_parser, parallel_intersection_test_params& input_params) 
function, and it contains the following parameters. 

 Required parameters:
 - `MeshA`, `-mA` or `--meshA` : path to the mesh A.
 - `MeshB`, `-mB` or `--meshB` : path to the mesh B.
 - `MeshC`, `-mC` or `--meshC` : path to the coupling mesh C.

Optional parameters:
 - `OutputBase`, `-mO` or `--output` : base of the output files (including folders). *Default*: `test_inter`.
 - `MeshingMethod` or `--meshingMethodType` : intersection meshing method. *Values*: `CGAL` or `LIBMESH_TETGEN`. *Default*: `CGAL`.

Boolean flags:
 - `StitchInterMeshes` : do not stich together the intersection meshes. 
 - `VerboseOutput` or `--verbose` : print some extra information, such as the coupling mesh partitioning.
 */

#include "CArl_build_intersections.h"

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
	if (command_line.search(2, "--inputfile", "-i")) {
		input_filename = command_line.next(input_filename);
		field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
	} else {
		field_parser = command_line;
	}

	carl::parallel_intersection_params input_params;
	carl::get_intersection_input_params(field_parser, input_params);

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


	// If set to do "verbose" output, print the current partitioning of the intersection search
	if(input_params.bVerbose)
	{	
		std::cout << " -> Nb. coupling elements before repartitioning (intersection partition) : "
		          << test_mesh_C.n_partitions() << " ( ";
		for(unsigned int iii = 0; iii < nodes; ++iii)
		{
			std::cout << test_mesh_C.n_elem_on_proc(iii) << " ";
		}
		std::cout << ")" << std::endl << std::endl;
	}

	// Preallocate the intersection data structures (intersection tables, etc ...)
	perf_log.push("Prepare intersection load");
	search_coupling_intersections.PreparePreallocationAndLoad(/* SearchMethod = BRUTE */);
	search_coupling_intersections.PreallocateAndPartitionCoupling();
	perf_log.push("Prepare intersection load");

	// Print the current partitioning of the intersection search
	if(input_params.bVerbose)
	{	
		std::cout << " -> Nb. coupling elements after repartitioning (intersection partition) : " 
		          << test_mesh_C.n_partitions() << " ( ";
		for(unsigned int iii = 0; iii < nodes; ++iii)
		{
				std::cout << test_mesh_C.n_elem_on_proc(iii) << " ";
		}
		std::cout << ")" << std::endl << std::endl;
	}

	// --- Do the intersection search!
	perf_log.push("Search intersection");
	search_coupling_intersections.BuildIntersections(/* SearchMethod = BRUTE */);
	search_coupling_intersections.CalculateGlobalVolume();
	perf_log.pop("Search intersection");

	// --- Join the intersection tables, stitch the meshes, build the restrictions!
	// - Declarations
	// LOCAL stitched mesh
	libMesh::Mesh test_mesh_full_I(LocalComm,3);

	// Object used to stitch the meshes and to join the local intersection tables into a global one
	carl::Stitch_Meshes	join_meshes(test_mesh_full_I,input_params.output_base + "_global");
	join_meshes.set_grid_constraints(test_mesh_A,test_mesh_B);

	// Set filenames
	if(rank == 0)
	{
		join_meshes.set_base_filenames(input_params.output_base,".e",nodes);

		// Join the intersection tables
		perf_log.push("Join intersection tables");
		join_meshes.join_tables();
		perf_log.pop("Join intersection tables");

		// Stitch the intersection meshes
		if(input_params.bStitchInterMeshes)
		{
			perf_log.push("Stitch intersection meshes");
			join_meshes.stitch_meshes();
			perf_log.pop("Stitch intersection meshes");
		}
	}

	// Build and save the mesh restrictions

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

	return 0;
}


