/*
 * main.cpp
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */
#include "parallel_intersection_test.h"

struct parallel_intersection_test_params {
	std::string mesh_A;
	std::string mesh_B;
	std::string mesh_C;
	std::string output_base;

	carl::SearchMethod search_type;

	bool bSkipIntersectionConstruction;
	bool bSkipIntersectionPartitioning;
	bool bSkipRestriction;
	bool bSkipMeshStitching;
	bool bExportScalingData;

	carl::IntersectionMeshingMethod inter_meshing_method;
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

	if(field_parser.search(1,"SkipIntersectionConstruction")) {
		input_params.bSkipIntersectionConstruction = true;
	}
	else
	{
		input_params.bSkipIntersectionConstruction = false;
	}

	if(field_parser.search(1,"SkipIntersectionPartitioning")) {
		input_params.bSkipIntersectionPartitioning = true;
	}
	else
	{
		input_params.bSkipIntersectionPartitioning = false;
	}

	if(field_parser.search(1,"SkipRestriction")) {
		input_params.bSkipRestriction = true;
	}
	else
	{
		input_params.bSkipRestriction = false;
	}

	if(field_parser.search(1,"ExportScalingData")) {
		input_params.bExportScalingData = true;
	}
	else
	{
		input_params.bExportScalingData = false;
	}

	if(field_parser.search(1,"SkipMeshStitching")) {
		input_params.bSkipMeshStitching = true;
	}
	else
	{
		input_params.bSkipMeshStitching = false;
	}


	std::string meshing_method;
	if (field_parser.search(2, "--meshingMethodType", "MeshingMethod")) {
		meshing_method = field_parser.next(
				search_type);
		if(meshing_method == "CGAL")
		{
			input_params.inter_meshing_method = carl::IntersectionMeshingMethod::CGAL;
		}
		else if(meshing_method == "TETGEN" )
		{
			input_params.inter_meshing_method = carl::IntersectionMeshingMethod::LIBMESH_TETGEN;
		}
		else
		{
			input_params.inter_meshing_method = carl::IntersectionMeshingMethod::LIBMESH_TETGEN;
		}
	}
	else
	{
		input_params.inter_meshing_method = carl::IntersectionMeshingMethod::LIBMESH_TETGEN;
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

	// Stitch the meshes!
	if(!input_params.bSkipIntersectionConstruction)
	{
		perf_log.push("Stitch intersection meshes");
		libMesh::Mesh test_mesh_full_I(LocalComm,3);
		carl::Stitch_Intersection_Meshes	join_meshes(test_mesh_full_I,input_params.output_base + "_stitched");
		join_meshes.set_grid_constraints(test_mesh_A,test_mesh_B);

		if(rank == 0)
		{
			join_meshes.set_base_filenames(input_params.output_base,".e",nodes);
			join_meshes.join_tables();
			if(!input_params.bSkipMeshStitching)
			{
				join_meshes.stitch_meshes();
			}
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


