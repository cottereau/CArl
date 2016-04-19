/*
 * main.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#include "main.h"

#include <tuple>

/*
 *
 * 		WARNING! BIG CAVEAT OF THE GANDER ALGORITHM!
 *
 * 		Essentially, the algorithm cannot guarantee the correct construction of
 * 	the intersections if one of the meshes is non-convex or non-connected.
 *
 * 		The basis of the algorithm is to define a "frontier" over which it can
 * 	build the intersections. It supposes that, for a given element B_i, all of
 * 	its intersecting elements {A_j} are in a connected ensemble - form a
 * 	"frontier". If the mesh A is non-convex (and moreso if it is non-connected),
 * 	we can find elements B_i intersecting two or more non-connected regions of
 * 	the mesh A, and since B_i will be marked as "already visited" when building
 * 	the intersections with the first region, the other ones will be ignored.
 *
 */

struct carl_intersection_input_params {
	std::string mesh_A_file;
	std::string mesh_B_file;
	std::string restriction_file;

	std::string search_method;

	std::string mesh_I_file;
	std::string intersection_table_full_file;
};

void get_input_params(GetPot& field_parser,
	carl_intersection_input_params& input_params) {

	// Set mesh files
	bool bHasMeshA = field_parser.search(3, "--meshA", "-mA", "MeshA");
	bool bHasMeshB = field_parser.search(3, "--meshB", "-mB", "MeshB");

	bool bUseDefaultValues = !bHasMeshB && !bHasMeshA;

	if (bUseDefaultValues) {
		// Then use default values for everything!
		std::cout << " > Using the default inputs parameters" << std::endl;
		input_params.mesh_A_file = "meshes/3D/test_restriction_cube.msh";
		input_params.mesh_B_file = "meshes/3D/test_restriction_micro.msh";
		input_params.restriction_file = "meshes/3D/restriction_data.txt";
	} else if (bHasMeshA && bHasMeshB) {
		field_parser.search(3, "--meshA", "-mA", "MeshA");
		input_params.mesh_A_file = field_parser.next(
				input_params.mesh_A_file);

		field_parser.search(3, "--meshB", "-mB", "MeshB");
		input_params.mesh_B_file = field_parser.next(
				input_params.mesh_B_file);
	} else {
		if (!bHasMeshA) {
			homemade_error_msg("Missing the probe mesh file (A)!");
		}
		if (!bHasMeshB) {
			homemade_error_msg("Missing the guide mesh file (B)!");
		}
	}

	// Set the restriction file
	bool bHasRestriction = field_parser.search(2, "--restriction",  "IntersectionRestriction");

	if (bHasRestriction)
	{
		field_parser.search(2, "--restriction",  "IntersectionRestriction");
		input_params.restriction_file = field_parser.next(
				input_params.restriction_file);
	} else if (!bUseDefaultValues) {
		homemade_error_msg("Missing the restriction file!");
	}

	// Set the search method
	bool bHasSearchMethod = field_parser.search(2, "--searchMethod", "SearchMethod");

	if (bHasSearchMethod) {
		field_parser.search(2, "--searchMethod", "SearchMethod");
		input_params.search_method = field_parser.next(
				input_params.search_method);

		if(!(input_params.search_method == "FRONT") &&
		   !(input_params.search_method == "TREE" ) &&
		   !(input_params.search_method == "FULL" ) )
		{
			std::cout << " > Unknown search method, using the default (FRONT)" << std::endl;
		}
	} else {
		std::cout << " > Using the default (TREE)" << std::endl;
		input_params.search_method = "TREE";
	}

	// Set the output file
	bool bHasMeshI = field_parser.search(3, "--meshI", "-mI", "MeshInter");
	bool bHasFullInterTable = field_parser.search(2, "--tableFullI", "FullIntersectionElementsTable");

	if (bHasMeshI) {
		field_parser.search(3, "--meshI", "-mI", "MeshInter");
		input_params.mesh_I_file = field_parser.next(
				input_params.mesh_I_file);
	} else {
		std::cout << " > Using the default output mesh file" << std::endl;
		input_params.mesh_I_file = "meshes/3D/output/test_restriction_inter_I_3D.msh";
	}

	if (bHasFullInterTable){
		field_parser.search(2, "--tableFullI", "FullIntersectionElementsTable");
		input_params.intersection_table_full_file = field_parser.next(
				input_params.intersection_table_full_file);
	} else {
		std::cout << " > Using the default output intersection table file" << std::endl;
		input_params.intersection_table_full_file = "meshes/equivalence_tables/intersection_carl_Full_I.dat";
	}
}

int main(int argc, char *argv[])
{
	const bool MASTER_bPerfLog_carl_intersection = true;
	libMesh::LibMeshInit init(argc, argv);

	libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_intersection);

	// - Set up inputs
	GetPot command_line(argc, argv);
	GetPot field_parser;
	std::string input_filename;

	if (command_line.search(1, "--inputfile")) {
		input_filename = command_line.next(input_filename);
		field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
	} else {
		field_parser = command_line;
	}

	carl_intersection_input_params input_params;
	get_input_params(field_parser, input_params);

	// --- Preamble - timing variables
	std::chrono::time_point<std::chrono::system_clock> 	timing_start,
														timing_end,
														timing_start_total,
														timing_end_total;

	std::chrono::duration<double> 	elapsed_seconds_input,
									elapsed_seconds_intersections,
									elapsed_seconds_export,
									elapsed_seconds_total;

	std::cout << " > Input parameters : " << std::endl
			  << "   Probe mesh : " << input_params.mesh_A_file << std::endl
			  << "   Guide mesh : " << input_params.mesh_B_file << std::endl
			  << "   Restriction file  : " << input_params.restriction_file << std::endl
			  << "   Search method     : " << input_params.search_method << std::endl << std::endl
			  << "   Intersection mesh : " << input_params.mesh_I_file << std::endl
			  << "   Intersection table : " << input_params.intersection_table_full_file << std::endl;

	// Import meshes
	timing_start 		= std::chrono::system_clock::now();
	timing_start_total 	= std::chrono::system_clock::now();

	std::cout << " ---> Importing mesh ... ";
	std::string sillyname;
	sillyname = "Testing A";
	Triangular_Mesh_3 dtA(sillyname,input_params.mesh_A_file);

	sillyname = "Testing B";
	Triangular_Mesh_3 dtB(sillyname,input_params.mesh_B_file);
	std::cout << "finished!" << std::endl;

	std::cout << " ---> Generating the Nef polyhedron ... " << std::endl;
	Nef_Polyhedron couplingRegion;
	if(argc == 3)
	{
		std::cout << " ---> Using the full space" << std::endl;
		couplingRegion.clear(Nef_Polyhedron::COMPLETE);
	}
	else
	{
		GenerateNefRestrictedRegion(input_params.restriction_file,couplingRegion);
	}
	std::cout << " ---> Generating the Nef polyhedron ... finished!" << std::endl;

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_input = timing_end-timing_start;

	// ****************************** //
	// Intersection                   //
	// ****************************** //
	timing_start = std::chrono::system_clock::now();

	std::cout << " ---> Building intersections ... ";
	std::cout.flush();
	TriangulationIntersectionVisitor_3 	tallyFastVector(std::min(dtA.LengthOrder(),dtB.LengthOrder()),
														dtA.get_nb_of_cells(),dtB.get_nb_of_cells());

	BuildMeshIntersections(dtA,dtB,tallyFastVector,couplingRegion);
	std::cout << "finished!" << std::endl;

	timing_end   = std::chrono::system_clock::now();
	elapsed_seconds_intersections = timing_end-timing_start;

	timing_start = std::chrono::system_clock::now();

	// Export the meshes
	std::cout << " ---> Exporting data ... ";
	tallyFastVector.IntersectionTDS_3.ExportMesh_3(input_params.mesh_I_file);
	tallyFastVector.PrintIntersectionTables(input_params.intersection_table_full_file);
	std::cout << "finished!" << std::endl;

	timing_end   		= std::chrono::system_clock::now();
	timing_end_total	= std::chrono::system_clock::now();

	elapsed_seconds_export = timing_end-timing_start;
	elapsed_seconds_total  = timing_end_total - timing_start_total;

	double volume  = tallyFastVector.IntersectionTDS_3.get_volume();

	std::cout << " ---> Timing : " << std::endl;
	std::cout 	<< dtA.get_nb_of_cells() << " "
				<< dtB.get_nb_of_cells() << " "
				<< tallyFastVector.IntersectionTDS_3.get_nb_of_cells() << " "
				<< volume << " "
				<< elapsed_seconds_input.count() << " "
				<< elapsed_seconds_intersections.count() << " "
				<< elapsed_seconds_export.count() << " "
				<< elapsed_seconds_total.count() << std::endl;

	return 0;
}
