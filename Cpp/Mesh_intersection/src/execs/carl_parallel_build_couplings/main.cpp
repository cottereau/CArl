#include "main.h"

/*
 * 		Sketch of the first version of the program -----------------------------
 *
 * 		-> Read the meshes A, B and the coupling region file.
 *
 * 		-> Read the intersection mesh I.
 *
 * 		-> Read the mediator mesh, M. (!)
 *
 * 		-> Read the mesh restrictions, R_A and R_B, and the equivalence tables
 * 		   between them and the original meshes, t_R_A->A and t_R_B->B.
 *
 * 		-> Copy R_A, R_B and I to all processors. (!)
 *
 * 		-> Run the coupling assemble program using the mediator space as the
 * 		   main loop index.
 *
 * 		-> Use the equivalence tables to fill the coupling matrices.
 *
 * 		-> Export them using PETSc MatView methods.
 *
 * 		Optimizations ----------------------------------------------------------
 *
 * 		-> Optimize the partitioning of M using MeshTools: build a weighting
 * 		   vector for the intersections, w_I, convert it to a w_M weighting, and
 * 		   partition M.
 *
 * 		-> Partition I by hand using this information.
 *
 * 		-> Partition R_A and R_B using this information.
 *
 * 		Functions I'll need ----------------------------------------------------
 *
 * 		TODO :	read meshes A, B and I in a way to make a single copy on each
 * 				processor - watch out for collisions!
 *
 * 		TODO :	read the intersection and the restriction equivalence tables.
 * 				Divide them on processor-by-processor tables.
 *
 * 		TODO :	adapt the coupling assemble methods in a way such that the main
 * 				loop is ran over the mediator mesh, and not the intersection
 * 				mesh.
 *
 * 		TODO :	export the matrix using Petsc interfaces.
 *
 * 		TODO :	build a matrix import function.
 *
 * 		Things to watch out for ------------------------------------------------
 *
 * 		-> Matrix renumerotations during the prepare_for_use step:
 * 		   This is done before the partitioning, and, according to the
 * 		   documentation, is done to guarantee that the elements are organized
 * 		   in contiguous blocks. From this, and from the current tests, it
 * 		   SHOULD not be a problem ... but still, watch out ...
 */

struct carl_coupling_generation_input_params {
	std::string physical_params_file;

	std::string mesh_BIG_file;
	std::string mesh_micro_file;

	std::string mesh_restrict_BIG_file;
	std::string mesh_restrict_micro_file;

	std::string mesh_mediator_file;
	std::string mesh_inter_file;

	std::string equivalence_table_restrict_BIG_file;
	std::string equivalence_table_restrict_micro_file;

	std::string intersection_table_BIG_micro_file;
	std::string intersection_table_I_file;

	bool b_UseMesh_BIG_AsMediator;
	bool b_UseMesh_micro_AsMediator;
	bool b_UseMesh_extra_AsMediator;

	double meanE;
	double meanMu;
	double coupling_const;
	double mean_distance;

	std::string output_coupling_matrix_med_BIG;
	std::string output_coupling_matrix_med_micro;
};

void get_input_params(GetPot& field_parser,
		carl_coupling_generation_input_params& input_params) {

	// Set mesh files
	if (field_parser.search(3, "--meshA", "-mA", "MeshA")) {
		input_params.mesh_BIG_file = field_parser.next(
				input_params.mesh_BIG_file);
	} else {
		homemade_error_msg("Missing the A mesh file!");
	}

	if (field_parser.search(3, "--meshB", "-mB", "MeshB")) {
		input_params.mesh_micro_file = field_parser.next(
				input_params.mesh_micro_file);
	} else {
		homemade_error_msg("Missing the B mesh file!");
	}

	if (field_parser.search(3, "--meshAR", "-mAR", "Mesh_A_Restriction")) {
		input_params.mesh_restrict_BIG_file = field_parser.next(
				input_params.mesh_restrict_BIG_file);
	} else {
		homemade_error_msg("Missing the restricted A mesh file!");
	}

	if (field_parser.search(3, "--meshBR", "-mBR", "Mesh_B_Restriction")) {
		input_params.mesh_restrict_micro_file = field_parser.next(
				input_params.mesh_restrict_micro_file);
	} else {
		homemade_error_msg("Missing the restricted B mesh file!");
	}

	if (field_parser.search(3, "--meshI", "-mI", "MeshInter")) {
		input_params.mesh_inter_file = field_parser.next(
				input_params.mesh_inter_file);
	} else {
		homemade_error_msg("Missing the intersection mesh file!");
	}

	// Set the equivalence and intersection tables
	if (field_parser.search(2, "--tableRA", "Mesh_A_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_BIG_file = field_parser.next(
				input_params.equivalence_table_restrict_BIG_file);
	} else {
		homemade_error_msg("Missing the equivalence table for the mesh A!");
	}

	if (field_parser.search(2, "--tableRB", "Mesh_B_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_micro_file = field_parser.next(
				input_params.equivalence_table_restrict_micro_file);
	} else {
		homemade_error_msg("Missing the equivalence table for the mesh B!");
	}

	if (field_parser.search(2, "--tablePairs", "IntersectionPairsTable")) {
		input_params.intersection_table_BIG_micro_file = field_parser.next(
				input_params.intersection_table_BIG_micro_file);
	} else {
		homemade_error_msg("Missing the intersection pairs file!");
	}

	if (field_parser.search(2, "--tableI", "IntersectionElementsTable")) {
		input_params.intersection_table_I_file = field_parser.next(
				input_params.intersection_table_I_file);
	} else {
		homemade_error_msg("Missing the intersection elements file!");
	}

	// Set the mediator mesh
	input_params.b_UseMesh_BIG_AsMediator = false;
	input_params.b_UseMesh_micro_AsMediator = false;
	input_params.b_UseMesh_extra_AsMediator = false;
	if (field_parser.search(1,"Use_A_AsMediator"))
	{
		input_params.b_UseMesh_BIG_AsMediator = true;
	}
	if (field_parser.search(1,"Use_B_AsMediator"))
	{
		input_params.b_UseMesh_micro_AsMediator = true;
	}
	if (field_parser.search(1,"Use_extra_AsMediator"))
	{
		input_params.b_UseMesh_extra_AsMediator = true;
	}
	if(input_params.b_UseMesh_BIG_AsMediator
			+ input_params.b_UseMesh_micro_AsMediator
			+ input_params.b_UseMesh_extra_AsMediator > 1)
	{
		homemade_error_msg("Choose only one mesh as mediator!");
	}

	if(input_params.b_UseMesh_BIG_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_BIG_file;
	}
	if(input_params.b_UseMesh_micro_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_micro_file;
	}
	if(input_params.b_UseMesh_extra_AsMediator)
	{
		if (field_parser.search(3, "--meshM", "-mM", "MeshMediator")) {
			input_params.mesh_mediator_file = field_parser.next(
					input_params.mesh_mediator_file);
		} else {
			homemade_error_msg("Missing the mediator mesh file!");
		}
	}

	// Set constant parameters
	if (field_parser.search(3, "-p", "--parameters", "PhysicalParameters")) {
		input_params.physical_params_file = field_parser.next(
				input_params.physical_params_file);
		calculate_average_params(input_params.physical_params_file,
				input_params.meanE,input_params.meanMu);
		input_params.coupling_const = eval_lambda_1(input_params.meanE,
				input_params.meanMu);
	}
	else if (field_parser.search(2, "--meanE", "MeanE") &&
			field_parser.search(2, "--meanMu", "MeanMu")) {
		field_parser.search(2, "--meanE", "MeanE");
		input_params.meanE = field_parser.next(input_params.meanE);
		field_parser.search(2, "--meanMu", "MeanMu");
		input_params.meanMu = field_parser.next(input_params.meanMu);
		input_params.coupling_const = eval_lambda_1(input_params.meanE,
				input_params.meanMu);
	}
	else if (field_parser.search(2, "--coupling", "CouplingConstant"))
	{
		input_params.coupling_const = field_parser.next(
				input_params.coupling_const);
	}
	else
	{
		homemade_error_msg("Missing the coupling constant");
	}

	input_params.mean_distance = 0.2;

	if (field_parser.search(2, "--dist", "CouplingMeshScale")) {
		input_params.mean_distance = field_parser.next(
				input_params.mean_distance);
	}

	// Set output files
	input_params.output_coupling_matrix_med_BIG =
			"meshes/parallel_test/output/coupling_matrix_mediator_A.dat";
	input_params.output_coupling_matrix_med_micro =
			"meshes/parallel_test/output/coupling_matrix_mediator_B.dat";
	if (field_parser.search(3, "-oA", "--outputA", "OutputMatrixA")) {
		input_params.output_coupling_matrix_med_BIG = field_parser.next(
				input_params.output_coupling_matrix_med_BIG);
	}
	if (field_parser.search(3, "-oB", "--outputB", "OutputMatrixB")) {
		input_params.output_coupling_matrix_med_micro = field_parser.next(
				input_params.output_coupling_matrix_med_micro);
	}

}
;

int main(int argc, char** argv) {

	// - Start libmesh --------------------------------------------------------
	const bool MASTER_bPerfLog_carl_libmesh = true;
	libMesh::LibMeshInit init(argc, argv);

	libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);

	// - Displacement conditions ----------------------------------------------
	boundary_displacement x_max_BIG(1.0, 0, 0);
	boundary_displacement x_min_BIG(-0.25, 0, 0);
	boundary_id_cube boundary_ids;

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

	carl_coupling_generation_input_params input_params;
	get_input_params(field_parser, input_params);

	return 0;
}
