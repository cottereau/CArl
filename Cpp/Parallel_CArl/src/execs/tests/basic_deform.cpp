#include "basic_deform.h"

/*
 * 		Sketch of the second version of the program ----------------------------
 *
 * 		-> Read the meshes A, B and the coupling region file.
 *
 * 		-> Read the intersection mesh I.
 *
 * 		-> Read the mediator mesh, M.
 *
 * 		-> Read the mesh restrictions, R_A and R_B, and the equivalence tables
 * 		   between them and the original meshes, t_R_A->A and t_R_B->B.
 *
 * 		-> A, B, M and I are partitioned over all processors, while each one
 * 		   will have a copy of R_A and R_B.
 *
 * 		TODO : For now, each processor will read the whole mesh R_A and R_B!
 * 		       This is not practical for large systems!
 *
 * 		-> Run the coupling assemble program using the intersection mesh as the
 * 		   main loop index.
 *
 *		-> The reduced meshes will be used to calculate the values, while the
 *		   equivalence tables will be used to convert them to the full meshes.
 *
 *		-> Solve the system using the LATIN method. Since the matrices are
 *		   associated to A, B and M (which were properly partitioned), no
 *		   further steps are needed to paralellize the solver.
 */

struct carl_basic_generation_input_params {
	std::string physical_params_file;

	std::string mesh_BIG_file;

	std::string output_file_BIG;
};

void get_input_params(GetPot& field_parser,
		carl_basic_generation_input_params& input_params) {

	// Set mesh files
	if (field_parser.search(3, "--meshA", "-mA", "MeshA")) {
		input_params.mesh_BIG_file = field_parser.next(
				input_params.mesh_BIG_file);
	} else {
		homemade_error_msg("Missing the A mesh file!");
	}

	// Set constant parameters
	if ( field_parser.search(3, "-p","--parameters","PhysicalParameters") )
	{
		input_params.physical_params_file = field_parser.next(input_params.physical_params_file);
	}
	else
	{
		homemade_error_msg("Missing the physical parameters file!");
	}

	// Set output files
	input_params.output_file_BIG = "meshes/3D/output/carl_multi_crystal_test_micro.exo";
	if ( field_parser.search(3, "-oA","--outputA", "OutputEXOFileA") )
	{
		input_params.output_file_BIG = field_parser.next(input_params.output_file_BIG);
	}
}
;

int main(int argc, char** argv) {

	// - Start libmesh ---------------------------------------------------------
	const bool MASTER_bPerfLog_carl_libmesh = true;
	libMesh::LibMeshInit init(argc, argv);

	libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);

	// - Displacement conditions -----------------------------------------------
	boundary_displacement y_max_BIG(0,1.0,0);
	boundary_displacement y_min_BIG(-0.25,0,0);
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

	carl_basic_generation_input_params input_params;
	get_input_params(field_parser, input_params);

	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// Set up the communicator and the rank variables
	libMesh::Parallel::Communicator& WorldComm = init.comm();

	int rank = WorldComm.rank();

//	WorldComm.split(rank,rank,LocalComm);

	// - Read the meshes -------------------------------------------------------

	// - Parallelized meshes: A, B, intersection, mediator and weight

	perf_log.push("Meshes - Parallel","Read files:");
	libMesh::Mesh mesh_BIG(WorldComm, dim);
	mesh_BIG.read(input_params.mesh_BIG_file);
	mesh_BIG.prepare_for_use();

//	// DEBUG - Test: print info per proc
//	{
//		std::ofstream mesh_info_ofstream;
//		mesh_info_ofstream.open("meshes/parallel_test/output/mesh_A_" + std::to_string(rank) + "_info.txt");
//		mesh_BIG.print_info(mesh_info_ofstream);
//		mesh_info_ofstream.close();
//
//		WorldComm.barrier();
//
//		mesh_info_ofstream.open("meshes/parallel_test/output/mesh_B_" + std::to_string(rank) + "_info.txt");
//		mesh_micro.print_info(mesh_info_ofstream);
//		mesh_info_ofstream.close();
//
//		WorldComm.barrier();
//
//		mesh_info_ofstream.open("meshes/parallel_test/output/mesh_inter_" + std::to_string(rank) + "_info.txt");
//		mesh_inter.print_info(mesh_info_ofstream);
//		mesh_info_ofstream.close();
//
//		WorldComm.barrier();
//	}

	perf_log.pop("Meshes - Parallel","Read files:");

	// - Build the BIG system --------------------------------------------------

	perf_log.push("Macro system","System initialization:");
	libMesh::EquationSystems equation_systems_BIG(mesh_BIG);

	// [MACRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_BIG
										= add_elasticity(equation_systems_BIG);

	elasticity_system_BIG.attach_assemble_function(assemble_elasticity);

	// [MACRO] Defining the boundaries with Dirichlet conditions
	set_displaced_border_translation(elasticity_system_BIG, y_max_BIG,boundary_ids.MAX_Y);
	set_clamped_border(elasticity_system_BIG, boundary_ids.MIN_Y);

	// [MACRO] Build stress system
	libMesh::ExplicitSystem& stress_system_BIG
										= add_stress(equation_systems_BIG);

	equation_systems_BIG.init();

	equation_systems_BIG.print_info();

	perf_log.pop("Macro system","System initialization:");

	perf_log.push("Physical properties","System initialization:");
	double BIG_E = 0;
	double BIG_Mu = 0;

	std::ifstream phys_params_file(input_params.physical_params_file);
	phys_params_file >> BIG_E >> BIG_Mu;
	phys_params_file.close();

	set_constant_physical_properties(equation_systems_BIG,BIG_E,BIG_Mu);

//	set_physical_properties(equation_systems_BIG,input_params.physical_params_file,BIG_E,BIG_Mu);
	perf_log.pop("Physical properties","System initialization:");

	// - Set the coupling matrix -----------------------------------------------
	std::cout << std::endl;
	std::cout << "| ---> Constants " << std::endl;
	std::cout << "| Macro :" << std::endl;
	std::cout << "|    E            : " << BIG_E << std::endl;
	std::cout << "|    Mu (lamba_2) : " << BIG_Mu << std::endl;
	std::cout << "|    lambda_1     : " << eval_lambda_1(BIG_E,BIG_Mu) << std::endl;

	// Solve !
	perf_log.push("Solve");
//	elasticity_system_BIG.solve();
	elasticity_system_BIG.assemble();
	elasticity_system_BIG.matrix->close();
	elasticity_system_BIG.rhs->close();

	libMesh::PetscMatrix<libMesh::Number> * sys_matrix =
			libMesh::libmesh_cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(elasticity_system_BIG.matrix);

	libMesh::PetscVector<libMesh::Number> * sys_vec =
			libMesh::libmesh_cast_ptr<libMesh::PetscVector<libMesh::Number> * >(elasticity_system_BIG.rhs);

	carl::base_CG_solver CG_test(WorldComm);
	CG_test.set_solver_matrix(*sys_matrix);
	CG_test.set_system_rhs(*sys_vec);

	CG_test.solve();

	*(elasticity_system_BIG.solution) = CG_test.get_solution();
	elasticity_system_BIG.solution->close();
	elasticity_system_BIG.update();
	perf_log.pop("Solve");

	perf_log.push("Compute stress - macro","Output:");
	compute_stresses(equation_systems_BIG);
	perf_log.pop("Compute stress - macro","Output:");

	// Export solution
#ifdef LIBMESH_HAVE_EXODUS_API
	perf_log.push("Save output","Output:");
	libMesh::ExodusII_IO exo_io_BIG(mesh_BIG, /*single_precision=*/true);

	std::set<std::string> system_names_BIG;
	system_names_BIG.insert("Elasticity");
	exo_io_BIG.write_equation_systems(input_params.output_file_BIG,equation_systems_BIG,&system_names_BIG);

	exo_io_BIG.write_element_data(equation_systems_BIG);
	perf_log.pop("Save output","Output:");
#endif



//	std::ofstream perf_log_file("meshes/parallel_test/output/perf_log_" + std::to_string(rank)  + ".txt");
//	perf_log_file << perf_log.get_log();
//	perf_log_file.close();
	return 0;
}
