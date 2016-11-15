#include "carl_libmesh_bricks.h"

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

struct carl_coupling_generation_input_params {
	std::string physical_params_file;

	std::string mesh_BIG_file;
	std::string mesh_micro_file;

	std::string mesh_restrict_BIG_file;
	std::string mesh_restrict_micro_file;

	std::string mesh_mediator_file;
	std::string mesh_inter_file;

	std::string mesh_weight_file;

	std::string equivalence_table_restrict_BIG_file;
	std::string equivalence_table_restrict_micro_file;
	std::string equivalence_table_mediator;

	std::string intersection_table_full;

	std::string weight_domain_idx_file;

	bool b_UseMesh_BIG_AsMediator;
	bool b_UseMesh_micro_AsMediator;
	bool b_UseMesh_extra_AsMediator;

	bool LATIN_b_UseRestartFiles;
	bool LATIN_b_PrintRestartFiles;

	double mean_distance;

	double k_dA;
	double k_dB;
	double k_cA;
	double k_cB;

	double LATIN_eps;
	int LATIN_conv_max;
	double LATIN_relax;

	std::string LATIN_convergence_output;
	std::string LATIN_restart_file_base;

	std::string output_file_BIG;
	std::string output_file_micro;

	std::string solver_type_string;
	carl::CoupledSolverType solver_type;
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

	if ( field_parser.search(3, "--meshWeight", "-mW", "MeshWeight") )
	{
		input_params.mesh_weight_file = field_parser.next(input_params.mesh_weight_file);
	}
	else
	{
		homemade_error_msg("Missing the weight mesh file!");
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

	if (field_parser.search(2, "--tableFullI", "FullIntersectionElementsTable")) {
		input_params.intersection_table_full = field_parser.next(
				input_params.intersection_table_full);
	} else {
		homemade_error_msg("Missing the full intersection elements file!");
	}

	// Set table files
	if( field_parser.search(2, "--weightIdx", "WeightIndexes") )
	{
		input_params.weight_domain_idx_file = field_parser.next(input_params.weight_domain_idx_file);
	}
	else
	{
		homemade_error_msg("Missing the weight value file!");
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
		input_params.equivalence_table_mediator = input_params.equivalence_table_restrict_BIG_file;
	}
	if(input_params.b_UseMesh_micro_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_micro_file;
		input_params.equivalence_table_mediator = input_params.equivalence_table_restrict_micro_file;
	}
	if(input_params.b_UseMesh_extra_AsMediator)
	{
//		if (field_parser.search(3, "--meshM", "-mM", "MeshMediator")) {
//			input_params.mesh_mediator_file = field_parser.next(
//					input_params.mesh_mediator_file);
//		} else {
//			homemade_error_msg("Missing the mediator mesh file!");
//		}
		libmesh_not_implemented_msg("Still implementing the external mesh case!");
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

	input_params.mean_distance = 0.2;

	if (field_parser.search(2, "--dist", "CouplingMeshScale")) {
		input_params.mean_distance = field_parser.next(
				input_params.mean_distance);
	}

	// Set coupling parameters
	bool bUseSameSearchCoeff = false;
	
	input_params.mean_distance = 0.2;

	input_params.k_dA = 2.5;
	input_params.k_dB = 2.5;
	input_params.k_cA = 2.5;
	input_params.k_cB = 2.5;

	if( field_parser.search(1,"LATINUseSameSearchParameters") )
	{
		bUseSameSearchCoeff = field_parser.next(bUseSameSearchCoeff);
	}

	if( bUseSameSearchCoeff )
	{
		if( field_parser.search(2, "-k","LATINk") )
		{
			input_params.k_dA = field_parser.next(input_params.k_dA);
			input_params.k_dB = input_params.k_dA;
			input_params.k_cA = input_params.k_dA;
			input_params.k_cB = input_params.k_dA;
		}
	}
	else
	{
		if( field_parser.search(2, "-kdA","LATINkDecoupledA") )
		{
			input_params.k_dA = field_parser.next(input_params.k_dA);
		}

		if( field_parser.search(2, "-kdB","LATINkDecoupledB") )
		{
			input_params.k_dB = field_parser.next(input_params.k_dB);
		}

		if( field_parser.search(2, "-kcA","LATINkCoupledA") )
		{
			input_params.k_cA = field_parser.next(input_params.k_cA);
		}

		if( field_parser.search(2, "-kcB","LATINkCoupledB") )
		{
			input_params.k_cB = field_parser.next(input_params.k_cB);
		}
	}

	if( field_parser.search(2, "--dist","LATINCouplingMeshScale") )
	{
		input_params.mean_distance = field_parser.next(input_params.mean_distance);
	}

	// Set LATIN parameters
	input_params.LATIN_b_UseRestartFiles = false;
	input_params.LATIN_b_PrintRestartFiles = false;
	input_params.LATIN_eps = 1E-2;
	input_params.LATIN_conv_max = 10000;
	input_params.LATIN_relax = 0.8;

	input_params.LATIN_convergence_output = "LATIN_convergence.dat";
	input_params.LATIN_restart_file_base = "LATIN_restart";

	if( field_parser.search(2, "--LATINeps","LATINEps") )
	{
		input_params.LATIN_eps = field_parser.next(input_params.LATIN_eps);
	}

	if( field_parser.search(2, "--LATINconv","LATINConvergneceLimit") )
	{
		input_params.LATIN_conv_max = field_parser.next(input_params.LATIN_conv_max);
	}

	if( field_parser.search(2, "--LATINrelax","LATINRelax") )
	{
		input_params.LATIN_relax = field_parser.next(input_params.LATIN_relax);
	}

	if( field_parser.search(2, "--LATINconvOutput","LATINConvergneceOutput") )
	{
		input_params.LATIN_convergence_output = field_parser.next(input_params.LATIN_convergence_output);
	}

	if( field_parser.search(1,"Use_restart_data") )
	{
		input_params.LATIN_b_UseRestartFiles = true;
	}
	if( field_parser.search(1,"Save_restart_data") )
	{
		input_params.LATIN_b_PrintRestartFiles = true;
	}

	if( field_parser.search(1,"LATINRestartDataFiles") )
	{
		input_params.LATIN_restart_file_base = field_parser.next(input_params.LATIN_restart_file_base);
	}

	// Set output files
	input_params.output_file_BIG = "meshes/3D/output/carl_multi_crystal_test_micro.exo";
	input_params.output_file_micro = "meshes/3D/output/carl_multi_crystal_test_macro.exo";
	if ( field_parser.search(3, "-oA","--outputA", "OutputEXOFileA") )
	{
		input_params.output_file_BIG = field_parser.next(input_params.output_file_BIG);
	}
	if ( field_parser.search(3, "-oB","--outputB", "OutputEXOFileB") )
	{
		input_params.output_file_micro = field_parser.next(input_params.output_file_micro);
	}

	input_params.solver_type = carl::LATIN_MODIFIED_STIFFNESS;
	if ( field_parser.search(1, "SolverType") )
	{
		input_params.solver_type_string = field_parser.next(input_params.solver_type_string);
		if(input_params.solver_type_string == "LATIN_Modified")
			input_params.solver_type = carl::LATIN_MODIFIED_STIFFNESS;
		if(input_params.solver_type_string == "LATIN_Original_Stiffness")
			input_params.solver_type = carl::LATIN_ORIGINAL_STIFFNESS;
	}
}
;

int main(int argc, char** argv) {

	// - Start libmesh ---------------------------------------------------------
	const bool MASTER_bPerfLog_carl_libmesh = true;
	libMesh::LibMeshInit init(argc, argv);

	libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);

	// - Displacement conditions -----------------------------------------------
	boundary_displacement x_max_BIG(1.0,0,0);
	boundary_displacement x_min_BIG(-1.0,0,0);
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

	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// Set up the communicator and the rank variables
	libMesh::Parallel::Communicator& WorldComm = init.comm();
	libMesh::Parallel::Communicator LocalComm;

	int rank = WorldComm.rank();
	int nodes = WorldComm.size();

	WorldComm.split(rank,rank,LocalComm);

	// - Read the meshes -------------------------------------------------------

	// - Parallelized meshes: A, B, mediator and weight

	perf_log.push("Meshes - Parallel","Read files:");
	libMesh::Mesh mesh_BIG(WorldComm, dim);
	mesh_BIG.read(input_params.mesh_BIG_file);
	mesh_BIG.prepare_for_use();

	libMesh::Mesh mesh_micro(WorldComm, dim);
	mesh_micro.read(input_params.mesh_micro_file);
	mesh_micro.prepare_for_use();

	libMesh::Mesh mesh_mediator(WorldComm, dim);
	mesh_mediator.allow_renumbering(false);
	mesh_mediator.read(input_params.mesh_mediator_file);
	mesh_mediator.prepare_for_use();

	libMesh::Mesh mesh_weight(WorldComm, dim);
	mesh_weight.allow_renumbering(false);
	mesh_weight.read(input_params.mesh_weight_file);
	mesh_weight.prepare_for_use();

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

	// - Semi-local meshes: restrict A and restrict B

	// -> libMesh's default mesh IS the SerialMesh, which creates a copy of
	//    itself on each processor, but partitions the iterators over each
	//    processor to allow the parallelization of the code. The usage of
	//    ParallelMesh is still a bit risky, and, performance-wise, is worse
	//    than SerialMesh for less than a few thousand processors.

	perf_log.push("Meshes - Serial","Read files:");
	libMesh::SerialMesh mesh_R_BIG(WorldComm, dim);
	mesh_R_BIG.allow_renumbering(false);
	mesh_R_BIG.read(input_params.mesh_restrict_BIG_file);
	mesh_R_BIG.prepare_for_use();

	libMesh::SerialMesh mesh_R_micro(WorldComm, dim);
	mesh_R_micro.allow_renumbering(false);
	mesh_R_micro.read(input_params.mesh_restrict_micro_file);
	mesh_R_micro.prepare_for_use();

//	// DEBUG - Test: print info per proc
//	{
//		std::ofstream mesh_info_ofstream;
//		mesh_info_ofstream.open("meshes/parallel_test/output/mesh_RA_" + std::to_string(rank) + "_info.txt");
//		mesh_R_BIG.print_info(mesh_info_ofstream);
//		mesh_info_ofstream.close();
//
//		WorldComm.barrier();
//
//		mesh_info_ofstream.open("meshes/parallel_test/output/mesh_RB_" + std::to_string(rank) + "_info.txt");
//		mesh_R_micro.print_info(mesh_info_ofstream);
//		mesh_info_ofstream.close();
//
//		WorldComm.barrier();
//
//		mesh_info_ofstream.open("meshes/parallel_test/output/mesh_mediator_" + std::to_string(rank) + "_info.txt");
//		mesh_mediator.print_info(mesh_info_ofstream);
//		mesh_info_ofstream.close();
//
//		std::ofstream mesh_data;
//		mesh_data.open("meshes/parallel_test/output/mesh_A_data_" + std::to_string(rank)  + ".dat");
//		libMesh::MeshBase::const_element_iterator       el     = mesh_BIG.active_local_elements_begin();
//		const libMesh::MeshBase::const_element_iterator end_el = mesh_BIG.active_local_elements_end();
//
//		for ( ; el != end_el; ++el)
//		{
//			const libMesh::Elem* elem = *el;
//			mesh_data << elem->id() << " " << elem->point(0) << " " << elem->point(1) << " " << elem->point(2)<< " " << elem->point(3)<< std::endl;
//		}
//		mesh_data.close();
//	}

	// - Local mesh: intersection mesh
	libMesh::Mesh mesh_inter(LocalComm, dim);
	mesh_inter.allow_renumbering(false);
	std::string local_inter_mesh_filename = input_params.mesh_inter_file + "_r_"
										+ std::to_string(rank) + "_n_" + std::to_string(nodes) + ".e";
	std::string local_inter_table_filename = input_params.intersection_table_full + "_r_"
										+ std::to_string(rank) + "_n_" + std::to_string(nodes) + "_inter_table_Full.dat";
	mesh_inter.read(local_inter_mesh_filename);
	mesh_inter.prepare_for_use();

	std::string global_inter_table_filename = input_params.intersection_table_full + "_stitched_inter_table_Full.dat";

	perf_log.pop("Meshes - Serial","Read files:");

	/*
	 * 		To do on the first proc
	 * 		- Read and build the equivalence tables between R_X and X, e_X - DONE
	 * 		- Read and build the intersection pairs table between A and B, p_AB - DONE
	 *		- Read and build the intersection indexes table, I_F - DONE
	 *
	 * 		Broadcast e_X, p_AB, I_F - DONE
	 *
	 * 		Build on each proc
	 * 		- The restricted intersection pairs table, p_R,AB- DONE
	 * 		- A local intersection indexes table, I_L - DONE
	 *
	 * 		Convert the pairs table to the libMesh indexing - DONE
	 */

	perf_log.push("Equivalence / intersection tables","Read files:");
	std::unordered_map<int,std::pair<int,int> > local_intersection_pairs_map;
	std::unordered_map<int,std::pair<int,int> > local_intersection_restricted_pairs_map;
	std::unordered_map<int,int> local_intersection_meshI_to_inter_map;

	std::unordered_map<int,int> equivalence_table_BIG_to_R_BIG;
	std::unordered_map<int,int> equivalence_table_micro_to_R_micro;
	std::unordered_map<int,int> equivalence_table_R_BIG_to_BIG;
	std::unordered_map<int,int> equivalence_table_R_micro_to_micro;

	//	Start by reading and broadcasting the equivalence tables
	carl::set_equivalence_tables(
			WorldComm,
			input_params.equivalence_table_restrict_BIG_file,
			input_params.equivalence_table_restrict_micro_file,

			equivalence_table_BIG_to_R_BIG,
			equivalence_table_micro_to_R_micro,
			equivalence_table_R_BIG_to_BIG,
			equivalence_table_R_micro_to_micro);

	if(input_params.b_UseMesh_BIG_AsMediator)
	{
		carl::set_local_intersection_tables(
				WorldComm,
				mesh_inter,
				local_inter_table_filename,
				input_params.equivalence_table_restrict_BIG_file,
				input_params.equivalence_table_restrict_micro_file,

				equivalence_table_BIG_to_R_BIG,
				equivalence_table_micro_to_R_micro,

				local_intersection_pairs_map,
				local_intersection_restricted_pairs_map,
				local_intersection_meshI_to_inter_map);
	}
	else if(input_params.b_UseMesh_micro_AsMediator)
	{
		carl::set_local_intersection_tables(
				WorldComm,
				mesh_inter,
				local_inter_table_filename,
				input_params.equivalence_table_restrict_micro_file,
				input_params.equivalence_table_restrict_BIG_file,

				equivalence_table_micro_to_R_micro,
				equivalence_table_BIG_to_R_BIG,

				local_intersection_pairs_map,
				local_intersection_restricted_pairs_map,
				local_intersection_meshI_to_inter_map);
	}

	std::unordered_multimap<int,int> inter_mediator_BIG;
	std::unordered_multimap<int,int> inter_mediator_micro;

	carl::set_global_mediator_system_intersection_lists(
			WorldComm,
			global_inter_table_filename,
			equivalence_table_BIG_to_R_BIG,
			equivalence_table_R_BIG_to_BIG,
			inter_mediator_BIG,
			inter_mediator_micro);

	perf_log.pop("Equivalence / intersection tables","Read files:");

	perf_log.push("Weight function domain","Read files:");

	// Set weight functions
	int domain_Idx_BIG = -1;
	int nb_of_domain_Idx = 1;
	std::vector<int> domain_Idx_micro;
	std::vector<int> domain_Idx_coupling;

	carl::set_weight_function_domain_idx(	input_params.weight_domain_idx_file,
											domain_Idx_BIG, nb_of_domain_Idx,
											domain_Idx_micro, domain_Idx_coupling
											);

	WorldComm.barrier();

	perf_log.pop("Weight function domain","Read files:");


	// - Generate the equation systems -----------------------------------------
	perf_log.push("Initialization","System initialization:");
	carl::coupled_system CoupledTest(WorldComm,input_params.solver_type);

	libMesh::EquationSystems& equation_systems_inter =
					CoupledTest.add_inter_EquationSystem("InterSys", mesh_inter);

	// Add the weight function mesh
	CoupledTest.add_alpha_mask("MicroSys",mesh_weight);
	CoupledTest.set_alpha_mask_parameters("MicroSys",domain_Idx_BIG,domain_Idx_micro[0],domain_Idx_coupling[0]);
	perf_log.pop("Initialization","System initialization:");

	// - Build the BIG system --------------------------------------------------

	perf_log.push("Macro system","System initialization:");
	libMesh::EquationSystems& equation_systems_BIG =
					CoupledTest.set_BIG_EquationSystem("BigSys", mesh_BIG);

	// [MACRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_BIG
										= add_elasticity(equation_systems_BIG);

	// [MACRO] Defining the boundaries with Dirichlet conditions
	//set_displaced_border_translation(elasticity_system_BIG, x_min_BIG,boundary_ids.MIN_X);
	set_clamped_border(elasticity_system_BIG, boundary_ids.MIN_X);

	// [MACRO] Build stress system
	libMesh::ExplicitSystem& stress_system_BIG
										= add_stress(equation_systems_BIG);

	equation_systems_BIG.init();

	perf_log.pop("Macro system","System initialization:");

	// - Build the micro system ------------------------------------------------

	perf_log.push("Micro system","System initialization:");

	libMesh::EquationSystems& equation_systems_micro =
					CoupledTest.add_micro_EquationSystem<libMesh::PetscMatrix<libMesh::Number> >("MicroSys", mesh_micro);

	// [MICRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_micro
										= add_elasticity(equation_systems_micro);

	// [MICRO] Defining the boundaries with Dirichlet conditions
	set_displaced_border_translation(elasticity_system_micro, x_max_BIG,boundary_ids.MAX_X);


	// [MICRO] Build stress system
	libMesh::ExplicitSystem& stress_system_micro
										= add_stress(equation_systems_micro);

	equation_systems_micro.init();
	perf_log.pop("Micro system","System initialization:");

	// - Build the RESTRICTED BIG system ---------------------------------------

	perf_log.push("RESTRICTED macro system","System initialization:");
	libMesh::EquationSystems& equation_systems_R_BIG =
					CoupledTest.set_Restricted_BIG_EquationSystem("BigSys", mesh_R_BIG);

	// [R. MACRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_R_BIG
										= add_elasticity(equation_systems_R_BIG);

	equation_systems_R_BIG.init();

	perf_log.pop("RESTRICTED macro system","System initialization:");

	// - Build the RESTRICTED micro system ------------------------------------------------

	perf_log.push("RESTRICTED micro system","System initialization:");

	libMesh::EquationSystems& equation_systems_R_micro =
					CoupledTest.add_Restricted_micro_EquationSystem("MicroSys", mesh_R_micro);

	// [R. MICRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_R_micro
										= add_elasticity(equation_systems_R_micro);

	equation_systems_R_micro.init();
	perf_log.pop("RESTRICTED micro system","System initialization:");

	// - Build the mediator system ---------------------------------------------

	perf_log.push("Mediator system","System initialization:");

	libMesh::EquationSystems& equation_systems_mediator =
					CoupledTest.add_mediator_EquationSystem("MediatorSys", mesh_mediator);

	libMesh::LinearImplicitSystem& elasticity_system_mediator
										= add_elasticity(equation_systems_mediator);

	equation_systems_mediator.init();

	WorldComm.barrier();

	perf_log.pop("Mediator system","System initialization:");

	// - Build the dummy inter system ------------------------------------------

	perf_log.push("Intersection system","System initialization:");

	libMesh::LinearImplicitSystem& elasticity_system_inter
										= add_elasticity(equation_systems_inter);

	equation_systems_inter.init();

	WorldComm.barrier();

	perf_log.pop("Intersection system","System initialization:");

	perf_log.push("Physical properties","System initialization:");
	double BIG_E = 0;
	double BIG_Mu = 0;

	double coupling_const = -1;
	std::ifstream phys_params_file(input_params.physical_params_file);
	carl::jump_lines(phys_params_file);
	phys_params_file >> BIG_E >> BIG_Mu;
	phys_params_file.close();

	set_constant_physical_properties(equation_systems_BIG,BIG_E,BIG_Mu);
	set_constant_physical_properties(equation_systems_micro,BIG_E,BIG_Mu);
	perf_log.pop("Physical properties","System initialization:");

	// - Set the coupling matrix -----------------------------------------------
	perf_log.push("Set coupling matrices");
	coupling_const = eval_lambda_1(BIG_E,BIG_Mu);
	CoupledTest.set_coupling_parameters("MicroSys",coupling_const,input_params.mean_distance);

	CoupledTest.use_H1_coupling("MicroSys");
	CoupledTest.assemble_coupling_elasticity_3D_parallel("BigSys","MicroSys",
			"InterSys","MediatorSys",
			mesh_R_BIG, mesh_R_micro,
			local_intersection_pairs_map,
			local_intersection_restricted_pairs_map,
			local_intersection_meshI_to_inter_map,
			inter_mediator_BIG,
			inter_mediator_micro);

	std::cout << std::endl;
	std::cout << "| ---> Constants " << std::endl;
	std::cout << "| Macro :" << std::endl;
	std::cout << "|    E            : " << BIG_E << std::endl;
	std::cout << "|    Mu (lamba_2) : " << BIG_Mu << std::endl;
	std::cout << "|    lambda_1     : " << eval_lambda_1(BIG_E,BIG_Mu) << std::endl;
	std::cout << "| LATIN :" << std::endl;
	std::cout << "|    k_dA, k_dB   : " << input_params.k_dA << " " << input_params.k_dB << std::endl;
	std::cout << "|    k_cA, k_cB   : " << input_params.k_cA << " " << input_params.k_cB << std::endl;
	std::cout << "|    kappa        : " << coupling_const << std::endl;
	std::cout << "|    e            : " << input_params.mean_distance << std::endl;
	std::cout << "|    restart file : " << input_params.LATIN_restart_file_base << "*" << std::endl;

	std::cout << std::endl << "| --> Testing the solver " << std::endl << std::endl;
	perf_log.push("Set up","LATIN Solver:");
	if(input_params.LATIN_b_PrintRestartFiles || input_params.LATIN_b_UseRestartFiles)
	{
		CoupledTest.set_restart(	input_params.LATIN_b_UseRestartFiles,
				input_params.LATIN_b_PrintRestartFiles,
				input_params.LATIN_restart_file_base);
	}

	CoupledTest.set_LATIN_solver(	"MicroSys","Elasticity",
									assemble_elasticity_with_weight,
									assemble_elasticity_heterogeneous_with_weight,
									input_params.k_dA, input_params.k_dB, input_params.k_cA, input_params.k_cB,
									input_params.LATIN_eps, input_params.LATIN_conv_max, input_params.LATIN_relax);
	perf_log.pop("Set up","LATIN Solver:");


	// Solve !
	perf_log.push("Solve","LATIN Solver:");
	CoupledTest.solve("MicroSys","Elasticity",input_params.LATIN_convergence_output);
	perf_log.pop("Solve","LATIN Solver:");

	// Calculate stress
	perf_log.push("Compute stress - micro","Output:");
	compute_stresses(equation_systems_micro);
	perf_log.pop("Compute stress - micro","Output:");

	perf_log.push("Compute stress - macro","Output:");
	compute_stresses(equation_systems_BIG);
	perf_log.pop("Compute stress - macro","Output:");

	// Export solution
#ifdef LIBMESH_HAVE_EXODUS_API
	perf_log.push("Save output","Output:");
	libMesh::ExodusII_IO exo_io_micro(mesh_micro, /*single_precision=*/true);

	std::set<std::string> system_names_micro;
	system_names_micro.insert("Elasticity");
	exo_io_micro.write_equation_systems(input_params.output_file_micro,equation_systems_micro,&system_names_micro);

	exo_io_micro.write_element_data(equation_systems_micro);

	libMesh::ExodusII_IO exo_io_BIG(mesh_BIG, /*single_precision=*/true);

	std::set<std::string> system_names_BIG;
	system_names_BIG.insert("Elasticity");
	exo_io_BIG.write_equation_systems(input_params.output_file_BIG,equation_systems_BIG,&system_names_BIG);

	exo_io_BIG.write_element_data(equation_systems_BIG);
	perf_log.pop("Save output","Output:");
#endif

	std::ofstream perf_log_file("perf_log_" + std::to_string(rank)  + ".txt");
	perf_log_file << perf_log.get_log();
	perf_log_file.close();
	return 0;
}
