#include "main.h"

struct carl_input_params
{
	std::string physical_params_file;

	std::string output_file_BIG;
	std::string output_file_micro;

	std::string mesh_BIG_file;
	std::string mesh_micro_file;
	std::string mesh_mediator_file;
	std::string mesh_inter_file;
	std::string mesh_weight_file;

	std::string equivalence_table_mediator_A_file;
	std::string intersection_table_mediator_B_file;
	std::string intersection_table_I_file;

	std::string weight_domain_idx_file;

	bool b_UseMeshAAsMediator;

	double mean_distance;

	double k_dA;
	double k_dB;
	double k_cA;
	double k_cB;

	double LATIN_eps;
	int LATIN_conv_max;
	double LATIN_relax;

	std::string LATIN_convergence_output;
};

void get_input_params(GetPot& field_parser, carl_input_params& input_params)
{
	// Set constant parameters
	if ( field_parser.search(3, "-p","--parameters","PhysicalParameters") )
	{
		input_params.physical_params_file = field_parser.next(input_params.physical_params_file);
	}
	else
	{
		homemade_error_msg("Missing the physical parameters file!");
	}

	// Set mesh files
	if ( field_parser.search(3, "--meshA", "-mA", "MeshA") )
	{
		input_params.mesh_BIG_file = field_parser.next(input_params.mesh_BIG_file);
	}
	else
	{
		homemade_error_msg("Missing the A mesh file!");
	}

	if ( field_parser.search(3, "--meshB", "-mB", "MeshB") )
	{
		input_params.mesh_micro_file = field_parser.next(input_params.mesh_micro_file);
	}
	else
	{
		homemade_error_msg("Missing the B mesh file!");
	}

	if ( field_parser.search(3, "--meshI", "-mI", "MeshInter") )
	{
		input_params.mesh_inter_file = field_parser.next(input_params.mesh_inter_file);
	}
	else
	{
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

	if ( field_parser.search(3, "--meshR", "-mR", "MeshMediator") && field_parser.search(2, "--tableA", "MediatorEquivalenceTableFile") )
	{
		// Then we have an external mediator mesh
		field_parser.search(3, "--meshR", "-mR", "MeshMediator");
		input_params.mesh_mediator_file = field_parser.next(input_params.mesh_mediator_file);

		field_parser.search(2, "--tableA", "MediatorEquivalenceTableFile");
		input_params.equivalence_table_mediator_A_file = field_parser.next(input_params.equivalence_table_mediator_A_file);

		input_params.b_UseMeshAAsMediator = false;
	}
	else if ( !field_parser.search(3, "--meshR", "-mR", "MeshMediator") && !field_parser.search(2, "--tableA", "MediatorEquivalenceTableFile") )
	{
		// Then use the mesh A as the mediator
		input_params.mesh_mediator_file = input_params.mesh_BIG_file;
		input_params.equivalence_table_mediator_A_file = "";

		input_params.b_UseMeshAAsMediator = true;
	}
	else if ( !field_parser.search(3, "--meshR", "-mR", "MeshMediator") )
	{
		homemade_error_msg("Missing the mediator mesh file!");
	}
	else if ( !field_parser.search(2, "--tableA", "MediatorEquivalenceTableFile") )
	{
		homemade_error_msg("Missing the mediator / A meshes equivalence file!");
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

	if( field_parser.search(2, "--tableB", "IntersectionPairsTable") )
	{
		input_params.intersection_table_mediator_B_file = field_parser.next(input_params.intersection_table_mediator_B_file);
	}
	else
	{
		homemade_error_msg("Missing the intersection pairs file!");
	}

	if(	field_parser.search(2, "--tableI", "IntersectionElementsTable") )
	{
		input_params.intersection_table_I_file = field_parser.next(input_params.intersection_table_I_file);
	}
	else
	{
		homemade_error_msg("Missing the intersection elements file!");
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
	input_params.LATIN_eps = 1E-2;
	input_params.LATIN_conv_max = 10000;
	input_params.LATIN_relax = 0.8;

	input_params.LATIN_convergence_output = "LATIN_convergence.dat";

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


};

int main (int argc, char** argv)
{
	/* To Do list --------------------------------------------------------------
	 *
	 * 		--->	LIBMESH DOCUMENTATION AND EXAMPLES
	 *
	 *		TODO :	Look at libMesh's implementation of coupled heat and fluid
	 *				examples
	 *
	 *				-> 	It's in the adjoint problem examples, and involves the
	 *				   	FEMSystem class. Take another look at Alvaro's articles
	 *
	 * -------------------------------------------------------------------------
	 *
	 *		--->	COUPLING MATRIX AND LATIN METHOD
	 *
	 * 		TODO :	Analyze the costs of inverting a matrix with PETSc and co.
	 * 				(hint: it's huge) and implement an "exact" projection matrix
	 * 				assembler
	 *
	 * ---------------------------------------------------------------------- */

	// - Start libmesh --------------------------------------------------------
	const bool MASTER_bPerfLog_carl_libmesh = true;
	libMesh::LibMeshInit init (argc, argv);

	libMesh::PerfLog perf_log ("Main program",MASTER_bPerfLog_carl_libmesh);

	// - Displacement conditions ----------------------------------------------
	boundary_displacement x_max_BIG(1.0,0,0);
	boundary_displacement x_min_BIG(-0.25,0,0);
	boundary_id_cube boundary_ids;

	// - Set up inputs
	GetPot command_line (argc, argv);
	GetPot field_parser;
	std::string input_filename;

	if( command_line.search(1, "--inputfile") )
	{
		input_filename = command_line.next(input_filename);
		field_parser.parse_input_file(input_filename,"#","\n"," \t\n");
	}
	else
	{
		field_parser = command_line;
	}

	carl_input_params input_params;
	get_input_params(field_parser, input_params);

	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// - Set meshes -----------------------------------------------------------
	libMesh::Mesh mesh_BIG(init.comm(), dim);
	std::unordered_map<int,int> mesh_BIG_NodeMap;
	std::unordered_map<int,int> mesh_BIG_ElemMap;
	carl::set_mesh_Gmsh(mesh_BIG,input_params.mesh_BIG_file,mesh_BIG_NodeMap,mesh_BIG_ElemMap);

	libMesh::Mesh mesh_micro(init.comm(), dim);
	std::unordered_map<int,int> mesh_micro_NodeMap;
	std::unordered_map<int,int> mesh_micro_ElemMap;
	carl::set_mesh_Gmsh(mesh_micro,input_params.mesh_micro_file,mesh_micro_NodeMap,mesh_micro_ElemMap);

	libMesh::Mesh mesh_inter(init.comm(), dim);
	std::unordered_map<int,int> mesh_inter_NodeMap;
	std::unordered_map<int,int> mesh_inter_ElemMap;
	carl::set_mesh_Gmsh(mesh_inter,input_params.mesh_inter_file,mesh_inter_NodeMap,mesh_inter_ElemMap);

	libMesh::Mesh mesh_weight(init.comm(), dim);
	carl::set_mesh_Gmsh(mesh_weight,input_params.mesh_weight_file);

	// - Set mediator mesh and index tables -----------------------------------
	std::unordered_map<int,int> equivalence_table_mediator_A;
	std::vector<std::pair<int,int> > intersection_table_mediator_B;
	std::unordered_multimap<int,int> intersection_table_I;

	libMesh::Mesh mesh_mediator(init.comm(), dim);
	std::unordered_map<int,int> mesh_mediator_NodeMap;
	std::unordered_map<int,int> mesh_mediator_ElemMap;
	if ( !input_params.b_UseMeshAAsMediator )
	{
		carl::set_mesh_Gmsh(mesh_mediator,input_params.mesh_mediator_file,mesh_mediator_NodeMap,mesh_mediator_ElemMap);

		carl::generate_intersection_tables_full(
				input_params.equivalence_table_mediator_A_file,
				input_params.intersection_table_mediator_B_file,
				input_params.intersection_table_I_file,
				mesh_mediator_ElemMap,
				mesh_micro_ElemMap,
				mesh_BIG_ElemMap,
				mesh_inter_ElemMap,
				equivalence_table_mediator_A,
				intersection_table_mediator_B,
				intersection_table_I);
	}
	else
	{
		mesh_mediator.copy_nodes_and_elements(mesh_BIG);
		input_params.mesh_mediator_file = input_params.mesh_BIG_file;
		mesh_mediator_NodeMap = mesh_BIG_NodeMap;
		mesh_mediator_ElemMap = mesh_BIG_ElemMap;

		carl::generate_intersection_tables_partial(
				input_params.intersection_table_mediator_B_file,
				input_params.intersection_table_I_file,
				mesh_mediator_ElemMap,
				mesh_micro_ElemMap,
				mesh_inter_ElemMap,
				intersection_table_mediator_B,
				intersection_table_I);
	}

	// Set weight functions
	int domain_Idx_BIG = -1;
	int nb_of_domain_Idx = 1;
	std::vector<int> domain_Idx_micro;
	std::vector<int> domain_Idx_coupling;

	carl::set_weight_function_domain_idx(	input_params.weight_domain_idx_file,
											domain_Idx_BIG, nb_of_domain_Idx,
											domain_Idx_micro, domain_Idx_coupling
											);

	// - Print info about the meshes and tables -------------------------------
	double vol = 0;
	libMesh::Elem* silly_elem;
	for(libMesh::MeshBase::element_iterator itBegin = mesh_BIG.elements_begin();
			itBegin != mesh_BIG.elements_end(); ++itBegin)
	{
		silly_elem = *itBegin;
		vol += silly_elem->volume();
	}
	std::cout << "| Mesh BIG info :" << std::endl;
	std::cout << "|    filename     " << input_params.mesh_BIG_file << std::endl;
	std::cout << "|    n_elem       " << mesh_BIG.n_elem() << std::endl;
	std::cout << "|    n_nodes      " << mesh_BIG.n_nodes() << std::endl;
	std::cout << "|    n_subdomains " << mesh_BIG.n_subdomains() << std::endl;
	std::cout << "|    volume       " << vol << std::endl << std::endl;

	vol = 0;
	for(libMesh::MeshBase::element_iterator itBegin = mesh_micro.elements_begin();
			itBegin != mesh_micro.elements_end(); ++itBegin)
	{
		silly_elem = *itBegin;
		vol += silly_elem->volume();
	}
	std::cout << "| Mesh micro info :" << std::endl;
	std::cout << "|    filename     " << input_params.mesh_micro_file << std::endl;
	std::cout << "|    n_elem       " << mesh_micro.n_elem() << std::endl;
	std::cout << "|    n_nodes      " << mesh_micro.n_nodes() << std::endl;
	std::cout << "|    n_subdomains " << mesh_micro.n_subdomains() << std::endl;
	std::cout << "|    volume       " << vol << std::endl << std::endl;

	vol = 0;
	for(libMesh::MeshBase::element_iterator itBegin = mesh_inter.elements_begin();
			itBegin != mesh_inter.elements_end(); ++itBegin)
	{
		silly_elem = *itBegin;
		vol += silly_elem->volume();
	}
	std::cout << "| Mesh inter info :" << std::endl;
	std::cout << "|    filename     " << input_params.mesh_inter_file << std::endl;
	std::cout << "|    n_elem       " << mesh_inter.n_elem() << std::endl;
	std::cout << "|    n_nodes      " << mesh_inter.n_nodes() << std::endl;
	std::cout << "|    n_subdomains " << mesh_inter.n_subdomains() << std::endl;
	std::cout << "|    volume       " << vol << std::endl << std::endl;

	vol = 0;
	for(libMesh::MeshBase::element_iterator itBegin = mesh_mediator.elements_begin();
			itBegin != mesh_mediator.elements_end(); ++itBegin)
	{
		silly_elem = *itBegin;
		vol += silly_elem->volume();
	}
	std::cout << "| Mesh mediatorion info :" << std::endl;
	std::cout << "|    filename     " << input_params.mesh_mediator_file << std::endl;
	std::cout << "|    n_elem       " << mesh_mediator.n_elem() << std::endl;
	std::cout << "|    n_nodes      " << mesh_mediator.n_nodes() << std::endl;
	std::cout << "|    n_subdomains " << mesh_mediator.n_subdomains() << std::endl;
	std::cout << "|    volume       " << vol << std::endl << std::endl;

	std::cout << "| Inter. table mediator / A :" << std::endl;
	std::cout << "|    filename     " << input_params.equivalence_table_mediator_A_file << std::endl << std::endl;

	std::cout << "| Inter. table mediator / B :" << std::endl;
	std::cout << "|    filename     " << input_params.intersection_table_mediator_B_file << std::endl << std::endl;

	std::cout << "| Inter. table I  :" << std::endl;
	std::cout << "|    filename     " << input_params.intersection_table_I_file << std::endl << std::endl;

	std::cout << "| Mesh weight info :" << std::endl;
	std::cout << "|    filename     " << input_params.mesh_weight_file << std::endl;
	std::cout << "|    n_elem       " << mesh_weight.n_elem() << std::endl;
	std::cout << "|    n_nodes      " << mesh_weight.n_nodes() << std::endl;
	std::cout << "|    n_subdomains " << mesh_weight.n_subdomains() << " " << nb_of_domain_Idx << std::endl;
	std::cout << "|    macro idx    " << domain_Idx_BIG << std::endl;
	std::cout << "|    micro idx    ";
	for(int iii = 0; iii < nb_of_domain_Idx; ++iii)
	{
		std::cout << domain_Idx_micro[iii] << " ";
	}
	std::cout << std::endl;
	std::cout << "|    coupling idx ";
	for(int iii = 0; iii < nb_of_domain_Idx; ++iii)
	{
		std::cout << domain_Idx_coupling[iii] << " ";
	}
	std::cout << std::endl << std::endl;

	// Generate the equation systems
	carl::coupled_system CoupledTest(mesh_micro.comm());

	libMesh::EquationSystems& equation_systems_inter =
					CoupledTest.add_inter_EquationSystem("InterSys", mesh_inter);


	// - Build the BIG system --------------------------------------------------

	perf_log.push("System initialization - BIG");

	libMesh::EquationSystems& equation_systems_BIG =
					CoupledTest.set_BIG_EquationSystem("BigSys", mesh_BIG);

	// [MACRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_BIG
										= add_elasticity(equation_systems_BIG);

	// [MACRO] Defining the boundaries with Dirichlet conditions
	set_displaced_border_translation(elasticity_system_BIG, x_max_BIG,boundary_ids.MAX_X);
//	set_displaced_border_translation(elasticity_system_BIG, x_min_BIG,boundary_ids.MIN_X);
	set_clamped_border(elasticity_system_BIG, boundary_ids.MIN_X);


	// [MACRO] Build stress system
	libMesh::ExplicitSystem& stress_system_BIG
										= add_stress(equation_systems_BIG);

	equation_systems_BIG.init();

	perf_log.pop("System initialization - BIG");

	// - Build the micro system ------------------------------------------------

	perf_log.push("System initialization - micro");

	libMesh::EquationSystems& equation_systems_micro =
					CoupledTest.add_micro_EquationSystem<libMesh::PetscMatrix<libMesh::Number> >("MicroSys", mesh_micro);

	// [MICRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_micro
										= add_elasticity(equation_systems_micro);

	// [MICRO] Build stress system
	libMesh::ExplicitSystem& stress_system_micro
										= add_stress(equation_systems_micro);

	equation_systems_micro.init();
	perf_log.pop("System initialization - micro");

	// [MICRO] Add the weight function mesh
	carl::weight_parameter_function& sillyalphas = CoupledTest.add_alpha_mask("MicroSys",mesh_weight);
	CoupledTest.set_alpha_mask_parameters("MicroSys",domain_Idx_BIG,domain_Idx_micro[0],domain_Idx_coupling[0]);

	// - Build the mediator system --------------------------------------------------

	perf_log.push("System initialization - mediator");

	libMesh::EquationSystems& equation_systems_mediator =
					CoupledTest.add_mediator_EquationSystem("MediatorSys", mesh_mediator);

	libMesh::LinearImplicitSystem& elasticity_system_mediator
										= add_elasticity(equation_systems_mediator);

	equation_systems_mediator.init();

	perf_log.pop("System initialization - mediator");

	// - Build the dummy inter system ------------------------------------------

	perf_log.push("System initialization - inter");

	libMesh::LinearImplicitSystem& elasticity_system_inter
										= add_elasticity(equation_systems_inter);

	equation_systems_inter.init();

	perf_log.pop("System initialization - inter");

	perf_log.push("Physical properties - micro");
	double BIG_E = 0;
	double BIG_Mu = 0;

//	double mean_distance = 0.2;
//
//	double k_dA = 2.5;
//	double k_dB = 2.5;
//	double k_cA = 2.5;
//	double k_cB = 2.5;

	double coupling_const = -1;
	set_physical_properties(equation_systems_micro,input_params.physical_params_file,BIG_E,BIG_Mu);
	coupling_const = eval_lambda_1(BIG_E,BIG_Mu);
	perf_log.pop("Physical properties - micro");

	perf_log.push("Physical properties - macro");
	set_constant_physical_properties(equation_systems_BIG,BIG_E,BIG_Mu);
	perf_log.pop("Physical properties - macro");

	perf_log.push("Build elasticity couplings");
	CoupledTest.set_coupling_parameters("MicroSys",coupling_const,input_params.mean_distance);

	CoupledTest.use_H1_coupling("MicroSys");

	CoupledTest.assemble_coupling_elasticity_3D(	"BigSys","MicroSys",
													"InterSys","MediatorSys",
													equivalence_table_mediator_A,
													intersection_table_mediator_B,
													intersection_table_I,
													input_params.b_UseMeshAAsMediator);
	perf_log.pop("Build elasticity couplings");

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

	CoupledTest.print_matrix_micro_info("MicroSys");
	CoupledTest.print_matrix_BIG_info("MicroSys");
	CoupledTest.print_matrix_mediator_info("MicroSys");

//	carl::check_coupling_matrix(CoupledTest.get_micro_coupling_matrix("MicroSys"),
//										mesh_inter,
//										coupling_const ,
//										"Mediator - Micro : coupling check");
//
//	carl::check_coupling_matrix(CoupledTest.get_BIG_coupling_matrix("MicroSys"),
//										mesh_inter,
//										coupling_const ,
//										"Mediator - BIG : coupling check");
//
//	carl::check_coupling_matrix(CoupledTest.get_mediator_coupling_matrix("MicroSys"),
//										mesh_inter,
//										coupling_const ,
//										"Mediator - Mediator : coupling check");

	std::cout << std::endl << "| --> Testing the solver " << std::endl << std::endl;
	perf_log.push("Set up","LATIN Solver:");
	CoupledTest.set_LATIN_solver(	"MicroSys","Elasticity",
									assemble_elasticity_with_weight,
									assemble_elasticity_heterogeneous_with_weight,
									input_params.k_dA, input_params.k_dB, input_params.k_cA, input_params.k_cB,
									input_params.LATIN_eps, input_params.LATIN_conv_max, input_params.LATIN_relax);
	perf_log.pop("Set up","LATIN Solver:");


	// Solve !
	perf_log.push("Solve","LATIN Solver:");
	CoupledTest.solve_LATIN("MicroSys","Elasticity",input_params.LATIN_convergence_output);

	perf_log.pop("Solve","LATIN Solver:");
//
//	// [MICRO] Solve
//	perf_log.push("Solve - micro");
//	elasticity_system_micro.assemble_before_solve = false;
//	elasticity_system_micro.solve();
//	perf_log.pop("Solve - micro");
//
//	// [MACRO] Solve
//	perf_log.push("Solve - macro");
//	elasticity_system_BIG.assemble_before_solve = false;
//	elasticity_system_BIG.solve();
//	perf_log.pop("Solve - macro");

	// Calculate stress
	perf_log.push("Compute stress - micro");
	compute_stresses(equation_systems_micro);
	perf_log.pop("Compute stress - micro");

	perf_log.push("Compute stress - macro");
	compute_stresses(equation_systems_BIG);
	perf_log.pop("Compute stress - macro");

	// Export solution
#ifdef LIBMESH_HAVE_EXODUS_API

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

#endif

	return 0;
}
