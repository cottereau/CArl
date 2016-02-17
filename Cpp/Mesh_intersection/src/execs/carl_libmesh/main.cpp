#include "main.h"

struct carl_input_params
{
	std::string physical_params_file;

	std::string output_file_BIG;
	std::string output_file_micro;

	std::string mesh_BIG_file;
	std::string mesh_micro_file;
	std::string mesh_restrict_file;
	std::string mesh_inter_file;
	std::string mesh_weight_file;

	std::string equivalence_table_restrict_A_file;
	std::string intersection_table_restrict_B_file;
	std::string intersection_table_I_file;

	std::string weight_domain_idx_file;

	bool b_UseMeshAAsMediator;
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
		input_params.mesh_restrict_file = field_parser.next(input_params.mesh_restrict_file);

		field_parser.search(2, "--tableA", "MediatorEquivalenceTableFile");
		input_params.equivalence_table_restrict_A_file = field_parser.next(input_params.equivalence_table_restrict_A_file);

		input_params.b_UseMeshAAsMediator = false;
	}
	else if ( !field_parser.search(3, "--meshR", "-mR", "MeshMediator") && !field_parser.search(2, "--tableA", "MediatorEquivalenceTableFile") )
	{
		// Then use the mesh A as the mediator
		input_params.mesh_restrict_file = input_params.mesh_BIG_file;
		input_params.equivalence_table_restrict_A_file = "";

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
		input_params.intersection_table_restrict_B_file = field_parser.next(input_params.intersection_table_restrict_B_file);
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

void set_mesh_Gmsh(	libMesh::Mesh& mesh, const std::string& mesh_file,
					std::unordered_map<int,int>& mesh_NodeMap, std::unordered_map<int,int>& mesh_ElemMap)
{
	libMesh::GmshIO meshBuffer(mesh);
	meshBuffer.read(mesh_file);
	mesh.prepare_for_use();
	carl::create_mesh_map(mesh_file,mesh_NodeMap,mesh_ElemMap);
};

void set_mesh_Gmsh(	libMesh::Mesh& mesh, const std::string& mesh_file)
{
	libMesh::GmshIO meshBuffer(mesh);
	meshBuffer.read(mesh_file);
	mesh.prepare_for_use();
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
	 *		--->	CGAL -> MESH TREATMENT
	 *
	 * 		DONE : 	Build with CGAL a BIG mesh restrictor, to build the
	 * 				interface space mesh
	 *
	 * 				>>>	It's an external program. Adapted the input so that it
	 * 				   	can use a "restricted" mesh input
	 *
	 * 		DONE : 	Verify if libMesh conserves the indexes of the input meshes.
	 * 				If so, build list of intersection relations between BIG and
	 * 				micro
	 *
	 * 				->	Gmsh file reading method :
	 * 					https://libmesh.github.io/doxygen/gmsh__io_8C_source.html#l00145
	 *
	 * 				->	Huh, Gmsh has a "$PhysicalNames" section, might be handy
	 *
	 * 				->	libMesh creates an INTERNAL map between its node indices
	 * 					and the Gmsh indices, and doesn't save the element
	 * 					indices ...
	 *
	 * 				->	Are the micro meshes going to be from the Abaqus format?
	 * 					If so, this might be a bit tricky. Abaqus support is
	 * 					preliminary (since 2011 ...) and only has a read
	 * 					function. Still, its code could be useful for rolling
	 * 					out the parser for the rest of the data ...
	 *
	 * 				>>>	Created mesh reading functions that build maps between
	 * 				    the input files and libmesh indexes. It's a local
	 * 				    unordered map (hash table) so no easy data scatter, but
	 * 				    it should work for now.
	 *
	 * -------------------------------------------------------------------------
	 *
	 *		--->	COUPLING MATRIX AND LATIN METHOD
	 *
	 * 		DONE :	Build a (for now) serial coupling matrix assembler
	 *
	 * 				-> The only "serial" part of it is defined by the element
	 * 				   index tables, so it should be easy to paralellize it.
	 *
	 *		DONE :  Implement the alpha mask reading
	 *
	 *		DONE :  Run examples with the multi-crystal case
	 *
	 * 		DONE :	Look for how libMesh deals with the matrix repartition
	 * 				between the processors
	 *
	 * 		TODO :	Analyze the costs of inverting a matrix with PETSc and co.
	 * 				(hint: it's huge) and implement an "exact" projection matrix
	 * 				assembler
	 *
	 * 		DONE :	Implement a "lumping" projection matrix assembler
	 *
	 * 		DONE :	Implement a parallel LATIN solver, capable of using any of
	 * 				the projection matrix assemblers
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
	set_mesh_Gmsh(mesh_BIG,input_params.mesh_BIG_file,mesh_BIG_NodeMap,mesh_BIG_ElemMap);

	libMesh::Mesh mesh_micro(init.comm(), dim);
	std::unordered_map<int,int> mesh_micro_NodeMap;
	std::unordered_map<int,int> mesh_micro_ElemMap;
	set_mesh_Gmsh(mesh_micro,input_params.mesh_micro_file,mesh_micro_NodeMap,mesh_micro_ElemMap);

	libMesh::Mesh mesh_inter(init.comm(), dim);
	std::unordered_map<int,int> mesh_inter_NodeMap;
	std::unordered_map<int,int> mesh_inter_ElemMap;
	set_mesh_Gmsh(mesh_inter,input_params.mesh_inter_file,mesh_inter_NodeMap,mesh_inter_ElemMap);

	libMesh::Mesh mesh_weight(init.comm(), dim);
	set_mesh_Gmsh(mesh_weight,input_params.mesh_weight_file);

	// - Set mediator mesh and index tables -----------------------------------
	std::unordered_map<int,int> equivalence_table_restrict_A;
	std::vector<std::pair<int,int> > intersection_table_restrict_B;
	std::unordered_multimap<int,int> intersection_table_I;

	libMesh::Mesh mesh_restrict(init.comm(), dim);
	std::unordered_map<int,int> mesh_restrict_NodeMap;
	std::unordered_map<int,int> mesh_restrict_ElemMap;
	if ( !input_params.b_UseMeshAAsMediator )
	{
		set_mesh_Gmsh(mesh_restrict,input_params.mesh_restrict_file,mesh_restrict_NodeMap,mesh_restrict_ElemMap);

		carl::generate_intersection_tables_full(
				input_params.equivalence_table_restrict_A_file,
				input_params.intersection_table_restrict_B_file,
				input_params.intersection_table_I_file,
				mesh_restrict_ElemMap,
				mesh_micro_ElemMap,
				mesh_BIG_ElemMap,
				mesh_inter_ElemMap,
				equivalence_table_restrict_A,
				intersection_table_restrict_B,
				intersection_table_I);
	}
	else
	{
		mesh_restrict.copy_nodes_and_elements(mesh_BIG);
		input_params.mesh_restrict_file = input_params.mesh_BIG_file;
		mesh_restrict_NodeMap = mesh_BIG_NodeMap;
		mesh_restrict_ElemMap = mesh_BIG_ElemMap;

		carl::generate_intersection_tables_partial(
				input_params.intersection_table_restrict_B_file,
				input_params.intersection_table_I_file,
				mesh_restrict_ElemMap,
				mesh_micro_ElemMap,
				mesh_inter_ElemMap,
				intersection_table_restrict_B,
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
	for(libMesh::MeshBase::element_iterator itBegin = mesh_restrict.elements_begin();
			itBegin != mesh_restrict.elements_end(); ++itBegin)
	{
		silly_elem = *itBegin;
		vol += silly_elem->volume();
	}
	std::cout << "| Mesh restriction info :" << std::endl;
	std::cout << "|    filename     " << input_params.mesh_restrict_file << std::endl;
	std::cout << "|    n_elem       " << mesh_restrict.n_elem() << std::endl;
	std::cout << "|    n_nodes      " << mesh_restrict.n_nodes() << std::endl;
	std::cout << "|    n_subdomains " << mesh_restrict.n_subdomains() << std::endl;
	std::cout << "|    volume       " << vol << std::endl << std::endl;

	std::cout << "| Inter. table restrict / A :" << std::endl;
	std::cout << "|    filename     " << input_params.equivalence_table_restrict_A_file << std::endl << std::endl;

	std::cout << "| Inter. table restrict / B :" << std::endl;
	std::cout << "|    filename     " << input_params.intersection_table_restrict_B_file << std::endl << std::endl;

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

	// [MICRO] Add the weight function mesh
	carl::weight_parameter_function& sillyalphas = CoupledTest.add_alpha_mask("MicroSys",mesh_weight);
	CoupledTest.set_alpha_mask_parameters("MicroSys",domain_Idx_BIG,domain_Idx_micro[0],domain_Idx_coupling[0]);

	perf_log.pop("System initialization - micro");

	// - Build the restrict system --------------------------------------------------

	perf_log.push("System initialization - restrict");

	libMesh::EquationSystems& equation_systems_restrict =
					CoupledTest.add_restrict_EquationSystem("RestrictSys", mesh_restrict);

	libMesh::LinearImplicitSystem& elasticity_system_restrict
										= add_elasticity(equation_systems_restrict);

	equation_systems_restrict.init();

	perf_log.pop("System initialization - restrict");

	// - Build the dummy inter system ------------------------------------------

	perf_log.push("System initialization - inter");

	libMesh::LinearImplicitSystem& elasticity_system_inter
										= add_elasticity(equation_systems_inter);

	equation_systems_inter.init();

	perf_log.pop("System initialization - inter");

	perf_log.push("Physical properties - micro");
	double BIG_E = 0;
	double BIG_Mu = 0;
	double mean_distance = 0.2;
	double k_dA = 2.5;
	double k_dB = 2.5;
	double k_cA = 2.5;
	double k_cB = 2.5;
	double coupling_const = -1;
	set_physical_properties(equation_systems_micro,input_params.physical_params_file,BIG_E,BIG_Mu);
	coupling_const = eval_lambda_1(BIG_E,BIG_Mu);
	perf_log.pop("Physical properties - micro");

	perf_log.push("Physical properties - macro");
	set_constant_physical_properties(equation_systems_BIG,BIG_E,BIG_Mu);
	perf_log.pop("Physical properties - macro");

	perf_log.push("Build elasticity couplings");
	CoupledTest.set_coupling_parameters("MicroSys",coupling_const,mean_distance);

	CoupledTest.use_H1_coupling("MicroSys");

	CoupledTest.assemble_coupling_elasticity_3D(	"BigSys","MicroSys",
													"InterSys","RestrictSys",
													equivalence_table_restrict_A,
													intersection_table_restrict_B,
													intersection_table_I,
													input_params.b_UseMeshAAsMediator);
	perf_log.pop("Build elasticity couplings");

	std::cout << "| ---> Constants " << std::endl;
	std::cout << "| Macro :" << std::endl;
	std::cout << "|    E            : " << BIG_E << std::endl;
	std::cout << "|    Mu (lamba_2) : " << BIG_Mu << std::endl;
	std::cout << "|    lambda_1     : " << eval_lambda_1(BIG_E,BIG_Mu) << std::endl;
	std::cout << "| LATIN :" << std::endl;
	std::cout << "|    k_dA, k_dB   : " << k_dA << " " << k_dB << std::endl;
	std::cout << "|    k_cA, k_cB   : " << k_cA << " " << k_cB << std::endl;
	std::cout << "|    kappa        : " << coupling_const << std::endl;
	std::cout << "|    e            : " << mean_distance << std::endl;

	CoupledTest.print_matrix_micro_info("MicroSys");
	CoupledTest.print_matrix_BIG_info("MicroSys");
	CoupledTest.print_matrix_restrict_info("MicroSys");

//	carl::check_coupling_matrix(CoupledTest.get_micro_coupling_matrix("MicroSys"),
//										mesh_inter,
//										coupling_const ,
//										"Restrict - Micro : coupling check");
//
//	carl::check_coupling_matrix(CoupledTest.get_BIG_coupling_matrix("MicroSys"),
//										mesh_inter,
//										coupling_const ,
//										"Restrict - BIG : coupling check");
//
//	carl::check_coupling_matrix(CoupledTest.get_restrict_coupling_matrix("MicroSys"),
//										mesh_inter,
//										coupling_const ,
//										"Restrict - Restrict : coupling check");

	std::cout << std::endl << "| --> Testing the solver " << std::endl << std::endl;
	perf_log.push("Set up","LATIN Solver:");
	CoupledTest.set_LATIN_solver(	"MicroSys","Elasticity",
									assemble_elasticity_with_weight,
									assemble_elasticity_heterogeneous_with_weight,
									k_dA, k_dB, k_cA, k_cB);
	perf_log.pop("Set up","LATIN Solver:");


	// Solve !
	perf_log.push("Solve","LATIN Solver:");
	CoupledTest.solve_LATIN("MicroSys","Elasticity");
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
