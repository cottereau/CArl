#include "main.h"

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
	boundary_displacement x_max_BIG(0.5,0,0);
	boundary_id_cube boundary_ids;

	// - Get inputs and set variables
	GetPot command_line (argc, argv);

	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// Set constant parameters
	std::string physicalParamsFile;
	if ( command_line.search(2, "-p","--parameters") )
	{
		physicalParamsFile = command_line.next(physicalParamsFile);
	}
	else
	{
		std::cerr << "Need a physical parameter file!" << std::endl;
		return 0;
	}

	// Set output
	std::string outputEXOFile_micro = "meshes/3D/output/carl_multi_crystal_test_micro.exo";
	std::string outputEXOFile_BIG  = "meshes/3D/output/carl_multi_crystal_test_macro.exo";

	if ( command_line.search(2, "-oA","--outputA") )
	{
		outputEXOFile_BIG = command_line.next(outputEXOFile_BIG);
	}
	if ( command_line.search(2, "-oB","--outputB") )
	{
		outputEXOFile_micro = command_line.next(outputEXOFile_micro);
	}

	// Set meshes
	libMesh::Mesh mesh_BIG(init.comm(), dim);
	libMesh::Mesh mesh_micro(init.comm(), dim);
	libMesh::Mesh mesh_restrict(init.comm(), dim);
	libMesh::Mesh mesh_inter(init.comm(), dim);
	libMesh::Mesh mesh_weight(init.comm(), dim);

	std::string mesh_BIG_Filename;
	std::string mesh_micro_Filename;
	std::string mesh_restrict_Filename;
	std::string mesh_inter_Filename;
	std::string mesh_weight_Filename;

	std::string equivalence_table_restrict_A_Filename;
	std::string intersection_table_restrict_B_Filename;
	std::string intersection_table_I_Filename;

	std::string weight_domain_idx_Filename;

	std::string extra_restrict_info = "";

	// Read meshes ...

	std::unordered_map<int,int> mesh_BIG_NodeMap;
	std::unordered_map<int,int> mesh_BIG_ElemMap;

	std::unordered_map<int,int> mesh_micro_NodeMap;
	std::unordered_map<int,int> mesh_micro_ElemMap;

	std::unordered_map<int,int> mesh_restrict_NodeMap;
	std::unordered_map<int,int> mesh_restrict_ElemMap;

	std::unordered_map<int,int> mesh_inter_NodeMap;
	std::unordered_map<int,int> mesh_inter_ElemMap;

	std::unordered_map<int,int> equivalence_table_restrict_A;
	std::vector<std::pair<int,int> > intersection_table_restrict_B;
	std::unordered_multimap<int,int> intersection_table_I;

	int domain_Idx_BIG = -1;
	int nb_of_domain_Idx = 1;
	std::vector<int> domain_Idx_micro;
	std::vector<int> domain_Idx_coupling;

	bool using_same_mesh_restrict_A = false;

	if ( 	command_line.search(2, "--meshA", "-mA") &&
			command_line.search(2, "--meshB", "-mB") &&
			command_line.search(2, "--meshI", "-mI") &&
			command_line.search(2, "--meshWeight", "-mW") &&
			command_line.search(1, "--weightIdx") &&
			command_line.search(1, "--tableB") &&
			command_line.search(1, "--tableI"))
	{
		command_line.search(2, "--meshA", "-mA");
		mesh_BIG_Filename = command_line.next(mesh_BIG_Filename);
		libMesh::GmshIO meshBuffer_BIG(mesh_BIG);
		meshBuffer_BIG.read(mesh_BIG_Filename);
		mesh_BIG.prepare_for_use();
		carl::create_mesh_map(mesh_BIG_Filename,mesh_BIG_NodeMap,mesh_BIG_ElemMap);

		command_line.search(2, "--meshB", "-mB");
		mesh_micro_Filename = command_line.next(mesh_micro_Filename);
		libMesh::GmshIO meshBuffer_micro(mesh_micro);
		meshBuffer_micro.read(mesh_micro_Filename);
		mesh_micro.prepare_for_use();
		carl::create_mesh_map(mesh_micro_Filename,mesh_micro_NodeMap,mesh_micro_ElemMap);

		command_line.search(2, "--meshI", "-mI");
		mesh_inter_Filename = command_line.next(mesh_inter_Filename);
		libMesh::GmshIO meshBuffer_inter(mesh_inter);
		meshBuffer_inter.read(mesh_inter_Filename);
		mesh_inter.prepare_for_use();
		carl::create_mesh_map(mesh_inter_Filename,mesh_inter_NodeMap,mesh_inter_ElemMap);

		command_line.search(2, "--meshW", "-mW");
		mesh_weight_Filename = command_line.next(mesh_weight_Filename);
		libMesh::GmshIO meshBuffer_weight(mesh_weight);
		meshBuffer_weight.read(mesh_weight_Filename);
		mesh_weight.prepare_for_use();

		command_line.search(1, "--tableB");
		intersection_table_restrict_B_Filename = command_line.next(intersection_table_restrict_B_Filename);

		command_line.search(1, "--tableI");
		intersection_table_I_Filename = command_line.next(intersection_table_I_Filename);

		command_line.search(1, "--weightIdx");
		weight_domain_idx_Filename = command_line.next(weight_domain_idx_Filename);
		carl::set_weight_function_domain_idx(	weight_domain_idx_Filename,
												domain_Idx_BIG, nb_of_domain_Idx,
												domain_Idx_micro, domain_Idx_coupling
												);

		if (command_line.search(2, "--meshR", "-mR") && command_line.search(1, "--tableA"))
		{
			command_line.search(2, "--meshR", "-mR");
			mesh_restrict_Filename = command_line.next(mesh_restrict_Filename);
			libMesh::GmshIO meshBuffer_restrict(mesh_restrict);
			meshBuffer_restrict.read(mesh_restrict_Filename);
			mesh_restrict.prepare_for_use();
			carl::create_mesh_map(mesh_restrict_Filename,mesh_restrict_NodeMap,mesh_restrict_ElemMap);

			command_line.search(1, "--tableA");
			equivalence_table_restrict_A_Filename = command_line.next(equivalence_table_restrict_A_Filename);

			carl::generate_intersection_tables_full(	equivalence_table_restrict_A_Filename,
												intersection_table_restrict_B_Filename,
												intersection_table_I_Filename,
												mesh_restrict_ElemMap,
												mesh_micro_ElemMap,
												mesh_BIG_ElemMap,
												mesh_inter_ElemMap,
												equivalence_table_restrict_A,
												intersection_table_restrict_B,
												intersection_table_I);
		}
		else if(!command_line.search(2, "--meshR", "-mR") && !command_line.search(1, "--tableA"))
		{
			mesh_restrict.copy_nodes_and_elements(mesh_BIG);
			mesh_restrict_Filename = mesh_BIG_Filename;
			mesh_restrict_NodeMap = mesh_BIG_NodeMap;
			mesh_restrict_ElemMap = mesh_BIG_ElemMap;

			using_same_mesh_restrict_A = true;

			extra_restrict_info = " (Same as A)";
			equivalence_table_restrict_A_Filename = " --- ";

			carl::generate_intersection_tables_partial(	intersection_table_restrict_B_Filename,
													intersection_table_I_Filename,
													mesh_restrict_ElemMap,
													mesh_micro_ElemMap,
													mesh_inter_ElemMap,
													intersection_table_restrict_B,
													intersection_table_I);
		}
		else if(!command_line.search(2, "--meshR", "-mR"))
		{
			libmesh_error_msg("\n---> ERROR: Missing restrition mesh file");
		}
		else if(!command_line.search(1, "--tableA"))
		{
			libmesh_error_msg("\n---> ERROR: Missing restrition table file");
		}

		double vol = 0;
		libMesh::Elem* silly_elem;
		for(libMesh::MeshBase::element_iterator itBegin = mesh_BIG.elements_begin();
				itBegin != mesh_BIG.elements_end(); ++itBegin)
		{
			silly_elem = *itBegin;
			vol += silly_elem->volume();
		}
		std::cout << "| Mesh BIG info :" << std::endl;
		std::cout << "|    filename     " << mesh_BIG_Filename << std::endl;
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
		std::cout << "|    filename     " << mesh_micro_Filename << std::endl;
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
		std::cout << "|    filename     " << mesh_inter_Filename << std::endl;
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
		std::cout << "| Mesh restriction info :" << extra_restrict_info << std::endl;
		std::cout << "|    filename     " << mesh_restrict_Filename << std::endl;
		std::cout << "|    n_elem       " << mesh_restrict.n_elem() << std::endl;
		std::cout << "|    n_nodes      " << mesh_restrict.n_nodes() << std::endl;
		std::cout << "|    n_subdomains " << mesh_restrict.n_subdomains() << std::endl;
		std::cout << "|    volume       " << vol << std::endl << std::endl;

		std::cout << "| Inter. table restrict / A :" << std::endl;
		std::cout << "|    filename     " << equivalence_table_restrict_A_Filename << std::endl << std::endl;

		std::cout << "| Inter. table restrict / B :" << std::endl;
		std::cout << "|    filename     " << intersection_table_restrict_B_Filename << std::endl << std::endl;

		std::cout << "| Inter. table I  :" << std::endl;
		std::cout << "|    filename     " << intersection_table_I_Filename << std::endl << std::endl;

		std::cout << "| Mesh weight info :" << std::endl;
		std::cout << "|    filename     " << mesh_weight_Filename << std::endl;
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
	}
	else
	{
		libmesh_error_msg("\n---> ERROR: Must give the four meshes, the intersection tables and the weight domains files!\n");
	}

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
	set_x_displacement(elasticity_system_BIG, x_max_BIG,boundary_ids);

	// [MACRO] Build stress system --------------------------------------------------
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
	set_physical_properties(equation_systems_micro,physicalParamsFile,BIG_E,BIG_Mu);
	perf_log.pop("Physical properties - micro");

	perf_log.push("Physical properties - macro");
	set_constant_physical_properties(equation_systems_BIG,BIG_E,BIG_Mu);
	perf_log.pop("Physical properties - macro");

	perf_log.push("Build elasticity couplings");
	CoupledTest.assemble_coupling_elasticity_3D(	"BigSys","MicroSys",
													"InterSys","RestrictSys",
													equivalence_table_restrict_A,
													intersection_table_restrict_B,
													intersection_table_I,
													eval_lambda_1(BIG_E,BIG_Mu)/(mean_distance*mean_distance),
													using_same_mesh_restrict_A);
	perf_log.pop("Build elasticity couplings");

	CoupledTest.print_matrix_micro_info("MicroSys");
	CoupledTest.print_matrix_BIG_info("MicroSys");
	CoupledTest.print_matrix_restrict_info("MicroSys");

	std::cout << std::endl << "| --> Testing the solver " << std::endl << std::endl;
	perf_log.push("Set up","LATIN Solver:");
	CoupledTest.set_LATIN_solver("MicroSys","Elasticity",assemble_elasticity_with_weight,assemble_elasticity_heterogeneous_with_weight);
	perf_log.pop("Set up","LATIN Solver:");

	// Solve !
	perf_log.push("Solve","LATIN Solver:");
	CoupledTest.solve_LATIN("MicroSys","Elasticity");
	perf_log.pop("Solve","LATIN Solver:");

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
	exo_io_micro.write_equation_systems(outputEXOFile_micro.c_str(),equation_systems_micro,&system_names_micro);

	exo_io_micro.write_element_data(equation_systems_micro);

	libMesh::ExodusII_IO exo_io_BIG(mesh_BIG, /*single_precision=*/true);

	std::set<std::string> system_names_BIG;
	system_names_BIG.insert("Elasticity");
	exo_io_BIG.write_equation_systems(outputEXOFile_BIG.c_str(),equation_systems_BIG,&system_names_BIG);

	exo_io_BIG.write_element_data(equation_systems_BIG);

#endif

	return 0;
}
