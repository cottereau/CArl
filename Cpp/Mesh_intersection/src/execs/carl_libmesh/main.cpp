#include "main.h"


//
//// Some functions to set up variables and co.
//libMesh::LinearImplicitSystem& add_elasticity(	libMesh::EquationSystems& input_systems,
//												libMesh::Order order = libMesh::FIRST,
//												libMesh::FEFamily family = libMesh::LAGRANGE)
//{
//	libMesh::ExplicitSystem& physical_variables =
//			input_systems.add_system<libMesh::ExplicitSystem> ("PhysicalConstants");
//
//	// Physical constants are set as constant, monomial
//	physical_variables.add_variable("E", libMesh::CONSTANT, libMesh::MONOMIAL);
//	physical_variables.add_variable("mu", libMesh::CONSTANT, libMesh::MONOMIAL);
//
//	libMesh::LinearImplicitSystem& elasticity_system =
//			input_systems.add_system<libMesh::LinearImplicitSystem> ("Elasticity");
//
//	elasticity_system.add_variable("u", order, family);
//	elasticity_system.add_variable("v", order, family);
//	elasticity_system.add_variable("w", order, family);
//
//	return elasticity_system;
//}
//
//// Some functions to set up variables and co.
//libMesh::LinearImplicitSystem& add_elasticity(	libMesh::EquationSystems& input_systems,
//						void fptr(	libMesh::EquationSystems& es,
//									const std::string& name),
//						libMesh::Order order = libMesh::FIRST,
//						libMesh::FEFamily family = libMesh::LAGRANGE)
//{
//	libMesh::ExplicitSystem& physical_variables =
//			input_systems.add_system<libMesh::ExplicitSystem> ("PhysicalConstants");
//
//	// Physical constants are set as constant, monomial
//	physical_variables.add_variable("E", libMesh::CONSTANT, libMesh::MONOMIAL);
//	physical_variables.add_variable("mu", libMesh::CONSTANT, libMesh::MONOMIAL);
//
//	libMesh::LinearImplicitSystem& elasticity_system =
//			input_systems.add_system<libMesh::LinearImplicitSystem> ("Elasticity");
//
//	elasticity_system.add_variable("u", order, family);
//	elasticity_system.add_variable("v", order, family);
//	elasticity_system.add_variable("w", order, family);
//
//	elasticity_system.attach_assemble_function(fptr);
//
//	return elasticity_system;
//}



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
	 *		TODO :  Implement the alpha mask reading
	 *
	 *		TODO :  Run examples with the multi-crystal case
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
	 * 		TODO :	Implement a parallel LATIN solver, capable of using any of
	 * 				the projection matrix assemblers
	 *
	 * ---------------------------------------------------------------------- */

	// - Start libmesh --------------------------------------------------------
	libMesh::LibMeshInit init (argc, argv);

	libMesh::PerfLog perf_log ("Main program");

	// - Displacement conditions ----------------------------------------------
	boundary_displacement x_max_BIG(0.5,0,0);
	boundary_displacement x_max_micro(0.5,0,0);

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

//	if ( command_line.search(1, "-o","--output") )
//	{
//		outputEXOFile_micro = command_line.next(outputEXOFile_micro);
//		outputEXOFile_micro = command_line.next(outputEXOFile_micro);
//	}
//	else
//	{
//		outputEXOFile_micro = "meshes/3D/output/carl_multi_crystal_test_micro.exo";
//	}

	// Set meshes
	libMesh::Mesh mesh_BIG(init.comm(), dim);
	libMesh::Mesh mesh_micro(init.comm(), dim);
	libMesh::Mesh mesh_restrict(init.comm(), dim);
	libMesh::Mesh mesh_inter(init.comm(), dim);

	std::string mesh_BIG_Filename;
	std::string mesh_micro_Filename;
	std::string mesh_restrict_Filename;
	std::string mesh_inter_Filename;

	std::string equivalence_table_restrict_A_Filename;
	std::string intersection_table_restrict_B_Filename;
	std::string intersection_table_I_Filename;

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

	bool using_same_mesh_restrict_A = false;


	if ( 	command_line.search(2, "--meshA", "-mA") &&
			command_line.search(2, "--meshB", "-mB") &&
			command_line.search(2, "--meshI", "-mI") &&
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

		command_line.search(1, "--tableB");
		intersection_table_restrict_B_Filename = command_line.next(intersection_table_restrict_B_Filename);

		command_line.search(1, "--tableI");
		intersection_table_I_Filename = command_line.next(intersection_table_I_Filename);

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
	}
	else
	{
		libmesh_error_msg("\n---> ERROR: Must give the three meshes and the intersection tables!\n");
	}

	// Generate the equation systems
	carl::coupled_system CoupledTest(mesh_micro.comm());


	libMesh::EquationSystems& equation_systems_inter =
					CoupledTest.add_inter_EquationSystem("InterSys", mesh_inter);

	// - Build the BIG system --------------------------------------------------

	perf_log.push("System initialization - BIG");

	libMesh::EquationSystems& equation_systems_BIG =
					CoupledTest.set_BIG_EquationSystem("BigSys", mesh_BIG);



	// [MACRO] Simple coupling
//	// TODO correct the coupling!!!
//	libMesh::LinearImplicitSystem& volume_BIG_system =
//			equation_systems_BIG.add_system<libMesh::LinearImplicitSystem> ("VolTest");
//
//	unsigned int sillyVar_BIG = volume_BIG_system.add_variable("SillyVar", libMesh::FIRST, libMesh::LAGRANGE);

	// [MACRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_BIG
										= add_elasticity_with_assemble(equation_systems_BIG,assemble_elasticity_heterogeneous);

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

	// [MICRO] Simple coupling
	// TODO correct the coupling!!!
//	libMesh::LinearImplicitSystem& volume_micro_system =
//			equation_systems_micro.add_system<libMesh::LinearImplicitSystem> ("VolTest");
//
//	unsigned int sillyVar_micro = volume_micro_system.add_variable("SillyVar", libMesh::FIRST, libMesh::LAGRANGE);

	// [MICRO] Set up the physical properties
	libMesh::LinearImplicitSystem& elasticity_system_micro
										= add_elasticity_with_assemble(equation_systems_micro,assemble_elasticity_heterogeneous);

	// [MICRO] Defining the boundaries with Dirichlet conditions
	set_x_displacement(elasticity_system_micro, x_max_micro,boundary_ids);

	// [MICRO] Build stress system
	libMesh::ExplicitSystem& stress_system_micro
										= add_stress(equation_systems_micro);

	equation_systems_micro.init();

	perf_log.pop("System initialization - micro");



	// - Build the restrict system --------------------------------------------------

	perf_log.push("System initialization - restrict");

	libMesh::EquationSystems& equation_systems_restrict =
					CoupledTest.add_restrict_EquationSystem("RestrictSys", mesh_restrict);

	// [RESTRICT] Simple coupling
	// TODO correct the coupling!!!
//	libMesh::LinearImplicitSystem& volume_restrict_system =
//			equation_systems_restrict.add_system<libMesh::LinearImplicitSystem> ("VolTest");
//
//	unsigned int sillyVar_restrict = volume_restrict_system.add_variable("SillyVar", libMesh::FIRST, libMesh::LAGRANGE);

	libMesh::LinearImplicitSystem& elasticity_system_restrict
										= add_elasticity(equation_systems_restrict);

	equation_systems_restrict.init();
	perf_log.pop("System initialization - restrict");



	// - Build the dummy inter system ------------------------------------------

	perf_log.push("System initialization - inter");
//	libMesh::LinearImplicitSystem& volume_inter_system =
//			equation_systems_inter.add_system<libMesh::LinearImplicitSystem> ("VolTest");
//
//	unsigned int sillyVar_inter = volume_inter_system.add_variable("SillyVar", libMesh::FIRST, libMesh::LAGRANGE);

	libMesh::LinearImplicitSystem& elasticity_system_inter
										= add_elasticity(equation_systems_inter);

	equation_systems_inter.init();
	perf_log.pop("System initialization - inter");

	perf_log.push("Build elasticity couplings");
	CoupledTest.assemble_coupling_elasticity_3D(	"BigSys","MicroSys",
													"InterSys","RestrictSys",
													equivalence_table_restrict_A,
													intersection_table_restrict_B,
													intersection_table_I,
													using_same_mesh_restrict_A);
	perf_log.pop("Build elasticity couplings");

	 perf_log.push("Physical properties - micro");
	 double BIG_E = 0;
	 double BIG_Mu = 0;
	 set_physical_properties(equation_systems_micro,physicalParamsFile,BIG_E,BIG_Mu);
	 perf_log.pop("Physical properties - micro");

	 perf_log.push("Physical properties - macro");
	 set_constant_physical_properties(equation_systems_BIG,BIG_E,BIG_Mu);
	 perf_log.pop("Physical properties - macro");

	// [MICRO] Solve
	perf_log.push("Solve - micro");
	elasticity_system_micro.solve();
	perf_log.pop("Solve - micro");

	// [MICRO] Calculate stress
	perf_log.push("Compute stress - micro");
	compute_stresses(equation_systems_micro);
	perf_log.pop("Compute stress - micro");

	// [MICRO] Export solution
#ifdef LIBMESH_HAVE_EXODUS_API

	libMesh::ExodusII_IO exo_io_micro(mesh_micro, /*single_precision=*/true);

	std::set<std::string> system_names_micro;
	system_names_micro.insert("Elasticity");
	exo_io_micro.write_equation_systems(outputEXOFile_micro.c_str(),equation_systems_micro,&system_names_micro);

	exo_io_micro.write_element_data(equation_systems_micro);

#endif

	// [MACRO] Solve
	perf_log.push("Solve - macro");
	elasticity_system_BIG.solve();
	perf_log.pop("Solve - macro");

	// [MACRO] Calculate stress
	perf_log.push("Compute stress - macro");
	compute_stresses(equation_systems_BIG);
	perf_log.pop("Compute stress - macro");

	// [MACRO] Export solution
#ifdef LIBMESH_HAVE_EXODUS_API

	libMesh::ExodusII_IO exo_io_BIG(mesh_BIG, /*single_precision=*/true);

	std::set<std::string> system_names_BIG;
	system_names_BIG.insert("Elasticity");
	exo_io_BIG.write_equation_systems(outputEXOFile_BIG.c_str(),equation_systems_BIG,&system_names_BIG);

	exo_io_BIG.write_element_data(equation_systems_BIG);

#endif

	CoupledTest.print_matrix_micro_info("MicroSys");
	CoupledTest.print_matrix_BIG_info("MicroSys");
	CoupledTest.print_matrix_restrict_info("MicroSys");

	libMesh::PetscMatrix<libMesh::Number>& LumpingTestInput = CoupledTest.get_restrict_coupling_matrix("MicroSys");
	libMesh::PetscMatrix<libMesh::Number> LumpingTestOutput(LumpingTestInput.comm());

	perf_log.push("Lumping");
	carl::lump_matrix(LumpingTestInput,LumpingTestOutput);
	perf_log.pop("Lumping");

	std::cout << std::endl << "| Testing the lumping " << std::endl;
	carl::print_matrix(LumpingTestOutput);

	return 0;
}
