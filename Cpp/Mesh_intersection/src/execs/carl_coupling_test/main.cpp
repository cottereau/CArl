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

	// - Displacement conditions ----------------------------------------------
	double x_max_x_displ = 0.5;
	double x_max_y_displ = 0;
	double x_max_z_displ = 0;

	// - Start libmesh --------------------------------------------------------
	libMesh::LibMeshInit init (argc, argv);

	// - Get inputs and set variables
	GetPot command_line (argc, argv);

	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

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

//	// Read mesh ...
//	int BOUNDARY_ID_MIN_Z = 0;
//	int BOUNDARY_ID_MIN_Y = 1;
//	int BOUNDARY_ID_MAX_X = 2;
//	int BOUNDARY_ID_MAX_Y = 3;
//	int BOUNDARY_ID_MIN_X = 4;
//	int BOUNDARY_ID_MAX_Z = 5;

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
//
//		++BOUNDARY_ID_MIN_Z;
//		++BOUNDARY_ID_MIN_Y;
//		++BOUNDARY_ID_MAX_X;
//		++BOUNDARY_ID_MAX_Y;
//		++BOUNDARY_ID_MIN_X;
//		++BOUNDARY_ID_MAX_Z;

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

	// Associate the boundary conditions -> no boundary conditions!

	// Generate the equation systems
	carl::coupled_system CoupledTest(mesh_micro.comm());

	libMesh::EquationSystems& equation_systems_BIG =
					CoupledTest.set_BIG_EquationSystem("BigSys", mesh_BIG);
	libMesh::EquationSystems& equation_systems_micro =
					CoupledTest.add_micro_EquationSystem<libMesh::PetscMatrix<libMesh::Number> >("MicroSys", mesh_micro);
	libMesh::EquationSystems& equation_systems_restrict =
					CoupledTest.add_restrict_EquationSystem("RestrictSys", mesh_restrict);
	libMesh::EquationSystems& equation_systems_inter =
					CoupledTest.add_inter_EquationSystem("InterSys", mesh_inter);

	// - Build the BIG system --------------------------------------------------

	libMesh::ImplicitSystem& elasticity_system_BIG
										= add_elasticity_with_assemble(equation_systems_BIG,assemble_elasticity);

	equation_systems_BIG.init();

	// - Build the micro system ------------------------------------------------

	libMesh::ImplicitSystem& elasticity_system_micro
										= add_elasticity_with_assemble(equation_systems_micro,assemble_elasticity);

	equation_systems_micro.init();

	// - Build the restrict system --------------------------------------------------

	libMesh::ImplicitSystem& elasticity_system_restrict
										= add_elasticity(equation_systems_restrict);

	equation_systems_restrict.init();

	// - Build the dummy inter system ------------------------------------------

	libMesh::ImplicitSystem& elasticity_system_inter
										= add_elasticity(equation_systems_inter);

	equation_systems_inter.init();

	CoupledTest.assemble_coupling_elasticity_3D(	"BigSys","MicroSys",
													"InterSys","RestrictSys",
													equivalence_table_restrict_A,
													intersection_table_restrict_B,
													intersection_table_I,
													using_same_mesh_restrict_A);

	CoupledTest.print_matrix_micro_info("MicroSys");
	CoupledTest.print_matrix_BIG_info("MicroSys");
	CoupledTest.print_matrix_restrict_info("MicroSys");

	libMesh::PetscMatrix<libMesh::Number>& LumpingTestInput = CoupledTest.get_restrict_coupling_matrix("MicroSys");
	libMesh::PetscMatrix<libMesh::Number> LumpingTestOutput(LumpingTestInput.comm());

	carl::lump_matrix(LumpingTestInput,LumpingTestOutput);

	std::cout << std::endl << "| Testing the lumping " << std::endl;
	carl::print_matrix(LumpingTestOutput);

	std::cout << std::endl << "| --> Testing the solver " << std::endl << std::endl;
	CoupledTest.set_LATIN_solver("MicroSys","Elasticity");

	return 0;
}
