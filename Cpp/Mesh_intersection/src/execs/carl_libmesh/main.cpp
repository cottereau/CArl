#include "main.h"

void generate_intersection_tables_partial(	std::string& intersection_table_restrict_B_Filename,
											std::string& intersection_table_I_Filename,
											std::unordered_map<int,int>& mesh_restrict_ElemMap,
											std::unordered_map<int,int>& mesh_micro_ElemMap,
											std::unordered_map<int,int>& mesh_inter_ElemMap,
											std::vector<std::pair<int,int> >& intersection_table_restrict_A,
											std::vector<std::pair<int,int> >& intersection_table_restrict_B,
											std::unordered_multimap<int,int>& intersection_table_I
											)
{
	std::ifstream table_restrict_B_file(intersection_table_restrict_B_Filename);
	std::ifstream table_I_file(intersection_table_I_Filename);

	int nbOfIntersections_restrict_B = -1;
	int nbOfIntersectionsI  = -1;
	int nbOfTotalTetrasI  = -1;
	int dummyInt = -1;
	int nbOfTetras = -1;
	int tetraIdx = -1;

	int extIdxA = -1;
	int extIdxB = -1;
	int extIdxI = -1;

	table_restrict_B_file >> nbOfIntersections_restrict_B;
	table_I_file >> nbOfIntersectionsI >> nbOfTotalTetrasI;

	homemade_assert_msg(nbOfIntersections_restrict_B == nbOfIntersectionsI, "Incompatible intersection table files!");

	intersection_table_restrict_B.resize(nbOfIntersections_restrict_B);
	intersection_table_I.reserve(2*nbOfTotalTetrasI);

	for(int iii = 0; iii < nbOfIntersections_restrict_B; ++iii)
	{
		table_restrict_B_file >> dummyInt >> extIdxA >> extIdxB;
		intersection_table_restrict_B[iii].first = mesh_restrict_ElemMap[extIdxA];
		intersection_table_restrict_B[iii].second = mesh_micro_ElemMap[extIdxB];

		table_I_file >> dummyInt >> nbOfTetras;
		for(int jjj = 0; jjj < nbOfTetras; ++jjj)
		{
			table_I_file >> extIdxI;
			tetraIdx = mesh_inter_ElemMap[extIdxI];
			intersection_table_I.emplace(dummyInt,tetraIdx);
		}
	}

	table_restrict_B_file.close();
	table_I_file.close();
};

void generate_intersection_tables_full(		std::string& intersection_table_restrict_A_Filename,
											std::string& intersection_table_restrict_B_Filename,
											std::string& intersection_table_I_Filename,
											std::unordered_map<int,int>& mesh_restrict_ElemMap,
											std::unordered_map<int,int>& mesh_micro_ElemMap,
											std::unordered_map<int,int>& mesh_inter_ElemMap,
											std::vector<std::pair<int,int> >& intersection_table_restrict_A,
											std::vector<std::pair<int,int> >& intersection_table_restrict_B,
											std::unordered_multimap<int,int>& intersection_table_I
											)
{
	libmesh_error_msg("\n---> ERROR: I shouldn't be here ...");
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
	 * 		TODO :	Build a (for now) serial coupling matrix assembler
	 *
	 * 		TODO :	Look for how libMesh deals with the matrix repartition
	 * 				between the processors
	 *
	 * 		TODO :	Analyze the costs of inverting a matrix with PETSc and co.
	 * 				(hint: it's huge) and implement an "exact" projection matrix
	 * 				assembler
	 *
	 * 		TODO :	Implement a "lumping" projection matrix assembler
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

	std::string intersection_table_restrict_A_Filename;
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

	std::vector<std::pair<int,int> > intersection_table_restrict_A;
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
		create_mesh_map(mesh_BIG_Filename,mesh_BIG_NodeMap,mesh_BIG_ElemMap);

		command_line.search(2, "--meshB", "-mB");
		mesh_micro_Filename = command_line.next(mesh_micro_Filename);
		libMesh::GmshIO meshBuffer_micro(mesh_micro);
		meshBuffer_micro.read(mesh_micro_Filename);
		mesh_micro.prepare_for_use();
		create_mesh_map(mesh_micro_Filename,mesh_micro_NodeMap,mesh_micro_ElemMap);

		command_line.search(2, "--meshI", "-mI");
		mesh_inter_Filename = command_line.next(mesh_inter_Filename);
		libMesh::GmshIO meshBuffer_inter(mesh_inter);
		meshBuffer_inter.read(mesh_inter_Filename);
		mesh_inter.prepare_for_use();
		create_mesh_map(mesh_inter_Filename,mesh_inter_NodeMap,mesh_inter_ElemMap);

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
			create_mesh_map(mesh_restrict_Filename,mesh_restrict_NodeMap,mesh_restrict_ElemMap);

			command_line.search(1, "--tableA");
			intersection_table_restrict_A_Filename = command_line.next(intersection_table_restrict_A_Filename);

			generate_intersection_tables_full(	intersection_table_restrict_A_Filename,
												intersection_table_restrict_B_Filename,
												intersection_table_I_Filename,
												mesh_restrict_ElemMap,
												mesh_micro_ElemMap,
												mesh_inter_ElemMap,
												intersection_table_restrict_A,
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
			intersection_table_restrict_A_Filename = " --- ";

			generate_intersection_tables_partial(	intersection_table_restrict_B_Filename,
													intersection_table_I_Filename,
													mesh_restrict_ElemMap,
													mesh_micro_ElemMap,
													mesh_inter_ElemMap,
													intersection_table_restrict_A,
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

		std::cout << "| Mesh BIG info :" << std::endl;
		std::cout << "|    filename     " << mesh_BIG_Filename << std::endl;
		std::cout << "|    n_elem       " << mesh_BIG.n_elem() << std::endl;
		std::cout << "|    n_nodes      " << mesh_BIG.n_nodes() << std::endl;
		std::cout << "|    n_subdomains " << mesh_BIG.n_subdomains() << std::endl << std::endl;

		std::cout << "| Mesh micro info :" << std::endl;
		std::cout << "|    filename     " << mesh_micro_Filename << std::endl;
		std::cout << "|    n_elem       " << mesh_micro.n_elem() << std::endl;
		std::cout << "|    n_nodes      " << mesh_micro.n_nodes() << std::endl;
		std::cout << "|    n_subdomains " << mesh_micro.n_subdomains() << std::endl << std::endl;

		std::cout << "| Mesh inter info :" << std::endl;
		std::cout << "|    filename     " << mesh_inter_Filename << std::endl;
		std::cout << "|    n_elem       " << mesh_inter.n_elem() << std::endl;
		std::cout << "|    n_nodes      " << mesh_inter.n_nodes() << std::endl;
		std::cout << "|    n_subdomains " << mesh_inter.n_subdomains() << std::endl << std::endl;

		std::cout << "| Mesh restriction info :" << extra_restrict_info << std::endl;
		std::cout << "|    filename     " << mesh_restrict_Filename << std::endl;
		std::cout << "|    n_elem       " << mesh_restrict.n_elem() << std::endl;
		std::cout << "|    n_nodes      " << mesh_restrict.n_nodes() << std::endl;
		std::cout << "|    n_subdomains " << mesh_restrict.n_subdomains() << std::endl << std::endl;

		std::cout << "| Inter. table restrict / A :" << std::endl;
		std::cout << "|    filename     " << intersection_table_restrict_A_Filename << std::endl << std::endl;

		std::cout << "| Inter. table restrict / B :" << std::endl;
		std::cout << "|    filename     " << intersection_table_restrict_B_Filename << std::endl << std::endl;

		std::cout << "| Inter. table I  :" << std::endl;
		std::cout << "|    filename     " << intersection_table_I_Filename << std::endl << std::endl;

////		DEBUG
//		int sillySize = intersection_table_AB.size();
//		for(int iii = 0; iii < sillySize; ++iii)
//		{
//			std::cout 	<< iii << " | "
//						<< intersection_table_AB[iii].first << " "
//						<< intersection_table_AB[iii].second << std::endl;
//		}
//
//		std::unordered_multimap<int,int>::iterator itInterTetras;
//
//		for(int iii = 0; iii < sillySize; ++iii)
//		{
//			auto sillyIrange = intersection_table_I.equal_range(iii);
//			std::cout 	<< iii << " | ";
//
//			for(	itInterTetras = sillyIrange.first;
//					itInterTetras != sillyIrange.second;
//					++itInterTetras)
//			{
//				std::cout 	<< itInterTetras->second << " ";
//			}
//			std::cout << std::endl;
//		}
//		std::cout << std::endl;
	}
	else
	{
		libmesh_error_msg("\n---> ERROR: Must give the three meshes and the intersection tables!\n");
	}

	// Associate the boundary conditions -> no boundary conditions!

	// Generate the equation systems
	carl::coupled_system CoupledTest;

	libMesh::EquationSystems& equation_systems_BIG =
					CoupledTest.set_BIG_EquationSystem("BigSys", mesh_BIG);
	libMesh::EquationSystems& equation_systems_micro =
					CoupledTest.add_micro_EquationSystem<libMesh::PetscMatrix<libMesh::Number> >("MicroSys", mesh_micro);
	libMesh::EquationSystems& equation_systems_restrict =
					CoupledTest.add_restrict_EquationSystem("RestrictSys", mesh_restrict);
	libMesh::EquationSystems& equation_systems_inter =
					CoupledTest.add_inter_EquationSystem("InterSys", mesh_inter);

	// - Build the BIG system --------------------------------------------------

	libMesh::LinearImplicitSystem& volume_BIG_system =
			equation_systems_BIG.add_system<libMesh::LinearImplicitSystem> ("VolTest");

	unsigned int sillyVar_BIG = volume_BIG_system.add_variable("SillyVar", libMesh::FIRST, libMesh::LAGRANGE);

	equation_systems_BIG.init();

	// - Build the micro system ------------------------------------------------

	libMesh::LinearImplicitSystem& volume_micro_system =
			equation_systems_micro.add_system<libMesh::LinearImplicitSystem> ("VolTest");

	unsigned int sillyVar_micro = volume_micro_system.add_variable("SillyVar", libMesh::FIRST, libMesh::LAGRANGE);

	equation_systems_micro.init();

	// - Build the restrict system --------------------------------------------------

	libMesh::LinearImplicitSystem& volume_restrict_system =
			equation_systems_restrict.add_system<libMesh::LinearImplicitSystem> ("VolTest");

	unsigned int sillyVar_restrict = volume_restrict_system.add_variable("SillyVar", libMesh::FIRST, libMesh::LAGRANGE);

	equation_systems_restrict.init();

	// - Build the dummy inter system ------------------------------------------

	libMesh::LinearImplicitSystem& volume_inter_system =
			equation_systems_inter.add_system<libMesh::LinearImplicitSystem> ("VolTest");

	unsigned int sillyVar_inter = volume_inter_system.add_variable("SillyVar", libMesh::FIRST, libMesh::LAGRANGE);

	equation_systems_inter.init();

	CoupledTest.assemble_coupling_matrices(	"BigSys","MicroSys",
											"InterSys","RestrictSys",
											intersection_table_restrict_B,
											intersection_table_I,
											using_same_mesh_restrict_A);

//	elasticity_system.attach_assemble_function(assemble_elasticity_heterogeneous);
//
//	// Defining the boundaries with Dirichlet conditions ...
//	std::set<boundary_id_type> boundary_ids_clamped;
//	std::set<boundary_id_type> boundary_ids_displaced;
//
//	boundary_ids_clamped.insert(BOUNDARY_ID_MIN_X);
//	boundary_ids_displaced.insert(BOUNDARY_ID_MAX_X);
//
//	std::vector<unsigned int> variables(3);
//	variables[0] = u_var; variables[1] = v_var; variables[2] = w_var;
//
//	ZeroFunction<> zf;
//	border_displacement right_border(u_var,v_var,w_var,
//									 x_max_x_displ,x_max_y_displ,x_max_z_displ);
//
//	// ... and set them
//	DirichletBoundary dirichlet_bc_clamped(	boundary_ids_clamped,
//											variables,
//											&zf);
//
//	elasticity_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc_clamped);
//
//	DirichletBoundary dirichlet_bc_displaced(	boundary_ids_displaced,
//												variables,
//												&right_border);
//
//	elasticity_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc_displaced);

	CoupledTest.print_matrix_micro_info("MicroSys");
	CoupledTest.print_matrix_BIG_info("MicroSys");

	return 0;
}
