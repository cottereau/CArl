#include "elasticity_3D_test.h"

int main (int argc, char** argv)
{
	libMesh::PerfLog perf_log ("Main program");

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

	std::string outputParamsFile;
	if ( command_line.search(1, "-o","--output") )
	{
		outputParamsFile = command_line.next(outputParamsFile);
	}
	else
	{
		outputParamsFile = "displacement_and_stress_hetero.exo";
	}

	// Set meshes
	libMesh::Mesh mesh(init.comm(), dim);
	std::string meshFilename;

	// Read mesh ...
	int BOUNDARY_ID_MIN_Z = 0;
	int BOUNDARY_ID_MIN_Y = 1;
	int BOUNDARY_ID_MAX_X = 2;
	int BOUNDARY_ID_MAX_Y = 3;
	int BOUNDARY_ID_MIN_X = 4;
	int BOUNDARY_ID_MAX_Z = 5;

	if ( command_line.search(2, "--mesh", "-m") )
	{
		meshFilename = command_line.next(meshFilename);
		libMesh::GmshIO meshBuffer(mesh);
		meshBuffer.read(meshFilename);
		mesh.prepare_for_use();

		++BOUNDARY_ID_MIN_Z;
		++BOUNDARY_ID_MIN_Y;
		++BOUNDARY_ID_MAX_X;
		++BOUNDARY_ID_MAX_Y;
		++BOUNDARY_ID_MIN_X;
		++BOUNDARY_ID_MAX_Z;
	}
	else
	{
		// ... or generate it
		libMesh::MeshTools::Generation::build_cube (	mesh,
											40,
											10,
											10,
											0., 3.,
											0., 1.,
											0., 1.,
											libMesh::HEX8);
	}

	// Associate the boundary conditions
	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	for ( ; el != end_el; ++el)
	{
		const libMesh::Elem* elem = *el;

		// Test if it is on max x/y/z and/or min y
		unsigned int side_max_x = 0, side_min_x = 0;

		bool found_side_max_x = false, found_side_min_x = false;

		for(unsigned int side=0; side<elem->n_sides(); side++)
		{
			if( mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_X))
			{
				side_max_x = side;
				found_side_max_x = true;
			}

			if( mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MIN_X))
			{
				side_min_x = side;
				found_side_min_x = true;
			}
		}
	}

	perf_log.push("System initialization");

	libMesh::EquationSystems equation_systems(mesh);

	// - Set up the physical properties ---------------------------------------
	libMesh::ExplicitSystem& physical_variables =
			equation_systems.add_system<libMesh::ExplicitSystem> ("PhysicalConstants");

	physical_variables.add_variable("E", libMesh::CONSTANT, libMesh::MONOMIAL);
	physical_variables.add_variable("mu", libMesh::CONSTANT, libMesh::MONOMIAL);

	// - Build elasticity equation system -------------------------------------

	libMesh::LinearImplicitSystem& elasticity_system =
			equation_systems.add_system<libMesh::LinearImplicitSystem> ("Elasticity");

	unsigned int u_var = elasticity_system.add_variable("u", libMesh::FIRST, libMesh::LAGRANGE);
	unsigned int v_var = elasticity_system.add_variable("v", libMesh::FIRST, libMesh::LAGRANGE);
	unsigned int w_var = elasticity_system.add_variable("w", libMesh::FIRST, libMesh::LAGRANGE);

	elasticity_system.attach_assemble_function(assemble_elasticity_heterogeneous);

	// Defining the boundaries with Dirichlet conditions ...
	std::set<libMesh::boundary_id_type> boundary_ids_clamped;
	std::set<libMesh::boundary_id_type> boundary_ids_displaced;

	boundary_ids_clamped.insert(BOUNDARY_ID_MIN_X);
	boundary_ids_displaced.insert(BOUNDARY_ID_MAX_X);

	std::vector<unsigned int> variables(3);
	variables[0] = u_var; variables[1] = v_var; variables[2] = w_var;

	libMesh::ZeroFunction<> zf;
	border_displacement right_border(u_var,v_var,w_var,
									 x_max_x_displ,x_max_y_displ,x_max_z_displ);

	// ... and set them
	libMesh::DirichletBoundary dirichlet_bc_clamped(	boundary_ids_clamped,
											variables,
											&zf);

	elasticity_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc_clamped);

	libMesh::DirichletBoundary dirichlet_bc_displaced(	boundary_ids_displaced,
												variables,
												&right_border);

	elasticity_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc_displaced);

	// - Build stress system --------------------------------------------------
	libMesh::ExplicitSystem& stress_system =
				equation_systems.add_system<libMesh::ExplicitSystem> ("StressSystem");

	stress_system.add_variable("sigma_00", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_01", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_02", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_10", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_11", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_12", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_20", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_21", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_22", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("vonMises", libMesh::CONSTANT, libMesh::MONOMIAL);

	equation_systems.init();
	perf_log.pop("System initialization");

	perf_log.push("Physical properties");
	set_physical_properties(equation_systems,physicalParamsFile);
	perf_log.pop("Physical properties");

	// Solve
	perf_log.push("Solve");
	elasticity_system.solve();
	perf_log.pop("Solve");

	// Calculate stress
	perf_log.push("Compute stress");
	compute_stresses(equation_systems);
	perf_log.pop("Compute stress");

	#ifdef LIBMESH_HAVE_EXODUS_API

	libMesh::ExodusII_IO exo_io(mesh, /*single_precision=*/true);

	std::set<std::string> system_names;
	system_names.insert("Elasticity");
	exo_io.write_equation_systems(outputParamsFile.c_str(),equation_systems,&system_names);

	exo_io.write_element_data(equation_systems);

	#endif // #ifdef LIBMESH_HAVE_EXODUS_API

	// Print short mesh info
	std::cout << "| Mesh info :" << std::endl;
	std::cout << "|    n_elem       " << mesh.n_elem() << std::endl;
	std::cout << "|    n_nodes      " << mesh.n_nodes() << std::endl;
	std::cout << "|    n_subdomains " << mesh.n_subdomains() << std::endl << std::endl;

	return 0;
}
