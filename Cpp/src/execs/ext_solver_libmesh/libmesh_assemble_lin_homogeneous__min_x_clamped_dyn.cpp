#include "libmesh_assemble_lin_homogeneous.h"
//#include "elasticity_system.h"
/**	\brief Program used to assemble the rigidity matrix and the vectors of a linear, homogeneous elasticity model with a clamped \f$x_{\mbox{Min}}\f$ face.
 * 
 *  Usage: `./libmesh_assemble_lin_homogeneous__min_x_clamped -i [input file]`
 *	
 * The input file is parsed by the get_input_params(GetPot& field_parser, libmesh_assemble_input_params& input_params) function, and it contains the following parameters. 
 *
 *	Required parameters:
 *	  - `Mesh` : path to the mesh.
 *    - `PhysicalParameters` : physical parameters.
 *    - `SystemType` : parameter used to tell the assembler which weight functions must be used. *Values*: `Micro` or `Macro`.
 *	  - `MeshWeight` : path to the mesh defining the domains of the Arlequin weight parameters.
 *    - `WeightIndexes` : path to the indices of the domains of the Arlequin weight parameters.
 *
 *  Optional parameter:
 *    - `OutputBase` or `--output` : base of the output files (including folders). *Default*: `test_system`.
 *
 *  Boolean flags:
 *    - `ExportRBVectors` : build and export the rigid body modes vectors.
 */

int main(int argc, char** argv) {

	// [USER] Fixed boudary
	boundary_id_cube fixed_bound_id = boundary_id_cube::MIN_X;
	
	// --- Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);

	// Do performance log?
	const bool MASTER_bPerfLog_carl_libmesh = true;
	libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);

	// libMesh's C++ / MPI communicator wrapper
	libMesh::Parallel::Communicator& WorldComm = init.comm();

	// Number of processors and processor rank.
	int rank = WorldComm.rank();
	int nodes = WorldComm.size();

	// --- Set up inputs

	// Command line parser
	GetPot command_line(argc, argv);

	// File parser
	GetPot field_parser;

	// If there is an input file, parse it to get the parameters. Else, parse the command line
	std::string input_filename;
	if (command_line.search(2, "--inputfile", "-i")) {
		input_filename = command_line.next(input_filename);
		field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
	} else {
		field_parser = command_line;
	}

	libmesh_assemble_input_params input_params;
	get_input_params(field_parser, input_params);

	// Check libMesh installation dimension
	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// --- Declare the three meshes to be intersected

	// - Parallelized meshes: A, B, mediator and weight
	perf_log.push("Meshes - Parallel","Read files:");
	libMesh::Mesh system_mesh(WorldComm, dim);
	system_mesh.read(input_params.mesh_file);
	system_mesh.prepare_for_use();

	libMesh::Mesh mesh_weight(WorldComm, dim); //mesh_weight: mesh object
	mesh_weight.allow_renumbering(false);
	mesh_weight.read(input_params.mesh_weight_file);
	mesh_weight.prepare_for_use();

	perf_log.pop("Meshes - Parallel","Read files:");

	// --- Generate the equation systems
	perf_log.push("System setup:");

    //=======================================================================
    const libMesh::Real deltat = .0000625;
    unsigned int n_time_steps = 300;
	// Set the equation systems object
	libMesh::EquationSystems equation_systems(system_mesh);
    // Add Newmark system with linear elasticity and physical parameters systems

    ElasticitySystem & elasticity_system =
        equation_systems.add_system<ElasticitySystem> ("Linear Elasticity");
    //ElasticitySystem& elasticity_system = add_dyn_newmark(equation_systems);
    
    //Create ExplicitSystem to help output velocity
    libMesh::ExplicitSystem * v_system;
    v_system = &add_vel_newmark(equation_systems);
    
    // Create ExplicitSystem to help output acceleration
    libMesh::ExplicitSystem * a_system;
    a_system = &add_acc_newmark(equation_systems);

    // Define and reset solver
    elasticity_system.time_solver.reset(new libMesh::NewmarkSolver(elasticity_system)); 

	// Set clamped border
	set_clamped_border(elasticity_system, fixed_bound_id);

	// Initialize the equation systems
	equation_systems.init();

    // Set the time stepping options
    elasticity_system.deltat = deltat;    

	// Homogeneous properties for the macro system
	set_homogeneous_physical_properties_dyn(equation_systems, input_params.physical_params_file);

    // And the nonlinear solver options
    libMesh::DiffSolver& solver = *(elasticity_system.time_solver->diff_solver().get());
    solver.quiet = false;
    solver.verbose = true;
    solver.max_nonlinear_iterations = 15;
    solver.relative_step_tolerance = 1.e-3;
    solver.relative_residual_tolerance = 0.0;
    solver.absolute_residual_tolerance = 0.0;

    // And the linear solver options
    solver.max_linear_iterations = 50000;
    solver.initial_linear_tolerance = 1.e-3;

	// Set the weight function object
	weight_parameter_function  system_weight(mesh_weight);
	system_weight.set_parameters(input_params.weight_domain_idx_file);

	perf_log.pop("System setup:");

	// Assemble!
	assemble_elasticity_with_weight_dyn(equation_systems,"Linear Elasticity",system_weight,
        input_params.system_type);

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
	elasticity_system.matrix->print_matlab(input_params.output_base + "_sys_mat.m");
	elasticity_system.rhs->print_matlab(input_params.output_base + "_sys_rhs_vec.m");
#endif

	// Export matrix and vector
	libMesh::PetscMatrix<libMesh::Number> * temp_mat_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(elasticity_system.matrix);
	libMesh::PetscVector<libMesh::Number> * temp_vec_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> * >(elasticity_system.rhs);

	carl::write_PETSC_matrix(*temp_mat_ptr, input_params.output_base + "_sys_mat.petscmat");
	carl::write_PETSC_vector(*temp_vec_ptr, input_params.output_base + "_sys_rhs_vec.petscvec");

	// If needed, print rigid body vectors
	if(input_params.bCalculateRBVectors)
	{
		MatNullSpace nullsp_sys;
		build_rigid_body_vectors(elasticity_system,nullsp_sys);
		write_rigid_body_vectors(nullsp_sys,input_params.output_base,WorldComm.rank());
		MatNullSpaceDestroy(&nullsp_sys);
	}

	return 0;
}

//ElasticitySystem& add_dyn_newmark(libMesh::EquationSystems& input_systems, 
//    libMesh::Order order, 
//    libMesh::FEFamily family)
//{
//    libMesh::ExplicitSystem& physical_variables =
//        input_systems.add_system<libMesh::ExplicitSystem> ("PhysicalConstants");
//
//    // Physical constants are set as constant, monomial
//    physical_variables.add_variable("E", libMesh::CONSTANT, libMesh::MONOMIAL);
//    physical_variables.add_variable("mu", libMesh::CONSTANT, libMesh::MONOMIAL);
//    physical_variables.add_variable("Rho", libMesh::CONSTANT, libMesh::MONOMIAL);
//    physical_variables.add_variable("Index", libMesh::CONSTANT, libMesh::MONOMIAL);
//
//    ElasticitySystem& elasticity_system =
//        input_systems.add_system<ElasticitySystem> ("Linear Elasticity");
//
//    return elasticity_system;
//}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=4 et tw=80 smartindent :                               */
