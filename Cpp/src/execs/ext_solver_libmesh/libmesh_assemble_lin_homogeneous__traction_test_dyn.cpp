#include "libmesh_assemble_lin_homogeneous.h"
int main(int argc, char** argv) {
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
    
    libMesh::Mesh mesh_weight(WorldComm, dim);
    mesh_weight.allow_renumbering(false);
    mesh_weight.read(input_params.mesh_weight_file);
    mesh_weight.prepare_for_use();
    
    perf_log.pop("Meshes - Parallel","Read files:");
    
    // --- Generate the equation systems
    perf_log.push("System setup:");
    
    //=======================================================================
    const libMesh::Real deltat = .0000625;
    unsigned int n_timesteps = 300;
    // [USER] Traction force density
    std::vector<double> traction_density(3,0);
    traction_density[0] = 100.;
    // Set the equation systems object
    libMesh::EquationSystems equation_systems(system_mesh);
    
    // Add Newmark system with linear elasticity and physical parameters systems
    libMesh::ExplicitSystem& physical_variables =
        equation_systems.add_system<libMesh::ExplicitSystem> ("PhysicalConstants");
    
    physical_variables.add_variable("E", libMesh::CONSTANT, libMesh::MONOMIAL);
    physical_variables.add_variable("mu", libMesh::CONSTANT, libMesh::MONOMIAL);
    physical_variables.add_variable("Index", libMesh::CONSTANT, libMesh::MONOMIAL);
    //
    ElasticitySystem & elasticity_system =
        equation_systems.add_system<ElasticitySystem> ("Linear Elasticity");
    
    //Create ExplicitSystem to help output velocity
    libMesh::ExplicitSystem * v_system;
    v_system = &add_vel_newmark(equation_systems);
    
    // Create ExplicitSystem to help output acceleration
    libMesh::ExplicitSystem * a_system;
    a_system = &add_acc_newmark(equation_systems);
    
    // Define and reset solver
    elasticity_system.time_solver.reset(new libMesh::NewmarkSolver(elasticity_system)); 
    
    // Set the time stepping options
    elasticity_system.deltat = deltat;    
    // Initialize the equation systems+Dirichlet
    // Homogeneous properties for the macro system
    equation_systems.init();
    
    set_homogeneous_physical_properties_dyn(equation_systems, input_params.physical_params_file);
    // Set the weight function object
    weight_parameter_function  system_weight(mesh_weight);
    system_weight.set_parameters(input_params.weight_domain_idx_file);
    
    perf_log.pop("System setup:");
    
    // Assemble!
    assemble_elasticity_with_weight_and_traction_dyn(equation_systems,"Elasticity",system_weight,
        input_params.system_type,boundary_id_cube::MAX_X,traction_density);
    
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
