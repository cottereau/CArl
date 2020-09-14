#include "libmesh_assemble_lin_homogeneous.h"
/** \brief Program used to assemble the rigidity matrix and the vectors of a linear, homogeneous elasticity model with a clamped \f$x_{\mbox{Min}}\f$ face.
 * 
 *  Usage: `./libmesh_assemble_lin_homogeneous__min_x_clamped_dyn -i [input file]`
 *  
 * The input file is parsed by the get_input_params(GetPot& field_parser, libmesh_assemble_input_params& input_params) function, and it contains the following parameters. 
 *
 *  Required parameters:
 *    - `Mesh` : path to the mesh.
 *    - `PhysicalParameters` : physical parameters.
 *    - `SystemType` : parameter used to tell the assembler which weight functions must be used. *Values*: `Micro` or `Macro`.
 *    - `MeshWeight` : path to the mesh defining the domains of the Arlequin weight parameters.
 *    - `WeightIndexes` : path to the indices of the domains of the Arlequin weight parameters.
 *
 *  Optional parameter:
 *    - `OutputBase` or `--output` : base of the output files (including folders). *Default*: `test_system`.
 *
 *  Boolean flags:
 *    - `ExportRBVectors` : build and export the rigid body modes vectors.
 */

using namespace std;
using namespace libMesh;

int main (int argc, char ** argv)
{

  const boundary_id_type bound_clamped_id = boundary_id_min_x;
  const boundary_id_type bound_force_id   = boundary_id_max_z;

  Real force = 10000.; 

  // Initialize libMesh.
  LibMeshInit init (argc, argv);

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

  // Declaration of a struct variable to catch all inputs from the input file
  libmesh_assemble_input_params input_params; // struct type 

  // Getting parsing argument from all file
  get_input_params(field_parser, input_params);

  // Check libMesh installation dimension
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
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
  perf_log.push("System setup:"); //add to stack of per_log 

  // Create an equation systems object
  EquationSystems equation_systems (system_mesh);

  ElasticitySystem & elasticity_system = equation_systems.add_system<ElasticitySystem> ("Linear Elasticity");
  //LinearImplicitSystem& elasticity_system
  //                  = add_elasticity(equation_systems); // define in common_assemble_functions_elasticity_3D.cpp
  //                  //by default order is FIRST and family is LAGRANGE (definition function in the file 
  //                  //common_assemble_functions_elasticity_3D.h). 
  //                  //See also http://libmesh.github.io/doxygen/namespacelibMesh.html#af3eb3e8751995b944fc135ea53b09da2

  // Set clamped boundary
  elasticity_system.set_dirichlet_bc(system_mesh,bound_clamped_id,bound_force_id,false);
  
  equation_systems.init();

  //elasticity_system.setForce(force);
  
  // Set material parameter
  elasticity_system.set_homogeneous_physical_properties(equation_systems, input_params.physical_params_file);
  
  // Set the weight function object
  weight_parameter_function  system_weight(mesh_weight);
  system_weight.set_parameters(input_params.weight_domain_idx_file);

  perf_log.pop("System setup:");

/*  // Solve this as a time-dependent or steady system
  std::string time_solver = std::string("newmark")/*infile("time_solver","DIE!");

  ExplicitSystem * v_system;
  ExplicitSystem * a_system;

  // Create ExplicitSystem to help output velocity
  v_system = &equation_systems.add_system<ExplicitSystem> ("Velocity");
  v_system->add_variable("u_vel", FIRST, LAGRANGE);
  v_system->add_variable("v_vel", FIRST, LAGRANGE);
  v_system->add_variable("w_vel", FIRST, LAGRANGE);

  // Create ExplicitSystem to help output acceleration
  a_system = &equation_systems.add_system<ExplicitSystem> ("Acceleration");
  a_system->add_variable("u_accel", FIRST, LAGRANGE);
  a_system->add_variable("v_accel", FIRST, LAGRANGE);
  a_system->add_variable("w_accel", FIRST, LAGRANGE);
  
  elasticity_system.time_solver = libmesh_make_unique<NewmarkSolver>(elasticity_system);
  
  // Initialize the system
  //->launch init_data() methode of ElasticitySystem
  equation_systems.init();

  // Set the time stepping options 
  // Read step time
  const Real deltat = static_cast<Real>(input_params.deltat);
  elasticity_system.deltat = deltat;

  // And the nonlinear solver options
  DiffSolver & solver = *(elasticity_system.time_solver->diff_solver().get());
  solver.quiet = input_params.solver_quiet;
  solver.verbose = !solver.quiet;
  solver.max_nonlinear_iterations = input_params.max_nonlinear_iterations;
  solver.relative_step_tolerance = static_cast<Real>(input_params.relative_step_tolerance);
  solver.relative_residual_tolerance = static_cast<Real>(input_params.relative_residual_tolerance);
  solver.absolute_residual_tolerance = static_cast<Real>(input_params.absolute_residual_tolerance);

  // And the linear solver options
  solver.max_linear_iterations = input_params.max_linear_iterations;
  solver.initial_linear_tolerance = static_cast<Real>(input_params.initial_linear_tolerance);
  //unsigned int n_timesteps = input_params.n_timesteps;
  // Print information about the elasticity_system to the screen.
  equation_systems.print_info();

  NewmarkSolver * newmark = cast_ptr<NewmarkSolver*>(elasticity_system.time_solver.get());
  newmark->compute_initial_accel();

  // Copy over initial velocity and acceleration for output.
  // Note we can do this because of the matching variables/FE spaces
  *(v_system->solution) = elasticity_system.get_vector("_old_solution_rate");
  *(a_system->solution) = elasticity_system.get_vector("_old_solution_accel");

  //libMesh::SparseMatrix< libMesh::Number > * MassTilde = 
  //  get_mass_tilde(equation_systems,"Linear Elasticity")
*/
  
  // Export matrix and vector
  libMesh::PetscMatrix<libMesh::Number> * temp_mat_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(elasticity_system.matrix);
  libMesh::PetscVector<libMesh::Number> * temp_vec_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> * >(elasticity_system.rhs);
  
  libMesh::PetscVector<libMesh::Number> * temp_vel_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&elasticity_system.get_vector("_old_solution_rate"));
  libMesh::PetscVector<libMesh::Number> * temp_acc_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&elasticity_system.get_vector("_old_solution_accel"));
  
  carl::write_PETSC_matrix(*temp_mat_ptr, input_params.output_base + "_sys_mat.petscmat");
  carl::write_PETSC_vector(*temp_vec_ptr, input_params.output_base + "_sys_rhs_vec.petscvec");
  carl::write_PETSC_vector(*temp_vel_ptr, input_params.output_base + "_sys_vel_n.petscvec");
  carl::write_PETSC_vector(*temp_acc_ptr, input_params.output_base + "_sys_acc_n.petscvec");
  
  // If needed, print rigid body vectors
  if(input_params.bCalculateRBVectors)
  {
    MatNullSpace nullsp_sys;
    build_rigid_body_vectors(elasticity_system,nullsp_sys);
    write_rigid_body_vectors(nullsp_sys,input_params.output_base,WorldComm.rank());
    MatNullSpaceDestroy(&nullsp_sys);
  }
  // All done.
  //
  return 0;
}
/* Local Variables:                                                        */
/* mode: c                                                              */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=4 et tw=80 smartindent :                               */